//! How the energy changes as the nuclei move: the forces that relax a
//! molecule into its shape.
//!
//! Analytic, term by term: the nuclei's repulsion; kinetic energy and
//! attraction to the nuclei, both through the basis functions moving with
//! their atoms and through each nucleus's own pull (Hellmann-Feynman); the
//! overlap term (Pulay) that a moving basis brings, weighted by the
//! energy-weighted density; the density-fitted Coulomb energy through its
//! three- and two-centre integrals; and exchange-correlation through the
//! moving basis functions, the grid's points moving with their atoms, and
//! Becke's partition changing as the atoms move. `tests/electrons.rs` checks
//! each against finite differences of its own energy.

use super::basis::{Basis, Shell};
use super::functional::evaluate;
use super::integrals::{eri_derivative, eri_three_gradient_with, one_electron_gradient, pairs as pairs_of, pairs_ext, Pair, ThreeScratch};
use super::linalg::{product_nt, Matrix};
use super::scf::{parallel_interleaved, Batches, Problem, Solution};

/// Which nucleus each shell sits on, by position.
fn owners(shells: &[Shell], nuclei: &[(f64, [f64; 3])]) -> Vec<usize> {
    shells
        .iter()
        .map(|s| {
            nuclei
                .iter()
                .position(|(_, p)| (0..3).all(|d| (p[d] - s.centre[d]).abs() < 1e-10))
                .expect("every shell sits on a nucleus")
        })
        .collect()
}

/// The energy gradient, hartree per bohr, one vector per nucleus.
pub fn gradient(problem: &Problem, solution: &Solution) -> Vec<[f64; 3]> {
    let parts = gradient_parts(problem, solution);
    let mut g = vec![[0.0; 3]; problem.nuclei.len()];
    for part in &parts {
        for (a, b) in g.iter_mut().zip(part) {
            for k in 0..3 {
                a[k] += b[k];
            }
        }
    }
    g
}

/// The gradient's terms separately: nuclear repulsion, one-electron (with the
/// overlap term), Coulomb, exchange-correlation, non-local correlation (zero
/// without it).
pub fn gradient_parts(problem: &Problem, solution: &Solution) -> [Vec<[f64; 3]>; 5] {
    let basis = &problem.basis;
    let nuc = &problem.nuclei;
    let n = basis.size;
    let owner = owners(&basis.shells, nuc);
    let mut d = solution.density_alpha.clone();
    let mut w = solution.weighted_alpha.clone();
    for k in 0..n * n {
        d.a[k] += solution.density_beta.a[k];
        w.a[k] += solution.weighted_beta.a[k];
    }
    let mut g = vec![[0.0; 3]; nuc.len()];
    // Nuclear repulsion.
    for i in 0..nuc.len() {
        for j in 0..nuc.len() {
            if i == j {
                continue;
            }
            let r: Vec<f64> = (0..3).map(|k| nuc[i].1[k] - nuc[j].1[k]).collect();
            let r3 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).powf(1.5);
            for k in 0..3 {
                g[i][k] -= nuc[i].0 * nuc[j].0 * r[k] / r3;
            }
        }
    }
    let nl = match problem.nonlocal {
        Some(_) => nonlocal_gradient(problem, &owner, &d),
        None => vec![[0.0; 3]; nuc.len()],
    };
    [g, one_electron_gradient(basis, nuc, &owner, &d, &w), coulomb_gradient(problem, &owner, &d, &solution.fitted), xc_gradient(problem, &owner, solution), nl]
}

/// The density-fitted Coulomb energy's gradient: `c . d(mn|P) D - 1/2 c . dV c`,
/// with `c` the coefficients the solve already fitted.
fn coulomb_gradient(problem: &Problem, owner: &[usize], d: &Matrix, c: &[f64]) -> Vec<[f64; 3]> {
    let basis = &problem.basis;
    let nuc = &problem.nuclei;
    let aux = problem.auxiliary.as_ref().expect("the gradient is for density-fitted Coulomb");
    assert_eq!(c.len(), aux.size, "the solution carries no fit for this auxiliary set");
    let aux_owner = owners(&aux.shells, nuc);
    let unit = Shell::unit();
    let shells = &basis.shells;
    let n = basis.size;
    let pairs: Vec<(usize, usize)> = (0..shells.len()).flat_map(|a| (0..=a).map(move |b| (a, b))).filter(|&(a, b)| {
        let (sa, sb) = (&shells[a], &shells[b]);
        let r2: f64 = (0..3).map(|k| (sa.centre[k] - sb.centre[k]).powi(2)).sum();
        sa.exponents.iter().any(|x| sb.exponents.iter().any(|y| (-(x * y) / (x + y) * r2).exp() > 1e-14))
    }).collect();
    let kets: Vec<Vec<Pair>> = aux.shells.iter().map(|p| pairs_of(p, &unit, 0)).collect();
    let job = |indices: &mut dyn Iterator<Item = usize>| {
        let mut g = vec![[0.0; 3]; nuc.len()];
        // One worker's working space, kept across its calls (see
        // `ThreeScratch`): allocating per call made the threads queue.
        let mut scratch = ThreeScratch::new();
        for (a, b) in indices.map(|i| pairs[i]) {
            let (sa, sb) = (&shells[a], &shells[b]);
            let f = if a == b { 1.0 } else { 2.0 };
            let (na_, nb_) = (sa.size(), sb.size());
            let mut dm = vec![0.0; na_ * nb_];
            for i in 0..na_ {
                for j in 0..nb_ {
                    dm[i * nb_ + j] = d.a[(basis.offsets[a] + i) * n + basis.offsets[b] + j];
                }
            }
            if dm.iter().all(|x| *x == 0.0) {
                continue;
            }
            let bra = pairs_ext(sa, sb, 1, 1);
            for (ip, p) in aux.shells.iter().enumerate() {
                // All three centres on one atom: the three derivatives sum to
                // zero (moving the atom moves all of them together), so the
                // atom's net is nothing. These are the tightest and most
                // expensive quartets there are.
                if owner[a] == owner[b] && owner[a] == aux_owner[ip] {
                    continue;
                }
                let cp = &c[aux.offsets[ip]..aux.offsets[ip] + p.size()];
                let (ta, tb) = eri_three_gradient_with(&bra, &kets[ip], sa.l, sb.l, p.l, &dm, cp, &mut scratch);
                for dir in 0..3 {
                    g[owner[a]][dir] += f * ta[dir];
                    g[owner[b]][dir] += f * tb[dir];
                    g[aux_owner[ip]][dir] -= f * (ta[dir] + tb[dir]);
                }
            }
        }
        g
    };
    let mut g = vec![[0.0; 3]; nuc.len()];
    for part in parallel_interleaved(pairs.len(), &job) {
        for (a, b) in g.iter_mut().zip(&part) {
            for k in 0..3 {
                a[k] += b[k];
            }
        }
    }
    // -1/2 c . dV c, with d(P|Q)/dQ = -d(P|Q)/dP.
    for (ip, p) in aux.shells.iter().enumerate() {
        for (iq, q) in aux.shells.iter().enumerate() {
            if aux_owner[ip] == aux_owner[iq] {
                continue;
            }
            let dp = eri_derivative(p, &unit, q, &unit);
            for dir in 0..3 {
                let mut t = 0.0;
                for i in 0..p.size() {
                    for j in 0..q.size() {
                        t += c[aux.offsets[ip] + i] * c[aux.offsets[iq] + j] * dp[dir][i * q.size() + j];
                    }
                }
                g[aux_owner[ip]][dir] -= 0.5 * t;
                g[aux_owner[iq]][dir] += 0.5 * t;
            }
        }
    }
    g
}

/// Exchange-correlation's gradient through the moving basis functions:
/// `-2 sum_points w sum_spin sum_(mu on A) [v_rho d_d phi_mu X_mu +
/// W . (grad d_d phi_mu X_mu + d_d phi_mu Y_mu)]`, with `X = D phi` and
/// `Y = D grad phi`.
fn xc_gradient(problem: &Problem, owner: &[usize], solution: &Solution) -> Vec<[f64; 3]> {
    let basis = &problem.basis;
    let nuc = &problem.nuclei;
    let n = basis.size;
    let atoms: Vec<([f64; 3], f64)> = nuc.iter().zip(&problem.sizes).map(|((_, p), r)| (*p, *r)).collect();
    let grid = super::grid::molecular_pruned(&atoms, problem.radial, problem.theta, problem.prune);
    let batches = Batches::new(basis, &grid);
    let f = problem.functional;
    let grad = f.needs_gradient();
    let (da, db) = (&solution.density_alpha, &solution.density_beta);
    // Function index -> owning nucleus.
    let mut fowner = vec![0usize; n];
    for (is, sh) in basis.shells.iter().enumerate() {
        for k in 0..sh.size() {
            fowner[basis.offsets[is] + k] = owner[is];
        }
    }
    let work = &batches.batches;
    let positions: Vec<[f64; 3]> = nuc.iter().map(|(_, p)| *p).collect();
    let same = da.a == db.a;
    let spins = if same { 1 } else { 2 };
    // hessian index of (e, dir)
    let hidx = |a: usize, b: usize| match (a.min(b), a.max(b)) {
        (0, 0) => 0,
        (1, 1) => 1,
        (2, 2) => 2,
        (0, 1) => 3,
        (0, 2) => 4,
        _ => 5,
    };
    let job = |batch_ids: &mut dyn Iterator<Item = usize>| {
        let mut g = vec![[0.0; 3]; nuc.len()];
        for b in batch_ids {
            let (pts, ws, funcs, shells) = &work[b];
            let idx = &batches.indices[b];
            let k = funcs.len();
            let np = pts.len();
            if k == 0 {
                continue;
            }
            // Values, gradients and second derivatives over the whole batch,
            // laid out point by function so that `X = Phi D` and
            // `Y_e = (d_e Phi) D` are four dense products.
            let mut val = vec![0.0; np * k];
            let mut gr = [vec![0.0; np * k], vec![0.0; np * k], vec![0.0; np * k]];
            let mut he: Vec<Vec<f64>> = vec![vec![0.0; np * k]; 6];
            let (mut v1, mut g1, mut h1) = (vec![0.0; k], vec![0.0; 3 * k], vec![0.0; 6 * k]);
            for (q, pt) in pts.iter().enumerate() {
                super::values::at_shells_hessian(basis, shells, *pt, &mut v1, &mut g1, &mut h1);
                for i in 0..k {
                    val[q * k + i] = v1[i];
                    for e in 0..3 {
                        gr[e][q * k + i] = g1[3 * i + e];
                    }
                    for e in 0..6 {
                        he[e][q * k + i] = h1[6 * i + e];
                    }
                }
            }
            let mut x: Vec<Vec<f64>> = Vec::new();
            let mut y: Vec<[Vec<f64>; 3]> = Vec::new();
            for dm in [da, db].iter().take(spins) {
                let mut sub = vec![0.0; k * k];
                for (i, &fi) in funcs.iter().enumerate() {
                    for (j, &fj) in funcs.iter().enumerate() {
                        sub[i * k + j] = dm.a[fi * n + fj];
                    }
                }
                let mut xs = vec![0.0; np * k];
                product_nt(&val, &sub, np, k, k, &mut xs);
                let mut ys = [vec![0.0; np * k], vec![0.0; np * k], vec![0.0; np * k]];
                for e in 0..3 {
                    product_nt(&gr[e], &sub, np, k, k, &mut ys[e]);
                }
                x.push(xs);
                y.push(ys);
            }
            for (q, (pt, &wt)) in pts.iter().zip(ws).enumerate() {
                let gi = idx[q];
                let r = q * k..q * k + k;
                let sp = |s: usize| s.min(spins - 1);
                let rho: Vec<f64> = (0..2).map(|s| val[r.clone()].iter().zip(&x[sp(s)][r.clone()]).map(|(a, b)| a * b).sum()).collect();
                if rho[0] + rho[1] < 1e-14 {
                    continue;
                }
                let gvec: Vec<[f64; 3]> = (0..2).map(|s| {
                    let mut v = [0.0; 3];
                    for e in 0..3 {
                        v[e] = 2.0 * gr[e][r.clone()].iter().zip(&x[sp(s)][r.clone()]).map(|(a, b)| a * b).sum::<f64>();
                    }
                    v
                }).collect();
                let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
                let (e_here, de) = evaluate(f, rho[0], rho[1], dot(gvec[0], gvec[0]), dot(gvec[0], gvec[1]), dot(gvec[1], gvec[1]));
                let wv = [
                    [2.0 * de[2] * gvec[0][0] + de[3] * gvec[1][0], 2.0 * de[2] * gvec[0][1] + de[3] * gvec[1][1], 2.0 * de[2] * gvec[0][2] + de[3] * gvec[1][2]],
                    [2.0 * de[4] * gvec[1][0] + de[3] * gvec[0][0], 2.0 * de[4] * gvec[1][1] + de[3] * gvec[0][1], 2.0 * de[4] * gvec[1][2] + de[3] * gvec[0][2]],
                ];
                let vr = [de[0], de[1]];
                // The basis functions moving with their atoms, at a fixed point.
                let mut moved = [0.0; 3];
                for i in 0..k {
                    let atom = fowner[funcs[i]];
                    let qi = q * k + i;
                    for dir in 0..3 {
                        let mut t = 0.0;
                        for s in 0..2 {
                            let xs = x[sp(s)][qi];
                            t += vr[s] * gr[dir][qi] * xs;
                            if grad {
                                for e in 0..3 {
                                    t += wv[s][e] * (he[hidx(e, dir)][qi] * xs + gr[dir][qi] * y[sp(s)][e][qi]);
                                }
                            }
                        }
                        g[atom][dir] -= 2.0 * wt * t;
                        moved[dir] -= 2.0 * wt * t;
                    }
                }
                // The point moving with its own atom: moving the point is moving
                // every function the other way, so it is minus the sum above.
                let own = batches_owner(&grid, gi);
                for dir in 0..3 {
                    g[own][dir] -= moved[dir];
                }
                // The partition itself moving: Becke's share of this point
                // changes as any atom moves.
                let share_grad = super::grid::becke_share_gradient(&positions, own, *pt);
                for (atom, sg) in share_grad.iter().enumerate() {
                    for dir in 0..3 {
                        g[atom][dir] += grid.raw[gi] * sg[dir] * e_here;
                    }
                }
            }
        }
        g
    };
    let mut g = vec![[0.0; 3]; nuc.len()];
    for part in parallel_interleaved(work.len(), &job) {
        for (a, b) in g.iter_mut().zip(&part) {
            for k in 0..3 {
                a[k] += b[k];
            }
        }
    }
    g
}

fn batches_owner(grid: &super::grid::Grid, i: usize) -> usize {
    grid.owner[i]
}

/// Basis functions of one nucleus — kept for the geometry module.
pub fn shells_on(basis: &Basis, nuclei: &[(f64, [f64; 3])], atom: usize) -> Vec<usize> {
    owners(&basis.shells, nuclei).iter().enumerate().filter(|(_, o)| **o == atom).map(|(i, _)| i).collect()
}

/// Non-local correlation's gradient, on the grid the field was solved with
/// (PLAY.md E7: forces come from the coarse grid). Four ways a nucleus moves
/// the energy: its basis functions move, changing the density and its
/// gradient at every fixed point (the semilocal formula, with the non-local
/// term's `dE/dn` and `dE/dgrad n`); its own points move with it, which at a
/// fixed density is minus that sum; its points move relative to every other
/// atom's, which changes the distances the kernel reads; and Becke's partition
/// of every point changes, through `dE/dw`.
fn nonlocal_gradient(problem: &Problem, owner: &[usize], d: &Matrix) -> Vec<[f64; 3]> {
    let spec = problem.nonlocal.expect("non-local correlation");
    let basis = &problem.basis;
    let nuc = &problem.nuclei;
    let n = basis.size;
    let atoms: Vec<([f64; 3], f64)> = nuc.iter().zip(&problem.sizes).map(|((_, p), r)| (*p, *r)).collect();
    let grid = super::grid::molecular_pruned(&atoms, spec.radial, spec.theta, false);
    let batches = Batches::new(basis, &grid);
    let dens = super::vdw::density_and_gradient_on(basis, &batches, d, grid.points.len());
    let forces = super::vdw::nonlocal_for_forces(&grid, &dens, spec.z_ab, spec.floor, super::vdw::kernel_table());
    let mut fowner = vec![0usize; n];
    for (is, sh) in basis.shells.iter().enumerate() {
        for k in 0..sh.size() {
            fowner[basis.offsets[is] + k] = owner[is];
        }
    }
    let positions: Vec<[f64; 3]> = nuc.iter().map(|(_, p)| *p).collect();
    let hidx = |a: usize, b: usize| match (a.min(b), a.max(b)) {
        (0, 0) => 0,
        (1, 1) => 1,
        (2, 2) => 2,
        (0, 1) => 3,
        (0, 2) => 4,
        _ => 5,
    };
    let work = &batches.batches;
    let job = |batch_ids: &mut dyn Iterator<Item = usize>| {
        let mut g = vec![[0.0; 3]; nuc.len()];
        for b in batch_ids {
            let (pts, ws, funcs, shells) = &work[b];
            let idx = &batches.indices[b];
            let k = funcs.len();
            if k == 0 {
                continue;
            }
            let mut sub = vec![0.0; k * k];
            for (i, &fi) in funcs.iter().enumerate() {
                for (j, &fj) in funcs.iter().enumerate() {
                    sub[i * k + j] = d.a[fi * n + fj];
                }
            }
            let (mut v1, mut g1, mut h1) = (vec![0.0; k], vec![0.0; 3 * k], vec![0.0; 6 * k]);
            let (mut x, mut y) = (vec![0.0; k], vec![[0.0; 3]; k]);
            for (q, (pt, &wt)) in pts.iter().zip(ws).enumerate() {
                let gi = idx[q];
                let (vn, vg) = (forces.nl.v_n[gi], forces.nl.v_g2[gi]);
                let own = grid.owner[gi];
                // The points moving apart, and the partition moving.
                for dir in 0..3 {
                    g[own][dir] += forces.de_dr[gi][dir];
                }
                if forces.de_dw[gi] != 0.0 {
                    let share = super::grid::becke_share_gradient(&positions, own, *pt);
                    for (atom, sg) in share.iter().enumerate() {
                        for dir in 0..3 {
                            g[atom][dir] += grid.raw[gi] * sg[dir] * forces.de_dw[gi];
                        }
                    }
                }
                if vn == 0.0 && vg == 0.0 {
                    continue;
                }
                super::values::at_shells_hessian(basis, shells, *pt, &mut v1, &mut g1, &mut h1);
                for i in 0..k {
                    x[i] = (0..k).map(|j| sub[i * k + j] * v1[j]).sum();
                    for e in 0..3 {
                        y[i][e] = (0..k).map(|j| sub[i * k + j] * g1[3 * j + e]).sum();
                    }
                }
                let gn = dens[gi].1;
                let wv = [2.0 * vg * gn[0], 2.0 * vg * gn[1], 2.0 * vg * gn[2]];
                // The basis functions moving with their atoms, at a fixed point.
                let mut moved = [0.0; 3];
                for i in 0..k {
                    let atom = fowner[funcs[i]];
                    for dir in 0..3 {
                        let mut t = vn * g1[3 * i + dir] * x[i];
                        for e in 0..3 {
                            t += wv[e] * (h1[6 * i + hidx(e, dir)] * x[i] + g1[3 * i + dir] * y[i][e]);
                        }
                        g[atom][dir] -= 2.0 * wt * t;
                        moved[dir] -= 2.0 * wt * t;
                    }
                }
                // The point moving with its own atom, at the density's values.
                for dir in 0..3 {
                    g[own][dir] -= moved[dir];
                }
            }
        }
        g
    };
    let mut g = vec![[0.0; 3]; nuc.len()];
    for part in parallel_interleaved(work.len(), &job) {
        for (a, b) in g.iter_mut().zip(&part) {
            for k in 0..3 {
                a[k] += b[k];
            }
        }
    }
    g
}
