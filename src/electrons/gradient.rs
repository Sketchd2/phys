//! How the energy changes as the nuclei move: the forces that relax a
//! molecule into its shape.
//!
//! Analytic, term by term: the nuclei's repulsion; kinetic energy and
//! attraction to the nuclei, both through the basis functions moving with
//! their atoms and through each nucleus's own pull (Hellmann-Feynman); the
//! overlap term (Pulay) that a moving basis brings, weighted by the
//! energy-weighted density; the density-fitted Coulomb energy through its
//! three- and two-centre integrals; and exchange-correlation through the
//! moving basis functions. What is left out is the grid's own motion — the
//! integration weights depend on where the nuclei are — and
//! `tests/electrons.rs` measures how much that leaves out against finite
//! differences of the energy, rather than assuming it is small.

use super::basis::{Basis, Shell};
use super::functional::evaluate;
use super::integrals::{eri_derivative, one_electron_gradient};
use super::linalg::Matrix;
use super::scf::{parallel, Batches, Fitted, Problem, Solution};

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
/// overlap term), Coulomb, exchange-correlation.
pub fn gradient_parts(problem: &Problem, solution: &Solution) -> [Vec<[f64; 3]>; 4] {
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
    [g, one_electron_gradient(basis, nuc, &owner, &d, &w), coulomb_gradient(problem, &owner, &d), xc_gradient(problem, &owner, solution)]
}

/// The density-fitted Coulomb energy's gradient: `c . d(mn|P) D - 1/2 c . dV c`.
fn coulomb_gradient(problem: &Problem, owner: &[usize], d: &Matrix) -> Vec<[f64; 3]> {
    let basis = &problem.basis;
    let nuc = &problem.nuclei;
    let aux = problem.auxiliary.as_ref().expect("the gradient is for density-fitted Coulomb");
    let aux_owner = owners(&aux.shells, nuc);
    let (s, _, _) = super::integrals::one_electron(basis, &[]);
    let c = Fitted::new(basis, aux, &s).coefficients(d);
    let unit = Shell::unit();
    let shells = &basis.shells;
    let n = basis.size;
    let pairs: Vec<(usize, usize)> = (0..shells.len()).flat_map(|a| (0..=a).map(move |b| (a, b))).filter(|&(a, b)| {
        let (sa, sb) = (&shells[a], &shells[b]);
        let r2: f64 = (0..3).map(|k| (sa.centre[k] - sb.centre[k]).powi(2)).sum();
        sa.exponents.iter().any(|x| sb.exponents.iter().any(|y| (-(x * y) / (x + y) * r2).exp() > 1e-14))
    }).collect();
    let job = |range: std::ops::Range<usize>| {
        let mut g = vec![[0.0; 3]; nuc.len()];
        for &(a, b) in &pairs[range] {
            let (sa, sb) = (&shells[a], &shells[b]);
            let f = if a == b { 1.0 } else { 2.0 };
            for (ip, p) in aux.shells.iter().enumerate() {
                let da = eri_derivative(sa, sb, p, &unit);
                let db = eri_derivative(sb, sa, p, &unit);
                let (na_, nb_, np_) = (sa.size(), sb.size(), p.size());
                for dir in 0..3 {
                    let (mut ta, mut tb) = (0.0, 0.0);
                    for i in 0..na_ {
                        for j in 0..nb_ {
                            let dm = d.a[(basis.offsets[a] + i) * n + basis.offsets[b] + j];
                            if dm == 0.0 {
                                continue;
                            }
                            for k in 0..np_ {
                                let ck = c[aux.offsets[ip] + k];
                                ta += dm * ck * da[dir][(i * nb_ + j) * np_ + k];
                                tb += dm * ck * db[dir][(j * na_ + i) * np_ + k];
                            }
                        }
                    }
                    g[owner[a]][dir] += f * ta;
                    g[owner[b]][dir] += f * tb;
                    g[aux_owner[ip]][dir] -= f * (ta + tb);
                }
            }
        }
        g
    };
    let mut g = vec![[0.0; 3]; nuc.len()];
    for part in parallel(pairs.len(), &job) {
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
    let job = |range: std::ops::Range<usize>| {
        let mut g = vec![[0.0; 3]; nuc.len()];
        for b in range {
            let (pts, ws, funcs, shells) = &work[b];
            let idx = &batches.indices[b];
            let k = funcs.len();
            if k == 0 {
                continue;
            }
            let mut val = vec![0.0; k];
            let mut gr = vec![0.0; 3 * k];
            let mut he = vec![0.0; 6 * k];
            for (q, (pt, &wt)) in pts.iter().zip(ws).enumerate() {
                let gi = idx[q];
                super::values::at_shells_hessian(basis, shells, *pt, &mut val, &mut gr, &mut he);
                let mut x = [vec![0.0; k], vec![0.0; k]];
                let mut y = [vec![0.0; 3 * k], vec![0.0; 3 * k]];
                for (s, dm) in [da, db].iter().enumerate() {
                    for (i, &fi) in funcs.iter().enumerate() {
                        let row = &dm.a[fi * n..fi * n + n];
                        let (mut xv, mut yv) = (0.0, [0.0; 3]);
                        for (j, &fj) in funcs.iter().enumerate() {
                            let dij = row[fj];
                            xv += dij * val[j];
                            for e in 0..3 {
                                yv[e] += dij * gr[3 * j + e];
                            }
                        }
                        x[s][i] = xv;
                        for e in 0..3 {
                            y[s][3 * i + e] = yv[e];
                        }
                    }
                }
                let rho: Vec<f64> = (0..2).map(|s| (0..k).map(|i| val[i] * x[s][i]).sum()).collect();
                if rho[0] + rho[1] < 1e-14 {
                    continue;
                }
                let gvec: Vec<[f64; 3]> = (0..2).map(|s| {
                    let mut v = [0.0; 3];
                    for i in 0..k {
                        for e in 0..3 {
                            v[e] += 2.0 * gr[3 * i + e] * x[s][i];
                        }
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
                // hessian index of (e, dir)
                let hidx = |a: usize, b: usize| match (a.min(b), a.max(b)) {
                    (0, 0) => 0,
                    (1, 1) => 1,
                    (2, 2) => 2,
                    (0, 1) => 3,
                    (0, 2) => 4,
                    _ => 5,
                };
                // The basis functions moving with their atoms, at a fixed point.
                let mut moved = [0.0; 3];
                for i in 0..k {
                    let atom = fowner[funcs[i]];
                    for dir in 0..3 {
                        let mut t = 0.0;
                        for s in 0..2 {
                            t += vr[s] * gr[3 * i + dir] * x[s][i];
                            if grad {
                                for e in 0..3 {
                                    t += wv[s][e] * (he[6 * i + hidx(e, dir)] * x[s][i] + gr[3 * i + dir] * y[s][3 * i + e]);
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
    for part in parallel(work.len(), &job) {
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
