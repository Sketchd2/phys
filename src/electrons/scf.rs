//! Kohn-Sham density-functional theory, solved self-consistently.
//!
//! The electrons of a fixed set of nuclei, as two spin densities, each the
//! occupied orbitals of an effective one-electron Hamiltonian that depends on
//! the densities themselves: kinetic energy, the nuclei, the electrons'
//! Coulomb field, and the exchange-correlation potential of
//! [`super::functional`]. Iterated until the densities reproduce themselves,
//! with Pulay's DIIS to get there.
//!
//! **Occupation is fractional across degenerate levels.** An atom with a
//! partly filled `p` shell has three equal `p` levels and fewer electrons than
//! places; putting them in whole would break the atom's spherical symmetry by
//! choice of axis. Shared equally, the density stays spherical — the
//! convention every published table of atomic energies uses, and the only one
//! that does not depend on which way the axes point.

use super::basis::Basis;
use super::functional::{evaluate, Functional};
use super::grid::{molecular, Grid};
use super::integrals::{eri_block, one_electron};
use super::linalg::{generalised, orthogonaliser, Matrix};
use super::values::at;

/// What to solve.
#[derive(Debug, Clone)]
pub struct Problem {
    pub basis: Basis,
    /// `(charge, position)`, atomic units.
    pub nuclei: Vec<(f64, [f64; 3])>,
    /// Size of each nucleus's grid region, bohr: the radial map's scale.
    pub sizes: Vec<f64>,
    pub alpha: f64,
    pub beta: f64,
    pub functional: Functional,
    pub radial: usize,
    pub theta: usize,
}

/// What came out.
#[derive(Debug, Clone)]
pub struct Solution {
    pub energy: f64,
    pub converged: bool,
    pub iterations: usize,
    pub levels_alpha: Vec<f64>,
    pub levels_beta: Vec<f64>,
    pub density_alpha: Matrix,
    pub density_beta: Matrix,
    /// The parts of the energy, hartree.
    pub kinetic_and_nuclear: f64,
    pub coulomb: f64,
    pub exchange_correlation: f64,
    pub nuclear_repulsion: f64,
}

/// Electron repulsion, packed by its eightfold symmetry: `(ij|kl)` with
/// `i >= j`, `k >= l`, `ij >= kl`.
pub struct Repulsion {
    n: usize,
    values: Vec<f64>,
}

#[inline]
fn pair(i: usize, j: usize) -> usize {
    if i >= j {
        i * (i + 1) / 2 + j
    } else {
        j * (j + 1) / 2 + i
    }
}

impl Repulsion {
    pub fn new(basis: &Basis) -> Repulsion {
        let n = basis.size;
        let np = n * (n + 1) / 2;
        let mut values = vec![0.0; np * (np + 1) / 2];
        let sh = &basis.shells;
        let off = &basis.offsets;
        for a in 0..sh.len() {
            for b in 0..=a {
                for c in 0..sh.len() {
                    for d in 0..=c {
                        // shell-pair ordering: (a,b) >= (c,d)
                        if (c, d) > (a, b) {
                            continue;
                        }
                        let block = eri_block(&sh[a], &sh[b], &sh[c], &sh[d]);
                        let (na, nb, nc, nd) = (sh[a].size(), sh[b].size(), sh[c].size(), sh[d].size());
                        let mut o = 0;
                        for i in 0..na {
                            for j in 0..nb {
                                for k in 0..nc {
                                    for l in 0..nd {
                                        let v = block[o];
                                        o += 1;
                                        let (fi, fj, fk, fl) = (off[a] + i, off[b] + j, off[c] + k, off[d] + l);
                                        let (p, q) = (pair(fi, fj), pair(fk, fl));
                                        values[pair(p, q)] = v;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Repulsion { n, values }
    }

    /// `J[D]_ij = sum_kl (ij|kl) D_kl`.
    pub fn coulomb(&self, d: &Matrix) -> Matrix {
        let n = self.n;
        let np = n * (n + 1) / 2;
        // Packed density with off-diagonal pairs doubled.
        let mut dp = vec![0.0; np];
        let mut idx = vec![(0usize, 0usize); np];
        for i in 0..n {
            for j in 0..=i {
                let p = pair(i, j);
                dp[p] = if i == j { d.get(i, i) } else { d.get(i, j) + d.get(j, i) };
                idx[p] = (i, j);
            }
        }
        let mut jp = vec![0.0; np];
        let mut k = 0;
        for p in 0..np {
            for q in 0..=p {
                let v = self.values[k];
                k += 1;
                jp[p] += v * dp[q];
                if p != q {
                    jp[q] += v * dp[p];
                }
            }
        }
        let mut j = Matrix::zeros(n);
        for (p, &(a, b)) in idx.iter().enumerate() {
            j.set(a, b, jp[p]);
            j.set(b, a, jp[p]);
        }
        j
    }
}

/// Occupation numbers for `count` electrons over ascending `levels`, shared
/// equally inside a degenerate group that is not completely filled.
pub fn occupy(levels: &[f64], count: f64) -> Vec<f64> {
    let mut occ = vec![0.0; levels.len()];
    let mut left = count;
    let mut i = 0;
    while i < levels.len() && left > 1e-12 {
        let tol = 1e-5 * levels[i].abs().max(1.0);
        let mut j = i + 1;
        while j < levels.len() && (levels[j] - levels[i]).abs() < tol {
            j += 1;
        }
        let group = (j - i) as f64;
        let each = (left / group).min(1.0);
        for o in &mut occ[i..j] {
            *o = each;
        }
        left -= each * group;
        i = j;
    }
    occ
}

fn density(c: &[f64], m: usize, n: usize, occ: &[f64]) -> Matrix {
    let mut d = Matrix::zeros(n);
    for (k, &o) in occ.iter().enumerate() {
        if o == 0.0 {
            continue;
        }
        for i in 0..n {
            let ci = c[i * m + k] * o;
            if ci == 0.0 {
                continue;
            }
            for j in 0..n {
                d.a[i * n + j] += ci * c[j * m + k];
            }
        }
    }
    d
}

/// Exchange-correlation energy and the two spin potentials' matrices.
pub fn exchange_correlation(basis: &Basis, grid: &Grid, f: Functional, da: &Matrix, db: &Matrix) -> (f64, Matrix, Matrix) {
    let n = basis.size;
    let grad = f.needs_gradient();
    let mut va = Matrix::zeros(n);
    let mut vb = Matrix::zeros(n);
    let mut exc = 0.0;
    let mut phi = vec![0.0; n];
    let mut dphi = vec![0.0; 3 * n];
    let mut xa = vec![0.0; n];
    let mut xb = vec![0.0; n];
    let mut aa = vec![0.0; n];
    let mut ab = vec![0.0; n];
    for (p, &w) in grid.points.iter().zip(&grid.weights) {
        at(basis, *p, &mut phi, if grad { Some(&mut dphi) } else { None });
        // X = D phi
        for i in 0..n {
            let (mut sa, mut sb) = (0.0, 0.0);
            for j in 0..n {
                sa += da.a[i * n + j] * phi[j];
                sb += db.a[i * n + j] * phi[j];
            }
            xa[i] = sa;
            xb[i] = sb;
        }
        let mut ra = 0.0;
        let mut rb = 0.0;
        let mut ga = [0.0; 3];
        let mut gb = [0.0; 3];
        for i in 0..n {
            ra += phi[i] * xa[i];
            rb += phi[i] * xb[i];
            if grad {
                for d in 0..3 {
                    ga[d] += 2.0 * dphi[3 * i + d] * xa[i];
                    gb[d] += 2.0 * dphi[3 * i + d] * xb[i];
                }
            }
        }
        if ra + rb < 1e-14 {
            continue;
        }
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        let (saa, sab, sbb) = (dot(ga, ga), dot(ga, gb), dot(gb, gb));
        let (e, de) = evaluate(f, ra, rb, saa, sab, sbb);
        exc += w * e;
        // A_mu = w (1/2 e_rho phi_mu + W . grad phi_mu); V += A phi^T + phi A^T.
        let wa = [2.0 * de[2] * ga[0] + de[3] * gb[0], 2.0 * de[2] * ga[1] + de[3] * gb[1], 2.0 * de[2] * ga[2] + de[3] * gb[2]];
        let wb = [2.0 * de[4] * gb[0] + de[3] * ga[0], 2.0 * de[4] * gb[1] + de[3] * ga[1], 2.0 * de[4] * gb[2] + de[3] * ga[2]];
        for i in 0..n {
            let mut a = 0.5 * de[0] * phi[i];
            let mut b = 0.5 * de[1] * phi[i];
            if grad {
                let g = [dphi[3 * i], dphi[3 * i + 1], dphi[3 * i + 2]];
                a += dot(wa, g);
                b += dot(wb, g);
            }
            aa[i] = w * a;
            ab[i] = w * b;
        }
        for i in 0..n {
            let (ai, bi, pi) = (aa[i], ab[i], phi[i]);
            let rowa = &mut va.a[i * n..i * n + n];
            for j in 0..n {
                rowa[j] += ai * phi[j] + pi * aa[j];
            }
            let rowb = &mut vb.a[i * n..i * n + n];
            for j in 0..n {
                rowb[j] += bi * phi[j] + pi * ab[j];
            }
        }
    }
    (exc, va, vb)
}

/// Solve.
pub fn solve(problem: &Problem, max_iterations: usize, tolerance: f64) -> Solution {
    let basis = &problem.basis;
    let n = basis.size;
    let (s, t, v) = one_electron(basis, &problem.nuclei);
    let mut h = t.clone();
    for k in 0..n * n {
        h.a[k] += v.a[k];
    }
    let (x, m) = orthogonaliser(&s, 1e-8);
    let eri = Repulsion::new(basis);
    let atoms: Vec<([f64; 3], f64)> = problem.nuclei.iter().zip(&problem.sizes).map(|((_, p), r)| (*p, *r)).collect();
    let grid = molecular(&atoms, problem.radial, problem.theta);
    let mut e_nn = 0.0;
    for (i, (zi, pi)) in problem.nuclei.iter().enumerate() {
        for (zj, pj) in problem.nuclei.iter().take(i) {
            let r = ((pi[0] - pj[0]).powi(2) + (pi[1] - pj[1]).powi(2) + (pi[2] - pj[2]).powi(2)).sqrt();
            e_nn += zi * zj / r;
        }
    }
    // Start from the bare-nucleus Hamiltonian.
    let (ea0, c0) = generalised(&h, &x, m);
    let mut da = density(&c0, m, n, &occupy(&ea0, problem.alpha));
    let mut db = density(&c0, m, n, &occupy(&ea0, problem.beta));
    let mut levels_a = ea0.clone();
    let mut levels_b = ea0;
    // DIIS history: Fock pair and error vector.
    let mut hist: Vec<(Matrix, Matrix, Vec<f64>)> = Vec::new();
    let mut energy = 0.0;
    let mut last = f64::INFINITY;
    let mut converged = false;
    let mut iterations = 0;
    let mut parts = (0.0, 0.0, 0.0);
    for it in 0..max_iterations {
        iterations = it + 1;
        let mut dt = da.clone();
        for k in 0..n * n {
            dt.a[k] += db.a[k];
        }
        let j = eri.coulomb(&dt);
        let (exc, vxa, vxb) = exchange_correlation(basis, &grid, problem.functional, &da, &db);
        let mut fa = h.clone();
        let mut fb = h.clone();
        for k in 0..n * n {
            fa.a[k] += j.a[k] + vxa.a[k];
            fb.a[k] += j.a[k] + vxb.a[k];
        }
        let e1 = h.dot(&dt);
        let ej = 0.5 * j.dot(&dt);
        energy = e1 + ej + exc + e_nn;
        parts = (e1, ej, exc);
        // Error FDS - SDF for each spin, together.
        let err = |f: &Matrix, d: &Matrix| {
            let fds = f.mul(d).mul(&s);
            let sdf = s.mul(d).mul(f);
            (0..n * n).map(|k| fds.a[k] - sdf.a[k]).collect::<Vec<f64>>()
        };
        let mut e = err(&fa, &da);
        e.extend(err(&fb, &db));
        let emax = e.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        if (energy - last).abs() < tolerance && emax < tolerance.sqrt() {
            converged = true;
            break;
        }
        last = energy;
        hist.push((fa.clone(), fb.clone(), e));
        if hist.len() > 10 {
            hist.remove(0);
        }
        let (fa_x, fb_x) = if hist.len() >= 2 {
            diis(&hist).unwrap_or((fa, fb))
        } else {
            (fa, fb)
        };
        let (ea, ca) = generalised(&fa_x, &x, m);
        let (eb, cb) = generalised(&fb_x, &x, m);
        da = density(&ca, m, n, &occupy(&ea, problem.alpha));
        db = density(&cb, m, n, &occupy(&eb, problem.beta));
        levels_a = ea;
        levels_b = eb;
    }
    Solution {
        energy,
        converged,
        iterations,
        levels_alpha: levels_a,
        levels_beta: levels_b,
        density_alpha: da,
        density_beta: db,
        kinetic_and_nuclear: parts.0,
        coulomb: parts.1,
        exchange_correlation: parts.2,
        nuclear_repulsion: e_nn,
    }
}

/// Pulay extrapolation of the Fock matrices.
fn diis(hist: &[(Matrix, Matrix, Vec<f64>)]) -> Option<(Matrix, Matrix)> {
    let k = hist.len();
    let dim = k + 1;
    let mut b = vec![0.0; dim * dim];
    for i in 0..k {
        for j in 0..k {
            b[i * dim + j] = hist[i].2.iter().zip(&hist[j].2).map(|(a, c)| a * c).sum();
        }
        b[i * dim + k] = -1.0;
        b[k * dim + i] = -1.0;
    }
    let mut rhs = vec![0.0; dim];
    rhs[k] = -1.0;
    let c = solve_linear(&mut b, &mut rhs, dim)?;
    let n = hist[0].0.n;
    let mut fa = Matrix::zeros(n);
    let mut fb = Matrix::zeros(n);
    for i in 0..k {
        for q in 0..n * n {
            fa.a[q] += c[i] * hist[i].0.a[q];
            fb.a[q] += c[i] * hist[i].1.a[q];
        }
    }
    Some((fa, fb))
}

/// Gaussian elimination with partial pivoting.
fn solve_linear(a: &mut [f64], b: &mut [f64], n: usize) -> Option<Vec<f64>> {
    for col in 0..n {
        let piv = (col..n).max_by(|&x, &y| a[x * n + col].abs().total_cmp(&a[y * n + col].abs()))?;
        if a[piv * n + col].abs() < 1e-300 {
            return None;
        }
        if piv != col {
            for k in 0..n {
                a.swap(piv * n + k, col * n + k);
            }
            b.swap(piv, col);
        }
        for r in (col + 1)..n {
            let f = a[r * n + col] / a[col * n + col];
            for k in col..n {
                a[r * n + k] -= f * a[col * n + k];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in (r + 1)..n {
            s -= a[r * n + k] * x[k];
        }
        x[r] = s / a[r * n + r];
    }
    Some(x)
}
