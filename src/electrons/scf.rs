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
use super::grid::Grid;
use super::integrals::{eri_block, one_electron};
use super::linalg::{generalised, orthogonaliser, Matrix};

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
    /// An auxiliary basis for fitting the Coulomb potential, or `None` for the
    /// exact four-index repulsion.
    pub auxiliary: Option<Basis>,
    /// Prune the grid's angular order near nuclei and far out.
    pub prune: bool,
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
    /// Energy-weighted densities `sum_i n_i e_i c_i c_i^T`, per spin: what the
    /// gradient's overlap (Pulay) term needs, because the basis moves with the
    /// atoms.
    pub weighted_alpha: Matrix,
    pub weighted_beta: Matrix,
    /// The parts of the energy, hartree.
    pub kinetic_and_nuclear: f64,
    pub coulomb: f64,
    pub exchange_correlation: f64,
    pub nuclear_repulsion: f64,
    /// Seconds spent: one-electron and repulsion integrals, building the grid,
    /// Coulomb builds, exchange-correlation builds, diagonalisation.
    pub timing: [f64; 5],
    /// Each iteration's energy and largest commutator error, for diagnosing
    /// convergence.
    pub history: Vec<(f64, f64)>,
    /// Basis functions kept by the orthogonaliser, and the overlap's smallest
    /// eigenvalue.
    pub kept: usize,
    pub smallest_overlap: f64,
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

/// The Coulomb potential by density fitting.
///
/// The electron density is fitted, in the Coulomb metric, by auxiliary
/// functions `P`: `d_P = sum (P|ls) D_ls`, `V c = d` with `V_PQ = (P|Q)`, and
/// `J_mn = sum_P (mn|P) c_P`. Three-centre integrals in place of four: for a
/// contracted basis this is what makes the cost per pair linear in the
/// auxiliary functions rather than in the square of the primitives. In the
/// Coulomb metric the fitted energy is a strict lower bound of the exact one,
/// and the gap is quadratic in the fitting error.
pub struct Fitted {
    n: usize,
    /// `(mn|P)` for each kept pair `m >= n`: `(m, n, values over P)`.
    three: Vec<(usize, usize, Vec<f64>)>,
    /// Pseudo-inverse of `V` by eigenvalues, dropping the near-dependent ones.
    vinv: Matrix,
    dropped: usize,
}

impl Fitted {
    pub fn new(basis: &Basis, aux: &Basis, s: &Matrix) -> Fitted {
        use super::basis::Shell;
        let unit = Shell::unit();
        let na = aux.size;
        let mut v = Matrix::zeros(na);
        for (ip, p) in aux.shells.iter().enumerate() {
            for (iq, q) in aux.shells.iter().enumerate().take(ip + 1) {
                let block = eri_block(p, &unit, q, &unit);
                for a in 0..p.size() {
                    for b in 0..q.size() {
                        let val = block[a * q.size() + b];
                        v.set(aux.offsets[ip] + a, aux.offsets[iq] + b, val);
                        v.set(aux.offsets[iq] + b, aux.offsets[ip] + a, val);
                    }
                }
            }
        }
        // The metric of a dense even-tempered set is nearly singular, and
        // inverting its smallest eigenvalues amplifies round-off past the fit
        // itself. Measured on water's Coulomb matrix, 1575 auxiliary
        // functions: dropping below 1e-9 of the largest, relative error 1.8e-8;
        // below 1e-14, 6e-4 and a "fitted" energy *above* the exact one, which
        // a Coulomb-metric fit cannot honestly give; below 1e-16, 6%.
        let (vals, vecs) = super::linalg::eigh(&v);
        let mut vinv = Matrix::zeros(na);
        let top = vals.last().cloned().unwrap_or(1.0);
        let mut dropped = 0;
        for k in 0..na {
            if vals[k] <= METRIC_CUTOFF * top {
                dropped += 1;
                continue;
            }
            for i in 0..na {
                let f = vecs.get(i, k) / vals[k];
                for j in 0..na {
                    vinv.a[i * na + j] += f * vecs.get(j, k);
                }
            }
        }
        // Three-centre integrals over shell pairs whose product is not
        // negligible anywhere. Not by overlap: an s and a p function on the
        // same atom overlap by exactly zero, and their product is a dipole
        // whose field is not — screening on the overlap dropped those and left
        // the Coulomb matrix 4% wrong whatever the auxiliary set. The product
        // of two Gaussians carries `exp(-a b / (a + b) R^2)`, which is what
        // decides it.
        let n = basis.size;
        let _ = s;
        let shells = &basis.shells;
        let pairs: Vec<(usize, usize)> = (0..shells.len()).flat_map(|a| (0..=a).map(move |b| (a, b))).filter(|&(a, b)| {
            let (sa, sb) = (&shells[a], &shells[b]);
            let r2: f64 = (0..3).map(|d| (sa.centre[d] - sb.centre[d]).powi(2)).sum();
            sa.exponents.iter().any(|x| sb.exponents.iter().any(|y| (-(x * y) / (x + y) * r2).exp() > 1e-14))
        }).collect();
        // Each auxiliary shell's pair with the unit function, built once.
        let kets: Vec<Vec<super::integrals::Pair>> = aux.shells.iter().map(|p| super::integrals::pairs(p, &unit, 0)).collect();
        let job = |range: std::ops::Range<usize>| {
            let mut out: Vec<(usize, usize, Vec<f64>)> = Vec::new();
            for &(a, b) in &pairs[range] {
                let (sa, sb) = (&shells[a], &shells[b]);
                let bra = super::integrals::pairs(sa, sb, 0);
                let mut vals = vec![vec![0.0; na]; sa.size() * sb.size()];
                for (ip, p) in aux.shells.iter().enumerate() {
                    let block = super::integrals::eri_from_pairs(&bra, &kets[ip], [sa.l, sb.l, p.l, 0]);
                    for i in 0..sa.size() {
                        for j in 0..sb.size() {
                            for k in 0..p.size() {
                                vals[i * sb.size() + j][aux.offsets[ip] + k] = block[(i * sb.size() + j) * p.size() + k];
                            }
                        }
                    }
                }
                for i in 0..sa.size() {
                    for j in 0..sb.size() {
                        let (m, nn) = (basis.offsets[a] + i, basis.offsets[b] + j);
                        if m >= nn {
                            out.push((m, nn, std::mem::take(&mut vals[i * sb.size() + j])));
                        }
                    }
                }
            }
            out
        };
        let three: Vec<(usize, usize, Vec<f64>)> = parallel(pairs.len(), &job).into_iter().flatten().collect();
        Fitted { n, three, vinv, dropped }
    }

    /// Auxiliary directions the metric's cutoff dropped.
    pub fn dropped(&self) -> usize {
        self.dropped
    }

    /// The fit's coefficients `c = V^-1 d` for density `d`.
    pub fn coefficients(&self, d: &Matrix) -> Vec<f64> {
        let na = self.vinv.n;
        let mut dp = vec![0.0; na];
        for (m, nn, vals) in &self.three {
            let w = if m == nn { d.get(*m, *m) } else { d.get(*m, *nn) + d.get(*nn, *m) };
            if w == 0.0 {
                continue;
            }
            for (x, v) in dp.iter_mut().zip(vals) {
                *x += w * v;
            }
        }
        (0..na).map(|i| self.vinv.a[i * na..i * na + na].iter().zip(&dp).map(|(a, b)| a * b).sum()).collect()
    }

    pub fn coulomb(&self, d: &Matrix) -> Matrix {
        let n = self.n;
        let na = self.vinv.n;
        let mut dp = vec![0.0; na];
        for (m, nn, vals) in &self.three {
            let w = if m == nn { d.get(*m, *m) } else { d.get(*m, *nn) + d.get(*nn, *m) };
            if w == 0.0 {
                continue;
            }
            for (x, v) in dp.iter_mut().zip(vals) {
                *x += w * v;
            }
        }
        let mut c = vec![0.0; na];
        for i in 0..na {
            let row = &self.vinv.a[i * na..i * na + na];
            c[i] = row.iter().zip(&dp).map(|(a, b)| a * b).sum();
        }
        let mut j = Matrix::zeros(n);
        for (m, nn, vals) in &self.three {
            let x: f64 = vals.iter().zip(&c).map(|(a, b)| a * b).sum();
            j.set(*m, *nn, x);
            j.set(*nn, *m, x);
        }
        j
    }
}

/// A commutator `FDS - SDF` smaller than this is a converged density.
pub const COMMUTATOR_CONVERGED: f64 = 1e-7;

/// Eigenvalues of the fitting metric below this fraction of the largest are
/// dropped rather than inverted. See [`Fitted::new`].
pub const METRIC_CUTOFF: f64 = 1e-9;

enum Coulomb {
    Exact(Repulsion),
    Fitted(Fitted),
}

impl Coulomb {
    fn build(&self, d: &Matrix) -> Matrix {
        match self {
            Coulomb::Exact(r) => r.coulomb(d),
            Coulomb::Fitted(f) => f.coulomb(d),
        }
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

fn weighted(c: &[f64], m: usize, n: usize, occ: &[f64], levels: &[f64]) -> Matrix {
    let w: Vec<f64> = occ.iter().zip(levels).map(|(o, e)| o * e).collect();
    density(c, m, n, &w)
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

/// Points grouped by where they are, and the basis functions that matter there.
pub struct Batches {
    /// `(points, weights, functions that are not negligible over them, and
    /// the shells those functions belong to)`.
    pub(crate) batches: Vec<(Vec<[f64; 3]>, Vec<f64>, Vec<usize>, Vec<usize>)>,
    /// For each batch, its points' indices in the grid it was built from.
    pub(crate) indices: Vec<Vec<usize>>,
}

/// Below this a function's (or its gradient's) magnitude is taken as zero.
const NEGLIGIBLE: f64 = 1e-12;

/// The largest `|r^l exp(-a r^2)|`, times `(1 + 2 a r)` for the gradient, can be
/// anywhere at distance `>= d` from its centre.
fn shell_bound(l: usize, a: f64, coef: f64, d: f64) -> f64 {
    let peak = if l == 0 { 0.0 } else { (l as f64 / (2.0 * a)).sqrt() };
    let r = d.max(peak);
    coef.abs() * r.powi(l as i32) * (-a * r * r).exp() * (1.0 + 2.0 * a * r + l as f64 / r.max(1e-12))
}

impl Batches {
    /// Box size for grouping points, bohr.
    const BOX: f64 = 1.5;
    const MAX_POINTS: usize = 256;

    pub fn new(basis: &Basis, grid: &Grid) -> Batches {
        let mut order: Vec<usize> = (0..grid.points.len()).collect();
        let cell = |p: [f64; 3]| [(p[0] / Self::BOX).floor() as i64, (p[1] / Self::BOX).floor() as i64, (p[2] / Self::BOX).floor() as i64];
        order.sort_by_key(|&i| cell(grid.points[i]));
        let mut batches = Vec::new();
        let mut indices = Vec::new();
        let mut start = 0;
        while start < order.len() {
            let c0 = cell(grid.points[order[start]]);
            let mut end = start + 1;
            while end < order.len() && end - start < Self::MAX_POINTS && cell(grid.points[order[end]]) == c0 {
                end += 1;
            }
            let pts: Vec<[f64; 3]> = order[start..end].iter().map(|&i| grid.points[i]).collect();
            let ws: Vec<f64> = order[start..end].iter().map(|&i| grid.weights[i]).collect();
            // Bounding box of the batch.
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for p in &pts {
                for d in 0..3 {
                    lo[d] = lo[d].min(p[d]);
                    hi[d] = hi[d].max(p[d]);
                }
            }
            let mut funcs = Vec::new();
            let mut shell_list = Vec::new();
            for (is, sh) in basis.shells.iter().enumerate() {
                let mut d2 = 0.0;
                for d in 0..3 {
                    let x = sh.centre[d];
                    let gap = if x < lo[d] { lo[d] - x } else if x > hi[d] { x - hi[d] } else { 0.0 };
                    d2 += gap * gap;
                }
                let dist = d2.sqrt();
                let bound: f64 = sh.exponents.iter().zip(&sh.coefficients).map(|(a, c)| shell_bound(sh.l, *a, *c, dist)).sum();
                if bound > NEGLIGIBLE {
                    shell_list.push(is);
                    for k in 0..sh.size() {
                        funcs.push(basis.offsets[is] + k);
                    }
                }
            }
            batches.push((pts, ws, funcs, shell_list));
            indices.push(order[start..end].to_vec());
            start = end;
        }
        Batches { batches, indices }
    }

    /// Mean number of functions a point sees, against the basis size.
    pub fn mean_functions(&self) -> f64 {
        let (mut n, mut w) = (0.0, 0.0);
        for (p, _, f, _) in &self.batches {
            n += (p.len() * f.len()) as f64;
            w += p.len() as f64;
        }
        n / w.max(1.0)
    }
}

/// `out[p][j] = sum_i a[p][i] m[i][j]`: `a` is `rows x k`, `m` is `k x k`.
fn rows_times(a: &[f64], m: &[f64], rows: usize, k: usize, out: &mut [f64]) {
    for o in out.iter_mut() {
        *o = 0.0;
    }
    for p in 0..rows {
        let ar = &a[p * k..p * k + k];
        let or = &mut out[p * k..p * k + k];
        for (i, &x) in ar.iter().enumerate() {
            if x == 0.0 {
                continue;
            }
            let mr = &m[i * k..i * k + k];
            for (o, &y) in or.iter_mut().zip(mr) {
                *o += x * y;
            }
        }
    }
}

/// `out[i][j] += sum_p a[p][i] b[p][j] + b[p][i] a[p][j]`, symmetric.
fn add_symmetric_products(a: &[f64], b: &[f64], rows: usize, k: usize, out: &mut [f64]) {
    for p in 0..rows {
        let ar = &a[p * k..p * k + k];
        let br = &b[p * k..p * k + k];
        for i in 0..k {
            let (ai, bi) = (ar[i], br[i]);
            let orow = &mut out[i * k..i * k + k];
            for j in 0..k {
                orow[j] += ai * br[j] + bi * ar[j];
            }
        }
    }
}

/// One batch's contribution: exchange-correlation energy, and the two
/// potential matrices restricted to the batch's functions.
///
/// Organised as dense products over the whole batch — basis values `Phi`
/// (`points x k`), `X = Phi D`, then `V += Phi^T A + A^T Phi` — rather than
/// point by point. The arithmetic is the same; done as matrix-matrix work the
/// inner loops are contiguous and vectorise.
fn xc_batch(basis: &Basis, f: Functional, da: &Matrix, db: &Matrix, pts: &[[f64; 3]], ws: &[f64], funcs: &[usize], shells: &[usize]) -> (f64, Vec<f64>, Vec<f64>) {
    let n = basis.size;
    let k = funcs.len();
    let b = pts.len();
    let grad = f.needs_gradient();
    let mut sa = vec![0.0; k * k];
    let mut sb = vec![0.0; k * k];
    for (i, &fi) in funcs.iter().enumerate() {
        for (j, &fj) in funcs.iter().enumerate() {
            sa[i * k + j] = da.a[fi * n + fj];
            sb[i * k + j] = db.a[fi * n + fj];
        }
    }
    // Basis values and gradients over the batch.
    let mut phi = vec![0.0; b * k];
    let mut gx = vec![0.0; b * k];
    let mut gy = vec![0.0; b * k];
    let mut gz = vec![0.0; b * k];
    let mut g3 = vec![0.0; 3 * k];
    for (p, pt) in pts.iter().enumerate() {
        super::values::at_shells(basis, shells, *pt, &mut phi[p * k..p * k + k], if grad { Some(&mut g3) } else { None });
        if grad {
            for i in 0..k {
                gx[p * k + i] = g3[3 * i];
                gy[p * k + i] = g3[3 * i + 1];
                gz[p * k + i] = g3[3 * i + 2];
            }
        }
    }
    let mut xa = vec![0.0; b * k];
    let mut xb = vec![0.0; b * k];
    rows_times(&phi, &sa, b, k, &mut xa);
    rows_times(&phi, &sb, b, k, &mut xb);
    let dotk = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(a, b)| a * b).sum::<f64>();
    let mut exc = 0.0;
    let mut aa = vec![0.0; b * k];
    let mut ab = vec![0.0; b * k];
    for p in 0..b {
        let r = p * k..p * k + k;
        let ra = dotk(&phi[r.clone()], &xa[r.clone()]);
        let rb = dotk(&phi[r.clone()], &xb[r.clone()]);
        if ra + rb < 1e-14 {
            continue;
        }
        let (mut ga, mut gb) = ([0.0; 3], [0.0; 3]);
        if grad {
            for (d, g) in [&gx, &gy, &gz].iter().enumerate() {
                ga[d] = 2.0 * dotk(&g[r.clone()], &xa[r.clone()]);
                gb[d] = 2.0 * dotk(&g[r.clone()], &xb[r.clone()]);
            }
        }
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        let (e, de) = evaluate(f, ra, rb, dot(ga, ga), dot(ga, gb), dot(gb, gb));
        let w = ws[p];
        exc += w * e;
        let wa = [2.0 * de[2] * ga[0] + de[3] * gb[0], 2.0 * de[2] * ga[1] + de[3] * gb[1], 2.0 * de[2] * ga[2] + de[3] * gb[2]];
        let wb = [2.0 * de[4] * gb[0] + de[3] * ga[0], 2.0 * de[4] * gb[1] + de[3] * ga[1], 2.0 * de[4] * gb[2] + de[3] * ga[2]];
        for i in 0..k {
            let q = p * k + i;
            let mut a = 0.5 * de[0] * phi[q];
            let mut c = 0.5 * de[1] * phi[q];
            if grad {
                a += wa[0] * gx[q] + wa[1] * gy[q] + wa[2] * gz[q];
                c += wb[0] * gx[q] + wb[1] * gy[q] + wb[2] * gz[q];
            }
            aa[q] = w * a;
            ab[q] = w * c;
        }
    }
    let mut va = vec![0.0; k * k];
    let mut vb = vec![0.0; k * k];
    add_symmetric_products(&phi, &aa, b, k, &mut va);
    add_symmetric_products(&phi, &ab, b, k, &mut vb);
    (exc, va, vb)
}

/// Exchange-correlation energy and the two spin potentials' matrices.
///
/// Batched by place and screened: a point only evaluates the functions that
/// are not negligible over its batch's box (oxygen's tightest `s` function is
/// zero a hundredth of a bohr from the nucleus and used to be computed at every
/// one of a hundred thousand points). Batches are spread over threads where
/// there are threads; on wasm32 they run in turn.
pub fn exchange_correlation(basis: &Basis, batches: &Batches, f: Functional, da: &Matrix, db: &Matrix) -> (f64, Matrix, Matrix) {
    let n = basis.size;
    let work = &batches.batches;
    let run = |range: std::ops::Range<usize>| {
        let mut va = Matrix::zeros(n);
        let mut vb = Matrix::zeros(n);
        let mut exc = 0.0;
        for (pts, ws, funcs, shells) in &work[range] {
            if funcs.is_empty() {
                continue;
            }
            let (e, a, b) = xc_batch(basis, f, da, db, pts, ws, funcs, shells);
            exc += e;
            let k = funcs.len();
            for (i, &fi) in funcs.iter().enumerate() {
                for (j, &fj) in funcs.iter().enumerate() {
                    va.a[fi * n + fj] += a[i * k + j];
                    vb.a[fi * n + fj] += b[i * k + j];
                }
            }
        }
        (exc, va, vb)
    };
    let parts = parallel(work.len(), &run);
    let mut va = Matrix::zeros(n);
    let mut vb = Matrix::zeros(n);
    let mut exc = 0.0;
    // Summed in a fixed order, so the result does not depend on scheduling.
    for (e, a, b) in parts {
        exc += e;
        for q in 0..n * n {
            va.a[q] += a.a[q];
            vb.a[q] += b.a[q];
        }
    }
    (exc, va, vb)
}

/// Split `0..len` into one contiguous range per worker and run `job` on each,
/// returning the results in range order: deterministic whatever the threads do.
pub fn parallel<T: Send>(len: usize, job: &(dyn Fn(std::ops::Range<usize>) -> T + Sync)) -> Vec<T> {
    #[cfg(not(target_arch = "wasm32"))]
    let workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).min(len.max(1));
    #[cfg(target_arch = "wasm32")]
    let workers = 1;
    let chunk = len.div_ceil(workers.max(1)).max(1);
    let ranges: Vec<std::ops::Range<usize>> = (0..workers).map(|w| (w * chunk).min(len)..((w + 1) * chunk).min(len)).collect();
    if workers <= 1 {
        return ranges.into_iter().map(job).collect();
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::thread::scope(|scope| {
            let handles: Vec<_> = ranges.into_iter().map(|r| scope.spawn(move || job(r))).collect();
            handles.into_iter().map(|h| h.join().expect("a worker panicked")).collect()
        })
    }
    #[cfg(target_arch = "wasm32")]
    {
        ranges.into_iter().map(job).collect()
    }
}

/// Solve.
pub fn solve(problem: &Problem, max_iterations: usize, tolerance: f64) -> Solution {
    let basis = &problem.basis;
    let n = basis.size;
    let mut timing = [0.0f64; 5];
    let clock = std::time::Instant::now();
    let (s, t, v) = one_electron(basis, &problem.nuclei);
    let mut h = t.clone();
    for k in 0..n * n {
        h.a[k] += v.a[k];
    }
    let (x, m) = orthogonaliser(&s, 1e-8);
    let smallest_overlap = super::linalg::eigh(&s).0.first().cloned().unwrap_or(0.0);
    let mut history = Vec::new();
    let eri = match &problem.auxiliary {
        None => Coulomb::Exact(Repulsion::new(basis)),
        Some(aux) => Coulomb::Fitted(Fitted::new(basis, aux, &s)),
    };
    timing[0] = clock.elapsed().as_secs_f64();
    let clock = std::time::Instant::now();
    let atoms: Vec<([f64; 3], f64)> = problem.nuclei.iter().zip(&problem.sizes).map(|((_, p), r)| (*p, *r)).collect();
    let grid = super::grid::molecular_pruned(&atoms, problem.radial, problem.theta, problem.prune);
    let batches = Batches::new(basis, &grid);
    timing[1] = clock.elapsed().as_secs_f64();
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
    let mut wa = weighted(&c0, m, n, &occupy(&ea0, problem.alpha), &ea0);
    let mut wb = weighted(&c0, m, n, &occupy(&ea0, problem.beta), &ea0);
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
        let clock = std::time::Instant::now();
        let j = eri.build(&dt);
        timing[2] += clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        let (exc, vxa, vxb) = exchange_correlation(basis, &batches, problem.functional, &da, &db);
        timing[3] += clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
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
        history.push((energy, emax));
        // Converged when the energy has stopped moving with the commutator
        // small, or when the commutator alone is negligible: the energy has a
        // numerical floor near 1e-7 hartree on a finite grid, and water once
        // spent 68 iterations after converging chasing it below 1e-9.
        if ((energy - last).abs() < tolerance && emax < tolerance.sqrt()) || emax < COMMUTATOR_CONVERGED {
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
        let (oa, ob) = (occupy(&ea, problem.alpha), occupy(&eb, problem.beta));
        da = density(&ca, m, n, &oa);
        db = density(&cb, m, n, &ob);
        wa = weighted(&ca, m, n, &oa, &ea);
        wb = weighted(&cb, m, n, &ob, &eb);
        timing[4] += clock.elapsed().as_secs_f64();
        levels_a = ea;
        levels_b = eb;
    }
    Solution {
        energy,
        converged,
        iterations,
        levels_alpha: levels_a,
        levels_beta: levels_b,
        weighted_alpha: wa,
        weighted_beta: wb,
        density_alpha: da,
        density_beta: db,
        kinetic_and_nuclear: parts.0,
        coulomb: parts.1,
        exchange_correlation: parts.2,
        nuclear_repulsion: e_nn,
        timing,
        history,
        kept: m,
        smallest_overlap,
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
