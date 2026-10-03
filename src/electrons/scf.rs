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
    /// A starting density per spin, or `None` for the bare nuclei's.
    pub guess: Option<(Matrix, Matrix)>,
    /// Non-local correlation solved self-consistently with the rest, or
    /// `None`. Its semilocal partner is `functional`.
    pub nonlocal: Option<super::vdw::NonlocalSpec>,
}

impl Problem {
    /// The same problem with every atom not in `keep` made a ghost: its basis
    /// functions (and auxiliary functions, and grid) stay, its nucleus's
    /// charge goes to zero, the electron count becomes `electrons`, and its
    /// share of the starting guess is removed — which is what a counterpoise
    /// correction solves each partner of a complex in. Keeping the ghosts'
    /// free-atom density in the guess started the field with the wrong number
    /// of electrons and cost a water monomer 24 iterations against the
    /// dimer's 11.
    pub fn with_ghosts(&self, keep: &[usize], electrons: f64) -> Problem {
        let mut q = self.clone();
        let ghost: Vec<bool> = (0..q.nuclei.len()).map(|i| !keep.contains(&i)).collect();
        for (n, g) in q.nuclei.iter_mut().zip(&ghost) {
            if *g {
                n.0 = 0.0;
            }
        }
        let spin = q.alpha - q.beta;
        q.alpha = 0.5 * (electrons + spin);
        q.beta = 0.5 * (electrons - spin);
        if let Some((da, db)) = &mut q.guess {
            let nb = q.basis.size;
            for (is, sh) in q.basis.shells.iter().enumerate() {
                let on = self.nuclei.iter().position(|(_, c)| (0..3).all(|k| (c[k] - sh.centre[k]).abs() < 1e-12));
                if on.is_some_and(|a| ghost[a]) {
                    for f in q.basis.offsets[is]..q.basis.offsets[is] + sh.size() {
                        for m in [&mut *da, &mut *db] {
                            for k in 0..nb {
                                m.a[f * nb + k] = 0.0;
                                m.a[k * nb + f] = 0.0;
                            }
                        }
                    }
                }
            }
        }
        q
    }
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
    /// The non-local correlation energy, included in `energy`; zero without
    /// one.
    pub nonlocal: f64,
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
    /// The fitted density's coefficients over the auxiliary set, for the
    /// final density; empty when Coulomb was exact. The gradient reads them
    /// rather than building the three-centre integrals a second time.
    pub fitted: Vec<f64>,
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
    /// The auxiliary functions the fit uses, and the Cholesky factor of the
    /// metric over them (see [`Fitted::new`]).
    kept: Vec<usize>,
    factor: Vec<f64>,
    /// `V` itself, which the robust energy needs (see [`Fitted::coulomb_and_energy`]).
    metric: Matrix,
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
        // The metric of a dense even-tempered set is nearly singular: of
        // water's 937 auxiliary functions, the others represent 230 to within
        // round-off. Factorised with pivoting, the factorisation takes the
        // most independent function left at each step and stops where what is
        // left is already represented, so the fit works in a well-conditioned
        // subset instead of inverting through the dependence. It replaced an
        // eigenvalue pseudo-inverse, which cost 3.1 s of a 4.8 s setup on
        // water against 0.3 s for this, and whose energy depended on where its
        // cut fell until the energy was put in the robust form
        // (`coulomb_and_energy`).
        let (kept, factor) = super::linalg::pivoted_cholesky(&v, METRIC_CUTOFF);
        let dropped = na - kept.len();
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
        let job = |indices: &mut dyn Iterator<Item = usize>| {
            let mut out: Vec<(usize, usize, Vec<f64>)> = Vec::new();
            // One worker's working space, kept across its calls: allocating
            // it per call made the threads queue on the allocator.
            let mut scratch = super::integrals::ThreeScratch::new();
            let mut block = Vec::new();
            for (a, b) in indices.map(|i| pairs[i]) {
                let (sa, sb) = (&shells[a], &shells[b]);
                let bra = super::integrals::pairs(sa, sb, 0);
                let mut vals = vec![vec![0.0; na]; sa.size() * sb.size()];
                for (ip, p) in aux.shells.iter().enumerate() {
                    super::integrals::eri_three_into(&bra, &kets[ip], sa.l, sb.l, p.l, &mut scratch, &mut block);
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
        let three: Vec<(usize, usize, Vec<f64>)> = parallel_interleaved(pairs.len(), &job).into_iter().flatten().collect();
        Fitted { n, three, kept, factor, metric: v, dropped }
    }

    /// Auxiliary directions the metric's cutoff dropped.
    pub fn dropped(&self) -> usize {
        self.dropped
    }

    /// The fit's coefficients `c = V^-1 d` for density `d`.
    pub fn coefficients(&self, d: &Matrix) -> Vec<f64> {
        let na = self.metric.n;
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
        let b: Vec<f64> = self.kept.iter().map(|&i| dp[i]).collect();
        let ck = super::linalg::cholesky_solve(&self.factor, self.kept.len(), &b);
        let mut c = vec![0.0; na];
        for (&i, x) in self.kept.iter().zip(ck) {
            c[i] = x;
        }
        c
    }

    /// The Coulomb matrix and the fitted Coulomb energy together.
    ///
    /// The energy is Dunlap's *robust* form, `b.c - 1/2 c.V.c` with
    /// `b = (mn|P) D`, rather than `1/2 b.c`: the two agree when `c` solves
    /// `V c = b` exactly, but the robust one is stationary in `c`, so an error
    /// in the coefficients enters only at second order and weighted by `V`.
    /// That matters because `V` is nearly singular — the coefficients along its
    /// smallest kept eigenvalues carry round-off amplified by `1/lambda`.
    /// Measured on a water molecule at a fixed density along a 2.4e-2 bohr
    /// path, `1/2 b.c` scattered by 2.2e-7 Ha about a smooth curve and the
    /// robust form by 1.1e-12, at the same distance from exact Coulomb; and
    /// the robust form no longer depends on where the metric is cut. The
    /// analytic gradient, `c.d(mn|P) D - 1/2 c.dV c`, was always the derivative
    /// of the robust form, and the Coulomb matrix is the same for both, since
    /// the extra term is stationary.
    pub fn coulomb_and_energy(&self, d: &Matrix) -> (Matrix, f64) {
        let c = self.coefficients(d);
        let j = self.coulomb_from(&c);
        let na = self.metric.n;
        let mut cvc = 0.0;
        for i in 0..na {
            let row = &self.metric.a[i * na..i * na + na];
            cvc += c[i] * row.iter().zip(&c).map(|(a, b)| a * b).sum::<f64>();
        }
        let e = j.dot(d) - 0.5 * cvc;
        (j, e)
    }

    /// The fitted Coulomb energy for density `d`, in the robust form.
    pub fn energy(&self, d: &Matrix) -> f64 {
        self.coulomb_and_energy(d).1
    }

    pub fn coulomb(&self, d: &Matrix) -> Matrix {
        self.coulomb_from(&self.coefficients(d))
    }

    fn coulomb_from(&self, c: &[f64]) -> Matrix {
        let n = self.n;
        let mut j = Matrix::zeros(n);
        for (m, nn, vals) in &self.three {
            let x: f64 = vals.iter().zip(c).map(|(a, b)| a * b).sum();
            j.set(*m, *nn, x);
            j.set(*nn, *m, x);
        }
        j
    }
}

/// A commutator `FDS - SDF` smaller than this is a converged density.
pub const COMMUTATOR_CONVERGED: f64 = 1e-7;

/// The pivoted factorisation of the fitting metric stops when the largest
/// diagonal left falls below this fraction of the largest there was. See
/// [`Fitted::new`]. Measured on water (937 auxiliary functions), the energy
/// at 1e-9 / 1e-11 / 1e-13 / nothing dropped: -76.387918596 / ...537 / ...535
/// / ...535 with 126 / 50 / 17 / 0 dropped — a fit in the Coulomb metric can
/// only rise towards the exact energy as functions are added, and 1e-11 is
/// 2e-9 from where it stops rising while still leaving out the fifty most
/// dependent.
pub const METRIC_CUTOFF: f64 = 1e-11;

enum Coulomb {
    Exact(Repulsion),
    Fitted(Fitted),
}

impl Coulomb {
    /// The Coulomb matrix and the Coulomb energy.
    fn build(&self, d: &Matrix) -> (Matrix, f64) {
        match self {
            Coulomb::Exact(r) => {
                let j = r.coulomb(d);
                let e = 0.5 * j.dot(d);
                (j, e)
            }
            Coulomb::Fitted(f) => f.coulomb_and_energy(d),
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
    pub batches: Vec<(Vec<[f64; 3]>, Vec<f64>, Vec<usize>, Vec<usize>)>,
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
    /// Most points in one batch.
    const MAX_POINTS: usize = 128;

    /// Points split by recursive bisection — each group halved at its median
    /// along its widest extent until it holds at most `MAX_POINTS` — and each
    /// batch given the shells that are not negligible anywhere in its box.
    ///
    /// Bisection rather than a fixed lattice of boxes: water's 72,986 points
    /// in 1.5-bohr boxes made 17,136 batches of 4.3 points on average, since
    /// a radial grid thins out away from its nucleus, and gathering and
    /// scattering a batch's 174 x 174 slice of the density for four points
    /// cost far more than the arithmetic on them.
    pub fn new(basis: &Basis, grid: &Grid) -> Batches {
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut stack: Vec<Vec<usize>> = vec![(0..grid.points.len()).collect()];
        while let Some(mut idx) = stack.pop() {
            if idx.len() <= Self::MAX_POINTS {
                if !idx.is_empty() {
                    groups.push(idx);
                }
                continue;
            }
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for &i in &idx {
                for d in 0..3 {
                    lo[d] = lo[d].min(grid.points[i][d]);
                    hi[d] = hi[d].max(grid.points[i][d]);
                }
            }
            let axis = (0..3).max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b]))).unwrap_or(0);
            idx.sort_by(|&a, &b| grid.points[a][axis].total_cmp(&grid.points[b][axis]).then(a.cmp(&b)));
            let upper = idx.split_off(idx.len() / 2);
            stack.push(upper);
            stack.push(idx);
        }
        let mut batches = Vec::with_capacity(groups.len());
        let mut indices = Vec::with_capacity(groups.len());
        for group in groups {
            let pts: Vec<[f64; 3]> = group.iter().map(|&i| grid.points[i]).collect();
            let ws: Vec<f64> = group.iter().map(|&i| grid.weights[i]).collect();
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
            indices.push(group);
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

/// One batch's contribution: exchange-correlation energy, and the two
/// potential matrices restricted to the batch's functions.
///
/// Organised as dense products over the whole batch — basis values `Phi`
/// (`points x k`), `X = Phi D`, then `V += Phi^T A + A^T Phi` — rather than
/// point by point. The arithmetic is the same; done as matrix-matrix work the
/// inner loops are contiguous and vectorise.
fn xc_batch(basis: &Basis, f: Functional, da: &Matrix, db: &Matrix, same: bool, pts: &[[f64; 3]], ws: &[f64], funcs: &[usize], shells: &[usize]) -> (f64, Vec<f64>, Vec<f64>) {
    use super::linalg::{product_nt, transpose};
    let n = basis.size;
    let k = funcs.len();
    let b = pts.len();
    let grad = f.needs_gradient();
    let spins = if same { 1 } else { 2 };
    let mut sub = [vec![0.0; k * k], vec![0.0; k * k]];
    for (s, dm) in [da, db].iter().enumerate().take(spins) {
        for (i, &fi) in funcs.iter().enumerate() {
            for (j, &fj) in funcs.iter().enumerate() {
                sub[s][i * k + j] = dm.a[fi * n + fj];
            }
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
    // X = Phi D, one per spin (the density is symmetric, so its rows serve).
    let mut x = [vec![0.0; b * k], vec![0.0; b * k]];
    for s in 0..spins {
        product_nt(&phi, &sub[s], b, k, k, &mut x[s]);
    }
    let dotk = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(a, b)| a * b).sum::<f64>();
    let mut exc = 0.0;
    let mut weights = [vec![0.0; b * k], vec![0.0; b * k]];
    for p in 0..b {
        let r = p * k..p * k + k;
        let ra = dotk(&phi[r.clone()], &x[0][r.clone()]);
        let rb = if same { ra } else { dotk(&phi[r.clone()], &x[1][r.clone()]) };
        if ra + rb < 1e-14 {
            continue;
        }
        let (mut ga, mut gb) = ([0.0; 3], [0.0; 3]);
        if grad {
            for (d, g) in [&gx, &gy, &gz].iter().enumerate() {
                ga[d] = 2.0 * dotk(&g[r.clone()], &x[0][r.clone()]);
                gb[d] = if same { ga[d] } else { 2.0 * dotk(&g[r.clone()], &x[1][r.clone()]) };
            }
        }
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        let (e, de) = evaluate(f, ra, rb, dot(ga, ga), dot(ga, gb), dot(gb, gb));
        let w = ws[p];
        exc += w * e;
        let wa = [2.0 * de[2] * ga[0] + de[3] * gb[0], 2.0 * de[2] * ga[1] + de[3] * gb[1], 2.0 * de[2] * ga[2] + de[3] * gb[2]];
        let wb = [2.0 * de[4] * gb[0] + de[3] * ga[0], 2.0 * de[4] * gb[1] + de[3] * ga[1], 2.0 * de[4] * gb[2] + de[3] * ga[2]];
        for (s, (vr, wv)) in [(de[0], wa), (de[1], wb)].iter().enumerate().take(spins) {
            for i in 0..k {
                let q = p * k + i;
                let mut a = 0.5 * vr * phi[q];
                if grad {
                    a += wv[0] * gx[q] + wv[1] * gy[q] + wv[2] * gz[q];
                }
                weights[s][q] = w * a;
            }
        }
    }
    // V = Phi^T A + A^T Phi, as M + M^T with M = Phi^T A.
    let mut phit = vec![0.0; k * b];
    transpose(&phi, b, k, &mut phit);
    let mut at = vec![0.0; k * b];
    let mut m = vec![0.0; k * k];
    let mut v = [vec![0.0; k * k], Vec::new()];
    for s in 0..spins {
        transpose(&weights[s], b, k, &mut at);
        product_nt(&phit, &at, k, k, b, &mut m);
        let mut out = vec![0.0; k * k];
        for i in 0..k {
            for j in 0..k {
                out[i * k + j] = m[i * k + j] + m[j * k + i];
            }
        }
        v[s] = out;
    }
    let [va, vb] = v;
    let vb = if same { va.clone() } else { vb };
    (exc, va, vb)
}


/// Exchange-correlation energy and the two spin potentials' matrices.
///
/// Batched by place and screened: a point only evaluates the functions that
/// are not negligible over its batch's box (oxygen's tightest `s` function is
/// zero a hundredth of a bohr from the nucleus and used to be computed at every
/// one of a hundred thousand points). Batches are spread over threads where
/// there are threads; on wasm32 they run in turn. When the two spins' densities
/// are the same to the last bit — a closed shell, which the solve keeps exactly
/// so — one spin's work is done once and serves both.
pub fn exchange_correlation(basis: &Basis, batches: &Batches, f: Functional, da: &Matrix, db: &Matrix) -> (f64, Matrix, Matrix) {
    let n = basis.size;
    let work = &batches.batches;
    let same = da.a == db.a;
    let run = |indices: &mut dyn Iterator<Item = usize>| {
        let mut va = Matrix::zeros(n);
        let mut vb = Matrix::zeros(n);
        let mut exc = 0.0;
        for bi in indices {
            let (pts, ws, funcs, shells) = &work[bi];
            if funcs.is_empty() {
                continue;
            }
            let (e, a, b) = xc_batch(basis, f, da, db, same, pts, ws, funcs, shells);
            exc += e;
            let k = funcs.len();
            for (i, &fi) in funcs.iter().enumerate() {
                let (row_a, row_b) = (&a[i * k..i * k + k], &b[i * k..i * k + k]);
                for (j, &fj) in funcs.iter().enumerate() {
                    va.a[fi * n + fj] += row_a[j];
                    vb.a[fi * n + fj] += row_b[j];
                }
            }
        }
        (exc, va, vb)
    };
    let parts = parallel_interleaved(work.len(), &run);
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

/// Run `job` on one share of `0..len` per worker, returning the results in
/// worker order. Each worker takes every `workers`-th item rather than a
/// contiguous block. Batches of grid points are sorted by place, so a block
/// holds the dense, expensive ones near a nucleus or the sparse cheap ones far
/// out, and contiguous blocks left three workers idle while the fourth
/// finished: 35 thread-seconds of work took 20 s on four threads. Which
/// worker takes which item depends only on the worker count, so the result is
/// the same however the threads are scheduled.
pub fn parallel_interleaved<T: Send>(len: usize, job: &(dyn Fn(&mut dyn Iterator<Item = usize>) -> T + Sync)) -> Vec<T> {
    #[cfg(not(target_arch = "wasm32"))]
    let workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).min(len.max(1));
    #[cfg(target_arch = "wasm32")]
    let workers = 1;
    if workers <= 1 {
        return vec![job(&mut (0..len))];
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers).map(|w| scope.spawn(move || job(&mut (w..len).step_by(workers)))).collect();
            handles.into_iter().map(|h| h.join().expect("a worker panicked")).collect()
        })
    }
    #[cfg(target_arch = "wasm32")]
    {
        vec![job(&mut (0..len))]
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
    // The non-local term's own grid, built once like the semilocal one.
    let nonlocal = problem.nonlocal.map(|spec| {
        let g = super::grid::molecular_pruned(&atoms, spec.radial, spec.theta, false);
        let b = Batches::new(basis, &g);
        (spec, g, b)
    });
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
    let (mut da, mut db) = match &problem.guess {
        Some((a, b)) => (a.clone(), b.clone()),
        None => (density(&c0, m, n, &occupy(&ea0, problem.alpha)), density(&c0, m, n, &occupy(&ea0, problem.beta))),
    };
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
    let mut last_nonlocal = 0.0;
    for it in 0..max_iterations {
        iterations = it + 1;
        let mut dt = da.clone();
        for k in 0..n * n {
            dt.a[k] += db.a[k];
        }
        let clock = std::time::Instant::now();
        let (j, ej) = eri.build(&dt);
        timing[2] += clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        let (exc, mut vxa, mut vxb) = exchange_correlation(basis, &batches, problem.functional, &da, &db);
        let mut enl = 0.0;
        if let Some((spec, g, b)) = &nonlocal {
            let dens = super::vdw::density_and_gradient_on(basis, b, &dt, g.points.len());
            let nl = super::vdw::nonlocal(g, &dens, spec.z_ab, spec.floor, super::vdw::kernel_table());
            let v = super::vdw::nonlocal_matrix(basis, b, &dens, &nl);
            for k in 0..n * n {
                vxa.a[k] += v.a[k];
                vxb.a[k] += v.a[k];
            }
            enl = nl.energy;
        }
        timing[3] += clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        let mut fa = h.clone();
        let mut fb = h.clone();
        for k in 0..n * n {
            fa.a[k] += j.a[k] + vxa.a[k];
            fb.a[k] += j.a[k] + vxb.a[k];
        }
        let e1 = h.dot(&dt);
        energy = e1 + ej + exc + enl + e_nn;
        last_nonlocal = enl;
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
    let fitted = match &eri {
        Coulomb::Fitted(f) => {
            let mut dt = da.clone();
            for k in 0..n * n {
                dt.a[k] += db.a[k];
            }
            f.coefficients(&dt)
        }
        Coulomb::Exact(_) => Vec::new(),
    };
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
        nonlocal: last_nonlocal,
        timing,
        history,
        kept: m,
        smallest_overlap,
        fitted,
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
