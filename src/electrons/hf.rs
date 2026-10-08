//! Hartree-Fock and second-order Moller-Plesset perturbation theory, by density
//! fitting — `PLAY.md` Phase 6, E9: what a molecule's electrons do to a
//! neighbour's, without a functional.
//!
//! **Why this is here.** Every density functional is an approximation to the
//! exchange-correlation energy that the Schrodinger equation does not need, and
//! the engine's first liquid showed what choosing one costs: PBE exchange with
//! the vdW-DF kernel overbinds water's dimer by 12% and its liquid by 4 kcal/mol
//! a molecule, revPBE underbinds the dimer, and neither builds a tetrahedral
//! network (`PLAY.md` E8). The literature has no functional that is right for
//! every molecule either. Perturbation theory has no such choice to make: its
//! inputs are the Hamiltonian, a basis and an auxiliary basis, all of which
//! the engine already derives. Its energy is not exact (MP2 overbinds
//! dispersion between large polarisable molecules, and a basis is finite), but
//! its faults are known and systematic, and unlike a functional's they come
//! with a path to better: the next term in the series.
//!
//! **What it is.** Closed-shell restricted Hartree-Fock with the Coulomb and
//! exchange matrices built from the same density-fitted three-centre table the
//! Kohn-Sham solver uses ([`Fitted`]), so the two share a table and a fit; then
//! the RI-MP2 correlation energy from the same table,
//! `E2 = sum_ijab (ia|jb) [2 (ia|jb) - (ib|ja)] / (e_i + e_j - e_a - e_b)`
//! with `(ia|jb)` the fitted integral. The core orbitals may be frozen (the
//! caller says how many), as is usual for a basis without core-valence
//! functions.
//!
//! **What is not here yet.** Open-shell references, forces, and anything past
//! second order. The memory is `count * n * nk` doubles for the half-transformed
//! integrals and the same again whitened, and `active * virtuals * nk` for the
//! transformed ones, which for six waters is a few gigabytes.

use super::integrals::one_electron;
use super::linalg::{generalised, orthogonaliser, product_nt, Matrix};
use super::scf::{diis, Fitted, Problem};
use std::sync::{Arc, Mutex};

/// A closed-shell Hartree-Fock solution.
pub struct HartreeFock {
    /// Total energy, hartree: one-electron, fitted Coulomb, exchange from the
    /// same fit, and the nuclei.
    pub energy: f64,
    pub converged: bool,
    pub iterations: usize,
    /// Orbitals as `c[ao * m + orbital]`, `m` of them, ascending in energy.
    pub c: Vec<f64>,
    pub m: usize,
    pub levels: Vec<f64>,
    /// Doubly occupied orbitals.
    pub occupied: usize,
}

/// The exchange matrix `K[P]_mn = sum_ls P_ls (ml|ns)` for `P = sum_i c_i c_i^T`
/// over the orbitals the whitened blocks `y` were made from.
fn exchange(y: &[f64], n: usize, count: usize, nk: usize) -> Matrix {
    let mut k = Matrix::zeros(n);
    let mut tmp = vec![0.0f64; n * n];
    for i in 0..count {
        let yi = &y[i * n * nk..(i + 1) * n * nk];
        product_cpu(yi, yi, n, n, nk, &mut tmp);
        for (a, b) in k.a.iter_mut().zip(&tmp) {
            *a += b;
        }
    }
    k
}

/// A dense-product engine the correlated methods may use in place of the CPU:
/// the work that grows fastest in these methods (the exchange matrix, the
/// transformation of the fitted integrals to virtual orbitals, the pair
/// integrals of MP2) is all one shape, `out = a b^T`, and a GPU does it in
/// single precision at a small fraction of the CPU's time. The core crate
/// carries no GPU code; `phys-gpu` implements this and a program that wants it
/// installs it, as it does the non-local correlation's rows.
///
/// **Single precision is not the default.** Without an engine every product is
/// the CPU's, in double. An engine's products carry single precision's relative
/// error (a few parts in a million over thousands of terms), which is harmless
/// where it enters as a sum of many small signed terms (the MP2 energy) and
/// worth checking where it does not (the exchange matrix, whose error in the
/// energy is second order for the converged orbitals); the precision is
/// measured in `PLAY.md` E8a before anything relies on it.
pub trait ProductEngine: Send + Sync {
    /// `out = a b^T`: `a` is `m x kd`, `b` is `n x kd`, rows contiguous, `out`
    /// is `m x n` and overwritten.
    fn product_nt(&self, a: &[f64], b: &[f64], m: usize, n: usize, kd: usize, out: &mut [f64]);
    fn name(&self) -> &str;
}

static PRODUCT_ENGINE: Mutex<Option<Arc<dyn ProductEngine>>> = Mutex::new(None);

/// Install an engine for this process (`None` removes it).
pub fn set_product_engine(engine: Option<Arc<dyn ProductEngine>>) {
    *PRODUCT_ENGINE.lock().expect("the product engine") = engine;
}

/// Below this many multiply-adds a product is left to the CPU: sending the
/// operands and fetching the result costs more than it saves.
const ENGINE_MIN_WORK: f64 = 2e8;

/// `out = a b^T` (`a` is `m x kd`, `b` is `n x kd`, rows contiguous): on the
/// installed engine if there is one and the product is large enough, otherwise
/// on the CPU across threads, which gives exactly what [`product_nt`] gives.
///
/// Used for the MP2 products, where single precision is harmless. **Not used
/// for the exchange matrix of the self-consistent field**, which stays in double
/// whatever is installed: measured, with it in single precision a water dimer's
/// Hartree-Fock did not converge, because the noise in the Fock matrix (a part
/// in ten million of its exchange) keeps the commutator and the energy change
/// from falling to the tolerance, and a converged energy is what the rest
/// depends on. The exchange is also the smaller part of an iteration at six
/// waters (about a sixth), so little is lost.
pub fn product_parallel(a: &[f64], b: &[f64], m: usize, n: usize, kd: usize, out: &mut [f64]) {
    if (m as f64) * (n as f64) * (kd as f64) >= ENGINE_MIN_WORK {
        let engine = PRODUCT_ENGINE.lock().expect("the product engine").clone();
        if let Some(e) = engine {
            e.product_nt(a, b, m, n, kd, out);
            return;
        }
    }
    product_cpu(a, b, m, n, kd, out);
}

/// [`product_parallel`] on the CPU in double precision, whatever is installed.
pub fn product_cpu(a: &[f64], b: &[f64], m: usize, n: usize, kd: usize, out: &mut [f64]) {
    let threads = std::thread::available_parallelism().map(|t| t.get()).unwrap_or(1).min(m.max(1));
    let rows = m.div_ceil(threads);
    std::thread::scope(|sc| {
        for (t, chunk) in out[..m * n].chunks_mut(rows * n).enumerate() {
            sc.spawn(move || {
                let mr = chunk.len() / n;
                product_nt(&a[t * rows * kd..(t * rows + mr) * kd], b, mr, n, kd, chunk);
            });
        }
    });
}

/// `sum_i c_i c_i^T` over the lowest `occupied` orbitals.
fn closed_shell_density(c: &[f64], n: usize, m: usize, occupied: usize) -> Matrix {
    let mut d = Matrix::zeros(n);
    for k in 0..occupied {
        for i in 0..n {
            let ci = c[i * m + k];
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

/// Closed-shell Hartree-Fock for `problem`, which must carry an auxiliary
/// basis and equal alpha and beta electrons. A ghost atom (a nucleus of zero
/// charge with its basis kept) is fine, which is what a counterpoise
/// correction needs.
pub fn hartree_fock(problem: &Problem, max_iterations: usize, tolerance: f64) -> HartreeFock {
    hartree_fock_in_field(problem, [0.0; 3], max_iterations, tolerance)
}

/// [`hartree_fock`] in a uniform electric field `field` (atomic units, the
/// force on a positive charge): the electrons feel `+ F . r` and the nuclei's
/// energy changes by `- sum Z F . R`, both about the coordinates' own origin,
/// so the total energy has a linear term in the field and a quadratic one
/// whose second difference is minus the polarisability. With no field the
/// arithmetic is exactly [`hartree_fock`]'s.
pub fn hartree_fock_in_field(problem: &Problem, field: [f64; 3], max_iterations: usize, tolerance: f64) -> HartreeFock {
    let basis = &problem.basis;
    let n = basis.size;
    let aux = problem.auxiliary.as_ref().expect("Hartree-Fock here is density fitted: the problem needs an auxiliary basis");
    assert!((problem.alpha - problem.beta).abs() < 1e-9 && (problem.alpha - problem.alpha.round()).abs() < 1e-9, "closed-shell Hartree-Fock only: {} alpha and {} beta electrons", problem.alpha, problem.beta);
    let occupied = problem.alpha.round() as usize;
    let (s, t, v) = one_electron(basis, &problem.nuclei);
    let mut h = t.clone();
    for k in 0..n * n {
        h.a[k] += v.a[k];
    }
    let mut e_field = 0.0;
    if field != [0.0; 3] {
        let dip = super::integrals::dipole(basis, [0.0; 3]);
        for d in 0..3 {
            for k in 0..n * n {
                h.a[k] += field[d] * dip[d].a[k];
            }
        }
        for (z, r) in &problem.nuclei {
            e_field -= z * (field[0] * r[0] + field[1] * r[1] + field[2] * r[2]);
        }
    }
    let (x, m) = orthogonaliser(&s, 1e-8);
    assert!(occupied <= m, "more occupied orbitals than the basis has directions");
    let fit = Fitted::new(basis, aux, &s);
    let mut e_nn = 0.0;
    for (i, (zi, pi)) in problem.nuclei.iter().enumerate() {
        for (zj, pj) in problem.nuclei.iter().take(i) {
            let r = ((pi[0] - pj[0]).powi(2) + (pi[1] - pj[1]).powi(2) + (pi[2] - pj[2]).powi(2)).sqrt();
            e_nn += zi * zj / r;
        }
    }
    // The first orbitals. From the bare nuclei, as the Kohn-Sham solver starts,
    // Hartree-Fock took 21 iterations for a water dimer; from the Fock matrix of
    // the free atoms' summed density, with its Coulomb potential from the fit
    // and a local exchange (LDA) standing in for the real one, which needs
    // orbitals it does not have yet, it takes fewer. The guess only has to put
    // the electrons roughly where they go: the energy it converges to is the
    // same.
    let (mut levels, mut c) = match &problem.guess {
        Some((da, db)) => {
            let mut dt = da.clone();
            for q in 0..n * n {
                dt.a[q] += db.a[q];
            }
            let (j0, _) = fit.coulomb_and_energy(&dt);
            let atoms: Vec<([f64; 3], f64)> = problem.nuclei.iter().zip(&problem.sizes).map(|((_, p), r)| (*p, *r)).collect();
            let grid = super::grid::molecular_pruned(&atoms, problem.radial, problem.theta, problem.prune);
            let batches = super::scf::Batches::new(basis, &grid);
            let (_, vxa, _) = super::scf::exchange_correlation(basis, &batches, super::functional::Functional::Lda, da, db);
            let mut f0 = h.clone();
            for q in 0..n * n {
                f0.a[q] += j0.a[q] + vxa.a[q];
            }
            generalised(&f0, &x, m)
        }
        None => generalised(&h, &x, m),
    };
    let mut hist: Vec<(Matrix, Matrix, Vec<f64>)> = Vec::new();
    let mut energy = 0.0;
    let mut last = f64::INFINITY;
    let mut converged = false;
    let mut iterations = 0;
    for it in 0..max_iterations {
        iterations = it + 1;
        let p = closed_shell_density(&c, n, m, occupied);
        let mut d = p.clone();
        for a in d.a.iter_mut() {
            *a *= 2.0;
        }
        let clock = std::time::Instant::now();
        let (j, ej) = fit.coulomb_and_energy(&d);
        let t_j = clock.elapsed().as_secs_f64();
        let (y, nk) = fit.whitened_half(&c, m, 0, occupied);
        let t_y = clock.elapsed().as_secs_f64() - t_j;
        let k = exchange(&y, n, occupied, nk);
        let t_k = clock.elapsed().as_secs_f64() - t_j - t_y;
        let mut f = h.clone();
        for q in 0..n * n {
            f.a[q] += j.a[q] - k.a[q];
        }
        energy = h.dot(&d) + ej - 0.5 * d.dot(&k) + e_nn + e_field;
        let err = {
            let fds = f.mul(&d).mul(&s);
            let sdf = s.mul(&d).mul(&f);
            (0..n * n).map(|q| fds.a[q] - sdf.a[q]).collect::<Vec<f64>>()
        };
        let emax = err.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        if ((energy - last).abs() < tolerance && emax < tolerance.sqrt()) || emax < 1e-9 {
            converged = true;
            break;
        }
        last = energy;
        hist.push((f.clone(), f.clone(), err));
        if hist.len() > 10 {
            hist.remove(0);
        }
        let fx = if hist.len() >= 2 { diis(&hist).map(|(a, _)| a).unwrap_or(f) } else { f };
        let clock = std::time::Instant::now();
        let (e, cc) = generalised(&fx, &x, m);
        levels = e;
        c = cc;
        if std::env::var_os("PHYS_PROFILE").is_some() {
            eprintln!("  HF iteration {}: Coulomb {t_j:.2} s, half transform {t_y:.2} s, exchange {t_k:.2} s, diagonalisation {:.2} s", it + 1, clock.elapsed().as_secs_f64());
        }
    }
    HartreeFock { energy, converged, iterations, c, m, levels, occupied }
}

/// The RI-MP2 correlation energy and its spin parts.
#[derive(Debug, Clone, Copy)]
pub struct Mp2 {
    /// `E_os + E_ss`, hartree.
    pub correlation: f64,
    /// Opposite-spin and same-spin parts, `sum (ia|jb)^2 / D` and
    /// `sum (ia|jb) [(ia|jb) - (ib|ja)] / D`. Reported so that anyone who
    /// wants the spin-component-scaled variants has what they need; the energy
    /// used is their plain sum.
    pub opposite_spin: f64,
    pub same_spin: f64,
}

/// RI-MP2 on a Hartree-Fock solution, the first `frozen` orbitals left out of
/// the correlation.
pub fn mp2(problem: &Problem, hf: &HartreeFock, frozen: usize) -> Mp2 {
    let basis = &problem.basis;
    let n = basis.size;
    let aux = problem.auxiliary.as_ref().expect("MP2 here is density fitted: the problem needs an auxiliary basis");
    let (s, _, _) = one_electron(basis, &problem.nuclei);
    let fit = Fitted::new(basis, aux, &s);
    let (m, nocc) = (hf.m, hf.occupied);
    assert!(frozen < nocc, "freezing {frozen} of {nocc} occupied orbitals leaves nothing to correlate");
    let active = nocc - frozen;
    let nv = m - nocc;
    let (y, nk) = fit.whitened_half(&hf.c, m, frozen, active);
    // Virtual orbitals as rows: cv[a][ao].
    let mut cv = vec![0.0f64; nv * n];
    for a in 0..nv {
        for i in 0..n {
            cv[a * n + i] = hf.c[i * m + nocc + a];
        }
    }
    // (ia|k) for each active occupied i: B[i][a][k] = sum_m cv[a][m] Y[i][m][k].
    let mut b = vec![0.0f64; active * nv * nk];
    for i in 0..active {
        let yi = &y[i * n * nk..(i + 1) * n * nk];
        // Y_i as [k][m] so the contraction over m is along rows.
        let mut yt = vec![0.0f64; nk * n];
        for mm in 0..n {
            for k in 0..nk {
                yt[k * n + mm] = yi[mm * nk + k];
            }
        }
        let bi = &mut b[i * nv * nk..(i + 1) * nv * nk];
        product_parallel(&cv, &yt, nv, nk, n, bi);
    }
    drop(y);
    let eo = &hf.levels[frozen..nocc];
    let ev = &hf.levels[nocc..m];
    let (mut os, mut ss) = (0.0f64, 0.0f64);
    let mut a_ij = vec![0.0f64; nv * nv];
    for i in 0..active {
        for j in i..active {
            let (bi, bj) = (&b[i * nv * nk..(i + 1) * nv * nk], &b[j * nv * nk..(j + 1) * nv * nk]);
            product_parallel(bi, bj, nv, nv, nk, &mut a_ij);
            let weight = if i == j { 1.0 } else { 2.0 };
            let (mut o, mut sm) = (0.0f64, 0.0f64);
            for a in 0..nv {
                for bb in 0..nv {
                    let iajb = a_ij[a * nv + bb];
                    let ibja = a_ij[bb * nv + a];
                    let d = eo[i] + eo[j] - ev[a] - ev[bb];
                    o += iajb * iajb / d;
                    sm += iajb * (iajb - ibja) / d;
                }
            }
            os += weight * o;
            ss += weight * sm;
        }
    }
    Mp2 { correlation: os + ss, opposite_spin: os, same_spin: ss }
}

/// The number of core orbitals to freeze for `z`: the 1s of each real atom
/// beyond helium. Atoms past neon are not handled.
pub fn frozen_core(z: &[u32], real: &[usize]) -> usize {
    real.iter()
        .map(|&i| match z[i] {
            0..=2 => 0,
            3..=10 => 1,
            other => panic!("no frozen-core rule for Z = {other}"),
        })
        .sum()
}

/// The static dipole polarisability tensors of the molecule in `problem`, by
/// finite field: minus the second derivative of the total energy with respect
/// to a uniform field, from central differences of `step` atomic units. Returns
/// the Hartree-Fock tensor and the MP2 one (the MP2 correlation energy at each
/// field's own Hartree-Fock orbitals added, which is the relaxed MP2
/// response). Diagonals from `E(+F) + E(-F) - 2 E(0)`, off-diagonals from the
/// four corners `E(+,+) - E(+,-) - E(-,+) + E(-,-)` over `4 F^2`: nineteen
/// solves, all in one basis and so one table. The axes are the box's; a
/// molecule not held in its principal axes gets the full tensor back. Bohr
/// cubed.
pub fn polarisabilities(problem: &Problem, frozen: usize, step: f64) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let energy = |f: [f64; 3]| -> (f64, f64) {
        let hf = hartree_fock_in_field(problem, f, 300, 1e-11);
        assert!(hf.converged, "Hartree-Fock in the field {f:?} did not converge");
        (hf.energy, hf.energy + mp2(problem, &hf, frozen).correlation)
    };
    let at = |a: f64, b: f64, c: f64| energy([a * step, b * step, c * step]);
    let e0 = at(0.0, 0.0, 0.0);
    let mut hf = [[0.0f64; 3]; 3];
    let mut mp = [[0.0f64; 3]; 3];
    for d in 0..3 {
        let mut plus = [0.0f64; 3];
        plus[d] = 1.0;
        let mut minus = [0.0f64; 3];
        minus[d] = -1.0;
        let (p, m) = (at(plus[0], plus[1], plus[2]), at(minus[0], minus[1], minus[2]));
        hf[d][d] = -(p.0 + m.0 - 2.0 * e0.0) / (step * step);
        mp[d][d] = -(p.1 + m.1 - 2.0 * e0.1) / (step * step);
    }
    for d in 0..3 {
        for e in d + 1..3 {
            let corner = |sd: f64, se: f64| {
                let mut f = [0.0; 3];
                f[d] = sd;
                f[e] = se;
                at(f[0], f[1], f[2])
            };
            let (pp, pm, mpl, mm) = (corner(1.0, 1.0), corner(1.0, -1.0), corner(-1.0, 1.0), corner(-1.0, -1.0));
            let vh = -(pp.0 - pm.0 - mpl.0 + mm.0) / (4.0 * step * step);
            let vm = -(pp.1 - pm.1 - mpl.1 + mm.1) / (4.0 * step * step);
            hf[d][e] = vh;
            hf[e][d] = vh;
            mp[d][e] = vm;
            mp[e][d] = vm;
        }
    }
    (hf, mp)
}
