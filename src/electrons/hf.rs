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
    let threads = std::thread::available_parallelism().map(|t| t.get()).unwrap_or(1).min(n.max(1));
    let rows = n.div_ceil(threads);
    let mut k = Matrix::zeros(n);
    std::thread::scope(|sc| {
        for (t, chunk) in k.a.chunks_mut(rows * n).enumerate() {
            sc.spawn(move || {
                let m = chunk.len() / n;
                let mut tmp = vec![0.0f64; m * n];
                for i in 0..count {
                    let yi = &y[i * n * nk..(i + 1) * n * nk];
                    product_nt(&yi[t * rows * nk..(t * rows + m) * nk], yi, m, n, nk, &mut tmp);
                    for (a, b) in chunk.iter_mut().zip(&tmp) {
                        *a += b;
                    }
                }
            });
        }
    });
    k
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
    // Start from the bare-nucleus Hamiltonian, as the Kohn-Sham solver does.
    let (mut levels, mut c) = generalised(&h, &x, m);
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
        let (j, ej) = fit.coulomb_and_energy(&d);
        let (y, nk) = fit.whitened_half(&c, m, 0, occupied);
        let k = exchange(&y, n, occupied, nk);
        let mut f = h.clone();
        for q in 0..n * n {
            f.a[q] += j.a[q] - k.a[q];
        }
        energy = h.dot(&d) + ej - 0.5 * d.dot(&k) + e_nn;
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
        let (e, cc) = generalised(&fx, &x, m);
        levels = e;
        c = cc;
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
    let threads = std::thread::available_parallelism().map(|t| t.get()).unwrap_or(1);
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
        let rows = nv.div_ceil(threads);
        let (cv_ref, yt_ref) = (&cv, &yt);
        std::thread::scope(|sc| {
            for (t, chunk) in bi.chunks_mut(rows * nk).enumerate() {
                sc.spawn(move || {
                    let mr = chunk.len() / nk;
                    product_nt(&cv_ref[t * rows * n..(t * rows + mr) * n], yt_ref, mr, nk, n, chunk);
                });
            }
        });
    }
    drop(y);
    let eo = &hf.levels[frozen..nocc];
    let ev = &hf.levels[nocc..m];
    let (mut os, mut ss) = (0.0f64, 0.0f64);
    let mut a_ij = vec![0.0f64; nv * nv];
    for i in 0..active {
        for j in i..active {
            let (bi, bj) = (&b[i * nv * nk..(i + 1) * nv * nk], &b[j * nv * nk..(j + 1) * nv * nk]);
            let rows = nv.div_ceil(threads);
            std::thread::scope(|sc| {
                for (t, chunk) in a_ij.chunks_mut(rows * nv).enumerate() {
                    sc.spawn(move || {
                        let mr = chunk.len() / nv;
                        product_nt(&bi[t * rows * nk..(t * rows + mr) * nk], bj, mr, nv, nk, chunk);
                    });
                }
            });
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
