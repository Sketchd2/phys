//! Hartree-Fock and MP2 by density fitting (`src/electrons/hf.rs`), checked
//! against the same theory done with exact four-index integrals on molecules
//! small enough to do it that way.

use phys::electrons::functional::Functional;
use phys::electrons::hf::{frozen_core, hartree_fock, mp2};
use phys::electrons::integrals::one_electron;
use phys::electrons::linalg::{generalised, orthogonaliser, Matrix};
use phys::electrons::molecule::Molecule;
use phys::electrons::scf::{Problem, Repulsion};

fn molecule(z: &[u32], pos: &[[f64; 3]]) -> Molecule {
    Molecule { z: z.to_vec(), positions: pos.to_vec(), charge: 0, unpaired: 0 }
}

/// Closed-shell Hartree-Fock with exact integrals, by plain damped iteration:
/// small systems only.
fn exact_hartree_fock(p: &Problem, eri: &Repulsion) -> (f64, Vec<f64>, Vec<f64>, usize) {
    let n = p.basis.size;
    let (s, t, v) = one_electron(&p.basis, &p.nuclei);
    let mut h = t.clone();
    for k in 0..n * n {
        h.a[k] += v.a[k];
    }
    let (x, m) = orthogonaliser(&s, 1e-8);
    let occ = p.alpha.round() as usize;
    let mut e_nn = 0.0;
    for (i, (zi, pi)) in p.nuclei.iter().enumerate() {
        for (zj, pj) in p.nuclei.iter().take(i) {
            e_nn += zi * zj / ((pi[0] - pj[0]).powi(2) + (pi[1] - pj[1]).powi(2) + (pi[2] - pj[2]).powi(2)).sqrt();
        }
    }
    let (mut levels, mut c) = generalised(&h, &x, m);
    let mut energy = 0.0;
    let mut dp = Matrix::zeros(n);
    for it in 0..300 {
        let mut pm = Matrix::zeros(n);
        for k in 0..occ {
            for i in 0..n {
                for j in 0..n {
                    pm.a[i * n + j] += c[i * m + k] * c[j * m + k];
                }
            }
        }
        if it > 0 {
            for q in 0..n * n {
                pm.a[q] = 0.5 * pm.a[q] + 0.5 * dp.a[q];
            }
        }
        dp = pm.clone();
        let mut d = pm.clone();
        for a in d.a.iter_mut() {
            *a *= 2.0;
        }
        let j = eri.coulomb(&d);
        let mut kx = Matrix::zeros(n);
        for a in 0..n {
            for b in 0..n {
                let mut sum = 0.0;
                for l in 0..n {
                    for sg in 0..n {
                        sum += pm.get(l, sg) * eri.get(a, l, b, sg);
                    }
                }
                kx.set(a, b, sum);
            }
        }
        let mut f = h.clone();
        for q in 0..n * n {
            f.a[q] += j.a[q] - kx.a[q];
        }
        let e_new = h.dot(&d) + 0.5 * j.dot(&d) - 0.5 * kx.dot(&d) + e_nn;
        let (l, cc) = generalised(&f, &x, m);
        levels = l;
        c = cc;
        if (e_new - energy).abs() < 1e-11 && it > 3 {
            energy = e_new;
            break;
        }
        energy = e_new;
    }
    (energy, c, levels, m)
}

/// MP2 with exact integrals on given orbitals: the full transformation.
fn exact_mp2(p: &Problem, eri: &Repulsion, c: &[f64], levels: &[f64], m: usize, frozen: usize) -> f64 {
    let n = p.basis.size;
    let nocc = p.alpha.round() as usize;
    // (ia|jb) for active i, j and virtual a, b, by four quarter transformations.
    let occ: Vec<usize> = (frozen..nocc).collect();
    let vir: Vec<usize> = (nocc..m).collect();
    let (no, nv) = (occ.len(), vir.len());
    // t1[i][ns-pair...]: do it plainly, n is tiny.
    let mut q1 = vec![0.0f64; no * n * n * n];
    for (ii, &i) in occ.iter().enumerate() {
        for nu in 0..n {
            for la in 0..n {
                for si in 0..n {
                    let mut sum = 0.0;
                    for mu in 0..n {
                        sum += c[mu * m + i] * eri.get(mu, nu, la, si);
                    }
                    q1[((ii * n + nu) * n + la) * n + si] = sum;
                }
            }
        }
    }
    let mut q2 = vec![0.0f64; no * nv * n * n];
    for ii in 0..no {
        for (aa, &a) in vir.iter().enumerate() {
            for la in 0..n {
                for si in 0..n {
                    let mut sum = 0.0;
                    for nu in 0..n {
                        sum += c[nu * m + a] * q1[((ii * n + nu) * n + la) * n + si];
                    }
                    q2[((ii * nv + aa) * n + la) * n + si] = sum;
                }
            }
        }
    }
    let mut q3 = vec![0.0f64; no * nv * no * n];
    for ii in 0..no {
        for aa in 0..nv {
            for (jj, &j) in occ.iter().enumerate() {
                for si in 0..n {
                    let mut sum = 0.0;
                    for la in 0..n {
                        sum += c[la * m + j] * q2[((ii * nv + aa) * n + la) * n + si];
                    }
                    q3[((ii * nv + aa) * no + jj) * n + si] = sum;
                }
            }
        }
    }
    let mut iajb = vec![0.0f64; no * nv * no * nv];
    for ii in 0..no {
        for aa in 0..nv {
            for jj in 0..no {
                for (bb, &b) in vir.iter().enumerate() {
                    let mut sum = 0.0;
                    for si in 0..n {
                        sum += c[si * m + b] * q3[((ii * nv + aa) * no + jj) * n + si];
                    }
                    iajb[((ii * nv + aa) * no + jj) * nv + bb] = sum;
                }
            }
        }
    }
    let mut e2 = 0.0;
    for ii in 0..no {
        for jj in 0..no {
            for aa in 0..nv {
                for bb in 0..nv {
                    let v1 = iajb[((ii * nv + aa) * no + jj) * nv + bb];
                    let v2 = iajb[((ii * nv + bb) * no + jj) * nv + aa];
                    let d = levels[occ[ii]] + levels[occ[jj]] - levels[vir[aa]] - levels[vir[bb]];
                    e2 += v1 * (2.0 * v1 - v2) / d;
                }
            }
        }
    }
    e2
}

/// The fitted Hartree-Fock energy is the exact one to the fit, and RI-MP2 is the
/// full four-index transformation's on the same orbitals to a small fraction of
/// the correlation energy: what the auxiliary basis costs. (The fitted Coulomb
/// energy is a lower bound in this metric, exchange has no such bound, so only
/// the size is tested.) One exact repulsion tensor per molecule, which is the
/// expensive part.
#[test]
fn fitted_theory_agrees_with_exact_integrals() {
    for (name, z, pos) in [("He", vec![2u32], vec![[0.0; 3]]), ("H2", vec![1, 1], vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])] {
        let p = molecule(&z, &pos).problem(Functional::Pbe);
        let eri = Repulsion::new(&p.basis);
        let fitted = hartree_fock(&p, 100, 1e-10);
        let (exact, ..) = exact_hartree_fock(&p, &eri);
        println!("  {name}: {} functions, Hartree-Fock fitted {:.8} ({} iterations), exact {:.8}, difference {:.2e} Ha", p.basis.size, fitted.energy, fitted.iterations, exact, fitted.energy - exact);
        assert!(fitted.converged);
        assert!((fitted.energy - exact).abs() < 2e-4, "{name}: the fit moved the energy by {:.2e}", fitted.energy - exact);
        let ri = mp2(&p, &fitted, 0);
        let four = exact_mp2(&p, &eri, &fitted.c, &fitted.levels, fitted.m, 0);
        println!("  {name}: RI-MP2 {:.8}, four-index {:.8} on the same orbitals; difference {:.2e} Ha ({:.2}%)", ri.correlation, four, ri.correlation - four, 100.0 * (ri.correlation - four) / four);
        assert!(ri.correlation < 0.0 && four < 0.0);
        assert!((ri.correlation - four).abs() < 0.02 * four.abs(), "{name}: RI-MP2 is off the exact by {:.1}%", 100.0 * (ri.correlation - four) / four);
    }
}

/// The frozen-core count.
#[test]
fn frozen_core_counts_the_first_shell_of_real_atoms() {
    let z = [8u32, 1, 1, 8, 1, 1];
    assert_eq!(frozen_core(&z, &[0, 1, 2]), 1, "one water of the dimer");
    assert_eq!(frozen_core(&z, &[0, 1, 2, 3, 4, 5]), 2);
    assert_eq!(frozen_core(&[2, 1], &[0, 1]), 0, "helium and hydrogen have no core to freeze");
}
