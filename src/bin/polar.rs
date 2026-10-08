//! A molecule's polarisability from its own electrons, and the atomic
//! polarisabilities that reproduce it — `docs/PLAY.md` Phase 6, E8a.
//!
//! ```sh
//! cargo run --release --bin phys-polar -- water [field step, default 0.002]
//! ```
//!
//! The molecule at the shape its growth relaxed it to (`grow-<name>.state`),
//! Hartree-Fock and MP2 energies in 19 uniform fields, the tensor from their
//! second differences (`electrons::hf::polarisabilities`). Then the induced-
//! dipole model's atomic polarisabilities, one per symmetry class of atom, found
//! by least squares so that the model's molecular tensor
//! (`induction::molecular_polarisability`: the atoms' dipoles coupled by the
//! screened dipole-dipole tensor) is the MP2 tensor. Writes
//! `polar-<name>.txt`: one `alpha <class> <value>` line per class, which a law's
//! text may carry.

use phys::electrons::functional::Functional;
use phys::electrons::grow::{equivalent_atoms, Resume};
use phys::electrons::hf::{frozen_core, polarisabilities};
use phys::electrons::molecule::Molecule;
use phys::induction::{molecular_polarisability, PolSite};
use phys::math::Vec3;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let step: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.002);
    let z: Vec<u32> = match name.as_str() {
        "water" => vec![8, 1, 1],
        "methane" => vec![6, 1, 1, 1, 1],
        "ammonia" => vec![7, 1, 1, 1],
        "methanol" => vec![6, 8, 1, 1, 1, 1],
        other => panic!("no atoms known for {other}"),
    };
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a growth state")).expect("readable");
    let mol = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let types = equivalent_atoms(&mol);
    let p = mol.problem(Functional::Pbe);
    let all: Vec<usize> = (0..z.len()).collect();
    let frozen = frozen_core(&z, &all);
    println!("{name}: {} functions, field step {step}, {frozen} frozen cores", p.basis.size);
    let t = std::time::Instant::now();
    let (hf, mp2) = polarisabilities(&p, frozen, step);
    println!("  19 fields in {:.0} s", t.elapsed().as_secs_f64());
    let mean = |a: &[[f64; 3]; 3]| (a[0][0] + a[1][1] + a[2][2]) / 3.0;
    for (label, a) in [("Hartree-Fock", &hf), ("MP2", &mp2)] {
        println!("  {label}: mean {:.4} bohr^3, diagonal {:.4} {:.4} {:.4}, off-diagonal {:.1e} {:.1e} {:.1e}", mean(a), a[0][0], a[1][1], a[2][2], a[0][1], a[0][2], a[1][2]);
    }
    // Atomic polarisabilities, one per class, by Levenberg-Marquardt on logs.
    let classes = types.iter().max().map(|m| m + 1).unwrap_or(1);
    let sites = |alpha: &[f64]| -> Vec<PolSite> { state.positions.iter().zip(&types).map(|(r, t)| PolSite { pos: Vec3 { x: r[0], y: r[1], z: r[2] }, alpha: alpha[*t] }).collect() };
    let residual = |logs: &[f64]| -> Vec<f64> {
        let alpha: Vec<f64> = logs.iter().map(|l| l.exp()).collect();
        let m = molecular_polarisability(&sites(&alpha));
        (0..3).flat_map(|i| (0..3).map(move |j| (i, j))).map(|(i, j)| m[i][j] - mp2[i][j]).collect()
    };
    let mut x = vec![(mean(&mp2) / z.len() as f64).ln(); classes];
    let cost = |r: &[f64]| r.iter().map(|v| v * v).sum::<f64>();
    let mut r = residual(&x);
    let mut c = cost(&r);
    let mut lambda = 1e-3;
    for _ in 0..200 {
        let mut jac = vec![vec![0.0; classes]; r.len()];
        for q in 0..classes {
            let (mut xp, mut xm) = (x.clone(), x.clone());
            xp[q] += 1e-6;
            xm[q] -= 1e-6;
            let (rp, rm) = (residual(&xp), residual(&xm));
            for i in 0..r.len() {
                jac[i][q] = (rp[i] - rm[i]) / 2e-6;
            }
        }
        let mut improved = false;
        for _ in 0..30 {
            // (J^T J + lambda diag) dx = J^T r, by Gaussian elimination.
            let mut a = vec![0.0; classes * classes];
            let mut b = vec![0.0; classes];
            for p_ in 0..classes {
                for q in 0..classes {
                    a[p_ * classes + q] = (0..r.len()).map(|i| jac[i][p_] * jac[i][q]).sum::<f64>();
                }
                a[p_ * classes + p_] *= 1.0 + lambda;
                b[p_] = (0..r.len()).map(|i| jac[i][p_] * r[i]).sum::<f64>();
            }
            for col in 0..classes {
                for row in col + 1..classes {
                    let f = a[row * classes + col] / a[col * classes + col];
                    for k in col..classes {
                        a[row * classes + k] -= f * a[col * classes + k];
                    }
                    b[row] -= f * b[col];
                }
            }
            let mut dx = vec![0.0; classes];
            for row in (0..classes).rev() {
                let mut s = b[row];
                for k in row + 1..classes {
                    s -= a[row * classes + k] * dx[k];
                }
                dx[row] = s / a[row * classes + row];
            }
            let trial: Vec<f64> = x.iter().zip(&dx).map(|(a, d)| a - d).collect();
            let rt = residual(&trial);
            if cost(&rt) < c {
                x = trial;
                r = rt;
                c = cost(&r);
                lambda = (lambda / 3.0).max(1e-12);
                improved = true;
                break;
            }
            lambda *= 4.0;
        }
        if !improved {
            break;
        }
    }
    let alpha: Vec<f64> = x.iter().map(|l| l.exp()).collect();
    let fitted = molecular_polarisability(&sites(&alpha));
    println!("  atomic polarisabilities by class {alpha:?} bohr^3 (classes {types:?}); the model's tensor: mean {:.4}, diagonal {:.4} {:.4} {:.4}; residual rms {:.2e}", mean(&fitted), fitted[0][0], fitted[1][1], fitted[2][2], (c / 9.0).sqrt());
    let mut text = String::new();
    for (t, a) in alpha.iter().enumerate() {
        text += &format!("alpha {t} {a:e}\n");
    }
    std::fs::write(format!("polar-{name}.txt"), text).expect("the polarisabilities written");
    println!("  written to polar-{name}.txt");
}
