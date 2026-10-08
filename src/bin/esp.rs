//! A molecule's charges from its own electrons, with no pair energy in them —
//! `docs/PLAY.md` E8b.
//!
//! ```sh
//! cargo run --release --bin phys-esp -- water
//! ```
//!
//! The molecule at the shape `phys-grow` relaxed it to, Hartree-Fock, and its
//! electrostatic potential (nuclei less electrons) on Merz-Kollman shells
//! outside the van der Waals surface: 1.4, 1.6, 1.8 and 2.0 times each atom's
//! radius, points kept where no atom's 1.4 radii reaches. One point charge
//! per symmetry class of atom and, if asked, one more on the bisector of the
//! first atom's two bonds at a distance searched for, are fitted to that
//! potential under neutrality. Nothing here knows what a dimer's energy is;
//! these are the charges the density says, and the law's repulsion and
//! dispersion are fitted to the pair energies around them (`phys-fit --esp`).
//! Writes `esp-<name>.txt`: `charges` and `bisector` lines as a law's text has.

use phys::electrons::functional::Functional;
use phys::electrons::grow::{equivalent_atoms, Resume};
use phys::electrons::hf::hartree_fock;
use phys::electrons::integrals::one_electron;
use phys::electrons::molecule::Molecule;
use phys::math::Vec3;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let z: Vec<u32> = match name.as_str() {
        "water" => vec![8, 1, 1],
        "methane" => vec![6, 1, 1, 1, 1],
        "ammonia" => vec![7, 1, 1, 1],
        "methanol" => vec![6, 8, 1, 1, 1, 1],
        other => panic!("no atoms known for {other}"),
    };
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a growth state")).expect("readable");
    let mol = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let classes = equivalent_atoms(&mol);
    let p = mol.problem(Functional::Pbe);
    let hf = hartree_fock(&p, 100, 1e-10);
    let n = p.basis.size;
    // The density matrix, 2 C C^T over the occupied orbitals.
    let mut d = phys::electrons::linalg::Matrix::zeros(n);
    for a in 0..n {
        for b in 0..n {
            let s: f64 = (0..hf.occupied).map(|i| hf.c[a * hf.m + i] * hf.c[b * hf.m + i]).sum();
            d.set(a, b, 2.0 * s);
        }
    }
    let pos: Vec<Vec3> = state.positions.iter().map(|r| Vec3 { x: r[0], y: r[1], z: r[2] }).collect();
    let radius = |zz: u32| phys::chem::elements::Element(zz as u8).vdw_radius().unwrap_or(1.5e-10) / 0.529177210903e-10;
    // Merz-Kollman points.
    let mut pts: Vec<Vec3> = Vec::new();
    let gauss = args.iter().any(|a| a == "--gauss");
    let shells: Vec<f64> = if gauss { vec![1.0, 1.15, 1.3, 1.5, 1.8, 2.2] } else { vec![1.4, 1.6, 1.8, 2.0] };
    let inner = if gauss { 1.0 } else { 1.4 };
    for scale in shells {
        for (i, c) in pos.iter().enumerate() {
            let r = radius(z[i]) * scale;
            let m = 400;
            for k in 0..m {
                let t = (k as f64 + 0.5) / m as f64;
                let phi = std::f64::consts::PI * (1.0 + 5f64.sqrt()) * k as f64;
                let (ct, st) = (1.0 - 2.0 * t, (1.0 - (1.0 - 2.0 * t).powi(2)).sqrt());
                let q = *c + Vec3 { x: st * phi.cos(), y: st * phi.sin(), z: ct }.scale(r);
                if pos.iter().enumerate().all(|(j, o)| (q - *o).norm() >= inner * radius(z[j])) {
                    pts.push(q);
                }
            }
        }
    }
    let esp: Vec<f64> = pts
        .iter()
        .map(|q| {
            let nuc: f64 = pos.iter().zip(&z).map(|(o, zz)| *zz as f64 / (*q - *o).norm()).sum();
            let (_, _, v) = one_electron(&p.basis, &[(1.0, [q.x, q.y, q.z])]);
            let mut e = 0.0;
            for a in 0..n {
                for b in 0..n {
                    e += d.get(a, b) * v.get(a, b);
                }
            }
            nuc + e
        })
        .collect();
    println!("{name}: {} functions, {} potential points", n, pts.len());
    let types = classes.iter().max().map(|m| m + 1).unwrap_or(1);
    let mut mult = vec![0usize; types];
    for c in &classes {
        mult[*c] += 1;
    }
    // For a bisector distance (None: no extra site), the least-squares charges
    // under neutrality, by the normal equations with a Lagrange multiplier.
    let solve = |dist: Option<f64>| -> (Vec<f64>, f64) {
        let sites = types + dist.is_some() as usize;
        let site_pos = |i: usize| -> Vec<Vec3> {
            if i < types {
                (0..z.len()).filter(|a| classes[*a] == i).map(|a| pos[a]).collect()
            } else {
                let o = pos[0];
                let (u, v) = ((pos[1] - o).unit(), (pos[2] - o).unit());
                vec![o + (u + v).unit().scale(dist.unwrap())]
            }
        };
        let rows: Vec<Vec<f64>> = pts.iter().map(|q| (0..sites).map(|i| site_pos(i).iter().map(|s| 1.0 / (*q - *s).norm()).sum()).collect()).collect();
        let m = sites + 1;
        let mut a = vec![vec![0.0; m]; m];
        let mut b = vec![0.0; m];
        for (row, v) in rows.iter().zip(&esp) {
            for i in 0..sites {
                b[i] += row[i] * v;
                for j in 0..sites {
                    a[i][j] += row[i] * row[j];
                }
            }
        }
        for i in 0..sites {
            let mi = if i < types { mult[i] as f64 } else { 1.0 };
            a[i][sites] = mi;
            a[sites][i] = mi;
        }
        // Gaussian elimination.
        for c in 0..m {
            let piv = (c..m).max_by(|&x, &y| a[x][c].abs().partial_cmp(&a[y][c].abs()).unwrap()).unwrap();
            a.swap(c, piv);
            b.swap(c, piv);
            for r in c + 1..m {
                let f = a[r][c] / a[c][c];
                for k in c..m {
                    a[r][k] -= f * a[c][k];
                }
                b[r] -= f * b[c];
            }
        }
        let mut x = vec![0.0; m];
        for c in (0..m).rev() {
            x[c] = (b[c] - (c + 1..m).map(|k| a[c][k] * x[k]).sum::<f64>()) / a[c][c];
        }
        let q: Vec<f64> = x[..sites].to_vec();
        let rms = (rows.iter().zip(&esp).map(|(row, v)| (row.iter().zip(&q).map(|(r, qq)| r * qq).sum::<f64>() - v).powi(2)).sum::<f64>() / esp.len() as f64).sqrt();
        (q, rms)
    };
    let (q0, r0) = solve(None);
    println!("  atom-centred charges {q0:?}: potential RMS {r0:.3e} hartree/e");
    let mut best = (f64::INFINITY, 0.0, vec![]);
    for k in 0..=24 {
        let dist = 0.02 * k as f64 + 0.02;
        let (q, r) = solve(Some(dist));
        println!("  bisector site at {dist:.2} bohr: charges {q:?}, RMS {r:.3e}");
        if r < best.0 {
            best = (r, dist, q);
        }
    }
    if gauss {
        // Gaussian charges: a width per type (the bisector site's own), the
        // charges linear for given widths and distance, so Nelder-Mead over
        // widths and distance and a constrained least squares inside.
        let site_pos = |i: usize, dist: f64| -> Vec<Vec3> {
            if i < types {
                (0..z.len()).filter(|a| classes[*a] == i).map(|a| pos[a]).collect()
            } else {
                let o = pos[0];
                let (u, v) = ((pos[1] - o).unit(), (pos[2] - o).unit());
                vec![o + (u + v).unit().scale(dist)]
            }
        };
        let sites = types + 1;
        let fit_at = |x: &[f64]| -> (Vec<f64>, f64) {
            // x: ln sigma per site (types + 1), then the distance.
            let sig: Vec<f64> = x[..sites].iter().map(|v| v.exp()).collect();
            let dist = x[sites].abs().max(1e-3);
            let rows: Vec<Vec<f64>> = pts.iter().map(|q| (0..sites).map(|i| site_pos(i, dist).iter().map(|s| { let r = (*q - *s).norm(); phys::induction::erf(r / sig[i]) / r }).sum()).collect()).collect();
            let m = sites + 1;
            let mut a = vec![vec![0.0; m]; m];
            let mut b = vec![0.0; m];
            for (row, v) in rows.iter().zip(&esp) {
                for i in 0..sites {
                    b[i] += row[i] * v;
                    for j in 0..sites {
                        a[i][j] += row[i] * row[j];
                    }
                }
            }
            for i in 0..sites {
                let mi = if i < types { mult[i] as f64 } else { 1.0 };
                a[i][sites] = mi;
                a[sites][i] = mi;
            }
            for cc in 0..m {
                let piv = (cc..m).max_by(|&u, &w| a[u][cc].abs().partial_cmp(&a[w][cc].abs()).unwrap()).unwrap();
                a.swap(cc, piv);
                b.swap(cc, piv);
                for r in cc + 1..m {
                    let f = a[r][cc] / a[cc][cc];
                    for k in cc..m {
                        a[r][k] -= f * a[cc][k];
                    }
                    b[r] -= f * b[cc];
                }
            }
            let mut sol = vec![0.0; m];
            for cc in (0..m).rev() {
                sol[cc] = (b[cc] - (cc + 1..m).map(|k| a[cc][k] * sol[k]).sum::<f64>()) / a[cc][cc];
            }
            let q: Vec<f64> = sol[..sites].to_vec();
            let rms = (rows.iter().zip(&esp).map(|(row, v)| (row.iter().zip(&q).map(|(r, qq)| r * qq).sum::<f64>() - v).powi(2)).sum::<f64>() / esp.len() as f64).sqrt();
            // A small pull of the widths toward 1 bohr and of charges toward
            // reasonable size keeps the near-degenerate directions tame.
            let pen = 1e-6 * q.iter().map(|v| v * v).sum::<f64>();
            (q, rms + pen)
        };
        // Nelder-Mead.
        let n = sites + 1;
        let mut simplex: Vec<Vec<f64>> = Vec::new();
        let base: Vec<f64> = (0..sites).map(|_| (0.8f64).ln()).chain(std::iter::once(0.4)).collect();
        simplex.push(base.clone());
        for k in 0..n {
            let mut v = base.clone();
            v[k] += if k < sites { 0.4 } else { 0.15 };
            simplex.push(v);
        }
        let mut vals: Vec<f64> = simplex.iter().map(|v| fit_at(v).1).collect();
        for _ in 0..600 {
            let mut idx: Vec<usize> = (0..=n).collect();
            idx.sort_by(|&a, &b| vals[a].partial_cmp(&vals[b]).unwrap());
            simplex = idx.iter().map(|&i| simplex[i].clone()).collect();
            vals = idx.iter().map(|&i| vals[i]).collect();
            let centroid: Vec<f64> = (0..n).map(|k| simplex[..n].iter().map(|v| v[k]).sum::<f64>() / n as f64).collect();
            let worst = simplex[n].clone();
            let along = |t: f64| -> Vec<f64> { (0..n).map(|k| centroid[k] + t * (worst[k] - centroid[k])).collect() };
            let refl = along(-1.0);
            let fr = fit_at(&refl).1;
            if fr < vals[0] {
                let exp = along(-2.0);
                let fe = fit_at(&exp).1;
                if fe < fr { simplex[n] = exp; vals[n] = fe; } else { simplex[n] = refl; vals[n] = fr; }
            } else if fr < vals[n - 1] {
                simplex[n] = refl;
                vals[n] = fr;
            } else {
                let con = along(0.5);
                let fc = fit_at(&con).1;
                if fc < vals[n] {
                    simplex[n] = con;
                    vals[n] = fc;
                } else {
                    for i in 1..=n {
                        for k in 0..n {
                            simplex[i][k] = simplex[0][k] + 0.5 * (simplex[i][k] - simplex[0][k]);
                        }
                        vals[i] = fit_at(&simplex[i]).1;
                    }
                }
            }
        }
        let best_x = simplex[0].clone();
        let (q, rms) = fit_at(&best_x);
        let sig: Vec<f64> = best_x[..sites].iter().map(|v| v.exp()).collect();
        let dist = best_x[sites].abs().max(1e-3);
        println!("  Gaussian charges: bisector {dist:.3} bohr, charges {q:?}, widths {sig:?} bohr, potential RMS {rms:.3e} hartree/e (point charges: {:.3e})", best.0);
        let text = format!("charges {}
sigma {}
bisector {} {:e}
", q.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "), sig.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "), types, dist);
        std::fs::write(format!("esp-{name}-gauss.txt"), text).expect("written");
        return;
    }
    let dip = |q: &[f64], dist: f64| -> f64 {
        let o = pos[0];
        let (u, v) = ((pos[1] - o).unit(), (pos[2] - o).unit());
        let site = o + (u + v).unit().scale(dist);
        let mut m = Vec3::ZERO;
        for (a, c) in classes.iter().enumerate() {
            m = m + pos[a].scale(q[*c]);
        }
        m = m + site.scale(q[types]);
        m.norm() / 0.393430307
    };
    println!("  best: bisector {:.2} bohr, charges {:?}, RMS {:.3e}, dipole {:.3} D", best.1, best.2, best.0, dip(&best.2, best.1));
    let text = format!("charges {}\nbisector {} {:e}\n", best.2.iter().map(|q| format!("{q:e}")).collect::<Vec<_>>().join(" "), types, best.1);
    std::fs::write(format!("esp-{name}.txt"), text).expect("written");
}
