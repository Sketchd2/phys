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
use phys::electrons::hf::{frozen_core, hartree_fock, hartree_fock_in_field, mp2};
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
    // The dipole Hartree-Fock's density gives, and the one MP2's energy gives
    // as minus its derivative with respect to a uniform field (the dipole of
    // the correlated, relaxed density). Hartree-Fock overestimates water's by
    // about 8%; the potential is scaled by the ratio when `--mp2-dipole` is
    // given, so the charges answer to the same method as the pair energies.
    let origin = [0.0; 3];
    let dm = phys::electrons::integrals::dipole(&p.basis, origin);
    let mut mu_hf = [0.0f64; 3];
    for k in 0..3 {
        let nuc: f64 = pos.iter().zip(&z).map(|(o, zz)| *zz as f64 * [o.x, o.y, o.z][k]).sum();
        let mut el = 0.0;
        for a in 0..n {
            for b in 0..n {
                el += d.get(a, b) * dm[k].get(a, b);
            }
        }
        mu_hf[k] = nuc - el;
    }
    let norm3 = |v: &[f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let scale = if args.iter().any(|a| a == "--mp2-dipole") {
        let all: Vec<usize> = (0..z.len()).collect();
        let fc = frozen_core(&z, &all);
        let step = 0.002;
        let mut mu_ff = [0.0f64; 3];
        let mut mu_ff_hf = [0.0f64; 3];
        for k in 0..3 {
            let mut e = [(0.0, 0.0); 2];
            for (s, sgn) in [1.0, -1.0].iter().enumerate() {
                let mut f = [0.0; 3];
                f[k] = sgn * step;
                let h = hartree_fock_in_field(&p, f, 100, 1e-11);
                let m = mp2(&p, &h, fc);
                e[s] = (h.energy, h.energy + m.correlation);
            }
            mu_ff_hf[k] = -(e[0].0 - e[1].0) / (2.0 * step);
            mu_ff[k] = -(e[0].1 - e[1].1) / (2.0 * step);
        }
        // The sign convention of the field is checked against the density's own
        // dipole, which it must reproduce.
        let sign = if (0..3).map(|k| mu_ff_hf[k] * mu_hf[k]).sum::<f64>() >= 0.0 { 1.0 } else { -1.0 };
        let mu_mp2 = [sign * mu_ff[0], sign * mu_ff[1], sign * mu_ff[2]];
        println!("  dipole: Hartree-Fock from the density {:.4} au ({:.3} D), from the field {:.4} au; MP2 from the field {:.4} au ({:.3} D)", norm3(&mu_hf), norm3(&mu_hf) / 0.393430307, norm3(&mu_ff_hf), norm3(&mu_mp2), norm3(&mu_mp2) / 0.393430307);
        norm3(&mu_mp2) / norm3(&mu_hf)
    } else {
        1.0
    };
    println!("  potential scaled by {scale:.4}");
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
            scale * (nuc + e)
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
    if args.iter().any(|a| a == "--sites") {
        // Sites where the molecule's own localised orbitals put a bond and a
        // lone pair (Boys), charges and widths fitted to the potential: the
        // places are derived, only the charge and size of each kind is fitted.
        let dip = phys::electrons::integrals::dipole(&p.basis, [0.0; 3]);
        let loc = phys::electrons::localise::boys(&hf.c, hf.m, n, 1..hf.occupied, &dip);
        let (o, p1, p2) = (pos[0], pos[1], pos[2]);
        let (u, v) = ((p1 - o).unit(), (p2 - o).unit());
        let nrm = u.cross(v).unit();
        // Frame coordinates (a, b, c) of a point: o + a u + b v + c n, by Cramer's rule.
        let det = |x: Vec3, y: Vec3, w: Vec3| x.dot(y.cross(w));
        let frame_of = |q: Vec3| -> [f64; 3] {
            let r = q - o;
            let d0 = det(u, v, nrm);
            [det(r, v, nrm) / d0, det(u, r, nrm) / d0, det(u, v, r) / d0]
        };
        let near_h = |q: Vec3| (q - p1).norm().min((q - p2).norm());
        let mut lone: Vec<[f64; 3]> = Vec::new();
        let mut bond: Vec<[f64; 3]> = Vec::new();
        for cen in &loc.centroids {
            let q = Vec3 { x: cen[0], y: cen[1], z: cen[2] };
            if near_h(q) < 1.6 {
                bond.push(frame_of(q))
            } else {
                lone.push(frame_of(q))
            }
        }
        println!("  localised in {} sweeps: {} bond centroids, {} lone pairs", loc.sweeps, bond.len(), lone.len());
        assert_eq!(lone.len(), 2, "water has two lone pairs");
        // Symmetrised: the two lone pairs are mirror images through the plane
        // that swaps the hydrogens (u and v exchange, the normal flips).
        let lp = {
            let (x, y) = (lone[0], lone[1]);
            let (pos_c, neg_c) = if x[2] >= y[2] { (x, y) } else { (y, x) };
            let a = (pos_c[0] + pos_c[1] + neg_c[0] + neg_c[1]) / 4.0;
            let cc = (pos_c[2] - neg_c[2]) / 2.0;
            [a, a, cc]
        };
        let want_bond = args.iter().any(|a| a == "--bond-sites");
        let bc_a = {
            let (x, y) = (bond[0], bond[1]);
            let (first, second) = if x[0] >= x[1] { (x, y) } else { (y, x) };
            (first[0] + second[1]) / 2.0
        };
        let site = |fr: [f64; 3]| -> Vec3 { o + u.scale(fr[0]) + v.scale(fr[1]) + nrm.scale(fr[2]) };
        // Classes: oxygen, hydrogens, lone pairs, optionally bond centres.
        let mut class_pos: Vec<Vec<Vec3>> = vec![vec![pos[0]], vec![pos[1], pos[2]], vec![site(lp), site([lp[1], lp[0], -lp[2]])]];
        let mut sites_text = format!("site 2 {:e} {:e} {:e}\nsite 2 {:e} {:e} {:e}\n", lp[0], lp[1], lp[2], lp[1], lp[0], -lp[2]);
        if want_bond {
            class_pos.push(vec![site([bc_a, 0.0, 0.0]), site([0.0, bc_a, 0.0])]);
            sites_text += &format!("site 3 {:e} 0e0 0e0\nsite 3 0e0 {:e} 0e0\n", bc_a, bc_a);
        }
        println!("  lone pairs at frame ({:.4}, {:.4}, +-{:.4}) bohr{}", lp[0], lp[1], lp[2], if want_bond { format!("; bond centres at {:.4} along each bond", bc_a) } else { String::new() });
        let k = class_pos.len();
        let mults: Vec<f64> = class_pos.iter().map(|c| c.len() as f64).collect();
        let fit_at = |x: &[f64]| -> (Vec<f64>, f64) {
            let sig: Vec<f64> = x.iter().map(|v| v.exp()).collect();
            let rows: Vec<Vec<f64>> = pts
                .iter()
                .map(|q| {
                    (0..k)
                        .map(|i| {
                            class_pos[i]
                                .iter()
                                .map(|s| {
                                    let r = (*q - *s).norm();
                                    phys::induction::erf(r / sig[i]) / r
                                })
                                .sum()
                        })
                        .collect()
                })
                .collect();
            let m = k + 1;
            let mut a = vec![vec![0.0; m]; m];
            let mut b = vec![0.0; m];
            for (row, val) in rows.iter().zip(&esp) {
                for i in 0..k {
                    b[i] += row[i] * val;
                    for j in 0..k {
                        a[i][j] += row[i] * row[j];
                    }
                }
            }
            for i in 0..k {
                a[i][k] = mults[i];
                a[k][i] = mults[i];
            }
            for cc in 0..m {
                let piv = (cc..m).max_by(|&u2, &w| a[u2][cc].abs().partial_cmp(&a[w][cc].abs()).unwrap()).unwrap();
                a.swap(cc, piv);
                b.swap(cc, piv);
                for r in cc + 1..m {
                    let f = a[r][cc] / a[cc][cc];
                    for kk in cc..m {
                        a[r][kk] -= f * a[cc][kk];
                    }
                    b[r] -= f * b[cc];
                }
            }
            let mut sol = vec![0.0; m];
            for cc in (0..m).rev() {
                sol[cc] = (b[cc] - (cc + 1..m).map(|kk| a[cc][kk] * sol[kk]).sum::<f64>()) / a[cc][cc];
            }
            let q: Vec<f64> = sol[..k].to_vec();
            let rms = (rows.iter().zip(&esp).map(|(row, val)| (row.iter().zip(&q).map(|(r, qq)| r * qq).sum::<f64>() - val).powi(2)).sum::<f64>() / esp.len() as f64).sqrt();
            (q.clone(), rms + 1e-6 * q.iter().map(|v| v * v).sum::<f64>())
        };
        // Nelder-Mead over the widths.
        let nn = k;
        let mut simplex: Vec<Vec<f64>> = vec![vec![(0.8f64).ln(); nn]];
        for i in 0..nn {
            let mut vtx = simplex[0].clone();
            vtx[i] += 0.4;
            simplex.push(vtx);
        }
        let mut vals: Vec<f64> = simplex.iter().map(|vtx| fit_at(vtx).1).collect();
        for _ in 0..500 {
            let mut idx: Vec<usize> = (0..=nn).collect();
            idx.sort_by(|&a, &b| vals[a].partial_cmp(&vals[b]).unwrap());
            simplex = idx.iter().map(|&i| simplex[i].clone()).collect();
            vals = idx.iter().map(|&i| vals[i]).collect();
            let centroid: Vec<f64> = (0..nn).map(|kk| simplex[..nn].iter().map(|vtx| vtx[kk]).sum::<f64>() / nn as f64).collect();
            let worst = simplex[nn].clone();
            let along = |t: f64| -> Vec<f64> { (0..nn).map(|kk| centroid[kk] + t * (worst[kk] - centroid[kk])).collect() };
            let refl = along(-1.0);
            let fr = fit_at(&refl).1;
            if fr < vals[0] {
                let ex = along(-2.0);
                let fe = fit_at(&ex).1;
                if fe < fr {
                    simplex[nn] = ex;
                    vals[nn] = fe;
                } else {
                    simplex[nn] = refl;
                    vals[nn] = fr;
                }
            } else if fr < vals[nn - 1] {
                simplex[nn] = refl;
                vals[nn] = fr;
            } else {
                let con = along(0.5);
                let fc = fit_at(&con).1;
                if fc < vals[nn] {
                    simplex[nn] = con;
                    vals[nn] = fc;
                } else {
                    for i in 1..=nn {
                        for kk in 0..nn {
                            simplex[i][kk] = simplex[0][kk] + 0.5 * (simplex[i][kk] - simplex[0][kk]);
                        }
                        vals[i] = fit_at(&simplex[i]).1;
                    }
                }
            }
        }
        let (q, rms) = fit_at(&simplex[0]);
        let sig: Vec<f64> = simplex[0].iter().map(|v| v.exp()).collect();
        println!("  charges {q:?}, widths {sig:?} bohr, potential RMS {rms:.3e} hartree/e (point charges on the atoms and one bisector site: {:.3e})", best.0);
        let text = format!("charges {}\nsigma {}\n{sites_text}", q.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "), sig.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "));
        std::fs::write(format!("esp-{name}-sites{}.txt", if scale != 1.0 { "-mp2" } else { "" }), text).expect("written");
        return;
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
        std::fs::write(format!("esp-{name}-gauss{}.txt", if scale != 1.0 { "-mp2" } else { "" }), text).expect("written");
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
