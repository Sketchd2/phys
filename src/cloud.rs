//! A molecule's charge clouds from its own electrons, for any molecule —
//! `docs/SOLUTIONS.md` S1, the law-building recipe generalised from water.
//!
//! The water law's electrostatics (`PLAY.md` E8b (l)) are cores and clouds: each nucleus less its
//! innermost electrons a point, the valence electrons Gaussian clouds on the atoms and on the
//! lone pair, with populations and widths found from the density. `phys-esp --cloud` did that for
//! water alone (a cloud on the oxygen, one on each hydrogen, one on the bisector). This does it for
//! a molecule it knows nothing about, with nothing typed in:
//!
//! * **cores** are the nuclei less their innermost electrons (one 1s orbital in the second row,
//!   1s2s2p in the third), counted from the molecule's own orbitals;
//! * **a cloud on every atom**, one population and one width for each symmetry class of atom;
//! * **a cloud on each atom's lone pairs**, where the molecule's own Boys-localised orbitals put
//!   them: an orbital whose centroid is near one nucleus and far (>1.6 bohr) from every other is a
//!   lone pair of that atom; the atom's lone pairs share one cloud, on the line from the nucleus to
//!   the mean of their centroids, at a distance that is fitted. For water that line is the
//!   bisector, for ammonia the three-fold axis, for methanol's oxygen the bisector of C-O-H's
//!   outside, and none of those is said anywhere;
//! * the **populations** are linear and are found under neutrality and non-negativity; the
//!   widths and the lone-pair distances are searched, bounded, by Nelder-Mead.
//!
//! The target is the electrostatic potential of the Hartree-Fock density on shells around the
//! molecule, with the electrons' part scaled so the dipole is the MP2 dipole (finite-field), as
//! `phys-esp` does. **This is a starting cloud and not the law's:** fitted to a potential outside
//! the molecule it cannot see penetration, and `PLAY.md` E8b (k)-(l) found it overshoots by
//! several kcal/mol at hydrogen-bond contact; the clouds are refitted to the pairs' exact
//! first-order electrostatic energies afterwards, starting from these (which are sound: positive,
//! bounded, and where the orbitals put the lone pairs, unlike the unconstrained potential fit,
//! whose output for water was degenerate).

use crate::electrons::functional::Functional;
use crate::electrons::grow::{equivalent_atoms, Resume};
use crate::electrons::hf::{frozen_core, hartree_fock, hartree_fock_in_field, mp2};
use crate::electrons::integrals::{dipole, one_electron};
use crate::electrons::localise::boys;
use crate::electrons::molecule::Molecule;
use crate::math::Vec3;

const A0_ANGSTROM: f64 = 0.529177210903;
const DEBYE_PER_AU: f64 = 1.0 / 0.393430307;

/// The cloud model found for a molecule.
#[derive(Debug, Clone)]
pub struct Cloud {
    pub name: String,
    /// The class of each atom (`equivalent_atoms`).
    pub classes: Vec<usize>,
    /// The valence charge of each atom class's core (nuclear charge less core electrons).
    pub core: Vec<f64>,
    /// One lone-pair cloud for each atom class that has any: the class, the number of lone pairs
    /// it stands for, its population (electrons), width (bohr) and distance from the nucleus
    /// along the line to its lone pairs.
    pub lone: Vec<LonePair>,
    /// Per atom class: the atom cloud's population and width.
    pub atom_pop: Vec<f64>,
    pub atom_width: Vec<f64>,
    /// The potential's RMS misfit, hartree/e, and what atom-centred point charges alone give.
    pub rms: f64,
    pub rms_points: f64,
    pub valence: f64,
    pub dipole_au: f64,
    /// The unit vector along Hartree-Fock's dipole, in the molecule's coordinates. Pair energies are
    /// invariant under reversing every charge in both molecules, so a fit to them alone can land on the
    /// mirror-charge molecule (oxygen positive, hydrogens negative); the dipole's sign pins it.
    pub dipole_dir: [f64; 3],
}

#[derive(Debug, Clone)]
pub struct LonePair {
    pub class: usize,
    /// The atoms of that class that carry one, with the unit vector from each nucleus toward its
    /// lone pairs' centre.
    pub directions: Vec<(usize, [f64; 3])>,
    pub pairs: usize,
    pub population: f64,
    pub width: f64,
    pub distance: f64,
    /// Where the localised orbitals put the centre, as a distance from the nucleus.
    pub centroid_distance: f64,
}

fn solve_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let m = b.len();
    for c in 0..m {
        let piv = (c..m).max_by(|&x, &y| a[x][c].abs().total_cmp(&a[y][c].abs())).unwrap();
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
    x
}

/// Nelder-Mead minimisation of `f` from `x0` with initial steps `step`.
pub fn nelder_mead(f: &dyn Fn(&[f64]) -> f64, x0: &[f64], step: &[f64], iterations: usize) -> Vec<f64> {
    let n = x0.len();
    let mut simplex: Vec<Vec<f64>> = vec![x0.to_vec()];
    for i in 0..n {
        let mut v = x0.to_vec();
        v[i] += step[i];
        simplex.push(v);
    }
    let mut vals: Vec<f64> = simplex.iter().map(|v| f(v)).collect();
    for _ in 0..iterations {
        let mut idx: Vec<usize> = (0..=n).collect();
        idx.sort_by(|&a, &b| vals[a].total_cmp(&vals[b]));
        simplex = idx.iter().map(|&i| simplex[i].clone()).collect();
        vals = idx.iter().map(|&i| vals[i]).collect();
        let centroid: Vec<f64> = (0..n).map(|k| simplex[..n].iter().map(|v| v[k]).sum::<f64>() / n as f64).collect();
        let worst = simplex[n].clone();
        let along = |t: f64| -> Vec<f64> { (0..n).map(|k| centroid[k] + t * (worst[k] - centroid[k])).collect() };
        let refl = along(-1.0);
        let fr = f(&refl);
        if fr < vals[0] {
            let ex = along(-2.0);
            let fe = f(&ex);
            if fe < fr {
                simplex[n] = ex;
                vals[n] = fe;
            } else {
                simplex[n] = refl;
                vals[n] = fr;
            }
        } else if fr < vals[n - 1] {
            simplex[n] = refl;
            vals[n] = fr;
        } else {
            let con = along(0.5);
            let fc = f(&con);
            if fc < vals[n] {
                simplex[n] = con;
                vals[n] = fc;
            } else {
                for i in 1..=n {
                    for k in 0..n {
                        simplex[i][k] = simplex[0][k] + 0.5 * (simplex[i][k] - simplex[0][k]);
                    }
                    vals[i] = f(&simplex[i]);
                }
            }
        }
    }
    let best = (0..=n).min_by(|&a, &b| vals[a].total_cmp(&vals[b])).unwrap();
    simplex[best].clone()
}

/// The cores' valence charge: the nuclear charge less the electrons of the innermost shells (1s in
/// the second row, 1s2s2p in the third), the same counts the frozen-core rule uses.
fn core_electrons(z: u32) -> f64 {
    match z {
        0..=2 => 0.0,
        3..=10 => 2.0,
        11..=18 => 10.0,
        other => panic!("no core count for Z = {other}"),
    }
}

/// Fit the cloud model to a molecule named as `phys-grow` knows it. `mp2_dipole` scales the
/// electrons' potential to the MP2 dipole (about eight per cent less than Hartree-Fock's for
/// water); the state `grow-<name>.state` must be in the working directory.
pub fn fit(name: &str, mp2_dipole: bool, log: &mut dyn FnMut(String)) -> Cloud {
    let z = crate::pairs::atoms_of(name);
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).unwrap_or_else(|_| panic!("no grow-{name}.state"))).expect("readable");
    let mol = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let classes = equivalent_atoms(&mol);
    let kclasses = classes.iter().max().map(|m| m + 1).unwrap_or(1);
    let p = mol.problem(Functional::Pbe);
    let hf = hartree_fock(&p, 100, 1e-10);
    assert!(hf.converged, "the Hartree-Fock of {name} did not converge");
    let n = p.basis.size;
    let mut d = crate::electrons::linalg::Matrix::zeros(n);
    for a in 0..n {
        for b in 0..n {
            let s: f64 = (0..hf.occupied).map(|i| hf.c[a * hf.m + i] * hf.c[b * hf.m + i]).sum();
            d.set(a, b, 2.0 * s);
        }
    }
    let pos: Vec<Vec3> = state.positions.iter().map(|r| Vec3 { x: r[0], y: r[1], z: r[2] }).collect();
    let radius = |zz: u32| crate::chem::elements::Element(zz as u8).vdw_radius().unwrap_or(1.5e-10) / (A0_ANGSTROM * 1e-10);
    // ---- the potential's points
    let mut pts: Vec<Vec3> = Vec::new();
    for scale in [1.0, 1.15, 1.3, 1.5, 1.8, 2.2] {
        for (i, c) in pos.iter().enumerate() {
            let r = radius(z[i]) * scale;
            let m = 400;
            for k in 0..m {
                let t = (k as f64 + 0.5) / m as f64;
                let phi = std::f64::consts::PI * (1.0 + 5f64.sqrt()) * k as f64;
                let (ct, st) = (1.0 - 2.0 * t, (1.0 - (1.0 - 2.0 * t).powi(2)).sqrt());
                let q = *c + Vec3 { x: st * phi.cos(), y: st * phi.sin(), z: ct }.scale(r);
                if pos.iter().enumerate().all(|(j, o)| (q - *o).norm() >= radius(z[j])) {
                    pts.push(q);
                }
            }
        }
    }
    // ---- the dipole, from the density and (if asked) from MP2's field derivative
    let dm = dipole(&p.basis, [0.0; 3]);
    let (mut mu_hf, mut mu_nuc) = ([0.0f64; 3], [0.0f64; 3]);
    for k in 0..3 {
        mu_nuc[k] = pos.iter().zip(&z).map(|(o, zz)| *zz as f64 * [o.x, o.y, o.z][k]).sum();
        let mut el = 0.0;
        for a in 0..n {
            for b in 0..n {
                el += d.get(a, b) * dm[k].get(a, b);
            }
        }
        mu_hf[k] = mu_nuc[k] - el;
    }
    let norm3 = |v: &[f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let mu_target = if mp2_dipole {
        let all: Vec<usize> = (0..z.len()).collect();
        let fc = frozen_core(&z, &all);
        let step = 0.002;
        let (mut mu_ff, mut mu_ff_hf) = ([0.0f64; 3], [0.0f64; 3]);
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
        let sign = if (0..3).map(|k| mu_ff_hf[k] * mu_hf[k]).sum::<f64>() >= 0.0 { 1.0 } else { -1.0 };
        norm3(&[sign * mu_ff[0], sign * mu_ff[1], sign * mu_ff[2]])
    } else {
        norm3(&mu_hf)
    };
    log(format!("{name}: dipole Hartree-Fock {:.4} au ({:.3} D), target {:.4} au ({:.3} D)", norm3(&mu_hf), norm3(&mu_hf) * DEBYE_PER_AU, mu_target, mu_target * DEBYE_PER_AU));
    // The potential is Hartree-Fock's, unscaled. Scaling the electrons' part to the MP2 dipole (as
    // `phys-esp` did) changes their monopole too, so the target is no longer neutral, which the fit
    // assumes: that mismatch is what made the old unconstrained fit degenerate. The MP2 dipole is
    // reported for the refit to hold (`phys-diag --dipole`), where it belongs.
    let s_e = 1.0;
    let cores: Vec<f64> = z.iter().map(|&zz| zz as f64 - core_electrons(zz)).collect();
    let valence: f64 = cores.iter().sum();
    // The potential of the electrons (scaled) plus all the nuclei, less the cores: what the clouds
    // (negative charge, `valence` electrons in all) must reproduce is minus this.
    let target: Vec<f64> = pts
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
            let total = nuc + s_e * e;
            total - pos.iter().zip(&cores).map(|(o, c)| c / (*q - *o).norm()).sum::<f64>()
        })
        .collect();
    log(format!("  {} potential points; the target's RMS {:.3e}, its mean {:.3e} hartree/e", pts.len(), (target.iter().map(|v| v * v).sum::<f64>() / target.len() as f64).sqrt(), target.iter().sum::<f64>() / target.len() as f64));
    // ---- lone pairs, from the localised valence orbitals
    let ncore: usize = z.iter().map(|&zz| match zz { 3..=10 => 1, 11..=18 => 5, _ => 0 }).sum();
    let loc = boys(&hf.c, hf.m, n, ncore..hf.occupied, &dm);
    let mut lone_by_atom: Vec<Vec<Vec3>> = vec![Vec::new(); z.len()];
    let mut bonds = 0;
    for cen in &loc.centroids {
        let q = Vec3 { x: cen[0], y: cen[1], z: cen[2] };
        let mut dist: Vec<(f64, usize)> = pos.iter().enumerate().map(|(i, o)| ((q - *o).norm(), i)).collect();
        dist.sort_by(|a, b| a.0.total_cmp(&b.0));
        if dist.len() > 1 && dist[1].0 < 1.6 {
            bonds += 1;
        } else {
            lone_by_atom[dist[0].1].push(q);
        }
    }
    let lone_total: usize = lone_by_atom.iter().map(|v| v.len()).sum();
    log(format!("  localised in {} sweeps: {} bond orbitals, {} lone pairs ({} core orbitals frozen out of the localisation)", loc.sweeps, bonds, lone_total, ncore));
    // One lone-pair cloud per atom class that has lone pairs.
    struct LpClass {
        class: usize,
        atoms: Vec<(usize, Vec3, usize)>, // (atom, unit direction, number of pairs)
        centroid_distance: f64,
    }
    let mut lp_classes: Vec<LpClass> = Vec::new();
    for a in 0..z.len() {
        if lone_by_atom[a].is_empty() {
            continue;
        }
        let m = lone_by_atom[a].iter().fold(Vec3::ZERO, |acc, q| acc + *q).scale(1.0 / lone_by_atom[a].len() as f64) - pos[a];
        let (dir, dist) = (m.unit(), m.norm());
        match lp_classes.iter_mut().find(|c| c.class == classes[a]) {
            Some(c) => c.atoms.push((a, dir, lone_by_atom[a].len())),
            None => lp_classes.push(LpClass { class: classes[a], atoms: vec![(a, dir, lone_by_atom[a].len())], centroid_distance: dist }),
        }
    }
    // ---- the fit: populations linear, widths and lone-pair distances searched
    let nl = lp_classes.len();
    let ncl = kclasses + nl;
    // Parameters: ln width of each atom class, then of each lone-pair class.
    let sigma_of = |x: f64| x.exp().clamp(0.25, 3.5);
    let sites_of = |x: &[f64]| -> (Vec<Vec<Vec3>>, Vec<f64>) {
        let mut sites: Vec<Vec<Vec3>> = Vec::new();
        let mut widths = Vec::new();
        for c in 0..kclasses {
            sites.push((0..z.len()).filter(|&a| classes[a] == c).map(|a| pos[a]).collect());
            widths.push(sigma_of(x[c]));
        }
        for (l, lc) in lp_classes.iter().enumerate() {
            // Where the orbitals put the lone pairs' centre: derived, not searched (a potential fit
            // outside the molecule cannot tell 0.03 bohr from 0.3, and ran to the bound).
            let dist = lc.centroid_distance;
            sites.push(lc.atoms.iter().map(|(a, dir, _)| pos[*a] + dir.scale(dist)).collect());
            widths.push(sigma_of(x[kclasses + l]));
        }
        (sites, widths)
    };
    let lp_pairs: Vec<f64> = lp_classes.iter().map(|lc| lc.atoms.iter().map(|a| a.2 as f64).sum::<f64>() / lc.atoms.len() as f64).collect();
    let solve_pops = |x: &[f64]| -> (Vec<f64>, f64, Vec<Vec<Vec3>>, Vec<f64>) {
        let (sites, widths) = sites_of(x);
        let rows: Vec<Vec<f64>> = pts.iter().map(|q| (0..ncl).map(|i| sites[i].iter().map(|s| -crate::induction::erf((*q - *s).norm() / widths[i]) / (*q - *s).norm()).sum()).collect()).collect();
        let mult: Vec<f64> = sites.iter().map(|s| s.len() as f64).collect();
        // The lone pairs' populations start from the orbitals' (two electrons a pair) and are held
        // there by a weak prior (5 % of the data's curvature in each): a potential outside the
        // molecule cannot tell a lone pair's electrons from the atom's own, and unconstrained, the
        // lone-pair cloud empties. Populations are non-negative and the total is the valence.
        let prior: Vec<f64> = (0..ncl).map(|i| if i >= kclasses { 2.0 * lp_pairs[i - kclasses] } else { 0.0 }).collect();
        let mut active: Vec<usize> = (0..ncl).collect();
        let mut pops = vec![0.0; ncl];
        for _ in 0..ncl {
            let m = active.len() + 1;
            let mut a = vec![vec![0.0; m]; m];
            let mut b = vec![0.0; m];
            for (row, val) in rows.iter().zip(&target) {
                for (ii, &i) in active.iter().enumerate() {
                    b[ii] += row[i] * val;
                    for (jj, &j) in active.iter().enumerate() {
                        a[ii][jj] += row[i] * row[j];
                    }
                }
            }
            for (ii, &i) in active.iter().enumerate() {
                if i >= kclasses {
                    let w = 0.05 * a[ii][ii];
                    a[ii][ii] += w;
                    b[ii] += w * prior[i];
                } else {
                    a[ii][ii] += 1e-9;
                }
                a[ii][active.len()] = mult[i];
                a[active.len()][ii] = mult[i];
            }
            b[active.len()] = valence;
            let sol = solve_linear(a, b);
            let worst = active.iter().enumerate().map(|(ii, &i)| (sol[ii], i)).min_by(|a, b| a.0.total_cmp(&b.0)).unwrap();
            if worst.0 >= 0.0 {
                for (ii, &i) in active.iter().enumerate() {
                    pops[i] = sol[ii];
                }
                break;
            }
            active.retain(|&i| i != worst.1);
            pops = vec![0.0; ncl];
        }
        let rms = (rows.iter().zip(&target).map(|(row, val)| (row.iter().zip(&pops).map(|(r, nn)| r * nn).sum::<f64>() - val).powi(2)).sum::<f64>() / target.len() as f64).sqrt();
        // A lone-pair cloud holds electrons for its pairs: no fewer than none, and keep it sane.
        (pops, rms, sites, widths)
    };
    let objective = |x: &[f64]| -> f64 {
        let (pops, rms, _, _) = solve_pops(x);
        let over = pops.iter().map(|p| (p - 8.0).max(0.0)).sum::<f64>();
        rms + 1e-3 * over
    };
    let mut x0: Vec<f64> = Vec::new();
    for c in 0..kclasses {
        let heavy = (0..z.len()).find(|&a| classes[a] == c).map(|a| z[a] > 1).unwrap_or(true);
        x0.push(if heavy { 0.8f64.ln() } else { 0.6f64.ln() });
    }
    for _ in 0..nl {
        x0.push(1.0f64.ln());
    }
    let step: Vec<f64> = vec![0.4; x0.len()];
    let best = nelder_mead(&objective, &x0, &step, 1500);
    let (pops, rms, sites, widths) = solve_pops(&best);
    // The potential atom-centred point charges alone reach, for the record.
    let rms_points = {
        let rows: Vec<Vec<f64>> = pts.iter().map(|q| (0..kclasses).map(|i| -(0..z.len()).filter(|&a| classes[a] == i).map(|a| 1.0 / (*q - pos[a]).norm()).sum::<f64>()).collect()).collect();
        let mult: Vec<f64> = (0..kclasses).map(|c| (0..z.len()).filter(|&a| classes[a] == c).count() as f64).collect();
        let m = kclasses + 1;
        let mut a = vec![vec![0.0; m]; m];
        let mut b = vec![0.0; m];
        for (row, val) in rows.iter().zip(&target) {
            for i in 0..kclasses {
                b[i] += row[i] * val;
                for j in 0..kclasses {
                    a[i][j] += row[i] * row[j];
                }
            }
        }
        for i in 0..kclasses {
            a[i][kclasses] = mult[i];
            a[kclasses][i] = mult[i];
        }
        b[kclasses] = valence;
        let sol = solve_linear(a, b);
        (rows.iter().zip(&target).map(|(row, val)| (row.iter().zip(&sol).map(|(r, q)| r * q).sum::<f64>() - val).powi(2)).sum::<f64>() / target.len() as f64).sqrt()
    };
    let core: Vec<f64> = (0..kclasses).map(|c| cores[(0..z.len()).find(|&a| classes[a] == c).unwrap()]).collect();
    let lone: Vec<LonePair> = lp_classes
        .iter()
        .enumerate()
        .map(|(l, lc)| LonePair {
            class: lc.class,
            directions: lc.atoms.iter().map(|(a, dir, _)| (*a, [dir.x, dir.y, dir.z])).collect(),
            pairs: lp_pairs[l].round() as usize,
            population: pops[kclasses + l],
            width: widths[kclasses + l],
            distance: lc.centroid_distance,
            centroid_distance: lc.centroid_distance,
        })
        .collect();
    let _ = sites;
    Cloud { name: name.to_string(), classes, core, lone, atom_pop: pops[..kclasses].to_vec(), atom_width: widths[..kclasses].to_vec(), rms, rms_points, valence, dipole_au: mu_target, dipole_dir: if norm3(&mu_hf) > 1e-8 { [mu_hf[0] / norm3(&mu_hf), mu_hf[1] / norm3(&mu_hf), mu_hf[2] / norm3(&mu_hf)] } else { [0.0; 3] } }
}

/// The frame coordinates `(a, b, c)` of `point` in the frame `ExtraSite` uses: `o + a u + b v + c n`
/// with `u`, `v` the unit vectors from atom 0 to atoms 1 and 2 and `n` their unit cross product.
pub fn frame_coordinates(positions: &[Vec3], point: Vec3) -> [f64; 3] {
    let o = positions[0];
    let (u, v) = ((positions[1] - o).unit(), (positions[2] - o).unit());
    let nrm = u.cross(v).unit();
    let det = |x: Vec3, y: Vec3, w: Vec3| x.dot(y.cross(w));
    let r = point - o;
    let d0 = det(u, v, nrm);
    [det(r, v, nrm) / d0, det(u, r, nrm) / d0, det(u, v, r) / d0]
}

impl Cloud {
    /// The law's text for this cloud: the site types are the atom classes (cores), then each
    /// lone-pair class, then each atom class's cloud, in the order `Kind::of_molecule_sites` is
    /// given them; the extra sites are `site` lines in the frame of atoms 0, 1, 2.
    pub fn law_text(&self, positions: &[Vec3]) -> String {
        let k = self.core.len();
        let nl = self.lone.len();
        let mut charges: Vec<f64> = self.core.clone();
        let mut sigma: Vec<f64> = vec![0.0; k];
        for lp in &self.lone {
            charges.push(-lp.population);
            sigma.push(lp.width);
        }
        for c in 0..k {
            charges.push(-self.atom_pop[c]);
            sigma.push(self.atom_width[c]);
        }
        let join = |v: &[f64]| v.iter().map(|x| format!("{x:e}")).collect::<Vec<_>>().join(" ");
        let mut text = format!("charges {}\nsigma {}\n", join(&charges), join(&sigma));
        for (l, lp) in self.lone.iter().enumerate() {
            for (a, dir) in &lp.directions {
                let site = positions[*a] + Vec3 { x: dir[0], y: dir[1], z: dir[2] }.scale(lp.distance);
                let f = frame_coordinates(positions, site);
                text += &format!("site {} {:e} {:e} {:e}\n", k + l, f[0], f[1], f[2]);
            }
        }
        for a in 0..positions.len() {
            let f = frame_coordinates(positions, positions[a]);
            text += &format!("site {} {:e} {:e} {:e}\n", k + nl + self.classes[a], f[0], f[1], f[2]);
        }
        text
    }
}

// ---- the structure file, and the refit to the pairs' exact electrostatic energies ---------------

impl Cloud {
    /// The structure as text, for the next stage to read: `core`, `atom` and `lone` lines.
    pub fn to_text(&self) -> String {
        let mut s = format!("name {}\nclasses {}\nvalence {}\ndipole {:e}\ndipoledir {:e} {:e} {:e}\n", self.name, self.classes.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(" "), self.valence, self.dipole_au, self.dipole_dir[0], self.dipole_dir[1], self.dipole_dir[2]);
        for (c, q) in self.core.iter().enumerate() {
            s += &format!("core {c} {q:e}\n");
        }
        for c in 0..self.atom_pop.len() {
            s += &format!("atom {c} {:e} {:e}\n", self.atom_pop[c], self.atom_width[c]);
        }
        for lp in &self.lone {
            for (a, d) in &lp.directions {
                s += &format!("lone {} {} {} {:e} {:e} {:e} {a} {:e} {:e} {:e}\n", lp.class, lp.pairs, lp.directions.len(), lp.population, lp.width, lp.distance, d[0], d[1], d[2]);
            }
        }
        s
    }

    pub fn from_text(text: &str) -> Option<Cloud> {
        let (mut name, mut classes, mut core, mut valence, mut dipole) = (String::new(), Vec::new(), Vec::new(), 0.0, 0.0);
        let mut dipole_dir = [0.0f64; 3];
        let (mut atom_pop, mut atom_width): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
        let mut lone: Vec<LonePair> = Vec::new();
        for line in text.lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            match w.first().copied() {
                Some("name") => name = w[1].to_string(),
                Some("classes") => classes = w[1..].iter().map(|x| x.parse().ok()).collect::<Option<Vec<usize>>>()?,
                Some("valence") => valence = w[1].parse().ok()?,
                Some("dipole") => dipole = w[1].parse().ok()?,
                Some("dipoledir") => dipole_dir = [w[1].parse().ok()?, w[2].parse().ok()?, w[3].parse().ok()?],
                Some("core") => core.push(w[2].parse().ok()?),
                Some("atom") => {
                    atom_pop.push(w[2].parse().ok()?);
                    atom_width.push(w[3].parse().ok()?);
                }
                Some("lone") if w.len() == 11 => {
                    let class: usize = w[1].parse().ok()?;
                    let dir = (w[7].parse().ok()?, [w[8].parse().ok()?, w[9].parse().ok()?, w[10].parse().ok()?]);
                    match lone.iter_mut().find(|l| l.class == class) {
                        Some(l) => l.directions.push(dir),
                        None => lone.push(LonePair { class, directions: vec![dir], pairs: w[2].parse().ok()?, population: w[4].parse().ok()?, width: w[5].parse().ok()?, distance: w[6].parse().ok()?, centroid_distance: w[6].parse().ok()? }),
                    }
                }
                _ => {}
            }
        }
        Some(Cloud { name, classes, core, lone, atom_pop, atom_width, rms: 0.0, rms_points: 0.0, valence, dipole_au: dipole, dipole_dir })
    }
}

/// One pair with its exact first-order electrostatic energy, kcal/mol.
pub struct EsRow {
    pub index: usize,
    pub a: Vec<(Vec3, usize)>,
    pub b: Vec<(Vec3, usize)>,
    pub reference: f64,
}

/// The pairs of `pairs-<name><tag>.txt` (comma-separated tags) that `es_file` (`k separation E`,
/// hartree) holds an energy for.
pub fn es_rows(name: &str, es_file: &str, tags: &str) -> Vec<EsRow> {
    let es: std::collections::HashMap<usize, f64> = std::fs::read_to_string(es_file)
        .unwrap_or_else(|_| panic!("no {es_file}"))
        .lines()
        .filter_map(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            Some((w.first()?.parse().ok()?, w.get(2)?.parse().ok()?))
        })
        .collect();
    let text: String = tags.split(',').map(|t| std::fs::read_to_string(format!("pairs-{name}{t}.txt")).unwrap_or_default()).collect::<Vec<_>>().join("\n");
    let mut rows = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#') && l.contains('|')) {
        let (head, tail) = line.split_once('|').expect("a |");
        let k: usize = head.split_whitespace().next().and_then(|x| x.parse().ok()).expect("an index");
        let Some(&reference) = es.get(&k) else { continue };
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<(Vec3, usize)> = t.chunks(4).map(|c| (Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }, c[3].parse().unwrap())).collect();
        let n = atoms.len() / 2;
        rows.push(EsRow { index: k, a: atoms[..n].to_vec(), b: atoms[n..].to_vec(), reference: reference * 627.509474 });
    }
    rows
}

/// The cloud's sites as `ExtraSite`s for a molecule whose atoms are `positions`, with the lone-pair
/// distances given and the types numbered cores, lone-pair clouds, atom clouds.
fn extra_sites(cloud: &Cloud, positions: &[Vec3], distances: &[f64]) -> Vec<crate::liquid::ExtraSite> {
    let k = cloud.core.len();
    let nl = cloud.lone.len();
    let mut out = Vec::new();
    for (l, lp) in cloud.lone.iter().enumerate() {
        for (a, dir) in &lp.directions {
            let site = positions[*a] + Vec3 { x: dir[0], y: dir[1], z: dir[2] }.scale(distances[l]);
            out.push(crate::liquid::ExtraSite { ty: k + l, frame: frame_coordinates(positions, site) });
        }
    }
    for a in 0..positions.len() {
        out.push(crate::liquid::ExtraSite { ty: k + nl + cloud.classes[a], frame: frame_coordinates(positions, positions[a]) });
    }
    out
}

/// Refit a cloud's populations, widths and lone-pair distances to the pairs' exact first-order
/// electrostatic energies (`phys-es-gpu`), half the pairs held out in alternation, the contact
/// region weighted as in a liquid, and the dipole held to `dipole_au` if given. Starts from
/// `start`. Returns the refitted cloud.
pub fn fit_es(start: &Cloud, reference_atoms: &[Vec3], rows: &[EsRow], dipole_au: Option<f64>, log: &mut dyn FnMut(String)) -> Cloud {
    use crate::liquid::{pair_energy, with_frame_sites, PairEnergy, SiteSite};
    let k = start.core.len();
    let nl = start.lone.len();
    let ncl = k + nl;
    let mult: Vec<f64> = (0..k).map(|c| start.classes.iter().filter(|&&x| x == c).count() as f64).collect();
    let lp_mult: Vec<f64> = start.lone.iter().map(|l| l.directions.len() as f64).collect();
    let cloud_mult = |i: usize| if i < k { mult[i] } else { lp_mult[i - k] };
    // Parameters: ln widths of the atom clouds then the lone-pair clouds, the lone-pair distances, then
    // the populations of every cloud but the last atom class (its follows from the valence electrons).
    let build = |x: &[f64]| -> (Cloud, SiteSite) {
        let mut c = start.clone();
        for i in 0..k {
            c.atom_width[i] = x[i].exp().clamp(0.15, 4.0);
        }
        for l in 0..nl {
            c.lone[l].width = x[k + l].exp().clamp(0.15, 4.0);
            c.lone[l].distance = x[ncl + l].abs().clamp(0.02, 1.6);
        }
        let mut pops = vec![0.0; ncl];
        let free: Vec<usize> = (0..k - 1).chain(k..ncl).collect();
        for (j, &i) in free.iter().enumerate() {
            pops[i] = x[ncl + nl + j];
        }
        let used: f64 = (0..ncl).filter(|&i| i != k - 1).map(|i| cloud_mult(i) * pops[i]).sum();
        pops[k - 1] = (start.valence - used) / mult[k - 1];
        for i in 0..k {
            c.atom_pop[i] = pops[i];
        }
        for l in 0..nl {
            c.lone[l].population = pops[k + l];
        }
        let types = ncl + k;
        let mut law = SiteSite { charge: Vec::new(), pair: vec![[0.0, 1.0, 0.0, 0.0]; types * types], sigma: Vec::new(), damp: Vec::new(), overlap: Vec::new() };
        law.charge = c.core.clone();
        law.sigma = vec![0.0; k];
        for lp in &c.lone {
            law.charge.push(-lp.population);
            law.sigma.push(lp.width);
        }
        for i in 0..k {
            law.charge.push(-c.atom_pop[i]);
            law.sigma.push(c.atom_width[i]);
        }
        (c, law)
    };
    let energy = |x: &[f64], r: &EsRow| -> f64 {
        let (c, law) = build(x);
        let d: Vec<f64> = c.lone.iter().map(|l| l.distance).collect();
        let place = |m: &[(Vec3, usize)]| {
            let pos: Vec<Vec3> = m.iter().map(|p| p.0).collect();
            with_frame_sites(m, &extra_sites(&c, &pos, &d))
        };
        pair_energy(&law, &PairEnergy { a: place(&r.a), b: place(&r.b), energy: 0.0 }) * 627.509474
    };
    let dipole_of = |x: &[f64]| -> f64 {
        let (c, law) = build(x);
        let d: Vec<f64> = c.lone.iter().map(|l| l.distance).collect();
        let m: Vec<(Vec3, usize)> = reference_atoms.iter().zip(&c.classes).map(|(p, &t)| (*p, t)).collect();
        let m = with_frame_sites(&m, &extra_sites(&c, reference_atoms, &d));
        let mut mu = Vec3::ZERO;
        for (p, t) in &m {
            mu = mu + p.scale(law.charge[*t]);
        }
        // Signed, along Hartree-Fock's dipole: the mirror-charge molecule has the same pair energies.
        mu.dot(Vec3 { x: start.dipole_dir[0], y: start.dipole_dir[1], z: start.dipole_dir[2] })
    };
    let (mut train, mut test): (Vec<&EsRow>, Vec<&EsRow>) = (Vec::new(), Vec::new());
    for (i, r) in rows.iter().enumerate() {
        if i % 2 == 0 {
            train.push(r);
        } else {
            test.push(r);
        }
    }
    let sep = |r: &EsRow| (r.a[0].0 - r.b[0].0).norm() * 0.529177210903;
    let weight = |r: &EsRow| 1.0 / (1.0 + (sep(r) - 2.7).max(0.0).powi(2));
    let cost = |x: &[f64]| -> f64 {
        let base = train.iter().map(|r| weight(r) * (energy(x, r) - r.reference).powi(2)).sum::<f64>() / train.len() as f64;
        let (c, _) = build(x);
        let negative: f64 = c.atom_pop.iter().chain(c.lone.iter().map(|l| &l.population)).map(|p| (-p).max(0.0)).sum();
        let dip = dipole_of(x);
        let dipole_term = match dipole_au {
            Some(t) => 3000.0 * ((dip - t) / t).powi(2),
            None => 3000.0 * (-dip).max(0.0).powi(2),
        };
        base + 1e3 * negative + dipole_term
    };
    let mut x0: Vec<f64> = Vec::new();
    for i in 0..k {
        x0.push(start.atom_width[i].ln());
    }
    for lp in &start.lone {
        x0.push(lp.width.ln());
    }
    for lp in &start.lone {
        x0.push(lp.distance);
    }
    for i in (0..k - 1).chain(k..ncl) {
        x0.push(if i < k { start.atom_pop[i] } else { start.lone[i - k].population });
    }
    let report = |label: &str, x: &[f64], log: &mut dyn FnMut(String)| {
        let stat = |set: &[&EsRow]| -> (f64, f64) {
            let n = set.len().max(1) as f64;
            (set.iter().map(|r| energy(x, r) - r.reference).sum::<f64>() / n, (set.iter().map(|r| (energy(x, r) - r.reference).powi(2)).sum::<f64>() / n).sqrt())
        };
        log(format!("  {label}: train mean/rms {:+.3}/{:.3}, held-out {:+.3}/{:.3} kcal/mol", stat(&train).0, stat(&train).1, stat(&test).0, stat(&test).1));
        for w in [0.0, 2.8, 3.1, 3.4, 4.7, 99.0].windows(2) {
            let sel: Vec<&EsRow> = rows.iter().filter(|r| sep(r) >= w[0] && sep(r) < w[1]).collect();
            if !sel.is_empty() {
                let (m, s) = stat(&sel);
                log(format!("      first-atom separation {:.1}-{:.1} A: {} pairs, mean {m:+.3}, rms {s:.3}", w[0], w[1].min(9.9), sel.len()));
            }
        }
    };
    report("before", &x0, log);
    let step: Vec<f64> = (0..x0.len()).map(|i| if i < ncl { 0.2 } else if i < ncl + nl { 0.1 } else { 0.3 }).collect();
    // The landscape has several minima (a lone-pair cloud sitting on its nucleus is one), so the given
    // start is joined by seeded random ones: widths log-uniform in 0.3-2.2 bohr, lone-pair distances
    // 0.1-0.8 bohr, populations spread over the valence electrons. The best fit on the training pairs wins.
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut rnd = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut best = nelder_mead(&cost, &x0, &step, 3000);
    let mut best_cost = cost(&best);
    let restarts = 24;
    for r in 0..restarts {
        let mut x = x0.clone();
        for i in 0..ncl {
            x[i] = (0.3f64.ln()) + rnd() * (2.2f64 / 0.3).ln();
        }
        for l in 0..nl {
            x[ncl + l] = 0.1 + rnd() * 0.7;
        }
        let free: Vec<usize> = (0..k - 1).chain(k..ncl).collect();
        let weights: Vec<f64> = (0..ncl).map(|_| 0.2 + rnd()).collect();
        let total: f64 = (0..ncl).map(|i| cloud_mult(i) * weights[i]).sum();
        for (j, &i) in free.iter().enumerate() {
            x[ncl + nl + j] = start.valence * weights[i] / total;
        }
        let cand = nelder_mead(&cost, &x, &step, 3000);
        let c = cost(&cand);
        if c < best_cost {
            log(format!("  restart {r}: cost {c:.4} beats {best_cost:.4}"));
            best_cost = c;
            best = cand;
        }
    }
    report("after ", &best, log);
    let (fitted, _) = build(&best);
    log(format!("  dipole {:.4} au ({:.3} D)", dipole_of(&best), dipole_of(&best) * DEBYE_PER_AU));
    fitted
}
