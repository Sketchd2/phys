//! Checks and measurements of the Ewald sum on a liquid snapshot
//! (`docs/PLAY.md` E8d).
//!
//! ```sh
//! phys-ewald law.txt grow-water.state snapshot.snap [s' bohr] [check]
//! ```
//!
//! Builds the liquid from the snapshot (the molecules' centres and orientations
//! fitted to their atoms), evaluates it with the law truncated as `phys-bulk`
//! does and with the Ewald sum at the splitting width(s) given, and prints the
//! energy per molecule and the pressure. With `check`, also: the energy's
//! independence of the splitting width, the force and the torque against finite
//! differences of the (re-minimised) energy, and the virial against minus the
//! derivative of the energy under rigid scaling of the box.

use crate::electrons::grow::{equivalent_atoms, Resume};
use crate::electrons::molecule::Molecule;
use crate::ewald::Ewald;
use crate::liquid::{EwaldLaw, Kind, Liquid, Polarisable, SiteLaw, SiteSite, AMU, BAR_PER_HARTREE_PER_BOHR3, K_B};
use crate::math::{Quat, Vec3};

const KCAL: f64 = 627.509474;

/// The rotation taking the principal-frame vectors `(b1, b2)` to `(a1, a2)`, as a quaternion.
fn rotation_between(b1: Vec3, b2: Vec3, a1: Vec3, a2: Vec3) -> Quat {
    let frame = |u1: Vec3, u2: Vec3| -> [Vec3; 3] {
        let e1 = u1.unit();
        let e3 = u1.cross(u2).unit();
        let e2 = e3.cross(e1);
        [e1, e2, e3]
    };
    let (a, b) = (frame(a1, a2), frame(b1, b2));
    // R = sum_k a_k b_k^T
    let mut r = [[0.0f64; 3]; 3];
    for k in 0..3 {
        let (ak, bk) = ([a[k].x, a[k].y, a[k].z], [b[k].x, b[k].y, b[k].z]);
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] += ak[i] * bk[j];
            }
        }
    }
    let tr = r[0][0] + r[1][1] + r[2][2];
    let (w, x, y, z);
    if tr > 0.0 {
        let s = (tr + 1.0).sqrt() * 2.0;
        w = 0.25 * s;
        x = (r[2][1] - r[1][2]) / s;
        y = (r[0][2] - r[2][0]) / s;
        z = (r[1][0] - r[0][1]) / s;
    } else if r[0][0] > r[1][1] && r[0][0] > r[2][2] {
        let s = (1.0 + r[0][0] - r[1][1] - r[2][2]).sqrt() * 2.0;
        w = (r[2][1] - r[1][2]) / s;
        x = 0.25 * s;
        y = (r[0][1] + r[1][0]) / s;
        z = (r[0][2] + r[2][0]) / s;
    } else if r[1][1] > r[2][2] {
        let s = (1.0 + r[1][1] - r[0][0] - r[2][2]).sqrt() * 2.0;
        w = (r[0][2] - r[2][0]) / s;
        x = (r[0][1] + r[1][0]) / s;
        y = 0.25 * s;
        z = (r[1][2] + r[2][1]) / s;
    } else {
        let s = (1.0 + r[2][2] - r[0][0] - r[1][1]).sqrt() * 2.0;
        w = (r[1][0] - r[0][1]) / s;
        x = (r[0][2] + r[2][0]) / s;
        y = (r[1][2] + r[2][1]) / s;
        z = 0.25 * s;
    }
    Quat { w, v: Vec3 { x, y, z } }.unit()
}

pub fn run(args: Vec<String>) {
    let law_path = args.first().expect("a law file");
    let state_path = args.get(1).expect("a growth state");
    let snap_path = args.get(2).expect("a snapshot");
    let widths: Vec<f64> = args.get(3).map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| vec![4.2]);
    let check = args.iter().any(|a| a == "check");
    let acc: f64 = std::env::var("PHYS_EWALD_ACC").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-7);
    let temperature = 298.0;

    let law_text = std::fs::read_to_string(law_path).unwrap_or_else(|_| panic!("no {law_path}"));
    let law = SiteSite::from_text(&law_text).expect("a readable law");
    let bisector = SiteSite::bisector_from_text(&law_text);
    let alpha = SiteSite::alpha_from_text(&law_text, law.charge.len());
    let polarised = Polarisable { law: &law, alpha };
    let state = Resume::from_text(&std::fs::read_to_string(state_path).expect("a state")).expect("readable");
    let z = vec![8u32, 1, 1];
    let mono = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let atom_types = equivalent_atoms(&mono);
    let frame = crate::liquid::extra_sites_from_text(&law_text);
    let kind = Kind::of_molecule_sites(&z, &state.positions, &atom_types, bisector, &frame);

    let snap = std::fs::read_to_string(snap_path).unwrap_or_else(|_| panic!("no {snap_path}"));
    let mut cell = [0.0; 3];
    let mut mols: Vec<Vec<Vec3>> = Vec::new();
    for line in snap.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w.first().copied() {
            Some("box") => (0..3).for_each(|k| cell[k] = w[1 + k].parse().unwrap()),
            Some("mol") => {
                let v: Vec<f64> = w[2..].iter().map(|x| x.parse().unwrap()).collect();
                mols.push(v.chunks(3).map(|c| Vec3 { x: c[0], y: c[1], z: c[2] }).collect());
            }
            _ => {}
        }
    }
    let n = mols.len();
    let volume = cell[0] * cell[1] * cell[2];
    let r_cut = 16.0f64.min(0.49 * cell[0]);
    let r_on = r_cut - 2.0;
    let molecule_g = kind.mass / AMU * 1.66053906660e-24;
    let _ = molecule_g;
    let mut liquid = Liquid::slab(kind, n, n as f64 / volume, cell[0], cell[2], temperature, 1, r_on, r_cut);
    liquid.cell = cell;
    let s = liquid.kind.sites.clone();
    for (i, m) in mols.iter().enumerate() {
        let q = rotation_between(s[1] - s[0], s[2] - s[0], m[1] - m[0], m[2] - m[0]);
        liquid.orientation[i] = q;
        liquid.com[i] = m[0] - q.rotate(s[0]);
    }
    let kt_per_volume = n as f64 * K_B * temperature / volume * BAR_PER_HARTREE_PER_BOHR3;
    println!("{n} molecules, box {:.3} bohr, cut-off {r_on}-{r_cut}, law {law_path}", cell[0]);

    let report = |label: &str, l: &Liquid, law: &dyn SiteLaw| -> (f64, f64) {
        let t = std::time::Instant::now();
        let f = l.forces(law);
        let p = kt_per_volume + f.virial / (3.0 * (l.cell[0] * l.cell[1] * l.cell[2])) * BAR_PER_HARTREE_PER_BOHR3;
        println!("  {label:28} U/N {:9.4} kcal/mol   P {:9.1} bar   ({:.1} s)", f.energy / n as f64 * KCAL, p, t.elapsed().as_secs_f64());
        (f.energy, p)
    };
    report("truncated (as phys-bulk)", &liquid, &polarised);
    let mut energies = Vec::new();
    for &sp in &widths {
        let ew = EwaldLaw { inner: &polarised, ewald: Ewald::new(sp, cell, acc) };
        let (e, _) = report(&format!("Ewald s' = {sp} ({} k)", ew.ewald.nk()), &liquid, &ew);
        energies.push(e);
    }
    if energies.len() > 1 {
        let spread = energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - energies.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("  spread of the energy over the widths: {:.3e} hartree ({:.3e} kcal/mol per molecule)", spread, spread / n as f64 * KCAL);
    }
    // `profile`: where a force evaluation spends its time, warm-started as in a run.
    if args.iter().any(|a| a == "profile") {
        let sp = widths[0];
        let plain = EwaldLaw { inner: &law, ewald: Ewald::new(sp, cell, acc) };
        let full = EwaldLaw { inner: &polarised, ewald: Ewald::new(sp, cell, acc) };
        let time = |label: &str, f: &dyn Fn(Option<&[Vec<Vec3>]>) -> Vec<Vec<Vec3>>| {
            let mut warm = f(None);
            let t = std::time::Instant::now();
            let reps = 6;
            for _ in 0..reps {
                warm = f(Some(&warm));
            }
            println!("  {label:46} {:7.1} ms", t.elapsed().as_secs_f64() * 1000.0 / reps as f64);
        };
        println!("profile at s' = {sp}, accuracy {acc} ({} k), per force evaluation:", full.ewald.nk());
        time("pair loop only (no induction, no Ewald)", &|_| {
            liquid.forces(&law);
            Vec::new()
        });
        time("pair loop + charge reciprocal sum (Ewald, no induction)", &|_| {
            liquid.forces(&plain);
            Vec::new()
        });
        time("truncated with induction (as phys-bulk)", &|w| liquid.forces_from(&polarised, w).dipoles);
        time("Ewald with induction (the production run)", &|w| liquid.forces_from(&full, w).dipoles);
        return;
    }
    if !check {
        return;
    }

    // Forces against finite differences of the energy.
    let sp = widths[0];
    let ew = EwaldLaw { inner: &polarised, ewald: Ewald::new(sp, cell, acc) };
    let f0 = liquid.forces(&ew);
    println!("force and torque of molecules against central differences of the energy (s' = {sp}):");
    for &i in &[0usize, 17, 101] {
        for axis in 0..3 {
            let delta = 2e-3;
            let e_of = |dx: f64| -> f64 {
                let mut l = liquid.clone();
                let mut d = Vec3::ZERO;
                match axis {
                    0 => d.x = dx,
                    1 => d.y = dx,
                    _ => d.z = dx,
                }
                l.com[i] += d;
                l.forces(&ew).energy
            };
            let fd = -(e_of(delta) - e_of(-delta)) / (2.0 * delta);
            let fa = [f0.force[i].x, f0.force[i].y, f0.force[i].z][axis];
            // Rotation about the body axis `axis`.
            let e_rot = |th: f64| -> f64 {
                let mut l = liquid.clone();
                let e = [Vec3 { x: 1.0, y: 0.0, z: 0.0 }, Vec3 { x: 0.0, y: 1.0, z: 0.0 }, Vec3 { x: 0.0, y: 0.0, z: 1.0 }][axis];
                l.orientation[i] = l.orientation[i].then(Quat::from_axis_angle(e, th)).unit();
                l.forces(&ew).energy
            };
            let td = -(e_rot(delta) - e_rot(-delta)) / (2.0 * delta);
            let tb = liquid.orientation[i].conjugate().rotate(f0.torque[i]);
            let ta = [tb.x, tb.y, tb.z][axis];
            println!("  molecule {i:3} axis {axis}: force {fa:+.6e} vs {fd:+.6e} (diff {:+.1e}); torque {ta:+.6e} vs {td:+.6e} (diff {:+.1e})", fa - fd, ta - td);
        }
    }

    // The virial against minus the energy's derivative under rigid scaling.
    println!("virial against -dU/dln(lambda) under rigid scaling of the box and the centres:");
    for (label, with_ewald) in [("truncated", false), ("Ewald", true)] {
        let e_at = |lam: f64| -> f64 {
            let mut l = liquid.clone();
            l.cell = [cell[0] * lam, cell[1] * lam, cell[2] * lam];
            for c in l.com.iter_mut() {
                *c = c.scale(lam);
            }
            if with_ewald {
                let e2 = EwaldLaw { inner: &polarised, ewald: Ewald::new(sp, l.cell, acc) };
                l.forces(&e2).energy
            } else {
                l.forces(&polarised).energy
            }
        };
        let eps = 2e-4;
        let w_fd = -(e_at(1.0 + eps) - e_at(1.0 - eps)) / ((1.0 + eps).ln() - (1.0 - eps).ln());
        let w = if with_ewald { f0.virial } else { liquid.forces(&polarised).virial };
        println!("  {label:10} virial {w:+.6e}  finite difference {w_fd:+.6e}  (diff {:+.2e} hartree = {:+.1} bar)", w - w_fd, (w - w_fd) / (3.0 * volume) * BAR_PER_HARTREE_PER_BOHR3);
    }
}
