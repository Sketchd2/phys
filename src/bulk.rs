//! A liquid in a periodic box at one temperature and density — `PLAY.md` E8's
//! first end-to-end look at what a fitted pair law does to a bulk liquid.
//!
//! ```sh
//! cargo run --release --bin phys-bulk -- water law-water.txt 298 [molecules] [g/cm3] [ps]
//! ```
//!
//! The molecule at the shape its growth relaxed it to (`grow-<name>.state`),
//! its atoms typed by symmetry class, moved by the law in the file (and its
//! off-atom site, if the law names one). The lattice's molecules are turned so
//! that no two atoms start closer than 2 A, then eased in with small steps and
//! a strong bath before the run proper. Reports the intermolecular energy per
//! molecule, the heat of vaporisation it implies (`-U/N + RT`, the vapour taken
//! as an ideal gas of molecules at rest in their potential), and the pressure
//! from the molecular virial, whose sign says which way the density wants to
//! move. Long-range electrostatics are cut off with the pair switch, not
//! summed: that is a known error of this prototype, not part of the law.
//!
//! It also accumulates the radial distribution function of the molecules' first
//! atoms (the oxygens, for water) over the production run, writes it to
//! `bulk-<name>-<T>.gOO.txt` (r in angstrom, g), and reports its first peak,
//! first minimum, second peak and the coordination number out to the minimum:
//! the quantities that tell an over-structured liquid (a tall first peak and a
//! deep minimum) from an under-structured one. Measured for water at 298 K:
//! first peak 2.57-2.75 at 2.80 A, first minimum 0.84, coordination about 4.3.

use crate::electrons::grow::{equivalent_atoms, Resume};
use crate::electrons::molecule::Molecule;
use crate::ewald::Ewald;
use crate::liquid::{EwaldLaw, Kind, Liquid, Polarisable, SiteLaw, SiteSite, AMU, BAR_PER_HARTREE_PER_BOHR3, K_B};
use crate::math::{Quat, Vec3};
use crate::rng::{Purpose, Stream};
use std::io::Write;
use std::time::Instant;

/// One atomic unit of time in picoseconds.
const PS: f64 = 2.418884326585747e-5;
/// Hartree in kcal/mol.
const KCAL: f64 = 627.509474;

fn random_rotation(s: &mut Stream) -> Quat {
    let (u1, u2, u3) = (s.uniform(), s.uniform(), s.uniform());
    let tau = std::f64::consts::TAU;
    let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
    Quat { w: a * (tau * u2).sin(), v: Vec3 { x: a * (tau * u2).cos(), y: b * (tau * u3).sin(), z: b * (tau * u3).cos() } }
}

/// The rotation taking the principal-frame vectors `(b1, b2)` to `(a1, a2)`, as a quaternion
/// (a three-atom molecule's orientation from its atoms; `PHYS_START_SNAP`).
fn rotation_between(b1: Vec3, b2: Vec3, a1: Vec3, a2: Vec3) -> Quat {
    let frame = |u1: Vec3, u2: Vec3| -> [Vec3; 3] {
        let e1 = u1.unit();
        let e3 = u1.cross(u2).unit();
        let e2 = e3.cross(e1);
        [e1, e2, e3]
    };
    let (a, b) = (frame(a1, a2), frame(b1, b2));
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
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let law_file = args.get(1).cloned().unwrap_or_else(|| format!("law-{name}.txt"));
    let temperature: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(298.0);
    let n: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(216);
    let density_g: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let total_ps: f64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(60.0);
    let law_text = std::fs::read_to_string(&law_file).unwrap_or_else(|_| panic!("no {law_file}"));
    let law = SiteSite::from_text(&law_text).expect("a readable law");
    let bisector = SiteSite::bisector_from_text(&law_text);
    // A law whose text carries `alpha` lines has induced dipoles.
    let alpha = SiteSite::alpha_from_text(&law_text, law.charge.len());
    let polarised = Polarisable { law: &law, alpha: alpha.clone() };
    let induced = alpha.iter().any(|a| *a > 0.0);
    let law_ref: &dyn SiteLaw = if induced { &polarised } else { &law };
    if induced {
        println!("induced dipoles on, polarisabilities by site type {alpha:?} bohr^3");
    }
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a growth state")).expect("readable");
    let z: Vec<u32> = match name.as_str() {
        "water" => vec![8, 1, 1],
        "methane" => vec![6, 1, 1, 1, 1],
        "ammonia" => vec![7, 1, 1, 1],
        "methanol" => vec![6, 8, 1, 1, 1, 1],
        other => panic!("no atoms known for {other}"),
    };
    let mono = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let atom_types = equivalent_atoms(&mono);
    let frame = crate::liquid::extra_sites_from_text(&law_text);
    let kind = Kind::of_molecule_sites(&z, &state.positions, &atom_types, bisector, &frame);
    let atoms = z.len();
    // The box: n molecules at the density asked.
    let molecule_g = kind.mass / AMU * 1.66053906660e-24;
    let per_bohr3 = density_g / molecule_g * 0.529177210903e-8f64.powi(3);
    let side = (n as f64 / per_bohr3).cbrt();
    // `PHYS_RCUT` (bohr) overrides the 16-bohr cap, for a study of what the
    // truncation does (PLAY.md E8d); the default is unchanged.
    let cap: f64 = std::env::var("PHYS_RCUT").ok().and_then(|v| v.parse().ok()).unwrap_or(16.0);
    let r_cut = (0.49 * side).min(cap);
    let r_on = r_cut - 2.0;
    let mut l = Liquid::slab(kind, n, per_bohr3, side, side, temperature, 0xb01c ^ temperature.to_bits(), r_on, r_cut);
    println!("{name} at {temperature} K: {n} molecules, {density_g} g/cm3, box {side:.2} bohr ({:.2} A), cut-off {r_cut:.1} bohr", side * 0.529177);
    // `PHYS_EWALD=<s' in bohr>` carries the long-range part of the charges and
    // the induced dipoles in reciprocal space (ewald.rs, PLAY.md E8d); the pair
    // terms keep only their short-range Coulomb remainder. Unset, the run is
    // the truncated one it has always been.
    let ewald_law;
    let law_ref: &dyn SiteLaw = match std::env::var("PHYS_EWALD").ok().and_then(|v| v.parse::<f64>().ok()) {
        Some(sp) => {
            ewald_law = EwaldLaw { inner: law_ref, ewald: Ewald::new(sp, l.cell, std::env::var("PHYS_EWALD_ACC").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-7)) };
            println!("  Ewald sum on: s' = {sp} bohr, {} reciprocal vectors", ewald_law.ewald.nk());
            &ewald_law
        }
        None => law_ref,
    };

    // Turn each molecule until no atom of it is within 2 A of an atom of any
    // before it.
    let mut rng = Stream::at(0xc1ea, 0, 0, Purpose::Positions);
    let min = 2.0 / 0.529177210903;
    let image = |l: &Liquid, mut d: Vec3| {
        d.x -= l.cell[0] * (d.x / l.cell[0]).round();
        d.y -= l.cell[1] * (d.y / l.cell[1]).round();
        d.z -= l.cell[2] * (d.z / l.cell[2]).round();
        d
    };
    let mut unplaced = 0;
    for i in 0..n {
        let mut ok = false;
        for _ in 0..2000 {
            l.orientation[i] = random_rotation(&mut rng);
            let mine: Vec<Vec3> = l.kind.sites.iter().take(atoms).map(|s| l.orientation[i].rotate(*s)).collect();
            let clear = (0..i).all(|j| {
                let d = image(&l, l.com[i] - l.com[j]);
                if d.norm() > 12.0 {
                    return true;
                }
                l.kind.sites.iter().take(atoms).all(|sj| {
                    let pj = l.orientation[j].rotate(*sj);
                    mine.iter().all(|pi| (d + *pi - pj).norm() > min)
                })
            });
            if clear {
                ok = true;
                break;
            }
        }
        if !ok {
            unplaced += 1;
        }
    }
    if unplaced > 0 {
        println!("  {unplaced} molecules could not be turned clear of their neighbours");
    }
    // `PHYS_START_SNAP=<snapshot>`: start from an earlier run's configuration (a three-atom
    // molecule's centre and orientation fitted to its atoms), the centres scaled to this box -
    // a test of whether a structure is the equilibrium one or the one a quench from a random
    // lattice selects. Unset, nothing changes.
    if let Ok(path) = std::env::var("PHYS_START_SNAP") {
        assert_eq!(atoms, 3, "PHYS_START_SNAP is for three-atom molecules");
        let text = std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("no {path}"));
        let mut old_side = 0.0;
        let mut mols: Vec<Vec<Vec3>> = Vec::new();
        for line in text.lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            match w.first().copied() {
                Some("box") => old_side = w[1].parse().unwrap(),
                Some("mol") => {
                    let v: Vec<f64> = w[2..].iter().map(|x| x.parse().unwrap()).collect();
                    mols.push(v.chunks(3).map(|c| Vec3 { x: c[0], y: c[1], z: c[2] }).collect());
                }
                _ => {}
            }
        }
        assert_eq!(mols.len(), n, "the snapshot holds {} molecules, the run {n}", mols.len());
        let sites = l.kind.sites.clone();
        let scale = l.cell[0] / old_side;
        for (i, m) in mols.iter().enumerate() {
            let q = rotation_between(sites[1] - sites[0], sites[2] - sites[0], m[1] - m[0], m[2] - m[0]);
            l.orientation[i] = q;
            l.com[i] = (m[0] - q.rotate(sites[0])).scale(scale);
        }
        println!("  started from {path}, centres scaled by {scale:.4}");
    }

    let out = format!("bulk-{name}-{temperature}.txt");
    let mut file = std::fs::File::create(&out).expect("the output file");
    writeln!(file, "# {name} at {temperature} K, {n} molecules, {density_g} g/cm3, cut-off {r_cut} bohr, law {law_file}").ok();
    writeln!(file, "# ps  U_per_molecule_kcal  pressure_bar  T_translation  T_rotation").ok();
    let mut stream = Liquid::noise(0xba7 ^ temperature.to_bits());
    let mut f = l.forces(law_ref);
    println!("  start: U/N {:.2} kcal/mol", f.energy / n as f64 * KCAL);
    // Easing in: tiny steps and a strong bath, then longer.
    let tau = |ps: f64| 1.0 / (ps / PS);
    let t0 = Instant::now();
    for (label, dt, steps, friction) in [("easing", 2.0, 4000usize, tau(0.02)), ("settling", 10.0, 6000, tau(0.1)), ("warming", 25.0, 8000, tau(0.5))] {
        for _ in 0..steps {
            f = l.step(law_ref, f, dt, Some((temperature, friction, &mut stream)));
        }
        let (kt, kr) = l.kinetic();
        println!("  {label}: U/N {:.2} kcal/mol, T {:.0}/{:.0} K, P {:.0} bar", f.energy / n as f64 * KCAL, 2.0 * kt / (3.0 * n as f64 * K_B), 2.0 * kr / (3.0 * n as f64 * K_B), l.pressure(&f) * BAR_PER_HARTREE_PER_BOHR3);
        if !f.energy.is_finite() || f.energy / n as f64 * KCAL > 1000.0 {
            println!("  the run blew up during {label}; the law's wall is too steep for these steps or the start is bad");
            return;
        }
    }
    let dt = 40.0;
    let steps = (total_ps / PS / dt) as usize;
    let settle = steps / 5;
    let every = ((0.1 / PS) / dt) as usize;
    let mut samples: Vec<(f64, f64, f64, f64)> = Vec::new();
    // g(r) of the first atoms, 0.1 A bins out to half the shortest box edge.
    let bin = 0.1 / 0.529177210903;
    let r_max = 0.5 * l.cell.iter().cloned().fold(f64::INFINITY, f64::min);
    let nbins = (r_max / bin) as usize;
    let mut hist = vec![0.0f64; nbins];
    let mut frames = 0usize;
    for step in 0..steps {
        f = l.step(law_ref, f, dt, Some((temperature, tau(1.0), &mut stream)));
        if step >= settle && step % every == 0 {
            let (kt, kr) = l.kinetic();
            let s = (f.energy / n as f64 * KCAL, l.pressure(&f) * BAR_PER_HARTREE_PER_BOHR3, 2.0 * kt / (3.0 * n as f64 * K_B), 2.0 * kr / (3.0 * n as f64 * K_B));
            writeln!(file, "{:.3} {:.4} {:.1} {:.1} {:.1}", step as f64 * dt * PS, s.0, s.1, s.2, s.3).ok();
            file.flush().ok();
            samples.push(s);
            let first: Vec<Vec3> = (0..n).map(|i| l.com[i] + l.orientation[i].rotate(l.kind.sites[0])).collect();
            for i in 0..n {
                for j in i + 1..n {
                    let d = image(&l, first[j] - first[i]).norm();
                    if d < nbins as f64 * bin {
                        hist[(d / bin) as usize] += 2.0;
                    }
                }
            }
            frames += 1;
        }
        if step % (steps / 12).max(1) == 0 {
            println!("  {:.0}% after {:.0} s: U/N {:.2} kcal/mol, P {:.0} bar", 100.0 * step as f64 / steps as f64, t0.elapsed().as_secs_f64(), f.energy / n as f64 * KCAL, l.pressure(&f) * BAR_PER_HARTREE_PER_BOHR3);
            std::io::stdout().flush().ok();
            if !f.energy.is_finite() {
                println!("  the run blew up");
                return;
            }
        }
    }
    // Block averages for an honest error bar.
    let mean = |k: usize| samples.iter().map(|s| [s.0, s.1, s.2, s.3][k]).sum::<f64>() / samples.len() as f64;
    let blocks = 5usize;
    let per = samples.len() / blocks;
    let spread = |k: usize| {
        let b: Vec<f64> = (0..blocks).map(|c| samples[c * per..(c + 1) * per].iter().map(|s| [s.0, s.1, s.2, s.3][k]).sum::<f64>() / per as f64).collect();
        let m = b.iter().sum::<f64>() / blocks as f64;
        (b.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (blocks - 1) as f64 / blocks as f64).sqrt()
    };
    let rt = K_B * temperature * KCAL;
    println!("{name} at {temperature} K ({} samples over {:.0} ps, block standard errors):", samples.len(), total_ps * 0.8);
    println!("  U/N {:.3} +- {:.3} kcal/mol; heat of vaporisation -U/N + RT = {:.3} kcal/mol", mean(0), spread(0), -mean(0) + rt);
    println!("  pressure {:.0} +- {:.0} bar at {density_g} g/cm3 (a liquid at the right density sits near 1 bar; positive says it wants to expand)", mean(1), spread(1));
    println!("  temperatures: translation {:.1} K, rotation {:.1} K (bath {temperature} K)", mean(2), mean(3));
    // The dispersion beyond the cut-off, from the law's own C6 (isotropic, g = 1):
    // P = -2 pi rho^2 C6 I, I = integral from r_on of (1 - S) / r^4; and the same for U/N.
    {
        let types = law.charge.len();
        let mut c6 = 0.0;
        for ta in &l.kind.types {
            for tb in &l.kind.types {
                c6 += law.pair[ta * types + tb][2];
            }
        }
        let rho = n as f64 / (l.cell[0] * l.cell[1] * l.cell[2]);
        let (mut integral, mut r, h) = (0.0, l.r_on, 0.001);
        while r < 2000.0 {
            let (sw, _) = crate::liquid::switch(r + h / 2.0, l.r_on, l.r_cut);
            integral += (1.0 - sw) / (r + h / 2.0).powi(4) * h;
            r += h;
        }
        let p_tail = -2.0 * std::f64::consts::PI * rho * rho * c6 * integral * BAR_PER_HARTREE_PER_BOHR3;
        let u_tail = -2.0 * std::f64::consts::PI * rho * c6 * integral * KCAL;
        println!("  dispersion beyond the cut-off (C6 {c6:.2}): {p_tail:.0} bar and {u_tail:.3} kcal/mol a molecule, so pressure with it {:.0} bar, U/N {:.3}", mean(1) + p_tail, mean(0) + u_tail);
    }
    // g(r): the pairs found in each shell against what an ideal gas of the
    // same density would put there.
    let volume = l.cell[0] * l.cell[1] * l.cell[2];
    let rho = n as f64 / volume;
    let g: Vec<f64> = (0..nbins)
        .map(|k| {
            let (r0, r1) = (k as f64 * bin, (k + 1) as f64 * bin);
            let shell = 4.0 / 3.0 * std::f64::consts::PI * (r1.powi(3) - r0.powi(3));
            hist[k] / (frames as f64 * n as f64 * rho * shell)
        })
        .collect();
    let a = 0.529177210903;
    let gfile = format!("bulk-{name}-{temperature}.gOO.txt");
    let mut gf = std::fs::File::create(&gfile).expect("the g(r) file");
    writeln!(gf, "# g(r) of the first atoms, {frames} frames, {n} molecules, 0.1 A bins: r_A g").ok();
    for (k, v) in g.iter().enumerate() {
        writeln!(gf, "{:.3} {:.4}", (k as f64 + 0.5) * 0.1, v).ok();
    }
    let at = |lo: f64, hi: f64, max: bool| -> (f64, f64) {
        let ks = (lo / 0.1) as usize..((hi / 0.1) as usize).min(nbins);
        let pick = ks.map(|k| ((k as f64 + 0.5) * 0.1, g[k]));
        if max { pick.fold((0.0, f64::MIN), |b, x| if x.1 > b.1 { x } else { b }) } else { pick.fold((0.0, f64::MAX), |b, x| if x.1 < b.1 { x } else { b }) }
    };
    let (r1, g1) = at(2.3, 3.4, true);
    let (rm, gm) = at(r1 + 0.2, 4.4, false);
    let (r2, g2) = at(rm, 5.8, true);
    let coordination: f64 = (0..nbins).filter(|&k| (k as f64 + 1.0) * 0.1 <= rm).map(|k| rho * 4.0 * std::f64::consts::PI * ((k as f64 + 0.5) * bin).powi(2) * bin * g[k]).sum();
    println!("  g(first atoms): first peak {g1:.2} at {r1:.2} A, first minimum {gm:.2} at {rm:.2} A, second peak {g2:.2} at {r2:.2} A, coordination to the minimum {coordination:.2} ({gfile})");
    let _ = a;
    // The last configuration, for taking pairs out of (phys-pairs-liquid-gpu).
    let snap = format!("bulk-{name}-{temperature}.snap");
    let mut sf = std::fs::File::create(&snap).expect("the snapshot file");
    writeln!(sf, "box {} {} {}", l.cell[0], l.cell[1], l.cell[2]).ok();
    for i in 0..n {
        let mut line = format!("mol {i}");
        for site in l.kind.sites.iter().take(atoms) {
            let p = l.com[i] + l.orientation[i].rotate(*site);
            line += &format!(" {:.8} {:.8} {:.8}", p.x, p.y, p.z);
        }
        writeln!(sf, "{line}").ok();
    }
    println!("  snapshot in {snap}");
}
