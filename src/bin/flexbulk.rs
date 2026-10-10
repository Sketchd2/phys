//! A liquid whose molecules move atom by atom — `PLAY.md` E8c, treatments (2)
//! and (3). The flexible counterpart of `phys-bulk`: the same law, the same
//! protocol and the same measurements, with the molecules' own energy surface
//! in place of a rigid body.
//!
//! ```sh
//! cargo run --release --bin phys-flexbulk -- law-water.txt 298 [molecules] [g/cm3] [ps] \
//!     [--bond ANGSTROM] [--dt AU] [--surface surface-water.txt] [--time STEPS]
//! ```
//!
//! `--bond L` holds both O-H lengths at `L` angstrom (SHAKE and RATTLE) and
//! leaves the angle free: treatment (2). Without it all three internal
//! coordinates move: treatment (3). `--dt` is the production step in atomic
//! units of time (default 40 with the bonds held, the rigid run's step, and 20
//! without, the O-H stretch's period being 9 fs). `--time N` builds the liquid,
//! eases it briefly and reports the cost of `N` production steps, nothing more.
//!
//! **What it reports beyond `phys-bulk`.** The molecules' own energy, the mean
//! bond length and angle, and a heat of vaporisation that includes the change
//! in a molecule's own energy between the gas and the liquid: the classical
//! canonical average of the surface over the gas molecule's internal
//! coordinates (`gas_intra`), not an assumed `kT/2` a mode. Everything is
//! classical: no zero-point motion (the path-integral step of E8c is the
//! comparison for that).

use phys::electrons::grow::{equivalent_atoms, Resume};
use phys::electrons::molecule::Molecule;
use phys::flex::{intramolecular, Flex, FlexSites, Surface};
use phys::liquid::{Kind, Liquid, Polarisable, SiteSite, AMU, BAR_PER_HARTREE_PER_BOHR3, K_B};
use phys::math::{Quat, Vec3};
use phys::rng::{Purpose, Stream};
use std::io::Write;
use std::time::Instant;

const PS: f64 = 2.418884326585747e-5;
const KCAL: f64 = 627.509474;
const BOHR: f64 = 0.529177210903;

fn random_rotation(s: &mut Stream) -> Quat {
    let (u1, u2, u3) = (s.uniform(), s.uniform(), s.uniform());
    let tau = std::f64::consts::TAU;
    let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
    Quat { w: a * (tau * u2).sin(), v: Vec3 { x: a * (tau * u2).cos(), y: b * (tau * u3).sin(), z: b * (tau * u3).cos() } }
}

/// A gas molecule's classical mean internal energy (hartree) on the surface at
/// `kt`, over the volume element `r1^2 r2^2 sin(theta)` of the internal
/// coordinates (`r^2 sin(theta)` less the held bond lengths when they are held),
/// by quadrature over the region the surface was fitted in.
fn gas_intra(surface: &Surface, kt: f64, bond: Option<f64>) -> f64 {
    let (r_lo, r_hi) = (0.86 / BOHR, 1.16 / BOHR);
    let (t_lo, t_hi) = (88.0f64.to_radians(), 124.0f64.to_radians());
    let (mut z, mut ez) = (0.0, 0.0);
    let nt = 400;
    match bond {
        Some(r) => {
            for k in 0..nt {
                let t = t_lo + (t_hi - t_lo) * (k as f64 + 0.5) / nt as f64;
                let e = surface.eval(r, r, t).0;
                let w = (-e / kt).exp() * t.sin();
                z += w;
                ez += w * e;
            }
        }
        None => {
            let nr = 60;
            let nt = 120;
            for i in 0..nr {
                let r1 = r_lo + (r_hi - r_lo) * (i as f64 + 0.5) / nr as f64;
                for j in 0..nr {
                    let r2 = r_lo + (r_hi - r_lo) * (j as f64 + 0.5) / nr as f64;
                    for k in 0..nt {
                        let t = t_lo + (t_hi - t_lo) * (k as f64 + 0.5) / nt as f64;
                        let e = surface.eval(r1, r2, t).0;
                        let w = (-e / kt).exp() * r1 * r1 * r2 * r2 * t.sin();
                        z += w;
                        ez += w * e;
                    }
                }
            }
        }
    }
    ez / z
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let flag = |n: &str| args.iter().position(|a| a == n).and_then(|i| args.get(i + 1)).cloned();
    let positional: Vec<&String> = {
        let mut v = Vec::new();
        let mut skip = false;
        for a in &args {
            if skip {
                skip = false;
            } else if a.starts_with("--") {
                skip = true;
            } else {
                v.push(a);
            }
        }
        v
    };
    let law_file = positional.first().map(|s| s.to_string()).expect("a law file");
    let temperature: f64 = positional.get(1).and_then(|s| s.parse().ok()).unwrap_or(298.0);
    let n: usize = positional.get(2).and_then(|s| s.parse().ok()).unwrap_or(216);
    let density_g: f64 = positional.get(3).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let total_ps: f64 = positional.get(4).and_then(|s| s.parse().ok()).unwrap_or(60.0);
    let bond: Option<f64> = flag("--bond").and_then(|v| v.parse::<f64>().ok()).map(|a| a / BOHR);
    let dt: f64 = flag("--dt").and_then(|v| v.parse().ok()).unwrap_or(if bond.is_some() { 40.0 } else { 20.0 });
    let surface_file = flag("--surface").unwrap_or_else(|| "surface-water.txt".into());
    let timing: Option<usize> = flag("--time").and_then(|v| v.parse().ok());
    let tag = flag("--tag").unwrap_or_else(|| if bond.is_some() { "flex2".into() } else { "flex3".into() });

    let law_text = std::fs::read_to_string(&law_file).unwrap_or_else(|_| panic!("no {law_file}"));
    let law = SiteSite::from_text(&law_text).expect("a readable law");
    let bisector = SiteSite::bisector_from_text(&law_text);
    let alpha = SiteSite::alpha_from_text(&law_text, law.charge.len());
    // `--no-induction` leaves the induced dipoles out: for timing only, to say how much of a step they are.
    let alpha = if args.iter().any(|a| a == "--no-induction") { vec![0.0; alpha.len()] } else { alpha };
    let law_ref = Polarisable { law: &law, alpha: alpha.clone() };
    let surface = Surface::from_text(&std::fs::read_to_string(&surface_file).unwrap_or_else(|_| panic!("no {surface_file}")));

    let state = Resume::from_text(&std::fs::read_to_string("grow-water.state").expect("a growth state")).expect("readable");
    let z: Vec<u32> = vec![8, 1, 1];
    let mono = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let atom_types = equivalent_atoms(&mono);
    let frame = phys::liquid::extra_sites_from_text(&law_text);
    let kind = Kind::of_molecule_sites(&z, &state.positions, &atom_types, bisector, &frame);
    let masses: Vec<f64> = z.iter().map(|&zz| phys::chem::elements::Element(zz as u8).mass_kg().expect("a mass") / 1.66053906660e-27 * AMU).collect();
    let molecule_g = kind.mass / AMU * 1.66053906660e-24;
    let per_bohr3 = density_g / molecule_g * 0.529177210903e-8f64.powi(3);
    let side = (n as f64 / per_bohr3).cbrt();
    let r_cut = (0.49 * side).min(16.0);
    let r_on = r_cut - 2.0;
    let mut lattice = Liquid::slab(kind, n, per_bohr3, side, side, temperature, 0xb01c ^ temperature.to_bits(), r_on, r_cut);
    println!(
        "water at {temperature} K: {n} molecules, {density_g} g/cm3, box {side:.2} bohr, cut-off {r_cut:.1} bohr; {} (production step {dt} au = {:.2} fs); {} threads",
        match bond {
            Some(b) => format!("bonds held at {:.5} A, angle free", b * BOHR),
            None => "all three internal coordinates free".into(),
        },
        dt * PS * 1000.0,
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    );

    // Turn each molecule until no atom is within 2 A of an earlier one's.
    let mut rng = Stream::at(0xc1ea, 0, 0, Purpose::Positions);
    let min = 2.0 / BOHR;
    let image = |cell: [f64; 3], mut d: Vec3| {
        d.x -= cell[0] * (d.x / cell[0]).round();
        d.y -= cell[1] * (d.y / cell[1]).round();
        d.z -= cell[2] * (d.z / cell[2]).round();
        d
    };
    for i in 0..n {
        for _ in 0..2000 {
            lattice.orientation[i] = random_rotation(&mut rng);
            let mine: Vec<Vec3> = lattice.kind.sites.iter().take(3).map(|s| lattice.orientation[i].rotate(*s)).collect();
            let clear = (0..i).all(|j| {
                let d = image(lattice.cell, lattice.com[i] - lattice.com[j]);
                d.norm() > 12.0 || lattice.kind.sites.iter().take(3).all(|sj| {
                    let pj = lattice.orientation[j].rotate(*sj);
                    mine.iter().all(|pi| (d + *pi - pj).norm() > min)
                })
            });
            if clear {
                break;
            }
        }
    }
    let mut pos: Vec<[Vec3; 3]> = (0..n)
        .map(|i| {
            let mut m = [Vec3::ZERO; 3];
            for a in 0..3 {
                m[a] = lattice.com[i] + lattice.orientation[i].rotate(lattice.kind.sites[a]);
            }
            if let Some(len) = bond {
                for h in [1usize, 2] {
                    m[h] = m[0] + (m[h] - m[0]).unit().scale(len);
                }
            }
            m
        })
        .collect();
    let _ = &mut pos;
    let mut s = Flex {
        sites: FlexSites::from_text(&law_text, [atom_types[0], atom_types[1], atom_types[2]]),
        masses: [masses[0], masses[1], masses[2]],
        cell: lattice.cell,
        pos,
        vel: vec![[Vec3::ZERO; 3]; n],
        r_on,
        r_cut,
        surface: surface.clone(),
        bond,
    };
    let kt = K_B * temperature;
    let tau = |ps: f64| 1.0 / (ps / PS);
    let mut stream = Liquid::noise(0xba7 ^ temperature.to_bits());
    let t_force = Instant::now();
    let mut f = s.forces(&law_ref, None);
    println!("  start: U/N {:.2} kcal/mol (one force evaluation took {:.2} s)", f.inter / n as f64 * KCAL, t_force.elapsed().as_secs_f64());
    let dof = s.degrees_of_freedom();

    // Easing in, the rigid run's schedule scaled to this production step.
    let k = dt / 40.0;
    let t0 = Instant::now();
    let stages: Vec<(&str, f64, usize, f64)> = if timing.is_some() {
        vec![("easing", 2.0 * k, 300, tau(0.02)), ("settling", 10.0 * k, 300, tau(0.1))]
    } else {
        vec![("easing", 2.0 * k, 4000, tau(0.02)), ("settling", 10.0 * k, 6000, tau(0.1)), ("warming", 25.0 * k, 8000, tau(0.5))]
    };
    for (label, h, steps, friction) in stages {
        let t = Instant::now();
        for _ in 0..steps {
            f = s.step(&law_ref, f, h, Some((kt, friction, &mut stream)));
        }
        let (kc, ka) = s.kinetic();
        println!(
            "  {label}: U/N {:.2} kcal/mol, intra {:.2}, T {:.0} (centres) / {:.0} (all) K, P {:.0} bar; {:.1} ms a step",
            f.inter / n as f64 * KCAL,
            f.intra / n as f64 * KCAL,
            2.0 * kc / (3.0 * n as f64 * K_B),
            2.0 * ka / (dof * K_B),
            s.pressure(&f) * BAR_PER_HARTREE_PER_BOHR3,
            t.elapsed().as_secs_f64() * 1000.0 / steps as f64
        );
        if !f.inter.is_finite() || f.inter / n as f64 * KCAL > 1000.0 {
            println!("  the run blew up during {label}");
            return;
        }
    }
    if let Some(steps) = timing {
        let t = Instant::now();
        for _ in 0..steps {
            f = s.step(&law_ref, f, dt, Some((kt, tau(1.0), &mut stream)));
        }
        let secs = t.elapsed().as_secs_f64();
        println!("  TIMING: {steps} production steps of {dt} au in {secs:.1} s = {:.1} ms a step; 100 ps would take {:.1} h", secs * 1000.0 / steps as f64, secs / steps as f64 * (100.0 / PS / dt) / 3600.0);
        return;
    }

    let out = format!("bulk-water-{temperature}-{tag}.txt");
    let mut file = std::fs::File::create(&out).expect("the output file");
    writeln!(file, "# water at {temperature} K, {n} molecules, {density_g} g/cm3, law {law_file}, {}", if bond.is_some() { "bonds held" } else { "fully flexible" }).ok();
    writeln!(file, "# ps  U_inter_per_molecule_kcal  E_intra_per_molecule_kcal  pressure_bar  T_centres  T_all  r_OH_A  theta_deg").ok();
    let steps = (total_ps / PS / dt) as usize;
    let settle = steps / 5;
    let every = (((0.1 / PS) / dt) as usize).max(1);
    // (U_inter, E_intra, P, T_centres, T_all, r_OH, theta) per sample.
    let mut samples: Vec<[f64; 7]> = Vec::new();
    let bin = 0.1 / BOHR;
    let r_max = 0.5 * s.cell.iter().cloned().fold(f64::INFINITY, f64::min);
    let nbins = (r_max / bin) as usize;
    let mut hist = vec![0.0f64; nbins];
    let mut frames = 0usize;
    let t1 = Instant::now();
    for step in 0..steps {
        f = s.step(&law_ref, f, dt, Some((kt, tau(1.0), &mut stream)));
        if step >= settle && step % every == 0 {
            let (kc, ka) = s.kinetic();
            let (mut r_sum, mut th_sum) = (0.0, 0.0);
            for m in &s.pos {
                let (a, b) = (m[1] - m[0], m[2] - m[0]);
                r_sum += a.norm() + b.norm();
                th_sum += (a.dot(b) / (a.norm() * b.norm())).clamp(-1.0, 1.0).acos();
            }
            let sample = [
                f.inter / n as f64 * KCAL,
                f.intra / n as f64 * KCAL,
                s.pressure(&f) * BAR_PER_HARTREE_PER_BOHR3,
                2.0 * kc / (3.0 * n as f64 * K_B),
                2.0 * ka / (dof * K_B),
                r_sum / (2.0 * n as f64) * BOHR,
                th_sum / n as f64,
            ];
            writeln!(file, "{:.3} {:.4} {:.4} {:.1} {:.1} {:.1} {:.5} {:.3}", step as f64 * dt * PS, sample[0], sample[1], sample[2], sample[3], sample[4], sample[5], sample[6].to_degrees()).ok();
            file.flush().ok();
            samples.push(sample);
            let first: Vec<Vec3> = (0..n).map(|i| s.pos[i][0]).collect();
            for i in 0..n {
                for j in i + 1..n {
                    let d = image(s.cell, first[j] - first[i]).norm();
                    if d < nbins as f64 * bin {
                        hist[(d / bin) as usize] += 2.0;
                    }
                }
            }
            frames += 1;
        }
        if step % (steps / 12).max(1) == 0 {
            println!("  {:.0}% after {:.0} s: U/N {:.2} kcal/mol, P {:.0} bar", 100.0 * step as f64 / steps as f64, t1.elapsed().as_secs_f64(), f.inter / n as f64 * KCAL, s.pressure(&f) * BAR_PER_HARTREE_PER_BOHR3);
            std::io::stdout().flush().ok();
            if !f.inter.is_finite() {
                println!("  the run blew up");
                return;
            }
        }
    }
    let ns = samples.len();
    let mean = |k: usize| samples.iter().map(|x| x[k]).sum::<f64>() / ns as f64;
    let blocks = 5usize;
    let per = ns / blocks;
    let spread = |k: usize| {
        let b: Vec<f64> = (0..blocks).map(|c| samples[c * per..(c + 1) * per].iter().map(|x| x[k]).sum::<f64>() / per as f64).collect();
        let m = b.iter().sum::<f64>() / blocks as f64;
        (b.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (blocks - 1) as f64 / blocks as f64).sqrt()
    };
    let rt = K_B * temperature * KCAL;
    let gas = gas_intra(&surface, kt, bond) * KCAL;
    let dh = -(mean(0) + mean(1) - gas) + rt;
    println!("water at {temperature} K ({ns} samples over {:.0} ps, block standard errors; {:.1} s of production):", total_ps * 0.8, t1.elapsed().as_secs_f64());
    println!("  U_inter/N {:.3} +- {:.3} kcal/mol; E_intra/N {:.3} +- {:.3} (a classical gas molecule's: {gas:.3})", mean(0), spread(0), mean(1), spread(1));
    println!("  heat of vaporisation -(U_inter/N + E_intra/N - E_intra_gas) + RT = {dh:.3} kcal/mol");
    println!("  pressure {:.0} +- {:.0} bar at {density_g} g/cm3", mean(2), spread(2));
    println!("  temperatures: centres {:.1} K, all atoms {:.1} K (bath {temperature} K)", mean(3), mean(4));
    println!("  mean r_OH {:.4} +- {:.4} A, mean H-O-H {:.2} degrees (the monomer's relaxed values: 0.9615 A, 104.15)", mean(5), spread(5), mean(6).to_degrees());
    let volume = s.cell[0] * s.cell[1] * s.cell[2];
    let rho = n as f64 / volume;
    let g: Vec<f64> = (0..nbins)
        .map(|k| {
            let (r0, r1) = (k as f64 * bin, (k + 1) as f64 * bin);
            let shell = 4.0 / 3.0 * std::f64::consts::PI * (r1.powi(3) - r0.powi(3));
            hist[k] / (frames as f64 * n as f64 * rho * shell)
        })
        .collect();
    let gfile = format!("bulk-water-{temperature}-{tag}.gOO.txt");
    let mut gf = std::fs::File::create(&gfile).expect("the g(r) file");
    writeln!(gf, "# g(r) of the oxygens, {frames} frames, {n} molecules, 0.1 A bins: r_A g").ok();
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
    println!("  g(OO): first peak {g1:.2} at {r1:.2} A, first minimum {gm:.2} at {rm:.2} A, second peak {g2:.2} at {r2:.2} A, coordination to the minimum {coordination:.2} ({gfile})");
    let snap = format!("bulk-water-{temperature}-{tag}.snap");
    let mut sf = std::fs::File::create(&snap).expect("the snapshot file");
    writeln!(sf, "box {} {} {}", s.cell[0], s.cell[1], s.cell[2]).ok();
    for i in 0..n {
        let mut line = format!("mol {i}");
        for p in &s.pos[i] {
            line += &format!(" {:.8} {:.8} {:.8}", p.x, p.y, p.z);
        }
        writeln!(sf, "{line}").ok();
    }
    println!("  snapshot in {snap}");
    let _ = intramolecular;
}
