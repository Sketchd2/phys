//! A liquid slab against its own vapour at one temperature — `PLAY.md` E8,
//! direct coexistence, by the owner's decision.
//!
//! ```sh
//! cargo run --release --bin phys-coexist -- water law-water.txt 400 [molecules] [cut-off bohr] [ps]
//! ```
//!
//! The molecule at the shape its growth relaxed it to (`grow-<name>.state`),
//! its atoms typed by symmetry class, moved by the law in the given file
//! (`SiteSite::to_text`). A slab is built at a starting density, brought to
//! the temperature by the thermostat for a fifth of the run, and then
//! sampled: every sample appends the liquid's and the vapour's densities and
//! the potential energy per molecule to `coexist-<name>-<T>.txt`, so a run
//! stopped part-way keeps what it measured.
//!
//! What is read off afterwards: the vapour's pressure, from its density by
//! the ideal-gas law (good at these pressures, poor near the critical point),
//! and the heat of vaporisation from the energies; Clausius-Clapeyron between
//! temperatures then gives the boiling point at any pressure.

use phys::electrons::grow::{equivalent_atoms, Resume};
use phys::electrons::molecule::Molecule;
use phys::liquid::{coexisting_densities, Kind, Liquid, SiteSite, AMU, K_B};
use std::io::Write;
use std::time::Instant;

/// One atomic unit of time in picoseconds.
const PS: f64 = 2.418884326585747e-5;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let law_file = args.get(1).cloned().unwrap_or_else(|| format!("law-{name}.txt"));
    let temperature: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(300.0);
    let n: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(512);
    let r_cut: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(14.0);
    let total_ps: f64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(500.0);
    let law = SiteSite::from_text(&std::fs::read_to_string(&law_file).unwrap_or_else(|_| panic!("no {law_file}"))).expect("a readable law");
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a growth state")).expect("readable");
    let z: Vec<u32> = match name.as_str() {
        "water" => vec![8, 1, 1],
        "methane" => vec![6, 1, 1, 1, 1],
        "ammonia" => vec![7, 1, 1, 1],
        "methanol" => vec![6, 8, 1, 1, 1, 1],
        other => panic!("no atoms known for {other}"),
    };
    let mono = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let types = equivalent_atoms(&mono);
    let masses: Vec<f64> = z.iter().map(|&zz| phys::chem::elements::Element(zz as u8).mass_kg().expect("a mass") / 1.66053906660e-27 * AMU).collect();
    let kind = Kind::from_atoms(&state.positions, &masses, &types);
    // A starting density from the molecule's mass at about a liquid's 0.8
    // g/cm^3; the slab finds its own as it settles.
    let molecule_g = kind.mass / AMU * 1.66053906660e-24;
    let start_density = 0.8 / molecule_g * 0.529177210903e-8f64.powi(3);
    let side = 2.5 * r_cut;
    let thickness = n as f64 / (start_density * side * side);
    let length = 3.0 * thickness + 2.0 * r_cut;
    let mut l = Liquid::slab(kind, n, start_density, side, length, temperature, 0x5eed ^ temperature.to_bits(), r_cut - 2.0, r_cut);
    let dt = 40.0;
    let steps = (total_ps / PS / dt) as usize;
    let settle = steps / 5;
    let every = ((0.5 / PS) / dt) as usize;
    let out = format!("coexist-{name}-{temperature}.txt");
    let mut file = std::fs::File::create(&out).expect("the output file");
    writeln!(file, "# {name} at {temperature} K: {n} molecules, box {side:.2} x {side:.2} x {length:.2} bohr, cut-off {r_cut} bohr, step {dt} a.u.; law {law_file}").ok();
    writeln!(file, "# ps liquid_density vapour_density (molecules/bohr^3) potential_per_molecule (hartree) temperature_K").ok();
    println!("{name} at {temperature} K: {n} molecules, {steps} steps ({total_ps} ps), settling for {settle}");
    let mut stream = Liquid::noise(0xba7 ^ temperature.to_bits());
    let mut f = l.forces(&law);
    let t0 = Instant::now();
    let friction = 1.0 / (1.0 / PS);
    for step in 0..steps {
        f = l.step(&law, f, dt, Some((temperature, friction, &mut stream)));
        if step >= settle && step % every == 0 {
            let (liquid, vapour) = coexisting_densities(&l.density_profile(120), 0.2, 0.25);
            let (kt, kr) = l.kinetic();
            let t_now = 2.0 * (kt + kr) / (6.0 * n as f64 * K_B);
            writeln!(file, "{:.3} {liquid:.6e} {vapour:.6e} {:.8e} {t_now:.2}", step as f64 * dt * PS, f.energy / n as f64).ok();
            file.flush().ok();
        }
        if step % (steps / 20).max(1) == 0 {
            println!("  {:.0}% after {:.0} s", 100.0 * step as f64 / steps as f64, t0.elapsed().as_secs_f64());
            std::io::stdout().flush().ok();
        }
    }
    println!("{name} at {temperature} K: done in {:.0} s; samples in {out}", t0.elapsed().as_secs_f64());
}
