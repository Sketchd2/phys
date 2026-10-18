//! A molecule's charge clouds from its own electrons, for any molecule — `docs/SOLUTIONS.md` S1.
//!
//! ```sh
//! phys-cloud ammonia [--mp2-dipole]
//! ```
//!
//! Reads `grow-<name>.state`, fits the core-and-cloud model of `phys::cloud` to the molecule's
//! electrostatic potential, prints it, and writes `cloud-<name>.txt`: the `charges`, `sigma` and
//! `site` lines of a law's text, a starting point for the refit to the pairs' exact first-order
//! energies. The widths and populations are found, none of them typed.

use phys::electrons::grow::Resume;
use phys::math::Vec3;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let mp2 = args.iter().any(|a| a == "--mp2-dipole");
    let cloud = phys::cloud::fit(&name, mp2, &mut |s| println!("{s}"));
    println!("  valence electrons {}, cores {:?}", cloud.valence, cloud.core);
    for (c, (p, w)) in cloud.atom_pop.iter().zip(&cloud.atom_width).enumerate() {
        println!("  atom class {c}: cloud {p:.3} electrons, width {w:.3} bohr");
    }
    for lp in &cloud.lone {
        println!("  lone-pair cloud of class {}: {} pair(s) on {} atom(s), {:.3} electrons, width {:.3} bohr, {:.3} bohr from the nucleus (the orbitals' centre is {:.3})", lp.class, lp.pairs, lp.directions.len(), lp.population, lp.width, lp.distance, lp.centroid_distance);
    }
    println!("  potential RMS {:.3e} hartree/e (atom-centred point charges alone: {:.3e})", cloud.rms, cloud.rms_points);
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a state")).expect("readable");
    let pos: Vec<Vec3> = state.positions.iter().map(|r| Vec3 { x: r[0], y: r[1], z: r[2] }).collect();
    std::fs::write(format!("cloud-{name}.txt"), cloud.law_text(&pos)).expect("written");
    std::fs::write(format!("cloud-{name}.cloud"), cloud.to_text()).expect("written");
    println!("  written to cloud-{name}.txt");
}
