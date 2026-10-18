//! Refit a molecule's charge clouds to the exact first-order electrostatic energies of pairs
//! (`phys-es-gpu`) — `docs/SOLUTIONS.md` S1, `phys::cloud::fit_es`.
//!
//! ```sh
//! phys-cloudfit ammonia cloud-ammonia.cloud es-ammonia.txt -gpu-mp2,-gpu-wall-mp2 [--dipole 0.57]
//! ```
//!
//! Starts from the structure file `phys-cloud` wrote, reads the pairs `pairs-<name><tag>.txt` and the
//! energies, fits populations, widths and lone-pair distances (neutral, non-negative, the dipole held
//! if asked), and writes `cloud-<name>-fit.cloud` and `cloud-<name>-fit.txt` (a law's electrostatic lines).

use phys::electrons::grow::Resume;
use phys::math::Vec3;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args[0].clone();
    let start = phys::cloud::Cloud::from_text(&std::fs::read_to_string(&args[1]).unwrap_or_else(|_| panic!("no {}", args[1]))).expect("a readable cloud");
    let rows = phys::cloud::es_rows(&name, &args[2], &args[3]);
    let dipole = args.iter().position(|a| a == "--dipole").and_then(|i| args.get(i + 1)).and_then(|v| v.parse::<f64>().ok());
    println!("{name}: {} pairs with an exact electrostatic energy", rows.len());
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a state")).expect("readable");
    let pos: Vec<Vec3> = state.positions.iter().map(|r| Vec3 { x: r[0], y: r[1], z: r[2] }).collect();
    let fitted = phys::cloud::fit_es(&start, &pos, &rows, dipole, &mut |s| println!("{s}"));
    std::fs::write(format!("cloud-{name}-fit.cloud"), fitted.to_text()).expect("written");
    std::fs::write(format!("cloud-{name}-fit.txt"), fitted.law_text(&pos)).expect("written");
    println!("  written to cloud-{name}-fit.cloud and cloud-{name}-fit.txt");
}
