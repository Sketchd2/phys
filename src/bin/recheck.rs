//! Check a pairs file against itself and against the calculation — see
//! `phys::pairs::recheck`.
//!
//! ```sh
//! phys-recheck pairs-water.txt water --sample 20 --seed 1
//! phys-recheck pairs-water.txt water --indices 124,125
//! ```
//!
//! For a file made with the GPU's final non-local energy use
//! `phys-recheck-gpu`, which solves it the same way.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(phys::pairs::recheck(&args));
}
