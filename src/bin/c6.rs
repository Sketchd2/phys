//! A molecule's atom-pair dispersion coefficients from its own electrons
//! (`docs/PLAY.md` E6, E8d): the non-local functional's `C6` between two molecules, split among
//! atom pairs by the density's partition (`electrons::partition`).
//!
//! ```sh
//! phys-c6 methane        # grow-methane.state in the working directory
//! ```
//!
//! Prints the class-pair values in the form `phys-fit --c6` takes (hartree bohr^6) and the
//! molecule-molecule total.

use phys::electrons::functional::Functional;
use phys::electrons::grow::{equivalent_atoms, Resume};
use phys::electrons::molecule::Molecule;
use phys::electrons::partition::{atom_pair_c6, partition, FreeAtom};
use phys::electrons::scf::solve;
use phys::electrons::vdw::Z_AB_DF1;

fn main() {
    let name = std::env::args().nth(1).unwrap_or_else(|| "water".into());
    let z: Vec<u32> = match name.as_str() {
        "water" => vec![8, 1, 1],
        "methane" => vec![6, 1, 1, 1, 1],
        "ammonia" => vec![7, 1, 1, 1],
        "methanol" => vec![6, 8, 1, 1, 1, 1],
        other => panic!("no atoms known for {other}"),
    };
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a growth state")).expect("readable");
    let mol = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let types = equivalent_atoms(&mol);
    let p = mol.problem(Functional::Pbe);
    let s = solve(&p, 200, 1e-10);
    let mut zs: Vec<u32> = z.clone();
    zs.sort();
    zs.dedup();
    let free: Vec<FreeAtom> = zs.iter().map(|&zz| FreeAtom::new(zz, Functional::Pbe)).collect();
    let part = partition(&mol, &p, &s, &free, 50, 12, Z_AB_DF1);
    let pairs = atom_pair_c6(&part, &part);
    let n = z.len();
    // One atom of each class stands for it.
    let classes = types.iter().max().map(|m| m + 1).unwrap_or(1);
    let rep: Vec<usize> = (0..classes).map(|c| types.iter().position(|t| *t == c).unwrap()).collect();
    let mut line = Vec::new();
    for a in 0..classes {
        for b in a..classes {
            line.push(format!("{a}-{b}={:.4}", pairs[rep[a]][rep[b]]));
        }
    }
    let total: f64 = (0..n).flat_map(|i| (0..n).map(move |j| (i, j))).map(|(i, j)| pairs[i][j]).sum();
    println!("{name}: atom classes {types:?}");
    println!("--c6 {}", line.join(","));
    println!("molecule-molecule C6 {total:.3} hartree bohr^6");
}
