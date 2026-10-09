//! A molecule's potential energy surface from its own electrons — `docs/PLAY.md`
//! E8c, the surface every treatment of a moving water molecule starts from.
//!
//! ```sh
//! cargo run --release --bin phys-pes -- water
//! ```
//!
//! Hartree-Fock and RI-MP2 (frozen core) energies of a water molecule over a
//! grid of its two bond lengths and its angle, each point in the engine's own
//! default basis and the same method the pair energies use. Writes
//! `pes-water.txt`, one line `r1 r2 theta E_hartree_fock E_mp2_correlation`
//! (bohr, radians, hartree) per point, and resumes: a point already in the file
//! is not done again. The grid has `r1 >= r2` only, the surface being symmetric
//! in the two hydrogens.

use phys::electrons::functional::Functional;
use phys::electrons::hf::{frozen_core, hartree_fock, mp2};
use phys::electrons::molecule::Molecule;

const ANGSTROM: f64 = 1.0 / 0.529177210903;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    assert_eq!(name, "water", "the surface is built for water first");
    let out = format!("pes-{name}.txt");
    let done: std::collections::HashSet<String> = std::fs::read_to_string(&out).unwrap_or_default().lines().filter(|l| !l.starts_with('#')).filter_map(|l| {
        let w: Vec<&str> = l.split_whitespace().collect();
        Some(format!("{}|{}|{}", w.first()?, w.get(1)?, w.get(2)?))
    }).collect();
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    use std::io::Write;
    if done.is_empty() {
        writeln!(file, "# water, Hartree-Fock + RI-MP2 (frozen core), default basis: r1 r2 theta (bohr, bohr, radians) E_HF E_MP2correlation (hartree)").ok();
    }
    // Bond lengths in angstrom, the angle in degrees.
    let r_grid = [0.86, 0.90, 0.94, 0.98, 1.02, 1.06, 1.10, 1.16];
    let t_grid: [f64; 7] = [88.0, 94.0, 100.0, 104.5, 110.0, 116.0, 124.0];
    let mut points = Vec::new();
    for (i, r1) in r_grid.iter().enumerate() {
        for r2 in &r_grid[..=i] {
            for t in &t_grid {
                points.push((r1 * ANGSTROM, r2 * ANGSTROM, t.to_radians()));
            }
        }
    }
    println!("{} points", points.len());
    for (r1, r2, t) in points {
        let key = format!("{r1}|{r2}|{t}");
        if done.contains(&key) {
            continue;
        }
        let start = std::time::Instant::now();
        let positions = vec![[0.0, 0.0, 0.0], [r1, 0.0, 0.0], [r2 * t.cos(), r2 * t.sin(), 0.0]];
        let mol = Molecule { z: vec![8, 1, 1], positions, charge: 0, unpaired: 0 };
        let p = mol.problem(Functional::Pbe);
        let hf = hartree_fock(&p, 200, 1e-10);
        assert!(hf.converged, "Hartree-Fock did not converge at {r1} {r2} {t}");
        let corr = mp2(&p, &hf, frozen_core(&mol.z, &[0, 1, 2]));
        writeln!(file, "{r1} {r2} {t} {:.12} {:.12}", hf.energy, corr.correlation).expect("written");
        file.flush().ok();
        println!("r1 {:.3} r2 {:.3} A, theta {:.1}: {:.8} hartree, {:.1} s", r1 / ANGSTROM, r2 / ANGSTROM, t.to_degrees(), hf.energy + corr.correlation, start.elapsed().as_secs_f64());
    }
}
