//! How much of a dimer's binding the basis leaves out (`docs/PLAY.md` E8f).
//!
//! ```sh
//! phys-shortfall water [scale ...]     # also ammonia, methane; scales default to 1
//! ```
//!
//! The S22 dimer of the molecule, the second monomer moved along the line of
//! the two heavy atoms so that their distance is `scale` times its S22 value,
//! Hartree-Fock + MP2 (frozen cores) twice over: each monomer solved in the
//! dimer's basis (the counterpoise correction, which approaches the complete
//! basis from below) and each in its own basis alone (which approaches it from
//! above, because the dimer borrows the partner's functions). The two bracket
//! the complete-basis answer and their gap is the measure of what the basis
//! is missing at that geometry. One line per scale: the distance, then
//! counterpoise HF, MP2 correlation and total, then the same uncorrected, in
//! kcal/mol.

use phys::electrons::functional::Functional;
use phys::electrons::hf::{frozen_core, hartree_fock, mp2};
use phys::electrons::molecule::Molecule;

type Dimer = (&'static str, &'static [u32], &'static [[f64; 3]], usize, f64);
const DIMERS: &[Dimer] = &[
    ("ammonia", &[7, 1, 1, 1, 7, 1, 1, 1], &[[-1.578718, -0.046611, 0.0], [-2.158621, 0.136396, -0.809565], [-2.158621, 0.136396, 0.809565], [-0.849471, 0.658193, 0.0], [1.578718, 0.046611, 0.0], [2.158621, -0.136396, -0.809565], [0.849471, -0.658193, 0.0], [2.158621, -0.136396, 0.809565]], 4, -3.1708),
    ("water", &[8, 1, 1, 8, 1, 1], &[[-1.551007, -0.11452, 0.0], [-1.934259, 0.762503, 0.0], [-0.599677, 0.040712, 0.0], [1.350625, 0.111469, 0.0], [1.680398, -0.373741, -0.758561], [1.680398, -0.373741, 0.758561]], 3, -5.0203),
    ("methane", &[6, 1, 1, 1, 1, 6, 1, 1, 1, 1], &[[0.0, -0.00014, 1.859161], [-0.888551, 0.51306, 1.494685], [0.888551, 0.51306, 1.494685], [0.0, -1.026339, 1.494868], [0.0, 8.9e-05, 2.948284], [0.0, 0.00014, -1.859161], [0.0, -8.9e-05, -2.948284], [-0.888551, -0.51306, -1.494685], [0.888551, -0.51306, -1.494685], [0.0, 1.026339, -1.494868]], 5, -0.5304),
];

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (name, z, pos, na, reference) = *DIMERS.iter().find(|d| Some(d.0) == args.first().map(|s| s.as_str())).expect("water, ammonia or methane");
    let scales: Vec<f64> = if args.len() > 1 { args[1..].iter().map(|s| s.parse().expect("a scale")).collect() } else { vec![1.0] };
    let a = 1.0 / 0.529177210903;
    let k = 627.509474;
    println!("# {name} dimer (S22; CCSD(T)/CBS {reference:.3} kcal/mol at scale 1): scale, heavy-atom distance (A), counterpoise HF, MP2 correlation, total; uncorrected HF, MP2 correlation, total; seconds");
    for scale in scales {
        let t0 = std::time::Instant::now();
        let axis: Vec<f64> = (0..3).map(|i| pos[na][i] - pos[0][i]).collect();
        let d0 = axis.iter().map(|x| x * x).sum::<f64>().sqrt();
        let p: Vec<[f64; 3]> = pos.iter().enumerate().map(|(i, r)| {
            let shift = if i >= na { scale - 1.0 } else { 0.0 };
            [(r[0] + shift * axis[0]) * a, (r[1] + shift * axis[1]) * a, (r[2] + shift * axis[2]) * a]
        }).collect();
        let ne = |r: std::ops::Range<usize>| -> f64 { z[r].iter().map(|&x| x as f64).sum() };
        let nat = z.len();
        let energy = |problem: &phys::electrons::scf::Problem, zs: &[u32], keep: &[usize]| -> (f64, f64) {
            let hf = hartree_fock(problem, 150, 1e-9);
            assert!(hf.converged, "{name}: Hartree-Fock did not converge");
            let c = mp2(problem, &hf, frozen_core(zs, keep));
            (hf.energy, c.correlation)
        };
        let dimer = Molecule { z: z.to_vec(), positions: p.clone(), charge: 0, unpaired: 0 }.problem(Functional::Pbe);
        let all: Vec<usize> = (0..nat).collect();
        let (ka, kb): (Vec<usize>, Vec<usize>) = ((0..na).collect(), (na..nat).collect());
        let ab = energy(&dimer, z, &all);
        let ga = energy(&dimer.with_ghosts(&ka, ne(0..na)), z, &ka);
        let gb = energy(&dimer.with_ghosts(&kb, ne(na..nat)), z, &kb);
        let alone = |r: std::ops::Range<usize>| -> (f64, f64) {
            let zs = z[r.clone()].to_vec();
            let keep: Vec<usize> = (0..zs.len()).collect();
            let m = Molecule { z: zs.clone(), positions: p[r].to_vec(), charge: 0, unpaired: 0 }.problem(Functional::Pbe);
            energy(&m, &zs, &keep)
        };
        let (ma, mb) = (alone(0..na), alone(na..nat));
        let (cp_hf, cp_c) = ((ab.0 - ga.0 - gb.0) * k, (ab.1 - ga.1 - gb.1) * k);
        let (un_hf, un_c) = ((ab.0 - ma.0 - mb.0) * k, (ab.1 - ma.1 - mb.1) * k);
        println!("{scale:.3} {:.4} {cp_hf:.4} {cp_c:.4} {:.4} {un_hf:.4} {un_c:.4} {:.4} {:.0}", d0 * scale, cp_hf + cp_c, un_hf + un_c, t0.elapsed().as_secs_f64());
    }
}
