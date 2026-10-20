use phys::electrons::functional::Functional;
use phys::electrons::hf::{frozen_core, hartree_fock, mp2};
use phys::electrons::molecule::Molecule;
use std::time::Instant;

fn run(more: Option<Vec<Vec<(usize, f64)>>>) -> (usize, f64, f64) {
    let a = 1.0 / 0.529177210903;
    // water at its S22 monomer geometry
    let pos: Vec<[f64; 3]> = [[-1.551007, -0.114520, 0.0], [-1.934259, 0.762503, 0.0], [-0.599677, 0.040712, 0.0]].iter().map(|p| [p[0] * a, p[1] * a, p[2] * a]).collect();
    let mol = Molecule { z: vec![8, 1, 1], positions: pos, charge: 0, unpaired: 0 };
    let p = mol.problem_plus(Functional::Pbe, None, more.as_deref());
    let hf = hartree_fock(&p, 100, 1e-9);
    let c = mp2(&p, &hf, frozen_core(&mol.z, &[0, 1, 2]));
    (p.basis.size, hf.energy, c.correlation)
}

#[test]
fn scan_higher_l() {
    let t = Instant::now();
    let (n0, hf0, c0) = run(None);
    println!("BASE {n0} functions: HF {hf0:.7} MP2 corr {c0:.7} ({:.0} s)", t.elapsed().as_secs_f64());
    for e in [0.4, 0.7, 1.1, 1.6, 2.4, 3.5] {
        let (n, hf, c) = run(Some(vec![vec![(4, e)], vec![], vec![]]));
        println!("g on O at {e}: {n} functions, HF {:+.2e}, MP2 corr {:+.3e} Ha", hf - hf0, c - c0);
    }
    for e in [0.3, 0.6, 1.0, 1.6, 2.5, 4.0] {
        let (n, hf, c) = run(Some(vec![vec![], vec![(3, e)], vec![(3, e)]]));
        println!("f on both H at {e}: {n} functions, HF {:+.2e}, MP2 corr {:+.3e} Ha", hf - hf0, c - c0);
    }
}
