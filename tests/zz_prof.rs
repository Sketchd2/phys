use phys::electrons::functional::Functional;
use phys::electrons::hf::{frozen_core, hartree_fock, mp2};
use phys::electrons::molecule::Molecule;
use std::time::Instant;

#[test]
fn profile_trimer() {
    let text = std::fs::read_to_string("cluster-water-3-c0.txt").unwrap();
    let mut pos: Vec<[f64; 3]> = Vec::new();
    for l in text.lines().filter(|l| l.starts_with("# mol ")) {
        let w: Vec<f64> = l.split_whitespace().skip(3).map(|x| x.parse().unwrap()).collect();
        for c in w.chunks(3) {
            pos.push([c[0], c[1], c[2]]);
        }
    }
    let mol = Molecule { z: vec![8, 1, 1, 8, 1, 1, 8, 1, 1], positions: pos, charge: 0, unpaired: 0 };
    let base = mol.problem(Functional::Pbe);
    println!("PROFILE trimer: {} functions", base.basis.size);
    let t = Instant::now();
    let hf = hartree_fock(&base, 100, 1e-9);
    let th = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let all: Vec<usize> = (0..9).collect();
    let c = mp2(&base, &hf, frozen_core(&mol.z, &all));
    println!("PROFILE trimer: HF {:.8} in {th:.1} s ({} iterations), MP2 {:.8} in {:.1} s", hf.energy, hf.iterations, c.correlation, t.elapsed().as_secs_f64());
}
