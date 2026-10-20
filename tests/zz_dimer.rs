use phys::electrons::functional::Functional;
use phys::electrons::hf::{frozen_core, hartree_fock, mp2};
use phys::electrons::molecule::Molecule;
use std::time::Instant;

fn dimer(label: &str, more: Option<Vec<Vec<(usize, f64)>>>) {
    let a = 1.0 / 0.529177210903;
    let pos: Vec<[f64; 3]> = [[-1.551007, -0.114520, 0.0], [-1.934259, 0.762503, 0.0], [-0.599677, 0.040712, 0.0], [1.350625, 0.111469, 0.0], [1.680398, -0.373741, -0.758561], [1.680398, -0.373741, 0.758561]].iter().map(|p| [p[0] * a, p[1] * a, p[2] * a]).collect();
    let mol = Molecule { z: vec![8, 1, 1, 8, 1, 1], positions: pos, charge: 0, unpaired: 0 };
    let base = mol.problem_plus(Functional::Pbe, None, more.as_deref());
    let t = Instant::now();
    let mut e = [(0.0, 0.0); 3];
    for (slot, (keep, ne)) in [(vec![0usize, 1, 2, 3, 4, 5], 20.0), (vec![0, 1, 2], 10.0), (vec![3, 4, 5], 10.0)].into_iter().enumerate() {
        let p = base.with_ghosts(&keep, ne);
        let hf = hartree_fock(&p, 100, 1e-9);
        let c = mp2(&p, &hf, frozen_core(&mol.z, &keep));
        e[slot] = (hf.energy, c.correlation);
    }
    let k = 627.509474;
    let hfi = (e[0].0 - e[1].0 - e[2].0) * k;
    let ci = (e[0].1 - e[1].1 - e[2].1) * k;
    println!("DIMER {label}: {} functions, HF {hfi:.3} + MP2 {ci:.3} = {:.3} kcal/mol (CCSD(T)/CBS -5.020); {:.0} s", base.basis.size, hfi + ci, t.elapsed().as_secs_f64());
}

#[test]
fn with_higher_l() {
    dimer("g(O) 2.0, f(H) 1.6", Some(vec![vec![(4, 2.0)], vec![(3, 1.6)], vec![(3, 1.6)], vec![(4, 2.0)], vec![(3, 1.6)], vec![(3, 1.6)]]));
}
