//! Ion-water pairs through the unlike-pair code (`docs/SOLUTIONS.md` S5): a sodium ion on the water
//! molecule's axis and a chloride ion on a straight O-H bond, each at a bond length typical of the
//! gas-phase complex, with the counterpoise MP2 interaction printed against the ranges the literature
//! gives (the author's recollection: Na+...OH2 about -24 kcal/mol at the CCSD(T) basis limit,
//! Cl-...HOH about -14.5 to -15; to be looked up and cited before they decide anything). This asserts
//! only that the energies are bound and of the right order, so a wrong charge or electron count, which
//! gives positive or wildly large numbers, fails it.

use phys::math::Vec3;

const ANGSTROM: f64 = 1.0 / 0.529177210903;
const KCAL: f64 = 627.509474;

#[test]
fn ions_bind_a_water_molecule_with_about_the_right_energy() {
    std::env::set_current_dir("..").expect("the repo root");
    if phys_gpu::install_fock().is_err() {
        eprintln!("no GPU Fock; the CPU does it (slow)");
    }
    let water = phys::pairs::Monomer::load("water", false);
    let text = std::fs::read_to_string("grow-water.state").expect("a state");
    let state = phys::electrons::grow::Resume::from_text(&text).expect("readable");
    let p: Vec<Vec3> = state.positions.iter().map(|v| Vec3 { x: v[0], y: v[1], z: v[2] }).collect();
    let (o, h1, h2) = (p[0], p[1], p[2]);
    let axis = (o.scale(2.0) - h1 - h2).unit();
    let na = phys::pairs::Monomer::load("sodium", false);
    let cl = phys::pairs::Monomer::load("chloride", false);
    assert_eq!((na.charge, cl.charge), (1, -1));
    let mut report = Vec::new();
    for (ion, pos, lo, hi, label) in [(&na, o + axis.scale(2.2 * ANGSTROM), -34.0, -18.0, "Na+ on the axis, Na-O 2.2 A"), (&cl, h1 + (h1 - o).unit().scale(2.2 * ANGSTROM), -22.0, -8.0, "Cl- on an O-H bond, Cl...H 2.2 A")] {
        let (total, hf) = phys::pairs::unlike_interaction(&water, ion, &p, &[pos]);
        println!("{label}: MP2 {:.2} kcal/mol (Hartree-Fock {:.2})", total * KCAL, hf * KCAL);
        report.push((label, total * KCAL, lo, hi));
    }
    for (label, e, lo, hi) in report {
        assert!(e > lo && e < hi, "{label}: {e:.2} kcal/mol is outside {lo} to {hi}");
    }
}
