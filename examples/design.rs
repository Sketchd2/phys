//! What the on-creation design pass does to a generated structure.
use phys::morph::{Morphology, Program};
use phys::prolong::prolong_structured;
use phys::state::Aggregate;
use phys::units::YEAR;

fn main() {
    for (label, prog, mass, budget) in [
        ("tree, 900 kg", Program::Tree, 900.0, 400usize),
        ("tree, 6 t", Program::Tree, 6000.0, 400),
        ("tree, 6 t, coarse", Program::Tree, 6000.0, 120),
        ("tower, 3000 t", Program::Tower, 3.0e6, 600),
    ] {
        let mut m = if prog.is_planned() {
            let mut m = Morphology::planned(prog, mass, 11, 0x77);
            m.progress = 1.0;
            m
        } else {
            let mut m = Morphology::new(prog, 0xACE, 0x1234, 0);
            m.age = 60.0 * YEAR;
            m
        };
        m.built = mass;
        let agg = Aggregate::neutral(mass, m.extent(), 291.0, prog.substrate());
        let (_, _, report) = prolong_structured(&agg, &m, budget, 7, 0x1234, 0);
        let d = report.design;
        println!(
            "{label:<20} peak {:.3} -> {:.3}   spread {:.3} -> {:.3}   {} passes   \
             volume error {:.2e}",
            d.peak_before, d.peak_after, d.spread_before, d.spread_after, d.passes,
            d.volume_error()
        );
    }
}
