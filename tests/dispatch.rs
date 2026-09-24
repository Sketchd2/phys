//! Dispatch reads state, not only size.
//!
//! `docs/PLAY.md` §3.3. `solvers::for_tier(Continuum)` is `Hydro`, and a
//! building, a wolf and a boulder are all `Continuum`. None of them is a fluid.
//!
//! A node holds both kinds at once — `sample_structured` lays a structure's
//! members out first and then "the unstructured remainder: litter, air,
//! rubble" — and until §3.3 the tier solver was handed the lot. The rule is
//! that **the tier says which regime the disordered contents are in, and the
//! node's own state says which contents are ordered**, in one node, in one
//! pass.
//!
//! Measured before the fix, on a forty-year-old tree standing on a planet and
//! advanced for one twentieth of a second: its members reached 1.9x10^8 m/s —
//! 64% of the speed of light — and travelled 9.6x10^6 m. The tree is 6.3 m
//! across.

use phys::engine::{default_spec, galaxy, World};
use phys::math::v3;
use phys::morph::{Environment, Program};
use phys::state::{Composition, Matter};
use phys::units::{Tier, YEAR};

const EARTH_MASS: f64 = 5.972e24;
const EARTH_RADIUS: f64 = 6.371e6;

/// A grown tree standing on a planet, materialised.
///
/// On a planet because a structure's weight is now derived from what it sits
/// inside (`PLAY.md` D6); materialised because an unmaterialised node has no
/// bodies for anything to dispatch on.
fn a_tree_on_a_planet() -> (World, phys::ids::NodeIdx) {
    let mut w = World::new(galaxy(0x3333, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let planet = w.tree.promote(root, 7, default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.spec.count = 8;
    }
    w.tree.refine(planet);
    let node = w.tree.promote(planet, 0, default_spec(Tier::Continuum));
    {
        let n = &mut w.tree.nodes[node.get()];
        n.matter = Matter::neutral(900.0, 6.0, 291.0, Composition::organic());
        n.spec.count = 600;
        n.motion.offset = v3(0.0, 0.0, EARTH_RADIUS);
    }
    w.plant(node, Program::Tree, Some(Environment::default()));
    for _ in 0..40 {
        w.grow_node(node, YEAR);
    }
    w.tree.refine(node);
    (w, node)
}

/// A node holds both, and the mask says which is which.
#[test]
fn a_node_holds_ordered_and_disordered_contents_at_once() {
    let (w, node) = a_tree_on_a_planet();
    let n = &w.tree.nodes[node.get()];
    assert_eq!(n.tier, Tier::Continuum, "a metre-scale tree is Continuum");
    assert_eq!(
        phys::solvers::for_tier(n.tier),
        phys::solvers::SolverKind::Hydro,
        "and Continuum dispatches to a fluid solver, which is the whole problem"
    );

    let mask = n.structural_mask().expect("a planted tree has ordered contents");
    let ordered = mask.iter().filter(|m| **m).count();
    let loose = mask.len() - ordered;
    assert_eq!(mask.len(), n.bodies.len(), "the mask is parallel to the bodies");
    assert!(ordered > 0 && loose > 0, "{ordered} ordered, {loose} loose — the test needs both");
    assert_eq!(
        ordered,
        n.last_report.structural_parts,
        "the mask should agree with what the sampler reported it built"
    );
}

/// The mask comes from the topology, which is persisted, and not from the
/// sample report, which is not.
///
/// A node reloaded from a file has to dispatch the same way as the one that
/// wrote it. `SampleReport` is documented as a diagnostic of the last
/// materialisation and is deliberately regenerated rather than stored, so
/// reading `structural_parts` would have made dispatch change across a save.
/// The node is **pinned** on purpose, and that is what makes this a test.
///
/// An unpinned node comes back with no bodies, and re-materialising it
/// regenerates the sample report too — so a report-based mask would agree and
/// the test would prove nothing. It was written that way first and passed with
/// the mask read from `structural_parts`. Pinned, the bodies come back from the
/// file and the report does not, which is the state a reloaded world is
/// actually in.
#[test]
fn the_mask_survives_what_the_sample_report_does_not() {
    let (mut w, node) = a_tree_on_a_planet();
    w.tree.pin(node);
    let before = w.tree.nodes[node.get()].structural_mask().unwrap();
    let ordered = before.iter().filter(|m| **m).count();

    let bytes = phys::persist::encode(w.view());
    let back = World::from_snapshot(phys::persist::decode(&bytes).expect("decode"), 20.0);

    let reloaded = &back.tree.nodes[node.get()];
    assert_eq!(
        reloaded.last_report.structural_parts, 0,
        "the sample report is not persisted, which is the point of this test"
    );
    assert_eq!(
        reloaded.bodies.len(),
        before.len(),
        "a pinned node's detail is written, so it should come back with it"
    );

    let after = reloaded.structural_mask().expect("a reloaded structure is still a structure");
    assert_eq!(
        after.iter().filter(|m| **m).count(),
        ordered,
        "a reloaded node dispatched differently from the one that wrote it"
    );
    assert_eq!(after, before, "and not merely the same count");
}

/// The headline: a structure's members are not the tier solver's to move.
#[test]
fn a_structures_members_are_not_handed_to_the_tier_solver() {
    let (mut w, node) = a_tree_on_a_planet();
    let mask = w.tree.nodes[node.get()].structural_mask().unwrap();
    let before: Vec<_> = w.tree.nodes[node.get()].bodies.iter().map(|b| b.pos).collect();
    let extent = w.tree.nodes[node.get()].matter.radius;

    for _ in 0..5 {
        w.advance_node(node, 0.05);
    }

    let n = &w.tree.nodes[node.get()];
    assert_eq!(n.bodies.len(), before.len(), "the node lost its detail mid-test");
    let mut worst_member = 0.0f64;
    let mut worst_loose = 0.0f64;
    for (i, b) in n.bodies.iter().enumerate() {
        let moved = (b.pos - before[i]).norm();
        if mask[i] {
            worst_member = worst_member.max(moved);
        } else {
            worst_loose = worst_loose.max(moved);
        }
    }
    assert_eq!(
        worst_member, 0.0,
        "a standing tree's members moved {worst_member} m in a quarter of a \
         second; the tree is {extent} m across"
    );
    // And the disordered contents are still the tier solver's, or the rule has
    // been applied by switching the solver off rather than by reading state.
    assert!(
        worst_loose > 0.0,
        "nothing at all moved, so this would pass with the solver removed"
    );
}

/// A node with no ordered contents is untouched by any of it.
#[test]
fn a_node_that_is_not_a_structure_is_not_partitioned() {
    let mut w = World::new(galaxy(0x9A11, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 256;
    let root = w.tree.root;
    w.pace_to(root);
    w.tree.refine(root);
    assert!(
        w.tree.nodes[root.get()].structural_mask().is_none(),
        "a galaxy has no ordered contents"
    );
    let before: Vec<_> = w.tree.nodes[root.get()].bodies.iter().map(|b| b.pos).collect();
    w.advance_node(root, w.pace);
    let moved = w.tree.nodes[root.get()]
        .bodies
        .iter()
        .enumerate()
        .map(|(i, b)| (b.pos - before[i]).norm())
        .fold(0.0f64, f64::max);
    assert!(moved > 0.0, "a galaxy's stars should still move: {moved}");
}

/// A fluid is substepped to what its own signal speed allows.
///
/// The scheduler's timestep answers to causality and to the tier; a pressure
/// wave answers to neither. `hydro::courant_dt` had been in the file the whole
/// time and nothing called it, while the identical problem in molecular
/// dynamics was handled three lines away. Measured on the loose contents of the
/// tree above: the stable step is 1.7x10^-5 s and a frame asked for 0.05 —
/// three thousand times over.
#[test]
fn a_fluid_is_substepped_to_its_courant_limit() {
    let (mut w, node) = a_tree_on_a_planet();
    let asked = 0.05;
    let report = w.advance_node(node, asked);
    assert!(
        report.dt_used > 0.0 && report.dt_used < asked,
        "a node that cannot integrate the whole span stably should cover the \
         part it can and say so: dt_used {} against {asked}",
        report.dt_used
    );
    assert!(
        report.steps > 1,
        "one step for a span three thousand times the stable one: {} steps",
        report.steps
    );
    // And the node's clock follows what was integrated, not what was asked,
    // so the shortfall becomes lateness the scheduler can see.
    assert!(
        w.tree.nodes[node.get()].time <= asked,
        "the clock ran past what was actually integrated"
    );
}

/// The equation of state says when it is being asked about something it does
/// not describe.
///
/// `docs/PLAY.md` §7's ninth Phase 2 item, on §3.7's own precedent: a node
/// crossed by its ensemble reports it rather than doing it quietly.
///
/// `Matter::pressure` is an ideal gas plus radiation and is the only equation
/// of state the engine has. There are two ways to be outside it. A node that
/// **is not a gas** — which D17 made answerable, because a node carrying a
/// mixture knows its own phase — and a node that is **mostly vacuum with solids
/// in it**, whose every body is the stand-in for a promoted child. The second
/// is what `docs/BACKLOG.md` measured detonating: a 12 m root holding a
/// 48-tonne box and a 10 kg ball, with zero collisions, took the ball from
/// 5.36 m/s to 566 over six frames.
///
/// This is Phase 2 making it **visible**. Water fixes it, with the liquid and
/// solid equation of state its second piece names.
#[test]
fn a_gas_law_asked_about_a_solid_says_so() {
    use phys::chem::{Mixture, Phase};
    use phys::engine::World;
    use phys::sampler::{MassSpectrum, Profile, SampleSpec};
    use phys::state::{BodyKind, Composition, Matter};
    use phys::tree::Tree;
    use phys::units::Tier;

    let build = |described: bool| {
        let spec = SampleSpec::new(32, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
        let matter = Matter::neutral(2650.0, 0.5, 290.0, Composition::primordial());
        let mut w = World::new(Tree::new(0xE05, matter, Tier::Continuum, spec), 20.0);
        let root = w.tree.root;
        if described {
            let silica = w
                .substances
                .intern(phys::material::substances::silica_arrangement())
                .expect("silica analyses");
            let mut mix = Mixture::new();
            mix.add(silica, Phase::Solid, 1.0);
            w.set_mixture(root, mix);
        }
        w.tree.refine(root);
        (w, root)
    };

    // Undescribed matter is *not* reported, and that is deliberate: "no
    // information" is not the same answer as "measured and wrong", and almost
    // every node in a galaxy genuinely is a gas.
    let (mut quiet, root) = build(false);
    for _ in 0..4 {
        quiet.advance_node(root, 1.0e-4);
    }
    println!("  undescribed: {} reports", quiet.stats.eos_outside_validity);
    assert_eq!(
        quiet.stats.eos_outside_validity, 0,
        "matter nobody has described must not be accused of being the wrong phase"
    );

    // The same node, told what it is made of. Nothing else changes.
    let (mut loud, root) = build(true);
    assert!(
        !loud.tree.nodes[root.get()].matter.gas_law_applies(),
        "a node of solid silicate is not a gas"
    );
    for _ in 0..4 {
        loud.advance_node(root, 1.0e-4);
    }
    println!(
        "  described as solid silicate: {} reports, at {:?}",
        loud.stats.eos_outside_validity, loud.stats.eos_outside_validity_at
    );
    assert!(
        loud.stats.eos_outside_validity > 0,
        "a solid priced through the gas law has to say so"
    );
    assert!(
        loud.stats.eos_outside_validity_at.is_some(),
        "a number worth chasing has to say where to look"
    );

    // And the other way to be outside it: a node that is **mostly vacuum with
    // solids in it**. Nothing here is described at all, so the phase reading
    // says nothing; what says something is that every body in the node is the
    // stand-in for a promoted child, so there are no contents of its own for a
    // fluid solver to be about.
    let spec = SampleSpec::new(2, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let matter = Matter::neutral(48_010.0, 12.0, 290.0, Composition::primordial());
    let mut w = World::new(Tree::new(0xB0F, matter, Tier::Continuum, spec), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let small = SampleSpec::new(2, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    w.tree.promote(root, 0, small);
    w.tree.promote(root, 1, small);
    assert!(
        w.tree.nodes[root.get()].matter.gas_law_applies(),
        "nothing has described it, so the phase reading has nothing to say"
    );
    let before = w.stats.eos_outside_validity;
    for _ in 0..4 {
        w.advance_node(root, 1.0e-4);
    }
    println!(
        "  every body a stand-in: {} reports",
        w.stats.eos_outside_validity - before
    );
    assert!(
        w.stats.eos_outside_validity > before,
        "a Continuum node that is two promoted solids and vacuum is being priced \
         as a hot dense gas, and has to say so"
    );
}
