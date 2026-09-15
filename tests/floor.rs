//! The resolution floor is reported, not applied in silence.
//!
//! `docs/PLAY.md` §3.7: the floor "is real, is derivable, and the engine should
//! *report* it rather than silently drop a node to its ensemble — the same
//! discipline `displacement_ratio` already applies to the small-displacement
//! regime."
//!
//! §3.4 derives it. A node is followed while `frame_span <= node_dt *
//! MAX_SUBSTEPS`; inside `Continuum` the binding term is `0.25 * h / c_signal`;
//! so `h >= 4 * frame_span * c_signal / MAX_SUBSTEPS`, which at the 50 ms frame
//! a world runs at is `c_signal / 1280`.

use phys::engine::{default_spec, galaxy, World, MAX_SUBSTEPS};
use phys::units::Tier;

/// The floor reproduces the table §3.4 states, from the formula §3.4 derives.
///
/// This is the test that found the pace defect. It read 5.312 m for air against
/// a stated 0.27 — exactly twenty times coarser — because the world was running
/// at one second per *frame* rather than per second. The formula was right and
/// the clock was not.
#[test]
fn the_floor_is_the_one_the_plan_derives() {
    let w = World::new(galaxy(0xF100, 1e9), 20.0);
    assert!(
        (w.frame_dt() - 0.05).abs() < 1e-12,
        "§3.4's table is computed at a 50 ms frame, and this world's is {}",
        w.frame_dt()
    );
    for (medium, signal, expected) in [("air", 340.0, 0.27), ("water", 1500.0, 1.2), ("rock", 5000.0, 3.9)] {
        let floor = w.resolution_floor(signal);
        assert!(
            (floor - expected).abs() < 0.05,
            "{medium} at {signal} m/s should be followable to about {expected} m, \
             and came out {floor}"
        );
        // And it is `c_signal / 1280` at this frame, which is the form the plan
        // quotes it in.
        assert!((floor - signal / 1280.0).abs() < 1e-9);
    }
    assert_eq!(w.resolution_floor(0.0), 0.0, "no signal, no floor");
    assert_eq!(w.resolution_floor(f64::NAN), 0.0);
}

/// It is a property of the material, not of the tier.
///
/// "Where exactly it falls off depends on the material, not on the tier" —
/// §3.4. Rock is followable to fifteen times the size air is.
#[test]
fn the_floor_moves_with_the_material_and_not_the_tier() {
    let w = World::new(galaxy(0xF102, 1e9), 20.0);
    let air = w.resolution_floor(340.0);
    let rock = w.resolution_floor(5000.0);
    assert!(
        rock > air * 14.0,
        "a stiffer medium is followable only at a coarser size: {rock} against {air}"
    );
}

/// Halving the frame span halves the floor: detail and the clock trade against
/// each other, which is the whole of §3.4's argument.
#[test]
fn a_shorter_frame_buys_finer_detail() {
    let mut w = World::new(galaxy(0xF103, 1e9), 20.0);
    let at_fifty_ms = w.resolution_floor(340.0);
    w.pace_fixed(0.025);
    let at_twenty_five = w.resolution_floor(340.0);
    assert!(
        (at_twenty_five / at_fifty_ms - 0.5).abs() < 1e-9,
        "halving the span should halve the floor: {at_twenty_five} against {at_fifty_ms}"
    );
}

/// A node asked for more substeps than it can have, and the engine says so.
///
/// Reported rather than clamped: a node wanting 289 million substeps and one
/// wanting 257 are both "capped" and are not the same situation.
#[test]
fn what_a_node_wanted_is_reported_uncapped() {
    let mut w = World::new(galaxy(0xF104, 1e9), 20.0);
    let root = w.tree.root;
    let path = w.drill_to(root, Tier::Continuum.max_radius(), &default_spec);
    assert!(path.len() >= 6, "expected a deep ladder, got {}", path.len());
    for _ in 0..8 {
        w.step_frame(50_000.0);
    }
    assert!(
        w.stats.worst_substeps > 0.0,
        "nothing reported what it needed"
    );
    assert!(
        w.stats.worst_substeps_at.is_some(),
        "a number worth chasing has to say where to look"
    );
    // The two counters between them account for every node the floor caught:
    // one that could be dropped to its ensemble, one that could not.
    let caught = w.stats.ensembled + w.stats.unreachable;
    if w.stats.worst_substeps > MAX_SUBSTEPS as f64 {
        assert!(
            caught > 0,
            "a node wanted {} substeps, past the cap of {MAX_SUBSTEPS}, and \
             nothing was reported as having been caught by the floor",
            w.stats.worst_substeps
        );
    }
}

/// Dropping a node to its ensemble is counted.
///
/// This is the silent path §3.7 names. `Stats::unreachable` already counted the
/// nodes that could *not* be dropped — pinned, bubbled, or with a promoted
/// child, so somebody is deliberately watching them run — and the ones that
/// *were* dropped went unrecorded, which is the wrong way round: falling behind
/// shows up in lateness anyway, and being replaced by a draw from your own
/// equilibrium shows up in nothing.
///
/// A **leaf** on purpose. `forgettable` requires no promoted children, so every
/// node of a drilled ladder but the last one is unreachable rather than
/// ensembled — which is how the first version of this test passed with the new
/// counter deleted, by asserting on the sum of the two.
#[test]
fn a_node_crossed_by_its_ensemble_is_counted() {
    use phys::sampler::{MassSpectrum, Profile, SampleSpec};
    use phys::state::{BodyKind, Composition, Matter};
    use phys::tree::Tree;

    // A centimetre of matter: far below the floor its own signal speed sets, so
    // its trajectory cannot be followed at any affordable substep count.
    let matter = Matter::neutral(1.0, 0.01, 290.0, Composition::primordial());
    let spec = SampleSpec::new(512, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let mut w = World::new(Tree::new(0xF105, matter, Tier::Continuum, spec), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);

    assert!(w.forgettable(root), "a leaf with nothing promoted out of it is forgettable");
    assert!(!w.can_resolve(root), "and it is below the floor, which is the setup");
    let floor = w.resolution_floor_of(root);
    assert!(
        floor > w.tree.nodes[root.get()].matter.radius * 100.0,
        "the node should be far below its own floor: {floor} m against a \
         {} m node",
        w.tree.nodes[root.get()].matter.radius
    );

    for _ in 0..6 {
        w.step_frame(50_000.0);
    }

    assert!(
        w.stats.ensembled > 0,
        "a node 376 times below its own resolution floor was crossed by its \
         ensemble and nothing said so; worst substeps {:.3e}",
        w.stats.worst_substeps
    );
    assert_eq!(
        w.stats.unreachable, 0,
        "and it is not the other counter: this node could be dropped, and was"
    );
    assert!(
        w.stats.worst_substeps > MAX_SUBSTEPS as f64,
        "the demand should be past the cap: {:.3e}",
        w.stats.worst_substeps
    );
}
