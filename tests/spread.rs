//! How far a node's contents actually extend.
//!
//! A node's radius is fixed when the node is created and nothing ever checked
//! that it still describes what the node holds. `docs/PLAY.md` Phase 1 asks for
//! the measurement on its own, before any of the three things that need it —
//! node splitting, D6's patch handoff, and promoting a fragment that has left
//! its parent — because it is to be built once rather than three times.
//!
//! It is also a detector, and it earns that immediately: run against worlds
//! that already exist it reproduces two defects `BACKLOG.md` records, from a
//! completely different direction. See
//! `it_finds_the_two_defects_already_on_the_list`.

use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::{v3, Vec3};
use phys::state::Spread;
use phys::units::Tier;

/// Six points on the axes at `+/-d`, each of radius `r` and mass `m`.
fn shell(d: f64, r: f64, m: f64) -> Vec<(Vec3, f64, f64)> {
    [
        v3(d, 0.0, 0.0), v3(-d, 0.0, 0.0),
        v3(0.0, d, 0.0), v3(0.0, -d, 0.0),
        v3(0.0, 0.0, d), v3(0.0, 0.0, -d),
    ]
    .into_iter()
    .map(|p| (p, m, r))
    .collect()
}

/// The arithmetic, against a configuration whose answer is known by hand.
#[test]
fn a_symmetric_shell_measures_its_own_geometry() {
    let s = Spread::of(shell(4.0, 0.5, 3.0));
    assert_eq!(s.count, 6);
    assert!(s.centre.norm() < 1e-12, "a symmetric shell is centred on the origin, not {:?}", s.centre);
    assert!((s.rms - 4.0).abs() < 1e-12, "rms of a shell at 4 m is 4 m, not {}", s.rms);
    // The furthest *surface*, not the furthest centre: a body is not inside a
    // volume its own bulk sticks out of.
    assert!(
        (s.furthest - 4.5).abs() < 1e-12,
        "the furthest surface of a 0.5 m body centred at 4 m is at 4.5 m, not {}",
        s.furthest
    );
    assert!((s.occupancy(4.5) - 1.0).abs() < 1e-12);
    assert!((s.occupancy(9.0) - 0.5).abs() < 1e-12);
}

/// The centre is found, not assumed. A node whose contents have drifted off its
/// own origin is one of the states worth noticing, and measuring the spread
/// about the origin would report the drift as size.
#[test]
fn the_centre_is_measured_rather_than_assumed() {
    let offset = v3(100.0, 0.0, 0.0);
    let moved: Vec<_> = shell(4.0, 0.5, 3.0)
        .into_iter()
        .map(|(p, m, r)| (p + offset, m, r))
        .collect();
    let s = Spread::of(moved);
    assert!((s.centre - offset).norm() < 1e-9, "centre came out {:?}", s.centre);
    assert!((s.rms - 4.0).abs() < 1e-9, "a shell that moved is not a bigger shell: {}", s.rms);
}

/// And it is the centre of *mass*. One heavy body pulls it.
#[test]
fn the_centre_is_mass_weighted() {
    let parts = vec![(v3(-1.0, 0.0, 0.0), 9.0, 0.0), (v3(9.0, 0.0, 0.0), 1.0, 0.0)];
    let s = Spread::of(parts);
    assert!((s.centre.x - 0.0).abs() < 1e-12, "centre of mass is at the origin, not {:?}", s.centre);
    // rms^2 = (9*1 + 1*81)/10 = 9, so rms = 3.
    assert!((s.rms - 3.0).abs() < 1e-12, "rms came out {}", s.rms);
}

/// Massless contents still have a geometry.
///
/// Answering with the origin would report a spread about a point nothing is
/// near, which is worse than saying nothing.
#[test]
fn contents_with_no_mass_still_have_an_extent() {
    let parts: Vec<_> = shell(4.0, 0.0, 0.0)
        .into_iter()
        .map(|(p, _, r)| (p + v3(50.0, 0.0, 0.0), 0.0, r))
        .collect();
    let s = Spread::of(parts);
    assert!((s.centre.x - 50.0).abs() < 1e-9, "centre came out {:?}", s.centre);
    assert!((s.rms - 4.0).abs() < 1e-9, "rms came out {}", s.rms);
}

/// Nothing has no extent, and says so rather than inventing one.
#[test]
fn empty_contents_have_no_spread() {
    let s = Spread::of(Vec::new());
    assert_eq!(s.count, 0);
    assert_eq!(s.rms, 0.0);
    assert_eq!(s.furthest, 0.0);
    assert_eq!(s.occupancy(10.0), 0.0);

    let w = World::new(galaxy(0xE0E0, 1e9), 20.0);
    assert_eq!(w.tree.spread(w.tree.root).count, 0, "an unmaterialised node holds nothing");
    assert_eq!(w.tree.spread(NodeIdx::NONE).count, 0);
}

/// A promoted child is contents at a finer resolution, not a second thing.
///
/// `children` runs parallel to `bodies`, so a slot holding a promoted child
/// must be counted once. Counting both would put the stand-in and the child it
/// stands for at slightly different places and report the difference as spread.
#[test]
fn a_promoted_child_is_counted_once() {
    let mut w = World::new(galaxy(0xADDA, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 64;
    let root = w.tree.root;
    w.tree.refine(root);
    let before = w.tree.spread(root);
    assert_eq!(before.count, 64);

    let tier = w.tree.nodes[root.get()].tier;
    let child = w.tree.promote(root, 3, default_spec(tier.finer()));
    assert!(!child.is_none());
    // Move the child well away from where its stand-in sits, so counting both
    // would be visible rather than a rounding difference.
    w.tree.nodes[child.get()].motion.offset = v3(1.0e21, 0.0, 0.0);

    let after = w.tree.spread(root);
    assert_eq!(
        after.count, 64,
        "promotion changed the resolution of one slot, not the number of things"
    );
    assert!(
        after.furthest >= 1.0e21,
        "the child's own position should be what its slot contributes, and the \
         furthest came out {}",
        after.furthest
    );
}

/// The measurement, run against worlds that already exist, independently found
/// two things `BACKLOG.md` records. One of them has since been fixed, and this
/// is now what holds the fix.
///
/// This is the case for having the measurement. Neither defect was looked for
/// here — both fall out of asking one question of every node, and both were
/// originally found by long and unrelated investigations.
///
/// Measured baselines, which matter to anyone who later wants a threshold: a
/// *healthy* node does **not** come out at or below one. A Plummer sphere's
/// tail legitimately reaches three to four radii, and the scenario shelf runs
/// 1.2 to 4.0. So "outgrew its radius" is not the signal; orders of magnitude
/// are.
#[test]
fn it_finds_the_two_defects_already_on_the_list() {
    // 1. **Was:** the sampler inflates anything bound by chemistry by about
    //    4.3e5, because its relaxation loop releases a *gravitational* binding
    //    and a granite block's is cohesive. Three scenarios came out here at
    //    over 1e5 and the other five at 3.06 to 3.98.
    //
    //    Fixed by `PLAY.md` §7's second Phase 2 item, which splits
    //    `binding_energy` into the half expansion releases and the half it does
    //    not. Now every scenario sits inside the healthy band, and this half of
    //    the test is what keeps it there — the assertion is inverted rather
    //    than deleted, because a measurement that caught something once is the
    //    cheapest guard against it coming back.
    let mut occupancies = Vec::new();
    for sc in phys::scenario::ALL {
        let mut w = World::new(sc.build(0xC0FFEE), 1.0);
        let root = w.tree.root;
        w.tree.nodes[0].spec.count = 256;
        w.tree.refine(root);
        if w.tree.nodes[root.get()].bodies.is_empty() {
            continue;
        }
        let occ = w.tree.spread(root).occupancy(w.tree.nodes[root.get()].matter.radius);
        occupancies.push((sc.name, occ));
    }
    assert!(
        occupancies.len() >= 8,
        "the whole shelf should refine: {occupancies:?}"
    );
    assert!(
        occupancies.iter().all(|(_, o)| *o < 10.0),
        "a sampled profile's tail reaches a few radii and no more, whatever \
         holds the node together: {occupancies:?}"
    );
    // And the four chemically bound ones are not merely under the ceiling but
    // *tighter* than the gravitationally bound ones, which is what a solid
    // ought to look like beside a Plummer sphere.
    let chemical: Vec<f64> = occupancies
        .iter()
        .filter(|(n, _)| {
            matches!(*n, "Granite block" | "Water vapour" | "Carbon atom" | "Iron nucleus")
        })
        .map(|(_, o)| *o)
        .collect();
    assert_eq!(chemical.len(), 4, "the four bound-by-chemistry scenarios: {occupancies:?}");
    assert!(
        chemical.iter().all(|o| *o < 2.0),
        "a solid's contents sit inside it: {occupancies:?}"
    );

    // 2. A node flung out of its parent by an unstable solver. The point of
    //    the detector is *where* it fires: in the node that caused it, on the
    //    frame it happened, rather than twenty tiers away inside a hash.
    //
    //    **Which detector reports it moved in Phase 3, and the reason is worth
    //    keeping.** `worst_occupancy` saw this fault only because nothing
    //    re-homed the victim: a promoted child drifted out of a node claiming a
    //    metre, stayed its child, and the ratio climbed without limit — 10^6
    //    here and 10^22 for a nucleus. `docs/PLAY.md` D16's crossing pass
    //    re-homes it on the frame it leaves, which fixes the tree and would
    //    have made the fault invisible. So the measurement moved to the
    //    crossing itself: a node that *steps* over a boundary crosses at a
    //    ratio a hair above one, and a node that is *flung* over it arrives
    //    with a number that says so. Measured here: 3.2x10^13 times what the
    //    parent's contents reach, on the first frame.
    let mut w = World::new(galaxy(0xABCD, 1e9), 20.0);
    let root = w.tree.root;
    // Paced to its subject. A world runs at one second per second (`PLAY.md`
    // D1), and the fault being detected is a solver flinging bodies out over a
    // span only a galactic pace supplies.
    w.pace_to(root);
    let path = w.drill_to(root, Tier::Continuum.max_radius(), &default_spec);
    assert!(path.len() >= 6, "expected a deep ladder, got {}", path.len());
    for _ in 0..20 {
        w.step_frame(50_000.0);
    }
    println!(
        "  the ladder: {} crossings, worst at {:.3e} times what the parent holds",
        w.stats.crossings, w.stats.worst_crossing
    );
    assert!(
        w.stats.crossings > 0,
        "the ladder flings nodes clean out of their parents and nothing crossed"
    );
    assert!(
        w.stats.worst_crossing > 1.0e6,
        "the ladder flings a node far outside its parent and the detector \
         reported only {:e}",
        w.stats.worst_crossing
    );
    assert!(
        w.stats.worst_crossing_at.is_some(),
        "a number worth chasing has to say where to look"
    );
}
