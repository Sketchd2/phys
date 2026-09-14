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

/// The measurement, run against worlds that already exist, independently finds
/// two things `BACKLOG.md` records.
///
/// This is the case for having it. Neither defect was looked for here — both
/// fall out of asking one question of every node, and both were originally
/// found by long and unrelated investigations.
///
/// Measured baselines, which matter to anyone who later wants a threshold: a
/// *healthy* node does **not** come out at or below one. A Plummer sphere's
/// tail legitimately reaches three to four radii, and the scenario shelf runs
/// 1.5 to 4.0. So "outgrew its radius" is not the signal; orders of magnitude
/// are.
#[test]
fn it_finds_the_two_defects_already_on_the_list() {
    // 1. The sampler inflates anything bound by chemistry by about 4.3e5,
    //    because its relaxation loop releases a *gravitational* binding and a
    //    granite block's is cohesive.
    let mut healthy = Vec::new();
    let mut inflated = Vec::new();
    for sc in phys::scenario::ALL {
        let mut w = World::new(sc.build(0xC0FFEE), 1.0);
        let root = w.tree.root;
        w.tree.nodes[0].spec.count = 256;
        w.tree.refine(root);
        if w.tree.nodes[root.get()].bodies.is_empty() {
            continue;
        }
        let occ = w.tree.spread(root).occupancy(w.tree.nodes[root.get()].matter.radius);
        if occ > 100.0 { &mut inflated } else { &mut healthy }.push((sc.name, occ));
    }
    assert!(
        healthy.iter().all(|(_, o)| *o < 10.0),
        "a sampled profile's tail reaches a few radii and no more: {healthy:?}"
    );
    assert_eq!(
        inflated.len(),
        3,
        "three scenarios set a chemical binding energy and should be the three \
         that come out inflated; got {inflated:?} against {healthy:?}"
    );
    assert!(
        inflated.iter().all(|(_, o)| *o > 1.0e5),
        "the inflation is 1.5^32, about 4.3e5: {inflated:?}"
    );

    // 2. A node whose bodies are flung out of it by an unstable solver. The
    //    point of the detector is *where* it fires: in the node that caused it,
    //    on the frame it happened, rather than twenty tiers away inside a hash.
    let mut w = World::new(galaxy(0xABCD, 1e9), 20.0);
    let root = w.tree.root;
    let path = w.drill_to(root, Tier::Continuum.max_radius(), &default_spec);
    assert!(path.len() >= 6, "expected a deep ladder, got {}", path.len());
    for _ in 0..20 {
        w.step_frame(50_000.0);
    }
    assert!(
        w.stats.worst_occupancy > 1.0e6,
        "the ladder flings a node's bodies far outside it and the detector \
         reported only {:e}",
        w.stats.worst_occupancy
    );
    assert!(
        w.stats.worst_occupancy_at.is_some(),
        "a number worth chasing has to say where to look"
    );
}
