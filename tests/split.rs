//! A node that has stopped being one region, and the inverse.
//!
//! `docs/BACKLOG.md`: "The tree can make a node's contents finer or coarser in
//! place, and it can push a child up a tier. It cannot say *these bodies are no
//! longer one neighbourhood* and hand them to two nodes." Phase 1 built the
//! measurement — `state::Spread` — and connected it to nothing. `docs/PLAY.md`
//! §7 Phase 3 is what needs it, because the same question drives a boundary
//! crossing one level up.
//!
//! The criterion is **connected components of the node's own contents, with the
//! split refused while one component holds the majority of the mass**. The
//! second half is not decoration: connectedness alone is a property of the
//! sample count rather than of the node, because the linking length is the mean
//! spacing and a centrally concentrated profile's outskirts are sparser than
//! the mean. Measured, on worlds nobody had touched:
//!
//! ```text
//! a rocky planet, 64 bodies      9 components, largest 86.5%
//! a rocky planet, 512 bodies    84 components, largest 77.8%
//! a planetary node, 4000        31 components, largest 99.05%, rest in pairs
//! a cloud pulled into two       20 components, largest 48.4%, next 47.7%
//! ```
//!
//! Only the last of those is two places.

use phys::engine::{default_spec, World};
use phys::ids::NodeIdx;
use phys::units::*;

/// A molecular cloud, materialised, with its draw pulled into two clumps a
/// radius and a half apart. Nothing else is touched: the bodies are the ones
/// the sampler made, moved.
fn a_cloud_in_two_minds(seed: u64) -> (World, NodeIdx, NodeIdx) {
    let sc = phys::scenario::ALL
        .iter()
        .find(|s| s.name == "Molecular cloud")
        .expect("the shelf has a molecular cloud");
    let mut w = World::new(sc.build(seed), 20.0);
    let root = w.tree.root;
    w.tree.nodes[root.get()].spec.count = 512;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    // A clump of the cloud, promoted so that the thing being split has a parent
    // to become a sibling in.
    let node = w.tree.promote(root, 0, default_spec(tier.finer()));
    assert!(!node.is_none());
    {
        let n = &mut w.tree.nodes[node.get()];
        n.spec.count = 512;
    }
    w.tree.refine(node);
    let r = w.tree.nodes[node.get()].matter.radius;
    // **Split by mass rather than by index.** The majority rule says a node is
    // one region while one component holds more than half of it, so a fixture
    // that halves the *count* of a draw with a mass spectrum lands wherever the
    // spectrum puts it — measured at 52.4 / 45.0 on this very cloud, which is a
    // majority and stays one node. That is the rule working as specified; it is
    // not what this test is about. Balancing the two sides by mass is what
    // makes this the two-places case rather than a test of where the boundary
    // falls.
    let mut order: Vec<usize> = (0..w.tree.nodes[node.get()].bodies.len()).collect();
    order.sort_by(|a, b| {
        let bodies = &w.tree.nodes[node.get()].bodies;
        bodies[*b].mass.total_cmp(&bodies[*a].mass)
    });
    let (mut left, mut right) = (0.0f64, 0.0f64);
    let mut side = vec![false; order.len()];
    for k in order {
        let m = w.tree.nodes[node.get()].bodies[k].mass;
        if left <= right {
            left += m;
        } else {
            side[k] = true;
            right += m;
        }
    }
    for (k, b) in w.tree.nodes[node.get()].bodies.iter_mut().enumerate() {
        let shift = if side[k] { 0.9 * r } else { -0.9 * r };
        b.pos = phys::math::v3(b.pos.x * 0.25 + shift, b.pos.y * 0.25, b.pos.z * 0.25);
    }
    // Somebody moved them, so the detail is theirs rather than the sampler's.
    w.tree.pin(node);
    (w, root, node)
}

/// The headline: a node whose contents have become two becomes two nodes.
#[test]
fn a_node_whose_contents_become_two_becomes_two_nodes() {
    let (mut w, root, node) = a_cloud_in_two_minds(0xBEEF);
    let before = w.conserved();
    let (_, _, components) = w.tree.components(node);
    let mass_before = w.tree.nodes[node.get()].matter.mass;
    let kids_before = w.tree.nodes[root.get()].children.iter().filter(|c| !c.is_none()).count();

    let new = w.tree.resolve_extent(node).expect("two clumps should become two nodes");

    let (a, b) = (&w.tree.nodes[node.get()], &w.tree.nodes[new.get()]);
    println!(
        "  {components} components -> two nodes: {} bodies at {:.3e} kg and {} at {:.3e} kg, \
         centres {:.3e} m apart",
        a.bodies.iter().filter(|x| x.mass > 0.0).count(),
        a.matter.mass,
        b.bodies.len(),
        b.matter.mass,
        (a.motion.offset - b.motion.offset).norm()
    );
    assert_eq!(w.tree.nodes[new.get()].parent, root, "the new node is a sibling, not a child");
    assert!(b.bodies.len() > 10, "the new node should hold the clump that left");
    assert!(
        (a.matter.mass + b.matter.mass - mass_before).abs() < mass_before * 1e-9,
        "the two halves do not add up to what there was"
    );
    assert_eq!(
        w.tree.nodes[root.get()].children.iter().filter(|c| !c.is_none()).count(),
        kids_before + 1,
        "the parent should hold one more thing than it did"
    );

    // Nothing is created or destroyed by dividing it.
    let after = w.conserved();
    let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(y.abs()).max(1e-300);
    println!(
        "  energy {:.6e} -> {:.6e} ({:.2e}), baryon {:.2e}",
        before.energy,
        after.energy,
        rel(before.energy, after.energy),
        rel(before.baryon, after.baryon)
    );
    assert!(rel(before.energy, after.energy) < 1e-9, "energy changed across a split");
    assert!(rel(before.baryon, after.baryon) < 1e-9, "baryon number changed across a split");
}

/// And the split settles: once it is two nodes, neither keeps dividing.
#[test]
fn a_split_settles() {
    let (mut w, _, node) = a_cloud_in_two_minds(0xBEEF);
    let first = w.tree.resolve_extent(node).expect("the first split");
    assert!(w.tree.resolve_extent(node).is_none(), "the node split twice");
    assert!(w.tree.resolve_extent(first).is_none(), "the piece that left split again");
    println!("  one split, then nothing: {} splits recorded", w.tree.stats.splits);
    assert_eq!(w.tree.stats.splits, 1);
}

/// A fresh draw is one region, at every sample count.
///
/// The regression that matters, and the one that cost the most: linking at the
/// mean spacing peels isolated bodies out of any centrally concentrated draw,
/// so without the majority rule a planet nobody has touched splits every few
/// frames for ever. Measured before the rule: three new nodes in ten frames.
#[test]
fn a_fresh_draw_is_one_region() {
    for count in [64usize, 512, 4000] {
        let sc = phys::scenario::ALL
            .iter()
            .find(|s| s.name == "Rocky planet")
            .expect("the shelf has a rocky planet");
        let mut w = World::new(sc.build(0xC0FFEE), 20.0);
        let root = w.tree.root;
        w.tree.nodes[root.get()].spec.count = count;
        w.tree.refine(root);
        let tier = w.tree.nodes[root.get()].tier;
        let node = w.tree.promote(root, 0, default_spec(tier.finer()));
        w.tree.nodes[node.get()].spec.count = count;
        w.tree.refine(node);
        let (_, _, components) = w.tree.components(node);
        let split = w.tree.resolve_extent(node);
        println!(
            "  a planet drawn with {count} bodies: {components} components, split {}",
            split.is_some()
        );
        assert!(
            split.is_none(),
            "a planet nobody has touched split into two at a count of {count}"
        );
    }
}

/// Two clumps that fall back together become one node again.
///
/// Without the inverse, `docs/BACKLOG.md` says, "two clumps that fall back
/// together stay two nodes forever" — and every event that ever separated
/// anything leaves a node behind.
#[test]
fn two_clumps_that_fall_back_together_become_one_node() {
    let (mut w, root, node) = a_cloud_in_two_minds(0xBEEF);
    let new = w.tree.resolve_extent(node).expect("split first");
    let before = w.conserved();
    let live_before = w.tree.live_count();
    let bodies_before = w.tree.nodes[node.get()].bodies.iter().filter(|b| b.mass > 0.0).count()
        + w.tree.nodes[new.get()].bodies.len();

    assert!(!w.tree.would_merge(node, new), "they have only just come apart");

    // Put them back together: the piece that left goes back over the centre of
    // mass of what stayed, which is where it was drawn from.
    w.tree.nodes[new.get()].motion.offset = w.tree.nodes[node.get()].motion.offset;
    assert!(w.tree.would_merge(node, new), "sitting on top of each other and still two");
    assert!(w.tree.merge(node, new), "the merge was refused");

    let n = &w.tree.nodes[node.get()];
    println!(
        "  merged: {} bodies, {:.3e} kg, live nodes {} -> {}",
        n.bodies.iter().filter(|b| b.mass > 0.0).count(),
        n.matter.mass,
        live_before,
        w.tree.live_count()
    );
    assert!(!w.tree.nodes[new.get()].alive, "the node that was folded in is still alive");
    assert_eq!(
        w.tree.nodes[root.get()].children.iter().filter(|c| **c == new).count(),
        0,
        "the parent still points at what it absorbed"
    );
    assert_eq!(
        w.tree.nodes[node.get()].bodies.iter().filter(|b| b.mass > 0.0).count(),
        bodies_before,
        "detail went missing in the merge"
    );
    let after = w.conserved();
    let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(y.abs()).max(1e-300);
    assert!(rel(before.energy, after.energy) < 1e-9, "energy changed across a merge");
    assert!(rel(before.baryon, after.baryon) < 1e-9, "baryon number changed across a merge");
}

/// The middle outcome: the radius follows contents that are the authority, and
/// does not follow a drawing of itself.
#[test]
fn the_radius_follows_contents_somebody_has_touched() {
    // Touched: pinned, so the bodies are what the node is.
    let (mut w, _, node) = a_cloud_in_two_minds(0xF00D);
    let claimed = w.tree.nodes[node.get()].matter.radius;
    // One clump rather than two, so the answer is the radius rather than a
    // split: pull everything outward instead of apart.
    for b in w.tree.nodes[node.get()].bodies.iter_mut() {
        b.pos = b.pos.scale(4.0);
    }
    w.tree.resolve_extent(node);
    let grown = w.tree.nodes[node.get()].matter.radius;
    println!("  touched: radius {claimed:.4e} -> {grown:.4e}");
    assert!(grown > claimed * 1.5, "the radius did not follow contents that had spread");

    // Untouched: the matter is the authority and the bodies are a drawing of
    // it, so the radius is an input to the detail rather than a measurement of
    // it. Measured on the biome world when this was the other way round: a
    // star's parcels dispersed, its radius followed them from 7e8 to 2.9e11 m,
    // its radiating area grew with it, and it cooled from 5800 K to 719 K.
    let sc = phys::scenario::ALL.iter().find(|s| s.name == "The Sun").unwrap();
    let mut w = World::new(sc.build(0x5A11), 20.0);
    let root = w.tree.root;
    w.tree.nodes[root.get()].spec.count = 64;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    let star = w.tree.promote(root, 0, default_spec(tier.finer()));
    w.tree.refine(star);
    let before = w.tree.nodes[star.get()].matter.radius;
    for b in w.tree.nodes[star.get()].bodies.iter_mut() {
        b.pos = b.pos.scale(4.0);
    }
    w.tree.resolve_extent(star);
    let after = w.tree.nodes[star.get()].matter.radius;
    println!("  untouched: radius {before:.4e} -> {after:.4e}");
    assert_eq!(
        before.to_bits(),
        after.to_bits(),
        "a regenerable node's radius followed a drawing of itself"
    );
}

/// End to end: the frame loop does it, on the cadence of the nodes it advanced.
#[test]
fn the_frame_loop_resolves_what_it_advanced() {
    let (mut w, root, node) = a_cloud_in_two_minds(0xBEEF);
    // Paced to the node, and somebody has to be looking. A world runs at one
    // second per second (`PLAY.md` D1) and a molecular cloud does nothing at
    // all in a twentieth of a second, so with neither of these the plan accepts
    // no work — which is correct, and is why the pass runs over what the plan
    // *did* advance rather than over every node in the world.
    w.pace_to(node);
    w.add_observer(phys::observe::Observer {
        anchor: root,
        offset: phys::math::v3(4.0 * w.tree.nodes[root.get()].matter.radius, 0.0, 0.0),
        look: phys::math::v3(-1.0, 0.0, 0.0),
        field: 1.2,
        angular_resolution: 1e-7,
        horizon: 1e6 * YEAR,
        priority: 100.0,
        ..Default::default()
    });
    let before = w.tree.live_count();
    for _ in 0..20 {
        w.step_frame(50_000.0);
    }
    println!(
        "  {} splits, {} merges, live nodes {before} -> {}",
        w.stats.splits,
        w.stats.merges,
        w.tree.live_count()
    );
    assert!(w.stats.splits > 0, "the frame loop never noticed");
    assert!(
        w.tree.nodes[node.get()].alive,
        "the node that split should still be the node it was"
    );
}
