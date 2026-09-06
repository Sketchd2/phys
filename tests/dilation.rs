//! How fast a node's own clock runs, and what happens when somebody overrides
//! it.
//!
//! Two claims, one mechanism.
//!
//! The physical one: relativity must actually *drive* evolution, not merely be
//! reported. Before this, `proper_time` was accumulated and never read, and
//! `gravitational_dilation` was written and never called, so a node moving at
//! 0.99c aged at exactly the same rate as one standing still.
//!
//! The unphysical one: an administrator can put a region in a time bubble and
//! run its interior faster than the universe around it, for testing and
//! balancing. That is the same multiplier in the same product, reported
//! separately so nothing can confuse the two.

use phys::coords::gamma;
use phys::dilation::{accept_bubble, physical_rate, TimeRate, MAX_BUBBLE, MIN_BUBBLE};
use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::observe::{Interaction, Property};
use phys::units::*;

fn a_world() -> World {
    let mut w = World::new(galaxy(0x71E5, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 512;
    w.time_rate = 0.05;
    w
}

/// A node deep enough to be worth stepping, and its parent.
fn deep(w: &mut World) -> NodeIdx {
    let root = w.tree.root;
    let d = *w.drill(root, Tier::Continuum, &default_spec).last().unwrap();
    w.tree.refine(d);
    w.pace_to(d);
    d
}

// ---------------------------------------------------------------------------
// the arithmetic
// ---------------------------------------------------------------------------

#[test]
fn the_rate_is_a_product_of_three_named_factors() {
    let r = TimeRate { kinematic: 0.5, gravitational: 0.8, bubble: 10.0 };
    assert_eq!(r.total(), 0.5 * 0.8 * 10.0);
    assert_eq!(r.physical(), 0.5 * 0.8);
    assert!(!r.is_physical());
    assert!(TimeRate::default().is_physical());
    assert_eq!(TimeRate::default().total(), 1.0);

    // Composing up the tree multiplies each factor with its own kind, so a
    // bubble can never be mistaken for gravity in a log.
    let c = r.compose(TimeRate { kinematic: 0.5, gravitational: 1.0, bubble: 2.0 });
    assert_eq!(c.kinematic, 0.25);
    assert_eq!(c.gravitational, 0.8);
    assert_eq!(c.bubble, 20.0);
}

/// A rate that has gone to zero or NaN would freeze or poison a node rather
/// than slow it, so it fails to "normal" instead.
#[test]
fn a_broken_rate_fails_to_one() {
    for bad in [0.0, f64::NAN, f64::INFINITY, -1.0] {
        let r = TimeRate { kinematic: bad, gravitational: 1.0, bubble: 1.0 };
        assert_eq!(r.total(), 1.0, "a rate of {bad} should fall back to 1");
    }
}

#[test]
fn speed_dilates_and_potential_dilates() {
    let still = physical_rate(v3(0.0, 0.0, 0.0), 0.0, 1.0);
    assert_eq!(still.kinematic, 1.0);
    assert_eq!(still.gravitational, 1.0);

    // 0.6c: gamma is exactly 1.25, so the clock runs at 0.8.
    let fast = physical_rate(v3(0.6 * C, 0.0, 0.0), 0.0, 1.0);
    println!("  0.6c gives 1/gamma = {:.6}", fast.kinematic);
    assert!((fast.kinematic - 0.8).abs() < 1e-9);

    // A potential well. `external_potential` is an energy, so it is divided by
    // the mass here — a fact easy to get wrong and worth pinning.
    let mass = 2.0;
    let phi = -0.18 * C * C;
    let low = physical_rate(v3(0.0, 0.0, 0.0), phi * mass, mass);
    println!("  phi = -0.18 c^2 gives sqrt(1+2phi/c^2) = {:.6}", low.gravitational);
    assert!((low.gravitational - (1.0f64 - 0.36).sqrt()).abs() < 1e-12);
    assert!(low.gravitational < 1.0);

    // Nothing there, no clock.
    let empty = physical_rate(v3(0.0, 0.0, 0.0), -1e30, 0.0);
    assert_eq!(empty.gravitational, 1.0);
}

#[test]
fn a_bubble_request_is_clamped_not_refused() {
    assert_eq!(accept_bubble(1.0), Some(1.0));
    assert_eq!(accept_bubble(1e12), Some(MAX_BUBBLE));
    assert_eq!(accept_bubble(1e-12), Some(MIN_BUBBLE));
    // Not an ambitious request — a caller bug.
    assert_eq!(accept_bubble(0.0), None);
    assert_eq!(accept_bubble(-2.0), None);
    assert_eq!(accept_bubble(f64::NAN), None);
}

// ---------------------------------------------------------------------------
// relativity actually drives evolution now
// ---------------------------------------------------------------------------

/// The bug this exists to prevent: proper time computed and thrown away.
#[test]
fn a_fast_node_evolves_on_its_own_clock() {
    let mut w = a_world();
    let d = deep(&mut w);

    // Not 1 — a node seven frames deep in a galaxy is already dilated, by its
    // ancestors' orbital speeds and by the halo it sits in. That is the point:
    // these numbers were always there and nothing used to read them.
    let resting = w.time_rate_of(d);
    println!(
        "  at rest in a galaxy: kinematic {:.9}, gravitational {:.9}",
        resting.kinematic, resting.gravitational
    );
    assert!(resting.kinematic < 1.0 && resting.kinematic > 0.99);
    assert!(resting.gravitational < 1.0, "sitting in a halo dilates a clock");
    assert!(resting.is_physical());

    // Now boost this node to 0.8c in its parent's frame. Gamma is 5/3, so its
    // own contribution should be exactly 0.6 — asserted as a ratio against the
    // ancestry it already had, which is the only part this test controls.
    w.tree.nodes[d.get()].frame.velocity = v3(0.8 * C, 0.0, 0.0);
    let moving = w.time_rate_of(d);
    let own = moving.kinematic / resting.kinematic;
    println!(
        "  boosted to 0.8c (gamma {:.4}): its own factor is {own:.9}, total rate {:.9}",
        gamma(v3(0.8 * C, 0.0, 0.0)),
        moving.total()
    );
    assert!((own - 0.6).abs() < 1e-6, "0.8c must contribute exactly 1/gamma = 0.6");
    assert!(moving.is_physical(), "relativity is not an administrative act");
    assert!(moving.total() < resting.total(), "going faster must slow the interior");
}

/// Dilation composes up the chain of frames, because each node's velocity is
/// measured in its parent's.
#[test]
fn dilation_composes_up_the_tree() {
    let mut w = a_world();
    let d = deep(&mut w);
    let parent = w.tree.nodes[d.get()].parent;
    assert!(!parent.is_none(), "need a parent for this to mean anything");

    w.tree.nodes[d.get()].frame.velocity = v3(0.6 * C, 0.0, 0.0);
    w.tree.nodes[parent.get()].frame.velocity = v3(0.0, 0.6 * C, 0.0);

    let child = w.time_rate_of(d);
    let above = w.time_rate_of(parent);
    println!(
        "  parent runs at {:.6}, child at {:.6} — child/parent = {:.6}",
        above.kinematic,
        child.kinematic,
        child.kinematic / above.kinematic
    );
    // The child carries its own 0.8 and everything its parent carries.
    assert!((child.kinematic / above.kinematic - 0.8).abs() < 1e-9);
    assert!(child.kinematic < above.kinematic);
}

/// The trajectory is not dilated. A fast object must still cross the world at
/// its own speed, or it would appear to crawl.
#[test]
fn a_dilated_node_still_travels_at_its_velocity() {
    let mut w = a_world();
    let d = deep(&mut w);
    let v = v3(0.0, 0.0, 0.9 * C);
    w.tree.nodes[d.get()].frame.velocity = v;
    let before = w.tree.nodes[d.get()].frame.offset;
    let t0 = w.tree.nodes[d.get()].time;

    w.advance_node(d, 1.0);

    let n = &w.tree.nodes[d.get()];
    let moved = (n.frame.offset - before).norm();
    let coordinate = n.time - t0;
    println!(
        "  {coordinate:.6} s of world time moved it {moved:.6e} m; v*dt would be {:.6e}",
        v.norm() * coordinate
    );
    assert!(
        (moved - v.norm() * coordinate).abs() / (v.norm() * coordinate) < 1e-9,
        "position must advance on coordinate time, not local time"
    );
    assert!((coordinate - 1.0).abs() < 1e-9, "one second asked, one second of world time");
}

// ---------------------------------------------------------------------------
// bubbles
// ---------------------------------------------------------------------------

#[test]
fn a_bubble_multiplies_the_interior_and_nothing_else() {
    let mut w = a_world();
    let d = deep(&mut w);

    assert_eq!(w.dilate(d, 100.0), Some(100.0));
    let r = w.time_rate_of(d);
    println!("  kinematic {:.6} gravitational {:.6} bubble {}", r.kinematic, r.gravitational, r.bubble);
    assert_eq!(r.bubble, 100.0);
    assert!(!r.is_physical(), "a bubble must never look like physics");
    assert!((r.total() / r.physical() - 100.0).abs() < 1e-9);

    // The trajectory is untouched.
    w.tree.nodes[d.get()].frame.velocity = v3(1000.0, 0.0, 0.0);
    let before = w.tree.nodes[d.get()].frame.offset;
    let t0 = w.tree.nodes[d.get()].time;
    w.advance_node(d, 1.0);
    let n = &w.tree.nodes[d.get()];
    let moved = (n.frame.offset - before).norm();
    let coordinate = n.time - t0;
    println!("  bubbled 100x, one second: moved {moved:.3} m over {coordinate:.6} s of world time");
    assert!(
        (moved - 1000.0 * coordinate).abs() < 1e-6,
        "a bubble must not move the object faster"
    );
}

/// A bubble is inherited by everything below it, because it is a property of
/// the subtree and not of one node.
#[test]
fn a_bubble_covers_the_subtree() {
    let mut w = a_world();
    let d = deep(&mut w);
    let parent = w.tree.nodes[d.get()].parent;

    w.dilate(parent, 50.0);
    assert_eq!(w.time_rate_of(parent).bubble, 50.0);
    assert_eq!(w.time_rate_of(d).bubble, 50.0, "a child inherits its ancestors' bubbles");

    // And they nest, multiplying, like every other factor in the product.
    w.dilate(d, 4.0);
    assert_eq!(w.time_rate_of(d).bubble, 200.0);
    assert_eq!(w.time_rate_of(parent).bubble, 50.0, "a child's bubble must not leak upward");
}

#[test]
fn a_bubble_is_lifted_by_setting_it_back_to_one() {
    let mut w = a_world();
    let d = deep(&mut w);
    w.dilate(d, 1000.0);
    assert!(!w.time_rate_of(d).is_physical());
    w.dilate(d, 1.0);
    assert!(w.time_rate_of(d).is_physical(), "rate 1 must remove the bubble");
    assert!(w.bubbles().is_empty());
}

#[test]
fn a_bubble_makes_a_node_proportionately_late() {
    let mut w = a_world();
    let d = deep(&mut w);
    let horizon = w.tree.nodes[d.get()].last_solved + 1.0;

    let plain = w.lateness(d, horizon);
    w.dilate(d, 100.0);
    let bubbled = w.lateness(d, horizon);
    println!("  lateness {plain:.6e} unbubbled, {bubbled:.6e} at 100x");
    assert!(
        (bubbled / plain - 100.0).abs() < 1e-6,
        "a node living 100x faster is 100x as late, so the scheduler prioritises it"
    );
}

/// Thermalising replaces a node's contents with a fresh draw from its
/// equilibrium. That is right for matter nobody is following and exactly wrong
/// for a region somebody has deliberately sped up to watch.
#[test]
fn a_bubbled_node_is_never_thermalised() {
    let mut w = a_world();
    let d = deep(&mut w);
    assert!(w.forgettable(d), "this node should be forgettable to begin with");
    w.dilate(d, 10.0);
    assert!(!w.forgettable(d), "a bubbled node must not be crossed by its ensemble");
    w.dilate(d, 1.0);
    assert!(w.forgettable(d));
}

#[test]
fn a_bubble_is_audited() {
    let mut w = a_world();
    let d = deep(&mut w);
    let before = w.audit.len();
    w.interact(Interaction::Dilate { target: d, rate: 250.0 });

    assert_eq!(w.audit.len(), before + 1, "an unphysical act must leave a record");
    let e = w.audit.last().unwrap();
    println!("  audit: {:?} at t={:.3e}, delta_energy {:.3e}", e.property, e.time, e.delta_energy);
    assert_eq!(e.property, Property::TimeRate);
    assert_eq!(e.key, w.tree.nodes[d.get()].key);
    // Setting a rate injects no energy. The divergence is a flow, not a step.
    assert_eq!(e.delta_energy, 0.0);
    assert_eq!(w.bubbles(), vec![(d, 250.0)]);
}

/// A bubble must not pin, bump the epoch, or disturb — see `World::dilate`.
#[test]
fn a_bubble_does_not_touch_the_detail() {
    let mut w = a_world();
    let d = deep(&mut w);
    let (epoch, pinned, disturbed) = {
        let n = &w.tree.nodes[d.get()];
        (n.epoch, n.pinned, n.last_disturbed)
    };
    w.dilate(d, 100.0);
    let n = &w.tree.nodes[d.get()];
    assert_eq!(n.epoch, epoch, "a bubble does not invalidate detail");
    assert_eq!(n.pinned, pinned, "a bubble does not make detail non-regenerable");
    assert_eq!(n.last_disturbed, disturbed, "a bubble is not an event in the node's history");
}

#[test]
fn a_bad_dilate_request_is_reported() {
    let mut w = a_world();
    let d = deep(&mut w);
    assert_eq!(w.dilate(d, 0.0), None);
    assert_eq!(w.dilate(d, f64::NAN), None);
    assert_eq!(w.dilate(NodeIdx::NONE, 2.0), None);
    assert_eq!(w.dilate(d, 1e30), Some(MAX_BUBBLE), "an ambitious request is clamped");
    assert!(w.audit.iter().all(|e| e.property != Property::TimeRate || e.time.is_finite()));
}

/// An untouched world has been given nothing it did not earn. This is the
/// assertion that makes `bubble_seconds` worth having.
#[test]
fn an_untouched_world_owes_nothing() {
    let mut w = a_world();
    let _ = deep(&mut w);
    for _ in 0..8 {
        w.step_frame(20_000.0);
    }
    println!(
        "  {} frames, bubble_seconds {:.3e}, bubbled nodes {}",
        w.stats.frames, w.stats.bubble_seconds, w.stats.bubbled
    );
    assert_eq!(w.stats.bubbled, 0);
    assert_eq!(w.stats.bubble_seconds, 0.0);
}

#[test]
fn a_bubble_is_metered() {
    let mut w = a_world();
    let d = deep(&mut w);
    w.dilate(d, 100.0);
    for _ in 0..8 {
        w.step_frame(20_000.0);
    }
    println!(
        "  after {} frames at 100x: {:.3e} interior seconds granted beyond the clock, {} bubbled",
        w.stats.frames, w.stats.bubble_seconds, w.stats.bubbled
    );
    assert_eq!(w.stats.bubbled, 1);
    assert!(
        w.stats.bubble_seconds > 0.0,
        "a bubble that ran must show up in the meter"
    );
}

#[test]
fn a_bubble_survives_a_save() {
    use phys::persist::{MemoryStore, WorldStore};
    let mut w = a_world();
    let d = deep(&mut w);
    w.dilate(d, 750.0);
    let key = w.tree.nodes[d.get()].key;

    let mut store = MemoryStore::default();
    store.save(w.view()).expect("save");
    let back = World::from_snapshot(store.load().expect("load"), 20.0);

    let found = back
        .tree
        .nodes
        .iter()
        .find(|n| n.alive && n.key == key)
        .expect("the node came back");
    println!("  bubble {} survived the round trip", found.bubble);
    assert_eq!(found.bubble, 750.0);
}
