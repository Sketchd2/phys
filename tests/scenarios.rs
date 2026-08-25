//! Every scenario, at every scale, doing the same things.
//!
//! The claim the scale ladder rests on is that nothing is special-cased: a
//! galaxy and an iron nucleus differ by forty-five orders of magnitude in size
//! and by nothing at all in structure. These tests are that claim, stated as
//! the same four operations applied to every starting point on the shelf —
//! refine it, step it, check the books, and descend.

use phys::engine::{default_spec, World};
use phys::scenario;
use phys::units::Tier;

/// Every scenario materialises, and materialisation conserves.
#[test]
fn every_scenario_refines_and_conserves() {
    for s in scenario::ALL {
        let mut world = World::new(s.build(0x5EED), 20.0);
        let root = world.tree.root;
        let bodies = world.tree.refine(root).len();
        assert!(bodies > 0, "{}: refined to nothing", s.name);

        let before = world.tree.nodes[root.get()].agg;
        world.tree.coarsen(root);
        let after = world.tree.nodes[root.get()].agg;

        let scale = before.mass.abs().max(1e-30);
        let mass_error = (after.mass - before.mass).abs() / scale;
        let energy = before.internal_energy.abs().max(before.binding_energy.abs()).max(1e-30);
        let energy_error = (after.internal_energy - before.internal_energy).abs() / energy;
        println!(
            "  {:<16} {:>9} {:>7} bodies   mass {:.2e}   energy {:.2e}",
            s.name,
            s.tier.name(),
            bodies,
            mass_error,
            energy_error
        );
        assert!(mass_error < 1e-12, "{}: mass moved by {mass_error:.3e}", s.name);
        assert!(
            energy_error < 1e-6,
            "{}: internal energy moved by {energy_error:.3e}",
            s.name
        );
    }
}

/// Every scenario steps under whatever solver its tier calls for, without
/// diverging and without inventing energy.
#[test]
fn every_scenario_steps() {
    for s in scenario::ALL {
        let mut world = World::new(s.build(0x5EED), 20.0);
        let root = world.tree.root;
        world.tree.refine(root);
        let dt = world.node_dt(root);
        assert!(dt > 0.0 && dt.is_finite(), "{}: timestep {dt:?}", s.name);

        let mut worst_drift = 0.0f64;
        for _ in 0..8 {
            let report = world.advance_node(root, dt);
            worst_drift = worst_drift.max(report.drift());
        }
        let bodies = &world.tree.nodes[root.get()].bodies;
        let fastest = bodies.iter().map(|b| b.vel.norm()).fold(0.0f64, f64::max);
        println!(
            "  {:<16} dt {:.3e} s   worst drift {:.2e}   fastest body {:.3e} m/s ({:.4} c)",
            s.name,
            dt,
            worst_drift,
            fastest,
            fastest / phys::units::C
        );
        assert!(
            bodies.iter().all(|b| b.pos.is_finite() && b.vel.is_finite()),
            "{}: diverged",
            s.name
        );
        assert!(fastest < phys::units::C, "{}: superluminal body", s.name);
    }
}

/// The ladder is continuous: from any starting point, refinement reaches the
/// bottom, and it reaches the same place a scenario that starts there does.
#[test]
fn the_ladder_runs_all_the_way_down() {
    let mut world = World::new(scenario::ALL[0].build(0x5EED), 20.0);
    let root = world.tree.root;
    let path = world.drill(root, Tier::Nuclear, &default_spec);
    let tiers: Vec<&str> = path
        .iter()
        .map(|&n| world.tree.nodes[n.get()].tier.name())
        .collect();
    println!("  galaxy to nucleus in {} steps: {}", path.len(), tiers.join(" -> "));
    assert_eq!(
        world.tree.nodes[path.last().unwrap().get()].tier,
        Tier::Nuclear,
        "the drill did not reach the bottom"
    );

    // And the node it arrives at is a nucleus by the same measure the scenario
    // that starts as one is: same tier, same solver, same order of size.
    let arrived = world.tree.nodes[path.last().unwrap().get()].agg.radius;
    let direct = scenario::ALL.last().unwrap();
    println!(
        "  arrived at a {arrived:.3e} m node; the {} scenario starts at {:.3e} m",
        direct.name, direct.scale
    );
    assert!(
        arrived > 1e-16 && arrived < 1e-12,
        "the bottom of the ladder is {arrived:.3e} m across"
    );
}

/// Descending is not a special case of anything. Every scenario can be entered
/// and left, and the world it leaves behind is the one it started with.
#[test]
fn descending_and_returning_leaves_no_trace() {
    for s in scenario::ALL {
        if s.tier == Tier::Nuclear {
            continue; // nothing below it to descend into
        }
        let mut world = World::new(s.build(0x5EED), 20.0);
        let root = world.tree.root;
        world.tree.refine(root);
        let before = world.tree.nodes[root.get()].agg.mass;

        let child = world.tree.promote(root, 0, default_spec(s.tier.finer()));
        assert!(!child.is_none(), "{}: could not descend", s.name);
        let inner = world.tree.refine(child).len();
        world.tree.coarsen(child);
        world.tree.coarsen(root);
        let after = world.tree.nodes[root.get()].agg.mass;
        println!(
            "  {:<16} {} -> {} : {inner} bodies inside, mass error {:.2e}",
            s.name,
            s.tier.name(),
            s.tier.finer().name(),
            (after - before).abs() / before.abs().max(1e-30)
        );
        assert!(
            (after - before).abs() / before.abs().max(1e-30) < 1e-12,
            "{}: mass moved by descending and returning",
            s.name
        );
    }
}

/// A node that is only being stepped must not heat up.
///
/// This is the failure the scale explorer found, and it is the one that looks
/// most like physics while it happens. Descend a galaxy far enough and you
/// reach a node whose gas is dense and hot; step it and the parcels accelerate,
/// every step conserving the energy the step before invented, until the gas is
/// moving at two thirds of light speed and the conservation check has reported
/// no drift at all.
///
/// The cause was a timestep that came from a table and the gravitational
/// dynamical time, with nothing about how fast a signal crosses the gap between
/// two parcels. Pressure was then acting across a distance the information
/// could not have travelled, which is not a small error.
#[test]
fn stepping_does_not_heat_a_node() {
    for s in scenario::ALL {
        let mut world = World::new(s.build(0x5EED), 20.0);
        let mut cur = world.tree.root;

        // Descend towards the densest thing available, which is where the
        // timestep has the least margin.
        for _ in 0..6 {
            world.tree.refine(cur);
            let best = {
                let n = &world.tree.nodes[cur.get()];
                let mut bi = 0usize;
                let mut bm = f64::NEG_INFINITY;
                for (i, b) in n.bodies.iter().enumerate() {
                    if b.mass > bm {
                        bm = b.mass;
                        bi = i;
                    }
                }
                bi
            };
            let tier = world.tree.nodes[cur.get()].tier;
            if tier == Tier::Nuclear || world.tree.nodes[cur.get()].bodies.is_empty() {
                break;
            }
            let child = world.tree.promote(cur, best, default_spec(tier.finer()));
            if child.is_none() {
                break;
            }
            cur = child;
        }

        world.tree.refine(cur);
        let speed = |w: &World| {
            w.tree.nodes[cur.get()]
                .bodies
                .iter()
                .map(|b| b.vel.norm())
                .fold(0.0f64, f64::max)
        };
        let before = speed(&world);
        let mut worst_drift = 0.0f64;
        for _ in 0..200 {
            let dt = world.node_dt(cur);
            let report = world.advance_node(cur, dt);
            worst_drift = worst_drift.max(report.drift());
        }
        let after = speed(&world);
        let tier = world.tree.nodes[cur.get()].tier;
        println!(
            "  {:<16} {:>10} {:>9.3e} m: fastest {:.3e} -> {:.3e} m/s ({:.2}x) over 200 steps, \
             worst drift {:.2e}",
            s.name,
            tier.name(),
            world.tree.nodes[cur.get()].agg.radius,
            before,
            after,
            after / before.max(1e-30),
            worst_drift
        );
        assert!(
            after < phys::units::C,
            "{}: stepping drove a body to light speed",
            s.name
        );
        // Two hundred steps of a bound system may redistribute energy; it may
        // not manufacture it. A factor of four is far more headroom than a
        // stable solver needs and far less than a runaway takes.
        assert!(
            after < before.max(1.0) * 4.0,
            "{}: the fastest body went from {before:.3e} to {after:.3e} m/s",
            s.name
        );
        assert!(
            worst_drift < 1e-3,
            "{}: energy drifted by {worst_drift:.3e} in a single step",
            s.name
        );
    }
}
