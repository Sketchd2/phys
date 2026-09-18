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

        let before = world.tree.nodes[root.get()].matter;
        world.tree.coarsen(root);
        let after = world.tree.nodes[root.get()].matter;

        let scale = before.mass.abs().max(1e-30);
        let mass_error = (after.mass - before.mass).abs() / scale;
        let energy = before
            .internal_energy
            .abs()
            .max(before.gravitational_binding.abs())
            .max(before.cohesive_binding.abs())
            .max(1e-30);
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
    let path = world.drill_to(root, Tier::Nuclear.max_radius(), &default_spec);
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
    let arrived = world.tree.nodes[path.last().unwrap().get()].matter.radius;
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
        let before = world.tree.nodes[root.get()].matter.mass;

        let child = world.tree.promote(root, 0, default_spec(s.tier.finer()));
        assert!(!child.is_none(), "{}: could not descend", s.name);
        let inner = world.tree.refine(child).len();
        world.tree.coarsen(child);
        world.tree.coarsen(root);
        let after = world.tree.nodes[root.get()].matter.mass;
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
            world.tree.nodes[cur.get()].matter.radius,
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

/// Every scenario's contents land inside the node they belong to.
///
/// `docs/PLAY.md` §7's second Phase 2 item, and the `docs/BACKLOG.md` entry
/// "The sampler inflates anything bound by chemistry by 4.3x10^5".
///
/// `sample` has a relaxation loop whose premise is that a configuration too
/// tightly bound to hold the energy it claims must be bigger, so it scales the
/// geometry up by 1.5 until the budget turns positive — up to thirty-two
/// times. Scaling moves the gravitational potential and nothing else, so
/// against a *chemical* deficit it cannot work: it ran all thirty-two
/// iterations, failed, fell through to a fallback and left `1.5^32 = 4.3x10^5`
/// of inflation in place that nothing undid.
///
/// Nothing geometric works on a node in that state — adjacency, contact, the
/// exchange pass and hydro's own neighbour finding all see contents scattered
/// hundreds of thousands of radii from a node they are supposed to be inside —
/// and the conserved-set tests never caught it, because the books close either
/// way: the energy budget absorbs whatever `phi` comes out to, which is the
/// sampler's design and is correct.
///
/// So this is a *geometric* assertion, deliberately, and the first one any
/// test has made about a Continuum node. Measured, on this exact probe, before
/// and after the split of `binding_energy` into a gravitational term and a
/// cohesive one:
///
/// ```text
///   scenario          relax before / after      max|pos|/R before / after
///   Spiral galaxy          0        0                3.354      3.354
///   Molecular cloud        0        0                3.574      3.574
///   The Sun                0        0                3.515      3.515
///   Rocky planet           0        0                3.515      3.515
///   Granite block         33        0            4.502e5        1.043
///   Water vapour          33        0            4.387e5        1.017
///   Carbon atom           33        0            4.738e5        1.098
///   Iron nucleus           0        0                1.080      1.080
/// ```
///
/// The four gravitationally bound scenarios are untouched to every digit,
/// which is the control: the loop is still there and still does its job for
/// the case it was written for. The nucleus is the other control — its binding
/// is nuclear and equally unmoved by expansion, and it never tripped the loop
/// only because its Fermi energy exceeds it.
///
/// The bound is 10 rather than 4, because a Plummer sphere genuinely has a
/// tail: the four gravitational rows sit at 3.4-3.6 and that is what a relaxed
/// self-gravitating configuration looks like. What is being caught is five
/// orders of magnitude, not a factor of two.
#[test]
fn every_scenario_samples_its_contents_inside_itself() {
    use phys::sampler::{sample, SampleSpec};

    for s in scenario::ALL {
        let tree = s.build(0x5EED);
        let node = &tree.nodes[tree.root.get()];
        let matter = node.matter;
        let spec = SampleSpec { count: 512, ..node.spec };
        let (bodies, report) = sample(&matter, spec, tree.world_seed, node.key.0, 0);
        assert!(!bodies.is_empty(), "{}: sampled nothing", s.name);

        let furthest =
            bodies.iter().map(|b| b.pos.norm()).fold(0.0f64, f64::max) / matter.radius;
        println!(
            "  {:<16} relax {:>2}   max|pos|/R {:>10.3e}   E_grav {:>10.3e}   E_cohesive {:>10.3e}",
            s.name, report.relaxations, furthest, matter.gravitational_binding,
            matter.cohesive_binding
        );
        assert!(
            furthest < 10.0,
            "{}: its furthest body is {furthest:.3e} radii out, after \
             {} relaxations — the sampler expanded the configuration trying to \
             release a binding that expansion does not release",
            s.name,
            report.relaxations
        );
    }
}
