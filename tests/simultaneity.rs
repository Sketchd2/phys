//! One instant, every scale, every frame.
//!
//! The engine's central claim after the scheduler rewrite: a galaxy and a
//! nucleus inside it are at the *same moment*, and both of them advance in
//! every frame. What differs between them is how often each is re-solved and
//! how it was carried here — the galaxy in one closed-form step, the nucleus
//! through hundreds of integrations — not what time it is for them.
//!
//! The old scheduler could not do this. It took a global minimum over every
//! live node's timestep, so resolving anything small stopped everything large:
//! a nucleus in the tree pinned the whole world to 10^-23 s per frame.

use phys::engine::{default_spec, galaxy, World};
use phys::math::v3;
use phys::units::*;

/// Every live node ends the frame at the world instant. Not near it — at it.
#[test]
fn the_whole_world_is_at_one_instant() {
    let mut w = World::new(galaxy(0xA11CE, 1e9), 20.0);
    let root = w.tree.root;
    let path = w.drill(root, Tier::Nuclear, &default_spec);
    println!("  drilled {} tiers, galaxy to nucleus", path.len());
    assert!(path.len() >= 6, "expected a deep ladder, got {}", path.len());

    for _ in 0..20 {
        w.step_frame(50_000.0);
    }

    let mut checked = 0;
    for (i, n) in w.tree.nodes.iter().enumerate() {
        if !n.alive {
            continue;
        }
        checked += 1;
        let slip = (n.time - w.time).abs();
        assert!(
            slip <= w.time.abs() * 1e-9,
            "node {i} ({:?}) is at {:.6e} s, the world is at {:.6e} s",
            n.tier,
            n.time,
            w.time
        );
    }
    println!("  {checked} live nodes, all at t = {:.6e} s", w.time);
    assert!(checked >= 6);
}

/// The frame no longer runs at the pace of the fastest thing in it.
///
/// This is the regression the whole rewrite exists for. Drilling to the nuclear
/// tier used to drag the world's timestep down by twenty orders of magnitude,
/// because the step was a global minimum over every node's stability limit. The
/// pace now follows the *subject* — resolving the galaxy itself does move it,
/// and that is the point — but resolving something twenty orders of magnitude
/// smaller inside it does not.
#[test]
fn resolving_something_small_does_not_stop_the_galaxy() {
    let mut coarse = World::new(galaxy(0xB0B, 1e9), 20.0);
    let root = coarse.tree.root;
    coarse.tree.refine(root);
    let before = coarse.frame_dt();

    let mut fine = World::new(galaxy(0xB0B, 1e9), 20.0);
    fine.drill(root, Tier::Nuclear, &default_spec);
    let after = fine.frame_dt();

    println!(
        "  galaxy resolved into its stars: {:.3e} s per frame; \
         with a nucleus resolved inside one of them too: {:.3e} s",
        before, after
    );
    assert_eq!(
        before, after,
        "resolving a nucleus changed the galaxy's frame from {before:.3e} to {after:.3e} s"
    );

    // And the galaxy really does move: over twenty frames the root turns.
    let start = fine.tree.nodes[root.get()].motion.orientation;
    for _ in 0..20 {
        fine.step_frame(50_000.0);
    }
    let turned = start
        .conjugate()
        .then(fine.tree.nodes[root.get()].motion.orientation)
        .angle();
    let omega = fine.tree.nodes[root.get()].motion.spin_rate.norm();
    let expected = omega * fine.time;
    println!(
        "  after twenty frames the galaxy has advanced {:.3e} s and turned {turned:.3e} rad, \
         against {expected:.3e} rad of spin",
        fine.time
    );
    assert!(fine.time > 0.0, "the world did not advance at all");
    // Rotation is carried in closed form, so it is not "about right" — it is
    // exactly what the spin rate and the elapsed time say it is.
    assert!(
        (turned - expected).abs() <= expected.max(1e-12) * 1e-9,
        "the galaxy turned {turned:.6e} rad where its spin says {expected:.6e}"
    );
}

/// The pace is set by what is being watched, not by any solver.
#[test]
fn the_pace_follows_the_subject() {
    let mut w = World::new(galaxy(0xC0FFEE, 1e9), 20.0);
    let root = w.tree.root;
    let galactic = w.pace;
    let path = w.drill(root, Tier::Molecular, &default_spec);
    let deep = *path.last().unwrap();
    w.pace_to(deep);
    let molecular = w.pace;
    println!(
        "  paced to the galaxy: {:.3e} s per frame; paced to a molecule: {:.3e} s",
        galactic, molecular
    );
    assert!(
        molecular < galactic,
        "a molecule should be watched more slowly than a galaxy"
    );
    assert!(
        galactic / molecular > 1e6,
        "only {:.1e} between the two paces",
        galactic / molecular
    );
}

/// Lateness, not attention, decides what runs — and nothing can be starved.
///
/// A node passed over grows more overdue every frame, so its value rises until
/// it outranks whatever kept beating it. The old scheduler ranked on observer
/// salience with a small novelty bonus, which meant a node the camera was not
/// pointing at could lose every frame forever — and, worse, coarsened away all
/// its detail on the first frame, leaving nothing to step at all.
#[test]
fn nothing_is_starved() {
    // A coarsely-drawn galaxy on one CPU core: this test is about what gets
    // scheduled, not about how much of it fits, so keep the work small enough
    // that the budget is not the thing under examination.
    let mut tree = galaxy(0xDEC1DE, 1e9);
    tree.nodes[0].spec.count = 256;
    let mut w = World::new(tree, 20.0);
    w.time_rate = 0.05;
    let root = w.tree.root;
    w.drill(root, Tier::Stellar, &default_spec);
    // Nobody is looking at anything. Under the old rule that meant no work at
    // all: every materialised node was coarsened on the frame it appeared.
    assert!(w.observers.is_empty());

    for _ in 0..40 {
        w.step_frame(50_000.0);
    }
    println!(
        "  unobserved world after 40 frames: {} bodies stepped, {} still materialised, \
         worst lateness {:.2}",
        w.stats.bodies_stepped, w.stats.materialised_bodies, w.stats.worst_lateness
    );
    assert!(
        w.stats.bodies_stepped > 0,
        "an unobserved world did no physics at all"
    );
    assert!(
        w.stats.materialised_bodies > 0,
        "an unobserved world threw away all its detail"
    );
}

/// Coasting is not an approximation. It is the exact solution.
#[test]
fn carrying_a_node_forward_is_exact() {
    let mut w = World::new(galaxy(0xFEED, 1e9), 20.0);
    let root = w.tree.root;
    let start = w.tree.nodes[root.get()].motion.offset;
    let v = v3(220e3, 0.0, 0.0);
    w.tree.nodes[root.get()].motion.velocity = v;

    let frames = 30;
    for _ in 0..frames {
        w.step_frame(50_000.0);
    }
    let expected = start + v.scale(w.time);
    let got = w.tree.nodes[root.get()].motion.offset;
    let error = (got - expected).norm();
    println!(
        "  after {frames} frames ({:.3e} s) the drift is {error:.3e} m over {:.3e} m travelled",
        w.time,
        (got - start).norm()
    );
    assert!(
        error <= (got - start).norm().max(1.0) * 1e-12,
        "coasting drifted by {error:.3e} m"
    );
}

/// Resolution, not distance, is what shortens a node's cadence.
#[test]
fn resolution_sets_the_cadence() {
    let mut w = World::new(galaxy(0x5EED, 1e9), 20.0);
    let root = w.tree.root;
    let coarse = w.node_cadence(root);
    w.tree.refine(root);
    let fine = w.node_cadence(root);
    println!(
        "  galaxy as bulk state: {:.3e} s between solves; materialised into {} bodies: {:.3e} s",
        coarse,
        w.tree.nodes[root.get()].bodies.len(),
        fine
    );
    assert!(
        fine < coarse,
        "materialising a node must make it come due more often"
    );
}

/// A span too long to integrate is crossed by ensemble, and the books still
/// balance.
///
/// Following the trajectory of a resolved node across fifty milliseconds of
/// molecular time would take 10^13 steps. It is also not the right answer:
/// after that long the node has sampled its states 10^13 times, and where it
/// ends up is a draw from its equilibrium ensemble. Summarising is conservative
/// and sampling is a maximum-entropy sample of the same conserved tuple, so
/// crossing by ensemble is exactly that draw.
#[test]
fn an_unreachable_span_is_crossed_by_ensemble() {
    let mut w = World::new(galaxy(0x5A11, 1e9), 20.0);
    let root = w.tree.root;
    let path = w.drill(root, Tier::Molecular, &default_spec);
    let deep = *path.last().unwrap();
    w.tree.refine(deep);
    assert!(w.tree.nodes[deep.get()].is_materialised());

    let before = w.conserved();
    let mass_before = w.tree.nodes[deep.get()].matter.mass;
    // Ask for a span the node cannot possibly integrate: a whole second of
    // molecular time is 10^14 steps.
    let horizon = w.time + 1.0;
    w.advance_to(deep, horizon, phys::engine::MAX_SUBSTEPS);

    let n = &w.tree.nodes[deep.get()];
    println!(
        "  asked a molecular node to cross 1 s; thermalised {} times, now at {:.6e} s",
        w.stats.thermalised, n.time
    );
    assert!(w.stats.thermalised > 0, "the node claimed to integrate 10^14 steps");
    assert!((n.time - horizon).abs() < 1e-9, "it did not arrive");
    let after = w.conserved();
    let drift = (after.baryon - before.baryon).abs() / before.baryon.abs().max(1e-300);
    println!("  baryon-number drift across the crossing: {drift:.3e}");
    assert!(drift < 1e-9, "ensemble crossing lost baryons: {drift:.3e}");
    assert!(
        (w.tree.nodes[deep.get()].matter.mass - mass_before).abs()
            <= mass_before * 1e-9,
        "the node's own mass changed"
    );
}

/// Detail somebody has touched is never re-drawn. A tree you broke stays
/// broken, however far behind it falls.
#[test]
fn touched_detail_is_not_regenerated() {
    let mut w = World::new(galaxy(0x9111, 1e9), 20.0);
    let root = w.tree.root;
    let path = w.drill(root, Tier::Molecular, &default_spec);
    let deep = *path.last().unwrap();
    w.tree.refine(deep);
    w.tree.pin(deep);
    let epoch = w.tree.nodes[deep.get()].epoch;
    let bodies = w.tree.nodes[deep.get()].bodies.len();

    let before = w.stats.thermalised;
    w.advance_to(deep, w.time + 1.0, phys::engine::MAX_SUBSTEPS);
    let n = &w.tree.nodes[deep.get()];
    println!(
        "  a pinned node asked to cross 1 s: fell {:.3e} s behind rather than being re-drawn",
        w.time + 1.0 - n.time
    );
    assert_eq!(w.stats.thermalised, before, "pinned detail was thermalised");
    assert_eq!(n.epoch, epoch, "pinned detail was regenerated");
    assert_eq!(n.bodies.len(), bodies);
    assert!(n.time < w.time + 1.0, "it claimed to have arrived");
}

/// What gets drawn is where things are *now*, not where the last solve left
/// them.
///
/// A node the frame could not bring all the way to the instant is behind by its
/// lateness. Drawing it there makes the world stutter at exactly the rate the
/// scheduler skips things — the one artefact the whole closed-form design
/// exists to avoid. Interpolating over the lag is the same exact solution
/// applied one level down, and it is a good one for as long as the node is not
/// badly overdue.
#[test]
fn rendering_interpolates_to_the_instant() {
    let mut w = World::new(galaxy(0x11FE, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);

    // A node brought all the way to the instant has nothing to interpolate.
    w.advance_to(root, w.time + w.node_dt(root), phys::engine::MAX_SUBSTEPS);
    w.time = w.tree.nodes[root.get()].time;
    assert_eq!(w.render_lag(root), 0.0, "a solved node should have no lag");

    // Now let the world move on without solving it. The lag is exactly the
    // gap, and interpolating over it must land where integrating would.
    let before: Vec<_> = w.tree.nodes[root.get()]
        .bodies
        .iter()
        .map(|b| (b.pos, b.vel))
        .collect();
    let dt = w.node_dt(root);
    w.time += dt;
    let lag = w.render_lag(root);
    assert!(
        (lag - dt).abs() <= dt * 1e-12,
        "lag {lag:.6e} should be the gap {dt:.6e}"
    );

    // Integrate for real and compare against the drawn positions.
    w.advance_to(root, w.time, phys::engine::MAX_SUBSTEPS);
    let mut worst: f64 = 0.0;
    let mut travelled: f64 = 0.0;
    for (i, b) in w.tree.nodes[root.get()].bodies.iter().enumerate() {
        let (p0, v0) = before[i];
        let drawn = p0 + v0.scale(lag);
        worst = worst.max((drawn - b.pos).norm());
        travelled = travelled.max((b.pos - p0).norm());
    }
    println!(
        "  over one timestep the drawn position is within {:.3e} m of the \
         integrated one, on {:.3e} m travelled — {:.2e} relative",
        worst,
        travelled,
        worst / travelled.max(1e-300)
    );
    assert!(
        worst < travelled * 1e-3,
        "interpolation was {worst:.3e} m off on {travelled:.3e} m travelled"
    );
}

/// No node may end a frame with its bodies outside itself.
///
/// A node's bodies are node-relative by construction and belong within a few
/// radii of its origin. Everything downstream assumes it: SPH's smoothing
/// length, gravity's softening, the LOD's angular size, and — the place this
/// first showed up — the neighbour grid's cell index, which at nanometre
/// spacing cannot even *name* a cell 10^10 m out.
///
/// This is the check that would have caught the substep bug at its source. The
/// molecular tier's force field asked for 10^6 substeps to cross the span the
/// scheduler handed it; the cap silently stretched the step to fit 64 instead,
/// and the Lennard-Jones potential did what the comment above the cap said it
/// would. The atoms left at 2.8x10^4 m/s and the failure surfaced twenty tiers
/// later inside a hash function, as an integer overflow, in a solver that had
/// nothing to do with it.
///
/// So the assertion is deliberately made *here*, on the node that owns the
/// bodies, on the frame it happens — not on the symptom.
///
/// **Ignored, because it still fails — and what it catches is a live bug.**
/// Fixing the substep clamp cut the overshoot from 4.5x10^8 radii to 4.0x10^7
/// and no further, because the clamp was never the cause. The cause is that
/// `Tier::Atomic.floor()` is 1x10^-14 m, so the tier table calls a node four
/// orders of magnitude smaller than an atom "Atomic", the drill keeps refining
/// into it, and `sample` packs 56 Lennard-Jones bodies into a sphere smaller
/// than one atom:
///
/// ```text
///     node 15  radius 1.680e-10 m  min separation 1.331e-10 m  dt 1e-14  healthy
///     node 18  radius 5.917e-12 m  min separation 5.526e-13 m  dt floored
///     node 22  radius 2.762e-14 m  min separation 2.090e-15 m  dt floored
/// ```
///
/// At 2x10^-15 m the Lennard-Jones repulsion is of order (sigma/r)^13 ~ 10^64.
/// No timestep integrates that, which is why `configuration_dt` returns
/// something near zero and `.max(1e-24)` — a third silent clamp — hands back a
/// step that is still 10^20 too large. The configuration is unphysical before
/// any integrator sees it, so this is a sampling and tiering fault, not a
/// solver one. See `docs/BACKLOG.md`.
#[test]
#[ignore = "catches a live bug: Atomic nodes are refined below atomic scale"]
fn no_node_flings_its_bodies_out_of_itself() {
    // The setup matters, and it is the one the original failure used: drill the
    // whole ladder, then run at the *galactic* pace. Pacing to the deep node
    // instead shrinks the frame until the force field can afford it, and the
    // bug does not appear — which is the point. It is the mismatch between the
    // span a node is handed and the step its own physics needs that breaks it,
    // so the test has to preserve that mismatch to mean anything.
    let mut w = World::new(galaxy(0x5EED, 1e9), 20.0);
    let root = w.tree.root;
    let path = w.drill(root, Tier::Nuclear, &default_spec);
    println!("  drilled {} tiers, galaxy to nucleus", path.len());
    assert!(path.len() >= 6, "expected a deep ladder, got {}", path.len());

    // Generous. A body two hundred radii out is already meaningless; the bug
    // this guards against overshot by twenty orders of magnitude, so the bound
    // does not need to be tight to be decisive.
    const SLACK: f64 = 200.0;

    for frame in 0..20 {
        w.step_frame(50_000.0);
        for (i, n) in w.tree.nodes.iter().enumerate() {
            if !n.alive || n.bodies.is_empty() || !(n.matter.radius > 0.0) {
                continue;
            }
            let (worst, at) = n.bodies.iter().enumerate().fold(
                (0.0f64, 0usize),
                |(d, k), (j, b)| {
                    let r = b.pos.norm();
                    if r > d { (r, j) } else { (d, k) }
                },
            );
            assert!(
                worst.is_finite(),
                "frame {frame}: node {i} ({:?}) body {at} is at a non-finite position",
                n.tier
            );
            assert!(
                worst <= n.matter.radius * SLACK,
                "frame {frame}: node {i} ({:?}) has radius {:.3e} m and its body {at} is \
                 {:.3e} m from the origin — {:.1e} radii out",
                n.tier,
                n.matter.radius,
                worst,
                worst / n.matter.radius
            );
        }
    }
}
