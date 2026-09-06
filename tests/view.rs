//! The seam between solving and looking.
//!
//! Phase 2's gate is "the viewer never touches a solver". That is not a thing a
//! comment can enforce, so these tests hold it three ways: `render` takes
//! `&self` and so cannot change anything; a `Scene` contains no engine type and
//! so cannot be stepped; and a scene survives a round trip through the wire, so
//! the split works across a process rather than only across a function call.

use phys::engine::{default_spec, galaxy, World};
use phys::units::*;
use phys::view::{decode, encode, kind_name, solver_name, tier_name, Scene, ViewRequest};
use phys::view::Speck;

fn ok(r: Result<Scene, phys::wire::WireError>) -> Scene {
    match r {
        Ok(s) => s,
        Err(e) => panic!("scene decode failed: {e}"),
    }
}

/// Mean spacing between bodies, in node radii.
///
/// Computed from the extent the bodies *actually* occupy rather than from the
/// node's nominal radius: a sampled profile has a tail, so a node whose radius
/// says one can easily have bodies out at forty, and the spacing that matters
/// is the one between the things being drawn.
fn spacing(s: &Scene) -> f32 {
    let extent = s
        .bodies
        .iter()
        .flat_map(|b| b.pos)
        .fold(0.0f32, |a, c| a.max(c.abs()));
    2.0 * extent.max(1e-6) / (s.bodies.len().max(1) as f32).cbrt()
}

fn a_world() -> World {
    let mut w = World::new(galaxy(0x5CE7E, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 2000;
    w.time_rate = 0.05;
    w
}

/// Rendering cannot change the world. This is the property the whole phase is
/// about, and `&self` is what enforces it — the test records that the state is
/// genuinely untouched, including the parts a lazy implementation would be
/// tempted to materialise on demand.
#[test]
fn rendering_changes_nothing() {
    let mut w = a_world();
    let root = w.tree.root;
    w.step_frame(50_000.0);

    let before_bytes = w.tree.detail_bytes();
    let before_time = w.time;
    let before_frames = w.stats.frames;
    let before_nodes = w.tree.nodes.len();
    let before_epoch = w.tree.nodes[root.get()].epoch;

    for _ in 0..50 {
        let _ = w.render(&ViewRequest::of(root));
    }

    assert_eq!(w.tree.detail_bytes(), before_bytes, "rendering materialised something");
    assert_eq!(w.time.to_bits(), before_time.to_bits(), "rendering advanced the clock");
    assert_eq!(w.stats.frames, before_frames, "rendering ran a frame");
    assert_eq!(w.tree.nodes.len(), before_nodes, "rendering added nodes");
    assert_eq!(w.tree.nodes[root.get()].epoch, before_epoch, "rendering bumped an epoch");
    println!("  fifty renders, world unchanged: {before_bytes} bytes of detail, t = {before_time:.3e}");
}

/// An unresolved node draws as nothing, rather than resolving itself to be
/// drawable. Detail is the solve side's decision.
#[test]
fn an_unresolved_node_renders_empty() {
    let w = a_world();
    let root = w.tree.root;
    assert!(!w.tree.nodes[root.get()].is_materialised());

    let scene = w.render(&ViewRequest::of(root));
    assert!(scene.bodies.is_empty(), "an unresolved node produced bodies");
    assert!(!scene.node.materialised);
    // But the node's own facts are still there — a client can say what it is
    // looking at without anything having been built.
    assert!(scene.node.mass > 0.0, "the aggregate should still describe itself");
    assert_eq!(scene.node.radius, w.tree.nodes[root.get()].agg.radius);
    println!(
        "  unresolved {} node: {:.3e} kg, {:.3e} m, 0 bodies",
        tier_name(scene.node.tier),
        scene.node.mass,
        scene.node.radius
    );
}

/// Positions come back in units of the node's radius, so one renderer draws a
/// galaxy and a nucleus with the same arithmetic.
#[test]
fn positions_are_node_relative_at_every_scale() {
    let mut w = a_world();
    let root = w.tree.root;
    let path = w.drill(root, Tier::Nuclear, &default_spec);

    let mut seen = 0;
    for &node in &path {
        w.tree.refine(node);
        let scene = w.render(&ViewRequest::of(node));
        if scene.bodies.is_empty() {
            continue;
        }
        seen += 1;
        let reach = scene
            .bodies
            .iter()
            .map(|b| (b.pos[0] * b.pos[0] + b.pos[1] * b.pos[1] + b.pos[2] * b.pos[2]).sqrt())
            .fold(0.0f32, f32::max);
        assert!(reach.is_finite(), "non-finite position at {}", tier_name(scene.node.tier));
        assert!(
            reach < 40.0,
            "{} bodies reach {reach} node radii — not node-relative",
            tier_name(scene.node.tier)
        );
    }
    println!("  {seen} tiers rendered, every one of order one in node radii");
    assert!(seen >= 4, "expected several tiers to have detail, saw {seen}");
}

/// A client with a smaller budget gets a *sample*, not a prefix. The first
/// thousand parcels of a disc are all in one place.
#[test]
fn a_body_budget_samples_rather_than_truncates() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    let all = w.render(&ViewRequest::of(root));
    assert!(all.bodies.len() > 500);

    let few = w.render(&ViewRequest { node: root, max_bodies: 100, trail: 4, since: f64::NEG_INFINITY });
    assert!(few.bodies.len() <= 110, "budget overshot: {}", few.bodies.len());
    assert!(few.bodies.len() >= 90, "budget undershot: {}", few.bodies.len());
    assert_eq!(few.node.body_count as usize, all.bodies.len(), "the true count should still travel");

    // Compare *medians*, not maxima. A stride sample legitimately misses the
    // rare outermost bodies of a centrally-concentrated profile, so the largest
    // radius is a noisy statistic; the median is what says whether the sample
    // has the same shape as the population. A prefix of a disc would not.
    let median = |s: &Scene| {
        let mut r: Vec<f32> = s
            .bodies
            .iter()
            .map(|b| (b.pos[0] * b.pos[0] + b.pos[1] * b.pos[1]).sqrt())
            .collect();
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        r[r.len() / 2]
    };
    let (whole, sampled) = (median(&all), median(&few));
    let ratio = sampled / whole.max(f32::MIN_POSITIVE);
    println!(
        "  median radius: {} bodies at {whole:.3}, {} sampled at {sampled:.3} — {ratio:.3}x",
        all.bodies.len(),
        few.bodies.len()
    );
    assert!(
        (0.8..1.25).contains(&ratio),
        "the sample has a different shape from the population: {ratio:.3}x"
    );
}

/// The scene is already carried to its instant, so a client never has to know
/// that some of the world was solved earlier than the rest.
#[test]
fn a_scene_is_interpolated_to_one_instant() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    for _ in 0..3 {
        w.step_frame(50_000.0);
    }
    let scene = w.render(&ViewRequest::of(root));
    assert_eq!(scene.instant.to_bits(), w.time.to_bits(), "the scene is at another instant");
    assert!(scene.node.lag >= 0.0, "negative lag");

    // The lag can be large in seconds and is still small in the units that
    // matter, and the two are tied together by the scheduler rather than by
    // luck. A node's cadence is *defined* as the time for its fastest body to
    // cross one resolution element, so a node that is `L` cadences late has
    // bodies that have moved at most `L` resolution elements since they were
    // solved. Lateness under one therefore bounds the extrapolation this render
    // performs to under one element — which is why drawing a coasting node
    // linearly is honest rather than a smear.
    let lateness = w.lateness(root, w.time);
    let elements = w.tree.nodes[root.get()].bodies.len() as f64;
    let element = w.tree.nodes[root.get()].agg.radius / elements.cbrt();
    let moved = w
        .tree
        .nodes[root.get()]
        .bodies
        .iter()
        .map(|b| b.vel.norm() * scene.node.lag)
        .fold(0.0f64, f64::max);
    println!(
        "  instant {:.6e} s, solved {:.3e} s ago = {lateness:.3} cadences; \
         the fastest body moved {:.3} resolution elements",
        scene.instant,
        scene.node.lag,
        moved / element
    );
    assert!(lateness < 1.0, "this world was meant to be keeping up");
    assert!(
        moved <= element * lateness.max(1e-9) * 1.5,
        "extrapolation ran past what the cadence permits: {:.3} elements at {lateness:.3} cadences",
        moved / element
    );
}

/// The split has to work across a process, not just across a function call.
#[test]
fn a_scene_survives_the_wire() {
    let mut w = a_world();
    let root = w.tree.root;
    let path = w.drill(root, Tier::Planetary, &default_spec);
    let deep = *path.last().unwrap();
    w.step_frame(50_000.0);
    // Refine *after* stepping: a frame may coarsen or thermalise a node, and a
    // wire test that compared two empty body lists would pass vacuously.
    w.tree.refine(deep);
    let scene = w.render(&ViewRequest::of(deep));
    assert!(!scene.bodies.is_empty(), "nothing to send, so nothing is being tested");
    let bytes = encode(&scene);
    let back = ok(decode(&bytes));

    assert_eq!(back.instant.to_bits(), scene.instant.to_bits());
    assert_eq!(back.node.key, scene.node.key);
    assert_eq!(back.node.tier, scene.node.tier);
    assert_eq!(back.node.solver, scene.node.solver);
    assert_eq!(back.node.mass.to_bits(), scene.node.mass.to_bits());
    assert_eq!(back.node.radius.to_bits(), scene.node.radius.to_bits());
    assert_eq!(back.node.cadence.to_bits(), scene.node.cadence.to_bits());
    assert_eq!(back.node.mixing_time.to_bits(), scene.node.mixing_time.to_bits());
    assert_eq!(back.node.body_count, scene.node.body_count);
    assert_eq!(back.world.time.to_bits(), scene.world.time.to_bits());
    assert_eq!(back.world.live_nodes, scene.world.live_nodes);
    assert_eq!(back.world.coasted, scene.world.coasted);
    assert_eq!(back.trail.len(), scene.trail.len());
    for (a, b) in scene.trail.iter().zip(back.trail.iter()) {
        assert_eq!(a.tier, b.tier);
        assert_eq!(a.radius.to_bits(), b.radius.to_bits());
    }
    assert_eq!(back.bodies.len(), scene.bodies.len());
    // Not equality: the wire is quantised on purpose, and
    // `quantisation_is_finer_than_the_physics` is where that error is bounded.
    // Here it is enough that every body came back as recognisably itself.
    let element = spacing(&scene);
    for (a, b) in scene.bodies.iter().zip(back.bodies.iter()) {
        assert_eq!(a.kind, b.kind, "a body changed kind crossing the wire");
        for k in 0..3 {
            assert!(
                (a.pos[k] - b.pos[k]).abs() < element * 0.01,
                "a body moved crossing the wire"
            );
        }
    }

    let per_body = bytes.len() as f64 / scene.bodies.len().max(1) as f64;
    println!(
        "  {} bodies in {} bytes — {per_body:.1} bytes each, the cost of one client frame",
        scene.bodies.len(),
        bytes.len()
    );
    // Re-encoding a *decoded* scene is not byte-identical, and asserting that
    // it were would be asserting the wrong thing: the ranges are measured from
    // the values present, so a decoded scene — whose values sit on the
    // quantisation grid — measures a marginally smaller span and lands the grid
    // somewhere marginally different. The server always encodes from the world,
    // never from a scene it decoded, and `rendering_is_deterministic` covers
    // that path.
    //
    // What a transport does need is that decoding is *stable*: sending a scene
    // back and forth must not let it wander.
    let twice = ok(decode(&encode(&back)));
    let mut wander = 0.0f32;
    for (a, b) in back.bodies.iter().zip(twice.bodies.iter()) {
        for k in 0..3 {
            wander = wander.max((a.pos[k] - b.pos[k]).abs());
        }
    }
    println!("  a second round trip moved bodies by at most {wander:.8} node radii");
    assert!(wander < element * 0.01, "the scene wandered on a second round trip");
}

/// Rendering the same unchanged world twice gives the same bytes. Without this
/// a client cannot tell "nothing happened" from "something happened that I
/// cannot see", and neither can a diff-based transport.
#[test]
fn rendering_is_deterministic() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    w.step_frame(50_000.0);
    let a = encode(&w.render(&ViewRequest::of(root)));
    let b = encode(&w.render(&ViewRequest::of(root)));
    println!("  {} bytes, twice", a.len());
    assert_eq!(a, b);
}

/// A malformed scene is a parse error, never a panic. Scenes arrive over a
/// network in every future this architecture has.
#[test]
fn a_malformed_scene_does_not_panic() {
    let mut w = a_world();
    w.tree.refine(w.tree.root);
    let good = encode(&w.render(&ViewRequest::of(w.tree.root)));

    assert!(decode(&[]).is_err());
    assert!(decode(b"NOPE\x01\x00").is_err());
    let mut cut = 0;
    for n in (0..good.len()).step_by(13) {
        if decode(&good[..n]).is_err() {
            cut += 1;
        }
    }
    for i in (0..good.len()).step_by(211) {
        let mut bad = good.clone();
        bad[i] ^= 0xA5;
        let _ = decode(&bad);
    }
    println!("  {cut} truncations and a scatter of flips of a {}-byte scene", good.len());
    assert!(cut > 0);
}

/// A client can name what it is looking at without linking the engine. These
/// tables are the whole vocabulary it needs.
#[test]
fn a_client_can_name_things_without_the_engine() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    let scene = w.render(&ViewRequest::of(root));

    assert_eq!(tier_name(scene.node.tier), "galactic");
    assert_eq!(solver_name(scene.node.solver), "Barnes-Hut");
    let kinds: std::collections::BTreeSet<&str> =
        scene.bodies.iter().map(|b| kind_name(b.kind)).collect();
    println!(
        "  a {} node solved by {}, holding: {}",
        tier_name(scene.node.tier),
        solver_name(scene.node.solver),
        kinds.iter().cloned().collect::<Vec<_>>().join(", ")
    );
    assert!(!kinds.contains(&"?"), "an unnamed body kind reached a client");
}


/// Quantising the wire is only safe if the step is small against *the node's
/// own resolution* — the engine never claims to know where anything is more
/// precisely than one resolution element. This is that claim, measured rather
/// than argued.
#[test]
fn quantisation_is_finer_than_the_physics() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    let scene = w.render(&ViewRequest::of(root));
    let back = ok(decode(&encode(&scene)));

    let n = scene.bodies.len();
    let element = spacing(&scene);
    let mut worst = 0.0f32;
    for (a, b) in scene.bodies.iter().zip(back.bodies.iter()) {
        for k in 0..3 {
            worst = worst.max((a.pos[k] - b.pos[k]).abs());
        }
    }
    // The ratio should come out at about `n^(1/3) / 65536` regardless of how
    // spread out the node is, because both the step and the spacing scale with
    // the measured extent. That cancellation is the reason an adaptive span is
    // the right design rather than merely a bigger one.
    let predicted = (n as f32).cbrt() / 65536.0;
    println!(
        "  {n} bodies: mean spacing {element:.5} node radii, worst quantisation error \
         {worst:.7} — {:.4}% of a spacing (predicted {:.4}%)",
        100.0 * worst / element,
        100.0 * predicted
    );
    assert!(
        worst < element * 0.01,
        "quantisation error {worst:.3e} is {:.2}% of the body spacing",
        100.0 * worst / element
    );
}

/// Velocity survives too, because the client needs it to carry a coasting node
/// forward itself.
#[test]
fn velocity_survives_the_wire_well_enough_to_extrapolate() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    let scene = w.render(&ViewRequest::of(root));
    let back = ok(decode(&encode(&scene)));

    // Carry both forward by one cadence and compare where they land. That is
    // the use the velocity is actually put to, so it is the error that matters.
    let dt = w.node_cadence(root) as f32;
    let mut worst = 0.0f32;
    for (a, b) in scene.bodies.iter().zip(back.bodies.iter()) {
        let (x, y) = (a.at(dt), b.at(dt));
        for k in 0..3 {
            worst = worst.max((x[k] - y[k]).abs());
        }
    }
    let element = spacing(&scene);
    println!(
        "  after carrying one whole cadence ({dt:.3e} s), the two disagree by \
         {worst:.6} node radii — {:.2}% of an element",
        100.0 * worst / element
    );
    assert!(
        worst < element * 0.25,
        "extrapolating a decoded velocity drifted {:.2}% of an element",
        100.0 * worst / element
    );
}

/// The thing that stops bandwidth scaling with how much world is on screen: a
/// node nobody re-solved has nothing new to say, so nothing is sent.
#[test]
fn an_unchanged_node_costs_almost_nothing() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    w.step_frame(50_000.0);

    let full = encode(&w.render(&ViewRequest::of(root)));

    // Ask again as a client that already has them, at the current instant.
    let quiet = w.render(&ViewRequest { node: root, max_bodies: 0, trail: 16, since: w.time });
    let quiet_bytes = encode(&quiet);

    println!(
        "  {} bodies: {} bytes when they are new, {} bytes when the client already \
         has them — {:.0}x less",
        full.len() / 17,
        full.len(),
        quiet_bytes.len(),
        full.len() as f64 / quiet_bytes.len() as f64
    );
    assert!(!quiet.bodies_included, "an unchanged node sent its bodies anyway");
    assert!(quiet.bodies.is_empty());
    assert!(quiet_bytes.len() < full.len() / 50, "the quiet answer is not small enough");

    // But the facts still come, so the client can still say what it is looking
    // at and how far behind it is.
    let back = ok(decode(&quiet_bytes));
    assert_eq!(back.node.body_count, w.tree.nodes[root.get()].bodies.len() as u32);
    assert!(back.node.mass > 0.0);
    assert!(!back.bodies_included);
}

/// And a node that *has* been re-solved still sends. A transport that never
/// sent anything would pass the test above.
#[test]
fn a_resolved_node_still_sends() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    let mark = w.time;
    // Step until the node actually comes due and is solved.
    for _ in 0..40 {
        w.step_frame(50_000.0);
        if w.tree.nodes[root.get()].last_solved > mark {
            break;
        }
    }
    assert!(
        w.tree.nodes[root.get()].last_solved > mark,
        "the node never got solved, so this proves nothing"
    );
    let scene = w.render(&ViewRequest { node: root, max_bodies: 0, trail: 16, since: mark });
    println!("  re-solved node sent {} bodies", scene.bodies.len());
    assert!(scene.bodies_included);
    assert!(!scene.bodies.is_empty());
}

/// The per-body cost on the wire, stated so a regression is visible.
#[test]
fn the_wire_cost_per_body_is_what_we_think() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    let scene = w.render(&ViewRequest::of(root));
    let bytes = encode(&scene);
    let overhead = encode(&Scene { bodies: Vec::new(), ..scene.clone() }).len();
    let per_body = (bytes.len() - overhead) as f64 / scene.bodies.len() as f64;
    println!(
        "  {} bodies, {} bytes total, {overhead} of header — {per_body:.1} bytes per body",
        scene.bodies.len(),
        bytes.len()
    );
    assert!(
        (per_body - 17.0).abs() < 0.5,
        "per-body cost moved to {per_body:.1}; update the estimate or find the regression"
    );
    let _ = Speck::default();
}


/// A node materialised long after it was last solved must not be flung across
/// the sky.
///
/// `refine` samples bodies from the aggregate *as it currently is*, so they are
/// valid now — but it has no clock, so `last_solved` stays wherever it was. A
/// renderer that carried them by `time - last_solved` extrapolated correct
/// positions at three and a half node radii out to 10^8, and nothing caught it
/// because a scene that absurd still round-tripped and still drew.
///
/// The cap in `render_lag` is one cadence, which is by definition one
/// resolution element of travel. This is that bug, kept.
#[test]
fn a_freshly_materialised_node_is_not_extrapolated() {
    let mut w = a_world();
    let root = w.tree.root;
    let deep = *w.drill(root, Tier::Planetary, &default_spec).last().unwrap();

    // Let a lot of world time pass with the node unresolved, so `last_solved`
    // is far behind, then materialise it.
    for _ in 0..5 {
        w.step_frame(50_000.0);
    }
    w.tree.refine(deep);

    let node = &w.tree.nodes[deep.get()];
    let truth = node
        .bodies
        .iter()
        .map(|b| b.pos.norm() / node.agg.radius)
        .fold(0.0f64, f64::max);
    let stale = w.time - node.last_solved;

    let scene = w.render(&ViewRequest::of(deep));
    let drawn = scene.bodies.iter().flat_map(|b| b.pos).fold(0.0f32, |a, c| a.max(c.abs()));

    println!(
        "  last solved {stale:.3e} s ago; bodies really reach {truth:.2} node radii, \
         drawn at {drawn:.2}"
    );
    assert!(
        (drawn as f64) < truth * 3.0 + 1.0,
        "a freshly materialised node was drawn at {drawn:.3e} radii when its bodies \
         are at {truth:.3e}"
    );
    // And the cap is the cadence, so the carry is at most one element.
    assert!(
        scene.node.lag <= w.node_cadence(deep) * 1.000001,
        "render lag {:.3e} exceeded the cadence {:.3e}",
        scene.node.lag,
        w.node_cadence(deep)
    );
}
