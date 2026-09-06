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

fn ok(r: Result<Scene, phys::wire::WireError>) -> Scene {
    match r {
        Ok(s) => s,
        Err(e) => panic!("scene decode failed: {e}"),
    }
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

    let few = w.render(&ViewRequest { node: root, max_bodies: 100, trail: 4 });
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
    assert_eq!(back.bodies, scene.bodies, "a body changed crossing the wire");

    let per_body = bytes.len() as f64 / scene.bodies.len().max(1) as f64;
    println!(
        "  {} bodies in {} bytes — {per_body:.1} bytes each, the cost of one client frame",
        scene.bodies.len(),
        bytes.len()
    );
    assert_eq!(encode(&back), bytes, "re-encoding gave different bytes");
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
