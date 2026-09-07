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
        .bodies()
        .iter()
        .flat_map(|b| b.pos)
        .fold(0.0f32, |a, c| a.max(c.abs()));
    2.0 * extent.max(1e-6) / (s.bodies().len().max(1) as f32).cbrt()
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
    assert!(scene.bodies().is_empty(), "an unresolved node produced bodies");
    assert!(!scene.node().materialised);
    // But the node's own facts are still there — a client can say what it is
    // looking at without anything having been built.
    assert!(scene.node().mass > 0.0, "the matter should still describe itself");
    assert_eq!(scene.node().radius, w.tree.nodes[root.get()].matter.radius);
    println!(
        "  unresolved {} node: {:.3e} kg, {:.3e} m, 0 bodies",
        tier_name(scene.node().tier),
        scene.node().mass,
        scene.node().radius
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
        if scene.bodies().is_empty() {
            continue;
        }
        seen += 1;
        let reach = scene
            .bodies()
            .iter()
            .map(|b| (b.pos[0] * b.pos[0] + b.pos[1] * b.pos[1] + b.pos[2] * b.pos[2]).sqrt())
            .fold(0.0f32, f32::max);
        assert!(reach.is_finite(), "non-finite position at {}", tier_name(scene.node().tier));
        assert!(
            reach < 40.0,
            "{} bodies reach {reach} node radii — not node-relative",
            tier_name(scene.node().tier)
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
    assert!(all.bodies().len() > 500);

    let few = w.render(&ViewRequest { node: root, max_bodies: 100, trail: 4, since: f64::NEG_INFINITY, ..Default::default() });
    assert!(few.bodies().len() <= 110, "budget overshot: {}", few.bodies().len());
    assert!(few.bodies().len() >= 90, "budget undershot: {}", few.bodies().len());
    assert_eq!(few.node().body_count as usize, all.bodies().len(), "the true count should still travel");

    // Compare *medians*, not maxima. A stride sample legitimately misses the
    // rare outermost bodies of a centrally-concentrated profile, so the largest
    // radius is a noisy statistic; the median is what says whether the sample
    // has the same shape as the population. A prefix of a disc would not.
    let median = |s: &Scene| {
        let mut r: Vec<f32> = s
            .bodies()
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
        all.bodies().len(),
        few.bodies().len()
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
    assert!(scene.node().lag >= 0.0, "negative lag");

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
    let element = w.tree.nodes[root.get()].matter.radius / elements.cbrt();
    let moved = w
        .tree
        .nodes[root.get()]
        .bodies
        .iter()
        .map(|b| b.vel.norm() * scene.node().lag)
        .fold(0.0f64, f64::max);
    println!(
        "  instant {:.6e} s, solved {:.3e} s ago = {lateness:.3} cadences; \
         the fastest body moved {:.3} resolution elements",
        scene.instant,
        scene.node().lag,
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
    assert!(!scene.bodies().is_empty(), "nothing to send, so nothing is being tested");
    let bytes = encode(&scene);
    let back = ok(decode(&bytes));

    assert_eq!(back.instant.to_bits(), scene.instant.to_bits());
    assert_eq!(back.node().key, scene.node().key);
    assert_eq!(back.node().tier, scene.node().tier);
    assert_eq!(back.node().solver, scene.node().solver);
    assert_eq!(back.node().mass.to_bits(), scene.node().mass.to_bits());
    assert_eq!(back.node().radius.to_bits(), scene.node().radius.to_bits());
    assert_eq!(back.node().cadence.to_bits(), scene.node().cadence.to_bits());
    assert_eq!(back.node().mixing_time.to_bits(), scene.node().mixing_time.to_bits());
    assert_eq!(back.node().body_count, scene.node().body_count);
    assert_eq!(back.world.time.to_bits(), scene.world.time.to_bits());
    assert_eq!(back.world.live_nodes, scene.world.live_nodes);
    assert_eq!(back.world.coasted, scene.world.coasted);
    assert_eq!(back.trail.len(), scene.trail.len());
    for (a, b) in scene.trail.iter().zip(back.trail.iter()) {
        assert_eq!(a.tier, b.tier);
        assert_eq!(a.radius.to_bits(), b.radius.to_bits());
    }
    assert_eq!(back.bodies().len(), scene.bodies().len());
    // Not equality: the wire is quantised on purpose, and
    // `quantisation_is_finer_than_the_physics` is where that error is bounded.
    // Here it is enough that every body came back as recognisably itself.
    let element = spacing(&scene);
    for (a, b) in scene.bodies().iter().zip(back.bodies().iter()) {
        assert_eq!(a.kind, b.kind, "a body changed kind crossing the wire");
        for k in 0..3 {
            assert!(
                (a.pos[k] - b.pos[k]).abs() < element * 0.01,
                "a body moved crossing the wire"
            );
        }
    }

    let per_body = bytes.len() as f64 / scene.bodies().len().max(1) as f64;
    println!(
        "  {} bodies in {} bytes — {per_body:.1} bytes each, the cost of one client frame",
        scene.bodies().len(),
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
    for (a, b) in back.bodies().iter().zip(twice.bodies().iter()) {
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

    assert_eq!(tier_name(scene.node().tier), "galactic");
    assert_eq!(solver_name(scene.node().solver), "Barnes-Hut");
    let kinds: std::collections::BTreeSet<&str> =
        scene.bodies().iter().map(|b| kind_name(b.kind)).collect();
    println!(
        "  a {} node solved by {}, holding: {}",
        tier_name(scene.node().tier),
        solver_name(scene.node().solver),
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

    let n = scene.bodies().len();
    let element = spacing(&scene);
    let mut worst = 0.0f32;
    for (a, b) in scene.bodies().iter().zip(back.bodies().iter()) {
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
    for (a, b) in scene.bodies().iter().zip(back.bodies().iter()) {
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
    let quiet = w.render(&ViewRequest { node: root, max_bodies: 0, trail: 16, since: w.time, ..Default::default() });
    let quiet_bytes = encode(&quiet);

    println!(
        "  {} bodies: {} bytes when they are new, {} bytes when the client already \
         has them — {:.0}x less",
        full.len() / 17,
        full.len(),
        quiet_bytes.len(),
        full.len() as f64 / quiet_bytes.len() as f64
    );
    assert!(!quiet.bodies_included(), "an unchanged node sent its bodies anyway");
    assert!(quiet.bodies().is_empty());
    assert!(quiet_bytes.len() < full.len() / 50, "the quiet answer is not small enough");

    // But the facts still come, so the client can still say what it is looking
    // at and how far behind it is.
    let back = ok(decode(&quiet_bytes));
    assert_eq!(back.node().body_count, w.tree.nodes[root.get()].bodies.len() as u32);
    assert!(back.node().mass > 0.0);
    assert!(!back.bodies_included());
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
    let scene = w.render(&ViewRequest { node: root, max_bodies: 0, trail: 16, since: mark, ..Default::default() });
    println!("  re-solved node sent {} bodies", scene.bodies().len());
    assert!(scene.bodies_included());
    assert!(!scene.bodies().is_empty());
}

/// A copy of a scene with every node's detail removed, so what is left is the
/// fixed cost of describing it.
fn strip_bodies(s: &Scene) -> Scene {
    let mut out = s.clone();
    for n in out.nodes.iter_mut() {
        n.detail = phys::view::Detail::Unchanged;
    }
    out
}

/// The per-body cost on the wire, stated so a regression is visible.
#[test]
fn the_wire_cost_per_body_is_what_we_think() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    let scene = w.render(&ViewRequest::of(root));
    let bytes = encode(&scene);
    let overhead = encode(&strip_bodies(&scene)).len();
    let per_body = (bytes.len() - overhead) as f64 / scene.bodies().len() as f64;
    println!(
        "  {} bodies, {} bytes total, {overhead} of header — {per_body:.1} bytes per body",
        scene.bodies().len(),
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
/// `refine` samples bodies from the matter *as it currently is*, so they are
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
        .map(|b| b.pos.norm() / node.matter.radius)
        .fold(0.0f64, f64::max);
    let stale = w.time - node.last_solved;

    let scene = w.render(&ViewRequest::of(deep));
    let drawn = scene.bodies().iter().flat_map(|b| b.pos).fold(0.0f32, |a, c| a.max(c.abs()));

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
        scene.node().lag <= w.node_cadence(deep) * 1.000001,
        "render lag {:.3e} exceeded the cadence {:.3e}",
        scene.node().lag,
        w.node_cadence(deep)
    );
}

// ---------------------------------------------------------------------------
// recipes: detail as instructions
// ---------------------------------------------------------------------------

/// The whole claim, checked at the byte level: what the client builds from a
/// recipe is what the server would have sent, and the recipe is orders of
/// magnitude smaller.
#[test]
fn a_recipe_builds_exactly_what_the_server_would_have_sent() {
    let mut w = a_world();
    let root = w.tree.root;
    let deep = *w.drill(root, Tier::Stellar, &default_spec).last().unwrap();
    assert!(
        !w.tree.nodes[deep.get()].is_materialised(),
        "the node is already materialised, so this proves nothing about recipes"
    );

    // What the client is told.
    let scene = w.render(&ViewRequest::of(deep));
    let phys::view::Detail::Recipe(r) = &scene.nodes[0].detail else {
        panic!("an unmaterialised node should come back as a recipe, got {:?}", scene.nodes[0].detail);
    };
    let recipe_bytes = encode(&scene).len();
    let from_recipe = r.build().expect("the client could not build its own recipe");

    // What the server would have sent, had it paid to materialise.
    w.tree.refine(deep);
    let explicit = w.render(&ViewRequest::of(deep));
    let explicit_bytes = encode(&explicit).len();
    let truth = explicit.bodies();

    println!(
        "  {} bodies: {recipe_bytes} bytes as a recipe, {explicit_bytes} as bodies — {:.0}x",
        truth.len(),
        explicit_bytes as f64 / recipe_bytes as f64
    );

    assert_eq!(from_recipe.len(), truth.len(), "the recipe built a different number of bodies");
    // Not "close": the same sampler, the same seed, the same epoch. The only
    // difference between the two paths is quantisation, and the recipe side
    // has not been through it.
    for (i, (a, b)) in from_recipe.iter().zip(truth).enumerate() {
        let d = (0..3).map(|k| (a.pos[k] - b.pos[k]).abs()).fold(0.0f32, f32::max);
        assert!(
            d < 1e-3,
            "body {i} differs by {d} node radii: {:?} vs {:?}",
            a.pos,
            b.pos
        );
        assert_eq!(a.kind, b.kind);
    }
    assert!(
        recipe_bytes * 20 < explicit_bytes,
        "a recipe should be at least twenty times smaller; got {recipe_bytes} vs {explicit_bytes}"
    );
}

/// A recipe survives the wire and is still buildable on the far side.
#[test]
fn a_recipe_crosses_the_wire() {
    let mut w = a_world();
    let root = w.tree.root;
    let deep = *w.drill(root, Tier::Stellar, &default_spec).last().unwrap();

    let sent = w.render(&ViewRequest::of(deep));
    let mut got = ok(decode(&encode(&sent)));
    assert!(matches!(got.nodes[0].detail, phys::view::Detail::Recipe(_)));

    let n = got.materialise().expect("a decoded recipe should still build");
    println!("  built {n} specks on the far side of the wire");
    assert!(n > 0);
    assert!(matches!(got.nodes[0].detail, phys::view::Detail::Explicit(_)));

    let here = match &sent.nodes[0].detail {
        phys::view::Detail::Recipe(r) => r.build().unwrap(),
        _ => unreachable!(),
    };
    assert_eq!(here, *got.nodes[0].detail.specks(), "the wire changed the recipe");
}

/// A mangled recipe is refused rather than drawn.
#[test]
fn a_corrupt_recipe_is_refused() {
    let mut w = a_world();
    let root = w.tree.root;
    let deep = *w.drill(root, Tier::Stellar, &default_spec).last().unwrap();
    let mut scene = w.render(&ViewRequest::of(deep));

    let phys::view::Detail::Recipe(r) = &mut scene.nodes[0].detail else {
        panic!("expected a recipe");
    };
    // The checksum is over the blob, so a lie about the epoch gets through it
    // and has to be caught by the mass check instead — which is the point of
    // having both.
    r.epoch = r.epoch.wrapping_add(7);
    let bumped = scene.clone().materialise();

    let phys::view::Detail::Recipe(r) = &mut scene.nodes[0].detail else { unreachable!() };
    r.checksum ^= 1;
    let err = scene.materialise().expect_err("a mangled blob must not build");
    println!("  epoch bumped: {bumped:?}\n  checksum flipped: {err}");
    assert!(matches!(err, phys::wire::WireError::BadRecipe { .. }));
    // The scene is left alone on failure: half-expanded is worse than not.
    assert!(matches!(scene.nodes[0].detail, phys::view::Detail::Recipe(_)));
}

/// A client that cannot run the sampler is not sent a recipe it cannot use.
#[test]
fn a_client_without_a_sampler_gets_no_recipes() {
    let mut w = a_world();
    let root = w.tree.root;
    let deep = *w.drill(root, Tier::Stellar, &default_spec).last().unwrap();
    let scene = w.render(&ViewRequest { allow_recipes: false, ..ViewRequest::of(deep) });
    assert!(scene.nodes[0].detail.is_unchanged());
    assert!(!scene.bodies_included());
}

/// Pinned detail was altered by an interaction, so it is not derivable and a
/// recipe for it would be a lie.
#[test]
fn pinned_detail_is_never_sent_as_a_recipe() {
    let mut w = a_world();
    let root = w.tree.root;
    let deep = *w.drill(root, Tier::Stellar, &default_spec).last().unwrap();
    w.tree.nodes[deep.get()].pinned = true;
    let scene = w.render(&ViewRequest::of(deep));
    assert!(
        !matches!(scene.nodes[0].detail, phys::view::Detail::Recipe(_)),
        "a pinned node must not be described by a recipe"
    );
}

// ---------------------------------------------------------------------------
// volume queries and level of detail
// ---------------------------------------------------------------------------

/// A world with a node and a crowd of neighbours around it, the shape a city
/// has: many things present, few of them being solved.
fn a_neighbourhood(n: usize) -> (World, phys::ids::NodeIdx) {
    let mut w = a_world();
    let root = w.tree.root;
    let here = *w.drill(root, Tier::Stellar, &default_spec).last().unwrap();
    w.tree.refine(here);
    let have = w.tree.nodes[here.get()].bodies.len();
    let spec = default_spec(w.tree.nodes[here.get()].tier.finer());
    let stride = (have / n.max(1)).max(1);
    for k in 0..n {
        let slot = k * stride;
        if slot >= have {
            break;
        }
        w.tree.promote(here, slot, spec);
    }
    // Pace the clock to this node, as any client watching it would. Without
    // this the world runs at its coarsest node's pace and one frame is a
    // hundred thousand years, which sends the neighbours 10^8 radii away
    // before the second frame — correct, and useless to look at.
    w.pace_to(here);
    (w, here)
}

#[test]
fn a_volume_query_returns_the_neighbourhood() {
    let (w, here) = a_neighbourhood(16);
    let one = w.render(&ViewRequest::of(here));
    assert_eq!(one.nodes.len(), 1, "without a volume, a request is about one node");

    let many = w.render(&ViewRequest::within(here, [0.0, 0.0, -6.0], 40.0));
    println!("  {} nodes gathered, frame node {:?}", many.nodes.len(), many.nodes[0].facts.key);
    assert!(many.nodes.len() > 8, "the neighbours should have been gathered");
    assert_eq!(many.nodes[0].facts.key, one.nodes[0].facts.key, "the frame node comes first");
    assert_eq!(many.nodes[0].offset, [0.0; 3], "the frame node is the origin of the frame");
    assert_eq!(many.nodes[0].scale, 1.0);

    // Everything is placed in the frame node's radii, and largest-on-screen
    // first, which is the order a cap has to cut from the back of.
    let angles: Vec<f32> = many.nodes[1..].iter().map(|n| n.angular_size).collect();
    assert!(
        angles.windows(2).all(|p| p[0] >= p[1]),
        "neighbours should come back largest first, got {angles:?}"
    );
}

#[test]
fn detail_follows_angular_size() {
    let (w, here) = a_neighbourhood(16);
    let eye = [0.0, 0.0, -6.0];

    let coarse = w.render(&ViewRequest {
        volume: Some(phys::view::Volume { eye, reach: 40.0, detail_angle: 1.0, max_nodes: 256 }),
        ..ViewRequest::of(here)
    });
    let fine = w.render(&ViewRequest {
        volume: Some(phys::view::Volume { eye, reach: 40.0, detail_angle: 0.0, max_nodes: 256 }),
        ..ViewRequest::of(here)
    });

    let detailed = |s: &Scene| s.nodes.iter().filter(|n| !n.detail.is_unchanged()).count();
    // What the *neighbourhood* costs. The frame node was asked for by name and
    // is detailed either way, so leaving it in would drown the thing being
    // measured — it holds four thousand bodies and the neighbours hold none.
    let around = |s: &Scene| {
        let mut c = s.clone();
        c.nodes.remove(0);
        encode(&c).len()
    };
    // What is left after LOD has cut the detail is the facts, and a stateless
    // client pays for those every frame. A `Client` stops paying twice, so the
    // third number is what the two mechanisms are worth together.
    let mut c = phys::view::Client::new();
    let req = ViewRequest {
        volume: Some(phys::view::Volume { eye, reach: 40.0, detail_angle: 1.0, max_nodes: 256 }),
        ..ViewRequest::of(here)
    };
    c.frame(&w, &req);
    let warm = c.frame(&w, &req);

    println!(
        "  at a 1-radian threshold {} of {} nodes are detailed; at zero, {} of {}\n  \
         the neighbourhood costs {} bytes with LOD against {} without,\n  \
         and {} once the client has already been told the facts",
        detailed(&coarse),
        coarse.nodes.len(),
        detailed(&fine),
        fine.nodes.len(),
        around(&coarse),
        around(&fine),
        around(&warm)
    );
    assert_eq!(coarse.nodes.len(), fine.nodes.len(), "LOD culls detail, not nodes");
    assert_eq!(detailed(&coarse), 1, "only the node asked for by name survives a 1-radian cut");
    assert_eq!(detailed(&fine), fine.nodes.len());
    assert!(
        around(&coarse) * 3 < around(&fine),
        "LOD must make the neighbourhood much cheaper, or it is not worth having"
    );
    assert!(
        around(&warm) * 3 < around(&coarse),
        "and once the facts are known there should be almost nothing left"
    );
}

#[test]
fn a_volume_query_is_capped() {
    let (w, here) = a_neighbourhood(24);
    let all = w.render(&ViewRequest::within(here, [0.0, 0.0, -6.0], 40.0));
    let capped = w.render(&ViewRequest {
        volume: Some(phys::view::Volume {
            eye: [0.0, 0.0, -6.0],
            reach: 40.0,
            detail_angle: 1e-3,
            max_nodes: 5,
        }),
        ..ViewRequest::of(here)
    });
    println!("  {} nodes uncapped, {} at a cap of 5", all.nodes.len(), capped.nodes.len());
    assert!(all.nodes.len() > 5);
    assert_eq!(capped.nodes.len(), 5);
    // The cap cuts from the back, so what survives is what is biggest on
    // screen — and the frame node keeps its place at the head regardless.
    assert_eq!(capped.nodes[0].facts.key, all.nodes[0].facts.key);
    for (a, b) in capped.nodes.iter().zip(all.nodes.iter()) {
        assert_eq!(a.facts.key, b.facts.key);
    }
}

#[test]
fn a_body_budget_is_shared_across_the_scene() {
    let (w, here) = a_neighbourhood(16);
    let budget = 500;
    let scene = w.render(&ViewRequest {
        max_bodies: budget,
        volume: Some(phys::view::Volume {
            eye: [0.0, 0.0, -6.0],
            reach: 40.0,
            detail_angle: 0.0,
            max_nodes: 256,
        }),
        ..ViewRequest::of(here)
    });
    let mut counts: Vec<usize> = scene
        .nodes
        .iter()
        .map(|n| match &n.detail {
            phys::view::Detail::Explicit(b) => b.len(),
            phys::view::Detail::Recipe(r) => (r.count as usize).div_ceil(r.stride.max(1) as usize),
            phys::view::Detail::Unchanged => 0,
        })
        .collect();
    let total: usize = counts.iter().sum();
    counts.sort_unstable();
    println!("  budget {budget} produced {total} specks across {} nodes", scene.nodes.len());
    assert!(counts[0] > 0, "every detailed node gets something");
    assert!(
        total < budget * 4,
        "a budget of {budget} produced {total} specks — the share-out is not bounding anything"
    );
}

// ---------------------------------------------------------------------------
// per-client deltas
// ---------------------------------------------------------------------------

#[test]
fn a_client_is_told_once_and_then_left_alone() {
    let (mut w, here) = a_neighbourhood(16);
    let req = ViewRequest::within(here, [0.0, 0.0, -6.0], 40.0);
    let mut c = phys::view::Client::new();

    let first = c.frame(&w, &req);
    let join = encode(&first).len();
    assert!(first.nodes.iter().all(|n| n.facts_included), "a new client needs every fact");
    assert!(first.nodes.iter().any(|n| !n.detail.is_unchanged()));

    // Nothing has been stepped, so nothing can have changed.
    let again = c.frame(&w, &req);
    let idle = encode(&again).len();
    println!("  joining cost {join} bytes; an idle frame costs {idle}");
    assert!(again.nodes.iter().all(|n| n.detail.is_unchanged()), "nothing changed, so nothing to send");
    assert!(again.nodes.iter().all(|n| !n.facts_included), "and the facts had already been sent");
    assert!(idle * 8 < join, "an idle frame should be a small fraction of joining");

    // A step re-solves something, and only that comes back.
    for _ in 0..6 {
        w.step_frame(50_000.0);
    }
    let moved = c.frame(&w, &req);
    let changed = moved.nodes.iter().filter(|n| !n.detail.is_unchanged()).count();
    println!("  after stepping, {changed} of {} nodes had something new", moved.nodes.len());
    assert!(changed < moved.nodes.len(), "a step should not invalidate the whole neighbourhood");
}

#[test]
fn a_client_is_told_what_to_forget() {
    let (w, here) = a_neighbourhood(16);
    let eye = [0.0, 0.0, -6.0];
    let mut c = phys::view::Client::new();

    let wide = c.frame(&w, &ViewRequest::within(here, eye, 40.0));
    let held = c.holding();
    assert_eq!(held, wide.nodes.len(), "the client holds what it was sent");

    // Walk away: a tighter reach drops most of the neighbourhood.
    let near = c.frame(&w, &ViewRequest::within(here, eye, 0.5));
    println!(
        "  held {held} nodes, then {} — {} dropped",
        c.holding(),
        near.dropped.len()
    );
    assert!(!near.dropped.is_empty(), "nodes that left the query must be named");
    assert_eq!(c.holding(), near.nodes.len());
    assert_eq!(held, near.nodes.len() + near.dropped.len());
    for d in &near.dropped {
        assert!(!near.nodes.iter().any(|n| n.facts.key == *d), "dropped and sent at once");
    }
}

/// Facts are 141 bytes and a placement is 21. For a crowd of things nobody is
/// touching, not re-sending the facts *is* the steady state.
#[test]
fn unchanged_facts_are_not_resent() {
    let (w, here) = a_neighbourhood(16);
    let req = ViewRequest::within(here, [0.0, 0.0, -6.0], 40.0);
    let mut c = phys::view::Client::new();

    let first = c.frame(&w, &req);
    let second = c.frame(&w, &req);
    println!(
        "  {} nodes: {} bytes with facts, {} without",
        second.nodes.len(),
        encode(&first).len(),
        encode(&second).len()
    );

    // And the client can put them back, which is the other half of the deal.
    let mut known = std::collections::HashMap::new();
    let mut a = ok(decode(&encode(&first)));
    let mut b = ok(decode(&encode(&second)));
    assert_eq!(a.merge_facts(&mut known), 0, "the first frame carried every fact");
    let filled = b.merge_facts(&mut known);
    println!("  {filled} of {} nodes had their facts filled in", b.nodes.len());
    assert_eq!(filled, b.nodes.len());
    for (x, y) in a.nodes.iter().zip(b.nodes.iter()) {
        assert_eq!(x.facts.key, y.facts.key);
        assert_eq!(x.facts.mass, y.facts.mass, "a filled-in fact must match the sent one");
        assert_eq!(x.facts.radius, y.facts.radius);
        assert_eq!(x.facts.epoch, y.facts.epoch);
    }
}

/// Detail generated at one epoch is *gone* at the next, not merely stale, so a
/// client holding it has to be told even though nothing was re-solved.
#[test]
fn a_new_epoch_forces_a_resend() {
    let (mut w, here) = a_neighbourhood(4);
    let req = ViewRequest::of(here);
    let mut c = phys::view::Client::new();

    assert!(!c.frame(&w, &req).nodes[0].detail.is_unchanged());
    assert!(c.frame(&w, &req).nodes[0].detail.is_unchanged());

    w.tree.nodes[here.get()].epoch += 1;
    let after = c.frame(&w, &req);
    println!("  epoch bumped; detail resent: {}", !after.nodes[0].detail.is_unchanged());
    assert!(!after.nodes[0].detail.is_unchanged(), "old-epoch detail is gone, not stale");
    assert!(after.nodes[0].facts_included, "and its facts travel with it");
}
