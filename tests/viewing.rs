//! Watching a test, rather than only asserting about it.
//!
//! `docs/VIEWING.md` measured what was missing and this is the first of it. The
//! assertions here are about the *instrument*: that a node with no topology
//! draws at all, that two promoted children land in two places, and that a
//! colour is a measurement rather than a label. The images themselves are never
//! asserted on — an image assertion is the archetype of a test that cannot
//! fail.

use phys::engine::{default_spec, galaxy, World};
use phys::film;
use phys::math::v3;
use phys::render::{Film, Paint, Shot};
use phys::state::{Composition, Matter};
use phys::units::Tier;

/// A node whose contents are loose draws its contents.
///
/// The first gap `VIEWING.md` measures: `draw_structure` iterated
/// `bodies.len().min(topo.base.len())`, so a node with no topology produced an
/// empty sky and everything Phase 1 built was in that class. Verified against
/// the defect by drawing the same scene through `draw_structure`, which takes
/// the topology path and is what the engine had.
#[test]
fn a_node_with_no_topology_draws_its_contents() {
    let mut w = World::new(galaxy(0x51DE, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.nodes[root.get()].spec.count = 256;
    w.tree.refine(root);
    assert!(
        w.tree.nodes[root.get()].topology.is_none(),
        "a sampled cloud has no topology, which is the whole point of this test"
    );

    let scene = film::of_node(&w, root);
    assert!(scene.bodies.len() > 100, "{} bodies", scene.bodies.len());
    assert!(scene.topology.is_none());

    let shot = film::framing(&w, root, 0.6, 0.35);
    let (_canvas, drawn) = film::shoot(&w, root, &shot);
    assert!(
        drawn.covered > 0.02 && drawn.covered < 0.90,
        "a cloud of {} bodies covered {:.3} of the frame",
        scene.bodies.len(),
        drawn.covered
    );
    assert_eq!(drawn.parts, scene.bodies.len(), "and every one of them was drawn");

    // The defect this catches, run on the same scene: the engine's own path
    // bounded the loop by the topology's length, so a node with none drew an
    // empty sky whatever it was holding.
    let mut old = shot.canvas();
    phys::render::draw_structure(
        &mut old,
        &shot.camera,
        &scene.bodies,
        &phys::topology::Topology::default(),
        &scene.intact,
        &shot.style,
    );
    assert_eq!(
        old.coverage(),
        0.0,
        "before this, {} bodies drew nothing at all",
        scene.bodies.len()
    );
}

/// Two promoted children are two frames, and they land in two places.
///
/// `VIEWING.md`'s second gap. Every position handed to the renderer is in one
/// node's own coordinates, and the offset between two of them is
/// `Tree::separation` — which also carries the error bound that says when a
/// composed scene is below its own precision.
#[test]
fn two_promoted_children_are_placed_in_one_scene() {
    let mut w = World::new(galaxy(0xC0FFEE, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.nodes[root.get()].spec.count = 8;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    let a = w.tree.promote(root, 0, default_spec(tier.finer()));
    let b = w.tree.promote(root, 1, default_spec(tier.finer()));
    for (n, off) in [(a, v3(-1.0e19, 0.0, 0.0)), (b, v3(1.0e19, 0.0, 0.0))] {
        w.tree.nodes[n.get()].motion.offset = off;
        w.tree.nodes[n.get()].matter =
            Matter::neutral(1.0e30, 2.0e18, 100.0, Composition::primordial());
        w.tree.nodes[n.get()].spec.count = 32;
        w.tree.refine(n);
    }

    let scene = film::of_node(&w, root);
    // Each child's bodies replaced its stand-in rather than joining it: 8 slots
    // with two of them expanded to 32.
    assert_eq!(scene.bodies.len(), 6 + 32 + 32, "{} bodies", scene.bodies.len());

    // Each child's contents are clustered about where the child is, rather
    // than all at the origin, which is what composing the frames buys.
    let near = |c: phys::math::Vec3| {
        scene.bodies.iter().filter(|b| (b.pos - c).norm() < 1.0e19).count()
    };
    assert!(
        near(v3(-1.0e19, 0.0, 0.0)) >= 32 && near(v3(1.0e19, 0.0, 0.0)) >= 32,
        "the two children landed apart: {} and {}",
        near(v3(-1.0e19, 0.0, 0.0)),
        near(v3(1.0e19, 0.0, 0.0))
    );
    assert!(
        scene.err < 1.0e14,
        "and the placement is far above its own round-off: {} m against 10^19",
        scene.err
    );
}

/// A colour is read off the thing, not supplied with it.
///
/// `VIEWING.md`: *a diagnostic image should be a measurement, not a legend.*
/// The range the ramp was read against comes back with the frame, so a colour
/// can be turned into a number — which is the whole difference between this and
/// a key in the corner of the picture. The scene is the one `VIEWING.md` names:
/// a hot node beside a cold one, which had nothing to show before this.
#[test]
fn a_paint_reads_a_measurement_off_the_bodies() {
    let mut w = World::new(galaxy(0x7E45, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.nodes[root.get()].spec.count = 8;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    for (slot, off, t) in [
        (0usize, v3(-1.0e19, 0.0, 0.0), 6000.0),
        (1, v3(1.0e19, 0.0, 0.0), 50.0),
    ] {
        let n = w.tree.promote(root, slot, default_spec(tier.finer()));
        w.tree.nodes[n.get()].matter =
            Matter::neutral(1.0e30, 2.0e18, t, Composition::primordial());
        w.tree.nodes[n.get()].motion.offset = off;
    }

    let scene = film::of_node(&w, root);
    let (lo, hi) = scene
        .bodies
        .iter()
        .fold((f64::INFINITY, 0.0f64), |(l, h), b| (l.min(b.temperature), h.max(b.temperature)));
    assert!(hi / lo > 10.0, "there is a gradient to see: {lo} K to {hi} K");

    let shot = film::framing(&w, root, 0.6, 0.35).painted(Paint::Temperature { range: None });
    let (_c, drawn) = film::shoot(&w, root, &shot);
    let range = drawn.range.expect("a temperature paint reports the range it used");
    assert!(
        (range.0 - lo).abs() < 1e-9 * lo.max(1.0) && (range.1 - hi).abs() < 1e-9 * hi.max(1.0),
        "the ramp should span the frame's own extremes, {range:?} against ({lo}, {hi})"
    );

    // A role paint has no measurement and says so.
    let plain = film::framing(&w, root, 0.6, 0.35);
    assert_eq!(film::shoot(&w, root, &plain).1.range, None);
}

/// Off unless asked, and a numbered sequence when it is.
#[test]
fn a_film_is_off_unless_the_environment_asks_for_it() {
    // Deliberately not unset-and-restored around a threaded test runner: the
    // variable is read once per call and the two cases are tested by two paths
    // rather than by mutating the process.
    let dir = std::env::temp_dir().join(format!("phys-film-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);

    let mut w = World::new(galaxy(0xF11, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.nodes[root.get()].spec.count = 32;
    w.tree.refine(root);

    if std::env::var_os("PHYS_FILM").is_none() {
        assert!(Film::open("nothing").is_none(), "no PHYS_FILM, no film");
        assert!(!Film::recording());
        let mut none = film::NodeFilm::open("nothing", root, 1.0e20, Paint::Role);
        assert!(none.take(&w).is_none(), "and taking a frame is a no-op");
        assert_eq!(none.frames(), 0);
    }

    // The writing half, with the directory supplied rather than the variable,
    // so the two halves of the test cannot interfere.
    std::fs::create_dir_all(&dir).unwrap();
    let shot = Shot::framing(phys::math::Vec3::ZERO, 1.0e20, 0.6, 0.35).sized(160, 120);
    let (canvas, _) = film::shoot(&w, root, &shot);
    let path = dir.join("frame-0000.png");
    phys::render::write_png(&canvas, &path.to_string_lossy()).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    assert_eq!(&bytes[..8], &[137, 80, 78, 71, 13, 10, 26, 10], "a PNG signature");
    assert!(bytes.len() > 160 * 120, "and a frame's worth of pixels: {}", bytes.len());
    let _ = std::fs::remove_dir_all(&dir);
}

/// The fog follows the subject, whatever size the subject is.
///
/// `VIEWING.md`'s third gap: `Style::fog` was the constant 60 m, so at six
/// metres everything was unfogged and at a kilometre everything was fog.
#[test]
fn the_depth_cue_scales_with_what_is_being_framed() {
    let small = Shot::framing(phys::math::Vec3::ZERO, 6.0, 0.6, 0.35);
    let large = Shot::framing(phys::math::Vec3::ZERO, 1000.0, 0.6, 0.35);
    assert!(
        (small.style.fog / 6.0 - large.style.fog / 1000.0).abs() < 1e-12,
        "fog is a multiple of the framing radius: {} and {}",
        small.style.fog,
        large.style.fog
    );
    assert!(small.style.fog > 6.0 && small.style.fog < 60.0);
}

/// A planet drawn from orbit is one disc, because nobody has resolved it.
///
/// The case that makes an accretion film watchable from its first frame: a node
/// holding no bodies is not an empty sky, it is one body of its own stated
/// size.
#[test]
fn an_unresolved_node_draws_as_itself() {
    let mut w = World::new(galaxy(0xBA11, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let planet = w.tree.promote(root, 0, default_spec(Tier::Planetary));
    w.tree.nodes[planet.get()].matter =
        Matter::neutral(5.972e24, 6.371e6, 290.0, Composition::primordial());

    let scene = film::of_node(&w, planet);
    assert_eq!(scene.bodies.len(), 1, "one disc, not an empty sky");
    assert_eq!(scene.bodies[0].radius, 6.371e6);
    let shot = film::framing(&w, planet, 0.6, 0.35);
    let (_c, drawn) = film::shoot(&w, planet, &shot);
    assert!(drawn.covered > 0.05, "and it covers the frame: {:.3}", drawn.covered);
}
