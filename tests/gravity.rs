//! Gravity is derived from what a thing is inside.
//!
//! `docs/PLAY.md` D6, and Phase 1's "`G_EARTH` deleted in favour of derived g".
//! `solvers::structure::G_EARTH` was `(0, 0, -9.80665)` and supplied the weight
//! for every structure and every falling piece — so debris fell at Earth's
//! surface gravity along its own structure's negative z on a moon, on a ship
//! under thrust, and in orbit, while the engine computed real gravitational
//! fields at every other tier and ignored them here.
//!
//! Nothing below knows what a planet is. `Tree::gravity_at` walks the ancestor
//! chain and adds what each one pulls with.

use phys::engine::{default_spec, galaxy, World};
use phys::math::v3;
use phys::state::{Composition, Matter};
use phys::units::{Tier, G};

const EARTH_MASS: f64 = 5.972e24;
const EARTH_RADIUS: f64 = 6.371e6;

/// A world holding one Earth-sized body, with a small node to place in it.
fn an_earth() -> (World, phys::ids::NodeIdx) {
    let mut w = World::new(galaxy(0xEA47, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let planet = w.tree.promote(root, 0, default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.spec.count = 8;
    }
    w.tree.refine(planet);
    let node = w.tree.promote(planet, 0, default_spec(Tier::Continuum));
    {
        let n = &mut w.tree.nodes[node.get()];
        n.matter = Matter::neutral(1.0, 1.0, 290.0, Composition::primordial());
    }
    (w, node)
}

fn at(w: &mut World, node: phys::ids::NodeIdx, offset: phys::math::Vec3) -> phys::math::Vec3 {
    w.tree.nodes[node.get()].motion.offset = offset;
    w.tree.gravity_at(node)
}

/// Surface gravity comes out of mass and radius, with no mention of Earth.
#[test]
fn a_body_of_earths_mass_and_radius_has_earths_surface_gravity() {
    let (mut w, node) = an_earth();
    let g = at(&mut w, node, v3(0.0, 0.0, EARTH_RADIUS));
    let expected = G * EARTH_MASS / (EARTH_RADIUS * EARTH_RADIUS);
    assert!(
        (g.norm() - expected).abs() / expected < 1e-9,
        "surface gravity came out {} m/s^2, not {expected}",
        g.norm()
    );
    assert!((expected - 9.82).abs() < 0.01, "and that should be about 9.82, not {expected}");
    // Pointing down, which here means at the centre.
    assert!(g.z < 0.0 && g.x.abs() < 1e-9 && g.y.abs() < 1e-9, "{g:?}");
}

/// Inside, the field is the mass *enclosed below*, which is Newton's shell
/// theorem rather than a guard against dividing by zero — and it is why this
/// needs no epsilon at the centre.
#[test]
fn inside_a_body_only_what_is_below_pulls() {
    let (mut w, node) = an_earth();
    let surface = at(&mut w, node, v3(0.0, 0.0, EARTH_RADIUS)).norm();
    let half = at(&mut w, node, v3(0.0, 0.0, 0.5 * EARTH_RADIUS)).norm();
    let tenth = at(&mut w, node, v3(0.0, 0.0, 0.1 * EARTH_RADIUS)).norm();
    let centre = at(&mut w, node, v3(0.0, 0.0, 0.0)).norm();

    assert!(
        (half / surface - 0.5).abs() < 1e-9,
        "at half a radius the field should be half the surface value, and is {}",
        half / surface
    );
    assert!(
        (tenth / surface - 0.1).abs() < 1e-9,
        "and linear in depth all the way down: {} at a tenth",
        tenth / surface
    );
    // At the centre the body itself pulls with nothing — but the *surroundings*
    // still do, and the honest statement is that what you feel there is what the
    // body's own centre feels. Asserting a flat zero measured 9.2e-15 m/s^2 and
    // was right to: that is the galaxy the planet is in.
    let ambient = w.tree.gravity_at(w.tree.nodes[node.get()].parent).norm();
    assert!(
        (centre - ambient).abs() <= ambient * 1e-9,
        "at the centre the field should be the surroundings alone ({ambient}), \
         and is {centre}"
    );
    assert!(
        centre < surface * 1e-12,
        "and the surroundings should be nothing beside the body: {centre} \
         against a surface {surface}"
    );
}

/// Outside, it falls off as the inverse square.
#[test]
fn outside_a_body_the_field_falls_off_as_the_inverse_square() {
    let (mut w, node) = an_earth();
    let near = at(&mut w, node, v3(0.0, 0.0, 2.0 * EARTH_RADIUS)).norm();
    let far = at(&mut w, node, v3(0.0, 0.0, 4.0 * EARTH_RADIUS)).norm();
    assert!(
        (near / far - 4.0).abs() < 1e-9,
        "doubling the distance changed the field by {}x, not 4x",
        near / far
    );
}

/// A thing does not pull on itself: it feels what is left over.
///
/// An ancestor's mass includes the node's own, because a promoted child's
/// stand-in body stays in its parent's list. For a speck on a planet that is a
/// rounding; for a node that is a serious fraction of what contains it, it is
/// the answer.
///
/// Stated at *half* the parent's mass on purpose. A node holding the whole of
/// its parent is refused by a separate guard and so passes whether or not the
/// subtraction happens, which is how the first version of this test passed with
/// the subtraction deleted.
#[test]
fn a_node_does_not_attract_itself() {
    let (mut w, node) = an_earth();
    let surface = at(&mut w, node, v3(0.0, 0.0, EARTH_RADIUS)).norm();

    w.tree.nodes[node.get()].matter.mass = 0.5 * EARTH_MASS;
    let g = at(&mut w, node, v3(0.0, 0.0, EARTH_RADIUS)).norm();
    assert!(
        (g / surface - 0.5).abs() < 1e-9,
        "a node that is half of what contains it should feel the other half — \
         {} of the surface field, not {}",
        0.5,
        g / surface
    );

    // And the whole of it feels nothing of it.
    w.tree.nodes[node.get()].matter.mass = EARTH_MASS;
    let g = at(&mut w, node, v3(0.0, 0.0, EARTH_RADIUS)).norm();
    let ambient = w.tree.gravity_at(w.tree.nodes[node.get()].parent).norm();
    assert!(
        (g - ambient).abs() <= ambient * 1e-9,
        "a node that is the whole of its parent should feel only what the \
         parent feels ({ambient}), and feels {g}"
    );
}

/// And in deep space there is essentially nothing, by the same arithmetic.
///
/// This is the case the old constant got most wrong: a structure adrift between
/// stars carried its weight at 9.80665 m/s^2.
#[test]
fn a_node_in_a_galaxy_is_in_free_fall() {
    let mut w = World::new(galaxy(0xDEE9, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.nodes[0].spec.count = 64;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    let node = w.tree.promote(root, 0, default_spec(tier.finer()));
    let g = w.tree.gravity_at(node).norm();
    assert!(
        g < 1e-6,
        "a node in interstellar space feels {g} m/s^2, where the constant it \
         replaces would have said 9.80665"
    );
    assert!(g > 0.0, "and not exactly nothing: a galaxy is still there");
}

/// The field a structure was proportioned in travels with it, because the
/// member radii depend on it and a client holds only part of the tree.
#[test]
fn the_field_a_structure_was_built_in_is_stored_on_the_node() {
    use phys::morph::{Environment, Program};
    let (mut w, node) = an_earth();
    w.tree.nodes[node.get()].motion.offset = v3(0.0, 0.0, EARTH_RADIUS);
    {
        let n = &mut w.tree.nodes[node.get()];
        n.matter = Matter::neutral(900.0, 6.0, 291.0, Program::Tree.substrate());
        n.spec.count = 600;
    }
    w.plant(node, Program::Tree, Some(Environment::default()));
    assert_eq!(
        w.tree.nodes[node.get()].gravity.norm(),
        0.0,
        "nothing has been sampled yet, so there is no field to have built it in"
    );

    w.tree.refine(node);
    let stored = w.tree.nodes[node.get()].gravity;
    assert!(
        (stored.norm() - 9.82).abs() < 0.01,
        "the node should have recorded the field it was proportioned in, and \
         recorded {stored:?}"
    );
    assert_eq!(stored, w.tree.gravity_at(node), "and it should be the derived one");
}
