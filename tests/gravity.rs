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
        n.matter = Matter::neutral(900.0, 6.0, 291.0, Composition::organic());
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

/// A node on the side of a planet carries its weight along its own down.
///
/// The `BACKLOG.md` entry this closes, whose stated trigger is "the first
/// oriented node, which is Phase 4's first terrain patch": `gravity_at` summed
/// each ancestor's pull in *that ancestor's* axes, so a node on the `+x` side
/// of a planet was told the field points along `-x`. True in the planet's
/// frame, and wrong for the node, which is generated with `+z` up — a structure
/// standing there would have carried its weight out sideways through its own
/// geometry.
///
/// A patch of ground is oriented by construction, so this is not an edge case
/// once there is terrain; it is every patch except the two on the axis.
#[test]
fn a_node_on_the_side_of_a_planet_is_pulled_along_its_own_down() {
    use phys::math::Quat;
    let (mut w, node) = an_earth();
    // Standing on the equator at longitude zero: the node's local +z points
    // away from the planet's centre, which is a quarter turn about +y.
    let up = Quat::from_axis_angle(v3(0.0, 1.0, 0.0), std::f64::consts::FRAC_PI_2);
    w.tree.nodes[node.get()].motion.orientation = up;
    let g = at(&mut w, node, v3(EARTH_RADIUS, 0.0, 0.0));

    let expected = G * EARTH_MASS / (EARTH_RADIUS * EARTH_RADIUS);
    assert!(
        (g.norm() - expected).abs() / expected < 1e-9,
        "the strength is unchanged by which way the node faces: {} against {expected}",
        g.norm()
    );
    assert!(
        g.z < 0.0 && g.x.abs() < 1e-9 * expected && g.y.abs() < 1e-9 * expected,
        "and it should point along the node's own -z, not the planet's -x: {g:?}"
    );

    // The unoriented node in the same place is the measurement of what this
    // was doing before: the same field, reported along the parent's -x.
    w.tree.nodes[node.get()].motion.orientation = Quat::IDENTITY;
    let parent_axes = at(&mut w, node, v3(EARTH_RADIUS, 0.0, 0.0));
    assert!(
        parent_axes.x < 0.0 && parent_axes.z.abs() < 1e-9 * expected,
        "an unoriented node still reads the parent's axes, which is right: {parent_axes:?}"
    );
}

/// Composition runs the whole chain, not one level of it.
///
/// Two quarter turns about different axes do not commute, so a node two levels
/// down a chain of oriented frames is the case that tells a real composition
/// apart from consulting the node's own orientation and stopping there.
#[test]
fn orientation_composes_all_the_way_up_the_chain() {
    use phys::math::Quat;
    let (mut w, node) = an_earth();
    let planet = w.tree.nodes[node.get()].parent;
    let a = Quat::from_axis_angle(v3(0.0, 1.0, 0.0), std::f64::consts::FRAC_PI_2);
    let b = Quat::from_axis_angle(v3(1.0, 0.0, 0.0), std::f64::consts::FRAC_PI_2);
    w.tree.nodes[planet.get()].motion.orientation = b;
    w.tree.nodes[node.get()].motion.orientation = a;

    let q = w.tree.axes_from(w.tree.root, node);
    let composed = a.then(b);
    assert!(
        (q.conjugate().then(composed).angle()).abs() < 1e-12,
        "axes_from should be the child's turn then the parent's: {q:?} against {composed:?}"
    );

    // And the same vector carried down lands where the composition says.
    let v = v3(1.0, 0.0, 0.0);
    let there = w.tree.into_axes_of(w.tree.root, node, v);
    let by_hand = composed.conjugate().rotate(v);
    assert!(
        (there - by_hand).norm() < 1e-12,
        "{there:?} against {by_hand:?}"
    );
}
