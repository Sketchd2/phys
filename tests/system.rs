//! What does not turn with its planet is not its child — `docs/PLAY.md` §7
//! Phase 5, the owner's rule.
//!
//! Since a child's place is stated in its parent's own turning axes, a thing
//! inside a spinning planet's node that is not going round with it — a parcel
//! of gas far out at rest in space, a satellite, a rocket that has left the
//! air — is carried correctly along its straight line, but in axes that turn
//! under it and as the child of a frame it does not share. The owner's rule:
//! it belongs to a node that does not turn, shared with the planet, and the
//! split puts it there. What decides it is motion — its velocity against the
//! frame over half the speed the frame's turn gives where it is — and the
//! system node is made on the first such split and folded back when it holds
//! the planet alone.

use phys::engine::{default_spec, World};
use phys::ids::NodeIdx;
use phys::math::{v3, Vec3};
use phys::state::{Composition, Matter};
use phys::tree::Tree;
use phys::units::*;

const EARTH_RADIUS: f64 = 6.371e6;
const EARTH_MASS: f64 = 5.972e24;
const AU: f64 = 1.496e11;
const DAY: f64 = 86_164.0;
/// Heavy enough for the world's books to see what happens to it: a planet
/// going round a star carries 1.8e29 kg m/s, and a million-kilogram parcel's
/// kilometre a second is below what a double can add to that.
const PARCEL: f64 = 1.0e20;

/// A star, and a planet at 1 AU turning once a sidereal day about its own z,
/// its contents going round with it, with one parcel of gas promoted inside it at `at` in the planet's axes,
/// moving at `velocity` against the planet's turning frame.
fn a_turning_planet(at: Vec3, velocity: Vec3) -> (World, NodeIdx, NodeIdx, NodeIdx) {
    let mass = 1.989e30;
    let radius = 6.957e8;
    let mut m = Matter::neutral(mass, radius, 5772.0, Composition::primordial());
    m.internal_energy = 0.3 * G * mass * mass / radius;
    m.gravitational_binding = -0.6 * G * mass * mass / radius;
    let mut spec = default_spec(Tier::Stellar);
    spec.count = 8;
    let mut w = World::new(Tree::new(0x5157, m, Tier::Stellar, spec), 20.0);
    let star = w.tree.root;
    w.tree.refine(star);
    let planet = w.tree.promote(star, 0, default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.matter.spin = v3(0.0, 0.0, std::f64::consts::TAU / DAY * n.matter.moment_of_inertia());
        n.spec.count = 8;
        n.motion.offset = v3(AU, 0.0, 0.0);
        n.motion.velocity = v3(0.0, 29_780.0, 0.0);
        n.sync_spin_rate();
    }
    w.tree.refine(planet);
    // Its contents going round with it, as the rock of a planet does: a turn
    // counts as a frame where it outweighs the random motion of what turns
    // (`World::turns_as_frame`), and drawn as hot gas this planet's lumps go
    // at kilometres a second against a turn of 465 m/s at its equator.
    for b in w.tree.nodes[planet.get()].bodies.iter_mut() {
        b.vel = Vec3::ZERO;
    }
    let parcel = w.tree.promote(planet, 0, default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[parcel.get()];
        n.matter = Matter::neutral(PARCEL, 2.0e4, 290.0, Composition::primordial());
        n.motion.offset = at;
        n.motion.velocity = velocity;
        n.motion.spin_rate = Vec3::ZERO;
    }
    (w, star, planet, parcel)
}

fn turn() -> Vec3 {
    v3(0.0, 0.0, std::f64::consts::TAU / DAY)
}

/// **A parcel at rest in space beside a turning planet leaves it**, for a
/// node that does not turn, made above the planet and holding both — and in
/// it the parcel is where it is in space, pulled by the planet it is beside.
#[test]
fn a_parcel_at_rest_in_space_leaves_its_turning_planet() {
    // Three planet radii out, and at rest in space: against the frame it is
    // going backwards at the whole of the frame's speed there, 1.39 km/s.
    let at = v3(3.0 * EARTH_RADIUS, 0.0, 0.0);
    let (mut w, star, planet, parcel) = a_turning_planet(at, Vec3::ZERO - turn().cross(at));
    let before = w.conserved();
    let key = w.tree.nodes[planet.get()].key;
    // From the planet's centre, root-aligned, whichever node each is in.
    let seen = |w: &World| {
        w.tree.offset_from(star, parcel, Vec3::ZERO).value - w.tree.offset_from(star, planet, Vec3::ZERO).value
    };
    let start = seen(&w);
    w.step_frame(2000.0);
    let system = w.tree.nodes[planet.get()].parent;
    // The books across the frame it parted in, where nothing else was solved.
    let parted = w.conserved();
    // Against the parcel's own share, which the books resolve to about 1e-10:
    // energy they do not, since it carries the star's rest mass.
    let p = PARCEL * turn().cross(at).norm();
    let moved_by = (parted.momentum - before.momentum).norm() / p;
    let turned_by = (parted.angular_momentum - before.angular_momentum).norm() / (p * at.norm());
    println!(
        "  across the frame it parted in: momentum moved {moved_by:.1e} and angular momentum {turned_by:.1e} of the parcel's against the frame"
    );
    assert!(moved_by < 1e-6, "parting moved the world's momentum by {moved_by:e} of the parcel's");
    assert!(turned_by < 1e-6, "parting moved the world's angular momentum by {turned_by:e} of the parcel's");
    println!(
        "  after one frame: {} parted, {} system nodes made; the planet's parent is {} (the star is {}), the parcel's {}",
        w.stats.parted,
        w.stats.systems_made,
        system.get(),
        star.get(),
        w.tree.nodes[parcel.get()].parent.get()
    );
    assert_ne!(system, star, "the planet is still the star's child: no system node was made");
    assert_eq!(w.tree.nodes[system.get()].parent, star, "the system node is not the star's child");
    assert_eq!(w.tree.nodes[parcel.get()].parent, system, "the parcel is still the turning planet's child");
    assert_eq!(w.tree.angular_velocity(system), Vec3::ZERO, "the system node turns");
    assert_eq!(w.tree.nodes[system.get()].key, key, "the system node did not take the planet's address");
    assert_eq!(w.stats.parted, 1);
    assert_eq!(w.stats.systems_made, 1);

    // Where it is in space, as the planet turns under it: a quarter of an
    // hour, 3.8 deg of the planet's turn.
    w.pace_fixed(30.0);
    for _ in 0..30 {
        w.step_frame(2000.0);
    }
    let end = seen(&w);
    let turned = start.cross(end).z.atan2(start.dot(end));
    // Root-aligned, against the planet's pull and the star's.
    let pull = w.tree.gravity_at_point(parcel, Vec3::ZERO);
    let from_star = w.tree.offset_from(star, parcel, Vec3::ZERO).value;
    let expected = end.scale(-G * EARTH_MASS / (end.norm2() * end.norm()))
        + from_star.scale(-G * 1.989e30 / (from_star.norm2() * from_star.norm()));
    let g = pull.norm();
    let closed = G * EARTH_MASS / end.norm2();
    println!(
        "  {:.0} s on: it moved {:.3e} rad round the planet against the planet's {:.3e}; fell {:.1} m; \
         gravity {g:.5} m/s2 against {:.5} from the planet and the star; falling freely is {:.0} m",
        w.time,
        turned,
        turn().z * w.time,
        start.norm() - end.norm(),
        expected.norm(),
        0.5 * closed * w.time * w.time
    );
    assert!(w.time > 800.0, "the world ran {} s, not a quarter of an hour", w.time);
    assert!(turned.abs() < 1e-3 * turn().z * w.time, "the parcel went round with the planet it left");
    // Towards the planet, and no further than falling freely would take it:
    // the system node is solved once a resolution element of its own, and a
    // node the size of a planet holding two things resolves 5,000 km, so it
    // is first solved at a minute and then coasts at what that gave it.
    let free = 0.5 * closed * w.time * w.time;
    assert!(start.norm() - end.norm() > 0.0, "the parcel did not fall towards the planet beside it");
    assert!(start.norm() - end.norm() < free, "the parcel fell further than falling freely takes it");
    assert!(
        (pull - expected).norm() < 1e-3 * expected.norm(),
        "beside the planet it feels {pull:?}, not the planet's and the star's {expected:?}"
    );
}

/// **A thing moving with the frame stays**, and nothing is made: a parcel
/// going round with the planet, and one a little off it — a ball thrown
/// across the ground, at a fifth of the frame's speed.
#[test]
fn a_parcel_moving_with_the_frame_stays() {
    let at = v3(0.0, 2.0 * EARTH_RADIUS, 0.0);
    let off = turn().cross(at).scale(0.2);
    let (mut w, star, planet, parcel) = a_turning_planet(at, off);
    for _ in 0..30 {
        w.step_frame(2000.0);
    }
    println!(
        "  moving at {:.1} m/s against a frame going {:.1} m/s there: {} parted, {} made",
        off.norm(),
        turn().cross(at).norm(),
        w.stats.parted,
        w.stats.systems_made
    );
    assert_eq!(w.tree.nodes[parcel.get()].parent, planet, "a parcel moving with its planet left it");
    assert_eq!(w.tree.nodes[planet.get()].parent, star);
    assert_eq!(w.stats.systems_made, 0);
}

/// **And back.** A parcel beside the planet in its system node that comes
/// inside what the planet holds, moving with its turn, is the planet's again,
/// and the system node — holding nothing but the planet — folds back into
/// the star, which gets the planet back at its own address.
///
/// The books are read across the two moves the pass makes, on a twin of the
/// world built and stepped the same way: the frame it happens in also solves
/// the parcel's own contents — a lump of gas drawn and solved every frame —
/// and that solve moves its own angular momentum by 1.9e-3 of what the parcel
/// carries against the frame, which is the parcel's business and not the
/// move's.
#[test]
fn a_parcel_that_turns_with_its_planet_again_rejoins_it_and_the_system_folds() {
    let at = v3(3.0 * EARTH_RADIUS, 0.0, 0.0);
    // Parted, and the system node's first solve spent — a new arrival is due
    // at once — and then the parcel put inside the planet's own contents,
    // going round with it.
    let world = || {
        let (mut w, star, planet, parcel) = a_turning_planet(at, Vec3::ZERO - turn().cross(at));
        let key = w.tree.nodes[planet.get()].key;
        for _ in 0..5 {
            w.step_frame(2000.0);
        }
        let system = w.tree.nodes[planet.get()].parent;
        assert_ne!(system, star, "nothing parted, so there is nothing to come back from");
        let inside = {
            let reach = w.tree.contents_reach(planet, parcel);
            let p = &w.tree.nodes[planet.get()].motion;
            let there = v3(0.0, 0.3 * reach, 0.0);
            let into = w.tree.facing(system).conjugate().then(w.tree.facing(planet)).unit();
            (p.offset + into.rotate(there), p.velocity + into.rotate(w.tree.own_turn(planet).cross(there)))
        };
        let n = &mut w.tree.nodes[parcel.get()];
        n.motion.offset = inside.0;
        n.motion.velocity = inside.1;
        (w, star, planet, parcel, system, key)
    };
    let (mut w, star, planet, parcel, system, key) = world();
    w.step_frame(2000.0);
    println!(
        "  {} rejoined, {} folded; the planet's parent is {} (the star {}), at its own address again: {}",
        w.stats.rejoined,
        w.stats.systems_folded,
        w.tree.nodes[planet.get()].parent.get(),
        star.get(),
        w.tree.nodes[planet.get()].key == key
    );
    assert_eq!(w.tree.nodes[parcel.get()].parent, planet, "the parcel did not rejoin its planet");
    assert_eq!(w.tree.nodes[planet.get()].parent, star, "the system node did not fold back");
    assert!(!w.tree.nodes[system.get()].alive, "the system node is still alive");
    assert_eq!(w.tree.nodes[planet.get()].key, key, "the planet did not get its address back");

    let (mut twin, _, planet, parcel, system, _) = world();
    let before = twin.conserved();
    assert!(twin.reparent(parcel, planet), "the twin's parcel could not rejoin");
    assert!(twin.tree.fold_into_parent(system).is_some(), "the twin's system node could not fold");
    let after = twin.conserved();
    let p = PARCEL * turn().cross(at).norm();
    let moved_by = (after.momentum - before.momentum).norm() / p;
    let turned_by = (after.angular_momentum - before.angular_momentum).norm() / (p * at.norm());
    println!("  the world's books across the two moves: momentum moved {moved_by:.1e} and angular momentum {turned_by:.1e} of the parcel's");
    assert!(moved_by < 1e-6, "rejoining moved the world's momentum by {moved_by:e} of the parcel's");
    assert!(turned_by < 1e-6, "rejoining moved the world's angular momentum by {turned_by:e} of the parcel's");
}
