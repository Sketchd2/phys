//! How often a node needs re-solving, and what it costs to carry it forward.
//!
//! The engine's scheduling rests on one quantity: how long before this node's
//! state has changed by more than the resolution it is being represented at.
//! `tau = resolution / v`, where `v` is the speed of the fastest-moving part of
//! the node. Everything about which node runs when comes out of that, so it had
//! better reproduce what the objects themselves do.

use phys::coords::Frame;
use phys::math::{v3, Quat, Vec3};
use phys::state::{Aggregate, Composition};
use phys::units::*;

/// Build an aggregate with a given rotation period about z.
fn spinning(mass: f64, radius: f64, temperature: f64, period: f64) -> Aggregate {
    let mut a = Aggregate::neutral(mass, radius, temperature, Composition::solar());
    if period > 0.0 {
        let omega = std::f64::consts::TAU / period;
        a.spin = v3(0.0, 0.0, omega * a.moment_of_inertia());
    }
    a
}

/// The characteristic time has to reproduce what the objects actually do.
///
/// Not "be of the right order" — a planet takes hours to change appreciably at
/// its own size and a swimming bacterium about a second, and one expression has
/// to give both from the bulk state alone.
#[test]
fn the_cadence_matches_what_the_object_does() {
    // Earth: internal energy booked as a tenth of the binding energy, as the
    // scenario does, and a sidereal day of rotation. Sitting still, so the only
    // thing changing is its own turning.
    let m = 5.972e24;
    let r = 6.371e6;
    let binding = -0.6 * G * m * m / r;
    let mut earth = spinning(m, r, 2000.0, 86164.0);
    earth.internal_energy = -0.1 * binding;
    earth.binding_energy = binding;

    let surface = earth.angular_velocity().norm() * r;
    let tau = earth.characteristic_time(r);
    println!(
        "  Earth alone: equator {surface:.0} m/s, stirring {:.0} m/s, tau {:.1} h",
        earth.stirring_speed(),
        tau / 3600.0
    );
    assert!(
        (surface - 465.0).abs() < 5.0,
        "equatorial speed {surface:.1} m/s, should be 465"
    );
    assert_eq!(
        earth.characteristic_speed(),
        surface,
        "rotation is the fastest thing an isolated Earth does"
    );
    // One update per radian of turn: 3.8 hours, and the orientation between
    // them is exact, so it turns smoothly rather than jumping.
    let radians = earth.angular_velocity().norm() * tau;
    println!("  and that is {radians:.3} radians of turn per update");
    assert!(
        (radians - 1.0).abs() < 1e-9,
        "should be one radian per update, got {radians:.3}"
    );

    // The same planet in orbit is a different question: 30 km/s of orbital
    // motion changes where it is far faster than rotation changes which way it
    // points, and the cadence follows.
    let mut orbiting = earth;
    orbiting.momentum = v3(0.0, 29_780.0 * m, 0.0);
    let tau_orbit = orbiting.characteristic_time(r);
    println!("  Earth in orbit: tau {:.1} min", tau_orbit / 60.0);
    assert!(
        tau_orbit < tau / 50.0,
        "orbital motion should dominate: {:.1} min against {:.1} h",
        tau_orbit / 60.0,
        tau / 3600.0
    );

    // A bacterium: ten microns, swimming at twenty microns a second. Nothing
    // about the expression changes — and its 300 K of thermal jiggle, which
    // would have demanded an update every five nanoseconds, does not count,
    // because a body in equilibrium is not changing.
    let bug = {
        let mut a = Aggregate::neutral(1e-15, 1.0e-5, 300.0, Composition::solar());
        a.internal_energy = 0.0;
        a.momentum = v3(2.0e-5 * 1e-15, 0.0, 0.0);
        a
    };
    let tau_bug = bug.characteristic_time(bug.radius);
    println!(
        "  bacterium: sound speed {:.0} m/s, stirring {:.0} m/s, tau {tau_bug:.2} s",
        bug.sound_speed(),
        bug.stirring_speed()
    );
    assert!(
        tau_bug > 0.1 && tau_bug < 5.0,
        "a bacterium's characteristic time is {tau_bug:.3} s"
    );

    // Stir it, and it does need re-solving sooner. The rule is not "ignore the
    // inside", it is "ignore the inside when the inside is settled".
    let mut stirred = bug;
    stirred.internal_energy = stirred.thermal_energy() + 0.5 * stirred.mass * 1.0e-6;
    let tau_stirred = stirred.characteristic_time(stirred.radius);
    println!("  the same bacterium, stirred: tau {tau_stirred:.3} s");
    assert!(
        tau_stirred < tau_bug,
        "coherent internal motion should shorten the cadence"
    );
}

/// Rotation must be in the speed, and for a gas giant it is the *dominant*
/// term. This is the case that shows why leaving it out is not a refinement.
#[test]
fn rotation_dominates_for_a_gas_giant() {
    // Jupiter: ten hours, and a cold hydrogen envelope whose sound speed is
    // about a kilometre a second.
    let mut jupiter = spinning(1.898e27, 6.99e7, 165.0, 9.925 * 3600.0);
    jupiter.internal_energy = 0.5 * jupiter.mass * 1.0e3 * 1.0e3;

    let surface = jupiter.angular_velocity().norm() * jupiter.radius;
    let internal = jupiter.sound_speed();
    println!(
        "  Jupiter: equator {:.1} km/s, internal {:.1} km/s — rotation wins by {:.1}x",
        surface / 1e3,
        internal / 1e3,
        surface / internal
    );
    assert!(
        surface > internal * 5.0,
        "rotation should dominate: {surface:.0} against {internal:.0} m/s"
    );
    assert_eq!(
        jupiter.characteristic_speed(),
        surface,
        "the characteristic speed should be the rotation"
    );

    // Without rotation the cadence would sample it once every two revolutions.
    let blind = jupiter.radius / internal;
    let turns = jupiter.angular_velocity().norm() * blind / std::f64::consts::TAU;
    println!("  ignoring rotation would update it every {turns:.2} revolutions");
    assert!(turns > 1.0, "the aliasing case should alias: {turns:.2} turns");
    // With it, one update per radian of turn or better.
    let good = jupiter.angular_velocity().norm() * jupiter.characteristic_time(jupiter.radius);
    assert!(good <= 1.001, "still under-sampled: {good:.3} rad per update");
}

/// Carrying a node forward between solves has to be exact, or a slow cadence
/// buys nothing: what it saves in solving it loses in drift.
#[test]
fn rotation_between_updates_is_exact() {
    let period = 86164.0;
    let omega = std::f64::consts::TAU / period;
    let mut frame = Frame::at_rest(Vec3::ZERO);
    frame.spin_rate = v3(0.0, 0.0, omega);

    // A point on the equator, carried forward a whole sidereal day in coarse
    // steps. A first-order angle-axis composition would drift; this must not.
    let point = v3(6.371e6, 0.0, 0.0);
    let steps = 24;
    let dt = period / steps as f64;
    for _ in 0..steps {
        frame.advance(dt);
    }
    let back = frame.body_to_parent(point);
    let error = (back - point).norm();
    println!(
        "  one sidereal day in {steps} steps of {:.0} s: point returned to within {error:.3e} m",
        dt
    );
    assert!(
        error < 1.0e-6,
        "a full rotation left the point {error:.3e} m from where it started"
    );
    assert!(
        (frame.orientation.norm() - 1.0).abs() < 1e-12,
        "the orientation drifted off the unit sphere"
    );

    // And the step size must not matter: the same day in one step.
    let mut coarse = Frame::at_rest(Vec3::ZERO);
    coarse.spin_rate = v3(0.0, 0.0, omega);
    coarse.advance(period);
    let one = coarse.body_to_parent(point);
    println!(
        "  the same day in one step: {:.3e} m from the twenty-four-step answer",
        (one - back).norm()
    );
    assert!(
        (one - back).norm() < 1.0e-6,
        "the answer depended on how it was stepped"
    );
}

/// A quarter turn is a quarter turn, and composition is not first-order.
#[test]
fn quaternions_compose_exactly() {
    let q = Quat::from_axis_angle(v3(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
    let p = q.rotate(v3(1.0, 0.0, 0.0));
    println!("  a quarter turn about z takes x to ({:.6}, {:.6}, {:.6})", p.x, p.y, p.z);
    assert!((p - v3(0.0, 1.0, 0.0)).norm() < 1e-12, "{p:?}");

    // Four quarter turns, the long way round.
    let full = q.then(q).then(q).then(q).unit();
    let r = full.rotate(v3(1.0, 0.0, 0.0));
    assert!((r - v3(1.0, 0.0, 0.0)).norm() < 1e-12, "{r:?}");

    // A large angle in one composition, which is where angle-axis addition
    // fails: two 120-degree turns about different axes do not commute, and
    // the result is not the sum of the vectors.
    let a = Quat::from_axis_angle(v3(0.0, 0.0, 1.0), 2.0);
    let b = Quat::from_axis_angle(v3(1.0, 0.0, 0.0), 2.0);
    let ab = a.then(b);
    let ba = b.then(a);
    let spread = (ab.rotate(v3(0.0, 1.0, 0.0)) - ba.rotate(v3(0.0, 1.0, 0.0))).norm();
    println!("  two 2-radian turns, order swapped: {spread:.3} apart");
    assert!(spread > 0.5, "rotation composition should not commute");
    assert!(ab.is_finite() && (ab.norm() - 1.0).abs() < 1e-12);
}

/// Resolution is what sets the cadence, not the node's own size. The same
/// planet needs updating far more often once you can see its surface.
#[test]
fn finer_resolution_demands_a_faster_cadence() {
    let m = 5.972e24;
    let r = 6.371e6;
    let binding = -0.6 * G * m * m / r;
    let mut earth = spinning(m, r, 2000.0, 86164.0);
    earth.internal_energy = -0.1 * binding;

    let as_a_point = earth.characteristic_time(r);
    let resolved = earth.characteristic_time(r / 4000f64.cbrt());
    let watched_at_a_km = earth.characteristic_time(1000.0);
    println!(
        "  Earth as a point: {:.0} min; resolved into 4000 parcels: {:.1} min; \
         watched at 1 km: {:.1} s",
        as_a_point / 60.0,
        resolved / 60.0,
        watched_at_a_km
    );
    assert!(resolved < as_a_point / 10.0, "resolving it changed nothing");
    assert!(
        watched_at_a_km < 10.0,
        "watching a kilometre of it should demand seconds, not {watched_at_a_km:.1} s"
    );
}
