//! Covalent bonds in the particle solver, against spectroscopy.
//!
//! A bond has three measured properties — length, dissociation energy and
//! vibrational frequency — and they are not independent: fixing any two fixes
//! the third. That makes them a real test rather than a self-consistency check,
//! because the solver has no freedom to get one right by getting another wrong.

use phys::math::{v3, Vec3};
use phys::solvers::md::{self, Bonded};
use phys::state::{Body, Composition};
use phys::units::*;

fn atom(species: Species, pos: Vec3, mass_amu: f64) -> Body {
    Body {
        pos,
        mass: mass_amu * AMU,
        composition: Composition::pure(species),
        ..Default::default()
    }
}

fn hydrogen_molecule(separation: f64) -> (Vec<Body>, Bonded) {
    let bodies = vec![
        atom(Species::Hydrogen, v3(0.0, 0.0, 0.0), 1.00794),
        atom(Species::Hydrogen, v3(separation, 0.0, 0.0), 1.00794),
    ];
    let mut bonded = Bonded::default();
    bonded.bond(&bodies, 0, 1);
    (bodies, bonded)
}

/// A diatomic must vibrate at `sqrt(k / mu)`.
///
/// This is the number a spectrometer measures. H2's fundamental is 4401 /cm;
/// the harmonic value from `k = 575 N/m` and the reduced mass of two protons is
/// close to it, and the solver has to land on the same figure by integrating
/// rather than by being told.
#[test]
fn a_diatomic_vibrates_at_its_spectroscopic_frequency() {
    let (r0, well, k) = md::covalent(Species::Hydrogen, Species::Hydrogen);
    // Displace by a hundredth of the bond length: small enough that the Morse
    // well is harmonic to well under a percent.
    let stretch = r0 * 0.01;
    let (mut bodies, bonded) = hydrogen_molecule(r0 + stretch);

    let reduced = bodies[0].mass * bodies[1].mass / (bodies[0].mass + bodies[1].mass);
    let expect = std::f64::consts::TAU * (reduced / k).sqrt();

    let dt = md::stable_dt_bonded(&bodies, &bonded);
    assert!(dt < expect / 20.0, "timestep {dt:.3e} s against a period of {expect:.3e} s");
    let params = md::MdParams::default();

    let mut crossings = Vec::new();
    let mut previous = stretch;
    for step in 1..4000u64 {
        md::step_bonded(&mut bodies, &bonded, dt, params, 1, 0x1, 0, step);
        let x = (bodies[1].pos - bodies[0].pos).norm() - r0;
        if (previous > 0.0) != (x > 0.0) {
            let frac = previous / (previous - x);
            crossings.push((step as f64 - 1.0 + frac) * dt);
        }
        previous = x;
        if crossings.len() >= 5 {
            break;
        }
    }
    assert!(crossings.len() >= 5, "only {} crossings", crossings.len());
    let measured = crossings[4] - crossings[2];
    let wavenumber = 1.0 / (measured * C * 100.0);
    println!(
        "  period {measured:.4e} s (theory {expect:.4e}), {wavenumber:.0} /cm \
         against a measured 4401 /cm for H2, well {:.2} eV",
        well / EV
    );
    assert!(
        (measured - expect).abs() / expect < 0.02,
        "period {measured:.4e} s against {expect:.4e} s"
    );
    // And the same number in the units a spectroscopist would recognise.
    assert!(
        (wavenumber - 4401.0).abs() / 4401.0 < 0.10,
        "{wavenumber:.0} /cm against a measured 4401 /cm"
    );
}

/// Dissociation must be emergent: more than the well depth and the molecule
/// comes apart, less and it does not.
///
/// Nothing in the solver tests an extension against a threshold. The bond stops
/// pulling because the Morse potential flattens out, which is what a real bond
/// does and what a harmonic one never does however large you make it.
#[test]
fn a_molecule_comes_apart_at_its_dissociation_energy() {
    let (r0, well, _) = md::covalent(Species::Hydrogen, Species::Hydrogen);
    let params = md::MdParams::default();

    for (label, fraction, should_break) in
        [("below", 0.75, false), ("above", 1.30, true)]
    {
        let (mut bodies, bonded) = hydrogen_molecule(r0);
        // Kinetic energy along the bond, in the centre-of-mass frame, so all of
        // it is available to the bond and none is bulk motion.
        let reduced = bodies[0].mass * bodies[1].mass / (bodies[0].mass + bodies[1].mass);
        let speed = (2.0 * fraction * well / reduced).sqrt() / 2.0;
        bodies[0].vel = v3(-speed, 0.0, 0.0);
        bodies[1].vel = v3(speed, 0.0, 0.0);

        let dt = md::stable_dt_bonded(&bodies, &bonded);
        let start = md::step_bonded(&mut bodies, &bonded, 0.0, params, 1, 0x1, 0, 0).before;
        let mut furthest: f64 = 0.0;
        for step in 1..12000u64 {
            md::step_bonded(&mut bodies, &bonded, dt, params, 1, 0x1, 0, step);
            furthest = furthest.max((bodies[1].pos - bodies[0].pos).norm());
        }
        let separation = (bodies[1].pos - bodies[0].pos).norm();
        let broke = !bonded.dissociating(&bodies).is_empty() && separation > 4.0 * r0;
        println!(
            "  {label:>5} the well ({:.0}% of {:.2} eV): reached {:.1} bond lengths, \
             energy drift {:.2e}",
            fraction * 100.0,
            well / EV,
            furthest / r0,
            drift(start, &bodies, &bonded)
        );
        assert_eq!(
            broke, should_break,
            "{label} the well: separation {:.2} bond lengths",
            separation / r0
        );
        if !should_break {
            // A bound molecule must turn around exactly where the Morse
            // potential says it will: `D_e (1 - e^{-a x})^2 = E` gives
            // `x = -ln(1 - sqrt(E/D_e)) / a`. Asserting against a round number
            // instead would only be testing that it stayed somewhere nearby.
            let alpha = bonded.bonds[0].alpha;
            let turning = r0 - (1.0 - fraction.sqrt()).ln() / alpha;
            assert!(
                (furthest - turning).abs() / turning < 0.02,
                "turned around at {:.3} A, the Morse potential says {:.3} A",
                furthest * 1e10,
                turning * 1e10
            );
        }
    }
}

fn drift(before: phys::state::Conserved, bodies: &[Body], bonded: &Bonded) -> f64 {
    let skip = bonded.exclusions(bodies.len());
    let u = bonded.energy(bodies)
        + md::potential_energy_excluding(bodies, md::MdParams::default(), &skip);
    let after = phys::solvers::measure(bodies, u);
    let scale = before.energy.abs().max(after.energy.abs()).max(1e-30);
    (after.energy - before.energy).abs() / scale
}

/// Bonded forces are internal, so they must not move the centre of mass or
/// change the total angular momentum.
///
/// The angle term is where this goes wrong in practice: it is easy to write the
/// two end forces correctly and forget that the centre atom takes their
/// reaction, and the result is a molecule that flies across the box under its
/// own bending.
#[test]
fn bonded_forces_are_internal() {
    let (r0, _, _) = md::covalent(Species::Hydrogen, Species::Oxygen);
    let theta: f64 = 104.5f64.to_radians();
    let mut bodies = vec![
        atom(Species::Oxygen, v3(0.0, 0.0, 0.0), 15.999),
        atom(Species::Hydrogen, v3(r0, 0.0, 0.0), 1.00794),
        atom(
            Species::Hydrogen,
            v3(r0 * theta.cos(), r0 * theta.sin(), 0.0),
            1.00794,
        ),
    ];
    let mut bonded = Bonded::default();
    bonded.bond(&bodies, 0, 1);
    bonded.bond(&bodies, 0, 2);
    bonded.bend(&bodies, 1, 0, 2, 70.0 * EV);

    // Distort it: both terms have to be away from their rest state or the test
    // passes on forces that are all zero.
    bodies[1].pos = v3(r0 * 1.3, 0.0, 0.0);
    bodies[2].pos = v3(r0 * 0.8 * 0.2f64.cos(), r0 * 0.8 * 0.2f64.sin(), 0.0);

    let f = bonded.forces(&bodies);
    let net: Vec3 = f.iter().fold(Vec3::ZERO, |a, b| a + *b);
    let scale: f64 = f.iter().map(|x| x.norm()).sum();
    let torque: Vec3 = bodies
        .iter()
        .zip(&f)
        .fold(Vec3::ZERO, |a, (b, force)| a + b.pos.cross(*force));
    let torque_scale: f64 = bodies
        .iter()
        .zip(&f)
        .map(|(b, force)| b.pos.norm() * force.norm())
        .sum();
    println!(
        "  net force {:.3e} N against {scale:.3e} N of individual forces, \
         net torque {:.3e} against {torque_scale:.3e}",
        net.norm(),
        torque.norm()
    );
    assert!(net.norm() < 1e-12 * scale, "bonded forces are not balanced");
    assert!(
        torque.norm() < 1e-12 * torque_scale,
        "bonded forces produce a net torque"
    );
}

/// A bend must actually restore the geometry it was given.
#[test]
fn an_angle_holds_a_molecule_in_shape() {
    let (r0, _, _) = md::covalent(Species::Hydrogen, Species::Oxygen);
    let rest: f64 = 104.5f64.to_radians();
    let place = |t: f64| v3(r0 * t.cos(), r0 * t.sin(), 0.0);
    let mut bodies = vec![
        atom(Species::Oxygen, v3(0.0, 0.0, 0.0), 15.999),
        atom(Species::Hydrogen, place(0.0), 1.00794),
        atom(Species::Hydrogen, place(rest), 1.00794),
    ];
    let mut bonded = Bonded::default();
    bonded.bond(&bodies, 0, 1);
    bonded.bond(&bodies, 0, 2);
    bonded.bend(&bodies, 1, 0, 2, 5.0 * EV);
    assert!(
        (bonded.angles[0].rest - rest).abs() < 1e-9,
        "the rest angle was not taken from the geometry"
    );

    // Squash it by 30 degrees and let go, with a little damping standing in for
    // the rest of the molecule's surroundings.
    bodies[2].pos = place(rest - 30f64.to_radians());
    let dt = md::stable_dt_bonded(&bodies, &bonded);
    let params = md::MdParams::default();
    let mut worst_after_settling: f64 = 0.0;
    for step in 1..30000u64 {
        md::step_bonded(&mut bodies, &bonded, dt, params, 1, 0x1, 0, step);
        for b in bodies.iter_mut() {
            b.vel = b.vel.scale(0.999);
        }
        if step > 20000 {
            let t = angle(&bodies);
            worst_after_settling = worst_after_settling.max((t - rest).abs());
        }
    }
    println!(
        "  settled to {:.2} degrees against a rest angle of {:.2}",
        angle(&bodies).to_degrees(),
        rest.to_degrees()
    );
    assert!(
        worst_after_settling.to_degrees() < 2.0,
        "settled {:.2} degrees off its rest angle",
        worst_after_settling.to_degrees()
    );
}

fn angle(b: &[Body]) -> f64 {
    let u = b[1].pos - b[0].pos;
    let v = b[2].pos - b[0].pos;
    (u.dot(v) / (u.norm() * v.norm())).clamp(-1.0, 1.0).acos()
}

/// The integrator must not pump the bonds.
///
/// A covalent bond oscillates a hundred times faster than anything else at this
/// tier, so it is the first thing to accumulate integration error. Velocity
/// Verlet is symplectic, which means the error stays bounded instead of growing
/// — and that is worth checking over enough steps for a drift to show.
#[test]
fn bonded_energy_stays_bounded() {
    let (r0, _, _) = md::covalent(Species::Hydrogen, Species::Hydrogen);
    let (mut bodies, bonded) = hydrogen_molecule(r0 * 1.15);
    let params = md::MdParams::default();
    let dt = md::stable_dt_bonded(&bodies, &bonded);
    let before = {
        let skip = bonded.exclusions(bodies.len());
        let u = bonded.energy(&bodies) + md::potential_energy_excluding(&bodies, params, &skip);
        phys::solvers::measure(&bodies, u)
    };

    let mut worst: f64 = 0.0;
    for step in 1..50000u64 {
        md::step_bonded(&mut bodies, &bonded, dt, params, 1, 0x1, 0, step);
        if step % 100 == 0 {
            worst = worst.max(drift(before, &bodies, &bonded));
        }
    }
    println!("  worst energy drift over 50000 steps: {worst:.3e}");
    assert!(worst < 1e-4, "energy drifted by {worst:.3e}");
}
