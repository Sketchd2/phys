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
            // A bound molecule must turn around exactly where its energy runs
            // out: potential at the furthest point equals the total it started
            // with. Asserting against a round number instead would only be
            // testing that it stayed somewhere nearby.
            //
            // Note this is *not* the bare Morse turning point,
            // `r0 - ln(1 - sqrt(E/D_e))/a`. The switch that lets bonds form and
            // break without an energy jump also lifts the potential towards
            // zero faster than Morse does, so the molecule turns around sooner.
            // The energy statement holds either way, which is why it is the one
            // worth asserting.
            let bond = bonded.bonds[0];
            let total = bond.morse(r0) + fraction * well;
            let at_turning = bond.energy(furthest);
            let bare = r0 - (1.0 - fraction.sqrt()).ln() / bond.alpha;
            println!(
                "        turned at {:.3} A where V = {:.4} D_e, total energy \
                 {:.4} D_e (bare Morse would turn at {:.3} A)",
                furthest * 1e10,
                at_turning / well,
                total / well,
                bare * 1e10
            );
            assert!(
                (at_turning - total).abs() < 0.01 * well,
                "turned around where V = {:.4} D_e, energy says {:.4} D_e",
                at_turning / well,
                total / well
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

// ---------------------------------------------------------------------------
// Reactions
// ---------------------------------------------------------------------------

/// Two atoms that meet must bond, and the meeting must cost nothing.
///
/// This is the whole difficulty of reactive chemistry in one assertion. Adding
/// a term to the potential ordinarily *changes* the potential, and a simulation
/// that gains energy every time two atoms find each other is worthless however
/// plausible its chemistry. Because a bond is only ever created where its
/// switch is zero, the moment of creation is energetically invisible.
#[test]
fn a_bond_forms_without_changing_the_energy() {
    let (_, well, _) = md::covalent(Species::Hydrogen, Species::Hydrogen);
    // Outside the capture radius, closing slowly. The capture radius is set by
    // where dispersion dies rather than by the bond length, so it is several
    // times the bond length — for hydrogen, 4.9 angstroms against a 0.74
    // angstrom bond.
    let capture = Bonded::capture_radius(Species::Hydrogen, Species::Hydrogen);
    let mut bodies = vec![
        atom(Species::Hydrogen, v3(0.0, 0.0, 0.0), 1.00794),
        atom(Species::Hydrogen, v3(1.3 * capture, 0.0, 0.0), 1.00794),
    ];
    let approach = 600.0;
    bodies[0].vel = v3(approach, 0.0, 0.0);
    bodies[1].vel = v3(-approach, 0.0, 0.0);

    let mut bonded = Bonded::default();
    assert_eq!(
        bonded.react(&bodies).formed,
        0,
        "bonded at {:.2} A, outside a capture radius of {:.2} A",
        1.3 * capture * 1e10,
        capture * 1e10
    );

    let params = md::MdParams::default();
    let dt = md::stable_dt_reactive(&bodies);
    let before = total_energy(&bodies, &bonded, params);
    let mut formed_at = None;
    let mut worst_jump: f64 = 0.0;
    for step in 1..40000u64 {
        let (_, reaction) =
            md::step_reactive(&mut bodies, &mut bonded, dt, params, 1, 0x1, 0, step);
        worst_jump = worst_jump.max(reaction.energy_change.abs());
        if reaction.formed > 0 && formed_at.is_none() {
            formed_at = Some(step);
        }
    }
    let step = formed_at.expect("the two atoms never bonded");
    let after = total_energy(&bodies, &bonded, params);
    let drift = (after - before).abs() / well;
    println!(
        "  bonded at step {step} ({:.2} fs); largest energy step the chemistry took: \
         {worst_jump:.3e} J; total drift over the run {drift:.3e} of a well depth; \
         separated again by the end",
        step as f64 * dt * 1e15
    );
    assert_eq!(worst_jump, 0.0, "forming a bond moved the potential energy");
    assert!(drift < 1e-4, "energy drifted by {drift:.3e} well depths");

    // And then it comes apart again, which is not a failure — it is the
    // physics. Two atoms approaching from infinity have positive total energy,
    // and nothing in a two-body collision can take any of it away, so they must
    // leave with what they arrived with. Real recombination needs a third body
    // to carry off the well depth, which is why `atoms_assemble_into_water`
    // has one and this test does not.
    assert!(
        bonded.bonds.is_empty(),
        "a two-body collision cannot leave a bound pair"
    );
}

/// Forming a bond releases its well depth as heat. That is what exothermic
/// means, and the conserved tuple has to show it coming from the potential.
#[test]
fn forming_a_bond_warms_the_gas() {
    let (r0, well, _) = md::covalent(Species::Hydrogen, Species::Hydrogen);
    let mut bodies = vec![
        atom(Species::Hydrogen, v3(0.0, 0.0, 0.0), 1.00794),
        atom(Species::Hydrogen, v3(3.0 * r0, 0.0, 0.0), 1.00794),
    ];
    bodies[0].vel = v3(200.0, 0.0, 0.0);
    bodies[1].vel = v3(-200.0, 0.0, 0.0);
    let mut bonded = Bonded::default();
    let params = md::MdParams::default();
    let dt = md::stable_dt_reactive(&bodies);

    let kinetic_before = phys::state::kinetic_energy_of(&bodies);
    let mut kinetic_peak: f64 = 0.0;
    for step in 1..20000u64 {
        md::step_reactive(&mut bodies, &mut bonded, dt, params, 1, 0x1, 0, step);
        kinetic_peak = kinetic_peak.max(phys::state::kinetic_energy_of(&bodies));
    }
    let released = kinetic_peak - kinetic_before;
    println!(
        "  kinetic energy {:.3} eV -> {:.3} eV at closest approach; the H–H well is {:.2} eV",
        kinetic_before / EV,
        kinetic_peak / EV,
        well / EV
    );
    assert!(bonded.bonds.len() == 1, "no bond formed");
    // The pair converts the whole well depth into motion on the way in. It does
    // not keep it — with nothing to carry the energy away this molecule
    // vibrates rather than settling, which is exactly why real synthesis needs
    // a third body.
    // Not the whole well: the pair starts a little way inside the switch, and
    // the switch has already lifted the potential slightly off the bare Morse
    // curve there. What matters is that it is the well's worth of energy and
    // not a rounding error.
    assert!(
        released > 0.8 * well,
        "released {:.3} eV of a {:.3} eV well",
        released / EV,
        well / EV
    );
}

/// Valence has to hold. A hydrogen atom that acquires a second neighbour is
/// not chemistry, it is a bug with a plausible-looking rendering.
#[test]
fn valence_limits_what_can_bond() {
    let (r0, _, _) = md::covalent(Species::Hydrogen, Species::Oxygen);
    // One oxygen with four hydrogens crowded around it. Oxygen takes two.
    let mut bodies = vec![atom(Species::Oxygen, v3(0.0, 0.0, 0.0), 15.999)];
    for k in 0..4 {
        let a = std::f64::consts::TAU * k as f64 / 4.0;
        bodies.push(atom(
            Species::Hydrogen,
            v3(1.2 * r0 * a.cos(), 1.2 * r0 * a.sin(), 0.0),
            1.00794,
        ));
    }
    let mut bonded = Bonded::default();
    let reaction = bonded.react(&bodies);
    let z = bonded.coordination(bodies.len());
    println!(
        "  {} bonds formed; oxygen holds {} (valence {}), hydrogens hold {:?}",
        reaction.formed,
        z[0],
        md::valence(Species::Oxygen),
        &z[1..]
    );
    assert_eq!(z[0], 2, "oxygen took {} bonds, valence is 2", z[0]);
    for (i, &held) in z[1..].iter().enumerate() {
        assert!(held <= 1, "hydrogen {i} took {held} bonds, valence is 1");
    }
    // And the angle between the two bonds it did take is water's.
    assert_eq!(bonded.angles.len(), 1, "one bend for a two-coordinate centre");
    assert!(
        (bonded.angles[0].rest.to_degrees() - 104.5).abs() < 0.01,
        "rest angle {:.2} degrees",
        bonded.angles[0].rest.to_degrees()
    );
}

/// Loose atoms must assemble themselves into the molecule their valences allow,
/// and it must be the same molecule every time.
#[test]
fn atoms_assemble_into_water() {
    let (r0, _, _) = md::covalent(Species::Hydrogen, Species::Oxygen);
    let build = || {
        let mut bodies = vec![
            atom(Species::Oxygen, v3(0.0, 0.0, 0.0), 15.999),
            atom(Species::Hydrogen, v3(2.1 * r0, 0.0, 0.0), 1.00794),
            atom(Species::Hydrogen, v3(-0.9 * r0, 1.9 * r0, 0.0), 1.00794),
        ];
        // A little inward drift so they meet rather than sitting still.
        bodies[1].vel = v3(-120.0, 0.0, 0.0);
        bodies[2].vel = v3(60.0, -110.0, 0.0);
        bodies
    };

    let params = md::MdParams::default();
    let run = || {
        let mut bodies = build();
        let mut bonded = Bonded::default();
        let dt = md::stable_dt_reactive(&bodies);
        for step in 1..60000u64 {
            md::step_reactive(&mut bodies, &mut bonded, dt, params, 1, 0x1, 0, step);
            // Bleed the heat of formation away, as a real solvent or a third
            // body would; otherwise the molecule keeps the 4.8 eV per bond it
            // just released and shakes itself apart.
            for b in bodies.iter_mut() {
                b.vel = b.vel.scale(0.9995);
            }
        }
        (bodies, bonded)
    };

    let (bodies, bonded) = run();
    let z = bonded.coordination(bodies.len());
    let theta = angle(&bodies).to_degrees();
    println!(
        "  assembled: oxygen holds {} bonds, H–O–H is {theta:.1} degrees, \
         O–H lengths {:.3} and {:.3} A (equilibrium {:.3})",
        z[0],
        (bodies[1].pos - bodies[0].pos).norm() * 1e10,
        (bodies[2].pos - bodies[0].pos).norm() * 1e10,
        r0 * 1e10
    );
    assert_eq!(bonded.bonds.len(), 2, "water has two bonds");
    assert_eq!(z[0], 2, "the oxygen should hold both");
    assert!(
        (theta - 104.5).abs() < 12.0,
        "H–O–H settled at {theta:.1} degrees, water is 104.5"
    );

    // Determinism: the same start must give the same molecule, bit for bit.
    let (again, _) = run();
    for (a, b) in bodies.iter().zip(&again) {
        assert_eq!(a.pos, b.pos, "the assembly did not replay");
    }
}

/// A bond that is stretched past its range must be released, and releasing it
/// must cost nothing either.
#[test]
fn a_bond_breaks_without_changing_the_energy() {
    let (r0, _, _) = md::covalent(Species::Hydrogen, Species::Hydrogen);
    let (mut bodies, _) = hydrogen_molecule(r0);
    let mut bonded = Bonded::default();
    assert_eq!(bonded.react(&bodies).formed, 1);

    // Pull the two apart at a speed nothing can hold.
    // Fast enough to clear the well: 4.75 eV over a reduced mass of half a
    // proton needs about 43 km/s of relative speed.
    bodies[0].vel = v3(-30000.0, 0.0, 0.0);
    bodies[1].vel = v3(30000.0, 0.0, 0.0);
    let params = md::MdParams::default();
    let dt = md::stable_dt_reactive(&bodies);

    let mut worst_jump: f64 = 0.0;
    let mut broke_at = None;
    for step in 1..20000u64 {
        let (_, reaction) =
            md::step_reactive(&mut bodies, &mut bonded, dt, params, 1, 0x1, 0, step);
        worst_jump = worst_jump.max(reaction.energy_change.abs());
        if reaction.broken > 0 && broke_at.is_none() {
            broke_at = Some(step);
        }
    }
    let step = broke_at.expect("the bond never let go");
    println!(
        "  released at step {step}, separation {:.2} bond lengths; largest energy step \
         the chemistry took: {worst_jump:.3e} J",
        (bodies[1].pos - bodies[0].pos).norm() / r0
    );
    assert_eq!(worst_jump, 0.0, "breaking a bond moved the potential energy");
    assert!(bonded.bonds.is_empty(), "the bond is still listed");
}

fn total_energy(bodies: &[Body], bonded: &Bonded, params: md::MdParams) -> f64 {
    let skip = bonded.exclusions(bodies.len());
    phys::state::kinetic_energy_of(bodies)
        + bonded.energy(bodies)
        + md::potential_energy_excluding(bodies, params, &skip)
}
