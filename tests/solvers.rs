//! Each tier's solver must be right on its own terms before the ladder means
//! anything.

use phys::math::{v3, Vec3};
use phys::solvers::*;
use phys::state::*;
use phys::units::*;

fn two_body(sep: f64, m1: f64, m2: f64) -> Vec<Body> {
    let v = (G * (m1 + m2) / sep).sqrt();
    vec![
        Body { pos: v3(-sep * m2 / (m1 + m2), 0.0, 0.0), vel: v3(0.0, -v * m2 / (m1 + m2), 0.0), mass: m1, ..Default::default() },
        Body { pos: v3(sep * m1 / (m1 + m2), 0.0, 0.0), vel: v3(0.0, v * m1 / (m1 + m2), 0.0), mass: m2, ..Default::default() },
    ]
}

/// Retarded gravity must not unbind a bound orbit.
///
/// This is Laplace's objection to finite-speed gravity, and a naive
/// implementation really does fail it: aberrating the source position alone
/// produces a tangential force that pumps the orbit. Fifty orbits is enough for
/// that failure to be unmistakable.
#[test]
fn retarded_gravity_keeps_orbits_bound() {
    let sep = AU;
    let period = std::f64::consts::TAU * (sep.powi(3) / (G * (M_SUN + M_EARTH))).sqrt();
    let params = gravity::GravityParams { theta: 0.0, softening: 0.0, retarded: true, post_newtonian: false };
    let mut b = two_body(sep, M_SUN, M_EARTH);
    let e0 = gravity::total_energy(&b, params);
    let l0 = total_spin(&b, Vec3::ZERO);
    let mut worst_radius = 0.0f64;
    let steps = 2000;
    for _ in 0..(steps * 50) {
        gravity::step_leapfrog(&mut b, period / steps as f64, params);
        let r = (b[1].pos - b[0].pos).norm();
        worst_radius = worst_radius.max((r - sep).abs() / sep);
    }
    let e1 = gravity::total_energy(&b, params);
    let l1 = total_spin(&b, Vec3::ZERO);
    let de = (e1 - e0).abs() / e0.abs();
    let dl = (l1 - l0).norm() / l0.norm();
    println!("50 orbits: dE/E={de:.3e} dL/L={dl:.3e} worst radius drift={worst_radius:.3e}");
    assert!(de < 1e-6, "energy drift {de:.3e}");
    assert!(dl < 1e-10, "angular momentum drift {dl:.3e}");
    assert!(worst_radius < 1e-3, "orbit radius drifted {worst_radius:.3e}");
}

/// Leapfrog is second order: halving the step must quarter the error.
#[test]
fn integrator_is_second_order() {
    let sep = AU;
    let period = std::f64::consts::TAU * (sep.powi(3) / (G * (M_SUN + M_EARTH))).sqrt();
    let params = gravity::GravityParams { theta: 0.0, softening: 0.0, retarded: false, post_newtonian: false };
    let mut errs = Vec::new();
    for div in [100usize, 200, 400] {
        let mut b = two_body(sep, M_SUN, M_EARTH);
        let dt = period / div as f64;
        for _ in 0..div {
            gravity::step_leapfrog(&mut b, dt, params);
        }
        // After exactly one period the bodies should be back where they started.
        let expected = two_body(sep, M_SUN, M_EARTH);
        errs.push((b[1].pos - expected[1].pos).norm() / sep);
    }
    println!("period-return error: {errs:?}");
    for w in errs.windows(2) {
        let ratio = w[0] / w[1];
        assert!(ratio > 3.0, "convergence ratio {ratio:.2}, expected ~4 (2nd order)");
    }
}

/// SPH forces are pairwise and antisymmetric, so momentum is conserved to
/// machine precision — not to truncation.
#[test]
fn sph_conserves_momentum_exactly() {
    let agg = Aggregate::neutral(1e30, 1e12, 1e4, Composition::solar());
    let spec = phys::prolong::ProlongSpec::new(600, phys::prolong::Profile::Uniform, phys::prolong::MassSpectrum::Equal, BodyKind::GasParcel);
    let (mut b, _) = phys::prolong::prolong(&agg, spec, 3, 0x77, 0);
    let params = hydro::HydroParams { h: 1e11, cooling: false, ..Default::default() };
    let p0 = total_momentum(&b);
    for _ in 0..20 {
        hydro::step(&mut b, 1.0, params);
    }
    let p1 = total_momentum(&b);
    let scale = Scales::of(&b).momentum;
    let drift = (p1 - p0).norm() / scale;
    println!("SPH momentum drift {drift:.3e} of total momentum content");
    assert!(drift < 1e-12, "SPH momentum drift {drift:.3e}");
}

/// Artificial viscosity must heat on compression and never cool.
#[test]
fn sph_viscosity_only_heats() {
    let agg = Aggregate::neutral(1e30, 1e12, 1e4, Composition::solar());
    let spec = phys::prolong::ProlongSpec::new(400, phys::prolong::Profile::Uniform, phys::prolong::MassSpectrum::Equal, BodyKind::GasParcel);
    let (mut b, _) = phys::prolong::prolong(&agg, spec, 3, 0x78, 0);
    // Drive a convergent flow.
    for body in b.iter_mut() {
        body.vel = body.pos.unit().scale(-3e4);
    }
    let params = hydro::HydroParams { h: 2e11, cooling: false, alpha: 1.0, beta: 2.0, ..Default::default() };
    let t0: f64 = b.iter().map(|x| x.temperature).sum::<f64>() / b.len() as f64;
    for _ in 0..10 {
        hydro::step(&mut b, 1e5, params);
    }
    let t1: f64 = b.iter().map(|x| x.temperature).sum::<f64>() / b.len() as f64;
    println!("compressive flow: T {t0:.1} -> {t1:.1} K");
    assert!(t1 >= t0, "compression must not cool: {t0} -> {t1}");
}

/// Radiative cooling must be positive, and stronger for metal-rich gas.
#[test]
fn cooling_curve_behaves() {
    for t in [100.0, 1e4, 1e5, 1e7, 1e8] {
        let low = hydro::cooling_rate(t, 1e-20, 0.0001);
        let high = hydro::cooling_rate(t, 1e-20, 0.02);
        assert!(low >= 0.0 && high >= 0.0, "cooling must not be negative at {t} K");
        if t < 1e7 {
            assert!(high > low, "metals must cool faster at {t} K");
        }
    }
    // The classic peak: line cooling around 10^5 K beats free-free at 10^7 K.
    let peak = hydro::cooling_rate(1e5, 1e-20, 0.02);
    let ff = hydro::cooling_rate(1e7, 1e-20, 0.02);
    assert!(peak > ff, "line-cooling peak should exceed bremsstrahlung");
}

/// The Lennard-Jones minimum is at 2^(1/6) sigma. If this moves, every
/// molecular structure the engine produces is wrong.
#[test]
fn lennard_jones_minimum() {
    let (sigma, epsilon) = md::lj_params(Species::Carbon);
    let r_min = 2f64.powf(1.0 / 6.0) * sigma;
    let f = md::lj_force(r_min, sigma, epsilon);
    assert!(f.abs() < 1e-3 * epsilon / sigma, "force at minimum: {f:.3e}");
    let u = md::lj_potential(r_min, sigma, epsilon);
    assert!((u + epsilon).abs() < 1e-9 * epsilon, "well depth {u:.3e} vs {epsilon:.3e}");
}

/// A molecular system with a thermostat must reach the target temperature and
/// stay there rather than exploding or freezing.
#[test]
fn molecular_dynamics_is_stable() {
    let agg = Aggregate::neutral(64.0 * 12.0 * AMU, 2e-9, 300.0, Composition::pure(Species::Carbon));
    let spec = phys::prolong::ProlongSpec::new(64, phys::prolong::Profile::Lattice, phys::prolong::MassSpectrum::Species, BodyKind::Atom);
    let (mut b, _) = phys::prolong::prolong(&agg, spec, 5, 0x1234, 0);
    let params = md::MdParams { thermostat: Some(300.0), friction: 1e12, ..Default::default() };
    let dt = md::stable_dt(&b);
    assert!(dt > 0.0 && dt < 1e-13, "implausible MD timestep {dt:.3e}");
    for tick in 0..2000u64 {
        md::step(&mut b, dt, params, 1, 0x1234, 0, tick);
    }
    for body in &b {
        assert!(body.pos.is_finite() && body.vel.is_finite(), "MD diverged");
        assert!(body.vel.norm() < C, "superluminal atom");
    }
    let t = md::temperature_of(&b);
    println!("thermostatted to {t:.1} K (target 300)");
    assert!(t > 100.0 && t < 1000.0, "thermostat produced {t} K, target 300");

    // And the noise must not repeat: two consecutive steps from the same
    // address must differ, or the "thermostat" is a constant force.
    let snapshot: Vec<Vec3> = b.iter().map(|x| x.vel).collect();
    md::step(&mut b, dt, params, 1, 0x1234, 0, 9000);
    let after: Vec<Vec3> = b.iter().map(|x| x.vel).collect();
    md::step(&mut b, dt, params, 1, 0x1234, 0, 9001);
    let after2: Vec<Vec3> = b.iter().map(|x| x.vel).collect();
    let d1: f64 = snapshot.iter().zip(&after).map(|(a, c)| (*a - *c).norm()).sum();
    let d2: f64 = after.iter().zip(&after2).map(|(a, c)| (*a - *c).norm()).sum();
    assert!(d1 > 0.0 && d2 > 0.0 && (d1 - d2).abs() / d1.max(d2) > 1e-6,
        "thermostat noise repeats between steps");
}

/// The pp chain at solar central conditions must produce roughly a solar
/// luminosity. This single number validates the whole nuclear rate path.
#[test]
fn stellar_burning_reproduces_the_sun() {
    let rho = 1.5e5; // kg/m^3, solar centre
    let t = 1.57e7; // K
    let x = 0.35;
    let eps = nuclear::pp_chain_rate(rho, t, x);
    // Integrating a realistic density/temperature profile gives ~2e-4 W/kg
    // averaged over the burning core; the central rate is a few times that.
    println!("central pp rate {eps:.3e} W/kg");
    assert!(eps > 1e-5 && eps < 1e-1, "pp rate {eps:.3e} W/kg is not solar");

    // Temperature sensitivity: the pp chain goes roughly as T^4.
    let e1 = nuclear::pp_chain_rate(rho, t, x);
    let e2 = nuclear::pp_chain_rate(rho, t * 1.1, x);
    let exponent = (e2 / e1).ln() / 1.1f64.ln();
    println!("d ln eps / d ln T = {exponent:.2}");
    assert!(exponent > 3.0 && exponent < 6.5, "pp temperature exponent {exponent:.2}");

    // CNO must be steeper still, and must take over at high temperature.
    let c1 = nuclear::cno_rate(rho, 2.5e7, x, 0.01);
    let c2 = nuclear::cno_rate(rho, 2.5e7 * 1.1, x, 0.01);
    let cno_exp = (c2 / c1).ln() / 1.1f64.ln();
    println!("CNO exponent {cno_exp:.2}");
    assert!(cno_exp > exponent, "CNO must be more temperature sensitive than pp");
}

/// Fusing hydrogen to helium releases 0.7% of the rest mass. Every stellar
/// energy budget in the engine depends on this number.
#[test]
fn fusion_energy_matches_the_binding_curve() {
    let mass = 1.0;
    let e = nuclear::fusion_energy(Species::Hydrogen, Species::Helium, mass);
    let fraction = e / (mass * C2);
    println!("H->He releases {:.4}% of rest mass", fraction * 100.0);
    assert!((fraction - 0.00712).abs() < 5e-4, "got {fraction:.5}");

    // Iron is the floor: nothing exothermic goes past it.
    for s in [Species::Carbon, Species::Oxygen, Species::Silicon] {
        assert!(
            nuclear::fusion_energy(s, Species::Iron, 1.0) > 0.0,
            "fusing {} to iron should release energy",
            s.name()
        );
    }
    assert!(
        nuclear::fusion_energy(Species::Iron, Species::Other, 1.0) < 0.0,
        "fusing past iron must cost energy"
    );
}

/// Quantum bookkeeping: uncertainty is enforced, spectra are right, and the
/// information bound is finite.
#[test]
fn quantum_limits_hold() {
    use phys::solvers::quantum::*;
    // Uncertainty principle.
    for dx in [1e-15, 1e-10, 1e-6] {
        let dp = enforce_uncertainty(dx, 0.0);
        assert!(dx * dp >= H_BAR / 2.0 * (1.0 - 1e-12), "uncertainty violated at dx={dx:.1e}");
    }
    // Lyman-alpha: n=2 -> n=1 in hydrogen is 121.6 nm.
    let h = Ensemble::hydrogenic(1.0, 5, 1e4);
    let de = h.transition_energy(1, 0);
    let lambda = H_PLANCK * C / de;
    println!("Lyman-alpha {:.2} nm", lambda * 1e9);
    assert!((lambda - 121.6e-9).abs() < 1e-9, "got {:.3} nm", lambda * 1e9);

    // Wien's law from the sampler: mean photon energy ~ 2.7 kT.
    let mut s = phys::rng::Stream::at(1, 2, 0, phys::rng::Purpose::PhotonEmission);
    let t = 5772.0;
    let mean: f64 = (0..20000).map(|_| sample_blackbody_photon(t, &mut s)).sum::<f64>() / 20000.0;
    let ratio = mean / (K_B * t);
    println!("mean photon energy = {ratio:.3} kT (Planck: 2.701)");
    assert!((ratio - 2.701).abs() < 0.15, "got {ratio:.3} kT");

    // The information bound is finite, which is what makes "down to subatomic"
    // a bounded request rather than an unbounded one.
    let states = detail_ceiling(1e-27, 300.0, M_PROTON);
    assert!(states.is_finite() && states > 0.0);
    println!("distinguishable states in a cubic nanometre of 300 K gas: {states:.3e}");
}

/// A decay sampled for a single nucleus is stochastic; averaged over many it
/// reproduces the half-life.
#[test]
fn decay_statistics_match_the_half_life() {
    use phys::solvers::nuclear::Isotope;
    let iso = Isotope::Neutron;
    let mut s = phys::rng::Stream::at(9, 9, 0, phys::rng::Purpose::Decay);
    let n = 40000;
    let mean: f64 = (0..n).map(|_| iso.sample_lifetime(&mut s)).sum::<f64>() / n as f64;
    let expected = 1.0 / iso.decay_constant();
    let err = (mean - expected).abs() / expected;
    println!("mean neutron lifetime {mean:.1} s (expected {expected:.1})");
    assert!(err < 0.02, "mean lifetime off by {:.1}%", err * 100.0);
}
