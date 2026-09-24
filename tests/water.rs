//! Phase 5 — Water. `docs/PLAY.md` §7 and §4.
//!
//! Each test here names the piece of the phase it holds, and the numbers it
//! prints are the ones the phase's commits quote.

use phys::chem::{Mixture, Phase};
use phys::engine::World;
use phys::eos::{Condensed, Eos};
use phys::material::substances;
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::Tree;
use phys::units::Tier;

/// A ball of one substance in one phase, at `Continuum`, materialised.
fn ball_of(
    seed: u64,
    arrangement: phys::chem::Arrangement,
    phase: Phase,
    radius: f64,
    count: usize,
) -> World {
    let props = phys::chem::analyse(&arrangement).expect("analyses");
    let rho = props.density;
    let mass = rho * 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    let spec = SampleSpec::new(count, Profile::Lattice, MassSpectrum::Equal, BodyKind::Grain);
    let matter = Matter::neutral(mass, radius, 290.0, Composition::primordial());
    let mut w = World::new(Tree::new(seed, matter, Tier::Continuum, spec), 20.0);
    let root = w.tree.root;
    let id = w.substances.intern(arrangement).expect("interns");
    let mut mix = Mixture::new();
    mix.add(id, phase, 1.0);
    w.set_mixture(root, mix);
    w.tree.refine(root);
    w
}

/// **A liquid's modulus is its vaporisation energy density, and a solid's is
/// its stiffness** — the owner's decision for Phase 5's second piece.
///
/// The number that makes the case: water's `cohesive_energy` is the O–H bonds,
/// which through the solid stiffness rule would price it a hundred times too
/// stiff. What Trouton prices is what holds the liquid together.
#[test]
fn water_is_priced_by_what_holds_it_together() {
    let water = substances::water();
    let liquid = Condensed::liquid(&water).expect("water has a liquid");
    let solid_rule = phys::material::dense_stiffness(&water) / 1.2;
    let c = liquid.sound_speed(liquid.rest_density);
    println!(
        "  water: rest {:.1} kg/m^3, K {:.3e} Pa (the bond rule would give {:.3e}), c {:.0} m/s",
        liquid.rest_density, liquid.bulk_modulus, solid_rule, c
    );
    // Real water: 2.2e9 Pa and 1500 m/s. The engine's water is 65% too dense
    // and its modulus high by the same factor, so the sound speed, which is
    // their ratio, is the number to hold.
    assert!((c - 1500.0).abs() / 1500.0 < 0.10, "water carries sound at {c:.0} m/s");
    assert!(solid_rule / liquid.bulk_modulus > 30.0, "the bond rule is the wrong energy");

    // At rest it presses with nothing; one per cent compressed, with K/7 of
    // (1.01^7 - 1) — about K/100.
    assert_eq!(liquid.pressure(liquid.rest_density), 0.0);
    let p1 = liquid.pressure(1.01 * liquid.rest_density);
    println!("  one per cent compressed: {p1:.3e} Pa");
    assert!((p1 / liquid.bulk_modulus - 0.0103).abs() < 0.001);
    // Below rest it is partly empty, not stretched.
    assert_eq!(liquid.pressure(0.5 * liquid.rest_density), 0.0);
    // And the inverse is the inverse.
    let back = liquid.density_at(p1);
    assert!((back / (1.01 * liquid.rest_density) - 1.0).abs() < 1e-12);

    let silica = Condensed::solid(&substances::silica()).expect("silica has a solid");
    println!(
        "  silica: rest {:.1} kg/m^3, K {:.3e} Pa, c {:.0} m/s",
        silica.rest_density,
        silica.bulk_modulus,
        silica.sound_speed(silica.rest_density)
    );
    assert!(silica.bulk_modulus > liquid.bulk_modulus);
}

/// **A bucket of water stops pricing at 4x10^8 Pa.** §4.2's measurement,
/// against the node's own equation of state.
#[test]
fn a_bucket_of_water_is_not_a_gas() {
    let w = ball_of(0xB0C, substances::water_arrangement(), Phase::Liquid, 0.15, 64);
    let root = w.tree.root;
    let m = w.tree.nodes[root.get()].matter;
    let eos = w.eos_of(root);
    let gas = m.pressure();
    let c_gas = m.sound_speed();
    let c = w.sound_speed_of(root);
    println!("  gas law: {gas:.3e} Pa, {c_gas:.0} m/s; its own: {:?}, {c:.0} m/s", eos);
    assert!(matches!(eos, Eos::Condensed(_)), "liquid water is condensed");
    assert!(gas > 1e8, "the gas law's answer is the one §4.2 measured");
    assert!((c - 1500.0).abs() / 1500.0 < 0.10);
}

fn fastest(w: &World) -> f64 {
    w.tree.nodes[w.tree.root.get()].bodies.iter().map(|b| b.vel.norm()).fold(0.0, f64::max)
}

/// **Water handed to SPH stays water.** Under the gas law the same ball of
/// liquid is a gas at 3x10^9 Pa and leaves; under its own equation of state it
/// carries no pressure at rest and moves only as far as its draw was
/// over-dense.
#[test]
fn water_handed_to_sph_is_not_blown_apart() {
    let run = |described: bool| {
        let mut w = ball_of(0xB0D, substances::water_arrangement(), Phase::Liquid, 0.15, 64);
        let root = w.tree.root;
        if !described {
            w.set_mixture(root, Mixture::new());
        }
        let books = |w: &World| -> f64 {
            w.tree.nodes[w.tree.root.get()]
                .bodies
                .iter()
                .map(|b| b.kinetic_energy() + b.internal_energy)
                .sum()
        };
        let before = books(&w);
        let mut covered = 0.0;
        for _ in 0..10 {
            covered += w.advance_node(root, 1.0e-3).dt_used;
        }
        let after = books(&w);
        (fastest(&w), covered, before, after, w.stats.eos_outside_validity)
    };
    let (gas, gcov, ..) = run(false);
    let (liq, lcov, b, a, reports) = run(true);
    let drift = (a - b) / b.abs().max(1e-30);
    println!("  as a gas: fastest {gas:.1} m/s over {gcov:.2e} s");
    println!("  as water: fastest {liq:.1} m/s over {lcov:.2e} s, energy moved {drift:.2e}, {reports} reports");
    assert_eq!(reports, 0, "water priced as water is inside its equation of state");
    assert!(liq < gas, "the liquid must not leave faster than the gas did");
}

/// **Two solids in a room stop detonating.** `docs/BACKLOG.md`'s entry, whose
/// trigger this phase is: a `Continuum` node whose every body is the stand-in
/// for a promoted child was priced as a hot dense gas, and the ball in it went
/// 5.36 -> 7.96 -> 20.9 -> 63.3 -> 194 -> 566 m/s with nothing touching it.
///
/// A stand-in now answers with its child's own equation of state. Told what
/// the two things are made of, neither is anywhere near its rest density at
/// its parent's smoothing length, and a solid below its rest density presses
/// with nothing.
#[test]
fn two_solids_in_a_room_do_not_detonate() {
    let run = |described: bool| {
        let spec = SampleSpec::new(2, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
        let matter = Matter::neutral(48_010.0, 12.0, 290.0, Composition::primordial());
        let mut w = World::new(Tree::new(0xB0F, matter, Tier::Continuum, spec), 20.0);
        let root = w.tree.root;
        w.tree.refine(root);
        {
            let nd = &mut w.tree.nodes[root.get()];
            nd.bodies[0].pos = phys::math::v3(-3.0, 0.0, 0.0);
            nd.bodies[0].vel = phys::math::Vec3::ZERO;
            nd.bodies[0].mass = 48_000.0;
            nd.bodies[1].pos = phys::math::v3(3.0, 0.0, 0.0);
            nd.bodies[1].vel = phys::math::v3(0.0, 5.3572, 0.0);
            nd.bodies[1].mass = 10.0;
        }
        let small = SampleSpec::new(2, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
        let a = w.tree.promote(root, 0, small);
        let b = w.tree.promote(root, 1, small);
        if described {
            let wood = w
                .substances
                .intern(substances::cellulose_arrangement())
                .expect("cellulose analyses");
            let mut mix = Mixture::new();
            mix.add(wood, Phase::Solid, 1.0);
            w.set_mixture(a, mix);
            w.set_mixture(b, mix);
        }
        let mut speeds = Vec::new();
        for _ in 0..6 {
            w.advance_node(root, 0.05);
            speeds.push(w.tree.nodes[b.get()].motion.velocity.norm());
        }
        (speeds, w.stats.eos_outside_validity)
    };
    let (gas, gas_reports) = run(false);
    let (solid, solid_reports) = run(true);
    println!("  undescribed: {:?}, {gas_reports} reports", gas.iter().map(|v| format!("{v:.2}")).collect::<Vec<_>>());
    println!("  two solids:  {:?}, {solid_reports} reports", solid.iter().map(|v| format!("{v:.4}")).collect::<Vec<_>>());
    assert!(gas_reports > 0, "undescribed stand-ins are still outside the gas law, and still say so");
    assert_eq!(solid_reports, 0);
    let worst = solid.iter().map(|v| (v - 5.3572).abs()).fold(0.0, f64::max);
    assert!(worst < 1e-3, "a solid in a room is not accelerated by its room: moved {worst:.2e} m/s");
}

/// The per-phase numbers, against the real ones they should be near.
#[test]
fn conductivity_is_a_law_per_phase() {
    let water = substances::water();
    let silica = substances::silica();
    let liquid = phys::transport::conductivity(&water, Phase::Liquid, 290.0);
    let steam = phys::transport::conductivity(&water, Phase::Gas, 400.0);
    let glass = phys::transport::conductivity(&silica, Phase::Solid, 290.0);
    println!("  water {liquid:.3} (real 0.60), steam {steam:.4} (real 0.027), silica {glass:.2} (fused 1.4, quartz 6-10) W/m/K");
    // Liquid and solid within a factor of a few of the real thing, and a gas
    // an order and more below either — which is the whole of why a duvet works.
    assert!(liquid > 0.2 && liquid < 3.0);
    assert!(glass > 0.5 && glass < 15.0);
    assert!(steam > 0.003 && steam < 0.1);
    assert!(steam * 10.0 < liquid);
}
