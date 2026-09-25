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
    // What it is made of, elementally, from what it is made of: a body's heat
    // capacity counts its particles off the composition, and a ball of
    // silicate declared primordial holds seventeen times the heat it should.
    // Rebuilt rather than patched, so the internal energy is the one this
    // composition holds at this temperature.
    let composition = mix.composition(&w.substances).0;
    let m = w.tree.nodes[root.get()].matter;
    w.tree.nodes[root.get()].matter = Matter::neutral(m.mass, m.radius, m.temperature, composition);
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

/// An open oak box — a floor and four sides — standing on an Earth's +z pole,
/// with `water` kg of water in it. Its own -z is down.
fn a_bucket(seed: u64, half: f64, water: f64, count: usize) -> World {
    use phys::assembly::{Assembly, Join};
    use phys::math::v3;
    const R: f64 = 6.371e6;
    let spec = SampleSpec::new(8, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let planet = Matter::neutral(5.97e24, R, 290.0, Composition::crustal());
    let mut w = World::new(Tree::new(seed, planet, Tier::Planetary, spec), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let oak = w.substances.intern(substances::cellulose_arrangement()).unwrap();
    let h2o = w.substances.intern(substances::water_arrangement()).unwrap();
    let t = 0.05 * half;
    let faces = [
        (v3(0.0, 0.0, -half), v3(half, half, t)),
        (v3(-half, 0.0, 0.0), v3(t, half, half)),
        (v3(half, 0.0, 0.0), v3(t, half, half)),
        (v3(0.0, -half, 0.0), v3(half, t, half)),
        (v3(0.0, half, 0.0), v3(half, t, half)),
    ];
    let mut parts = Vec::new();
    let mut joins = Vec::new();
    for (i, (at, h)) in faces.iter().enumerate() {
        parts.push(phys::state::Body::solid(*at, *h, 8.0 * h.x * h.y * h.z * 700.0, oak, i as u32));
        joins.push(if i == 0 { Join::NONE } else { Join::new(0, oak, 2.0 * half * 2.0 * t) });
    }
    let parts = Assembly::new(parts, joins);
    let wood = parts.mass();
    let patch_spec = SampleSpec::new(count, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let bucket = w.tree.promote(root, 0, patch_spec);
    {
        let n = &mut w.tree.nodes[bucket.get()];
        n.matter = Matter::neutral(wood + water, half, 290.0, Composition::organic());
        n.motion.offset = v3(0.0, 0.0, R + 2.0 * half);
        n.motion.velocity = phys::math::Vec3::ZERO;
        n.spec = patch_spec;
    }
    let mut mix = Mixture::new();
    mix.add(oak, Phase::Solid, wood / (wood + water));
    mix.add(h2o, Phase::Liquid, water / (wood + water));
    let composition = mix.composition(&w.substances).0;
    w.tree.nodes[bucket.get()].matter = Matter::neutral(wood + water, half, 290.0, composition);
    w.set_mixture(bucket, mix);
    w.assemble(bucket, phys::morph::Program::Tree, parts, None);
    w.tree.refine(bucket);
    w
}

/// **A free surface**: `docs/PLAY.md` §4.3's first piece, which "everything
/// else here waits on".
///
/// Nothing stores a water level. The node's matter says how much water there
/// is, its recipe says what shape the bucket is, and the draw stands the water
/// on the floor in columns at the spacing its own density gives it until it
/// runs out — so the top is level and its height is a consequence. Then it is
/// left alone in the field it is in, held by the bucket, and has to stay: not
/// sinking into the floor, not leaving, and settling rather than ringing.
#[test]
fn water_lies_level_in_a_bucket_and_stays_there() {
    let mut w = a_bucket(0x3E7, 0.5, 250.0, 800);
    let key = w.tree.nodes.iter().find(|n| n.alive && n.morphology.is_some()).unwrap().key;
    let node = |w: &World| {
        phys::ids::NodeIdx(w.tree.nodes.iter().position(|n| n.alive && n.key == key).unwrap() as u32)
    };
    let split = |w: &World| {
        let n = &w.tree.nodes[node(w).get()];
        let mask = n.structural_mask().expect("the bucket is a structure");
        let pick = |want: bool| -> Vec<phys::state::Body> {
            n.bodies.iter().zip(mask.iter()).filter(|(_, o)| **o == want).map(|(b, _)| *b).collect()
        };
        (pick(false), pick(true))
    };
    let (water, wood) = split(&w);
    let spacing = (water[0].mass / w.tree.nodes[node(&w).get()].rest_density).cbrt();
    // The floor is the first part; its top is where the water stands.
    let floor = wood[0].pos.z + wood[0].half.z;
    let heights = |water: &[phys::state::Body]| {
        let lo = water.iter().map(|b| b.pos.z).fold(f64::INFINITY, f64::min);
        let hi = water.iter().map(|b| b.pos.z).fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    };
    let (lo, hi) = heights(&water);
    let top: Vec<f64> = water
        .iter()
        .filter(|a| {
            !water.iter().any(|b| {
                let d = b.pos - a.pos;
                d.z > 0.5 * spacing && (d.x * d.x + d.y * d.y).sqrt() < 0.5 * spacing
            })
        })
        .map(|b| b.pos.z)
        .collect();
    let spread = top.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - top.iter().cloned().fold(f64::INFINITY, f64::min);
    println!(
        "  {} parcels {spacing:.4} m apart, standing {:.4} spacings above the floor; \
         the top is level to {spread:.2e} m, {:.4} m above the floor",
        water.len(),
        (lo - floor) / spacing,
        hi - floor
    );
    assert!(water.len() > 50, "the water is drawn as parcels");
    assert!(((lo - floor) / spacing - 0.5).abs() < 1e-6, "the bottom layer stands on the floor");
    assert!(spread < 1e-9, "the free surface is not level");

    // Left alone for three seconds of world time, in the field it is in.
    let (mut covered, mut solves, mut peak) = (0.0, 0, 0.0f64);
    while covered < 3.0 && solves < 1200 {
        covered += w.advance_node(node(&w), 0.05).dt_used;
        solves += 1;
        let (water, _) = split(&w);
        peak = peak.max(water.iter().map(|b| b.vel.norm()).fold(0.0, f64::max));
    }
    let (water, _) = split(&w);
    let fastest = water.iter().map(|b| b.vel.norm()).fold(0.0, f64::max);
    let (lo_after, hi_after) = heights(&water);
    let out = water.iter().filter(|b| b.pos.x.abs() > 0.5 || b.pos.y.abs() > 0.5).count();
    println!(
        "  after {covered:.3} s in {solves} solves: fastest {fastest:.4} m/s (peak {peak:.4}), \
         bottom {:.4} -> {:.4} spacings above the floor, top {:.4} -> {:.4} m, {out} out",
        (lo - floor) / spacing,
        (lo_after - floor) / spacing,
        hi - floor,
        hi_after - floor
    );
    assert!(covered >= 3.0, "the solver could not cover the span: {covered}");
    assert_eq!(out, 0, "water left the bucket");
    assert!(lo_after - floor > 0.25 * spacing, "the water sank into the floor");
    assert!(fastest < 0.1 * peak.max(0.01), "the water is still ringing at {fastest} m/s");
}
/// **Two touching things reach the same temperature** — `docs/BACKLOG.md`'s
/// trigger for a conductive coefficient, and the reason Phase 5 carries one.
///
/// Two neighbouring parcels of one block of silicate, one at 900 K and one at
/// 300 K. With radiation alone they share heat at
/// `sigma A (T_a + T_b)(T_a^2 + T_b^2)`; with conduction as well, at
/// Einstein–Cahill–Pohl's rate through the face they share. The ratio of the
/// two conductances is the claim, and the equilibration is the consequence.
#[test]
fn two_touching_things_reach_the_same_temperature() {
    let mut w = ball_of(0xC0D, substances::silica_arrangement(), Phase::Solid, 0.10, 64);
    let root = w.tree.root;
    w.tree.pin(root);
    let nb = w.tree.neighbourhood(root);
    let pairs = nb.pairs(nb.resolution()).expect("in reach");
    assert!(!pairs.is_empty(), "the parcels of a solid block have neighbours");
    let (i, j) = pairs[0];
    let (_, pi, ri) = nb.at(i).unwrap();
    let (_, pj, rj) = nb.at(j).unwrap();
    let d = (pj - pi).norm();

    let props = substances::silica();
    let k = phys::transport::conductivity(&props, Phase::Solid, 300.0);
    let conductive = phys::transport::continuum_conductance(k, k, d);
    let area = phys::neighbourhood::radiative_area(ri, rj, d);
    let radiative = phys::neighbourhood::radiative_conductance(900.0, 300.0, area);
    println!(
        "  silica conducts at {k:.2} W/m/K; neighbours {d:.4} m apart: conduction {conductive:.3e} W/K, \
         radiation {radiative:.3e} W/K, {:.0}x",
        conductive / radiative
    );
    assert!(conductive > 10.0 * radiative, "touching, conduction is what carries the heat");

    // The solver heats every parcel by viscosity as it moves, and that is not
    // the exchange. A condensed pressure does not depend on temperature, so a
    // control with the hot parcel at 300 K follows the same trajectories to
    // the bit, and the difference between the two is the heat that crossed.
    let cold_after = |hot: f64| {
        let mut w = ball_of(0xC0D, substances::silica_arrangement(), Phase::Solid, 0.10, 64);
        let root = w.tree.root;
        w.tree.pin(root);
        {
            let b = &mut w.tree.nodes[root.get()].bodies;
            for x in b.iter_mut() {
                x.temperature = 300.0;
            }
            b[i].temperature = hot;
        }
        let span = w.advance_node(root, 1.0e-3).dt_used;
        let b = &w.tree.nodes[root.get()].bodies;

        (b[j].temperature, span, b[i].heat_capacity(), b[j].heat_capacity())
    };
    let (warmed, span, c_i, c_j) = cold_after(900.0);
    let (control, _, _, _) = cold_after(300.0);
    let gained = (warmed - control) * c_j;
    let by_radiation = radiative * 600.0 * span;
    let by_both = (radiative + conductive) * 600.0 * span;
    println!(
        "  over {span:.3e} s the cold parcel took {gained:.3e} J: radiation alone \
         would give {by_radiation:.3e}, with conduction {by_both:.3e}"
    );
    assert!(
        (gained / by_both - 1.0).abs() < 0.05,
        "the heat that crossed is not what the conductance says"
    );

    // And the time they take to meet, which is the diffusion time of the rock
    // and not a property of the frame: `C_h / G` for the pair.
    let tau = c_i * c_j / (c_i + c_j) / (radiative + conductive);
    let reached = phys::neighbourhood::exchange(
        phys::neighbourhood::Reservoir::new(900.0, c_i),
        phys::neighbourhood::Reservoir::new(300.0, c_j),
        radiative + conductive,
        5.0 * tau,
    );
    let (hot, cold) = (900.0 - reached / c_i, 300.0 + reached / c_j);
    println!("  they meet with a time constant of {tau:.0} s; after five: {hot:.2} K and {cold:.2} K");
    assert!((hot - cold).abs() < 0.01 * 600.0);
    // A rock's diffusivity is about 1e-6 m^2/s, so 3 cm takes minutes.
    assert!(tau > 30.0 && tau < 3000.0, "{tau} s is not a rock's diffusion time");
}

