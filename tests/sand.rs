//! Phase 5 — what lets sand be sand. `docs/PLAY.md` §7's granular strength,
//! which `docs/BACKLOG.md` measured missing: a poured pile of dry silica came
//! out at 3.26x10^7 Pa in tension against the 4.5x10^3 Pa water at 3 m/s can
//! press with, so a footprint in sand was as permanent as a scar in granite.

use phys::chem::{Mixture, Phase};
use phys::engine::World;
use phys::erode::{capillary_cohesion, loose_grain_threshold, pocket_friction, relaxation_rate};
use phys::ids::NodeIdx;
use phys::material::substances;
use phys::math::v3;
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::Tree;
use phys::units::Tier;

/// A 1 m patch of silicate ground on an Earth's +z pole, holding `water` kg of
/// water in its pores (or over it) beside 1500 kg of sand. Emplaced, so
/// nothing ever saw it freeze.
fn ground(seed: u64, water: f64) -> (World, NodeIdx) {
    const R: f64 = 6.371e6;
    let spec = SampleSpec::new(8, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let planet = Matter::neutral(5.97e24, R, 290.0, Composition::crustal());
    let mut w = World::new(Tree::new(seed, planet, Tier::Planetary, spec), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let patch_spec = SampleSpec::new(400, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let patch = w.tree.promote(root, 0, patch_spec);
    let silica = w.substances.intern(substances::silica_arrangement()).unwrap();
    let h2o = w.substances.intern(substances::water_arrangement()).unwrap();
    let sand = 1500.0;
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, sand / (sand + water));
    if water > 0.0 {
        mix.add(h2o, Phase::Liquid, water / (sand + water));
    }
    let composition = mix.composition(&w.substances).0;
    {
        let n = &mut w.tree.nodes[patch.get()];
        n.matter = Matter::neutral(sand + water, 1.0, 290.0, composition);
        n.motion.offset = v3(0.0, 0.0, R + 1.0);
        n.spec = patch_spec;
    }
    w.set_mixture(patch, mix);
    w.emplace(patch, phys::morph::Program::Terrain, sand, None);
    w.tree.refine(patch);
    (w, patch)
}

/// The pocket a grain sits in, from packing geometry and nothing else.
#[test]
fn a_loose_grain_sits_in_a_pocket_its_packing_sets() {
    let close = pocket_friction(phys::erode::CLOSE_PACKING);
    let loose = pocket_friction(phys::sampler::RANDOM_LOOSE_PACKING);
    println!(
        "  close-packed pocket {:.4} ({:.2} deg), random loose {:.4} ({:.2} deg)",
        close,
        close.atan().to_degrees(),
        loose,
        loose.atan().to_degrees()
    );
    assert!((close - 1.0 / (2.0 * 2.0f64.sqrt())).abs() < 1e-12, "three touching spheres");
    assert!(loose > close, "a looser pile's pockets are deeper");
}

/// **Sand is held by its weight, and damp sand by its menisci**; neither by a
/// bond. What a flux has to press with to take a grain, on the engine's own
/// silicate ground, dry in air, damp in air and drowned — and the flow speed
/// each one takes, against granite's.
#[test]
fn what_holds_a_grain_of_sand_is_not_a_bond() {
    let (dry_w, dry) = ground(0x5A4D, 0.0);
    let (damp_w, damp) = ground(0x5A4D, 60.0);
    let (wet_w, wet) = ground(0x5A4D, 400.0);
    let dry_m = dry_w.material_of(dry).unwrap();
    let grain = dry_m.flaw_size;
    let bonded = dry_m.strength_of(grain, 290.0);
    let air = 1.225;
    let water = substances::water();
    let water_rho = phys::eos::Condensed::liquid(&water).unwrap().rest_density;
    let gamma = phys::erode::surface_tension(&water);
    let dry_s = dry_w.loose_grain_strength(dry, &dry_m, air).expect("emplaced ground is loose");
    let damp_s = damp_w.loose_grain_strength(damp, &damp_w.material_of(damp).unwrap(), air).unwrap();
    let wet_s = wet_w.loose_grain_strength(wet, &wet_w.material_of(wet).unwrap(), water_rho).unwrap();
    let speed = |rho: f64, sigma: f64| (2.0 * sigma / rho).sqrt();
    println!("  grain {grain:.3e} m (the deposited flaw scale); water's surface tension {gamma:.3} J/m^2 (real 0.072)");
    println!("  bonded, as before       {bonded:.3e} Pa");
    println!("  dry, in air             {dry_s:.3e} Pa, lifted by wind at {:.2} m/s", speed(air, dry_s));
    println!("  damp, in air            {damp_s:.3e} Pa, lifted by wind at {:.1} m/s", speed(air, damp_s));
    println!("  drowned, under water    {wet_s:.3e} Pa, lifted by a current of {:.3} m/s", speed(water_rho, wet_s));
    assert!(bonded > 1e6 * dry_s, "a loose grain is not held by the bond Griffith prices");
    assert!(damp_s > 10.0 * dry_s, "menisci hold damp sand far harder than its weight does");
    assert!(wet_s < dry_s, "under water a grain weighs less and has no menisci");
    assert!(speed(water_rho, wet_s) < 0.3, "a tidal current moves drowned sand");
    assert!(speed(air, damp_s) > 10.0, "a breeze does not strip a damp beach");
    let _ = capillary_cohesion(gamma, grain, 0.6);
    let _ = loose_grain_threshold(2644.0, air, 9.8, grain, 0.6);
}

/// §5.2's honesty test, finally: **one expression gives the wet-sand case and
/// the granite case** with only material and flux differing. The same mark in
/// the same current relaxes in hours on drowned sand and never on rock the
/// engine saw freeze.
#[test]
fn a_squiggle_in_drowned_sand_goes_and_one_in_granite_stays() {
    let (w, wet) = ground(0x5A4E, 400.0);
    let water_rho = phys::eos::Condensed::liquid(&substances::water()).unwrap().rest_density;
    let m = w.material_of(wet).unwrap();
    let sand = w.loose_grain_strength(wet, &m, water_rho).unwrap();
    let granite = m.strength_of(m.flaw_size, 290.0);
    let (current, span) = (0.3, 0.05);
    let on_sand = relaxation_rate(sand, water_rho, current, span);
    let on_rock = relaxation_rate(granite, water_rho, current, span);
    println!(
        "  a 5 cm squiggle in a 0.3 m/s current: sand relaxes at {on_sand:.3e} /s (half-life {:.1} s), \
         bonded rock at {on_rock:.3e} /s",
        std::f64::consts::LN_2 / on_sand
    );
    assert!(on_sand > 0.0 && std::f64::consts::LN_2 / on_sand < 6.0 * 3600.0, "gone well inside a tide");
    assert_eq!(on_rock, 0.0, "and the rock is not touched");
}
