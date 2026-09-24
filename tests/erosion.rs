//! A deviation decays, and forgetting is what the decay reaches.
//!
//! `docs/PLAY.md` §5. The engine's fourth axiom is implemented as a binary:
//! anything touched is pinned and persisted, outright and for ever. §5.1 says
//! that at play resolution the binary is wrong, because it makes a footprint as
//! permanent as a felled trunk, and that the fix is to make the exemption a
//! *decay* — a deviation carries an amplitude, the amplitude falls at a rate
//! the physics sets, and forgetting becomes the deviation reaching zero.
//!
//! These are the two halves of Phase 4's second done-when: a mark drawn in a
//! soft surface goes, a cut that carried material away stays, and neither of
//! them was tagged as important when it was made. What separates them is a
//! measurement over the conserved set and nothing else.

use phys::chem::{Mixture, Phase};
use phys::engine::{default_spec, World};
use phys::erode::{relaxation_rate, Deviation};
use phys::ids::NodeIdx;
use phys::morph::Program;
use phys::observe::Interaction;
use phys::recipe::Recipe;
use phys::state::{Composition, Matter};
use phys::tree::Tree;
use phys::units::*;

/// A patch of ground that a flux is running over.
///
/// Ice a whisker under its own melting point, with meltwater over it: a surface
/// the engine can honestly say is soft, which is what stands in for wet sand
/// until an uncemented aggregate has a derived strength of its own. See
/// `phys::erode`'s module documentation for the measurement that says why.
fn a_soft_flat(seed: u64, softness: f64, pace: f64) -> (World, NodeIdx) {
    let spec = default_spec(Tier::Continuum);
    let ground = Matter::neutral(2.0e9, 5.0, 200.0, Composition::primordial());
    let mut w = World::new(Tree::new(seed, ground, Tier::Continuum, spec), 20.0);
    w.pace_fixed(pace);
    let patch = w.tree.root;
    let ice = w
        .substances
        .intern(phys::material::substances::water_arrangement())
        .expect("water analyses");
    let mut mix = Mixture::new();
    mix.add(ice, Phase::Solid, 1.0);
    w.set_mixture(patch, mix);
    // The tide running over the flat. Stated by the scenario, the way a
    // mixture is: this is a place with water moving across it, and the engine
    // measures everything else from there.
    w.plant(
        patch,
        Program::Terrain,
        Some(phys::morph::Environment {
            fluid_density: 1000.0,
            flow_speed: FLOW,
            ..Default::default()
        }),
    );

    // **How warm "soft" is, solved for rather than stated.** The engine has no
    // uncemented aggregate, so the only surface it can honestly call soft is
    // one near its own softening point — and where that is is a property of
    // what the surface is made of, which is exactly the thing a test may not
    // hard-code. So: take the material the node measures for itself, ask what
    // stress the flux it is in presses with, and put the temperature where the
    // material's own strength is `softness` times that.
    let material = w.material_of(patch).expect("the ground measures a material");
    let press = 0.5 * w.environment_at(patch).fluid_density * FLOW * FLOW;
    let nominal = material.strength_of(material.flaw_size, 0.0).max(1e-30);
    let want = (softness * press / nominal).clamp(0.0, 1.0);
    let t = material.thermal_onset
        + (1.0 - want) * (material.thermal_gone - material.thermal_onset);
    w.tree.nodes[patch.get()].matter.set_temperature(t);
    (w, patch)
}

/// How fast the water is moving over it, m/s.
const FLOW: f64 = 3.0;

fn field_of(w: &World, idx: NodeIdx) -> Vec<Deviation> {
    w.tree.nodes[idx.get()]
        .morphology
        .as_ref()
        .map(|m| m.field.clone())
        .unwrap_or_default()
}

/// One expression, three inputs, and the geometry is one of them.
#[test]
fn the_rate_is_the_material_the_flux_and_the_feature() {
    let soft = 2.0e2;
    let water = (1000.0, 1.0);
    // Geometry: a narrow notch has a steeper gradient than a broad dish and
    // goes faster, at the same material and the same flux.
    let narrow = relaxation_rate(soft, water.0, water.1, 0.01);
    let broad = relaxation_rate(soft, water.0, water.1, 1.0);
    println!("  a 1 cm feature relaxes at {narrow:.4e} /s, a 1 m one at {broad:.4e}");
    assert!(narrow > broad * 50.0, "geometry has to matter: {narrow:.3e} against {broad:.3e}");

    // Flux: the same feature in the same material, under still air, does not
    // move at all, because nothing reaches the threshold.
    let still = relaxation_rate(soft, 1.225, 0.5, 0.01);
    println!("  and under still air, {still:.4e} /s");
    assert_eq!(still, 0.0, "a flux that cannot lift a grain lifts none");

    // Material: rock under the same water does not move either, and that is
    // the granite half of the honesty test.
    let granite = relaxation_rate(9.72e7, water.0, water.1, 0.01);
    println!("  and a scar in rock at 9.72e7 Pa, under the same water, {granite:.4e} /s");
    assert_eq!(granite, 0.0);
}

/// A mark that moved nothing decays away, and the node forgets it.
#[test]
fn a_mark_that_took_nothing_is_forgotten() {
    // A whisker under melting: the surface is soft, and the engine says so
    // from the substance's own softening range rather than being told.
    // A tenth as strong as the flux presses: soft ground with meltwater on it.
    let (mut w, patch) = a_soft_flat(0xE401, 0.1, 1.0);
    let before = w.tree.nodes[patch.get()].matter.mass;
    w.interact(Interaction::Mark {
        target: patch,
        deviation: Deviation::reshaped(0.2, -0.1, 0.05, -0.02),
    });
    assert_eq!(field_of(&w, patch).len(), 1, "the mark is in the field");
    let drawn = field_of(&w, patch)[0].amplitude;

    let mut gone_at = None;
    for f in 0..4000 {
        w.step_frame(20_000.0);
        if field_of(&w, patch).is_empty() {
            gone_at = Some((f, w.time));
            break;
        }
    }
    let (frame, when) = gone_at.expect("the mark never went");
    let after = w.tree.nodes[patch.get()].matter.mass;
    let env = w.environment_at(patch);
    let material = w.material_of(patch).expect("a material");
    let cohesion =
        material.strength_of(material.flaw_size, w.tree.nodes[patch.get()].matter.temperature);
    let rate = relaxation_rate(cohesion, env.fluid_density, env.flow_speed, 0.05);
    println!(
        "  a 0.10 m mark drawn 2 cm deep, in ground of {cohesion:.4e} Pa under {:.0} kg/m^3 \
         at {:.1} m/s: {rate:.4e} /s, half-life {:.4} s",
        env.fluid_density,
        env.flow_speed,
        std::f64::consts::LN_2 / rate
    );
    println!(
        "  gone after {:.1} s ({frame} frames), and the patch weighs {after:.6e} kg \
         against {before:.6e}",
        when
    );
    assert!(w.tree.stats.deviations_forgotten >= 1, "nothing was forgotten");
    assert!(
        (after - before).abs() <= 1e-12 * before,
        "forgetting is not allowed to be a source or a sink"
    );
    let _ = drawn;
}

/// A cut that carried material away is still there a year later.
///
/// Nothing tagged it. The only difference between this and the mark above is
/// that the node weighs less than it did, which is a measurement over the
/// conserved set — §5.3 and §5.8.
#[test]
fn a_cut_that_took_mass_is_never_forgotten() {
    let (mut w, patch) = a_soft_flat(0xE402, 0.1, YEAR / 500.0);
    let before = w.tree.nodes[patch.get()].matter.mass;
    // A channel: 0.4 m across, 0.3 m deep, and the spoil is gone from the node.
    let volume = std::f64::consts::PI * 0.2 * 0.2 * 0.3;
    let carried = volume * 917.0;
    w.interact(Interaction::Mark {
        target: patch,
        deviation: Deviation::cut(-0.3, 0.4, 0.2, -0.3, carried),
    });
    let after_cut = w.tree.nodes[patch.get()].matter.mass;
    println!(
        "  cutting it took {carried:.3} kg out of the patch: {before:.4e} -> {after_cut:.4e} kg"
    );
    assert!(carried > 0.0);
    assert!((before - after_cut - carried).abs() <= 1e-9 * before);

    // A year of the same flux that removed the mark above in seconds.
    let mut years = 0.0;
    while w.time < YEAR {
        w.step_frame(20_000.0);
        years = w.time / YEAR;
    }
    let field = field_of(&w, patch);
    println!(
        "  and after {years:.2} years of the same water, the channel is still {:.4} m deep",
        -field.first().map(|d| d.amplitude).unwrap_or(0.0)
    );
    assert_eq!(field.len(), 1, "a deviation that took mass away may not be dropped");
    assert!(
        (field[0].amplitude + 0.3).abs() < 1e-12,
        "and it may not decay either: {:.6} m",
        field[0].amplitude
    );
    let now = w.tree.nodes[patch.get()].matter.mass;
    assert!(
        (now - after_cut).abs() <= 1e-9 * after_cut,
        "and the mass it took never came back: {now:.6e} against {after_cut:.6e}"
    );
}

/// The field is in the ground, not beside it: what the recipe draws moves.
#[test]
fn a_mark_changes_what_the_ground_draws() {
    // Two worlds from one seed: the same ground, and the same ground marked.
    // `refine` returns what a node already holds, so asking one world twice
    // would measure nothing — and this way the comparison also says that a
    // patch nobody marked draws exactly as it always did.
    let (mut plain_world, plain_patch) = a_soft_flat(0xE403, 1.0e-6, 1.0);
    let (mut w, patch) = a_soft_flat(0xE403, 1.0e-6, 1.0);
    plain_world.tree.refine(plain_patch);
    let plain: Vec<f64> = plain_world.tree.nodes[plain_patch.get()]
        .bodies
        .iter()
        .map(|b| b.pos.z)
        .collect();
    assert!(!plain.is_empty(), "the patch drew nothing");

    w.interact(Interaction::Mark {
        target: patch,
        deviation: Deviation::reshaped(0.0, 0.0, 1.0, -0.5),
    });
    w.tree.refine(patch);
    let marked: Vec<f64> = w.tree.nodes[patch.get()].bodies.iter().map(|b| b.pos.z).collect();
    assert_eq!(plain.len(), marked.len(), "the mark changed the cell count");
    let worst = plain
        .iter()
        .zip(marked.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    println!("  a half-metre mark moved the ground by {worst:.4} m where it was drawn");
    assert!(worst > 0.1, "the mark is not in what the ground draws: {worst:.6} m");

    // And the recipe is still a rule: the mark costs five numbers, not a
    // heightmap.
    let bytes = w.tree.nodes[patch.get()].morphology.as_ref().unwrap().state_bytes();
    println!("  and the patch is {bytes} B of rule with the mark in it");
    assert!(bytes < 512);
    assert!(matches!(
        w.tree.nodes[patch.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()),
        Some(Recipe::Tiled(_))
    ));
}
