//! Forests, deserts and winters — none of which are written down anywhere.
//!
//! There is no `Biome` type in this engine and there is not going to be one.
//! A forest is what a patch of ground becomes when its environment supports
//! trees; a desert is the same laws not firing. The test that matters is
//! therefore not "does a forest appear" but "do two patches differing only in
//! their contents end up different", because that is the difference between a
//! world with rules and a world with a lookup table.
//!
//! Everything here goes through machinery written for something else:
//! `chem::react` was written for salt dissolving in a beaker, `environment_at`
//! for shading a tree under its parent, and `evolve_matter` knows only how to
//! radiate and absorb. Nothing in any of them has heard of snow.

use phys::chem::{Arrangement, Bond, Element, Mixture, Order, Phase};
use phys::engine::World;
use phys::ids::NodeIdx;
use phys::morph::{Environment, Program};
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Matter};
use phys::tree::Tree;
use phys::units::*;

const PATCH_MASS: f64 = 2.4e9;

fn ground() -> Tree {
    let matter = Matter::neutral(PATCH_MASS, 500.0, 288.0, Program::Terrain.substrate());
    Tree::new(
        0xB10E,
        matter,
        Tier::Continuum,
        SampleSpec {
            count: 64,
            profile: Profile::Uniform,
            spectrum: MassSpectrum::Equal,
            kind: BodyKind::Grain,
            composition_scatter: 0.0,
            turbulent_fraction: 0.0,
        },
    )
}

/// Water, built from first principles like everything else in the chemistry
/// layer. The engine is never told which substance this is.
fn water(w: &mut World) -> phys::chem::SubstanceId {
    w.substances
        .intern(Arrangement::molecule(
            vec![Element(8), Element(1), Element(1)],
            vec![Bond::new(0, 1, Order::Single), Bond::new(0, 2, Order::Single)],
        ))
        .expect("water analyses")
}

fn a_patch(temperature: f64) -> World {
    let mut w = World::new(ground(), 20.0);
    let root = w.tree.root;
    w.tree.nodes[root.get()].matter.temperature = temperature;
    // `step_frame` takes a wall-clock budget in microseconds, not a span: how
    // much *world* time a frame covers is the pace, and by default that follows
    // whatever is being watched. A 500 m patch of ground is watched at about a
    // second and a half a frame, which is the right answer for looking at it
    // and useless for asking what it does over a season. So the clock is driven
    // by hand here — a day a frame — which is what `PaceMode::Fixed` is for.
    w.pace_fixed(86_400.0);
    w
}

/// Two patches, identical but for what they are made of.
///
/// The wet one sustains more than the dry one. Nothing in the engine branches on a
/// biome — `environment_at` measures the *liquid phase* of whatever the node
/// contains, because what a growing thing needs is a mobile solvent, and the
/// engine deliberately does not know which substance is water.
#[test]
fn a_desert_is_a_forest_that_had_no_water() {
    let mut wet = a_patch(288.0);
    let root = wet.tree.root;
    let h2o = water(&mut wet);
    let mut mix = Mixture::new();
    mix.add(h2o, Phase::Liquid, 0.2);
    wet.set_mixture(root, mix);

    let mut dry = a_patch(288.0);
    let dry_root = dry.tree.root;
    let h2o_dry = water(&mut dry);
    let mut dust = Mixture::new();
    // The same substance, in the phase a desert holds it in: locked in the
    // rock, not available to anything.
    dust.add(h2o_dry, Phase::Solid, 0.2);
    dry.set_mixture(dry_root, dust);

    let grown = |w: &World, r: NodeIdx| w.tree.nodes[r.get()].morphology.as_ref().unwrap().built;

    let we = wet.environment_at(root);
    let de = dry.environment_at(dry_root);
    println!(
        "  wet patch: water {:.3}, thermal factor {:.3}\n  dry patch: water {:.3}, thermal factor {:.3}",
        we.water,
        we.thermal_factor(),
        de.water,
        de.thermal_factor()
    );
    assert!(we.water > 0.1, "the wet patch measured no solvent");
    assert!(de.water < 1e-9, "the dry patch measured solvent it does not have");

    // And it shows in what grows. Same program, same seed, same span, and
    // `None` for the environment so both feel the physics rather than an
    // authored situation — which is the whole point of the comparison.
    //
    // Emplaced as an established tree rather than planted as a seed: a
    // one-kilogram seedling has almost no crown to catch light with, so over a
    // few months neither patch would move enough to tell them apart. This is
    // asking whether an existing tree keeps growing, which is what a forest and
    // a desert actually differ by.
    wet.emplace(root, Program::Tree, 1.0e3, None);
    dry.emplace(dry_root, Program::Tree, 1.0e3, None);
    let start = grown(&wet, root);
    // Three hundred days of weather.
    for _ in 0..300 {
        wet.step_frame(50_000.0);
        dry.step_frame(50_000.0);
    }
    println!(
        "  from {start:.4e} kg of wood after 300 days: wet {:.6e} kg, dry {:.6e} kg",
        grown(&wet, root),
        grown(&dry, dry_root)
    );
    // Both decline, and that is not a failure — at 200 W/m^2 a standing tonne
    // of wood costs more in respiration than its crown can capture, so both
    // patches are above their carrying capacity and shrinking toward it. The
    // carrying capacity is emergent: capture scales with area and upkeep with
    // mass, so the two balance at a finite size that nothing imposes. What the
    // water changes is *where* that balance sits.
    assert!(
        grown(&wet, root) > grown(&dry, dry_root),
        "the dry patch sustains as much wood as the wet one: {:.6e} against {:.6e}",
        grown(&dry, dry_root),
        grown(&wet, root)
    );
}

/// A node radiates, and the energy it radiates is gone.
///
/// The one law `evolve_matter` applies. Before it existed, `matter.luminosity`
/// was computed from Stefan-Boltzmann and *read* — for illumination, for flux
/// at an observer — and never spent, so every star in the world shone for free.
#[test]
fn a_node_left_in_the_dark_cools() {
    let mut w = a_patch(320.0);
    let root = w.tree.root;
    // Nothing shining on it: an unparented root gets the default flux, so make
    // the darkness explicit.
    w.environments.insert(
        w.tree.nodes[root.get()].key,
        Environment { light_flux: 0.0, ..Default::default() },
    );

    let before = w.tree.nodes[root.get()].matter.temperature;
    for _ in 0..50 {
        w.step_frame(50_000.0);
    }
    let after = w.tree.nodes[root.get()].matter.temperature;
    println!(
        "  {before:.2} K -> {after:.2} K, radiating {:.3e} J away in total",
        w.stats.radiated
    );
    assert!(after < before, "a node in the dark did not cool");
    assert!(w.stats.radiated > 0.0, "nothing was radiated");
    assert!(
        w.stats.radiation_deficit == 0.0,
        "a node radiated energy it did not hold: {:.3e} J",
        w.stats.radiation_deficit
    );
    assert!(after >= 2.725, "cooled below the microwave background");
}

/// Winter and spring, neither of which is written anywhere.
///
/// Driven by *light*, because with `evolve_matter` running nothing else would
/// be honest: setting a node's temperature by hand no longer holds, since the
/// node radiates the heat away again within a few frames. So the only input is
/// how brightly the parent shines, and everything downstream of that —
/// temperature, phase, whether there is any water to be had — follows.
///
/// The chain, none of which mentions a season:
///   `evolve_matter` absorbs less than it radiates, so the patch cools
///   → `chem::react` moves its liquid into the solid phase as it passes the
///     melting point it derived for itself
///   → `environment_at` measures no mobile solvent
///   → growth stops.
/// Turn the light back up and every step runs the other way.
#[test]
fn a_patch_freezes_when_the_light_goes_and_thaws_when_it_returns() {
    // A star with a patch of ground in orbit, so the light is derived from
    // something shining rather than authored. `environment_at` turns the
    // parent's luminosity into a flux at the child's distance.
    let matter = Matter::neutral(2.0e30, 7.0e8, 5800.0, Program::Terrain.substrate());
    let tree = Tree::new(
        0x5EA_50,
        matter,
        Tier::Stellar,
        SampleSpec {
            count: 8,
            profile: Profile::Uniform,
            spectrum: MassSpectrum::Equal,
            kind: BodyKind::Grain,
            composition_scatter: 0.0,
            turbulent_fraction: 0.0,
        },
    );
    let mut w = World::new(tree, 20.0);
    w.pace_fixed(86_400.0);
    let star = w.tree.root;
    w.tree.refine(star);
    let patch = w.tree.promote(star, 0, w.tree.nodes[star.get()].spec);
    assert!(!patch.is_none());
    // Promoting a body out of a star gives a 10^29 kg lump of star, which is
    // not a patch of ground. Test setup, stated plainly rather than smuggled:
    // the node is overwritten with the patch this test is about, and put an
    // astronomical unit away so the flux at it is a flux and not an immersion.
    {
        let n = &mut w.tree.nodes[patch.get()];
        n.matter = Matter::neutral(PATCH_MASS, 500.0, 280.0, Program::Terrain.substrate());
        n.motion.offset = phys::math::v3(1.496e11, 0.0, 0.0);
        n.motion.velocity = phys::math::Vec3::ZERO;
    }

    let h2o = water(&mut w);
    let mut mix = Mixture::new();
    mix.add(h2o, Phase::Liquid, 0.3);
    w.set_mixture(patch, mix);

    let liquid = |w: &World| w.mixture_of(patch).in_phase(Phase::Liquid);
    let solid = |w: &World| w.mixture_of(patch).in_phase(Phase::Solid);
    let temp = |w: &World| w.tree.nodes[patch.get()].matter.temperature;
    // The star is dimmed by *cooling* it, not by assigning a luminosity:
    // `evolve_matter` derives luminosity from temperature every frame, so a
    // hand-set value is overwritten before it is ever read. That is the right
    // behaviour and it is worth stating, because assigning the luminosity is
    // the obvious thing to try and it silently does nothing.
    let sun = |w: &mut World, t: f64| {
        let n = &mut w.tree.nodes[star.get()];
        n.matter.temperature = t;
        n.matter.luminosity = phys::state::stefan_boltzmann(n.matter.radius, t);
    };

    // Summer: bright enough that the patch settles above its melting point.
    // A sunlike surface temperature, which at this radius and distance puts
    // about 1.4 kW/m^2 on the patch and settles it near 280 K — above the
    // 241.9 K this chemistry derives for water's melting point.
    let bright = 5800.0;
    sun(&mut w, bright);
    for _ in 0..120 {
        w.step_frame(50_000.0);
    }
    let summer = (temp(&w), liquid(&w), solid(&w), w.environment_at(patch).water);

    // Winter: the star dims. Nothing else is touched.
    // Cool enough that almost nothing reaches the patch.
    sun(&mut w, 1500.0);
    for _ in 0..120 {
        w.step_frame(50_000.0);
    }
    let winter = (temp(&w), liquid(&w), solid(&w), w.environment_at(patch).water);

    // Spring: it brightens again.
    sun(&mut w, bright);
    for _ in 0..240 {
        w.step_frame(50_000.0);
    }
    let spring = (temp(&w), liquid(&w), solid(&w), w.environment_at(patch).water);

    for (label, s) in [("summer", summer), ("winter", winter), ("spring", spring)] {
        println!(
            "  {label:<7} {:>7.1} K   liquid {:.3}   solid {:.3}   water available {:.3}",
            s.0, s.1, s.2, s.3
        );
    }

    assert!(winter.0 < summer.0, "the patch did not cool when the light went");
    assert!(winter.2 > summer.2, "nothing froze");
    assert!(winter.3 < summer.3, "the frozen patch measures as much water as the thawed one");
    assert!(spring.0 > winter.0, "the patch did not warm when the light returned");
    assert!(spring.1 > winter.1, "the thaw melted nothing");
    assert!(spring.3 > winter.3, "the thawed patch is no wetter than the frozen one");
}
