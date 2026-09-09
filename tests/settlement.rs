//! A small moon, a town on it, and a few buildings in the town.
//!
//! The point of the case is the *ladder*, not any one level of it. Each level
//! is a morphology program, and the bodies one program emits are the nodes the
//! next program runs in: terrain emits ground columns, a settlement emits
//! plots, a plot promoted becomes a tower that knows about floors. No level
//! knows about the ones above or below it, and nothing anywhere is a generator
//! that knows what a town is made of all the way down.

use phys::engine::World;
use phys::ids::NodeIdx;
use phys::math::{v3, Vec3};
use phys::morph::{Environment, Program};
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Matter};
use phys::tree::Tree;
use phys::units::*;

/// A moon 400 km across: small enough to be a moon, big enough to stand on.
const MOON_RADIUS: f64 = 2.0e5;
const MOON_DENSITY: f64 = 2400.0;

/// The patch of ground the town sits on, and the town itself.
const PATCH_SIDE: f64 = 2.0e3;
/// A small town: at the settlement's areal density this is about 300 m across,
/// which is a few streets rather than a city.
const TOWN_MASS: f64 = 6.3e7;

fn moon() -> Tree {
    let volume = 4.0 / 3.0 * std::f64::consts::PI * MOON_RADIUS.powi(3);
    let matter = Matter::neutral(
        volume * MOON_DENSITY,
        MOON_RADIUS,
        220.0,
        Program::Terrain.substrate(),
    );
    Tree::new(
        0x70_04_11,
        matter,
        Tier::containing(MOON_RADIUS),
        SampleSpec {
            count: 512,
            profile: Profile::Uniform,
            spectrum: MassSpectrum::Equal,
            kind: BodyKind::Grain,
            composition_scatter: 0.0,
            turbulent_fraction: 0.0,
        },
    )
}

/// Build moon → ground → town → building and return the chain.
fn a_moon_with_a_town(w: &mut World) -> (NodeIdx, NodeIdx, NodeIdx, NodeIdx) {
    let moon = w.tree.root;

    // A patch of that moon's surface, promoted out of it and given terrain.
    w.tree.refine(moon);
    let ground = w.tree.promote(moon, 0, w.tree.nodes[moon.get()].spec);
    assert!(!ground.is_none(), "the moon should yield a surface patch");
    let patch_mass = PATCH_SIDE * PATCH_SIDE * (PATCH_SIDE * 0.125) * MOON_DENSITY;
    w.emplace(ground, Program::Terrain, patch_mass, Environment::default());

    // The town, standing on that patch.
    w.tree.refine(ground);
    let mut town_spec = w.tree.nodes[ground.get()].spec;
    // A few buildings, which is what the case asks for. The count is the
    // resolution the town is drawn at, so this is "show me the plots", not
    // "there are only this many".
    town_spec.count = 24;
    let town = w.tree.promote(ground, 0, town_spec);
    assert!(!town.is_none(), "the ground should yield a plot for the town");
    w.emplace(town, Program::Settlement, TOWN_MASS, Environment::default());

    // And one building out of the town.
    w.tree.refine(town);
    let building = w.tree.promote(town, 0, w.tree.nodes[town.get()].spec);
    assert!(!building.is_none(), "the town should yield a plot");
    let plot_mass = w.tree.nodes[building.get()].matter.mass;
    w.emplace(building, Program::Tower, plot_mass, Environment::default());

    (moon, ground, town, building)
}

/// The ladder builds, and each level is the kind of thing it should be.
#[test]
fn a_moon_carries_a_town_carries_buildings() {
    let mut w = World::new(moon(), 20.0);
    let (moon, ground, town, building) = a_moon_with_a_town(&mut w);

    for (label, n, program) in [
        ("ground", ground, Program::Terrain),
        ("town", town, Program::Settlement),
        ("building", building, Program::Tower),
    ] {
        let (tier, radius, built, got) = {
            let node = &w.tree.nodes[n.get()];
            let m = node.morphology.as_ref().expect("should carry a program");
            (node.tier, node.matter.radius, m.built, m.program)
        };
        assert_eq!(got, program, "{label} has the wrong program");
        let parts = w.tree.refine(n).len();
        println!(
            "  {label:<9} {:>10} tier   radius {radius:>10.3e} m   {parts:>4} parts   {built:.3e} kg",
            format!("{tier:?}")
        );
        assert!(parts > 0, "{label} materialised into nothing");
        assert!(radius.is_finite() && radius > 0.0);
    }

    // Each level is smaller than the one holding it. Not a tautology: the
    // radius comes from the program's own `extent`, and a program that got its
    // arithmetic wrong would produce a town larger than its moon.
    let r = |n: NodeIdx| w.tree.nodes[n.get()].matter.radius;
    assert!(r(moon) > r(ground), "the ground is not inside the moon");
    assert!(r(ground) > r(town), "the town is not inside its patch");
    assert!(r(town) > r(building), "the building is not inside the town");
}

/// The ground is a surface: its parts vary in height, and they all reach down
/// to the same bedrock.
///
/// A terrain program that returned a constant height would still produce parts,
/// still conserve mass, and still pass every other test here — and would be a
/// car park rather than a moon.
#[test]
fn terrain_has_relief() {
    let mut w = World::new(moon(), 20.0);
    let (_, ground, _, _) = a_moon_with_a_town(&mut w);

    let bodies = w.tree.refine(ground).to_vec();
    let topo = w.tree.nodes[ground.get()].topology.clone().expect("terrain has a topology");
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..topo.tip.len() {
        lo = lo.min(topo.tip[i].z);
        hi = hi.max(topo.tip[i].z);
    }
    let radius = w.tree.nodes[ground.get()].matter.radius;
    println!(
        "  {} columns, surface from {:.2} m to {:.2} m — relief {:.2} m over a {:.0} m patch",
        bodies.len(),
        lo,
        hi,
        hi - lo,
        radius * 2.0
    );
    assert!(bodies.len() >= 4, "a patch of ground is not four columns");
    assert!(
        hi - lo > radius * 1e-3,
        "the ground is perfectly flat: relief {:.3e} m across {radius:.3e} m",
        hi - lo
    );

    // Every column is anchored. Terrain is what load paths terminate in, so a
    // column resting on another column would mean the ground had a ceiling.
    let anchored = topo.support.iter().filter(|s| **s == phys::morph::NO_SUPPORT).count();
    assert_eq!(anchored, topo.support.len(), "some ground is not anchored to bedrock");
}

/// A half-built town is a core with its edges missing, not a scatter.
#[test]
fn a_town_fills_in_from_its_centre() {
    let mut w = World::new(moon(), 20.0);
    let (_, _, town, _) = a_moon_with_a_town(&mut w);

    let full = w.tree.refine(town).len();
    let spread = |w: &mut World| {
        let bodies = w.tree.refine(town).to_vec();
        bodies.iter().map(|b| b.pos.norm()).fold(0.0f64, f64::max)
    };
    let spread_full = spread(&mut w);

    // Wind it back to a third built and regenerate.
    w.tree.nodes[town.get()].morphology.as_mut().unwrap().progress = 0.33;
    w.tree.nodes[town.get()].bodies.clear();
    let third = w.tree.refine(town).len();
    let spread_third = spread(&mut w);

    println!(
        "  fully built: {full} plots out to {spread_full:.1} m; a third built: {third} plots out to {spread_third:.1} m"
    );
    assert!(third < full, "a third-built town has as many plots as a finished one");
    assert!(
        spread_third < spread_full,
        "the partly built town reaches as far as the finished one, so it is filling in at random rather than from the centre"
    );
}

/// Promoting a plot gives a building that stands on the town's ground.
#[test]
fn a_plot_becomes_a_building_that_knows_about_floors() {
    let mut w = World::new(moon(), 20.0);
    let (_, _, town, building) = a_moon_with_a_town(&mut w);

    let bodies = w.tree.refine(building).to_vec();
    let topo = w.tree.nodes[building.get()].topology.clone().expect("a tower has a topology");
    // A tower is columns and beams stacked, so unlike terrain most of its parts
    // rest on another part rather than on the ground.
    let anchored = topo.support.iter().filter(|s| **s == phys::morph::NO_SUPPORT).count();
    println!(
        "  {} parts, {} of them founded on the ground, {} carried by the frame",
        bodies.len(),
        anchored,
        topo.support.len() - anchored
    );
    assert!(bodies.len() > 4, "a building is more than four parts");
    assert!(
        anchored < topo.support.len(),
        "every part of the building is on the ground, so it has no floors above the first"
    );
    // And it is a child of the town, so it moves when the town does.
    assert_eq!(w.tree.nodes[building.get()].parent, town);
}

/// Building the whole ladder neither creates nor destroys anything.
#[test]
fn generating_a_world_conserves_it() {
    let mut w = World::new(moon(), 20.0);
    let before = w.conserved();
    let (_, ground, town, building) = a_moon_with_a_town(&mut w);
    for n in [ground, town, building] {
        w.tree.refine(n);
    }
    let after = w.conserved();
    let rel = (before.baryon - after.baryon).abs() / before.baryon.abs().max(1e-300);
    println!(
        "  baryon {:.6e} -> {:.6e}, relative change {rel:.3e}",
        before.baryon, after.baryon
    );
    assert!(rel < 1e-9, "generating the world changed its baryon number by {rel:.3e}");
    let _ = v3(0.0, 0.0, 0.0);
    let _: Vec3 = Vec3::ZERO;
}
