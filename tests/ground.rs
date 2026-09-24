//! Ground: a planet's surface, parameterised rather than typed.
//!
//! `docs/PLAY.md` D6 and Phase 4. A planetary surface is not a new kind of
//! thing — it is a node whose own gravity has rounded it, divided into patches
//! by a cubed sphere, and a patch is a node like any other. Its cells *are* its
//! children, so the surface tree is the scale tree with no adapter between
//! them, and nothing is generated until somebody approaches it.

use phys::chem::{Mixture, Phase};
use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::{v3, Vec3};
use phys::observe::Observer;
use phys::recipe::Recipe;
use phys::state::{Composition, Matter};
use phys::units::*;

const EARTH_MASS: f64 = 5.972e24;
const EARTH_RADIUS: f64 = 6.371e6;

/// A world with one Earth-sized body in it, said to be made of silicate.
///
/// The mixture is the scenario speaking, exactly as the composition already is:
/// a node nobody has said what it is made of has no fact of the matter about
/// its surface, and `World::assess_surface` gives it none.
fn an_earth(seed: u64, cells: usize) -> (World, NodeIdx) {
    let mut w = World::new(galaxy(seed, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let planet = w.tree.promote(root, 0, default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.spec.count = cells * cells + 1;
    }
    let silica = w
        .substances
        .intern(phys::material::substances::silica_arrangement())
        .expect("silica analyses");
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0);
    w.set_mixture(planet, mix);
    (w, planet)
}

/// Make a node a ball of rock of this radius, keeping what it is made of.
///
/// Overwriting `matter` wholesale would drop the mixture with it, and a node
/// with no mixture has not been said what it is made of — which
/// `assess_surface` answers by giving it no surface, correctly and for a reason
/// that has nothing to do with its size.
fn resize(w: &mut World, node: NodeIdx, radius: f64) {
    let mix = w.mixture_of(node);
    let m = 4.0 / 3.0 * std::f64::consts::PI * radius * radius * radius * 2600.0;
    let n = &mut w.tree.nodes[node.get()];
    n.matter = Matter::neutral(m, radius, 290.0, Composition::primordial());
    n.matter.mixture = mix;
}

/// Descend towards a direction from the planet's centre, one level at a time.
///
/// By the parameterisation rather than by searching for the nearest cell: a
/// face is nine thousand kilometres across, and the nearest cell *centre* to a
/// point above the middle of one is not the cell under that point. This is the
/// same question `World::approach` asks, asked by hand.
fn descend(w: &mut World, from: NodeIdx, towards: Vec3, target: f64) -> Vec<NodeIdx> {
    let mut path = vec![from];
    let mut here = from;
    for _ in 0..40 {
        let Some((cell, side)) = w.tree.nodes[here.get()].morphology.as_ref().and_then(|m| {
            let Some(Recipe::Tiled(t)) = m.recipe.as_ref() else { return None };
            Some((t.cell_of_direction(towards)?, t.side))
        }) else {
            break;
        };
        if side <= target {
            break;
        }
        w.tree.refine(here);
        if cell >= w.tree.nodes[here.get()].bodies.len() {
            break;
        }
        let spec = w.tree.nodes[here.get()].spec;
        let child = w.tree.promote(here, cell, spec);
        if child.is_none() {
            break;
        }
        path.push(child);
        here = child;
    }
    path
}

/// The side of the patch a node holds, metres.
fn side_of(w: &World, node: NodeIdx) -> f64 {
    match w.tree.nodes[node.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
        Some(Recipe::Tiled(t)) => t.side,
        _ => f64::INFINITY,
    }
}

/// A body its own gravity has rounded acquires a surface; one it has not does
/// not.
///
/// The gate is measured on both sides: the central pressure of a
/// self-gravitating sphere against the strength of what it is made of. Nothing
/// here knows what a planet is.
#[test]
fn a_body_its_own_gravity_has_rounded_has_a_surface() {
    let (mut w, planet) = an_earth(0x6204, 8);
    assert!(w.assess_surface(planet), "an Earth is round");
    assert!(
        matches!(
            w.tree.nodes[planet.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()),
            Some(Recipe::Tiled(_))
        ),
        "and what it got is a tiling"
    );
    assert!(!w.assess_surface(planet), "and it is written down once");

    // The same rock, small enough that its own weight is nothing beside its
    // strength. This is the potato radius and it is derived, not stated.
    let (mut w, rock) = an_earth(0x6205, 8);
    resize(&mut w, rock, 500.0);
    assert!(!w.assess_surface(rock), "a 500 m asteroid is a potato");

    // And somewhere between them it turns over. Measured rather than asserted
    // at a number: the transition is wherever the two curves cross.
    let mut smallest = f64::INFINITY;
    for km in [0.5f64, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0] {
        let (mut w, body) = an_earth(0x6206, 8);
        resize(&mut w, body, km * 1000.0);
        if w.assess_surface(body) {
            smallest = smallest.min(km);
        }
    }
    println!("  the smallest rock that rounds itself is {smallest:.1} km across");
    // Tens of kilometres rather than the two hundred a real potato radius sits
    // at, and the difference is not a fudge: the criterion here is the
    // *fracture* strength at the body's own size, which D14's size effect makes
    // small for something this big, where a real rock at depth is held by its
    // confined yield strength. The direction is the safe one — the engine
    // rounds things a little sooner than nature does — and it is one law rather
    // than a radius somebody chose.
    assert!(
        (1.0..=200.0).contains(&smallest),
        "the transition should be kilometres to hundreds of kilometres, and is {smallest}"
    );
}

/// A planet nobody has visited costs its matter and nothing else.
#[test]
fn a_planet_nobody_has_visited_holds_no_detail() {
    let (mut w, planet) = an_earth(0x6207, 8);
    w.assess_surface(planet);
    let n = &w.tree.nodes[planet.get()];
    assert!(n.bodies.is_empty(), "nothing is materialised until it is asked for");
    let bytes = n.detail_bytes();
    let state = n.morphology.as_ref().unwrap().state_bytes();
    println!("  a whole planet's surface: {state} B of rule, {bytes} B of detail");
    assert_eq!(bytes, 0);
    assert!(state < 512, "and the rule is a handful of numbers: {state} B");
}

/// From orbit to a square metre, and the mass is conserved every step.
///
/// The phase's first done-when. Each level's cells *are* its children: the
/// bodies a patch holds are the sub-patches it divides into, so descending is
/// promoting, and there is no second structure to keep in step.
#[test]
fn an_observer_descends_from_orbit_to_a_square_metre() {
    let (mut w, planet) = an_earth(0x6208, 8);
    assert!(w.assess_surface(planet));
    let path = descend(&mut w, planet, v3(1.0, 0.2, 0.1).unit(), 1.0);
    let deepest = *path.last().unwrap();
    let side = side_of(&w, deepest);
    println!(
        "  {} levels from {:.3e} m across to {:.4} m across",
        path.len() - 1,
        EARTH_RADIUS,
        side
    );
    assert!(side <= 1.0, "the descent arrived at a patch {side:.4} m across");
    assert!(path.len() < 20, "and it took {} levels", path.len() - 1);

    // Mass conserves at every level: a patch's cells and what is under them are
    // a partition of it.
    for &node in &path {
        let n = &w.tree.nodes[node.get()];
        if n.bodies.is_empty() {
            continue;
        }
        let held: f64 = n.bodies.iter().map(|b| b.mass).sum();
        assert!(
            (held - n.matter.mass).abs() <= 1e-12 * n.matter.mass,
            "a patch of {:.4e} kg holds {held:.4e}",
            n.matter.mass
        );
    }

    // And gravity arrives along the standing thing's own down, at the strength
    // the planet's mass and radius say. This is the `BACKLOG.md` entry whose
    // trigger was "the first terrain patch".
    let g = w.tree.gravity_at(deepest);
    let closed = G * EARTH_MASS / (EARTH_RADIUS * EARTH_RADIUS);
    println!(
        "  gravity there {:.4} m/s^2 against {closed:.4}, tilted {:.3} degrees off its own down",
        g.norm(),
        (g.x.hypot(g.y) / g.norm().max(1e-30)).asin().to_degrees()
    );
    assert!(g.z < 0.0, "down is down: {g:?}");
    assert!(
        g.x.hypot(g.y) < 0.02 * g.norm(),
        "and it is not sideways: {g:?}"
    );
    assert!(
        (g.norm() / closed - 1.0).abs() < 0.05,
        "and it is the field the planet's own mass makes: {:.4} against {closed:.4}",
        g.norm()
    );
}

/// The ground behind an observer regenerates bit-identically.
///
/// The invariant `tests/consistency.rs` guards, applied to a surface: a patch
/// is an address and a rule, so releasing it and asking again has to give back
/// the same ground rather than a plausible one.
#[test]
fn terrain_behind_an_observer_regenerates_bit_for_bit() {
    let (mut w, planet) = an_earth(0x6209, 8);
    assert!(w.assess_surface(planet));
    let path = descend(&mut w, planet, v3(0.3, 1.0, 0.2).unit(), 50.0);
    let patch = *path.last().unwrap();

    let before: Vec<phys::state::Body> = w.tree.refine(patch).to_vec();
    let recipe = w.tree.nodes[patch.get()].morphology.as_ref().unwrap().recipe.clone();
    assert!(!before.is_empty());

    // Look away: the detail is released, because a patch is a structure with a
    // recipe and `collapsible` says so.
    assert!(w.collapsible(patch), "a patch must be able to release its detail");
    w.tree.coarsen(patch);
    assert!(w.tree.nodes[patch.get()].bodies.is_empty());

    let after: Vec<phys::state::Body> = w.tree.refine(patch).to_vec();
    assert_eq!(after.len(), before.len(), "a different number of cells came back");
    let mut worst = 0.0f64;
    for (a, b) in before.iter().zip(&after) {
        worst = worst.max((a.pos - b.pos).norm());
        assert_eq!(a.pos, b.pos, "a cell moved");
        assert_eq!(a.mass, b.mass, "a cell changed mass");
        assert_eq!(a.half, b.half, "a cell changed shape");
        assert_eq!(a.orientation, b.orientation, "a cell turned");
    }
    assert_eq!(
        recipe,
        w.tree.nodes[patch.get()].morphology.as_ref().unwrap().recipe,
        "and the rule that draws it is the same rule"
    );
    println!("  {} cells regenerated, worst displacement {worst:.1e} m", after.len());
}

/// Ten kilometres across patch boundaries, and the ground is there all the way.
///
/// The done-when's "travels ten kilometres across patch boundaries". What it
/// measures is that a walk along the surface keeps finding ground at the same
/// height above the planet's centre — the patches meet, rather than leaving a
/// gap or a step at every boundary.
#[test]
fn a_walk_of_ten_kilometres_crosses_patch_boundaries() {
    let (mut w, planet) = an_earth(0x620A, 8);
    assert!(w.assess_surface(planet));

    // Ten kilometres is an angle on a planet this size.
    let start = v3(1.0, 0.0, 0.0);
    let across = v3(0.0, 1.0, 0.0);
    let total = 10_000.0;
    let steps = 40;
    let mut patches: Vec<NodeIdx> = Vec::new();
    let mut heights: Vec<f64> = Vec::new();
    for k in 0..=steps {
        let t = total * k as f64 / steps as f64 / EARTH_RADIUS;
        let dir = (start + across.scale(t)).unit();
        let path = descend(&mut w, planet, dir, 200.0);
        let patch = *path.last().unwrap();
        if patches.last() != Some(&patch) {
            patches.push(patch);
        }
        // How far this patch's own surface is from the planet's centre.
        let r = w.tree.offset_from(planet, patch, Vec3::ZERO).value.norm();
        heights.push(r + 0.5 * w.tree.nodes[patch.get()].matter.radius);
    }
    let distinct = {
        let mut u = patches.clone();
        u.sort();
        u.dedup();
        u.len()
    };
    let hi = heights.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let lo = heights.iter().cloned().fold(f64::INFINITY, f64::min);
    println!(
        "  {steps} steps over 10 km crossed into {distinct} patches; the ground stays \
         within {:.1} m of {:.4e} m from the centre",
        hi - lo,
        0.5 * (hi + lo)
    );
    assert!(distinct > 1, "ten kilometres should cross a patch boundary");
    assert!(
        hi - lo < 0.02 * (hi + lo) * 0.5,
        "the ground steps by {:.1} m between patches",
        hi - lo
    );
}

/// An observer approaching generates the surface, and nothing else does.
#[test]
fn nothing_is_generated_until_it_is_approached() {
    let (mut w, planet) = an_earth(0x620B, 8);
    assert!(w.assess_surface(planet));
    assert_eq!(w.tree.live_count(), 2, "a root and a planet");

    // An observer a metre above the ground, looking down.
    let obs = Observer {
        anchor: planet,
        offset: v3(EARTH_RADIUS + 1.0, 0.0, 0.0),
        look: v3(-1.0, 0.0, 0.0),
        angular_resolution: 1e-3,
        ..Default::default()
    };
    w.add_observer(obs);
    for _ in 0..40 {
        w.step_frame(20_000.0);
    }
    let live = w.tree.live_count();
    let approaches = w.stats.approaches;
    let deepest = w
        .tree
        .nodes
        .iter()
        .filter(|n| n.alive)
        .map(|n| n.matter.radius)
        .fold(f64::INFINITY, f64::min);
    println!(
        "  an observer a metre up drew {live} nodes out of a planet in {approaches} \
         approaches; the finest is {deepest:.3} m across"
    );
    assert!(approaches > 0, "the observer generated nothing");
    assert!(deepest < 1.0e5, "and it got somewhere: {deepest:.3e} m");
}
