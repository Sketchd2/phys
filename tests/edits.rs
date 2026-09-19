//! A change is an edit to the recipe, not a pile of stored bodies.
//!
//! `docs/PLAY.md` D19. Before it there was one flag and one question: `pinned`
//! meant "this node's detail cannot be regenerated", `Tree::pin` set it on the
//! whole ancestry, and the scheduler refused to coarsen anything carrying it.
//! Three separate things were riding on that:
//!
//! * **whether the detail can be drawn again** — a genuine property of the
//!   node, and the only one `pinned` is actually about;
//! * **whether an ancestor's materialisation still matches what `sample` would
//!   produce** — true of every ancestor of an edited node, and *not* a reason
//!   to keep a body list, because the ancestor's own bodies are still a sample;
//! * **whether the node may be collapsed at all** — which the scheduler was
//!   answering from the mixing time, so a grown or built thing, which mixes in
//!   infinite time and correctly so, could never be collapsed.
//!
//! These tests hold the two claims that separating them buys: a broken
//! structure stays broken while costing a recipe rather than a body list, and
//! editing one leaf does not file its whole ancestry in the store.

use phys::engine::{galaxy, World};
use phys::math::v3;
use phys::morph::{Environment, Program};
use phys::observe::{Interaction, Property};
use phys::solvers::structure::weather;
use phys::state::{Body, Composition, Matter};
use phys::units::*;

/// A tree-sized node on the surface of an Earth. Same construction as
/// `tests/topology.rs`, and for the same reason: a structure is loaded by the
/// field it is actually in, so the scene has to contain something to fall
/// towards.
fn on_an_earth(seed: u64, mass: f64, radius: f64, count: usize) -> (World, phys::ids::NodeIdx) {
    const EARTH_MASS: f64 = 5.972e24;
    const EARTH_RADIUS: f64 = 6.371e6;
    let mut w = World::new(galaxy(seed, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let planet = w.tree.promote(root, 7, phys::engine::default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.spec.count = 8;
    }
    w.tree.refine(planet);
    let node = w.tree.promote(planet, 0, phys::engine::default_spec(Tier::Continuum));
    {
        let n = &mut w.tree.nodes[node.get()];
        n.matter = Matter::neutral(mass, radius, 291.0, Program::Tree.substrate());
        n.spec.count = count;
        n.motion.offset = v3(0.0, 0.0, EARTH_RADIUS);
    }
    (w, node)
}

/// A storm takes a limb off. The break is in the recipe, so the tree may be
/// collapsed to it and must never be redrawn without it.
///
/// The two halves pull in opposite directions and that is the whole point.
/// *Forgetting* — a fresh draw from the equilibrium ensemble — would mend the
/// break, so a damaged structure must refuse it forever. *Collapsing* —
/// throwing the bodies away and running the description again — reproduces the
/// break exactly, so a damaged structure must accept it as readily as anything
/// else. One flag could not say both, which is why there are now two.
#[test]
fn a_broken_tree_collapses_to_its_recipe_and_comes_back_broken() {
    let (mut w, node) = on_an_earth(0x0D19, 900.0, 6.0, 3000);
    // Emplaced rather than planted: a seedling under snow is not a measurement
    // of anything, and this test is about what happens to the *record* of a
    // break rather than about what breaks.
    w.emplace(node, Program::Tree, 850.0, Some(Environment::default()));

    let intact = w.tree.refine(node).iter().filter(|b| b.radius > 0.0).count();
    let out = w.damage(node, &[weather::snow(0.25, 450.0, 30.0)]);
    assert!(out.broken_joints > 0, "the storm did nothing, so there is no edit to test");

    // The break is expressible, so nothing was filed. A `pinned` here would be
    // the old behaviour: a body list in the store for a change the genome and
    // eight bytes of event log already describe.
    assert!(
        !w.tree.nodes[node.get()].pinned,
        "a break the recipe can express must not pin a body list"
    );
    assert!(w.tree.nodes[node.get()].contains_edit, "the break was not recorded as an edit");
    assert!(
        w.tree.nodes[node.get()].parent.is_none()
            || w.tree.nodes[w.tree.nodes[node.get()].parent.get()].contains_edit,
        "the ancestry does not know there is an edit below it"
    );

    // Never forgettable...
    assert!(
        !w.forgettable(node),
        "a broken tree redrawn from its ensemble would come back mended"
    );
    // ...and still collapsible.
    assert!(
        w.collapsible(node),
        "a structure that carries its own description must still be able to release its detail"
    );

    let broken = w.tree.refine(node).iter().filter(|b| b.radius > 0.0).count();
    let detail_bytes = w.tree.nodes[node.get()].bodies.len() * std::mem::size_of::<Body>();
    let recipe_bytes = w.tree.nodes[node.get()].morphology.as_ref().unwrap().state_bytes();

    let filed_before = w.tree.persisted.len();
    w.tree.coarsen(node);
    assert!(
        w.tree.nodes[node.get()].bodies.is_empty(),
        "the node kept its bodies through a coarsen"
    );
    assert_eq!(
        w.tree.persisted.len(),
        filed_before,
        "collapsing a described structure filed a body list it can regenerate"
    );

    // And it comes back with the same members missing, not the ones it was
    // generated with.
    let again = w.tree.refine(node).iter().filter(|b| b.radius > 0.0).count();
    println!(
        "  {intact} structural parts -> {broken} after {} joints broke; regenerated {again}",
        out.broken_joints
    );
    println!(
        "  detail {detail_bytes} B -> recipe {recipe_bytes} B ({:.0}x)",
        detail_bytes as f64 / recipe_bytes as f64
    );
    assert!(broken < intact, "the storm broke joints but the structure did not lose parts");
    assert_eq!(again, broken, "the regenerated tree is not the tree that was broken");
    assert!(
        recipe_bytes < detail_bytes,
        "the description is not smaller than the detail it replaces"
    );
}

/// Editing one leaf files one node's detail, not its whole ancestry's.
///
/// This is the measurement D19 was written for. `Tree::pin` walked to the root,
/// so authoring a temperature on one node marked every node above it
/// unregenerable; each of those then wrote its entire sampled body list into
/// the store the moment it collapsed, and refused to collapse again afterwards.
/// A planet does not become unregenerable because somebody felled a tree on it.
#[test]
fn editing_a_leaf_does_not_file_its_ancestry() {
    let mut w = World::new(galaxy(0xED17, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 64;
    let root = w.tree.root;
    let path = w.drill_to(root, Tier::Continuum.max_radius(), &phys::engine::default_spec);
    let leaf = *path.last().unwrap();

    // One authored property on one node, which is the smallest edit there is.
    w.interact(Interaction::Author { target: leaf, property: Property::Temperature, value: 350.0 });

    let pinned = w.tree.nodes.iter().filter(|n| n.alive && n.pinned).count();
    let edited = w.tree.nodes.iter().filter(|n| n.alive && n.contains_edit).count();
    let resident: usize = w.tree.nodes.iter().filter(|n| n.alive).map(|n| n.bodies.len()).sum();

    // Now release everything above the leaf, bottom-up. A pinned ancestor files
    // every body it holds; an ancestor that merely contains an edit files none.
    for &idx in path.iter().rev().skip(1) {
        w.tree.coarsen(idx);
    }
    let filed: usize = w.tree.persisted.values().map(|b| b.len()).sum();

    println!(
        "  one authored leaf in a {}-deep ladder: {pinned} pinned, {edited} containing an edit",
        path.len()
    );
    println!("  {resident} resident bodies above it, {filed} filed on release");

    assert_eq!(pinned, 1, "authoring one node pinned {pinned} nodes");
    assert!(edited >= path.len() - 1, "the ancestry was not marked as containing an edit");
    assert!(resident > 10_000, "the ladder is too small for the measurement to mean anything");
    assert_eq!(filed, 0, "an ancestor filed a body list the sampler regenerates for free");
}
