//! A node's tier follows its size.
//!
//! A tier is a physics regime, and `Tier::containing` derives it from a radius.
//! Until `docs/PLAY.md` Phase 1's "tier revisited on size change", that
//! derivation ran exactly once — at `promote` — and every later thing that
//! changed a node's size left it behind. `plant`, `emplace`, growth and damage
//! all set `matter.radius` from the program's extent, so a node could sit two
//! tiers from what its own radius said.
//!
//! `BACKLOG.md` measured the visible case (a terrain patch 1.4 km across still
//! marked `Planetary`) and called it benign. It was not: the *spec* travels with
//! the tier, and a node that keeps a coarse spec keeps a coarse **sampler**. A
//! tree seed planted in a galactic node kept `kind: Star, profile: Plummer,
//! spectrum: Kroupa` — so materialising a 1.45 m seedling would have sampled
//! stars inside it.

use phys::engine::{default_spec, galaxy, World};
use phys::morph::Program;
use phys::units::Tier;

/// A galaxy refined once, with slot 0 promoted — a node whose tier came from a
/// body 10^20 m across.
fn a_coarse_node() -> (World, phys::ids::NodeIdx) {
    let mut w = World::new(galaxy(0xB0A7, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 6;
    let root = w.tree.root;
    w.tree.refine(root);
    let child = w.tree.promote(root, 0, default_spec(Tier::Stellar));
    assert!(!child.is_none());
    w.tree.nodes[child.get()].spec.count = 3;
    w.tree.refine(child);
    (w, child)
}

/// Planting a tree in it makes it a tree-sized node, and the tier says so.
#[test]
fn planting_a_seed_retiers_the_node_that_holds_it() {
    let (mut w, child) = a_coarse_node();
    let before = w.tree.nodes[child.get()].tier;
    let before_radius = w.tree.nodes[child.get()].matter.radius;
    assert_eq!(before, Tier::Galactic, "the setup should start coarse");

    w.plant(child, Program::Tree, Some(Default::default()));

    let n = &w.tree.nodes[child.get()];
    assert!(
        n.matter.radius < before_radius * 1e-9,
        "planting should have made it a structure's size: {:e} -> {:e}",
        before_radius,
        n.matter.radius
    );
    assert_eq!(
        n.tier,
        Tier::containing(n.matter.radius),
        "a {:e} m node is {:?} by the size table and is marked {:?}",
        n.matter.radius,
        Tier::containing(n.matter.radius),
        n.tier
    );
    assert_eq!(w.tree.stats.retiers, 1, "the move should have been counted");
}

/// And the refinement policy moves with it, which is the half that bites.
///
/// A tier that changes without its spec is the failure `Tree::promote`
/// documents at length: the node materialises under a policy meant for a
/// different scale. Arriving through `plant` rather than through `promote` made
/// it no less wrong — a seedling whose sampler draws from a stellar initial
/// mass function produces stars.
#[test]
fn the_refinement_policy_moves_with_the_tier() {
    let (mut w, child) = a_coarse_node();
    let before = w.tree.nodes[child.get()].spec;
    assert_eq!(
        before.kind,
        phys::state::BodyKind::Star,
        "the setup should start with a stellar policy, and has {:?}",
        before.kind
    );

    w.plant(child, Program::Tree, Some(Default::default()));

    let after = w.tree.nodes[child.get()].spec;
    assert_ne!(
        after.kind,
        phys::state::BodyKind::Star,
        "a 1.4 m tree seedling still samples stars"
    );
    assert!(
        phys::sampler::tier_of(after.kind) >= w.tree.nodes[child.get()].tier,
        "the policy is still coarser than the node: {:?} for a {:?} node",
        after.kind,
        w.tree.nodes[child.get()].tier
    );
}

/// `emplace` is the other way a node becomes a structure, and takes the same
/// path. This is `BACKLOG.md`'s own case: a patch that stayed `Planetary` while
/// being a kilometre across.
#[test]
fn emplacing_a_structure_retiers_the_node() {
    let (mut w, child) = a_coarse_node();
    w.emplace(child, Program::Terrain, 1.0e9, None);
    let n = &w.tree.nodes[child.get()];
    assert_eq!(
        n.tier,
        Tier::containing(n.matter.radius),
        "a {:e} m terrain patch is marked {:?}",
        n.matter.radius,
        n.tier
    );
}

/// A node is inside its parent, so it cannot be coarser than one.
///
/// Structural rather than physical: a rounding that said otherwise would put a
/// child on a coarser solver than the thing containing it. The clamp has to be
/// made to *bind* to be under test — a node whose own radius already lands at
/// or below its parent's tier would pass with the clamp deleted, which is how
/// the first version of this passed with exactly that mutation applied. So the
/// node is authored to a size its parent's band does not contain.
#[test]
fn a_child_is_never_coarser_than_its_parent() {
    use phys::observe::{Interaction, Property};

    let (mut w, child) = a_coarse_node();
    // Make the parent fine, so that a large child would be coarser than it.
    w.plant(child, Program::Tree, Some(Default::default()));
    assert_eq!(w.tree.nodes[child.get()].tier, Tier::Continuum);
    w.tree.refine(child);
    let tier = w.tree.nodes[child.get()].tier;
    let grandchild = w.tree.promote(child, 0, default_spec(tier));
    assert!(!grandchild.is_none(), "the planted node should have parts to promote");

    // A size squarely in the planetary band, inside a `Continuum` parent.
    w.interact(Interaction::Author {
        target: grandchild,
        property: Property::Radius,
        value: 1.0e6,
    });
    let n = &w.tree.nodes[grandchild.get()];
    assert_eq!(
        Tier::containing(n.matter.radius),
        Tier::Planetary,
        "the setup should ask for a tier coarser than the parent's, and asks for {:?}",
        Tier::containing(n.matter.radius)
    );
    assert_eq!(
        n.tier,
        Tier::Continuum,
        "a {:e} m node inside a Continuum parent came out {:?}",
        n.matter.radius,
        n.tier
    );

    for n in w.tree.nodes.iter() {
        if !n.alive || n.parent.is_none() {
            continue;
        }
        let parent = w.tree.nodes[n.parent.get()].tier;
        assert!(n.tier >= parent, "a {:?} node sits inside a {:?} one", n.tier, parent);
    }
}

/// Coarsening does not retier, and that is deliberate.
///
/// There the radius is being *restored* rather than changed, and a node's tier,
/// spec and identity are preserved across a refine/coarsen round trip on
/// purpose — it is what makes leaving and coming back idempotent.
///
/// To be under test at all the node has to be one a retier *would* move: a node
/// whose tier already agrees with its radius round-trips unchanged whether or
/// not coarsen retiers, which is how the first version of this passed with a
/// `retier` call added to `coarsen`. So the node is given a tier its own radius
/// disagrees with, and the round trip has to preserve the disagreement.
///
/// Note which path this exercises: an undisturbed node coarsens through
/// `coarsen`'s *idempotent* branch — the one that keeps the coarse state as the
/// authority when the detail turned out to say nothing new, and the one that
/// makes leaving and coming back bit-for-bit identical rather than merely
/// accurate. That is the branch a retier would spoil, and it is a different
/// early return from the one a disturbed node takes.
#[test]
fn a_round_trip_through_coarsen_leaves_the_tier_alone() {
    let mut w = World::new(galaxy(0xC0A25E, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 64;
    let root = w.tree.root;
    w.tree.refine(root);
    let child = w.tree.promote(root, 0, default_spec(Tier::Stellar));
    assert!(!child.is_none());
    w.tree.nodes[child.get()].spec.count = 8;
    w.tree.refine(child);

    // A tier its radius does not agree with. Nothing legitimate sets one by
    // hand; this stands in for a node that has drifted, which is the state a
    // round trip must not quietly tidy up.
    w.tree.nodes[child.get()].tier = Tier::Molecular;
    let spec = w.tree.nodes[child.get()].spec;
    let retiers = w.tree.stats.retiers;
    assert_ne!(
        Tier::containing(w.tree.nodes[child.get()].matter.radius),
        Tier::Molecular,
        "the setup should be inconsistent, and is not"
    );

    w.tree.coarsen(child);
    w.tree.refine(child);

    assert_eq!(
        w.tree.nodes[child.get()].tier,
        Tier::Molecular,
        "a round trip moved the tier, and coarsening is where a node's tier, \
         spec and identity are preserved on purpose"
    );
    assert_eq!(
        w.tree.nodes[child.get()].spec.kind,
        spec.kind,
        "a round trip moved the refinement policy"
    );
    assert_eq!(
        w.tree.stats.retiers, retiers,
        "coarsening retiered something, and it must not"
    );
}

/// Growth changes a node's size every step, and the tier must not chatter.
///
/// A tier is a physics regime: changing one changes the node's solver, its
/// timestep and its cadence, so a node that flickered across a boundary would
/// take its solver with it. Measured over 200 frames of a growing tree: it
/// moves once, when the seed is planted, and never again.
#[test]
fn growth_does_not_make_the_tier_chatter() {
    let (mut w, child) = a_coarse_node();
    w.plant(child, Program::Tree, Some(Default::default()));
    let settled = w.tree.stats.retiers;
    let tier = w.tree.nodes[child.get()].tier;

    for _ in 0..200 {
        w.step_frame(50_000.0);
    }

    assert_eq!(
        w.tree.stats.retiers, settled,
        "growing moved a tier {} times after the plant settled it",
        w.tree.stats.retiers - settled
    );
    assert_eq!(w.tree.nodes[child.get()].tier, tier);
}
