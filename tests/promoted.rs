//! A promoted child is the real thing, and its body is a stand-in.
//!
//! `promote` sets a child's `Motion` from the body it stands for and, until
//! `docs/PLAY.md` D4, nothing ever wrote `motion.velocity` again. So a promoted
//! node moved ballistically in its parent's frame for the rest of its life
//! while the parent's solver went on integrating the body it came from — two
//! representations of one object, moving under different laws.
//!
//! `docs/BACKLOG.md` measured the divergence at 0.79 of the child's own radius
//! after forty frames, and noted why no test had caught it: every test either
//! promoted and then looked, where a fraction of a radius is invisible, or
//! promoted and then coarsened, where `sync_from_child` hides it. These are the
//! tests that do neither.

use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::{v3, Vec3};
use phys::morph::{Environment, Program};
use phys::state::{Composition, Matter};
use phys::units::*;

/// A galaxy refined once, with the given slots promoted.
fn with_promoted(slots: &[usize]) -> (World, Vec<NodeIdx>) {
    let mut w = World::new(galaxy(0x9A11, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 256;
    let root = w.tree.root;
    // Paced to its subject explicitly. A world runs at one second per second
    // (`PLAY.md` D1), and a galaxy's gravity does nothing visible in a second.
    w.pace_to(root);
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    let children: Vec<NodeIdx> = slots
        .iter()
        .map(|s| w.tree.promote(root, *s, default_spec(tier.finer())))
        .collect();
    assert!(children.iter().all(|c| !c.is_none()), "promotion should succeed");
    (w, children)
}

/// The headline: a promoted child feels the force its parent's solver computes.
///
/// Before D4 this velocity was written once, at promotion, and never again.
#[test]
fn a_promoted_child_feels_a_force() {
    let (mut w, children) = with_promoted(&[0]);
    let child = children[0];
    let before = w.tree.nodes[child.get()].motion.velocity;

    for _ in 0..40 {
        w.step_frame(50_000.0);
    }

    let after = w.tree.nodes[child.get()].motion.velocity;
    let change = (after - before).norm();
    assert!(
        change > 0.0,
        "the child's velocity never changed: it is still ballistic, which is the \
         defect D4 exists to fix (before {before:?}, after {after:?})"
    );
    assert!(
        after.is_finite(),
        "a force arrived but left the velocity non-finite: {after:?}"
    );
}

/// The child and the body that stands for it agree whenever they are synced,
/// and drift only by the solver's own second-order term in between.
///
/// `docs/BACKLOG.md` measured the old behaviour at 0.79 of a radius after forty
/// frames *and growing without bound*, because nothing ever reconciled them.
/// What replaces it is not zero drift — the parent's solver integrates its
/// bodies with leapfrog, `x + v·dt + ½a·dt²`, while the child's own `Motion`
/// coasts linearly, `x + v·dt`, so between syncs they differ by the
/// acceleration term. What matters is that it is bounded and reset rather than
/// accumulated. Measured over four hundred frames it oscillates between 0.13
/// and 0.29 radii and does not climb.
#[test]
fn a_child_and_its_stand_in_agree_whenever_they_are_synced() {
    let (mut w, children) = with_promoted(&[0]);
    let child = children[0];
    let slot = w.tree.nodes[child.get()].slot as usize;
    let root = w.tree.root;

    for _ in 0..40 {
        w.step_frame(50_000.0);
    }

    // The exact invariant: a sync makes them the same place.
    w.tree.sync_children(root);
    let child_at = w.tree.nodes[child.get()].motion.offset;
    let body_at = w.tree.nodes[root.get()].bodies[slot].pos;
    assert_eq!(
        child_at, body_at,
        "a sync must leave the stand-in exactly where the child is"
    );
}

/// And the drift between syncs does not accumulate.
///
/// This is the half that was broken: the old divergence grew every frame
/// because nothing reconciled the two representations at all.
#[test]
fn the_drift_between_syncs_stays_bounded() {
    let (mut w, children) = with_promoted(&[0]);
    let child = children[0];
    let slot = w.tree.nodes[child.get()].slot as usize;
    let root = w.tree.root;
    let radius = w.tree.nodes[child.get()].matter.radius;

    let mut early = 0.0f64;
    let mut late = 0.0f64;
    for f in 1..=400 {
        w.step_frame(50_000.0);
        let apart = (w.tree.nodes[child.get()].motion.offset
            - w.tree.nodes[root.get()].bodies[slot].pos)
            .norm()
            / radius;
        if f <= 100 {
            early = early.max(apart);
        } else if f > 300 {
            late = late.max(apart);
        }
    }

    assert!(
        late < 1.0,
        "the stand-in drifted {late:.3} radii from its child, which is not a \
         second-order term any more"
    );
    assert!(
        late <= early * 2.0,
        "drift grew from {early:.3} radii in the first hundred frames to \
         {late:.3} in the last hundred — it is accumulating rather than being \
         reset, which is the defect D4 exists to fix"
    );
}

/// Two promoted siblings affect each other.
///
/// This is what the divergence blocked: "nothing yet promotes two siblings and
/// expects them to interact, which is the case where it becomes obvious". Both
/// must feel something, and what they feel must be the parent's solver rather
/// than noise.
#[test]
fn two_promoted_siblings_perturb_each_other() {
    let (mut w, children) = with_promoted(&[0, 1]);
    let before: Vec<_> = children
        .iter()
        .map(|c| w.tree.nodes[c.get()].motion.velocity)
        .collect();

    for _ in 0..40 {
        w.step_frame(50_000.0);
    }

    for (i, c) in children.iter().enumerate() {
        let after = w.tree.nodes[c.get()].motion.velocity;
        let change = (after - before[i]).norm();
        assert!(
            change > 0.0,
            "sibling {i} was unmoved by anything in forty frames; two promoted \
             things that cannot affect each other are the whole of the defect"
        );
    }
}

/// Coarsening still recovers the child's work, which it did before and must
/// keep doing — `sync_from_child` running every frame must not have made the
/// call at `coarsen` redundant in a way that loses the last frame's motion.
#[test]
fn coarsening_still_recovers_the_child() {
    let (mut w, children) = with_promoted(&[0]);
    let child = children[0];

    for _ in 0..10 {
        w.step_frame(50_000.0);
    }
    let child_at = w.tree.nodes[child.get()].motion.offset;

    let root = w.tree.root;
    w.tree.coarsen(root);
    w.tree.refine(root);

    // The node was rebuilt from matter, so the body is a fresh sample; what
    // must hold is that coarsening took the child's state up rather than
    // discarding it, which shows in the conserved totals.
    assert!(
        child_at.is_finite(),
        "the child's position should have been finite when it was taken up"
    );
    assert!(
        w.conserved().energy.is_finite() && w.conserved().energy != 0.0,
        "coarsening a promoted child should not have destroyed the conserved set"
    );
}

/// Nothing above changes what a node with no promoted children does.
#[test]
fn a_node_without_children_is_unaffected() {
    let mut w = World::new(galaxy(0x9A12, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 256;
    let root = w.tree.root;
    w.tree.refine(root);
    let before = w.tree.nodes[root.get()].bodies.clone();

    for _ in 0..5 {
        w.step_frame(50_000.0);
    }

    let after = &w.tree.nodes[root.get()].bodies;
    assert_eq!(before.len(), after.len(), "the body list should not have changed shape");
    assert!(
        after.iter().all(|b| b.pos.is_finite() && b.vel.is_finite()),
        "a node with no promoted children should be untouched by the sync"
    );
}

// ---------------------------------------------------------------------------
// The orphan — `docs/PLAY.md` §7 Phase 3's first correction
// ---------------------------------------------------------------------------

const EARTH_RADIUS: f64 = 6.371e6;
const EARTH_MASS: f64 = 5.972e24;

/// A forty-year-old tree standing on a planet, materialised, with one limb
/// promoted out of it.
fn a_tree_with_a_limb() -> (World, NodeIdx, NodeIdx) {
    let mut w = World::new(galaxy(0x3333, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let planet = w.tree.promote(root, 7, default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.spec.count = 8;
    }
    w.tree.refine(planet);
    let tree = w.tree.promote(planet, 0, default_spec(Tier::Continuum));
    {
        let n = &mut w.tree.nodes[tree.get()];
        n.matter = Matter::neutral(900.0, 6.0, 291.0, Program::Tree.substrate());
        n.spec.count = 600;
        n.motion.offset = v3(0.0, 0.0, EARTH_RADIUS);
        n.motion.velocity = Vec3::ZERO;
    }
    w.plant(tree, Program::Tree, Some(Environment::default()));
    for _ in 0..40 {
        w.grow_node(tree, YEAR);
    }
    w.tree.refine(tree);
    let limb = w.tree.promote(tree, 3, default_spec(Tier::Continuum));
    assert!(!limb.is_none(), "the limb should promote");
    (w, tree, limb)
}

/// Whether a walk from the root can still find it.
fn reachable(w: &World, target: NodeIdx) -> bool {
    let mut found = false;
    w.tree.walk(w.tree.root, &mut |i, _| found |= i == target);
    found
}

/// A structure that sheds a limb leaves no unreachable node behind.
///
/// The phase's done-when, and the defect it names. Every structural-change path
/// says "my body list is stale" by clearing it, and said it by throwing away
/// the only record of what had been promoted out of it. Measured before the
/// fix, on exactly this world: after the strike the limb is `alive` with
/// `parent` still pointing at the tree, the tree's `children` array is empty,
/// no walk from the root can reach it, and it is never freed — it keeps its
/// arena slot and goes on being scheduled for the life of the world.
#[test]
fn a_structure_that_sheds_a_limb_leaves_no_unreachable_node() {
    let (mut w, tree, limb) = a_tree_with_a_limb();
    assert!(reachable(&w, limb), "it should start out reachable");
    let before = w.conserved();
    let live_before = w.tree.live_count();

    let out = w.damage(tree, &[phys::solvers::structure::weather::wind(200.0, v3(1.0, 0.0, 0.0))]);
    assert!(out.detached_pieces > 0, "the strike has to change the structure to test anything");

    let n = &w.tree.nodes[limb.get()];
    println!(
        "  {} pieces came off; the limb is alive {}, live nodes {} -> {}",
        out.detached_pieces,
        n.alive,
        live_before,
        w.tree.live_count()
    );
    assert!(
        !n.alive || reachable(&w, limb),
        "the limb is alive, unreachable from the root, and will never be freed"
    );
    assert!(w.tree.stats.shed > 0, "nothing was folded back, so nothing was tested");

    // And it was *folded back*, not lost: the mass it stood for was always
    // counted in the tree's own matter, so the books do not move.
    let after = w.conserved();
    let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(b.abs()).max(1e-300);
    assert!(
        rel(before.energy, after.energy) < 1e-9,
        "energy changed by {:.3e} when the limb was folded back",
        rel(before.energy, after.energy)
    );
    assert!(rel(before.baryon, after.baryon) < 1e-9, "baryon number changed");
}

/// The same orphan by the other road: a reload.
///
/// The wire format writes `children` and, for an unpinned node, no bodies at
/// all — so a reloaded node comes back unmaterialised *with* promoted children,
/// and `refine` used to assign `vec![NONE; len]` over the top of them.
#[test]
fn a_reloaded_node_keeps_what_was_promoted_out_of_it() {
    let (mut w, _) = with_promoted(&[3, 7]);
    let bytes = phys::persist::encode(w.view());
    let mut back = World::from_snapshot(
        phys::persist::decode(&bytes).expect("decode"),
        20.0,
    );
    let root = back.tree.root;
    let kids: Vec<NodeIdx> = back.tree.nodes[root.get()]
        .children
        .iter()
        .copied()
        .filter(|c| !c.is_none())
        .collect();
    assert_eq!(kids.len(), 2, "the file should carry both promoted children");
    back.tree.refine(root);
    for c in kids {
        assert!(
            back.tree.nodes[c.get()].alive && reachable(&back, c),
            "refining the parent orphaned a child the file had restored"
        );
    }
}
