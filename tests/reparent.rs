//! Moving a node between frames.
//!
//! Until this existed a node's place in the hierarchy was fixed at creation for
//! life, which is right for containment that does not change — a star does not
//! leave its cluster during play — and wrong for everything at the scale the
//! play space is scoped to, where a thing picked up enters your frame and a
//! thing thrown enters the ground's.

use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::Vec3;


fn a_world() -> World {
    let mut w = World::new(galaxy(0x9E77, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 400;
    w
}

/// Promote two siblings out of the same parent, and hand one to the other.
fn two_siblings(w: &mut World) -> (NodeIdx, NodeIdx, NodeIdx) {
    let root = w.tree.root;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    let a = w.tree.promote(root, 3, default_spec(tier.finer()));
    let b = w.tree.promote(root, 7, default_spec(tier.finer()));
    assert!(!a.is_none() && !b.is_none());
    (root, a, b)
}

/// Nothing is created or destroyed by moving it.
///
/// The one that would be easiest to get wrong and hardest to notice. A promoted
/// body is a *stand-in*: `sum_conserved` counts the child and skips the body
/// wherever a slot is promoted. Leave the vacated body behind and the old
/// parent silently keeps the mass that just left it, and the world gains a
/// duplicate of whatever was moved.
#[test]
fn a_move_conserves_everything() {
    let mut w = a_world();
    let (_, a, b) = two_siblings(&mut w);

    let before = w.conserved();
    assert!(w.reparent(a, b), "the move should be accepted");
    let after = w.conserved();

    let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(y.abs()).max(1e-300);
    println!(
        "  mass-energy {:.6e} -> {:.6e}   baryon {:.6e} -> {:.6e}",
        before.energy, after.energy, before.baryon, after.baryon
    );
    assert!(
        rel(before.energy, after.energy) < 1e-9,
        "energy changed by {:.3e}",
        rel(before.energy, after.energy)
    );
    assert!(
        rel(before.baryon, after.baryon) < 1e-9,
        "baryon number changed by {:.3e}",
        rel(before.baryon, after.baryon)
    );
    assert!(
        (before.charge - after.charge).abs() <= before.charge.abs() * 1e-9 + 1e-30,
        "charge changed"
    );
}

/// The moved node is where it was, seen from its new parent.
///
/// A `Motion` is relative to the parent, so a move that does not re-express it
/// teleports the object. The check is that the separation between the moved
/// node and a third node — measured through the tree, which is how everything
/// in the engine measures — is the same before and after.
#[test]
fn a_move_does_not_teleport() {
    let mut w = a_world();
    let (root, a, b) = two_siblings(&mut w);

    let before = w.tree.separation(root, Vec3::ZERO, a, Vec3::ZERO).value;
    assert!(w.reparent(a, b));
    let after = w.tree.separation(root, Vec3::ZERO, a, Vec3::ZERO).value;

    let drift = (after - before).norm();
    let scale = before.norm().max(1e-30);
    println!(
        "  from the root: {:.6e} m before, {:.6e} m after, drift {:.3e} m ({:.2e} relative)",
        before.norm(),
        after.norm(),
        drift,
        drift / scale
    );
    assert!(
        drift / scale < 1e-9,
        "moving it shifted it by {drift:.3e} m, {:.3e} of its distance",
        drift / scale
    );
    assert_eq!(w.tree.nodes[a.get()].parent, b, "it should now be under b");
}

/// A node's key is its path, so a move rekeys it and everything beneath it.
#[test]
fn a_move_rekeys_the_subtree() {
    let mut w = a_world();
    let (_, a, b) = two_siblings(&mut w);

    // Give `a` a grandchild, so there is a subtree to carry.
    let tier = w.tree.nodes[a.get()].tier;
    w.tree.refine(a);
    let grand = w.tree.promote(a, 1, default_spec(tier.finer()));
    assert!(!grand.is_none());

    let a_key = w.tree.nodes[a.get()].key;
    let g_key = w.tree.nodes[grand.get()].key;
    assert!(w.reparent(a, b));

    assert_ne!(w.tree.nodes[a.get()].key, a_key, "the moved node kept its key");
    assert_ne!(w.tree.nodes[grand.get()].key, g_key, "a descendant kept its key");
    // And the derivation still holds: a child's key comes from its parent's.
    let slot = w.tree.nodes[grand.get()].slot as u64;
    assert_eq!(
        w.tree.nodes[grand.get()].key,
        w.tree.nodes[a.get()].key.child(slot),
        "the descendant's key is no longer derivable from its parent's"
    );
    assert_eq!(w.tree.nodes[grand.get()].depth, w.tree.nodes[a.get()].depth + 1);
}

/// Chemistry travels with the object.
///
/// The migration is the reason `World::reparent` exists at all rather than
/// callers using `Tree::reparent` directly: a thrown beaker that arrives
/// without its contents is the symptom of forgetting one side table.
#[test]
fn a_move_carries_what_the_node_is_made_of() {
    use phys::chem::{Arrangement, Bond, Element, Lattice, Mixture, Order, Phase};

    let mut w = a_world();
    let (_, a, b) = two_siblings(&mut w);

    let salt = w
        .substances
        .intern(Arrangement::crystal(
            vec![Element(11), Element(17)],
            vec![Bond::new(0, 1, Order::Ionic)],
            Lattice::Cubic { a: 3.55e-10 },
        ))
        .expect("salt analyses");
    let mut mix = Mixture::new();
    mix.add(salt, Phase::Solid, 0.25);
    w.set_mixture(a, mix);
    assert!(!w.mixture_of(a).is_empty(), "precondition");

    assert!(w.reparent(a, b));
    assert!(
        !w.mixture_of(a).is_empty(),
        "the node arrived without its chemistry: a side table was not migrated"
    );
}

/// Identity is issued, and a move does not change it.
///
/// This is the property the whole of D2 exists for: what a thing *is* stops
/// being a function of where it happens to be. Before it, the node's address
/// was also its name, so moving it renamed it and every side table keyed by
/// the old name was left pointing at a stranger.
#[test]
fn a_moved_node_keeps_its_name_while_its_address_changes() {
    let mut w = a_world();
    let (_, a, b) = two_siblings(&mut w);

    let id = w.identify(a);
    let key_before = w.tree.nodes[a.get()].key;
    assert!(!id.is_none(), "an identified node has a name");

    assert!(w.reparent(a, b));

    let key_after = w.tree.nodes[a.get()].key;
    assert_ne!(key_before, key_after, "the address must have changed, or this proves nothing");
    assert_eq!(
        w.identity(a),
        Some(id),
        "the name must not have changed: identity is issued, not derived from position"
    );
}

/// An identity is never handed out twice, and asking does not create one.
#[test]
fn identities_are_issued_once_and_never_reused() {
    let mut w = a_world();
    let (_, a, b) = two_siblings(&mut w);

    assert_eq!(w.identity(a), None, "asking must not issue");
    let first = w.identify(a);
    assert_eq!(w.identify(a), first, "identifying twice gives the same name");

    let other = w.identify(b);
    assert_ne!(first, other, "two things cannot share a name");
    assert!(w.next_entity > other.0, "the counter moves past what it issued");
}

/// The side tables stay where they are, and what a node is made of moves with
/// the node.
///
/// `World::reparent` used to hold the only enumeration of the tables, so adding
/// one meant remembering to add a line there. Now only the address-to-identity
/// index moves, and this asserts the tables themselves are untouched — which is
/// what makes a *new* side table safe by default rather than safe if somebody
/// remembered.
///
/// **Speciation is no longer one of them.** `docs/PLAY.md` D17 puts the
/// `Mixture` on the node's own `Matter`, so it needs no protection from a move
/// at all: what a thing is made of travels with the thing, which is what the
/// side table was contriving to imitate by keying on a name. This test used to
/// read the table directly and assert the entry had not moved; it now asserts
/// the stronger thing, that the answer to "what is this made of" is unchanged
/// across a move that rekeys the node.
#[test]
fn a_move_does_not_touch_the_tables_keyed_by_identity() {
    use phys::chem::{Arrangement, Bond, Element, Lattice, Mixture, Order, Phase};

    let mut w = a_world();
    let (_, a, b) = two_siblings(&mut w);

    let salt = w
        .substances
        .intern(Arrangement::crystal(
            vec![Element(11), Element(17)],
            vec![Bond::new(0, 1, Order::Ionic)],
            Lattice::Cubic { a: 3.55e-10 },
        ))
        .expect("salt analyses");
    let mut mix = Mixture::new();
    mix.add(salt, Phase::Solid, 0.25);
    w.set_mixture(a, mix);

    // An environment *is* still keyed by identity, so it is what holds the
    // original claim up now that chemistry does not.
    let id = w.identify(a);
    w.environments
        .insert(id, phys::morph::Environment { light_flux: 340.0, ..Default::default() });
    let key_before = w.tree.nodes[a.get()].key;
    let before = w.mixture_of(a).fraction_of(salt);
    assert!(before > 0.0, "precondition: the node is a quarter salt");

    assert!(w.reparent(a, b));

    assert_ne!(
        w.tree.nodes[a.get()].key, key_before,
        "a move rekeys the node, or this test is measuring nothing"
    );
    assert_eq!(
        w.mixture_of(a).fraction_of(salt),
        before,
        "what a node is made of is part of its matter and a move cannot touch it"
    );
    assert_eq!(w.identity(a), Some(id), "and the node still answers to its name");
    assert!(
        w.environments.contains_key(&id),
        "the environment table is keyed by that name and should not have moved"
    );
}

/// Naming a node is deterministic, because only *events* name one.
///
/// A clock is bookkeeping the frame budget created, and which nodes the budget
/// advances depends on a wall-clock allowance — so if a clock named its node,
/// identity would depend on how fast the machine is. `next_entity` is
/// persisted, so that divergence would be durable, and `docs/PLAY.md` D10 makes
/// world state a function of the seed and the input log. An input naming an
/// entity would then bind to a different node on replay.
///
/// Measured before the fix: 50,000 us of budget gave `next_entity` 2 and one
/// checksum, 5,000 us and below gave 5 and another. Same length, different
/// world.
#[test]
fn identity_does_not_depend_on_how_fast_the_machine_is() {
    fn run(budget_us: f64) -> (u64, usize, Vec<phys::ids::PathKey>) {
        let mut w = a_world();
        w.tree.nodes[0].spec.count = 4096;
        let root = w.tree.root;
        // Paced to its subject, which is what makes this a test about the
        // *budget*: a world runs at one second per second (`PLAY.md` D1), and a
        // fixed span would hide the starvation this is looking for.
        w.pace_to(root);
        for idx in w.drill_to(root, 1.0e2, &default_spec) {
            w.tree.refine(idx);
        }
        for _ in 0..40 {
            w.step_frame(budget_us);
        }
        let mut advanced: Vec<_> = w.clocks.keys().copied().collect();
        advanced.sort();
        (w.next_entity, w.identities.len(), advanced)
    }

    let generous = run(50_000.0);
    for budget in [5_000.0, 200.0, 1.0] {
        let starved = run(budget);
        assert_eq!(
            starved.0, generous.0,
            "a {budget} us budget issued a different number of names than a generous one"
        );
        assert_eq!(starved.1, generous.1, "and a different index");
    }

    // The control: the budget must actually have changed what the engine did,
    // or this test is asserting nothing. Clocks are keyed by address precisely
    // so that they may vary with the budget without naming anything. *Which*
    // nodes, not how many: measured at Phase 5, a generous budget and a
    // starved one each advanced two nodes over forty frames, and not the same
    // two.
    let starved = run(1.0);
    assert_ne!(
        starved.2, generous.2,
        "the budget did not change which nodes were advanced, so this proves nothing"
    );
}

/// A move that is not a move is refused, and a cycle above all.
///
/// A cycle would not merely be wrong. `lca`, `offset_from` and `disturb` all
/// walk parents with `while !cur.is_none()`, so a node that is its own ancestor
/// hangs the engine rather than producing a wrong answer.
#[test]
fn impossible_moves_are_refused() {
    let mut w = a_world();
    let (root, a, b) = two_siblings(&mut w);
    let tier = w.tree.nodes[a.get()].tier;
    w.tree.refine(a);
    let grand = w.tree.promote(a, 2, default_spec(tier.finer()));

    assert!(!w.reparent(root, a), "the root has no outside to move to");
    assert!(!w.reparent(a, a), "a node cannot contain itself");
    assert!(!w.reparent(a, grand), "a node cannot move into its own descendant");
    assert!(!w.reparent(a, root), "a is already under the root");
    assert!(!w.reparent(NodeIdx::NONE, b), "nothing is not a node");

    // And after all that refusal the tree is still walkable, which is the
    // property a cycle would have destroyed.
    assert_eq!(w.tree.lca(grand, b), root);
    let _ = w.tree.separation(grand, Vec3::ZERO, b, Vec3::ZERO);
}
