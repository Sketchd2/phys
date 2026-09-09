//! Something changes hands.
//!
//! ```sh
//! cargo run --release --bin phys-rehome
//! ```
//!
//! Until re-parenting existed, a node's place in the tree was fixed when it was
//! created and never changed again. That is right for the containment the
//! engine was built on — a star does not leave its cluster during play, and a
//! nucleus does not leave its atom — and wrong for everything at the scale the
//! play space is scoped to, where containment changes constantly: pick a thing
//! up and it enters your frame, throw it and it enters the ground's, walk into
//! a ship and you enter the ship's.
//!
//! This walks one object through three frames and checks, at every step, the
//! three things a move must not break. It is meant to be read as well as run:
//! each section prints what it did and the number that says whether it worked.

use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::Vec3;

fn rule(title: &str) {
    println!("\n\x1b[1m{title}\x1b[0m");
    println!("{}", "─".repeat(title.chars().count()));
}

fn si(v: f64, unit: &str) -> String {
    if !v.is_finite() {
        return format!("∞ {unit}");
    }
    format!("{v:.4e} {unit}")
}

/// Where a node is, measured from the root the way everything else in the
/// engine measures: through the tree, to the common ancestor and back down.
fn world_position(w: &World, n: NodeIdx) -> Vec3 {
    w.tree.separation(w.tree.root, Vec3::ZERO, n, Vec3::ZERO).value
}

fn describe(w: &World, label: &str, n: NodeIdx) {
    let node = &w.tree.nodes[n.get()];
    println!(
        "  {label:<10} parent {:>3}   depth {}   key {}   {} from the root",
        node.parent.get(),
        node.depth,
        node.key,
        si(world_position(w, n).norm(), "m"),
    );
}

fn main() {
    let mut w = World::new(galaxy(0x9E77, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 400;
    let root = w.tree.root;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;

    // Three siblings out of the same parent. `carrier` and `destination` are
    // the two frames; `cargo` is the thing that moves between them.
    let cargo = w.tree.promote(root, 3, default_spec(tier.finer()));
    let carrier = w.tree.promote(root, 7, default_spec(tier.finer()));
    let destination = w.tree.promote(root, 11, default_spec(tier.finer()));

    rule("Three objects, all children of the root");
    describe(&w, "cargo", cargo);
    describe(&w, "carrier", carrier);
    describe(&w, "destination", destination);

    // The invariants, sampled before anything moves.
    let conserved0 = w.conserved();
    let where0 = world_position(&w, cargo);

    rule("Picked up: cargo moves into the carrier's frame");
    let ok = w.reparent(cargo, carrier);
    println!("  accepted: {ok}");
    describe(&w, "cargo", cargo);
    report(&w, conserved0, where0, cargo);

    rule("Handed over: cargo moves again, carrier to destination");
    let ok = w.reparent(cargo, destination);
    println!("  accepted: {ok}");
    describe(&w, "cargo", cargo);
    report(&w, conserved0, where0, cargo);

    rule("Put down: cargo returns to the root");
    let ok = w.reparent(cargo, root);
    println!("  accepted: {ok}");
    describe(&w, "cargo", cargo);
    report(&w, conserved0, where0, cargo);

    rule("Moves that are refused, and why they must be");
    println!("  the root, which has no outside to move to ...... {}", w.reparent(root, cargo));
    println!("  a node into itself ............................. {}", w.reparent(cargo, cargo));
    println!("  a node into its own descendant ................. {}", {
        let t = w.tree.nodes[cargo.get()].tier;
        w.tree.refine(cargo);
        let inner = w.tree.promote(cargo, 1, default_spec(t.finer()));
        w.reparent(cargo, inner)
    });
    println!(
        "\n  A cycle is the one that matters. `lca`, `offset_from` and `disturb`\n  \
         all walk parents with `while !cur.is_none()`, so a node that became its\n  \
         own ancestor would not give a wrong answer — it would hang. After all\n  \
         three refusals the tree is still walkable:"
    );
    println!("    lca(cargo, carrier) = {:?}", w.tree.lca(cargo, carrier).get());

    rule("And the world still runs");
    for _ in 0..20 {
        w.step_frame(50_000.0);
    }
    let after = w.conserved();
    println!(
        "  20 frames later: {} live nodes, {} moves recorded, energy drift {:.3e}",
        w.tree.live_count(),
        w.tree.stats.reparents,
        rel(conserved0.energy, after.energy),
    );
    println!(
        "\n  Re-parenting is what turns the tree from a record of what owns what\n  \
         into an index of what is where. Everything at play scale needs it: a\n  \
         thrown branch cannot leave the tree it broke off without it, and two\n  \
         objects cannot change hands."
    );
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(1e-300)
}

/// The three things a move must not break, checked after every one.
fn report(w: &World, conserved0: phys::state::Conserved, where0: Vec3, n: NodeIdx) {
    let now = w.conserved();
    let here = world_position(w, n);
    let drift = (here - where0).norm();
    println!(
        "  conserved: energy {:.3e}, baryon {:.3e}   |   did not teleport: {} ({:.2e} of its distance)",
        rel(conserved0.energy, now.energy),
        rel(conserved0.baryon, now.baryon),
        si(drift, "m"),
        drift / where0.norm().max(1e-300),
    );
}
