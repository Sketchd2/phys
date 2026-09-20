//! What Phase 3's two per-frame passes cost, and what the phase measured.
//!
//! `cargo run --release --example phase3`

use phys::engine::{default_spec, galaxy, World};
use std::time::Instant;

fn rule(title: &str) {
    println!("\n\x1b[1m{}\x1b[0m", title);
    println!("{}", "-".repeat(title.len()));
}

/// A world with `n` live nodes, built by promoting siblings out of one parent.
fn a_world_of(n: usize) -> World {
    let mut w = World::new(galaxy(0x9E77, 1e9), 20.0);
    w.tree.nodes[0].spec.count = (n + 8).max(64);
    let root = w.tree.root;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    for slot in 0..n {
        w.tree.promote(root, slot, default_spec(tier.finer()));
    }
    w
}

fn main() {
    rule("the crossing pass, over a world nothing is leaving");
    println!("  {:>8}  {:>12}  {:>14}", "nodes", "per pass (us)", "per node (ns)");
    for n in [16usize, 128, 1024, 8192] {
        let mut w = a_world_of(n);
        // Warm: the first call pays for whatever the pass touches first.
        w.cross_boundaries();
        let reps = 200;
        let t = Instant::now();
        for _ in 0..reps {
            w.cross_boundaries();
        }
        let us = t.elapsed().as_secs_f64() * 1e6 / reps as f64;
        let live = w.tree.live_count();
        println!(
            "  {live:>8}  {us:>12.2}  {:>14.1}",
            us * 1e3 / live as f64
        );
    }

    rule("a crossing, end to end");
    let mut w = a_world_of(64);
    let root = w.tree.root;
    // Push one node out of its parent's reach and time the frame that notices.
    let child = w.tree.nodes[root.get()].children.iter().copied().find(|c| !c.is_none()).unwrap();
    let reach = w.tree.contents_reach(root, child);
    w.tree.nodes[child.get()].motion.offset =
        w.tree.nodes[child.get()].motion.offset.unit().scale(reach * 4.0);
    let before = w.stats.crossings_refused;
    let t = Instant::now();
    w.cross_boundaries();
    let us = t.elapsed().as_secs_f64() * 1e6;
    println!(
        "  one node past its parent's reach: {us:.2} us for the pass, {} refusals \
         (it is a child of the root, which has no outside)",
        w.stats.crossings_refused - before
    );

    rule("the extent pass, per node advanced");
    println!("  {:>8}  {:>10}  {:>14}", "bodies", "components (us)", "components");
    for count in [256usize, 1024, 4096, 16384] {
        let mut w = World::new(galaxy(0x1234, 1e9), 20.0);
        w.tree.nodes[0].spec.count = count;
        let root = w.tree.root;
        w.tree.refine(root);
        let tier = w.tree.nodes[root.get()].tier;
        let node = w.tree.promote(root, 0, default_spec(tier.finer()));
        w.tree.nodes[node.get()].spec.count = count;
        w.tree.refine(node);
        let (_, _, components) = w.tree.components(node);
        let reps = 20;
        let t = Instant::now();
        for _ in 0..reps {
            let _ = w.tree.components(node);
        }
        let us = t.elapsed().as_secs_f64() * 1e6 / reps as f64;
        // And what advancing the same node costs, which is what the pass
        // rides on: it only ever runs for a node the frame advanced.
        let dt = w.node_dt(node);
        let t = Instant::now();
        for _ in 0..reps {
            w.advance_node(node, dt);
        }
        let solve = t.elapsed().as_secs_f64() * 1e6 / reps as f64;
        println!(
            "  {count:>8}  {us:>10.1}  {components:>14}   solve {solve:>9.1} us, \
             the check is {:>5.1}% of it",
            100.0 * us / solve.max(1e-9)
        );
    }
}
