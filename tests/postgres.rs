//! The world in a database.
//!
//! These need a live PostgreSQL and are skipped without one. Point `PHYS_PG` at
//! a database the test may freely truncate:
//!
//! ```sh
//! PHYS_PG='host=127.0.0.1 port=5432 user=phys dbname=phys' \
//!   cargo test --release --features postgres --test postgres
//! ```
//!
//! Skipping rather than failing is deliberate: a contributor without a database
//! should still get a green suite, and the property these tests check — that a
//! store is swappable — is only meaningful when there is something to swap to.

#![cfg(feature = "postgres")]

use phys::engine::{default_spec, galaxy, World};
use phys::math::v3;
use phys::observe::{Interaction, Quantity};
use phys::persist::{dirty_nodes, WorldStore};
use phys::store_pg::PostgresStore;
use phys::units::*;

fn store() -> Option<PostgresStore> {
    let url = std::env::var("PHYS_PG").ok()?;
    match PostgresStore::connect(&url) {
        Ok(mut s) => {
            s.clear().expect("truncate");
            Some(s)
        }
        Err(e) => {
            eprintln!("PHYS_PG is set but unusable ({e}); skipping");
            None
        }
    }
}

macro_rules! db {
    () => {
        match store() {
            Some(s) => s,
            None => {
                println!("  no PHYS_PG, skipped");
                return;
            }
        }
    };
}

fn a_world() -> World {
    let mut w = World::new(galaxy(0xDB, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 256;
    w.time_rate = 0.05;
    w
}

/// The same world, through a completely different store, comes back the same.
/// That is the swappability claim, and it is worth stating as a test rather
/// than as a trait.
#[test]
fn a_world_round_trips_through_postgres() {
    let mut pg = db!();
    let mut w = a_world();
    let root = w.tree.root;
    let path = w.drill(root, Tier::Continuum, &default_spec);
    let deep = *path.last().unwrap();
    w.plant(deep, phys::morph::Program::Tree, phys::morph::Environment::default());
    for _ in 0..6 {
        w.step_frame(50_000.0);
    }
    w.interact(Interaction::Impulse { target: deep, dp: v3(1.0, 2.0, 3.0) });
    w.interact(Interaction::Pin { target: deep });
    w.ledger.commit(w.tree.nodes[root.get()].key, Quantity::DecayTime, 99.5, 7.0);

    let before = w.conserved();
    pg.save(w.view()).expect("save");
    let rows = pg.node_count().expect("count");

    let back = World::from_snapshot(pg.load().expect("load"), 20.0);
    let after = back.conserved();

    println!("  {rows} node rows; baryons {:.6e} -> {:.6e}", before.baryon, after.baryon);
    assert_eq!(rows as usize, w.tree.nodes.len(), "a node did not become a row");
    assert_eq!(back.tree.nodes.len(), w.tree.nodes.len());
    assert_eq!(back.time.to_bits(), w.time.to_bits(), "world clock");
    assert_eq!(back.tree.world_seed, w.tree.world_seed);
    assert_eq!(back.tree.root, w.tree.root);
    assert!(back.tree.nodes[deep.get()].morphology.is_some(), "the structure was lost");
    assert!(back.tree.nodes[deep.get()].pinned, "the pin was lost");
    assert_eq!(back.ledger.len(), w.ledger.len(), "ledger");

    let drift = (after.baryon - before.baryon).abs() / before.baryon.abs().max(1e-300);
    assert!(drift < phys::tree::IDEMPOTENT_TOLERANCE, "baryon drift {drift:.3e}");
}

/// A file store and a database store must produce the same world. If they do
/// not, one of them is lying and the trait is not a swap point.
#[test]
fn postgres_and_a_file_agree() {
    let mut pg = db!();
    let mut w = a_world();
    let root = w.tree.root;
    w.drill(root, Tier::Planetary, &default_spec);
    for _ in 0..4 {
        w.step_frame(50_000.0);
    }

    let file = std::env::temp_dir().join(format!("phys-pg-cmp-{}.world", std::process::id()));
    let mut fs = phys::persist::FileStore::new(&file);
    fs.save(w.view()).expect("file save");
    pg.save(w.view()).expect("pg save");

    let a = World::from_snapshot(fs.load().expect("file load"), 20.0);
    let b = World::from_snapshot(pg.load().expect("pg load"), 20.0);

    assert_eq!(a.tree.nodes.len(), b.tree.nodes.len());
    let mut worst = 0u64;
    for (i, (x, y)) in a.tree.nodes.iter().zip(b.tree.nodes.iter()).enumerate() {
        assert_eq!(x.key, y.key, "node {i} key");
        assert_eq!(x.alive, y.alive, "node {i} alive");
        assert_eq!(x.epoch, y.epoch, "node {i} epoch");
        worst = worst.max(x.agg.mass.to_bits() ^ y.agg.mass.to_bits());
        worst = worst.max(x.frame.offset.x.to_bits() ^ y.frame.offset.x.to_bits());
        worst = worst.max(x.last_solved.to_bits() ^ y.last_solved.to_bits());
    }
    println!("  {} nodes, worst bit difference between backends: {worst}", a.tree.nodes.len());
    assert_eq!(worst, 0, "the two stores disagree about the world");
    let _ = std::fs::remove_file(&file);
}

/// Writes track *events*, not how much world there is. This is the property the
/// whole scheduler design was for, and the database is where it finally pays.
#[test]
fn writes_are_proportional_to_events() {
    let mut pg = db!();
    let mut w = a_world();
    let root = w.tree.root;
    w.drill(root, Tier::Planetary, &default_spec);
    for _ in 0..3 {
        w.step_frame(50_000.0);
    }
    pg.save(w.view()).expect("initial save");
    let total = w.tree.nodes.len();

    // Advance, then flush only what changed.
    let mark = w.time;
    for _ in 0..3 {
        w.step_frame(50_000.0);
    }
    let flushed = pg.save_since(w.view(), mark).expect("incremental save");

    println!(
        "  {total} nodes in the world; {} written, {} skipped as unchanged",
        flushed.nodes_written, flushed.nodes_skipped
    );
    assert!(flushed.incremental, "postgres should not fall back to a full write");
    assert_eq!(flushed.nodes_written + flushed.nodes_skipped, total);
    assert!(
        flushed.nodes_written < total,
        "every node was written; nothing was coasted, so the claim is untested here"
    );

    // The other direction matters as much: something that *did* happen must be
    // written. A store that skipped everything would pass the check above.
    let mark = w.time;
    let deep = *w.drill(root, Tier::Continuum, &default_spec).last().unwrap();
    w.interact(Interaction::Impulse { target: deep, dp: v3(9.0, 0.0, 0.0) });
    for _ in 0..4 {
        w.step_frame(50_000.0);
    }
    let after_event = pg.save_since(w.view(), mark).expect("save after an event");
    println!(
        "  after an impulse landed: {} written, {} skipped",
        after_event.nodes_written, after_event.nodes_skipped
    );
    assert!(
        after_event.nodes_written > 0,
        "something happened and the store wrote nothing"
    );

    // And skipping those writes lost nothing: the reader carries them forward.
    let back = World::from_snapshot(pg.load().expect("load"), 20.0);
    assert_eq!(back.time.to_bits(), w.time.to_bits(), "clock");
    for (i, (x, y)) in w.tree.nodes.iter().zip(back.tree.nodes.iter()).enumerate() {
        if !x.alive {
            continue;
        }
        assert_eq!(y.time.to_bits(), w.time.to_bits(), "node {i} was not settled to the instant");
        let slip = (x.frame.offset - y.frame.offset).norm();
        let scale = x.frame.offset.norm().max(x.agg.radius).max(1.0);
        assert!(
            slip <= scale * 1e-9,
            "node {i} was reconstructed {slip:.3e} m from where it should be"
        );
    }
}

/// The dirty set is computed from the node's own timestamps, so it needs no
/// extra bookkeeping and cannot drift out of step with the scheduler.
#[test]
fn the_dirty_set_is_derived_not_tracked() {
    let mut w = a_world();
    let root = w.tree.root;
    w.drill(root, Tier::Planetary, &default_spec);
    for _ in 0..3 {
        w.step_frame(50_000.0);
    }

    // Nothing has happened since now, so nothing is dirty.
    let quiet = dirty_nodes(&w.tree, w.time);
    assert!(quiet.is_empty(), "{} nodes claimed to be dirty with no time passed", quiet.len());

    // Touching one node makes it — and its ancestry, because a changed child
    // means the parent's sample no longer describes it — dirty. Note the
    // stepping: an impulse is *posted*, not applied, and arrives after the
    // light delay. Nothing is dirty until it lands, which is correct and is
    // exactly why the mailbox itself has to be durable.
    let mark = w.time;
    let deep = *w.drill(root, Tier::Continuum, &default_spec).last().unwrap();
    w.interact(Interaction::Impulse { target: deep, dp: v3(5.0, 0.0, 0.0) });
    assert!(w.mailbox.pending() > 0, "the impulse should be in flight, not applied");
    assert!(
        dirty_nodes(&w.tree, mark).is_empty(),
        "an influence marked nodes dirty before it arrived"
    );

    for _ in 0..4 {
        w.step_frame(50_000.0);
    }
    let touched = dirty_nodes(&w.tree, mark);
    println!(
        "  impulse landed after {} frames; {} of {} nodes dirty",
        4,
        touched.len(),
        w.tree.nodes.len()
    );
    assert!(!touched.is_empty(), "the impulse landed and marked nothing dirty");
}

/// An influence posted and not yet arrived is an action somebody took. A save
/// in the light-delay between the act and its landing must not lose it.
#[test]
fn an_impulse_in_flight_survives_a_save() {
    let mut pg = db!();
    let mut w = a_world();
    let root = w.tree.root;
    let deep = *w.drill(root, Tier::Planetary, &default_spec).last().unwrap();
    w.step_frame(50_000.0);

    w.interact(Interaction::Impulse { target: deep, dp: v3(11.0, -4.0, 2.0) });
    let flying = w.mailbox.pending();
    assert!(flying > 0, "nothing in flight to test");

    pg.save(w.view()).expect("save");
    let mut back = World::from_snapshot(pg.load().expect("load"), 20.0);
    assert_eq!(back.mailbox.pending(), flying, "an in-flight influence was dropped");

    // And it still lands, on the reloaded world, with the same effect.
    let before = back.tree.nodes[deep.get()].agg.momentum;
    for _ in 0..6 {
        w.step_frame(50_000.0);
        back.step_frame(50_000.0);
    }
    let after = back.tree.nodes[deep.get()].agg.momentum;
    println!(
        "  {flying} in flight across the save; momentum {:.4e} -> {:.4e} after it landed",
        before.norm(),
        after.norm()
    );
    assert_eq!(
        w.tree.nodes[deep.get()].agg.momentum.x.to_bits(),
        after.x.to_bits(),
        "the reloaded world applied the impulse differently"
    );
}

/// The node table is queryable, which is the entire reason for not storing a
/// blob. A shard will be a `WHERE` clause over these columns.
#[test]
fn nodes_are_rows_you_can_ask_questions_about() {
    let mut pg = db!();
    let mut w = a_world();
    let root = w.tree.root;
    w.drill(root, Tier::Molecular, &default_spec);
    w.step_frame(50_000.0);
    pg.save(w.view()).expect("save");

    let rows = pg
        .client()
        .query(
            "SELECT tier, count(*) AS n FROM node WHERE alive GROUP BY tier ORDER BY tier",
            &[],
        )
        .expect("group by tier");
    println!("  live nodes by tier:");
    let mut total = 0i64;
    for r in &rows {
        let tier: i16 = r.get("tier");
        let n: i64 = r.get("n");
        total += n;
        println!("    {:<12} {n}", phys::view::tier_name(tier as u8));
    }
    assert!(rows.len() >= 4, "expected several tiers, got {}", rows.len());
    assert_eq!(total, w.tree.live_count() as i64, "the query disagrees with the tree");

    // The children of a node, without reading any other node.
    let kids: i64 = pg
        .client()
        .query_one("SELECT count(*) FROM node WHERE parent = $1", &[&(root.0 as i32)])
        .expect("children query")
        .get(0);
    println!("  the root has {kids} children, fetched without reading the tree");
    assert!(kids >= 1);
}
