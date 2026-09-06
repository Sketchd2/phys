//! A world that outlives its process.
//!
//! Run it twice. The first run grows a tree, breaks it in a storm, and saves.
//! The second run loads that file and finds the same tree, still broken in the
//! same places — without the file having stored the tree.
//!
//! ```sh
//! cargo run --release --bin phys-persist          # grow, damage, save
//! cargo run --release --bin phys-persist          # load, verify, continue
//! cargo run --release --bin phys-persist -- reset # start again
//! ```

use phys::engine::{default_spec, World};
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::morph::{Environment, Program};
use phys::persist::{FileStore, WorldStore};
use phys::solvers::structure::weather;
use phys::state::Aggregate;
use phys::tree::Tree;
use phys::units::*;

const UPS: f64 = 20.0;

fn rule(title: &str) {
    println!("\n\x1b[1m{title}\x1b[0m");
    println!("{}", "─".repeat(title.len()));
}

fn bytes(n: usize) -> String {
    if n >= 1_000_000 {
        format!("{:.1} MB", n as f64 / 1e6)
    } else if n >= 1_000 {
        format!("{:.1} kB", n as f64 / 1e3)
    } else {
        format!("{n} B")
    }
}

/// A patch of ground with enough carbon in it to grow something.
fn a_meadow() -> Tree {
    let mut agg = Aggregate::neutral(4.0e4, 12.0, 291.0, Program::Tree.substrate());
    agg.internal_energy = agg.thermal_energy();
    Tree::new(0x5011, agg, Tier::Continuum, default_spec(Tier::Continuum))
}

fn main() {
    let path = std::path::PathBuf::from("world.phys");
    let mut store = FileStore::new(&path);

    if std::env::args().nth(1).as_deref() == Some("reset") {
        let _ = std::fs::remove_file(&path);
        println!("removed {}", path.display());
        return;
    }

    if path.exists() {
        second_run(&store, &path);
    } else {
        first_run(&mut store, &path);
    }
}

// ---------------------------------------------------------------------------

fn first_run(store: &mut FileStore, path: &std::path::Path) {
    rule("1. A world with something growing in it");

    let mut w = World::new(a_meadow(), UPS);
    let root = w.tree.root;
    w.plant(
        root,
        Program::Tree,
        Environment {
            light_flux: 340.0,
            temperature: 291.0,
            water: 0.75,
            crowding: 0.1,
            reservoir_mass: 4.0e4,
            labour: 0.0,
        },
    );
    println!("  planted a tree in a 12 m patch of ground");

    // Twenty years of growth, on the aggregate. The tree does not exist yet:
    // growth runs on bulk state, so this costs one ODE step per frame however
    // elaborate the thing being grown is.
    //
    // `paced_to` is cleared because the pace normally follows whatever is being
    // watched, and here we are driving it by hand.
    w.paced_to = NodeIdx::NONE;
    w.pace = 90.0 * 24.0 * 3600.0;
    w.time_rate = 1.0;
    for _ in 0..80 {
        w.step_frame(50_000.0);
    }
    let m = w.tree.nodes[root.get()].morphology.as_ref().unwrap();
    println!(
        "  grew it for {:.0} years: {:.1} m tall, {:.0} kg of wood, still {} bodies materialised",
        m.age / YEAR,
        w.tree.nodes[root.get()].agg.radius,
        m.built,
        w.tree.nodes[root.get()].bodies.len()
    );

    rule("2. Now make it real, and break it");

    let whole = w.tree.refine(root).len();
    let structural = w.tree.nodes[root.get()].last_report.structural_parts;
    let live_detail = w.tree.detail_bytes();
    println!("  materialised {whole} parts, {structural} of them structural");
    println!("  that is {} of geometry, resident", bytes(live_detail));

    let gale = weather::wind(38.0, v3(1.0, 0.0, 0.15));
    let out = w.damage(root, &[weather::gravity(), gale]);
    println!(
        "\n  a 38 m/s gale: {} joints broken, {} pieces came away, {:.2} kg of it",
        out.broken_joints, out.detached_pieces, out.detached_mass
    );

    // `damage` clears the body list on purpose. What it leaves behind is a set
    // of events on the morphology — this limb, severed, at this time — and the
    // geometry is regenerated from the program *plus* those events. So the
    // damage is a hundred bytes of history, not a million vertices of wreckage.
    let m = w.tree.nodes[root.get()].morphology.as_ref().unwrap();
    println!(
        "  recorded as {} events on the morphology; the geometry was dropped",
        m.events.len()
    );

    let before = w.conserved();

    rule("3. Save");

    store.save(w.view()).expect("save failed");
    let on_disk = store.size();
    println!("  wrote {} to {}", bytes(on_disk), path.display());
    println!(
        "  against {} of geometry the live world had been holding — {:.0}x smaller",
        bytes(live_detail),
        live_detail as f64 / on_disk.max(1) as f64
    );
    println!("\n  energy on the books: {:.6e} J", before.energy);
    println!("  baryons:             {:.6e}", before.baryon);
    println!("\n  Run it again — the process ends here, the world does not.");
}

// ---------------------------------------------------------------------------

fn second_run(store: &FileStore, path: &std::path::Path) {
    rule("1. Load");

    let snapshot = match store.load() {
        Ok(s) => s,
        Err(e) => {
            eprintln!("could not read {}: {e}", path.display());
            eprintln!("run with `reset` to start again");
            std::process::exit(1);
        }
    };
    let mut w = World::from_snapshot(snapshot, UPS);
    let root = w.tree.root;
    println!("  read {} from {}", bytes(store.size()), path.display());
    println!("  {} live nodes, world clock at {:.3} years", w.tree.live_count(), w.time / YEAR);

    rule("2. Is it the same world?");

    let after = w.conserved();
    println!("  energy on the books: {:.6e} J", after.energy);
    println!("  baryons:             {:.6e}", after.baryon);

    match w.tree.nodes[root.get()].morphology.as_ref() {
        Some(m) => {
            let severed = m
                .events
                .iter()
                .filter(|e| e.kind == phys::morph::EventKind::Severed)
                .count();
            println!(
                "\n  the tree is {:.0} years old, {:.0} kg of wood, {} events in its history",
                m.age / YEAR,
                m.built,
                m.events.len()
            );
            if severed > 0 {
                println!(
                    "  \x1b[1m{severed} of them are limbs it lost\x1b[0m — and it is still missing them"
                );
            }
        }
        None => println!("  no structure here (did the first run get that far?)"),
    }

    rule("3. What the file did not contain");

    println!("  {} of geometry on load", bytes(w.tree.detail_bytes()));
    let n = w.tree.refine(root).len();
    println!("  asked for it back: {n} parts, rebuilt from the program and its event log");
    println!("  {} now resident", bytes(w.tree.detail_bytes()));
    println!(
        "\n  The file is small because most of a world can be worked out again.\n  \
         It is not empty because the part somebody changed cannot — and the\n  \
         tree that comes back is the broken one, not a fresh copy of the program."
    );

    rule("4. Carry on");

    let t0 = w.time;
    w.pace = 30.0 * 24.0 * 3600.0;
    for _ in 0..20 {
        w.step_frame(50_000.0);
    }
    println!("  ran {:.1} more years of world time", (w.time - t0) / YEAR);
    println!("  save again to keep it, or `reset` to start over");
}
