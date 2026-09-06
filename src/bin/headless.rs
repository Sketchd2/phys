//! The authoritative loop, and a client that never sees it.
//!
//! Two halves of what used to be one program, run separately:
//!
//! ```sh
//! cargo run --release --bin phys-headless -- serve   # solve, write scenes
//! cargo run --release --bin phys-headless -- watch   # read scenes, draw
//! ```
//!
//! `serve` runs the world and writes a scene per frame. `watch` reads those
//! bytes back and reports what a renderer would draw — with no `World`, no
//! `Tree`, and no solver anywhere in its call graph. It could as easily be
//! reading from a socket; a file is just the cheapest process boundary that
//! proves the point.
//!
//! What this measures is the number Phase 3 needs: what one client actually
//! costs per frame.

use phys::engine::{default_spec, galaxy, World};
use phys::units::*;
use phys::view::{self, kind_name, solver_name, tier_name, Scene, ViewRequest};

const FRAMES: usize = 60;
const STREAM: &str = "scenes.phys";

fn rule(title: &str) {
    println!("\n\x1b[1m{title}\x1b[0m");
    println!("{}", "─".repeat(title.chars().count()));
}

fn si(v: f64, unit: &str) -> String {
    const P: [(f64, &str); 11] = [
        (1e24, "Y"),
        (1e21, "Z"),
        (1e18, "E"),
        (1e15, "P"),
        (1e12, "T"),
        (1e9, "G"),
        (1e6, "M"),
        (1e3, "k"),
        (1.0, ""),
        (1e-3, "m"),
        (1e-6, "µ"),
    ];
    if !v.is_finite() {
        return format!("∞ {unit}");
    }
    let a = v.abs();
    // Past yotta there is no prefix and a galaxy weighs a billion of them, so
    // fall back to an exponent rather than printing "6285263048.58 Ykg".
    if a >= 1e27 || (a > 0.0 && a < 1e-6) {
        return format!("{v:.3e} {unit}");
    }
    for (s, p) in P {
        if a >= s {
            return format!("{:.2} {p}{unit}", v / s);
        }
    }
    format!("{v:.3e} {unit}")
}

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("watch") => watch(),
        Some("serve") | None => serve(),
        Some(other) => {
            eprintln!("unknown command {other:?}; expected `serve` or `watch`");
            std::process::exit(2);
        }
    }
}

// ---------------------------------------------------------------------------
// the solve side
// ---------------------------------------------------------------------------

fn serve() {
    rule("Solving");

    let mut w = World::new(galaxy(0xA11A5, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 4000;
    let root = w.tree.root;
    let path = w.drill(root, Tier::Stellar, &default_spec);
    let watched = *path.last().unwrap();
    w.tree.refine(watched);
    w.pace_to(watched);
    println!("  a galaxy, drilled to a {} node", tier_name(w.tree.nodes[watched.get()].tier.index() as u8));
    println!("  writing {FRAMES} frames to {STREAM}");

    // One scene per frame, length-prefixed so the reader can walk them.
    let mut stream: Vec<u8> = Vec::new();
    let mut solve_us = 0.0;
    let mut render_us = 0.0;

    // What a client that re-asks for everything every frame would cost, against
    // what one that carries its own bodies forward costs. The difference is the
    // entire answer to "this scales badly with scene complexity".
    let mut naive_bytes = 0usize;
    let mut sent_since = f64::NEG_INFINITY;
    let mut updates = 0usize;

    for _ in 0..FRAMES {
        let t = std::time::Instant::now();
        w.step_frame(50_000.0);
        solve_us += t.elapsed().as_secs_f64() * 1e6;

        let t = std::time::Instant::now();
        let everything =
            w.render(&ViewRequest { node: watched, max_bodies: 0, trail: 16, since: f64::NEG_INFINITY });
        naive_bytes += view::encode(&everything).len();

        // What we actually send: bodies only when the node has been re-solved.
        let scene =
            w.render(&ViewRequest { node: watched, max_bodies: 0, trail: 16, since: sent_since });
        if scene.bodies_included {
            sent_since = w.time;
            updates += 1;
        }
        let bytes = view::encode(&scene);
        render_us += t.elapsed().as_secs_f64() * 1e6;

        stream.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
        stream.extend_from_slice(&bytes);
    }

    if let Err(e) = std::fs::write(STREAM, &stream) {
        eprintln!("could not write {STREAM}: {e}");
        std::process::exit(1);
    }

    rule("What that cost");
    let per_frame = stream.len() as f64 / FRAMES as f64;
    let naive_per_frame = naive_bytes as f64 / FRAMES as f64;
    println!("  solve:  {:>8.2} ms/frame", solve_us / FRAMES as f64 / 1e3);
    println!(
        "  render: {:>8.2} ms/frame  ({:.1}% of the solve)",
        render_us / FRAMES as f64 / 1e3,
        100.0 * render_us / solve_us.max(1e-9)
    );

    rule("Bandwidth");
    println!("  re-sending every body every frame:");
    println!(
        "    {:>8.1} kB/frame   {:.2} MB/s at 20 fps",
        naive_per_frame / 1e3,
        naive_per_frame * 20.0 / 1e6
    );
    println!("  sending only when the node is re-solved ({updates} of {FRAMES} frames):");
    println!(
        "    {:>8.1} kB/frame   {:.2} MB/s at 20 fps   \x1b[1m{:.1}x less\x1b[0m",
        per_frame / 1e3,
        per_frame * 20.0 / 1e6,
        naive_per_frame / per_frame.max(1.0)
    );
    println!(
        "\n  The saving is not compression. A node nobody re-solved has bodies the\n  \
         client can carry forward itself, so there is nothing to say about it —\n  \
         which makes traffic track what is *happening* rather than what is *visible*.\n  \
         A city of still buildings costs nothing after the first frame."
    );
    println!("\n  Now run `watch`. Nothing it does will touch this engine.");
}

// ---------------------------------------------------------------------------
// the view side — no engine below this line
// ---------------------------------------------------------------------------

fn watch() {
    let bytes = match std::fs::read(STREAM) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("could not read {STREAM}: {e}");
            eprintln!("run `serve` first");
            std::process::exit(1);
        }
    };

    rule("Watching");
    println!("  {} of scenes, and no world to ask", si(bytes.len() as f64, "B"));

    let mut at = 0usize;
    let mut frames = 0usize;
    let mut first: Option<Scene> = None;
    let mut last: Option<Scene> = None;
    let mut drawn = 0u64;
    let mut held: Vec<phys::view::Speck> = Vec::new();
    let mut held_at = 0.0f64;
    let mut updates = 0usize;
    let mut carried = 0usize;

    while at + 4 <= bytes.len() {
        let n = u32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]]) as usize;
        at += 4;
        if at + n > bytes.len() {
            eprintln!("  stream truncated at frame {frames}");
            break;
        }
        let scene = match view::decode(&bytes[at..at + n]) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("  frame {frames} did not decode: {e}");
                break;
            }
        };
        at += n;
        frames += 1;

        // The client's half of the bargain. A frame that carries no bodies is
        // not an empty frame — it means nothing was re-solved, so what we
        // already hold is still right, carried forward by its own velocities.
        // This is the same closed-form step the engine performs internally, and
        // doing it here is what makes the server's silence affordable.
        if scene.bodies_included {
            held = scene.bodies.clone();
            updates += 1;
            held_at = scene.instant;
        } else {
            let dt = (scene.instant - held_at) as f32;
            for b in held.iter_mut() {
                b.pos = b.at(dt);
            }
            held_at = scene.instant;
            carried += 1;
        }
        drawn += held.len() as u64;

        let mut shown = scene.clone();
        shown.bodies = held.clone();
        if first.is_none() {
            first = Some(shown.clone());
        }
        last = Some(shown);
    }

    let (Some(first), Some(last)) = (first, last) else {
        eprintln!("  no complete frames in the stream");
        std::process::exit(1);
    };

    rule("What a client can say about a world it cannot touch");
    println!(
        "  {frames} frames, {drawn} bodies drawn — {updates} frames brought new bodies,\n  \
         {carried} were carried forward by the client from what it already had"
    );
    println!(
        "\n  node:    {} tier, solved by {}",
        tier_name(last.node.tier),
        solver_name(last.node.solver)
    );
    println!("  mass:    {}", si(last.node.mass, "kg"));
    println!("  radius:  {}", si(last.node.radius, "m"));
    println!("  cadence: {} between solves", si(last.node.cadence, "s"));
    println!(
        "  forgets: {}",
        if last.node.mixing_time.is_finite() {
            si(last.node.mixing_time, "s")
        } else {
            "never — structured matter".to_string()
        }
    );

    let kinds: std::collections::BTreeSet<&str> =
        last.bodies.iter().map(|b| kind_name(b.kind)).collect();
    println!("  holding: {}", kinds.into_iter().collect::<Vec<_>>().join(", "));

    // Node radii per second: everything a client sees is scaled to the node,
    // which is what lets one renderer draw a galaxy and a nucleus.
    let (lo, hi) = last.channel_range(|b| b.speed());
    println!(
        "  speeds:  {:.3e} to {:.3e} node radii/s  ({} to {} at this scale)",
        lo,
        hi,
        si(lo as f64 * last.node.radius, "m/s"),
        si(hi as f64 * last.node.radius, "m/s")
    );

    rule("And that the world moved");
    println!(
        "  world clock: {} to {}",
        si(first.instant, "s"),
        si(last.instant, "s")
    );
    println!("  {} of world time across {frames} frames", si(last.instant - first.instant, "s"));
    println!(
        "  live nodes {} → {}, coasted {} per frame, worst lateness {:.3}",
        first.world.live_nodes, last.world.live_nodes, last.world.coasted, last.world.worst_lateness
    );

    // The bodies moved between the first frame and the last, which is the only
    // evidence a client ever has that anything is happening.
    let moved = first
        .bodies
        .iter()
        .zip(last.bodies.iter())
        .map(|(a, b)| {
            let d = [b.pos[0] - a.pos[0], b.pos[1] - a.pos[1], b.pos[2] - a.pos[2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .fold(0.0f32, f32::max);
    println!("  furthest a body travelled: {moved:.4} node radii");

    println!(
        "\n  Every number above came out of a byte stream. This half of the\n  \
         program has no World, no Tree and no solver — it could not advance\n  \
         the simulation if it wanted to."
    );
}
