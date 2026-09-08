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
//! costs per frame, under each of the four things the seam can do about it.

use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::units::*;
use phys::view::{self, kind_name, solver_name, tier_name, Client, Scene, ViewRequest, Volume};

const FRAMES: usize = 60;
const STREAM: &str = "scenes.phys";
/// Child nodes to hang off the watched node. The "city": many things nearby,
/// almost none of them doing anything.
const NEIGHBOURS: usize = 200;
/// The client's body budget for the whole scene. Every row of the table gets
/// the same one, so the table measures mechanism rather than generosity.
const BUDGET: usize = 20_000;

/// One measured way of answering the same request.
struct Row {
    name: &'static str,
    join: usize,
    total: usize,
    frames: usize,
    /// Node-frames of detail actually delivered, so a row cannot look good by
    /// saying nothing.
    updates: usize,
}

impl Row {
    fn new(name: &'static str) -> Row {
        Row { name, join: 0, total: 0, frames: 0, updates: 0 }
    }

    fn take(&mut self, s: &Scene) {
        let n = view::encode(s).len();
        if self.frames == 0 {
            self.join = n;
        }
        self.frames += 1;
        self.total += n;
        self.updates += s.nodes.iter().filter(|v| !v.detail.is_unchanged()).count();
    }

    /// Bytes per frame after the first.
    fn steady(&self) -> f64 {
        if self.frames < 2 {
            return self.join as f64;
        }
        (self.total - self.join) as f64 / (self.frames - 1) as f64
    }

    fn print(&self, against: f64) {
        println!(
            "  {:<26} {:>7.1} kB {:>7.1} kB {:>7.2} MB/s {:>9}   {:>5.1}x",
            self.name,
            self.join as f64 / 1e3,
            self.steady() / 1e3,
            self.steady() * 20.0 / 1e6,
            self.updates,
            against / self.steady().max(1.0),
        );
    }
}

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

/// A galaxy, drilled to a stellar node, with neighbours promoted around it.
///
/// The neighbours are the point. Each is a real node with its own matter
/// that the engine has *not* materialised, which is the shape a city has: many
/// things present, few of them being solved. What a client can be told about
/// them is the entire question this program measures.
fn build() -> (World, NodeIdx) {
    let mut w = World::new(galaxy(0xA11A5, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 4000;
    let root = w.tree.root;
    let path = w.drill_to(root, Tier::Stellar.max_radius(), &default_spec);
    let watched = *path.last().unwrap();
    w.tree.refine(watched);

    // Promote a spread of the watched node's own bodies into nodes. Every
    // k-th rather than the first k: a prefix of a sampled population all sits
    // in one place, and a city with every building in one corner would flatter
    // the LOD numbers.
    let have = w.tree.nodes[watched.get()].bodies.len();
    let spec = default_spec(w.tree.nodes[watched.get()].tier.finer());
    let stride = (have / NEIGHBOURS.max(1)).max(1);
    for k in 0..NEIGHBOURS {
        let slot = k * stride;
        if slot >= have {
            break;
        }
        w.tree.promote(watched, slot, spec);
    }
    w.pace_to(watched);
    (w, watched)
}

fn serve() {
    rule("Solving");
    let (mut w, watched) = build();
    let kids = w.tree.nodes[watched.get()].children.iter().filter(|c| !c.is_none()).count();
    println!(
        "  a galaxy, drilled to a {} node with {kids} neighbours promoted around it",
        tier_name(w.tree.nodes[watched.get()].tier.index() as u8)
    );
    println!("  writing {FRAMES} frames to {STREAM}");

    // Four ways to answer the *same* question, measured against each other.
    // Same volume, same body budget, same 60 frames — only the mechanism
    // differs, which is the only way the comparison means anything.
    //
    // Two numbers per row, because the average of the two hides both. The
    // first frame is what joining costs: the client holds nothing and must be
    // told everything. Every frame after that is what *running* costs, and it
    // is the one a server sizes its uplink against.
    let mut naive = Row::new("every body, every frame");
    let mut on_solve = Row::new("only when re-solved");
    let mut delta = Row::new("per-node deltas");
    let mut full = Row::new("+ LOD and recipes");

    let mut since = f64::NEG_INFINITY;
    let mut plain = Client::new();
    let mut smart = Client::new();
    let mut solve_us = 0.0;
    let mut render_us = 0.0;

    // Standing back a few node radii, looking at the whole neighbourhood. A
    // real client would move this every frame; holding it still keeps the
    // numbers comparable.
    let eye = [0.0, 0.0, -6.0];
    // `detail_angle: 0` gathers the same nodes and details every one of them,
    // so the first three rows differ from the fourth only in LOD and recipes.
    let wide = Volume { eye, reach: 40.0, detail_angle: 0.0, max_nodes: 256 };
    let lod = Volume { eye, reach: 40.0, detail_angle: 2e-3, max_nodes: 256 };
    let base = ViewRequest {
        node: watched,
        max_bodies: BUDGET,
        volume: Some(wide),
        allow_recipes: false,
        ..Default::default()
    };

    let mut stream: Vec<u8> = Vec::new();
    for _ in 0..FRAMES {
        let t = std::time::Instant::now();
        w.step_frame(50_000.0);
        solve_us += t.elapsed().as_secs_f64() * 1e6;

        let t = std::time::Instant::now();

        naive.take(&w.render(&base));

        let s = w.render(&ViewRequest { since, ..base.clone() });
        if s.nodes.iter().any(|n| !n.detail.is_unchanged()) {
            since = w.time;
        }
        on_solve.take(&s);

        delta.take(&plain.frame(&w, &base));

        let scene = smart.frame(
            &w,
            &ViewRequest { volume: Some(lod), allow_recipes: true, ..base.clone() },
        );
        let bytes = view::encode(&scene);
        full.take(&scene);
        render_us += t.elapsed().as_secs_f64() * 1e6;

        stream.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
        stream.extend_from_slice(&bytes);
    }

    if let Err(e) = std::fs::write(STREAM, &stream) {
        eprintln!("could not write {STREAM}: {e}");
        std::process::exit(1);
    }

    rule("What that cost");
    println!("  solve:  {:>8.2} ms/frame", solve_us / FRAMES as f64 / 1e3);
    println!(
        "  render: {:>8.2} ms/frame  ({:.1}% of the solve, for four scenes a frame)",
        render_us / FRAMES as f64 / 1e3,
        100.0 * render_us / solve_us.max(1e-9)
    );

    rule("Bandwidth");
    println!(
        "  {:<26} {:>9} {:>10} {:>10} {:>9}",
        "", "join", "steady", "at 20fps", "updates"
    );
    for r in [&naive, &on_solve, &delta, &full] {
        r.print(naive.steady());
    }
    println!(
        "\n  \x1b[1mupdates\x1b[0m is node-frames of detail actually delivered, because a\n  \
         scheme can always be cheap by being silent. The second row is exactly\n  \
         that: one `since` instant gates every node in the query, so a node that\n  \
         *was* re-solved goes unmentioned whenever the one the client asked\n  \
         about was not. It buys its bytes with detail the client never gets.\n  \
         Per-node tracking is both cheaper and right, because it also stops\n  \
         re-sending the facts of things nobody is touching — which is what the\n  \
         steady state of a crowded scene is almost entirely made of.\n\n  \
         The last row costs *more* to join than the first. That is not a\n  \
         regression: it is delivering two hundred nodes the first row cannot\n  \
         describe at all, because the engine never materialised them and has no\n  \
         bodies to send. Three hundred bytes of recipe each is what that costs,\n  \
         once."
    );

    println!("\n  Now run `watch`. Nothing it does will touch this engine.");
}

// ---------------------------------------------------------------------------
// the view side — no engine below this line
// ---------------------------------------------------------------------------

/// What the client holds for one node between frames.
struct Holding {
    specks: Vec<phys::view::Speck>,
    at: f64,
    /// True when the client built these itself out of a recipe.
    grown: bool,
}

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
    let mut held: std::collections::HashMap<u128, Holding> = std::collections::HashMap::new();
    // The client's own copy of every node's facts. The server sends them once
    // and then stops, so a client that wants to know what it is drawing has to
    // remember. `merge_facts` is both halves of that.
    let mut known: std::collections::HashMap<phys::ids::PathKey, phys::view::NodeFacts> =
        std::collections::HashMap::new();
    let mut filled = 0usize;
    let mut sent = 0usize;
    let mut grown = 0usize;
    let mut carried = 0usize;
    let mut dropped = 0usize;
    let mut refused = 0usize;

    while at + 4 <= bytes.len() {
        let n = u32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]]) as usize;
        at += 4;
        if at + n > bytes.len() {
            eprintln!("  stream truncated at frame {frames}");
            break;
        }
        let mut scene = match view::decode(&bytes[at..at + n]) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("  frame {frames} did not decode: {e}");
                break;
            }
        };
        at += n;
        frames += 1;
        filled += scene.merge_facts(&mut known);

        // The client's half of the bargain, in four parts.
        //
        // A recipe is built here, by the same sampler the engine uses, from an
        // matter that arrived in a few hundred bytes. If this build cannot
        // reproduce it the scene says so rather than drawing something wrong,
        // and a real client would re-ask with `allow_recipes: false`.
        let from_recipe: std::collections::HashSet<u128> = scene
            .nodes
            .iter()
            .filter(|n| matches!(n.detail, view::Detail::Recipe(_)))
            .map(|n| n.facts.key.0)
            .collect();
        if let Err(e) = scene.materialise() {
            refused += from_recipe.len();
            eprintln!("  frame {frames}: {e} — would re-ask without recipes");
        }

        // A node that left the query is released. This is the only reason the
        // client's memory does not grow without bound as it moves.
        dropped += scene.dropped.len();
        for k in &scene.dropped {
            held.remove(&k.0);
        }
        scene.forget_dropped(&mut known);

        // And a node with nothing new is carried forward at its own velocities:
        // the same closed-form step the engine performs internally, which is
        // what makes the server's silence affordable.
        for v in &scene.nodes {
            match &v.detail {
                view::Detail::Unchanged => {
                    if let Some(h) = held.get_mut(&v.facts.key.0) {
                        let dt = (scene.instant - h.at) as f32;
                        for b in h.specks.iter_mut() {
                            b.pos = b.at(dt);
                        }
                        h.at = scene.instant;
                        carried += 1;
                    }
                }
                view::Detail::Explicit(b) => {
                    let from_recipe = from_recipe.contains(&v.facts.key.0);
                    if from_recipe {
                        grown += 1;
                    } else {
                        sent += 1;
                    }
                    held.insert(
                        v.facts.key.0,
                        Holding { specks: b.clone(), at: scene.instant, grown: from_recipe },
                    );
                }
                // Turned into `Explicit` by `materialise` above, unless this
                // build could not reproduce it — in which case there is
                // nothing to hold and the node draws as matter alone.
                view::Detail::Recipe(_) => {}
            }
        }

        drawn += held.values().map(|h| h.specks.len() as u64).sum::<u64>();
        if first.is_none() {
            first = Some(scene.clone());
        }
        last = Some(scene);
    }

    let (Some(first), Some(last)) = (first, last) else {
        eprintln!("  no complete frames in the stream");
        std::process::exit(1);
    };

    rule("What a client can say about a world it cannot touch");
    println!(
        "  {frames} frames, {drawn} specks drawn across {} nodes it is holding",
        held.len()
    );
    println!(
        "  {sent} node updates arrived as bodies, {grown} as recipes the client built\n  \
         itself, {carried} were carried forward from what it already had, {dropped} released.\n  \
         {filled} node-frames arrived as a bare key, their facts filled in from what\n  \
         the client already knew — which is most of what the steady state saves."
    );
    if refused > 0 {
        println!("  \x1b[1m{refused} recipes did not check out\x1b[0m — this build samples differently");
    }

    let n = last.node();
    println!("\n  node:    {} tier, solved by {}", tier_name(n.tier), solver_name(n.solver));
    println!("  mass:    {}", si(n.mass, "kg"));
    println!("  radius:  {}", si(n.radius, "m"));
    println!("  cadence: {} between solves", si(n.cadence, "s"));
    println!(
        "  forgets: {}",
        if n.mixing_time.is_finite() {
            si(n.mixing_time, "s")
        } else {
            "never — structured matter".to_string()
        }
    );

    let hold = held.get(&n.key.0);
    let specks: &[phys::view::Speck] = hold.map(|h| &h.specks[..]).unwrap_or(&[]);
    if let Some(h) = hold {
        println!(
            "  source:  {}",
            if h.grown { "built here, from a recipe" } else { "sent as bodies" }
        );
    }
    let kinds: std::collections::BTreeSet<&str> = specks.iter().map(|b| kind_name(b.kind)).collect();
    println!("  holding: {}", kinds.into_iter().collect::<Vec<_>>().join(", "));

    // Node radii per second: everything a client sees is scaled to the node,
    // which is what lets one renderer draw a galaxy and a nucleus.
    let (lo, hi) = specks
        .iter()
        .map(|b| b.speed())
        .fold((f32::INFINITY, 0.0f32), |(l, h), v| (l.min(v), h.max(v)));
    if lo.is_finite() {
        println!(
            "  speeds:  {:.3e} to {:.3e} node radii/s  ({} to {} at this scale)",
            lo,
            hi,
            si(lo as f64 * n.radius, "m/s"),
            si(hi as f64 * n.radius, "m/s")
        );
    }

    rule("And what else is nearby");
    let mut sizes: Vec<(f32, &phys::view::NodeView)> =
        last.nodes.iter().map(|v| (v.angular_size, v)).collect();
    sizes.sort_by(|a, b| b.0.total_cmp(&a.0));
    println!(
        "  {} nodes in the query, {} of them detailed",
        last.nodes.len(),
        last.nodes.iter().filter(|v| !v.detail.is_unchanged()).count()
    );
    for (a, v) in sizes.iter().take(6) {
        println!(
            "    {:>10} at {:>7.2} radii, {:.4} rad across — {}",
            tier_name(v.facts.tier),
            (v.offset[0].powi(2) + v.offset[1].powi(2) + v.offset[2].powi(2)).sqrt(),
            a,
            if v.detail.is_unchanged() { "nothing new".to_string() } else { format!("{} specks", v.detail.specks().len()) }
        );
    }

    rule("And that the world moved");
    println!("  world clock: {} to {}", si(first.instant, "s"), si(last.instant, "s"));
    println!("  {} of world time across {frames} frames", si(last.instant - first.instant, "s"));
    println!(
        "  live nodes {} → {}, coasted {} per frame, worst lateness {:.3}",
        first.world.live_nodes, last.world.live_nodes, last.world.coasted, last.world.worst_lateness
    );
    if last.world.unreachable > 0 {
        println!(
            "  {} node-frames could not be integrated at this pace and could not be
               crossed by ensemble either — something is watching them. Their lateness
               does not recover on its own.",
            last.world.unreachable
        );
    }

    println!(
        "\n  Every number above came out of a byte stream. This half of the\n  \
         program has no World, no Tree and no solver — it could not advance\n  \
         the simulation if it wanted to. It *can* draw scenery the server\n  \
         never built, because a recipe is instructions rather than a picture."
    );
}
