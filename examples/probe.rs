//! Phase 0 probes: the measurements `docs/PLAY.md` §6 says must exist before
//! any of it is designed further.
//!
//! Each probe answers one row of that table and prints what it measured. The
//! numbers here are what `docs/PERFORMANCE.md` and `docs/PLAY.md` cite; running
//! this is how they are checked rather than remembered. `tests/probes.rs` pins
//! the same behaviours as assertions so a change has to walk past them.
//!
//! Probes that need machinery the plan has not built yet — a neighbour index, a
//! gait solver, a derived erosion rate — say so and report nothing, rather than
//! reporting a number that came from somewhere else.

use phys::engine::{default_spec, galaxy, World, MAX_SUBSTEPS};
use phys::morph::{Environment, Event, EventKind, Morphology, Program};
use phys::math::v3;
use phys::rng::{Purpose, Stream};
use phys::units::*;

const PROGRAMS: [Program; 6] = [
    Program::Tree,
    Program::Coral,
    Program::Tower,
    Program::Wall,
    Program::Terrain,
    Program::Settlement,
];

fn rule(title: &str) {
    println!("\n\x1b[1m{}\x1b[0m", title);
    println!("{}", "-".repeat(title.len()));
}

/// A morphology grown or built to a usable size, so a probe has something real
/// to measure. Growth programs eat light and water; planned ones eat labour.
fn matured(program: Program, seed: u64) -> Morphology {
    // A planned program builds nothing without a design mass: `advance` returns
    // early on `design_mass <= 0.0`. Growth programs leave it zero by contract.
    let mut m = if program.is_planned() {
        Morphology::planned(program, 5.0e4, seed, 0x1234)
    } else {
        Morphology::new(program, seed, 0x1234, 0)
    };
    let mut env = Environment::default();
    env.labour = 0.02;
    env.reservoir_mass = 1.0e9;
    for _ in 0..400 {
        m.advance(YEAR * 0.25, &env);
    }
    m
}

// ---------------------------------------------------------------------------
// §3.4 — where does the trajectory path run out, at one second per second?
// ---------------------------------------------------------------------------

fn substeps_per_tier() {
    rule("§3.4  substeps per frame by tier, at 1 s/s");
    let mut w = World::new(galaxy(0xAC70, 1e9), 20.0);
    w.pace_fixed(1.0 / 20.0);
    let span = w.frame_dt();
    println!("frame covers {:.4} s of world time; MAX_SUBSTEPS = {}\n", span, MAX_SUBSTEPS);
    println!("{:<12} {:>10} {:>8} {:>13} {:>13}", "tier", "radius", "parts", "node_dt (s)", "substeps");

    let root = w.tree.root;
    for idx in w.drill_to(root, 1e-3, &default_spec) {
        let (tier, r) = {
            let n = &w.tree.nodes[idx.get()];
            (n.tier, n.matter.radius)
        };
        w.tree.refine(idx);
        let parts = {
            let n = &w.tree.nodes[idx.get()];
            if n.bodies.is_empty() { n.spec.count } else { n.bodies.len() }
        };
        let dt = w.node_dt(idx);
        let need = if dt > 0.0 { span / dt } else { f64::INFINITY };
        println!(
            "{:<12} {:>10.3e} {:>8} {:>13.3e} {:>13.1}{}",
            tier.name(), r, parts, dt, need,
            if need > MAX_SUBSTEPS as f64 { "  <- ensemble" } else { "" }
        );
    }

    println!("\nfinest h still followed, h = 4 * span * c / MAX_SUBSTEPS:");
    for (name, c) in [("air 340 m/s", 340.0), ("water 1500 m/s", 1500.0), ("rock 5000 m/s", 5000.0)] {
        println!("  {:<16} {:.4} m", name, 4.0 * span * c / MAX_SUBSTEPS as f64);
    }
}

// ---------------------------------------------------------------------------
// §5.9 — at what edit count does a structure regrow what was removed?
// ---------------------------------------------------------------------------

fn event_cap_regrowth() {
    rule("§5.9  severed parts that return, by program");
    println!("{:<12} {:>8} {:>12} {:>12} {:>14}", "program", "intact", "at 64 sev", "at 65 sev", "back at 65");
    for program in PROGRAMS {
        let mut m = matured(program, 0xBEEF);
        let budget = 512;
        let base = m.render(budget);
        let sites: Vec<u32> = base.site.iter().copied().take(70).collect();
        if sites.len() < 66 {
            println!("{:<12} {:>8} {:>12} {:>12} {:>14}", program.name(), base.site.len(), "-", "-", "too few parts");
            continue;
        }
        let mut at64 = 0usize;
        let mut at65 = 0usize;
        let mut back = 0usize;
        for (i, s) in sites.iter().enumerate().take(65) {
            m.record(Event { at: m.age, kind: EventKind::Severed, site: *s, magnitude: 0.001 }, 290.0);
            if i + 1 == 64 {
                at64 = m.render(budget).site.len();
            }
            if i + 1 == 65 {
                let sk = m.render(budget);
                at65 = sk.site.len();
                back = sites[..65].iter().filter(|s| sk.site.contains(s)).count();
            }
        }
        println!("{:<12} {:>8} {:>12} {:>12} {:>14}", program.name(), base.site.len(), at64, at65, back);
    }
    println!("\n(`back at 65` counts severed sites present again. 65 means the log");
    println!(" never suppressed any of them; 33 means the cap dropped the first 32.)");
}

// ---------------------------------------------------------------------------
// §5.8 — does every program book the embodied energy of what it builds?
// ---------------------------------------------------------------------------

fn embodied_energy_booked() {
    rule("§5.8  embodied energy booked at construction");
    println!("{:<12} {:>12} {:>14} {:>14} {:>10}", "program", "built (kg)", "stored (J)", "J/kg", "booked?");
    for program in PROGRAMS {
        let m = matured(program, 0xB00C);
        let stored = m.stored_energy();
        let per_kg = if m.built > 0.0 { stored / m.built } else { 0.0 };
        println!(
            "{:<12} {:>12.4e} {:>14.4e} {:>14.4e} {:>10}",
            program.name(), m.built, stored, per_kg,
            if stored > 0.0 { "yes" } else { "NO" }
        );
    }
    println!("\n(Terrain books zero by design — `energy_density()` returns 0.0 for it,");
    println!(" because rock was not manufactured. That is the contrast §5.8 relies on:");
    println!(" sand carries no embodied energy and a built thing does.)");
}

// ---------------------------------------------------------------------------
// §5.8 — does the embodied-energy check have the digits?
// ---------------------------------------------------------------------------

fn embodied_energy_digits() {
    rule("§5.8  one member's embodied energy against the whole structure");
    println!("{:<12} {:>7} {:>13} {:>13} {:>13} {:>9}", "program", "parts", "whole (J)", "median part", "smallest part", "digits");
    for program in PROGRAMS {
        let m = matured(program, 0xD161);
        let sk = m.render(512);
        let n = sk.site.len().max(1);
        let whole = m.stored_energy();
        let total_mass: f64 = sk.mass.iter().sum();
        if whole <= 0.0 || total_mass <= 0.0 {
            println!("{:<12} {:>7} {:>13.4e} {:>13} {:>13} {:>9}", program.name(), n, whole, "-", "-", "no energy");
            continue;
        }
        let mut shares: Vec<f64> = sk.mass.iter().map(|m0| whole * m0 / total_mass).collect();
        shares.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let smallest = shares[0];
        let median = shares[shares.len() / 2];
        // The worst case is detecting the smallest member's removal.
        let digits = 15.95 - (whole / smallest.max(f64::MIN_POSITIVE)).log10();
        println!(
            "{:<12} {:>7} {:>13.4e} {:>13.4e} {:>13.4e} {:>9.1}",
            program.name(), n, whole, median, smallest, digits
        );
    }
    println!("\n(`digits` is what f64's 15.95 leaves after differencing a whole");
    println!(" structure to find its *smallest* member — the worst case. Below");
    println!(" about 6, the check §5.8 depends on is not reliable.)");
}

// ---------------------------------------------------------------------------
// D11 — can an undesigned genome make something coherent?
// ---------------------------------------------------------------------------

fn undesigned_genomes() {
    rule("D11  a hundred undesigned genomes per program");
    println!("{:<12} {:>10} {:>10} {:>10} {:>12}", "program", "coherent", "degenerate", "impossible", "median parts");
    for program in PROGRAMS {
        let mut st = Stream::at(0x6E0E, 0x99, 0, Purpose::Structure);
        let (mut ok, mut degen, mut bad) = (0, 0, 0);
        let mut parts = Vec::new();
        for _ in 0..100 {
            let mut m = matured(program, 0x6E0E);
            for g in m.genome.iter_mut() {
                *g = st.uniform() as f32;
            }
            let sk = m.render(512);
            let finite = sk.pos.iter().all(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite())
                && sk.mass.iter().all(|x| x.is_finite() && *x >= 0.0)
                && sk.radius.iter().all(|x| x.is_finite() && *x >= 0.0);
            if !finite {
                bad += 1;
            } else if sk.site.len() < 2 {
                degen += 1;
            } else {
                ok += 1;
                parts.push(sk.site.len());
            }
        }
        parts.sort_unstable();
        let median = parts.get(parts.len() / 2).copied().unwrap_or(0);
        println!("{:<12} {:>10} {:>10} {:>10} {:>12}", program.name(), ok, degen, bad, median);
    }
    println!("\n(`impossible` is a non-finite skeleton — the failure D11 says must not");
    println!(" happen. `degenerate` is under two parts: ugly, which is allowed.)");
}

// ---------------------------------------------------------------------------
// D12 — does a grown shape respond to the conditions it grew in?
// ---------------------------------------------------------------------------

fn shape_versus_conditions() {
    rule("D12  same genome, three light fields, normalised to one size");
    let fields = [("dim 40 W/m2", 40.0), ("open 400 W/m2", 400.0), ("bright 900 W/m2", 900.0)];
    let mut grown = Vec::new();
    for (name, light) in fields {
        let mut m = Morphology::new(Program::Tree, 0xF0BE, 0x1234, 0);
        let mut env = Environment::default();
        env.light_flux = light;
        env.reservoir_mass = 1.0e6;
        for _ in 0..400 {
            m.advance(YEAR * 0.25, &env);
        }
        grown.push((name, m));
    }
    println!("{:<16} {:>12} {:>12} {:>12}", "light", "age (yr)", "built (kg)", "extent (m)");
    for (name, m) in &grown {
        println!("{:<16} {:>12.1} {:>12.4e} {:>12.4e}", name, m.age / YEAR, m.built, m.extent());
    }

    // Now hold age and built equal across all three and re-render. Anything that
    // still differs is shape carrying information about the conditions.
    let (age, built) = (grown[1].1.age, grown[1].1.built);
    let mut skeletons = Vec::new();
    for (name, m) in grown.iter() {
        let mut c = m.clone();
        c.age = age;
        c.built = built;
        skeletons.push((*name, c.render(256)));
    }
    let base = &skeletons[0].1;
    println!("\nat equal age and mass:");
    for (name, sk) in &skeletons {
        let same = sk.site.len() == base.site.len()
            && sk.pos.iter().zip(base.pos.iter()).all(|(a, b)| a.x == b.x && a.y == b.y && a.z == b.z);
        println!("  {:<16} {} parts, identical to the first: {}", name, sk.site.len(), same);
    }
    println!("\n(D12's claim is that this is impossible today: `render` is pure in");
    println!(" (genome, age, built, progress, events) and never sees Environment,");
    println!(" so conditions reach shape only through how much mass they grew.)");
}

// ---------------------------------------------------------------------------
// D12 — is growth Markovian in its conditions?
// ---------------------------------------------------------------------------

fn growth_markovian() {
    rule("D12  two condition histories, equal integrals, opposite order");
    let steps = 200;
    let (lo, hi) = (80.0, 720.0); // mean 400, same as the flat case
    let orders: [(&str, fn(usize, usize) -> f64); 3] = [
        ("flat 400", |_, _| 400.0),
        ("rising", |i, n| 80.0 + (720.0 - 80.0) * (i as f64) / (n as f64 - 1.0)),
        ("falling", |i, n| 720.0 - (720.0 - 80.0) * (i as f64) / (n as f64 - 1.0)),
    ];
    println!("light integral is identical in all three (mean {:.0} W/m2)\n", (lo + hi) / 2.0);
    println!("{:<12} {:>12} {:>12} {:>12} {:>10}", "history", "built (kg)", "extent (m)", "stored (J)", "parts");
    let mut first: Option<f64> = None;
    for (name, f) in orders {
        let mut m = Morphology::new(Program::Tree, 0x3A12, 0x1234, 0);
        let mut env = Environment::default();
        env.reservoir_mass = 1.0e6;
        for i in 0..steps {
            env.light_flux = f(i, steps);
            m.advance(YEAR * 0.5, &env);
        }
        let sk = m.render(256);
        println!(
            "{:<12} {:>12.5e} {:>12.5e} {:>12.5e} {:>10}",
            name, m.built, m.extent(), m.stored_energy(), sk.site.len()
        );
        if first.is_none() {
            first = Some(m.built);
        } else if let Some(b) = first {
            let rel = ((m.built - b) / b).abs();
            println!("{:<12} {:>12} relative difference in mass against flat: {:.3e}", "", "", rel);
        }
    }
    println!("\n(If ordering changes the outcome, the response is a function of a");
    println!(" history and D12's per-species surface does not tabulate.)");
}

// ---------------------------------------------------------------------------
// D5 — does a swinging limb leave the small-displacement regime?
// ---------------------------------------------------------------------------

fn limb_chord_rotation() {
    use phys::solvers::dynamics::Dynamics;
    use phys::solvers::frame::{Dof, Framework, Member};
    use phys::topology::Material;

    rule("D5  a two-segment limb swung through a stride");

    // Hip, knee, foot: two 0.4 m segments, 30 mm radius, hip built in.
    let build = || Framework {
        joints: vec![v3(0.0, 0.0, 0.0), v3(0.0, 0.0, -0.4), v3(0.0, 0.0, -0.8)],
        members: vec![
            Member { a: 0, b: 1, radius: 0.03, truss: false, integrity: 1.0 },
            Member { a: 1, b: 2, radius: 0.03, truss: false, integrity: 1.0 },
        ],
        fixed: vec![true, false, false],
        material: Material::GREEN_WOOD,
        lumped: Vec::new(),
        mass_scale: 0.0,
        stiff_scale: 1.0,
    };

    println!("{:>10} {:>13} {:>12} {:>12} {:>8} {:>10}", "force (N)", "travel (m)", "travel / L", "disp_ratio", "broken", "verdict");
    let limb_length = 0.8;
    for force in [1.0, 5.0, 20.0, 60.0, 150.0, 400.0, 700.0, 1000.0] {
        let mut d = Dynamics::new(build());
        let mut load = vec![Dof::default(); 3];
        load[2].t = v3(force, 0.0, 0.0);
        let mut broken = 0usize;
        let mut diverged = false;
        for _ in 0..400 {
            let r = d.step(&load, 1.0e-3);
            broken += r.broken.len();
            if !r.displacement_ratio.is_finite() || r.displacement_ratio > 10.0 {
                diverged = true;
                break;
            }
        }
        let travel = d.displacement[2].t.norm();
        let ratio = d.displacement_ratio();
        let verdict = if diverged || broken > 0 {
            "FAILED"
        } else if ratio > 0.3 {
            "wrong"
        } else if ratio > 0.1 {
            "past regime"
        } else {
            "linear ok"
        };
        if diverged {
            println!("{:>10.1} {:>13} {:>12} {:>12} {:>8} {:>10}", force, "-", "-", "-", broken, verdict);
        } else {
            println!(
                "{:>10.1} {:>13.5} {:>12.5} {:>12.5} {:>8} {:>10}",
                force, travel, travel / limb_length, ratio, broken, verdict
            );
        }
    }
    println!("\n(A stride swings a leg through most of a radian, so the case that");
    println!(" matters is `travel / L` of order 0.5. A row marked FAILED is the");
    println!(" member rupturing or the solve diverging, not a measurement — the");
    println!(" reading is how far the limb gets before one of those happens.)");
}

// ---------------------------------------------------------------------------
// D1/D10 — how large is a checkpoint of an interactive subtree?
// ---------------------------------------------------------------------------

fn checkpoint_size() {
    rule("D1/D10  snapshot size, pinned against regenerable");

    // Only pinned nodes write their bodies; everything else persists as the
    // address and epoch that regenerate it. An interactive subtree is largely
    // pinned by definition -- anything an actor touched is -- so that is the
    // case a replay checkpoint actually pays for.
    {
        let w = World::new(galaxy(0xC4EC, 1e9), 20.0);
        let live = w.tree.nodes.iter().filter(|n| n.alive).count();
        println!("unrefined world: {} live node(s), {} bytes\n", live, phys::persist::encode(w.view()).len());
    }
    println!("{:>8} {:>10} {:>13} {:>13} {:>13}", "nodes", "bodies", "regen (B)", "pinned (B)", "B/body pinned");
    for count in [512usize, 4096, 16384] {
        let mut sizes = [0usize; 2];
        let mut bodies = 0usize;
        let mut live = 0usize;
        for (which, pin) in [(0usize, false), (1usize, true)] {
            let mut w = World::new(galaxy(0xC4EC, 1e9), 20.0);
            w.tree.nodes[0].spec.count = count;
            let root = w.tree.root;
            let path = w.drill_to(root, 1.0e4, &default_spec);
            for idx in &path {
                w.tree.refine(*idx);
                if pin {
                    w.tree.nodes[idx.get()].pinned = true;
                }
            }
            bodies = w.tree.nodes.iter().filter(|n| n.alive).map(|n| n.bodies.len()).sum();
            live = w.tree.nodes.iter().filter(|n| n.alive).count();
            sizes[which] = phys::persist::encode(w.view()).len();
        }
        let per_body = if bodies > 0 {
            (sizes[1] as f64 - sizes[0] as f64) / bodies as f64
        } else {
            0.0
        };
        println!(
            "{:>8} {:>10} {:>13} {:>13} {:>13.1}",
            live, bodies, sizes[0], sizes[1], per_body
        );
    }
    println!("\n(`regen` stores no bodies at all -- the node's address and epoch");
    println!(" regenerate them bit-for-bit. `pinned` is what a checkpoint of a");
    println!(" region an actor has touched actually costs, and the last column is");
    println!(" the marginal cost of one pinned body.)");
}

// ---------------------------------------------------------------------------
// Probes blocked on machinery the plan has not built
// ---------------------------------------------------------------------------

fn blocked() {
    rule("not runnable yet, and why");
    let rows = [
        ("neighbour query cost", "D3's index does not exist; a prototype is Phase 1 work"),

        ("promoted-sibling conservation", "needs D4's per-frame sync to measure against"),
        ("gait optimisation time", "no gait solver exists; D7 is Phase 4"),

        ("busy-square deviation count", "needs a derived erosion rate (§5.2), which is Phase 2"),
        ("erosion: granite vs wet sand", "same — the expression is what Phase 2 writes"),
        ("edit list summarises losslessly", "needs the deviation transform (§5.7), Phase 2"),
        ("embodied energy from topology", "needs the derivation §5.8 proposes; today it is carried"),
        ("regenerable from summarised history", "needs D12's response surface"),
        ("§3.4 crossover on terrestrial matter", "§6 already defers this to Phase 2"),
    ];
    for (what, why) in rows {
        println!("  {:<38} {}", what, why);
    }
}

fn main() {
    let only = std::env::args().nth(1);
    let run = |name: &str| only.as_deref().map_or(true, |o| o == name);

    if run("substeps") { substeps_per_tier(); }
    if run("events") { event_cap_regrowth(); }
    if run("embodied") { embodied_energy_booked(); }
    if run("digits") { embodied_energy_digits(); }
    if run("genomes") { undesigned_genomes(); }
    if run("shape") { shape_versus_conditions(); }
    if run("markov") { growth_markovian(); }
    if run("limb") { limb_chord_rotation(); }
    if run("checkpoint") { checkpoint_size(); }
    if run("blocked") { blocked(); }
    println!();
}
