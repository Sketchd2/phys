//! Two identical trees. One of them is in a hurry.
//!
//! ```sh
//! cargo run --release --bin phys-bubble
//! ```
//!
//! A time bubble is an administrator's tool: encase a node and run its interior
//! faster than the universe around it, so a century of growth or a reaction
//! going to completion can be watched, measured and balanced in an afternoon.
//! It is deliberately unphysical, and everything here is arranged to show that
//! it is unphysical in *exactly* the intended way and no other.
//!
//! The two trees start from the same program in the same world under the same
//! clock. The only difference between them is one `Dilate` command. What should
//! follow is that the bubbled one ages faster and does nothing else differently
//! — same world clock, same position, same physics, and an audit trail saying
//! who did it and a meter saying how much.

use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::morph::{Environment, Program};
use phys::observe::Interaction;
use phys::units::*;

const RATE: f64 = 500.0;
const FRAMES: usize = 60;

fn rule(title: &str) {
    println!("\n\x1b[1m{title}\x1b[0m");
    println!("{}", "─".repeat(title.chars().count()));
}

fn si(v: f64, unit: &str) -> String {
    const P: [(f64, &str); 8] =
        [(1e18, "E"), (1e15, "P"), (1e12, "T"), (1e9, "G"), (1e6, "M"), (1e3, "k"), (1.0, ""), (1e-3, "m")];
    if !v.is_finite() {
        return format!("∞ {unit}");
    }
    let a = v.abs();
    if a >= 1e21 || (a > 0.0 && a < 1e-3) {
        return format!("{v:.3e} {unit}");
    }
    for (s, p) in P {
        if a >= s {
            return format!("{:.2} {p}{unit}", v / s);
        }
    }
    format!("{v:.3e} {unit}")
}

const YEAR: f64 = 31_557_600.0;

/// Ages here span four orders of magnitude between the two trees, so a fixed
/// two decimal places would print one of them as zero.
fn years(seconds: f64) -> String {
    let y = seconds / YEAR;
    if y >= 1.0 {
        format!("{y:.3} yr")
    } else if y > 0.0 {
        format!("{y:.3e} yr")
    } else {
        "0".to_string()
    }
}

/// Age, built mass and stored chemical energy — what a grower has to show for
/// itself.
fn state(w: &World, i: NodeIdx) -> (f64, f64, f64) {
    let n = &w.tree.nodes[i.get()];
    match &n.morphology {
        Some(m) => (m.age, m.built, n.matter.chemical_energy),
        None => (0.0, 0.0, 0.0),
    }
}

fn main() {
    // One world, one clock, two siblings of the same parent so that neither can
    // be said to be sitting somewhere more convenient than the other.
    let mut w = World::new(galaxy(0xB0BB1E, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 512;
    let root = w.tree.root;
    let here = *w.drill(root, Tier::Continuum, &default_spec).last().unwrap();
    w.tree.refine(here);

    let spec = default_spec(w.tree.nodes[here.get()].tier.finer());
    let slow = w.tree.promote(here, 0, spec);
    let fast = w.tree.promote(here, 1, spec);
    for i in [slow, fast] {
        w.plant(i, Program::Tree, Environment::default());
    }
    // Watching a tree grow means a clock measured in days, not the
    // milliseconds a resolved continuum node's dynamics ask for.
    w.pace_fixed(86_400.0);
    w.time_rate = 1.0;

    rule("Two trees");
    println!("  planted as siblings in one {} node, under one clock",
        w.tree.nodes[here.get()].tier.name());
    println!("  the only difference between them is the next line");

    rule("The command");
    // Through the command path, not the method, because that is how an
    // administrator would reach it — and unlike an impulse it applies at once,
    // since a bubble is not an influence crossing the world to get there.
    w.interact(Interaction::Dilate { target: fast, rate: RATE });
    let applied = w.tree.nodes[fast.get()].bubble;
    if applied == 1.0 {
        eprintln!("  the world refused the bubble");
        std::process::exit(1);
    }
    println!("  Dilate {{ target: <tree B>, rate: {RATE} }}  ->  accepted at {applied}x");
    println!("  audit now holds {} entry(s); {} node(s) bubbled",
        w.audit.len(), w.bubbles().len());

    let r = w.time_rate_of(fast);
    println!(
        "\n  tree A runs at {:.9} s/s   (kinematic {:.9}, gravitational {:.9}, bubble {})",
        w.time_rate_of(slow).total(),
        w.time_rate_of(slow).kinematic,
        w.time_rate_of(slow).gravitational,
        w.time_rate_of(slow).bubble
    );
    println!(
        "  tree B runs at {:.9} s/s   (kinematic {:.9}, gravitational {:.9}, bubble {})",
        r.total(),
        r.kinematic,
        r.gravitational,
        r.bubble
    );
    println!(
        "\n  Neither is exactly one even before the bubble: a node seven frames\n  \
         deep in a galaxy is already dilated by its ancestors' orbital speeds\n  \
         and by the halo it sits in. Those numbers were always computed. Until\n  \
         now nothing read them."
    );

    let t0 = w.time;
    for _ in 0..FRAMES {
        w.step_frame(30_000.0);
    }
    let elapsed = w.time - t0;

    rule("After one run of the world");
    let (age_a, built_a, chem_a) = state(&w, slow);
    let (age_b, built_b, chem_b) = state(&w, fast);
    println!(
        "  world clock advanced {} ({}) across {FRAMES} frames",
        si(elapsed, "s"),
        years(elapsed)
    );
    println!(
        "  — less than the {} asked for, because the throttle reports what the \n  \
         machine actually sustained rather than what was requested",
        si(86_400.0 * FRAMES as f64, "s")
    );
    println!("\n  {:<10} {:>16} {:>16} {:>18}", "", "age", "built", "stored energy");
    println!("  {:<10} {:>16} {:>16} {:>18}", "tree A", years(age_a), si(built_a, "g"), si(chem_a, "J"));
    println!("  {:<10} {:>16} {:>16} {:>18}", "tree B", years(age_b), si(built_b, "g"), si(chem_b, "J"));
    if age_a > 0.0 {
        println!("\n  tree B has lived {:.1}x as long as tree A", age_b / age_a);
    }

    rule("And nothing else moved");
    let (pa, pb) = (
        w.tree.nodes[slow.get()].motion.offset,
        w.tree.nodes[fast.get()].motion.offset,
    );
    println!("  world instant is one number: {} for both", si(w.time, "s"));
    println!("  A is at {} from its parent, B at {}", si(pa.norm(), "m"), si(pb.norm(), "m"));
    println!(
        "  node clocks: A {}, B {}  — identical, because a bubble scales what a\n  \
         node *does*, never where it *is*. An object that travelled faster is an\n  \
         object with more velocity, which the engine already models.",
        si(w.tree.nodes[slow.get()].time, "s"),
        si(w.tree.nodes[fast.get()].time, "s")
    );

    rule("What it cost, stated");
    println!(
        "  {} of interior evolution handed out beyond what the clock paid for",
        years(w.stats.bubble_seconds)
    );
    println!(
        "\n  Energy is still conserved inside the bubble: tree B really did burn\n  \
         {} of its own reserves, and every joule is accounted for. What it did\n  \
         *not* receive is {} of sunlight, because only {} passed outside. A\n  \
         bubble does not create energy — it breaks the balance between what a\n  \
         node takes in and what it spends, and the meter above is the size of\n  \
         that break in the one unit that covers every process at once.",
        years(age_b),
        years(age_b),
        years(elapsed)
    );
    println!(
        "\n  worst lateness {:.3e}. A day a frame is far coarser than a resolved\n  \
         continuum node's dynamics can follow, and the engine says so rather\n  \
         than pretending: growth runs on the matter and costs O(1), so it\n  \
         keeps up on nodes whose trajectories cannot. Tree B is additionally\n  \
         {}x as late as its own clock says, which is how the scheduler knows to\n  \
         prioritise it — and why a bubbled node is never thermalised away.",
        w.stats.worst_lateness,
        RATE
    );
}
