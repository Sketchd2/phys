//! A guided tour of the engine: galaxy to nucleus and back, with the
//! consistency checks running the whole way.

use phys::budget::FrameBudget;
use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::observe::*;
use phys::state::*;
use phys::tree::Residency;
use phys::units::*;

fn rule(title: &str) {
    println!("\n\x1b[1m{}\x1b[0m", title);
    println!("{}", "-".repeat(title.len().max(60)));
}

fn main() {
    let seed = 0x5EED_1234_ABCD_0001;
    let stars = 1e9;

    rule("1. The problem");
    let baryons = stars * 0.8 * M_SUN / M_PROTON;
    let budget = FrameBudget::ups(20.0);
    let affordable = budget.sim_budget_us() * phys::budget::cost::GPU_SPEEDUP
        / phys::budget::cost::GRAVITY_STEP_US;
    println!("  galaxy of {:.0e} stars contains       {:.3e} baryons", stars, baryons);
    println!("  RTX 2060 can step per 50 ms frame     {:.3e} bodies", affordable);
    println!("  ratio the architecture must supply    {:.2e}", baryons / affordable);
    println!("\n  Nothing is going to close that gap. The engine's job is not to");
    println!("  simulate the galaxy, but to be indistinguishable from one.");

    let mut w = World::new(galaxy(seed, stars), 20.0);
    // The cost model is told the truth about what it is running on. This demo
    // executes the CPU reference implementation, so it budgets for one Ryzen
    // core; the RTX 2060 projection is reported separately rather than assumed.
    w.gpu = false;
    let root = w.tree.root;
    w.tree.nodes[root.get()].residency = Residency::Observed;

    let obs = w.add_observer(Observer {
        anchor: root,
        offset: v3(8.0 * KPC, 0.0, 0.0),
        look: v3(-1.0, 0.0, 0.0),
        field: 3.14,
        angular_resolution: 1e-6,
        integration_time: 3600.0,
        horizon: 1e4 * YEAR,
        priority: 1.0,
        ..Default::default()
    });
    w.gate = phys::causal::CausalGate::new(1e4 * YEAR);

    rule("2. The scale ladder");
    println!("  {:<11} {:>12} {:>13} {:>14} {:>13}", "tier", "length", "timestep", "light cross", "solver");
    for t in Tier::ALL {
        println!(
            "  {:<11} {:>12} {:>13} {:>14} {:>13?}",
            t.name(),
            fmt_length(t.length()),
            fmt_time(t.dt()),
            fmt_time(t.light_crossing()),
            phys::solvers::for_tier(t)
        );
    }

    rule("3. Descending: galaxy to nucleus");
    println!("  Each step materialises one node and promotes its largest child.");
    println!("  Only the chain is built - never the siblings' interiors.\n");
    println!("  {:<11} {:>14} {:>12} {:>9} {:>10} {:>11}", "tier", "mass", "radius", "children", "T (K)", "cons. err");

    let path = w.drill(root, Tier::Nuclear, &default_spec);
    for &idx in &path {
        let n = &w.tree.nodes[idx.get()];
        println!(
            "  {:<11} {:>14} {:>12} {:>9} {:>10.3e} {:>11.2e}",
            n.tier.name(),
            fmt_mass(n.agg.mass),
            fmt_length(n.agg.radius),
            n.bodies.len(),
            n.agg.temperature,
            n.last_report.conservation_error
        );
    }
    let leaf = *path.last().unwrap();
    println!("\n  depth reached: {} tiers, {} live nodes, {} materialised bodies",
        path.len(), w.tree.live_count(), w.tree.materialised_bodies());
    println!("  memory for the whole chain: {:.2} MB", w.tree.detail_bytes() as f64 / 1e6);
    println!("  a galaxy's worth of baryons would need: {:.2e} MB",
        baryons * std::mem::size_of::<Body>() as f64 / 1e6);

    rule("4. Regeneration is exact");
    println!("  Coarsen a node, throw the detail away, rebuild it. The engine");
    println!("  claims the rebuild is bit-identical - so check every bit.\n");
    // Deliberately on a *sibling* branch: coarsening a node releases everything
    // beneath it, and the drilled chain is still in use.
    let sibling = w.tree.promote(root, 7, default_spec(Tier::Galactic));
    w.tree.refine(sibling);
    let target = sibling;
    let before: Vec<Body> = w.tree.nodes[target.get()].bodies.clone();
    let agg_before = w.tree.nodes[target.get()].agg;
    let err = w.tree.coarsen(target);
    let agg_mid = w.tree.nodes[target.get()].agg;
    let after: Vec<Body> = w.tree.refine(target).to_vec();
    let _ = NodeIdx::NONE;
    let identical = before.len() == after.len()
        && before.iter().zip(&after).all(|(a, b)| {
            a.pos == b.pos && a.vel == b.vel && a.mass == b.mass && a.charge == b.charge
                && a.spin == b.spin && a.temperature == b.temperature
        });
    let worst = before
        .iter()
        .zip(&after)
        .map(|(a, b)| {
            let dp = (a.pos - b.pos).norm() / a.pos.norm().max(1e-300);
            let dv = (a.vel - b.vel).norm() / a.vel.norm().max(1e-300);
            dp.max(dv)
        })
        .fold(0.0f64, f64::max);
    println!("  bodies before / after            {} / {}", before.len(), after.len());
    println!("  every field bit-identical        {}", if identical { "yes" } else { "NO" });
    println!("  worst per-body deviation         {:.3e}", worst);
    println!("  conservation error, round trip   {:.3e}", err);
    println!("  mass    {:.17e} -> {:.17e}", agg_before.mass, agg_mid.mass);
    println!("  energy  {:.17e} -> {:.17e}", agg_before.total_energy(), agg_mid.total_energy());
    println!("\n  The coarse state was left untouched because the fine detail");
    println!("  agreed with it to {:.0e}, so there was nothing new to record.", err);

    rule("5. Observation is retarded");
    let sightings = w.look(obs, Instrument::Imager);
    println!("  Nothing is seen as it is now; everything is seen as it was.\n");
    println!("  {:<20} {:>14} {:>16} {:>12}", "target", "distance", "light delay", "doppler");
    for s in sightings.iter().take(6) {
        let n = &w.tree.nodes[s.node.get()];
        println!(
            "  {:<20} {:>14} {:>16} {:>12.9}",
            format!("{} {}", n.tier.name(), s.key.short()),
            fmt_length(s.view.distance),
            fmt_time(s.view.delay),
            s.doppler
        );
    }

    rule("6. Measurement commits");
    println!("  The first measurement of an undetermined quantity samples it.");
    println!("  Every later query returns the same value - forever.\n");
    for i in 0..3 {
        let r = w.measure(leaf, Instrument::Thermometer, Quantity::DecayTime);
        let fact = w.ledger.peek(w.tree.nodes[leaf.get()].key, Quantity::DecayTime);
        if let (Some(Reading::Temperature { kelvin, uncertainty }), Some(f)) = (r, fact) {
            println!(
                "  query {}: T = {:.6e} +/- {:.2e} K   committed decay time = {:.9e} s (seq {})",
                i + 1, kelvin, uncertainty, f.value, f.sequence
            );
        }
    }
    println!("\n  ledger holds {} facts in {:.2} kB - the only part of the world",
        w.ledger.len(), w.ledger.bytes() as f64 / 1e3);
    println!("  that costs permanent memory.");

    rule("7. Interaction propagates at c");
    let d = w.tree.separation(root, v3(8.0 * KPC, 0.0, 0.0), leaf, v3(0.0, 0.0, 0.0)).value.norm();
    println!("  Firing a 10^40 J pulse at a target {} away.", fmt_length(d));
    w.interact(Interaction::Deposit { target: leaf, joules: 1e40, radius: 1.0 });
    println!("  influences in flight: {}", w.mailbox.pending());
    if let Some(t) = w.mailbox.next_arrival() {
        println!("  arrives at t = {} (delay {})", fmt_time(t), fmt_time(t - w.time));
    }
    println!("  It cannot arrive sooner, and the scheduler will not advance the");
    println!("  target past that instant until it does.");

    rule("8. Frames");
    // A second observer, sitting inside the deepest node with a femtometre
    // aperture. This is the expensive case: it forces the whole chain from the
    // galaxy down to the nucleus to stay resolved, every frame.
    w.add_observer(Observer {
        anchor: leaf,
        offset: v3(0.0, 0.0, 0.0),
        look: v3(0.0, 0.0, 1.0),
        field: 3.14,
        angular_resolution: 1e-3,
        integration_time: 1e-21,
        horizon: 1e-18,
        priority: 50.0,
        ..Default::default()
    });
    println!("  Fixed 50 ms budget. What varies is how much simulated time");
    println!("  passes and how much of the world is resolved, never frame rate.\n");
    println!("  {:>5} {:>9} {:>8} {:>8} {:>7} {:>7} {:>11} {:>7}",
        "frame", "wall (ms)", "solved", "coasted", "defer", "bodies", "sim step", "late");
    for f in 0..8 {
        let plan = w.step_frame(50_000.0);
        println!(
            "  {:>5} {:>9.2} {:>8} {:>8} {:>7} {:>7} {:>11} {:>7.0e}",
            f,
            w.stats.last_frame_us / 1e3,
            plan.accepted.len(),
            w.stats.coasted,
            plan.deferred,
            w.stats.materialised_bodies,
            fmt_time(w.frame_dt()),
            w.stats.worst_lateness
        );
    }
    println!("\n  Coasted nodes were carried to the same instant in closed form:");
    println!("  exact, one add each, and nearly the whole world nearly every frame.");
    println!("  Lateness is how many of its own characteristic times the stalest");
    println!("  node has gone unsolved - one number for a nucleus and a galaxy.");
    println!("\n  cost-model calibration after 8 frames: {:.3}x (learned from measurement)",
        w.budget.calibration());
    println!("  measured frame time (EMA):            {:.1} ms", w.budget.measured_us / 1e3);
    println!("  same work on the GPU path (x{:.0}):      {:.1} ms",
        phys::budget::cost::GPU_SPEEDUP,
        w.budget.measured_us / phys::budget::cost::GPU_SPEEDUP / 1e3);
    println!("\n  Deferred work is not lost work - it is detail the frame chose");
    println!("  not to buy. The frame rate is the invariant; fidelity is not.");

    rule("9. Structures: matter with a history");
    println!("  A tree is not a sample from a distribution - it is the specific");
    println!("  record of how it grew. So the generator changes, and the state");
    println!("  it is generated from is a few dozen bytes.\n");

    // Plant a tree on a node of its own and let it grow, unobserved.
    let forest = w.tree.promote(root, 11, default_spec(Tier::Stellar));
    {
        let n = &mut w.tree.nodes[forest.get()];
        n.agg = Aggregate::neutral(2.0, 0.4, 291.0, phys::morph::Program::Tree.substrate());
        n.spec.count = 6000;
    }
    w.plant(forest, phys::morph::Program::Tree, phys::morph::Environment::default());

    let entropy0 = w.tree.nodes[forest.get()].agg.entropy;
    println!("  {:>6} {:>12} {:>10} {:>14} {:>14}", "year", "biomass", "height", "absorbed (J)", "total dS (J/K)");
    let mut absorbed = 0.0;
    let mut entropy = 0.0;
    for year in 0..=160 {
        if year > 0 {
            for _ in 0..12 {
                if let Some(t) = w.grow_node(forest, YEAR / 12.0) {
                    absorbed += t.net_boundary_flux();
                    entropy += t.total_entropy_change();
                }
            }
        }
        if year % 40 == 0 {
            let m = w.tree.nodes[forest.get()].morphology.as_ref().unwrap();
            println!(
                "  {:>6} {:>12} {:>9.1} m {:>14.4e} {:>14.4e}",
                year, fmt_mass(m.built), m.tree_height(), absorbed, entropy
            );
        }
    }
    println!("\n  It never materialised once - {} bodies exist for it.",
        w.tree.nodes[forest.get()].bodies.len());

    let bodies = w.tree.refine(forest).to_vec();
    let m = w.tree.nodes[forest.get()].morphology.as_ref().unwrap();
    println!("  Rendered on demand:      {} parts, {:.2} MB",
        bodies.len(), bodies.len() as f64 * std::mem::size_of::<Body>() as f64 / 1e6);
    println!("  Developmental state:     {} bytes  ({}x smaller)",
        m.state_bytes(),
        (bodies.len() * std::mem::size_of::<Body>()) / m.state_bytes().max(1));
    println!("  Free energy stored:      {:.4e} J", w.tree.nodes[forest.get()].agg.chemical_energy);
    let agg_now = w.tree.nodes[forest.get()].agg;
    println!("  Local entropy:           {:+.4e} J/K   (it ordered itself)",
        agg_now.entropy - entropy0);
    println!("  Exported to surroundings:{:+.4e} J/K   (so the total rose)",
        agg_now.entropy_exported);
    println!("\n  Growth ran on the aggregate: {} steps, no fine structure touched.",
        w.tree.stats.growth_steps);
    println!("  Transactions refused for breaking the books: {}", w.rejected_transactions);

    rule("10. Where the memory went");
    let procedural = w.tree.detail_bytes();
    let persisted: usize = w.tree.persisted.values().map(|v| v.len() * std::mem::size_of::<Body>()).sum();
    println!("  materialised, regenerable   {:>10.2} MB", (procedural - persisted) as f64 / 1e6);
    println!("  pinned, must be stored      {:>10.2} MB", persisted as f64 / 1e6);
    println!("  ledger of committed facts   {:>10.2} kB", w.ledger.bytes() as f64 / 1e3);
    println!("  nodes                       {:>10}", w.tree.live_count());
    println!("\n  {}", w.summary());

    rule("11. Invariants");
    println!("  worst conservation error across every scale transition:  {:.3e}", w.tree.stats.worst_conservation_error);
    println!("  worst causality violation between any two nodes:        {:.3e} s", w.check_causality());
    println!("  materialisations {}, coarsenings {}, bodies created {}, discarded {}",
        w.tree.stats.materialisations, w.tree.stats.coarsenings,
        w.tree.stats.bodies_created, w.tree.stats.bodies_discarded);
    let c = w.conserved();
    println!("  live totals: E = {:.6e} J   |P| = {:.6e}   B = {:.6e}", c.energy, c.momentum.norm(), c.baryon);
    println!();
}
