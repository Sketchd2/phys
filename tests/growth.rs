//! Structures that have a history: growth, construction, and the books that
//! have to balance for either to be legitimate.

use phys::engine::{galaxy, World};
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::morph::*;
use phys::prolong::*;
use phys::state::*;
use phys::units::*;

fn oak(mass: f64) -> (Aggregate, Morphology) {
    let mut m = Morphology::new(Program::Tree, 0xACE, 0x1234, 0);
    m.built = mass;
    m.age = 40.0 * YEAR;
    let mut agg = Aggregate::neutral(mass, m.extent(), 291.0, Program::Tree.substrate());
    agg.chemical_energy = m.stored_energy();
    (agg, m)
}

/// A generated structure is held to exactly the conservation standard a
/// sampled gas cloud is, because it goes through the same projection.
#[test]
fn structures_conserve_like_everything_else() {
    let mut worst = 0.0f64;
    for program in [Program::Tree, Program::Coral, Program::Tower, Program::Wall] {
        for mass in [1.0, 1e3, 1e6] {
            let mut m = Morphology::new(program, 0xACE, 0x99, 0);
            m.built = mass;
            m.design_mass = mass;
            m.progress = 0.6;
            m.age = 10.0 * YEAR;
            let mut agg =
                Aggregate::neutral(mass, m.extent(), 290.0, program.substrate());
            agg.chemical_energy = m.stored_energy();
            agg.momentum = v3(mass * 0.5, -mass * 0.2, 0.0);
            agg.spin = v3(0.0, 0.0, mass * 1e-2);

            for budget in [16usize, 256, 4096] {
                let (bodies, r) = prolong_structured(&agg, &m, budget, 7, 0x99, 0);
                assert!(!bodies.is_empty(), "{:?} produced no geometry", program);
                let mut back = restrict(&bodies, r.potential);
                back.chemical_energy = agg.chemical_energy;
                back.entropy_exported = agg.entropy_exported;
                back.external_potential = agg.external_potential;
                let scales = Scales::of(&bodies);
                let err = back.conserved().error_against(&agg.conserved(), &scales);
                worst = worst.max(err);
                assert!(
                    err < 1e-9,
                    "{:?} mass {mass} budget {budget}: conservation error {err:.3e}",
                    program
                );
                // The mass has to land in the structure, not near it.
                assert!((back.mass - agg.mass).abs() / agg.mass < 1e-12);
            }
        }
    }
    println!("worst structural conservation error {worst:.3e}");
}

/// The whole point: the same tree comes back, not a different one.
#[test]
fn the_same_tree_comes_back() {
    let (agg, m) = oak(900.0);
    let a = prolong_structured(&agg, &m, 2000, 7, 0x1234, 0).0;
    let b = prolong_structured(&agg, &m, 2000, 7, 0x1234, 0).0;
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(x.pos, y.pos, "the tree regrew differently");
        assert_eq!(x.mass, y.mass);
        assert_eq!(x.radius, y.radius);
    }
    // A different tree in the same forest must genuinely differ.
    let mut other = Morphology::new(Program::Tree, 0xACE, 0x5678, 0);
    other.built = m.built;
    other.age = m.age;
    let c = prolong_structured(&agg, &other, 2000, 7, 0x5678, 0).0;
    let differing = a.iter().zip(&c).filter(|(p, q)| p.pos != q.pos).count();
    assert!(
        differing > a.len() / 2,
        "two trees in the same forest are identical"
    );
}

/// The state that stands in for the structure is tiny — that is the trick.
#[test]
fn developmental_state_is_small() {
    let (agg, m) = oak(900.0);
    let bodies = prolong_structured(&agg, &m, 20_000, 7, 0x1234, 0).0;
    let rendered = bodies.len() * std::mem::size_of::<Body>();
    let state = m.state_bytes();
    println!(
        "{} parts = {} bytes rendered, {} bytes of state ({}x)",
        bodies.len(),
        rendered,
        state,
        rendered / state.max(1)
    );
    assert!(state < 512, "developmental state is {state} bytes");
    assert!(rendered / state.max(1) > 100, "not enough compression to matter");
}

/// Growth is a transaction. It cannot mint free energy or order.
#[test]
fn growth_obeys_both_laws() {
    let mut m = Morphology::new(Program::Tree, 1, 2, 0);
    m.built = 50.0;
    let env = Environment::default();
    let mut total_absorbed = 0.0;
    let mut total_entropy = 0.0;
    for _ in 0..200 {
        let txn = m.advance(YEAR / 12.0, &env);
        txn.validate().expect("growth transaction must balance");
        // First law, stated exactly.
        let inflow = txn.energy_absorbed + txn.energy_released;
        let outflow = txn.energy_stored + txn.heat_released + txn.energy_radiated;
        assert!((inflow - outflow).abs() <= 1e-9 * inflow.max(1e-30));
        // Second law: local entropy may fall, the total may not.
        assert!(txn.total_entropy_change() >= -1e-30, "entropy fell");
        total_absorbed += txn.energy_absorbed;
        total_entropy += txn.total_entropy_change();
    }
    println!(
        "16 yr: {:.1} kg, {:.3e} J absorbed, total dS = +{:.3e} J/K",
        m.built, total_absorbed, total_entropy
    );
    assert!(m.built > 50.0, "the tree did not grow");
    assert!(total_entropy > 0.0, "growth exported no entropy at all");
    // And the local entropy really does go the other way — that is the point.
    let sample = m.advance(YEAR, &env);
    assert!(sample.entropy_local < 0.0, "ordering had no entropy cost");
}

/// A program that stored more energy than it absorbed would be rejected.
#[test]
fn impossible_transactions_are_refused() {
    let bad = Transaction {
        mass_incorporated: 1.0,
        composition: Composition::solar(),
        energy_absorbed: 100.0,
        energy_stored: 150.0,
        heat_released: -50.0,
        energy_radiated: 0.0,
        energy_released: 0.0,
        entropy_local: -1.0,
        entropy_exported: 0.0,
    };
    assert!(bad.validate().is_err(), "a free-energy machine was accepted");

    let refrigerator = Transaction {
        mass_incorporated: 1.0,
        composition: Composition::solar(),
        energy_absorbed: 100.0,
        energy_stored: 100.0,
        heat_released: 0.0,
        energy_radiated: 0.0,
        energy_released: 0.0,
        entropy_local: -1.0,
        entropy_exported: 0.0,
    };
    assert!(
        refrigerator.validate().is_err(),
        "ordering for free was accepted"
    );
}

/// Carrying capacity is not a constant anyone typed in — it falls out of the
/// allometry, because capture scales with area and upkeep with mass.
#[test]
fn growth_saturates_from_allometry() {
    let mut m = Morphology::new(Program::Tree, 3, 4, 0);
    m.built = 1.0;
    let env = Environment::default();
    let mut history = Vec::new();
    for year in 0..600 {
        for _ in 0..12 {
            m.advance(YEAR / 12.0, &env);
        }
        if year % 100 == 0 {
            history.push((year, m.built, m.tree_height()));
        }
    }
    for (y, mass, h) in &history {
        println!("  year {y:>3}: {mass:>10.1} kg, {h:>5.1} m tall");
    }
    let early = history[1].1 - history[0].1;
    let late = history[history.len() - 1].1 - history[history.len() - 2].1;
    assert!(m.built > 100.0, "no meaningful growth: {} kg", m.built);
    assert!(late < early, "growth never slowed: {early} then {late}");
    let h = m.tree_height();
    assert!((5.0..90.0).contains(&h), "implausible tree height {h:.1} m");
}

/// Growth stops in winter, which is what makes tree rings.
#[test]
fn growth_responds_to_conditions() {
    let base = Environment::default();
    let winter = Environment { temperature: 268.0, ..base };
    let drought = Environment { water: 0.05, ..base };
    let shade = Environment { light_flux: 4.0, ..base };

    let grow = |env: &Environment| {
        let mut m = Morphology::new(Program::Tree, 5, 6, 0);
        m.built = 100.0;
        for _ in 0..24 {
            m.advance(YEAR / 12.0, env);
        }
        m.built
    };
    let (b, w, d, s) = (grow(&base), grow(&winter), grow(&drought), grow(&shade));
    println!("2 yr biomass — base {b:.1}, winter {w:.1}, drought {d:.1}, shade {s:.1}");
    assert!(b > w && b > d && b > s, "conditions had no effect");
    assert!(w < 100.0, "a tree frozen for two years should lose mass to upkeep");
}

/// Planned construction: progress-driven, and it finishes.
#[test]
fn construction_completes_and_reports_it() {
    let mut m = Morphology::planned(Program::Tower, 4.0e6, 11, 0x77);
    let env = Environment { labour: 1.0 / (200.0 * 86400.0), ..Default::default() };
    let mut done_at = None;
    for day in 0..500 {
        m.advance(86400.0, &env);
        if m.progress >= 1.0 && done_at.is_none() {
            done_at = Some(day);
        }
    }
    let day = done_at.expect("the tower was never finished");
    println!("tower topped out on day {day}, {:.0} t placed", m.built / 1000.0);
    assert!((150..=260).contains(&day), "finished on day {day}");
    assert!((m.built - 4.0e6).abs() / 4.0e6 < 1e-6, "wrong final mass");
    assert!(
        m.events.iter().any(|e| e.kind == EventKind::Completed),
        "completion was not recorded"
    );
}

/// A half-built tower renders as the design masked by its progress.
#[test]
fn partial_construction_renders_partially() {
    let mut m = Morphology::planned(Program::Tower, 4.0e6, 11, 0x77);
    let mut counts = Vec::new();
    for p in [0.1, 0.35, 0.7, 1.0] {
        m.progress = p;
        m.built = 4.0e6 * p;
        counts.push(m.render(2000).len());
    }
    println!("parts placed at 10/35/70/100% completion: {counts:?}");
    for w in counts.windows(2) {
        assert!(w[1] > w[0], "the tower did not grow: {counts:?}");
    }
}

/// Severing a limb must survive coarsening — the whole point of a structure is
/// that its history is visible.
#[test]
fn damage_persists_through_regeneration() {
    let (agg, mut m) = oak(900.0);
    // A budget large enough that the whole tree fits, so the count is the
    // structure's own size rather than the budget's.
    let before = prolong_structured(&agg, &m, 20_000, 7, 0x1234, 0).0;
    let mass_before: f64 = before.iter().map(|b| b.mass).sum();

    let built_before = m.built;
    let txn = m.record(
        Event { at: m.age, kind: EventKind::Severed, site: 2, magnitude: 0.25 },
        291.0,
    );
    txn.validate().expect("severing must balance");
    let after = prolong_structured(&agg, &m, 20_000, 7, 0x1234, 0).0;

    // The structure lost both geometry and mass...
    assert!(m.built < built_before, "severing removed no mass from the structure");
    assert!(
        txn.energy_released > 0.0,
        "the free energy in the severed limb went nowhere"
    );
    // ...but the node still holds every kilogram: the limb is litter now, not
    // gone. This is the assertion that catches a tree getting heavier when you
    // prune it.
    let mass_after: f64 = after.iter().map(|b| b.mass).sum();
    assert!(
        (mass_after - mass_before).abs() / mass_before < 1e-9,
        "node mass changed when a limb was severed"
    );
    let structural: f64 = after
        .iter()
        .filter(|b| b.radius > 0.0)
        .map(|b| b.mass)
        .sum();
    println!(
        "structure {:.0} -> {:.0} kg, node total unchanged at {:.0} kg",
        mass_before, structural, mass_after
    );
    assert!(structural < mass_before * 0.9, "structure did not get lighter");

    // And it is still deterministic afterwards.
    let again = prolong_structured(&agg, &m, 20_000, 7, 0x1234, 0).0;
    for (x, y) in after.iter().zip(&again) {
        assert_eq!(x.pos, y.pos);
    }
    let struct_parts = after.iter().filter(|b| b.radius > 0.0).count();
    println!(
        "severing one limb: {} structural parts -> {}, plus {} of litter",
        before.len(),
        struct_parts,
        after.len() - struct_parts
    );
}

/// The event log cannot grow without bound, however long someone keeps poking.
#[test]
fn event_log_is_bounded() {
    let (_, mut m) = oak(900.0);
    for i in 0..5000 {
        m.record(
            Event {
                at: i as f64,
                kind: EventKind::Damaged,
                site: i as u32,
                magnitude: 0.001,
            },
            290.0,
        );
    }
    println!("after 5000 events the log holds {} and the state is {} bytes",
        m.events.len(), m.state_bytes());
    assert!(m.events.len() <= 64, "event log grew to {}", m.events.len());
    assert!(m.state_bytes() < 2048);
}

/// End to end in the engine: a planted tree grows while nobody looks at it,
/// and the world's energy changes by exactly what crossed the boundary.
#[test]
fn engine_grows_unobserved_structures() {
    let mut w = World::new(galaxy(0x1111, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let node = w.tree.promote(root, 3, phys::engine::default_spec(Tier::Stellar));
    assert!(!node.is_none());

    // Give it a plausible tree-sized aggregate, then plant.
    {
        let n = &mut w.tree.nodes[node.get()];
        n.agg = Aggregate::neutral(500.0, 8.0, 291.0, Program::Tree.substrate());
    }
    w.plant(node, Program::Tree, Environment::default());
    let start_mass = w.tree.nodes[node.get()].agg.mass;
    // Audited on the non-rest energy: rest mass is nine orders larger and would
    // swamp the entire growth budget in round-off.
    let e0 = w.tree.nodes[node.get()].agg.non_rest_energy();
    let s0 = w.tree.nodes[node.get()].agg.total_entropy();
    let absorbed0 = w.tree.stats.external_energy_absorbed;

    // Grow for a simulated decade without ever materialising the tree.
    for _ in 0..120 {
        w.grow_node(node, YEAR / 12.0);
    }
    let n = &w.tree.nodes[node.get()];
    let e1 = n.agg.non_rest_energy();
    let s1 = n.agg.total_entropy();
    let absorbed = w.tree.stats.external_energy_absorbed - absorbed0;

    println!(
        "10 yr unobserved: {:.1} -> {:.1} kg, height {:.2} m, {} materialisations",
        start_mass,
        n.morphology.as_ref().unwrap().built,
        n.morphology.as_ref().unwrap().tree_height(),
        w.tree.stats.materialisations
    );

    assert_eq!(w.rejected_transactions, 0, "a transaction failed to balance");
    assert!(n.morphology.as_ref().unwrap().built > start_mass, "no growth");
    assert!(n.bodies.is_empty(), "growth materialised the structure");

    // The first law across the node boundary, exactly.
    let delta = e1 - e0;
    let err = (delta - absorbed).abs() / absorbed.abs().max(1e-30);
    println!("dE = {delta:.6e} J, net boundary flux = {absorbed:.6e} J, mismatch {err:.2e}");
    assert!(err < 1e-12, "node energy changed by {delta:.6e} but net flux was {absorbed:.6e}");
    // Mass, composition and baryon number are untouched: growth moves carbon
    // from the node's air into its wood, both of which are inside the node.
    assert_eq!(n.agg.mass, start_mass, "growth changed the node's mass");

    // The second law, across the same decade.
    assert!(s1 > s0, "total entropy did not increase");
    // ...and the local entropy went the other way, which is the interesting bit.
    assert!(
        n.agg.entropy < s0,
        "local entropy did not fall while the structure ordered itself"
    );
}

/// A structure survives the refine/coarsen cycle that everything else does.
#[test]
fn structures_round_trip_through_the_tree() {
    let mut w = World::new(galaxy(0x2222, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let node = w.tree.promote(root, 5, phys::engine::default_spec(Tier::Stellar));
    {
        let n = &mut w.tree.nodes[node.get()];
        n.agg = Aggregate::neutral(1200.0, 10.0, 291.0, Program::Tree.substrate());
        n.spec.count = 1500;
    }
    w.plant(node, Program::Tree, Environment::default());
    for _ in 0..60 {
        w.grow_node(node, YEAR / 12.0);
    }

    let e0 = w.tree.nodes[node.get()].agg.total_energy();
    let chem0 = w.tree.nodes[node.get()].agg.chemical_energy;
    let ent0 = w.tree.nodes[node.get()].agg.entropy;
    let first = w.tree.refine(node).to_vec();
    assert!(!first.is_empty(), "structure did not materialise");
    let err = w.tree.coarsen(node);
    let n = &w.tree.nodes[node.get()];

    println!("{} parts, round-trip error {err:.3e}", first.len());
    assert!(err < 1e-9, "structural round trip error {err:.3e}");
    assert_eq!(n.agg.chemical_energy, chem0, "stored free energy was lost");
    assert_eq!(n.agg.entropy, ent0, "restriction overwrote structural entropy");
    assert!((n.agg.total_energy() - e0).abs() / e0.abs() < 1e-12);

    let second = w.tree.refine(node).to_vec();
    for (a, b) in first.iter().zip(&second) {
        assert_eq!(a.pos, b.pos, "the structure regenerated differently");
    }
}
