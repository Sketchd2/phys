//! The user-facing surface: observing, measuring, perturbing, authoring.

use phys::engine::{default_spec, galaxy, World};
use phys::math::v3;
use phys::observe::*;
use phys::units::*;

/// Index of the deepest live node, which is where the interesting physics is.
fn deepest(w: &World) -> phys::ids::NodeIdx {
    let mut best = w.tree.root;
    let mut d = 0;
    for (i, n) in w.tree.nodes.iter().enumerate() {
        if n.alive && n.depth >= d {
            d = n.depth;
            best = phys::ids::NodeIdx(i as u32);
        }
    }
    best
}

fn world() -> (World, usize) {
    let mut w = World::new(galaxy(0xAB1E, 1e9), 20.0);
    let root = w.tree.root;
    let obs = w.add_observer(Observer {
        anchor: root,
        offset: v3(8.0 * KPC, 0.0, 0.0),
        look: v3(-1.0, 0.0, 0.0),
        field: 3.2,
        angular_resolution: 1e-8,
        integration_time: 3600.0,
        horizon: 1e4 * YEAR,
        priority: 1.0,
        ..Default::default()
    });
    w.gate = phys::causal::CausalGate::new(1e4 * YEAR);
    let path = w.drill(root, Tier::Planetary, &default_spec);
    for &n in &path {
        w.tree.pin(n);
    }
    (w, obs)
}

/// A measured quantity is committed: ask again, get the same answer, forever.
#[test]
fn measurement_commits_permanently() {
    let (mut w, _) = world();
    let target = w.tree.nodes.iter().position(|n| n.alive && n.depth == 2).unwrap();
    let target = phys::ids::NodeIdx(target as u32);
    let key = w.tree.nodes[target.get()].key;

    let first = w.measure(target, Instrument::Thermometer, Quantity::DecayTime);
    assert!(first.is_some());
    let fact = w.ledger.peek(key, Quantity::DecayTime).expect("must commit");
    for _ in 0..25 {
        w.measure(target, Instrument::Thermometer, Quantity::DecayTime);
        let again = w.ledger.peek(key, Quantity::DecayTime).unwrap();
        assert_eq!(again.value, fact.value, "committed value changed");
        assert_eq!(again.sequence, fact.sequence, "fact was re-sampled");
    }
    assert_eq!(w.ledger.len(), 1, "one quantity, one fact");
}

/// Measurement disturbs. An interferometer that reports a position without
/// depositing momentum is claiming a free lunch.
#[test]
fn measurement_disturbs() {
    let (mut w, _) = world();
    let target = deepest(&w);
    let before = w.tree.nodes[target.get()].agg.internal_energy;
    let r = w.measure(target, Instrument::Interferometer, Quantity::Position);
    match r {
        Some(Reading::Position { uncertainty, disturbance, .. }) => {
            assert!(uncertainty > 0.0);
            assert!(disturbance > 0.0, "position measurement must cost momentum");
            // dx dp >= hbar/2
            let dp = (2.0 * w.tree.nodes[target.get()].agg.mass * disturbance).sqrt();
            assert!(
                uncertainty * dp >= H_BAR / 2.0 * 0.99,
                "uncertainty principle violated"
            );
        }
        other => panic!("expected a position reading, got {other:?}"),
    }
    let after = w.tree.nodes[target.get()].agg.internal_energy;
    assert!(after >= before, "measurement removed energy");
    assert!(w.tree.nodes[target.get()].pinned, "a measured node must be pinned");

    // At this scale the disturbance is real but utterly negligible — locating a
    // planet-sized body to 10^12 m deposits ~10^-90 J, which is the correct
    // answer and is why nobody notices quantum mechanics on a planet. The
    // effect only becomes significant when the target is small and the
    // precision is high, so check the sharp case directly.
    let mut s = phys::rng::Stream::at(1, 2, 0, phys::rng::Purpose::QuantumMeasure);
    let m = phys::solvers::quantum::measure_position(
        M_PROTON,
        v3(0.0, 0.0, 0.0),
        v3(1.0, 0.0, 0.0),
        1e-15,
        &mut s,
    );
    // Confining a proton to a femtometre costs ~20 MeV. That is not a rounding
    // error, it is why you cannot watch a nucleus without changing it.
    let mev = m.disturbance / MEV;
    println!("locating a proton to 1 fm deposits {mev:.1} MeV");
    assert!(mev > 1.0 && mev < 1000.0, "got {mev:.3} MeV");
}

/// A user impulse cannot take effect before its light arrives.
#[test]
fn interactions_respect_light_delay() {
    let (mut w, _) = world();
    let target = phys::ids::NodeIdx(
        w.tree.nodes.iter().position(|n| n.alive && n.depth == 2).unwrap() as u32,
    );
    let d = w
        .tree
        .separation(w.tree.root, v3(8.0 * KPC, 0.0, 0.0), target, v3(0.0, 0.0, 0.0))
        .value
        .norm();
    let before = w.tree.nodes[target.get()].agg.momentum;
    w.interact(Interaction::Impulse { target, dp: v3(1e40, 0.0, 0.0) });

    let arrival = w.mailbox.next_arrival().expect("influence must be queued");
    let expected = w.time + d / C;
    assert!((arrival - expected).abs() / expected < 1e-9, "arrival {arrival} vs {expected}");

    // Stepping short of the arrival must not change the target.
    for _ in 0..3 {
        w.step_frame(50_000.0);
        if w.time < arrival {
            assert_eq!(
                w.tree.nodes[target.get()].agg.momentum, before,
                "influence arrived early"
            );
        }
    }
}

/// Pinned detail survives coarsening: once a user has touched something, it is
/// no longer derivable and must be stored.
#[test]
fn touched_detail_is_persisted() {
    let (mut w, _) = world();
    let target = phys::ids::NodeIdx(
        w.tree.nodes.iter().position(|n| n.alive && n.depth == 1).unwrap() as u32,
    );
    w.tree.refine(target);
    let key = w.tree.nodes[target.get()].key;
    // Alter the detail, as an interaction would.
    w.tree.nodes[target.get()].bodies[0].vel += v3(1e5, 0.0, 0.0);
    let marker = w.tree.nodes[target.get()].bodies[0].vel;
    w.tree.pin(target);
    w.tree.coarsen(target);
    assert!(w.tree.persisted.contains_key(&key), "pinned detail was thrown away");
    let back = w.tree.refine(target).to_vec();
    assert_eq!(back[0].vel, marker, "persisted detail came back changed");
}

/// Authoring is allowed but audited: the engine reports exactly what the user
/// broke rather than absorbing it silently.
#[test]
fn authoring_is_audited() {
    let (mut w, _) = world();
    let target = phys::ids::NodeIdx(
        w.tree.nodes.iter().position(|n| n.alive && n.depth == 2).unwrap() as u32,
    );
    let before = w.tree.nodes[target.get()].agg.total_energy();
    w.interact(Interaction::Author {
        target,
        property: Property::Temperature,
        value: 1e8,
    });
    let after = w.tree.nodes[target.get()].agg.total_energy();
    assert_eq!(w.audit.len(), 1, "authoring must be recorded");
    let ev = w.audit[0];
    assert!((ev.delta_energy - (after - before)).abs() <= (after - before).abs() * 1e-9);
    assert!(ev.delta_energy != 0.0, "heating to 10^8 K should cost energy");
    println!("user injected {:.3e} J, recorded in the audit log", ev.delta_energy);
}

/// Different instruments give different, physically appropriate answers about
/// the same object.
#[test]
fn instruments_disagree_appropriately() {
    let (mut w, obs) = world();
    for instrument in [
        Instrument::Imager,
        Instrument::Spectrometer,
        Instrument::Thermometer,
        Instrument::ParticleDetector,
        Instrument::MassSpectrometer,
    ] {
        let sightings = w.look(obs, instrument);
        assert!(!sightings.is_empty(), "{instrument:?} saw nothing");
        for s in sightings.iter().take(3) {
            match &s.reading {
                Reading::Spectrum { bins, .. } => assert!(bins.iter().all(|(l, f)| *l > 0.0 && f.is_finite())),
                Reading::Temperature { kelvin, .. } => assert!(*kelvin >= 2.725),
                Reading::Composition { fractions } => {
                    let sum: f64 = fractions.iter().sum();
                    assert!((sum - 1.0).abs() < 1e-9, "composition sums to {sum}");
                }
                Reading::Flux { value, .. } => assert!(value.is_finite() && *value >= 0.0),
                Reading::ParticleCount { expected, .. } => assert!(expected.is_finite()),
                Reading::Unresolved { angular_size, needed } => assert!(angular_size < needed),
                Reading::Position { .. } | Reading::Bulk { .. } => {}
            }
        }
    }
}

/// Zooming in and out repeatedly must not drift. This is the interaction
/// pattern a real user generates constantly.
#[test]
fn repeated_zoom_does_not_drift() {
    let (mut w, _) = world();
    let target = phys::ids::NodeIdx(
        w.tree.nodes.iter().position(|n| n.alive && n.depth == 1).unwrap() as u32,
    );
    w.tree.nodes[target.get()].pinned = false;
    w.tree.nodes[target.get()].residency = phys::tree::Residency::Speculative;
    // One cycle to settle: the first coarsening after a drill legitimately
    // updates the aggregate, because it folds back state the promoted children
    // evolved on their own. What must not drift is everything after that.
    w.tree.coarsen(target);
    let first = w.tree.refine(target).to_vec();
    let e0 = w.tree.nodes[target.get()].agg.total_energy();
    for _ in 0..50 {
        w.tree.coarsen(target);
        w.tree.refine(target);
    }
    let last = w.tree.nodes[target.get()].bodies.clone();
    let e1 = w.tree.nodes[target.get()].agg.total_energy();
    assert_eq!(e0, e1, "50 zoom cycles changed the energy");
    for (a, b) in first.iter().zip(&last) {
        assert_eq!(a.pos, b.pos, "50 zoom cycles moved a body");
        assert_eq!(a.vel, b.vel);
    }
    println!("50 zoom cycles: {} idempotent", w.tree.stats.idempotent_coarsenings);
}
