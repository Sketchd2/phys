//! Nothing outruns light — asserted, not assumed.

use phys::causal::*;
use phys::coords::*;
use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::observe::*;
use phys::units::*;

/// An influence cannot be delivered before its light has had time to arrive.
#[test]
fn influences_never_arrive_early() {
    let mut mb = Mailbox::new();
    let distances = [1.0, 1e6, AU, PARSEC, 8.0 * KPC];
    for (i, d) in distances.iter().enumerate() {
        mb.post(NodeIdx(i as u32), 0.0, *d, InfluenceKind::Radiation, 1e30, v3(0.0, 0.0, 0.0));
    }
    for d in distances {
        let just_before = d / C * (1.0 - 1e-12);
        let arrivals = mb.drain_until(just_before);
        for a in &arrivals {
            assert!(
                a.source_distance / C <= just_before,
                "influence from {} m delivered after {} s",
                a.source_distance,
                just_before
            );
        }
    }
    let all = mb.drain_until(1e30);
    assert!(mb.pending() == 0);
    assert!(!all.is_empty() || mb.delivered > 0);
}

/// Delivery order is total and deterministic, so replay is exact.
#[test]
fn delivery_order_is_deterministic() {
    let order = |()| {
        let mut mb = Mailbox::new();
        for i in 0..200u32 {
            // Deliberately many identical arrival times, to exercise tie-breaks.
            mb.post(NodeIdx(i), 0.0, (i % 5) as f64 * AU, InfluenceKind::Blast, 1.0, v3(0.0, 0.0, 0.0));
        }
        mb.drain_until(1e12)
            .iter()
            .map(|i| (i.target.0, i.arrives))
            .collect::<Vec<_>>()
    };
    assert_eq!(order(()), order(()));
}

/// The retarded-time solver must converge, including at relativistic speeds
/// where the naive iteration is closest to failing.
#[test]
fn retarded_time_converges() {
    for beta in [0.0, 0.1, 0.5, 0.9, 0.99] {
        let v = v3(beta * C, 0.0, 0.0);
        let start = v3(-1e12, 0.0, 0.0);
        let src = |t: f64| start + v.scale(t);
        let observer = v3(0.0, 1e12, 0.0);
        let t_obs = 1e4;
        let (t_ret, pos, dist) = retarded_time(observer, t_obs, &src);
        // The defining condition: |x_obs - x_src(t_ret)| = c (t_obs - t_ret).
        let residual = ((observer - pos).norm() - C * (t_obs - t_ret)).abs() / dist.max(1.0);
        assert!(
            residual < 1e-9,
            "beta={beta}: retardation residual {residual:.3e}"
        );
        assert!(t_ret <= t_obs, "retarded time must precede observation");
    }
}

/// A light cone gate admits exactly what light can reach.
#[test]
fn causal_gate_matches_light_travel() {
    let gate = CausalGate::new(1000.0 * YEAR);
    let r = gate.radius();
    assert!(gate.reaches(r * 0.999));
    assert!(!gate.reaches(r * 1.001));
    assert!((gate.urgency(0.0) - 1.0).abs() < 1e-12);
    assert!(gate.urgency(r) < 1e-12);
}

/// Velocity composition never exceeds c, however many frames are nested.
#[test]
fn nested_boosts_stay_subluminal() {
    let mut v = v3(0.0, 0.0, 0.0);
    for _ in 0..200 {
        v = velocity_add(v, v3(0.9 * C, 0.0, 0.0));
        assert!(v.norm() < C, "composed velocity reached {} c", v.norm() / C);
    }
    assert!(v.norm() > 0.999 * C, "should approach c asymptotically");
}

/// The engine's own end-to-end check: no two disjoint regions may be skewed in
/// time by more than the light travel between them.
#[test]
fn engine_preserves_causality() {
    let mut w = World::new(galaxy(0xBEEF, 1e9), 20.0);
    let root = w.tree.root;
    w.add_observer(Observer {
        anchor: root,
        offset: v3(8.0 * KPC, 0.0, 0.0),
        look: v3(-1.0, 0.0, 0.0),
        field: 3.2,
        angular_resolution: 1e-7,
        horizon: 1e4 * YEAR,
        priority: 1.0,
        ..Default::default()
    });
    w.gate = CausalGate::new(1e4 * YEAR);
    let path = w.drill(root, Tier::Continuum, &default_spec);
    assert!(path.len() > 3);
    for _ in 0..12 {
        w.step_frame(50_000.0);
        let violation = w.check_causality();
        assert_eq!(
            violation, 0.0,
            "causality violated by {violation:.3e} s"
        );
    }
}

/// History is interpolated, not fabricated: a query outside the retained window
/// must announce itself.
#[test]
fn history_reports_its_limits() {
    let mut h = History::new(16);
    for i in 0..16 {
        h.push(Snapshot {
            t: i as f64,
            offset: v3(i as f64 * 1000.0, 0.0, 0.0),
            velocity: v3(1000.0, 0.0, 0.0),
            mass: 1.0,
            luminosity: 0.0,
            temperature: 100.0,
        });
    }
    assert!(h.sample(7.5).is_ok(), "inside the window");
    assert!(h.sample(-100.0).is_err(), "before the window must be flagged");
    assert!(h.sample(1e6).is_err(), "far beyond the window must be flagged");
    // Interpolation must be exact for uniform motion.
    let s = h.sample(7.5).unwrap();
    assert!((s.offset.x - 7500.0).abs() < 1e-6, "got {}", s.offset.x);
}

/// An observer sees the past, and the delay is the light travel time.
#[test]
fn observation_is_retarded() {
    let mut w = World::new(galaxy(0xF00D, 1e9), 20.0);
    let root = w.tree.root;
    let obs = w.add_observer(Observer {
        anchor: root,
        offset: v3(8.0 * KPC, 0.0, 0.0),
        look: v3(-1.0, 0.0, 0.0),
        field: 3.2,
        angular_resolution: 1e-7,
        ..Default::default()
    });
    let path = w.drill(root, Tier::Planetary, &default_spec);
    // Pin the chain: without it the survey correctly coarsens away detail the
    // observer's resolution does not demand, and there is nothing left to see.
    for &n in &path {
        w.tree.pin(n);
    }
    for _ in 0..3 {
        w.step_frame(50_000.0);
    }
    let sightings = w.look(obs, Instrument::Imager);
    assert!(!sightings.is_empty(), "observer should see something");
    for s in &sightings {
        let expected = s.view.distance / C;
        let err = (s.view.delay - expected).abs() / expected.max(1e-30);
        assert!(err < 1e-6, "delay {} vs light travel {}", s.view.delay, expected);
        assert!(s.view.t_retarded <= w.time + 1e-9);
    }
}
