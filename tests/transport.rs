//! Phase 5 — multi-resolution transport: water crossing between a coarse
//! ocean, which is a description on its planet, and a fine region, where it is
//! drawn as parcels. `docs/PLAY.md` §4.3's fifth piece.
//!
//! What crosses is the sea the ocean knows, as the one wave it knows how to be
//! (`ocean::Train`) — the owner's decision.

use phys::math::{v3, Vec3};
use phys::ocean::{dispersion, Train, BREAKING_STEEPNESS, MCCOWAN};

const G: f64 = 9.81;
const WATER: f64 = 1000.0;
const TENSION: f64 = 0.072;

/// **The train is the sea's own numbers.** Its length is the ocean's limiting
/// steepness, its energy the sea's, and its period what the dispersion
/// relation gives that length in deep water: for a 2 m sea, 14 m long, and
/// `omega^2 = g k` to the capillary term.
#[test]
fn a_sea_is_one_wave_with_the_seas_own_numbers() {
    let h = 2.0;
    let t = Train::of(h, v3(1.0, 0.0, 0.0), 4000.0, G, TENSION, WATER).unwrap();
    let length = std::f64::consts::TAU / t.k;
    let sea = WATER * G * h * h / 16.0;
    let deep = (G * t.k).sqrt();
    println!(
        "  a {h} m sea: {length:.3} m long, period {:.3} s, {:.4} J/m^2 against the sea's {sea:.4}; omega {:.6} against sqrt(g k) {deep:.6}",
        t.period(),
        t.energy(),
        t.omega
    );
    assert!((length - h / BREAKING_STEEPNESS).abs() < 1e-9 * length);
    assert!((t.energy() - sea).abs() < 1e-9 * sea, "the train does not carry the sea's energy");
    assert!((t.omega / deep - 1.0).abs() < 1e-3, "deep water does not disperse as deep water");
}

/// **Into shallower water it keeps its period and its energy flux**, grows,
/// and breaks at McCowan's limit.
#[test]
fn a_wave_coming_ashore_keeps_its_period_and_breaks_at_its_depth() {
    let deep = Train::of(1.0, v3(1.0, 0.0, 0.0), 4000.0, G, TENSION, WATER).unwrap();
    let flux = deep.energy() * deep.group_speed();
    let mut rows = Vec::new();
    for d in [20.0, 5.0, 2.0, 1.0, 0.5] {
        let t = deep.in_depth(d).unwrap();
        rows.push((d, 2.0 * t.amplitude, t.period(), t.energy() * t.group_speed() / flux));
    }
    for (d, h, p, f) in &rows {
        println!("  {d:>4} m deep: height {h:.4} m, period {p:.4} s, energy flux {f:.6} of the deep water's");
    }
    for &(d, h, p, f) in &rows {
        assert!((p / deep.period() - 1.0).abs() < 1e-9, "the period changed at {d} m");
        if h < MCCOWAN * d - 1e-12 {
            assert!((f - 1.0).abs() < 1e-6, "energy flux not conserved at {d} m: {f}");
        } else {
            assert!((h - MCCOWAN * d).abs() < 1e-12, "not broken at McCowan's limit at {d} m");
        }
    }
    // A steep sea breaks before it reaches water shallow enough to grow in:
    // the energy flux is kept, and the shoaling coefficient dips below one in
    // between. A long swell is what grows — 100 m long, reaching 1 m of water.
    let k = std::f64::consts::TAU / 100.0;
    let swell = Train { amplitude: 0.1, k, omega: dispersion(k, 4000.0, G, TENSION, WATER), depth: 4000.0, ..deep };
    let ashore = swell.in_depth(1.0).unwrap();
    // Against the shallow-water limit, `sqrt(c_g deep / sqrt(g d))`.
    let limit = (0.5 * (G / k).sqrt() / (G * 1.0f64).sqrt()).sqrt();
    let grew = ashore.amplitude / swell.amplitude;
    println!("  a 100 m swell, 0.2 m high: {:.4} m high in a metre of water, {grew:.4}x against the shallow limit's {limit:.4}", 2.0 * ashore.amplitude);
    assert!(grew > 1.0 && (grew / limit - 1.0).abs() < 0.02, "a swell did not shoal as shallow water says");
}

/// **The water moves as the surface does**: at the surface its vertical
/// velocity is the surface's own rate of rise, and at the bed it has none.
#[test]
fn the_water_under_a_wave_moves_as_its_surface_does() {
    let t = Train::of(0.8, v3(1.0, 0.0, 0.0), 3.0, G, TENSION, WATER).unwrap().in_depth(3.0).unwrap();
    let up = v3(0.0, 0.0, 1.0);
    let x = v3(1.3, 0.0, 0.0);
    let dt = 1e-6;
    let (mut worst_top, mut worst_bed): (f64, f64) = (0.0, 0.0);
    for i in 0..20 {
        let at = 0.37 * i as f64;
        let rise = (t.surface(x, at + dt) - t.surface(x, at - dt)) / (2.0 * dt);
        worst_top = worst_top.max((t.velocity(x, 0.0, up, at).dot(up) - rise).abs());
        worst_bed = worst_bed.max(t.velocity(x, -t.depth, up, at).dot(up).abs());
    }
    let scale = t.amplitude * t.omega;
    println!("  surface: {:.2e} of a omega off its rise; bed: {:.2e} of a omega through it", worst_top / scale, worst_bed / scale);
    assert!(worst_top < 1e-6 * scale && worst_bed < 1e-12 * scale);
    let _ = (dispersion, Vec3::ZERO);
}
