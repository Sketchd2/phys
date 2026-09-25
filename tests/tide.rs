//! Phase 5 — the tide. `docs/PLAY.md` §7: "a squiggle drawn in sand is gone by
//! the next tide", and the owner's decision that the tide is dynamic — the
//! ocean's response through the shallow-water equations to the tidal
//! potential of the bodies outside the planet, not `-potential / g` stood
//! still.

use phys::chem::{Mixture, Phase};
use phys::engine::World;
use phys::ids::NodeIdx;
use phys::material::substances;
use phys::math::v3;
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::Tree;
use phys::units::Tier;

const EARTH: f64 = 5.972e24;
const RADIUS: f64 = 6.371e6;
const MOON: f64 = 7.342e22;
const DISTANCE: f64 = 3.844e8;
const DAY: f64 = 86_164.0;

/// An Earth with an ocean and its moon, in one node. `water` is the ocean's
/// mass: 1.4x10^21 kg is the Earth's own.
fn earth_and_moon(seed: u64, water: f64) -> (World, NodeIdx) {
    let spec = SampleSpec::new(2, Profile::Uniform, MassSpectrum::Equal, BodyKind::Planet);
    let pair = Matter::neutral(EARTH + MOON, 4.0e8, 290.0, Composition::crustal());
    let mut w = World::new(Tree::new(seed, pair, Tier::Planetary, spec), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    {
        let b = &mut w.tree.nodes[root.get()].bodies;
        // The barycentre at the origin, the two on a circular orbit about it.
        let (me, mm) = (EARTH, MOON);
        let omega = (phys::units::G * (me + mm) / DISTANCE.powi(3)).sqrt();
        let (re, rm) = (DISTANCE * mm / (me + mm), DISTANCE * me / (me + mm));
        b[0].pos = v3(-re, 0.0, 0.0);
        b[0].vel = v3(0.0, -omega * re, 0.0);
        b[0].mass = me;
        b[0].radius = RADIUS;
        b[1].pos = v3(rm, 0.0, 0.0);
        b[1].vel = v3(0.0, omega * rm, 0.0);
        b[1].mass = mm;
        b[1].radius = 1.737e6;
    }
    let earth = w.tree.promote(root, 0, SampleSpec::new(8, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain));
    let silica = w.substances.intern(substances::silica_arrangement()).unwrap();
    let h2o = w.substances.intern(substances::water_arrangement()).unwrap();
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0 - water / EARTH);
    mix.add(h2o, Phase::Liquid, water / EARTH);
    let composition = mix.composition(&w.substances).0;
    {
        let n = &mut w.tree.nodes[earth.get()];
        let offset = n.motion.offset;
        n.matter = Matter::neutral(EARTH, RADIUS, 290.0, composition);
        n.matter.spin = v3(0.0, 0.0, 2.0 * std::f64::consts::PI / DAY * n.matter.moment_of_inertia());
        n.motion.offset = offset;
        n.sync_spin_rate();
    }
    w.set_mixture(earth, mix);
    w.pace_fixed(300.0);
    (w, earth)
}

/// The strongest period in a record, scanned from `lo` to `hi` seconds.
fn peak_period(trace: &[(f64, f64)], lo: f64, hi: f64) -> (f64, f64) {
    let mean = trace.iter().map(|p| p.1).sum::<f64>() / trace.len() as f64;
    let mut best = (0.0, 0.0);
    let mut period = lo;
    while period <= hi {
        let w = std::f64::consts::TAU / period;
        let (mut c, mut s) = (0.0, 0.0);
        for (t, y) in trace {
            c += (y - mean) * (w * t).cos();
            s += (y - mean) * (w * t).sin();
        }
        let amp = 2.0 * (c * c + s * s).sqrt() / trace.len() as f64;
        if amp > best.1 {
            best = (period, amp);
        }
        period *= 1.002;
    }
    best
}

/// **The tide comes in twice a lunar day, and the water is all still there.**
///
/// A point on the equator is watched for sixty days of an Earth turning under
/// its orbiting moon — long enough for the spectrum to tell the tide from the
/// ocean's own free oscillations beside it, which a week cannot. The strongest period in its elevation is the lunar
/// semidiurnal one — half of the time between the moon passing overhead,
/// 12.42 h — and the ocean's volume does not move beyond the rounding of a
/// scheme that gives every edge's flux to both sides.
///
/// The range beats from day to day, and that is recorded rather than hidden:
/// switching the tide on from rest rings the ocean's own free oscillations,
/// and on a flat seabed nothing damps them quickly — the log-law drag's time
/// scale at these speeds is months. A real ocean's shelves and ridges take a
/// tide's energy in about a day; this one has none yet.
#[test]
fn the_tide_comes_in_twice_a_lunar_day() {
    let (mut w, earth) = earth_and_moon(0x71DE, 1.4e21);
    let mut trace: Vec<(f64, f64)> = Vec::new();
    let point = v3(1.0, 0.0, 0.0);
    for _ in 0..(60.0 * DAY / 300.0) as usize {
        w.step_frame(1_000_000.0);
        let Some(o) = w.tree.nodes[earth.get()].ocean.as_ref() else { continue };
        trace.push((w.time, o.cells[o.cell_of(point)].eta));
    }
    let o = w.tree.nodes[earth.get()].ocean.as_ref().expect("an Earth with water has an ocean");
    let volume = o.depth * 4.0 * std::f64::consts::PI * o.radius * o.radius;
    let lost = o.excess_volume() / volume;
    let late: Vec<(f64, f64)> = trace.iter().copied().filter(|(t, _)| *t > 2.0 * DAY).collect();
    let (period, amplitude) = peak_period(&late, 4.0 * 3600.0, 48.0 * 3600.0);
    let orbit = std::f64::consts::TAU / (phys::units::G * (EARTH + MOON) / DISTANCE.powi(3)).sqrt();
    let lunar_day = 1.0 / (1.0 / DAY - 1.0 / orbit);
    let equilibrium = 1.5 * MOON / EARTH * RADIUS.powi(4) / DISTANCE.powi(3);
    println!(
        "  an ocean {:.0} m deep, {} cells, stepped every {:.0} s, seabed drag {:.2e}",
        o.depth,
        o.cells.len(),
        o.stable_step(),
        o.drag
    );
    println!(
        "  at the equator the strongest period is {:.2} h (half a lunar day {:.2} h), amplitude {:.3} m \
         (equilibrium {:.3} m); volume moved {lost:.1e} of the ocean",
        period / 3600.0,
        lunar_day / 7200.0,
        amplitude,
        equilibrium
    );
    for day in [1, 2, 5, 10, 20, 40, 59] {
        let d: Vec<f64> = trace.iter().filter(|p| p.0 > day as f64 * DAY && p.0 <= (day + 1) as f64 * DAY).map(|p| p.1).collect();
        let range = d.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - d.iter().cloned().fold(f64::INFINITY, f64::min);
        print!(" day {day} {range:.2} m;");
    }
    println!();
    assert!(lost.abs() < 1e-12, "water was made or lost");
    assert!((period - lunar_day / 2.0).abs() < 0.03 * lunar_day / 2.0, "the tide is not semidiurnal: {period} s");
    assert!(amplitude > 0.1 * equilibrium && amplitude < 10.0 * equilibrium);
}
