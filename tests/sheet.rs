//! Phase 5 — water over a patch of ground as a sheet: the shallow-water
//! equations on the patch's own columns (`shallow::Sheet`), the owner's
//! decision over SPH for a patch of shore.

use phys::math::{v3, Vec3};
use phys::ocean::{Sea, Train};
use phys::shallow::{Floor, SeaAtEdge, Sheet};
use phys::solvers::hydro::Wall;

const G: f64 = 9.81;
const WATER: f64 = 1000.0;

/// A slab of ground `half` across each way, its top at `top`, with columns
/// `dx` wide.
fn slab(half: [f64; 2], top: f64) -> Wall {
    Wall {
        centre: v3(0.0, 0.0, top - 0.5),
        radius: 0.0,
        half: v3(half[0], half[1], 0.5),
        orientation: phys::math::Quat::IDENTITY,
        axis: Vec3::ZERO,
        owner: None,
    }
}

/// A box of ground standing on it.
fn block(centre: Vec3, half: Vec3) -> Wall {
    Wall { centre, radius: 0.0, half, orientation: phys::math::Quat::IDENTITY, axis: Vec3::ZERO, owner: None }
}

/// The sea at `level` above the region's centre, with `train` on it.
fn sea(level: f64, train: Option<Train>) -> SeaAtEdge {
    SeaAtEdge {
        sea: Sea { up: v3(0.0, 0.0, 1.0), level, current: Vec3::ZERO, train, density: WATER, g: G },
        centre: v3(0.0, 0.0, 6.371e6),
        height: 0.0,
    }
}

fn run(sheet: &mut Sheet, sea: &SeaAtEdge, span: f64) -> (f64, phys::shallow::Crossing) {
    let mut t = 0.0;
    let mut total = phys::shallow::Crossing::default();
    while t < span {
        let dt = sheet.stable_step(G).min(span - t);
        let c = sheet.step(dt, G, sea, t);
        total.mass += c.mass;
        total.momentum += c.momentum;
        total.bed += c.bed;
        total.heat += c.heat;
        t += dt;
    }
    (t, total)
}

/// **Still water over uneven ground stays still**, under a still sea at the
/// same level: the reconstruction balances the bed's steps exactly, so a lake
/// at rest is a state of the scheme and not an approximation to one.
#[test]
fn still_water_over_uneven_ground_stays_still() {
    let dx = 0.05;
    let walls = vec![
        slab([1.0, 1.0], 0.0),
        block(v3(0.2, 0.1, 0.05), v3(0.15, 0.3, 0.1)),
        block(v3(-0.4, -0.3, 0.02), v3(0.2, 0.1, 0.07)),
    ];
    let floor = Floor::of(&walls, v3(0.0, 0.0, -G), dx, Vec3::ZERO);
    let mut sheet = Sheet::on(&floor, 0.3, WATER, 1e-3).expect("a floor");
    let m0 = sheet.mass();
    let edge = sea(0.3, None);
    let (t, crossed) = run(&mut sheet, &edge, 5.0);
    let fastest = (0..sheet.depth.len()).map(|k| sheet.velocity(k).norm()).fold(0.0, f64::max);
    println!(
        "  {} columns {dx} m wide over two blocks, {t:.1} s: fastest {fastest:.2e} m/s, mass moved {:.2e} of {m0:.1} kg",
        sheet.nx * sheet.ny,
        (sheet.mass() - m0) / m0
    );
    assert!(fastest < 1e-10, "still water moved at {fastest} m/s");
    assert!(((sheet.mass() - m0) / m0).abs() < 1e-12 && crossed.mass.abs() < 1e-9 * m0);
}

/// **A dam breaks as Ritter's solution says.** Water 0.4 m deep behind a line
/// across a dry flat bed, let go: through the rarefaction the depth is
/// `(2 sqrt(g h0) - x / t)^2 / 9 g`, which is `4/9 h0` at the dam's old place
/// whatever the time. Read where the water is deep — its thin front is held
/// by the bed's drag, which Ritter leaves out. Read along
/// the middle of a bed wide enough that what runs off its sides has not
/// reached the middle by then.
#[test]
fn a_dam_breaks_as_ritter_says() {
    let dx = 0.02;
    let (h0, span) = (0.4, 0.4);
    let walls = vec![slab([3.0, 1.5], 0.0)];
    let floor = Floor::of(&walls, v3(0.0, 0.0, -G), dx, Vec3::ZERO);
    // The sea is below the ground, so the edge is dry land.
    let edge = sea(-1.0, None);
    let mut sheet = Sheet::on(&floor, -1.0, WATER, 1e-9).expect("a floor");
    for k in 0..sheet.depth.len() {
        if sheet.holds(k) && sheet.foot(k).dot(sheet.u) < 0.0 {
            sheet.depth[k] = h0;
        }
    }
    let (along, across) = (sheet.u, sheet.v);
    run(&mut sheet, &edge, span);
    let line: Vec<usize> = (0..sheet.depth.len()).filter(|&k| sheet.holds(k) && sheet.foot(k).dot(across).abs() < 0.5 * dx).collect();
    let middle = line.iter().copied().min_by(|&a, &b| {
        sheet.foot(a).dot(along).abs().total_cmp(&sheet.foot(b).dot(along).abs())
    });
    let _ = middle;
    // Ritter's profile through the rarefaction: `h = (2 c0 - x / t)^2 / 9 g`.
    let c0 = (G * h0).sqrt();
    let depth_at = |x: f64| {
        let k = line.iter().copied().min_by(|&a, &b| (sheet.foot(a).dot(along) - x).abs().total_cmp(&(sheet.foot(b).dot(along) - x).abs())).unwrap();
        (sheet.foot(k).dot(along), sheet.depth[k])
    };
    let mut worst: f64 = 0.0;
    for x in [-0.5 * c0 * span, 0.0, 0.5 * c0 * span] {
        let (at, h) = depth_at(x);
        let ritter = (2.0 * c0 - at / span).powi(2) / (9.0 * G);
        println!("  after {span} s, {at:+.3} m from the dam: {h:.4} m deep against Ritter's {ritter:.4}");
        worst = worst.max((h / ritter - 1.0).abs());
    }
    assert!(worst < 0.05, "off Ritter's profile by {worst}");
}

/// **The sea's wave comes in through the edge**, at the sea's period, and
/// what the sheet sends back goes out: over a flat bed half a metre down, the
/// surface in the middle of the patch rises and falls with the train's period
/// and near its height, and settles into that rather than ringing up.
#[test]
fn the_seas_wave_comes_in_through_the_edge() {
    let dx = 0.05;
    let depth = 0.5;
    let walls = vec![slab([1.0, 1.0], -depth)];
    let floor = Floor::of(&walls, v3(0.0, 0.0, -G), dx, Vec3::ZERO);
    // A long swell: 0.1 m high and 10 m long in deep water, carried into half
    // a metre.
    let k = std::f64::consts::TAU / 10.0;
    let deep = Train { amplitude: 0.05, k, omega: phys::ocean::dispersion(k, 4000.0, G, 0.072, WATER), heading: v3(1.0, 0.0, 0.0), depth: 4000.0, g: G, tension: 0.072, density: WATER };
    let edge = sea(0.0, Some(deep));
    let there = deep.in_depth(depth).unwrap();
    let mut sheet = Sheet::on(&floor, 0.0, WATER, 1e-4).expect("a floor");
    let middle = (0..sheet.depth.len())
        .filter(|&k| sheet.holds(k))
        .min_by(|&a, &b| sheet.foot(a).norm().total_cmp(&sheet.foot(b).norm()))
        .unwrap();
    let mut trace: Vec<(f64, f64)> = Vec::new();
    let mut t = 0.0;
    let span = 6.0 * there.period();
    while t < span {
        let dt = sheet.stable_step(G).min(0.01);
        sheet.step(dt, G, &edge, t);
        t += dt;
        trace.push((t, sheet.depth[middle] - depth));
    }
    let late: Vec<&(f64, f64)> = trace.iter().filter(|p| p.0 > 3.0 * there.period()).collect();
    let hi = late.iter().map(|p| p.1).fold(f64::MIN, f64::max);
    let lo = late.iter().map(|p| p.1).fold(f64::MAX, f64::min);
    // Upward crossings of the mean, for the period.
    let mean = late.iter().map(|p| p.1).sum::<f64>() / late.len() as f64;
    let ups: Vec<f64> = late.windows(2).filter(|w| w[0].1 < mean && w[1].1 >= mean).map(|w| w[1].0).collect();
    let period = if ups.len() > 1 { (ups[ups.len() - 1] - ups[0]) / (ups.len() - 1) as f64 } else { 0.0 };
    let height = hi - lo;
    println!(
        "  the train in {depth} m: {:.4} m high, period {:.3} s; in the middle of the patch {height:.4} m high, period {period:.3} s",
        2.0 * there.amplitude,
        there.period()
    );
    assert!((period / there.period() - 1.0).abs() < 0.03, "the patch does not move at the sea's period: {period}");
    assert!((height / (2.0 * there.amplitude) - 1.0).abs() < 0.2, "the patch's wave is {height} m high");
}
