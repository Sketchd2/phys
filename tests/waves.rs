//! Phase 5 — the sea the wind raises. `docs/PLAY.md` §7's beach is "the waves
//! flow through it", and the owner's decisions for this phase say where waves
//! come from and where they stop: the wind's stress is the log-law's, with the
//! roughness measured off the sea's own height; the sea takes energy at
//! `tau c` until its phase speed reaches the wind's; and it breaks at
//! Michell's steepness in deep water and at McCowan's depth near a shore. The
//! wind itself is stated by the scene until a planet's climate is derived.

use phys::material::substances;
use phys::math::{v3, Vec3};
use phys::ocean::{self, Ocean, Wind};

const RADIUS: f64 = 6.371e6;
const DAY: f64 = 86_164.0;

/// The engine's own water: its rest density and its surface tension.
fn water() -> (f64, f64) {
    let props = substances::water();
    let rho = phys::eos::Condensed::liquid(&props).expect("water is a liquid").rest_density;
    (rho, phys::erode::surface_tension(&props))
}

/// An ocean on its own, `depth` deep, with nothing but what a test does to it.
fn an_ocean(depth: f64) -> Ocean {
    let (rho, sigma) = water();
    Ocean::new(8, RADIUS, depth, 9.81, 1e-3, rho, sigma)
}

/// A wind of `speed` over every cell, eastward, at ten metres.
fn everywhere(o: &Ocean, speed: f64) -> Vec<Wind> {
    o.cells
        .iter()
        .enumerate()
        .map(|(k, c)| {
            let east = v3(0.0, 0.0, 1.0).cross(c.up);
            let east = if east.norm() > 1e-6 { east.unit() } else { v3(1.0, 0.0, 0.0) };
            Wind { cell: k, velocity: east.scale(speed), height: 10.0, air_density: 1.225, area: c.area }
        })
        .collect()
}

/// **A wind raises a sea until the sea outruns it.** Ten metres a second over
/// a deep ocean for four days: the sea grows, and stops where its phase speed
/// at the limiting steepness is the wind's — `2 pi u^2 / 7 g`, 9.16 m — and
/// never passes it.
#[test]
fn a_wind_raises_a_sea_until_it_outruns_it() {
    let mut o = an_ocean(4000.0);
    let u = 10.0;
    let winds = everywhere(&o, u);
    let cell = o.cell_of(v3(1.0, 0.0, 0.0));
    let limit = ocean::saturated_height(u, o.g);
    let (rho, sigma) = (o.density, o.tension);
    println!(
        "  water {rho:.0} kg/m^3, {sigma:.4} N/m; the slowest wave on it goes {:.3} m/s; \
         a {u} m/s wind outruns a sea of {limit:.2} m",
        ocean::slowest_wave(o.g, sigma, rho)
    );
    let mut tallest: f64 = 0.0;
    let mut reached = None;
    let dt = 60.0;
    for step in 0..(4.0 * DAY / dt) as usize {
        o.raise(dt, &winds);
        let h = o.wave_height(cell);
        tallest = tallest.max(h);
        if reached.is_none() && h > 0.9 * limit {
            reached = Some(step as f64 * dt);
        }
    }
    let h = o.wave_height(cell);
    let (c, _) = ocean::phase_speed(h, o.depth, o.g, sigma, rho);
    println!(
        "  after four days {h:.3} m, phase speed {c:.3} m/s; 90% of the limit at {:.1} h",
        reached.map(|t| t / 3600.0).unwrap_or(f64::NAN)
    );
    assert!(reached.is_some(), "the sea never came near outrunning the wind");
    assert!(tallest <= limit * 1.001, "the sea grew past the wind: {tallest:.3} m against {limit:.3} m");
    assert!((h - limit).abs() < 0.02 * limit, "the sea stopped short of the wind: {h:.3} m");
}

/// A wind slower than the slowest wave raises nothing, which is surface
/// tension's doing and not a threshold anybody stated.
#[test]
fn a_breath_of_air_raises_no_sea() {
    let mut o = an_ocean(4000.0);
    let floor = ocean::slowest_wave(o.g, o.tension, o.density);
    let winds = everywhere(&o, 0.9 * floor);
    for _ in 0..1000 {
        o.raise(60.0, &winds);
    }
    let h = o.wave_height(o.cell_of(v3(1.0, 0.0, 0.0)));
    println!("  {:.3} m/s of wind under a {floor:.3} m/s slowest wave: {h:.3e} m of sea", 0.9 * floor);
    assert_eq!(h, 0.0);
}

/// **In shallow water a sea breaks at McCowan's depth.** Two metres of water
/// under twenty metres a second of wind: the sea stops at 0.78 of the depth,
/// and what it would have held beyond that is reported as broken.
#[test]
fn a_shallow_sea_breaks_at_its_depth() {
    let mut o = an_ocean(2.0);
    let winds = everywhere(&o, 20.0);
    let mut broken = 0.0;
    for _ in 0..(DAY / 60.0) as usize {
        broken += o.raise(60.0, &winds).iter().map(|b| b.broken).sum::<f64>();
    }
    let h = o.wave_height(o.cell_of(v3(1.0, 0.0, 0.0)));
    println!("  2 m of water, 20 m/s of wind: {h:.4} m of sea, {broken:.3e} J broken");
    assert!((h - ocean::MCCOWAN * 2.0).abs() < 1e-9, "not at the depth limit: {h} m");
    assert!(broken > 0.0);
}

/// **Swell outlives the wind that raised it.** A sea raised over one cell and
/// left alone travels the way it was heading, at its group speed, and none of
/// its energy is made or lost on the way.
#[test]
fn swell_travels_at_its_group_speed_and_keeps_its_energy() {
    let mut o = an_ocean(4000.0);
    let start = o.cell_of(v3(1.0, 0.0, 0.0));
    let east = v3(0.0, 1.0, 0.0);
    let h = 3.0;
    o.cells[start].waves = ocean::energy_of_height(h, o.density, o.g);
    o.cells[start].heading = east;
    let total = |o: &Ocean| o.cells.iter().map(|c| c.waves * c.area).sum::<f64>();
    let centroid = |o: &Ocean| {
        let mut m = Vec3::ZERO;
        for c in &o.cells {
            m += c.up.scale(c.waves * c.area);
        }
        m.unit()
    };
    let e0 = total(&o);
    // Measured from the centre of the cell it started in, which is not the
    // point that picked it.
    let from = o.cells[start].up;
    let (_, group) = ocean::phase_speed(h, o.depth, o.g, o.tension, o.density);
    let span = 2.0 * DAY;
    let mut t = 0.0;
    while t < span {
        let step = o.wave_step().min(span - t);
        o.propagate(step);
        t += step;
    }
    let moved = centroid(&o).dot(from).clamp(-1.0, 1.0).acos() * RADIUS;
    let heading = (centroid(&o) - from).dot(east);
    let e1 = total(&o);
    println!(
        "  a {h} m sea left alone for two days: its centre moved {:.0} km, {:.0} km at its group speed \
         of {group:.2} m/s; energy {:.2e} of itself",
        moved / 1e3,
        group * span / 1e3,
        (e1 - e0) / e0
    );
    assert!(((e1 - e0) / e0).abs() < 1e-12, "swell made or lost energy");
    assert!(heading > 0.0, "the swell did not go the way it was heading");
    assert!((moved / (group * span) - 1.0).abs() < 0.15, "swell did not travel at its group speed");
}
