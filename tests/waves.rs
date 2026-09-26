//! Phase 5 — the sea the wind raises. `docs/PLAY.md` §7's beach is "the waves
//! flow through it", and the owner's decisions for this phase say where waves
//! come from and where they stop: the wind's stress is the log-law's, with the
//! roughness measured off the sea's own height; the sea takes energy at
//! `tau c` until its phase speed reaches the wind's; and it breaks at
//! Michell's steepness in deep water and at McCowan's depth near a shore. The
//! wind itself is stated by the scene until a planet's climate is derived.

use phys::chem::{Arrangement, Bond, Element, Mixture, Order, Phase};
use phys::engine::World;
use phys::ids::NodeIdx;
use phys::recipe::Recipe;
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::Tree;
use phys::units::Tier;
use phys::material::substances;
use phys::math::{v3, Vec3};
use phys::ocean::{self, Ocean, Wind};

const RADIUS: f64 = 6.371e6;
const EARTH: f64 = 5.972e24;
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

/// Nitrogen, which is what the air is here.
fn nitrogen() -> Arrangement {
    Arrangement::molecule(vec![Element(7), Element(7)], vec![Bond::new(0, 1, Order::Triple)])
}

/// An Earth with an ocean and an atmosphere, turning, and a mass of air over
/// its sea that the scene states moving at `wind` m/s eastward over the ground
/// at ten metres.
fn an_earth_with_wind(wind: f64, r: f64) -> (World, NodeIdx, NodeIdx) {
    // Eight cells a face, the way `tests/ground.rs` states an Earth.
    let spec = SampleSpec::new(65, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let mut w = World::new(Tree::new(0xA1B, Matter::neutral(EARTH, RADIUS, 290.0, Composition::crustal()), Tier::Planetary, spec), 20.0);
    let earth = w.tree.root;
    let silica = w.substances.intern(substances::silica_arrangement()).unwrap();
    let h2o = w.substances.intern(substances::water_arrangement()).unwrap();
    let n2 = w.substances.intern(nitrogen()).unwrap();
    // The Earth's own ocean and the Earth's own atmosphere, by mass.
    let (water, air) = (1.4e21, 5.1e18);
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0 - (water + air) / EARTH);
    mix.add(h2o, Phase::Liquid, water / EARTH);
    mix.add(n2, Phase::Gas, air / EARTH);
    let composition = mix.composition(&w.substances).0;
    {
        let n = &mut w.tree.nodes[earth.get()];
        n.matter = Matter::neutral(EARTH, RADIUS, 290.0, composition);
        n.matter.spin = v3(0.0, 0.0, std::f64::consts::TAU / DAY * n.matter.moment_of_inertia());
        n.sync_spin_rate();
    }
    w.set_mixture(earth, mix);
    assert!(w.assess_ocean(earth), "an Earth with water has an ocean");
    // A planet already, with its surface written, so that nothing about it is
    // still to be decided once the air is over it.
    assert!(w.assess_surface(earth), "an Earth is a body its gravity has rounded");
    w.tree.refine(earth);


    // The air: a thousand kilometres of nitrogen, its centre ten metres over
    // the sea on the equator, at the density the Earth's own atmosphere has
    // there — which is what makes it weigh nothing in it.
    let mut air_mix = Mixture::new();
    air_mix.add(n2, Phase::Gas, 1.0);
    let sea = w.tree.nodes[earth.get()].ocean.as_ref().unwrap().radius;
    let at = v3(sea + 10.0, 0.0, 0.0);
    let ambient = w.ambient_density(earth, at.norm());
    let (_, rho0, height) = w.atmosphere_of(earth).expect("an Earth with nitrogen has an atmosphere");
    println!("  the atmosphere: {rho0:.3} kg/m^3 at the sea, scale height {:.2} km", height / 1e3);
    let mass = ambient * 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
    // Over the sea, held by the face of the ground it stands over — the way
    // anything near a planet's surface is — and placed there as something
    // loose among the face's cells rather than drawn from one of them, which
    // would take the cell out of the ground (`Tree::place`).
    let face = {
        let dir = v3(1.0, 0.0, 0.0);
        let f = match w.tree.nodes[earth.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
            Some(Recipe::Tiled(t)) => t.cell_of_direction(dir).expect("the Earth has a face under the equator"),
            _ => panic!("the Earth has no surface"),
        };
        let face = w.tree.promote(earth, f, w.tree.nodes[earth.get()].spec);
        w.tree.refine(face);
        face
    };
    let at = at - w.tree.offset_from(earth, face, Vec3::ZERO).value;
    let composition = air_mix.composition(&w.substances).0;
    let body = phys::state::Body {
        pos: at,
        mass,
        radius: r,
        temperature: 290.0,
        composition,
        ..Default::default()
    };
    let air = w.tree.place(face, body, SampleSpec::new(1, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain));
    // What it is made of: the face's own mixture came down with the slot, and
    // this is air.
    w.tree.nodes[air.get()].matter = Matter::neutral(mass, r, 290.0, composition);
    w.set_mixture(air, air_mix);
    // Relative to the face, which goes round with the Earth: the ground under
    // the air moves at the Earth's spin crossed with where the air is, less
    // what the face itself is doing there.
    let world_at = w.tree.offset_from(earth, air, Vec3::ZERO).value;
    let ground = w.tree.nodes[earth.get()].motion.spin_rate.cross(world_at) - w.tree.velocity_from(earth, face);
    w.tree.nodes[air.get()].motion.velocity = ground + v3(0.0, wind, 0.0);
    w.pace_fixed(60.0);
    (w, earth, air)
}

/// The air's speed over the water, m/s: its velocity less the ground's
/// turning where it is.
fn wind_of(w: &World, earth: NodeIdx, air: NodeIdx) -> f64 {
    let at = w.tree.offset_from(earth, air, Vec3::ZERO).value;
    let rel = w.tree.velocity_from(earth, air) - w.tree.nodes[earth.get()].motion.spin_rate.cross(at);
    (rel - at.unit().scale(rel.dot(at.unit()))).norm()
}

/// **A wind the scene states raises the sea it outruns, and pays for it.** A
/// thousand kilometres of air ten metres over the equator of a turning Earth,
/// blowing east at 10 m/s, for a day. The sea under it grows until its phase
/// speed reaches the wind's — the owner's limit — and stands there; the far
/// side of the planet stays calm; and the air slows by the stress it put into
/// the water, against the same air with no wind, which drifts by nothing.
///
/// Measured before, each with its fix in `World`: the ocean of a planet at the
/// root of its world stood still while the ground went round, and the air
/// swept 428 of 1536 cells and raised 0.95 m of sea on the far side; the wind
/// was read as one velocity, 14.8 m/s over the cells at the footprint's edge;
/// the face holding the air was solved every 240 to 1200 s against a 188 s
/// bob; the pull toward the axis was held for a step, and air at rest went
/// east at 0.0135 m/s an hour; the field was read where the face would be
/// rather than where its contents were; and the air sped up to 11.6 m/s.
#[test]
fn a_stated_wind_raises_the_sea_it_outruns_and_pays_for_it() {
    let run = |wind: f64| {
        let (mut w, earth, air) = an_earth_with_wind(wind, 1.0e6);
        let u0 = wind_of(&w, earth, air);
        // Held from its first frame on: `World::hold` is what says so.
        w.step_frame(1.0e6);
        if wind > 0.0 {
            let winds = w.winds_over(earth);
            assert!(!winds.is_empty(), "the air is over the sea and blows over nothing");
            let spread = winds.iter().map(|x| x.0.velocity.norm()).fold((f64::INFINITY, 0.0f64), |(a, b), u| (a.min(u), b.max(u)));
            println!("  the air covers {} cells at {:.3} to {:.3} m/s, from {u0:.3} m/s at its centre", winds.len(), spread.0, spread.1);
            // Not exactly even: a thing going round rigidly keeps its velocity
            // along a straight line, and the sea a footprint's width away
            // tilts under it by the arc between — cos 8 deg, 0.990.
            assert!(spread.1 - spread.0 < 0.02 * u0, "a held mass of air blows unevenly over its own footprint");
        }
        let mut lowest = f64::INFINITY;
        let mut highest = f64::NEG_INFINITY;
        // The world's books, read where the planet, the ground holding the
        // air and the world are at one instant.
        let face = w.tree.nodes[air.get()].parent;
        let mut first: Option<phys::state::Conserved> = None;
        let (mut turned, mut moved) = (0.0f64, 0.0f64);
        for _ in 1..(DAY / 60.0) as usize {
            w.step_frame(1.0e6);
            let (te, tf) = (w.tree.nodes[earth.get()].time, w.tree.nodes[face.get()].time);
            if (te - tf).abs() < 1e-6 && (te - w.time).abs() < 1e-6 {
                let c = w.conserved();
                let c0 = *first.get_or_insert(c);
                turned = turned.max((c.angular_momentum - c0.angular_momentum).norm() / c0.angular_momentum.norm());
                moved = moved.max((c.momentum - c0.momentum).norm());
            }
            let h = w.tree.offset_from(earth, air, Vec3::ZERO).value.norm() - w.tree.nodes[earth.get()].ocean.as_ref().unwrap().radius;
            lowest = lowest.min(h);
            highest = highest.max(h);
        }
        (w, earth, air, u0, lowest, highest, turned, moved)
    };
    let (w, earth, air, u0, lowest, highest, turned, moved) = run(10.0);
    let carried = w.tree.nodes[air.get()].matter.mass * w.tree.velocity_from(earth, air).norm();
    let u1 = wind_of(&w, earth, air);
    // The ground the air is held by, solved every frame of it: its pieces go
    // round with the planet. Drawn at the rate a uniform sphere of the face's
    // radius gives its angular momentum, they slipped at 5.6 to 8.2 km/s; with
    // their pull on each other through a Barnes-Hut tree, the ringing grew to
    // 22.8 m/s in a day; carrying their weight as a stress, they buckled; and
    // with a support turned to where the planet had been solved to rather
    // than where they were, they slipped at 48 m/s.
    let slip = {
        let face = w.tree.nodes[air.get()].parent;
        let fc = &w.tree.nodes[face.get()];
        let spin = w.tree.nodes[earth.get()].motion.spin_rate;
        let centre = w.tree.offset_at(earth, face, fc.time);
        let fv = w.tree.velocity_at(face, fc.time);
        let pieces = fc.ground.as_ref().map(|g| g.pieces).unwrap_or(0);
        assert!(pieces > 0, "the face holding the air was never solved as ground");
        fc.bodies[..pieces].iter().map(|b| (fv + b.vel - spin.cross(centre + b.pos)).norm()).fold(0.0f64, f64::max)
    };
    let o = w.tree.nodes[earth.get()].ocean.as_ref().unwrap();
    let up = w.tree.offset_from(earth, air, Vec3::ZERO).value.unit();
    let into = w.tree.facing(earth).conjugate();
    let under = o.wave_height(o.cell_of(into.rotate(up)));
    let far = o.wave_height(o.cell_of(into.rotate(Vec3::ZERO - up)));
    let lit = (0..o.cells.len()).filter(|&k| o.wave_height(k) > 0.01).count();
    let limit = |u: f64| phys::ocean::outrun_height(u, o.depth, o.g, o.tension, o.density);
    let (still, _, still_air, s0, _, _, _, _) = run(0.0);
    let drift = wind_of(&still, still.tree.root, still_air) - s0;
    println!(
        "  a day of it: {under:.3} m of sea under the air against an outrun height of {:.3} to {:.3} m; \
         {far:.1e} m on the far side; {lit} of {} cells over 1 cm; the air between {lowest:.2} and {highest:.2} m",
        limit(u1),
        limit(u0),
        o.cells.len()
    );
    println!("  the wind {u0:.4} -> {u1:.4} m/s; the same air with none drifted {drift:.1e} m/s");
    println!("  the ground under it: its pieces slip over the turning planet at {slip:.1e} m/s at worst");
    println!(
        "  the world's books over the day: angular momentum moved {turned:.1e} of itself, momentum {moved:.1e} kg m/s ({:.1e} of the air's)",
        moved / carried
    );
    assert!(under > 0.97 * limit(u1) && under < 1.01 * limit(u0), "the sea is not the one the wind outruns: {under} m");
    assert!(far < 1e-3, "the far side of the planet has a sea: {far} m");
    assert!(slip < 0.05, "the ground holding the air slips over its planet at {slip} m/s");
    // Measured before: 1.5e-7 of the angular momentum and 5e19 kg m/s, 2% of
    // the air's own, with the air's push on its planet kept on the face
    // rather than passed to the planet's own bodies.
    assert!(turned < 1e-7, "the world's angular momentum moved by {turned:.2e} of itself");
    assert!(moved < 1e-2 * carried, "the world's momentum moved by {moved:.2e} kg m/s");
    assert!(lowest > 0.0 && highest < 50.0, "the air left its level: {lowest} to {highest} m");
    assert!(drift.abs() < 1e-3, "held air with no wind moved at {drift} m/s");
    assert!(u0 - u1 > 10.0 * drift.abs().max(1e-3), "the air paid nothing for the sea it raised: {u0} -> {u1}");
}

