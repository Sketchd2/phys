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

const EARTH: f64 = 5.97e24;
const RADIUS: f64 = 6.371e6;

/// An Earth with its ground written and its sea on it, not turning, and on
/// its +z pole a patch of silicate ground under `water` kg of water whose
/// surface stands at the sea's. The sea over it runs `sea` m high towards +x.
fn a_shore(seed: u64, water: f64, count: usize, sea: f64) -> (phys::engine::World, phys::ids::NodeIdx, phys::ids::NodeIdx) {
    use phys::chem::{Mixture, Phase};
    use phys::material::substances;
    use phys::sampler::{MassSpectrum, Profile, SampleSpec};
    use phys::state::{BodyKind, Composition, Matter};
    let spec = SampleSpec::new(65, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let planet = Matter::neutral(EARTH, RADIUS, 290.0, Composition::crustal());
    let mut w = phys::engine::World::new(phys::tree::Tree::new(seed, planet, phys::units::Tier::Planetary, spec), 20.0);
    let earth = w.tree.root;
    let silica = w.substances.intern(substances::silica_arrangement()).unwrap();
    let h2o = w.substances.intern(substances::water_arrangement()).unwrap();
    let ocean = 1.4e21;
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0 - ocean / EARTH);
    mix.add(h2o, Phase::Liquid, ocean / EARTH);
    let composition = mix.composition(&w.substances).0;
    w.tree.nodes[earth.get()].matter = Matter::neutral(EARTH, RADIUS, 290.0, composition);
    w.set_mixture(earth, mix);
    assert!(w.assess_surface(earth), "an Earth has ground");
    assert!(w.assess_ocean(earth), "and a sea on it");
    let sand = 1500.0;
    let patch_spec = SampleSpec::new(count, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let body = phys::state::Body { pos: v3(0.0, 0.0, RADIUS), mass: sand + water, radius: 1.0, temperature: 290.0, ..Default::default() };
    let patch = w.tree.place(earth, body, patch_spec);
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, sand / (sand + water));
    mix.add(h2o, Phase::Liquid, water / (sand + water));
    let composition = mix.composition(&w.substances).0;
    {
        let n = &mut w.tree.nodes[patch.get()];
        n.matter = Matter::neutral(sand + water, 1.0, 290.0, composition);
        n.spec = patch_spec;
    }
    w.set_mixture(patch, mix);
    w.emplace(patch, phys::morph::Program::Terrain, sand, None);
    w.tree.refine(patch);
    // Its water's surface at the sea's: the top layer's centres are half a
    // spacing under it.
    let top = {
        let parcels = water_of(&w, patch);
        let s = (parcels[0].mass / w.tree.nodes[patch.get()].rest_density).cbrt();
        parcels.iter().map(|b| b.pos.z).fold(f64::MIN, f64::max) + 0.5 * s
    };
    let level = w.ocean_of(earth).unwrap().radius;
    w.tree.nodes[patch.get()].motion.offset = v3(0.0, 0.0, level - top);
    {
        let o = w.ocean_of_mut(earth).unwrap();
        let c = o.cell_of(v3(0.0, 0.0, 1.0));
        o.cells[c].waves = phys::ocean::energy_of_height(sea, o.density, o.g);
        o.cells[c].heading = v3(1.0, 0.0, 0.0);
    }
    (w, earth, patch)
}

/// A patch's water: its loose bodies that hold anything.
fn water_of(w: &phys::engine::World, patch: phys::ids::NodeIdx) -> Vec<phys::state::Body> {
    let n = &w.tree.nodes[patch.get()];
    let mask = n.structural_mask().expect("a patch of ground is a structure");
    n.bodies
        .iter()
        .enumerate()
        .filter(|(i, b)| b.mass > 0.0 && !mask.get(*i).copied().unwrap_or(false))
        .map(|(_, b)| *b)
        .collect()
}

/// **What crosses a patch of shore's open edge comes out of the sea's own
/// books**, exactly: the sea is a node (`World::assess_ocean`) whose matter is
/// far too large to hold a parcel of water, so what it gives and takes is
/// written in its account (`ocean::Account`). Over a patch of ground whose
/// water stands at the sea's level, for ten solves: the water the patch gained
/// is the mass the sea's account lost; the momentum the water gained, less
/// what the field and the floor gave it (`SolveReport::outside`), is the
/// momentum the sea's account lost; and the world's momentum moves by exactly
/// the outside push. Left out of the sea's books, all three move by what
/// crossed.
#[test]
fn what_crosses_the_edge_comes_out_of_the_seas_books() {
    let (mut w, earth, patch) = a_shore(0x5EA, 1500.0, 3000, 0.3);
    let account = |w: &phys::engine::World| w.ocean_of(earth).unwrap().account;
    let water = |w: &phys::engine::World| water_of(w, patch).iter().map(|b| b.mass).sum::<f64>();
    let momentum = |w: &phys::engine::World| w.tree.nodes[patch.get()].bodies.iter().fold(Vec3::ZERO, |a, b| a + b.momentum());
    let (a0, m0, p0, c0) = (account(&w), water(&w), momentum(&w), w.conserved());
    let mut outside = Vec3::ZERO;
    for _ in 0..10 {
        outside += w.advance_node(patch, 0.02).outside;
    }
    let (a1, m1, p1, c1) = (account(&w), water(&w), momentum(&w), w.conserved());
    let (came, went) = (w.stats.parcels_from_sea, w.stats.parcels_to_sea);
    let crossed = p1 - p0 - outside;
    let world = c1.momentum - c0.momentum - outside;
    println!(
        "  {came} parcels in and {went} out: the water {:+.6} kg, the sea's account {:+.6} kg",
        m1 - m0,
        a1.mass - a0.mass
    );
    println!(
        "  momentum across the edge {:.4} kg m/s, the sea's account off it by {:.2e}; the world off the outside push by {:.2e}",
        crossed.norm(),
        (crossed + (a1.momentum - a0.momentum)).norm(),
        world.norm()
    );
    assert!(came > 0 && went > 0, "nothing crossed");
    assert!(((m1 - m0) + (a1.mass - a0.mass)).abs() < 1e-9 * m0, "mass was made or lost at the edge");
    assert!((crossed + (a1.momentum - a0.momentum)).norm() < 1e-9 * crossed.norm(), "the sea's books do not hold what crossed");
    assert!(world.norm() < 1e-6 * crossed.norm(), "the world's momentum moved by more than the outside push");
}
