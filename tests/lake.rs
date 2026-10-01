//! Phase 5's lakes and runoff: water standing on ground that no sea reaches.
//!
//! > A basin is cut in a patch of dry ground and water arrives on it — as a
//! > film over the whole patch, which the slope runs downhill — and it stands
//! > as a lake: a flat surface, the mass that arrived, the books closed, and
//! > still there a year later.
//!
//! **Watch it.** `PHYS_FILM=<dir> cargo test --test lake -- --nocapture`.

use phys::chem::{Mixture, Phase};
use phys::engine::World;
use phys::erode::Deviation;
use phys::film;
use phys::ids::NodeIdx;
use phys::material::substances;
use phys::math::{v3, Vec3};
use phys::observe::{Interaction, Observer};
use phys::render::{Film, Paint, Shot};
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::Tree;
use phys::units::Tier;

const EARTH: f64 = 5.97e24;
const RADIUS: f64 = 6.371e6;
const CELLS: usize = 54;
const DAY: f64 = 86_164.0;

fn site() -> Vec3 {
    v3(0.0123, 0.0217, 1.0).unit()
}

/// A dry Earth — its ground and no sea — and a patch of it reached by an
/// observer standing a metre over it, with a bowl in the patch: a broad
/// reshaping of its surface a hand deep. Returns the world, the patch, its
/// half-side and the water's mixture.
fn a_basin(seed: u64) -> (World, NodeIdx, f64, Mixture) {
    let spec = SampleSpec::new(CELLS * CELLS + 1, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let planet = Matter::neutral(EARTH, RADIUS, 290.0, Composition::crustal());
    let mut w = World::new(Tree::new(seed, planet, Tier::Planetary, spec), 20.0);
    let earth = w.tree.root;
    let silica = w.substances.intern(substances::silica_arrangement()).unwrap();
    let h2o = w.substances.intern(substances::water_arrangement()).unwrap();
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0);
    let composition = mix.composition(&w.substances).0;
    w.tree.nodes[earth.get()].matter = Matter::neutral(EARTH, RADIUS, 290.0, composition);
    w.set_mixture(earth, mix);
    assert!(w.assess_surface(earth), "an Earth has ground");
    w.add_observer(Observer { anchor: earth, offset: site().scale(RADIUS + 1.0), look: site().scale(-1.0), angular_resolution: 1.5, ..Default::default() });
    w.approach();
    let (patch, half) = w
        .tree
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.alive)
        .filter_map(|(i, n)| match n.morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
            Some(phys::recipe::Recipe::Tiled(t)) if !t.is_ball() => Some((NodeIdx(i as u32), t.field_half())),
            _ => None,
        })
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("the observer reached the ground");
    // The camera, linked to the patch it is looking at: a metre over it along
    // its own up, so that it goes round with it and keeps it drawn.
    let up = Vec3::ZERO - w.tree.nodes[patch.get()].gravity.unit();
    w.observers[0] = Observer { anchor: patch, offset: up, look: Vec3::ZERO - up, angular_resolution: 1.5, ..Default::default() };
    // The bowl: one broad reshaping, down.
    if let Some(m) = w.tree.nodes[patch.get()].morphology.as_mut() {
        m.field.push(Deviation::reshaped(0.0, 0.0, 0.8 * half, -0.03));
    }
    w.tree.redraw_members(patch);
    let mut water = Mixture::new();
    water.add(h2o, Phase::Liquid, 1.0);
    (w, patch, half, water)
}

/// The lake over a patch: its node, if there is one, and the surface its wet
/// columns stand at, as `(lowest, highest, wet columns, mass)`.
fn the_lake(w: &World, patch: NodeIdx) -> Option<(f64, f64, usize, f64)> {
    let node = w.sheet_of(patch)?;
    let s = w.tree.nodes[node.get()].sheet.as_ref()?;
    let wet: Vec<f64> = (0..s.depth.len()).filter(|&k| s.holds(k) && s.depth[k] > 1.0e-3).map(|k| s.bed[k] + s.depth[k]).collect();
    let lo = wet.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = wet.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    Some((lo, hi, wet.len(), s.mass()))
}

#[test]
fn rain_on_a_basin_stands_as_a_lake() {
    let (mut w, patch, half, water) = a_basin(11);
    let mut film = Film::open("lake");
    w.pace_fixed(1.0);
    for _ in 0..3 {
        w.step_frame(1.0e6);
    }
    assert!(the_lake(&w, patch).is_none(), "water stands on a dry patch");
    let before = w.tree.nodes[patch.get()].matter.mass;
    let mass = 382.0;
    w.interact(Interaction::Pour { target: patch, mass, liquid: water });
    w.step_frame(1.0e6);
    let after = w.tree.nodes[patch.get()].matter.mass;
    println!("  poured {mass} kg: the world's mass moved {:.3} kg", after - before);
    assert!((after - before - mass).abs() < 1e-6 * mass, "the world's books did not take the water");
    let mut shown = 0;
    for f in 0..600 {
        w.step_frame(1.0e6);
        if f % 20 == 0 {
            if let Some(film) = film.as_mut() {
                let scene = film::of_node(&w, patch);
                let values = film::fluidity(&w, patch, &scene.bodies);
                let shot = Shot::framing(Vec3::ZERO, half * 1.4, 0.6, 0.35).painted(Paint::Measured { values, range: Some((0.0, 1.0)) }).sized(320, 240);
                let (canvas, _) = film::shoot(&w, patch, &shot);
                film.shoot(&canvas);
                shown += 1;
            }
        }
    }
    let (lo, hi, wet, held) = the_lake(&w, patch).expect("water stands");
    let node = w.sheet_of(patch).expect("a sheet");
    let (covered, all, volume, bottom) = {
        let s = w.tree.nodes[node.get()].sheet.as_ref().unwrap();
        let covered = (0..s.bed.len()).filter(|&k| s.holds(k)).count();
        let volume: f64 = (0..s.bed.len()).filter(|&k| s.holds(k)).map(|k| s.depth[k]).sum::<f64>() * s.dx * s.dx;
        let bottom = (0..s.bed.len()).filter(|&k| s.holds(k)).map(|k| s.bed[k]).fold(f64::INFINITY, f64::min);
        (covered, s.density, volume, bottom)
    };
    println!(
        "  after ten minutes: {wet} of {covered} columns wet, surface {lo:.5} to {hi:.5} m (flat to {:.1e} m), lowest bed {bottom:.4} m, {held:.2} kg held of {mass} poured, volume {volume:.5} m^3 against {:.5} ({shown} frames filmed)",
        hi - lo,
        held / all
    );
    // Runoff: the film that landed on the whole patch is in the bowl.
    assert!(wet < covered / 2, "the water did not run downhill: {wet} of {covered} columns are wet");
    // A lake: flat, and exactly the water that stood.
    assert!(hi - lo < 1.0e-3, "the surface is not flat: {lo} to {hi} m");
    assert!(hi > bottom, "the water is below the bed's lowest point");
    assert!((volume - held / all).abs() < 1e-9 * volume.max(1e-30), "the volume is not the mass");
    assert!(held > 0.0 && held < mass, "the lake holds {held} kg of {mass}: the ground's pores take the rest");
    let patch_mass = w.tree.nodes[patch.get()].matter.mass;
    // The year, a day a frame.
    let level = (lo + hi) / 2.0;
    w.pace_fixed(DAY);
    for d in 0..365 {
        w.step_frame(1.0e6);
        if d % 30 == 0 {
            if let Some(film) = film.as_mut() {
                let scene = film::of_node(&w, patch);
                let values = film::fluidity(&w, patch, &scene.bodies);
                let shot = Shot::framing(Vec3::ZERO, half * 1.4, 0.6, 0.35).painted(Paint::Measured { values, range: Some((0.0, 1.0)) }).sized(320, 240);
                let (canvas, _) = film::shoot(&w, patch, &shot);
                film.shoot(&canvas);
            }
        }
    }
    let (lo2, hi2, wet2, held2) = the_lake(&w, patch).expect("the lake is still there a year later");
    println!(
        "  a year on: surface {lo2:.5} to {hi2:.5} m, {wet2} columns, {held2:.3} kg; the patch's mass {patch_mass:.3} -> {:.3}",
        w.tree.nodes[patch.get()].matter.mass
    );
    assert!((held2 - held).abs() < 1e-6 * held, "the lake gained or lost water: {held} -> {held2} kg");
    assert!(((lo2 + hi2) / 2.0 - level).abs() < 1e-3, "the level moved: {level} -> {}", (lo2 + hi2) / 2.0);
    assert!(hi2 - lo2 < 1.0e-3, "the surface is not flat a year on");
}
