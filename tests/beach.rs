//! Phase 5's done-when: `docs/PLAY.md` §4's beach test, and the tidal half of
//! Phase 4's deviation clause.
//!
//! > An actor is standing on a beach. They carve a 5 cm channel in the sand.
//! > Do the waves flow through it, in real time?

use phys::chem::{Mixture, Phase};
use phys::engine::World;
use phys::erode::Deviation;
use phys::ids::NodeIdx;
use phys::material::substances;
use phys::math::v3;
use phys::observe::Interaction;
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::Tree;
use phys::units::Tier;

const EARTH: f64 = 5.97e24;
const RADIUS: f64 = 6.371e6;

/// An Earth with its ground written and its sea on it, and on its +z pole a
/// patch of sand about a metre across placed with its mean surface at the
/// sea's level — **the shore, stated by the scene**, the owner's decision
/// while ground on a sphere has no relief: where the patch's own relief is
/// below the sea it is under water, and where it is above, it is beach. The
/// sea over it is `sea` m high, running towards +x. Returns the world, the
/// Earth, the patch and the patch's half-side.
fn a_beach(seed: u64, sea: f64) -> (World, NodeIdx, NodeIdx, f64) {
    let spec = SampleSpec::new(65, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let planet = Matter::neutral(EARTH, RADIUS, 290.0, Composition::crustal());
    let mut w = World::new(Tree::new(seed, planet, Tier::Planetary, spec), 20.0);
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
    let level = w.ocean_of(earth).unwrap().radius;
    // Sand: 64 columns a side, so a column is under two centimetres and a
    // 5 cm channel is three of them across.
    let sand = 800.0;
    let patch_spec = SampleSpec::new(64 * 64, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let body = phys::state::Body { pos: v3(0.0, 0.0, level), mass: sand, radius: 0.6, temperature: 290.0, ..Default::default() };
    let patch = w.tree.place(earth, body, patch_spec);
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0);
    let composition = mix.composition(&w.substances).0;
    {
        let n = &mut w.tree.nodes[patch.get()];
        n.matter = Matter::neutral(sand, 0.6, 290.0, composition);
        n.spec = patch_spec;
    }
    w.set_mixture(patch, mix);
    w.emplace(patch, phys::morph::Program::Terrain, sand, None);
    w.tree.refine(patch);
    let half = match w.tree.nodes[patch.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
        Some(phys::recipe::Recipe::Tiled(t)) => 0.5 * t.side,
        _ => panic!("ground is tiled"),
    };
    // Its mean surface at the sea's level.
    let mean_top = {
        let n = &w.tree.nodes[patch.get()];
        let mask = n.structural_mask().unwrap();
        let tops: Vec<f64> = n.bodies.iter().zip(&mask).filter(|(_, o)| **o).map(|(b, _)| b.pos.z + b.half.z).collect();
        tops.iter().sum::<f64>() / tops.len() as f64
    };
    w.tree.nodes[patch.get()].motion.offset = v3(0.0, 0.0, level - mean_top);
    {
        let o = w.ocean_of_mut(earth).unwrap();
        let c = o.cell_of(v3(0.0, 0.0, 1.0));
        o.cells[c].waves = phys::ocean::energy_of_height(sea, o.density, o.g);
        o.cells[c].heading = v3(1.0, 0.0, 0.0);
    }
    (w, earth, patch, half)
}

/// The water that went along the channel's line, across it at `across`
/// metres up the beach, over `span` seconds of world: m^3, both ways.
struct Run {
    passed: f64,
    deepest: f64,
    /// The period the water at the line rises and falls at, s.
    period: f64,
    /// The sea's wave's own period there, s.
    train: f64,
    world: f64,
    computing: f64,
}

/// Run a beach for `span` s of world, carved first if `carve` is given, and
/// measure the water crossing the line `y = across` within the channel's
/// width about `x = along`.
fn run(sea: f64, carve: Option<&dyn Fn(&mut World, NodeIdx, f64)>, along: f64, across: f64, span: f64) -> Run {
    let (mut w, _earth, patch, half) = a_beach(0xBEAC, sea);
    {
        let o = w.ocean_of_mut(_earth).unwrap();
        let c = o.cell_of(v3(0.0, 0.0, 1.0));
        o.cells[c].heading = v3(0.0, 1.0, 0.0);
    }
    w.advance_node(patch, 0.05);
    let sheet = w.sheet_of(patch).expect("the sea lies over the shore");
    if let Some(cut) = carve {
        cut(&mut w, patch, half);
    }
    let started = std::time::Instant::now();
    let (mut t, mut passed, mut deepest) = (0.0, 0.0, 0.0f64);
    let mut trace: Vec<(f64, f64)> = Vec::new();
    while t < span {
        let dt = w.advance_node(patch, 0.05).dt_used;
        t += dt;
        let s = w.tree.nodes[sheet.get()].sheet.as_ref().unwrap();
        let mut depth = 0.0;
        for k in 0..s.depth.len() {
            let at = s.foot(k);
            if s.holds(k) && (at.x - along).abs() < 0.025 && (at.y - across).abs() < 0.5 * s.dx {
                // After the first five seconds: whatever the still sea fills
                // a cut with has mostly filled it by then.
                if t > 5.0 {
                    passed += s.flow[k][0].abs() * s.dx * dt;
                }
                deepest = deepest.max(s.depth[k]);
                depth += s.depth[k];
            }
        }
        trace.push((t, depth));
    }
    let computing = started.elapsed().as_secs_f64();
    let late: Vec<&(f64, f64)> = trace.iter().filter(|p| p.0 > 5.0).collect();
    let mean = late.iter().map(|p| p.1).sum::<f64>() / late.len().max(1) as f64;
    let ups: Vec<f64> = late.windows(2).filter(|w| w[0].1 < mean && w[1].1 >= mean).map(|w| w[1].0).collect();
    let period = if ups.len() > 1 { (ups[ups.len() - 1] - ups[0]) / (ups.len() - 1) as f64 } else { 0.0 };
    let train = {
        let o = w.ocean_of(_earth).unwrap();
        let c = o.cell_of(v3(0.0, 0.0, 1.0));
        phys::ocean::Train::of(o.wave_height(c), v3(0.0, 1.0, 0.0), o.depth, o.g, o.tension, o.density).map(|t| t.period()).unwrap_or(0.0)
    };
    Run { passed, deepest, period, train, world: t, computing }
}

/// **An actor carves a 5 cm channel in the sand, and the waves flow through
/// it, in real time.** `docs/PLAY.md` §4's beach test, which is Phase 5's
/// done-when.
///
/// The same shore twice from one seed, with the sea's 0.2 m sea running up
/// it: once as it is, and once after a channel 5 cm deep and 5 cm across at
/// half its depth is cut from the water up across the beach — a line of marks
/// (`Interaction::Mark`), each taking its sand away, so the ground changes
/// where it was cut and the sheet over it lies on the ground as it is — and
/// the same channel once more under a still sea, which floods a cut below its
/// level and sloshes in it for a while without any waves. Across the
/// channel's line a quarter of a metre up the dry sand, from five seconds on,
/// the water that passes is counted and its rise and fall timed; and the
/// carved run's computing is timed against the world it covered.
#[test]
fn waves_flow_through_a_channel_carved_in_the_sand() {
    let (along, from, to, across) = (0.33, -0.2, 0.6, 0.25);
    let cut = |w: &mut World, patch: NodeIdx, half: f64| {
        let density = match w.tree.nodes[patch.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
            Some(phys::recipe::Recipe::Tiled(t)) => t.density,
            _ => panic!("ground is tiled"),
        };
        // Gaussian marks a half-width apart sum to a groove 1.7725 times one
        // of them deep; each is scaled so the groove is 5 cm, and its width
        // at half that depth, `2 sqrt(ln 2)` half-widths, is 5 cm.
        let span = 0.05 / (2.0 * std::f64::consts::LN_2.sqrt());
        let amplitude = -0.05 / 1.7725;
        let mut y = from;
        while y <= to {
            let mark = Deviation::cut(along / half, y / half, span, amplitude, 0.0);
            let moved = density * mark.volume().abs();
            w.interact(Interaction::Mark { target: patch, deviation: Deviation { moved, ..mark } });
            y += span;
        }
    };
    let plain = run(0.2, None, along, across, 15.0);
    let still = run(1e-9, Some(&cut), along, across, 15.0);
    let carved = run(0.2, Some(&cut), along, across, 15.0);
    println!(
        "  across the channel's line {across} m up the beach, from 5 s to {:.1} s: {:.5} m^3 through the channel under the sea's waves, {:.5} under a still sea, {:.5} over uncut sand",
        carved.world, carved.passed, still.passed, plain.passed
    );
    println!(
        "  the water there rose and fell every {:.3} s against the waves' {:.3} s, and stood up to {:.4} m deep; the uncut sand {:.4} m",
        carved.period, carved.train, carved.deepest, plain.deepest
    );
    println!("  {:.3} s of computing for {:.2} s of world", carved.computing, carved.world);
    assert!(plain.deepest == 0.0 && plain.passed == 0.0, "the uncut sand there was not dry");
    assert!(carved.passed > 3.0 * still.passed, "what went through the channel was not mostly the waves");
    assert!((carved.period / carved.train - 1.0).abs() < 0.02, "the channel's water does not move with the waves");
    assert!(carved.computing < carved.world, "not in real time: {:.3} s for {:.2} s", carved.computing, carved.world);
}
