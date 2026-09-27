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
use phys::math::Vec3;
use phys::observe::{Interaction, Observer};
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::Tree;
use phys::units::Tier;

const EARTH: f64 = 5.97e24;
const RADIUS: f64 = 6.371e6;

/// Where on the Earth the beach is: a little off the middle of a face of its
/// cube, so that the descent to it is over one cell at every level rather than
/// on the corner four of them share.
fn shore() -> Vec3 {
    v3(0.0123, 0.0217, 1.0).unit()
}

/// How the Earth divides its ground: 54 cells a side, so that four levels down
/// from a face a patch is 1.08 m across on the face's mean — 1.50 m here, near
/// the face's middle, where the cube map's cells are widest — and its cells
/// 2.8 cm: a 5 cm channel is two of them across.
const CELLS: usize = 54;

/// An Earth with its ground written and its sea on it, and **the beach the
/// Earth's own ground** — the owner's decision for Phase 5: a patch of it,
/// reached by `World::approach` under an observer standing a metre over the
/// sea at `shore()`, held and turned as ground is. The Earth's ground is a
/// smooth sphere at the sea's own level, so the scene states the one thing it
/// does not have — a slope, rising 0.15 towards +y, the side the sea's waves
/// run to — as two broad reshapings of the patch's surface, one up and one
/// down, either side of it. The sea over it is `sea` m high, running towards
/// +y. Returns the world, the Earth, the patch and the patch's half-side.
fn a_beach(seed: u64, sea: f64) -> (World, NodeIdx, NodeIdx, f64) {
    let spec = SampleSpec::new(CELLS * CELLS + 1, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
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
    // Somebody standing at the water's edge, who wants to see the ground at
    // their feet to a metre and a half.
    w.add_observer(Observer {
        anchor: earth,
        offset: shore().scale(level + 1.0),
        look: shore().scale(-1.0),
        angular_resolution: 1.5,
        ..Default::default()
    });
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
    // The slope, stated by the scene.
    let (y0, span) = (2.0, 2.0 * half * 0.92);
    let amplitude = 0.15 / (4.0 * y0 * half / (span * span) * (-(y0 * half / span).powi(2)).exp());
    if let Some(m) = w.tree.nodes[patch.get()].morphology.as_mut() {
        m.field.push(Deviation::reshaped(0.0, y0, span, amplitude));
        m.field.push(Deviation::reshaped(0.0, -y0, span, -amplitude));
    }
    w.tree.refine(patch);
    // And its mean surface at the sea's: the scene states whatever difference
    // there is between where a sphere's ground stands and where its sea does
    // as a lift of the whole patch. Measured, there is none to a tenth of a
    // millimetre — the sea's radius is the ground's.
    let lift = {
        let n = &w.tree.nodes[patch.get()];
        let centre = w.tree.offset_from(earth, patch, Vec3::ZERO).value;
        let up = centre.unit();
        let cells = CELLS * CELLS;
        let tops: Vec<f64> = n.bodies[..cells].iter().map(|b| (b.pos + b.orientation.rotate(v3(0.0, 0.0, b.half.z))).dot(up)).collect();
        let mean = tops.iter().sum::<f64>() / tops.len() as f64;
        let o = w.ocean_of(earth).unwrap();
        let sea = o.radius + o.cells[o.cell_of(shore())].eta - centre.norm();
        sea - mean
    };
    println!("  the scene lifts the beach {lift:.4} m to the sea's level");
    if let Some(m) = w.tree.nodes[patch.get()].morphology.as_mut() {
        m.field.push(Deviation::reshaped(0.0, 0.0, 1.0e4, lift));
    }
    w.tree.redraw_members(patch);
    {
        let o = w.ocean_of_mut(earth).unwrap();
        let c = o.cell_of(shore());
        o.cells[c].waves = phys::ocean::energy_of_height(sea, o.density, o.g);
        o.cells[c].heading = v3(0.0, 1.0, 0.0);
    }
    (w, earth, patch, half)
}

/// Where a point of the patch's surface is, and its two directions across,
/// in the axes its water is in: `(centre, x, y)`.
fn axes_of(w: &World, patch: NodeIdx) -> (Vec3, Vec3, Vec3) {
    let p = |x: f64, y: f64| w.field_point(patch, x, y).unwrap();
    let centre = p(0.0, 0.0);
    ((centre), (p(0.1, 0.0) - p(-0.1, 0.0)).unit(), (p(0.0, 0.1) - p(0.0, -0.1)).unit())
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
    w.advance_node(patch, 0.05);
    let (centre, ex, ey) = axes_of(&w, patch);
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
            let at = s.foot(k) - centre;
            if s.holds(k) && (at.dot(ex) - along).abs() < 0.025 && (at.dot(ey) - across).abs() < 0.5 * s.dx {
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
    // The period, from the times of the water's peaks there, each placed by
    // the parabola through its sample and the two either side of it. Not by
    // up-crossings of the mean: the channel drains between five and eight
    // seconds as the still water in it settles, the depth stays under its
    // own mean through two waves of it, and on the same run the up-crossings
    // read 0.967 s where the peaks read 0.9475.
    let late: Vec<&(f64, f64)> = trace.iter().filter(|p| p.0 > 5.0).collect();
    let peaks: Vec<f64> = late
        .windows(3)
        .filter(|w| w[1].1 > w[0].1 && w[1].1 >= w[2].1)
        .map(|w| {
            let bend = w[0].1 - 2.0 * w[1].1 + w[2].1;
            let off = if bend != 0.0 { 0.5 * (w[0].1 - w[2].1) / bend } else { 0.0 };
            w[1].0 + off * (w[2].0 - w[1].0)
        })
        .collect();
    let period = if peaks.len() > 1 { (peaks[peaks.len() - 1] - peaks[0]) / (peaks.len() - 1) as f64 } else { 0.0 };
    let train = {
        let o = w.ocean_of(_earth).unwrap();
        let c = o.cell_of(shore());
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
    let (along, from, to, across) = (0.33, -0.2, 0.5, 0.25);
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





