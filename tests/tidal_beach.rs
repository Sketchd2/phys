//! Phase 5's tidal done-when, from Phase 4's deviation clause:
//!
//! > a squiggle drawn in sand is gone by the next tide while a channel that
//! > redirects drainage is still there a year later.
//!
//! The beach is a patch of an Earth's own ground on its equator, under its
//! moon, with a light wind over the sea in front of it: the moon's tide alone
//! runs at a centimetre a second over the sand, below what moves a grain, and
//! it is the swell the wind raises that the tide brings up and down the beach.
//! The wind is the one thing the scene states (the owner's decision: a stated
//! wind, raised by the engine into the sea it outruns).
//!
//! **Watch it.** `PHYS_FILM=<dir> cargo test --test tidal_beach -- --nocapture`
//! writes the beach as a film: the patch's pieces and the water over it,
//! painted by how much of each is water, a frame an hour for the first tide
//! and a frame a month for the year.

use phys::chem::{Arrangement, Bond, Element, Mixture, Order, Phase};
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
const MOON: f64 = 7.342e22;
const DISTANCE: f64 = 3.844e8;
const DAY: f64 = 86_164.0;
/// The lunar semidiurnal tide: half of the time between the moon passing
/// overhead, 12.42 h.
const TIDE: f64 = 44_712.0;
/// A light breeze, m/s: the sea it outruns is 0.09 m at a metre a second and
/// goes as the square, so 0.21 m at this.
const WIND: f64 = 1.5;
/// A frame of the first tide, s.
const FRAME: f64 = 300.0;

/// Where the tidal beach is: on the equator, where the moon's tide is
/// semidiurnal, a little off the middle of the face of the cube over it.
fn equator() -> Vec3 {
    v3(1.0, 0.0123, 0.0217).unit()
}

/// Nitrogen, which is what the air is here.
fn nitrogen() -> Arrangement {
    Arrangement::molecule(vec![Element(7), Element(7)], vec![Bond::new(0, 1, Order::Triple)])
}

/// An Earth turning under its orbiting moon, with its ground written 54 cells
/// a side, its sea and its air, and a beach of its own ground on the equator:
/// reached as `beach.rs`'s is, and then watched from a camera linked to the
/// beach itself — the owner's decision, that a camera is linked to a node — so
/// that it goes round with the ground. The scene states the beach's slope, up
/// towards the Earth's north, and a wind of `wind` m/s eastward over the sea
/// in front of it, ten metres up. Returns the world, the Earth, the patch, its
/// half-side and the air.
fn a_tidal_beach(seed: u64, wind: f64) -> (World, NodeIdx, NodeIdx, f64, NodeIdx) {
    let spec = SampleSpec::new(2, Profile::Uniform, MassSpectrum::Equal, BodyKind::Planet);
    let pair = Matter::neutral(EARTH + MOON, 4.0e8, 290.0, Composition::crustal());
    let mut w = World::new(Tree::new(seed, pair, Tier::Planetary, spec), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    {
        let b = &mut w.tree.nodes[root.get()].bodies;
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
    let earth = w.tree.promote(root, 0, SampleSpec::new(CELLS * CELLS + 1, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain));
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
        let offset = n.motion.offset;
        n.matter = Matter::neutral(EARTH, RADIUS, 290.0, composition);
        n.matter.spin = v3(0.0, 0.0, std::f64::consts::TAU / DAY * n.matter.moment_of_inertia());
        n.motion.offset = offset;
        n.sync_spin_rate();
    }
    w.set_mixture(earth, mix);
    assert!(w.assess_surface(earth), "an Earth has ground");
    assert!(w.assess_ocean(earth), "and a sea on it");
    assert!(w.assess_atmosphere(earth), "and nitrogen over it, air");
    let level = w.ocean_of(earth).unwrap().radius;
    w.add_observer(Observer {
        anchor: earth,
        offset: equator().scale(level + 1.0),
        look: equator().scale(-1.0),
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
    // The camera, now linked to the beach it is looking at: a metre over it
    // along its own up, in its own axes, so that it goes round with it.
    let up = Vec3::ZERO - w.tree.nodes[patch.get()].gravity.unit();
    w.observers[0] = Observer { anchor: patch, offset: up, look: Vec3::ZERO - up, angular_resolution: 1.5, ..Default::default() };
    let (y0, span) = (2.0, 2.0 * half * 0.92);
    let amplitude = 0.15 / (4.0 * y0 * half / (span * span) * (-(y0 * half / span).powi(2)).exp());
    if let Some(m) = w.tree.nodes[patch.get()].morphology.as_mut() {
        m.field.push(Deviation::reshaped(0.0, y0, span, amplitude));
        m.field.push(Deviation::reshaped(0.0, -y0, span, -amplitude));
    }
    w.tree.refine(patch);

    // The air: three hundred kilometres of nitrogen ten metres over the sea in
    // front of the beach, at the density the Earth's own atmosphere has there
    // (which is what makes it weigh nothing in it), held by the face of the
    // ground it stands over as anything near a planet's surface is, and
    // placed among the face's cells rather than drawn from one.
    let mut air_mix = Mixture::new();
    air_mix.add(n2, Phase::Gas, 1.0);
    let at = equator().scale(level + 10.0);
    let ambient = w.ambient_density(earth, at.norm());
    let r: f64 = 3.0e5;
    let mass = ambient * 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
    let mut face = patch;
    while w.tree.nodes[face.get()].parent != earth {
        face = w.tree.nodes[face.get()].parent;
    }
    let at = w.tree.facing(face).conjugate().rotate(at - w.tree.offset_from(earth, face, Vec3::ZERO).value);
    let composition = air_mix.composition(&w.substances).0;
    let body = phys::state::Body { pos: at, mass, radius: r, temperature: 290.0, composition, ..Default::default() };
    let air = w.tree.place(face, body, SampleSpec::new(1, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain));
    w.tree.nodes[air.get()].matter = Matter::neutral(mass, r, 290.0, composition);
    w.set_mixture(air, air_mix);
    // Relative to the face, in its own turning axes: the ground under the air
    // goes round with it, so what the air is doing there is the wind and
    // nothing else. East is the way the ground turns.
    let east = v3(0.0, 0.0, 1.0).cross(equator()).unit();
    let into = w.tree.facing(face).conjugate();
    w.tree.nodes[air.get()].motion.velocity = into.rotate(east.scale(wind));
    w.pace_fixed(FRAME);
    (w, earth, patch, half, air)
}

/// The deviation of a given place, summed: a squiggle is one mark.
fn amplitude_of(w: &World, patch: NodeIdx, x: f64, y: f64) -> f64 {
    let field = &w.tree.nodes[patch.get()].morphology.as_ref().unwrap().field;
    field.iter().filter(|d| d.moved == 0.0 && (d.x - x).abs() < 1e-9 && (d.y - y).abs() < 1e-9).map(|d| d.amplitude).sum()
}

/// What is left of a channel cut down the beach: the depth its marks sum to,
/// and how much they took out.
fn channel_of(w: &World, patch: NodeIdx) -> (f64, f64) {
    let field = &w.tree.nodes[patch.get()].morphology.as_ref().unwrap().field;
    let cut = field.iter().filter(|d| d.moved > 0.0);
    (cut.clone().map(|d| d.amplitude).sum(), cut.map(|d| d.moved).sum())
}

/// How deep the channel is across its line at a point along it, m: the
/// surface the whole field makes beside the channel less the surface on its
/// line — what the channel's own marks take out.
fn depth_at(w: &World, patch: NodeIdx, half: f64, y: f64) -> f64 {
    let field = &w.tree.nodes[patch.get()].morphology.as_ref().unwrap().field;
    let surface = |x: f64| field.iter().map(|d| d.height_at(x, y, half)).sum::<f64>();
    surface(0.4 + 0.3) - surface(0.4)
}

/// The patch of ground the channel is in, wherever the tree has put it: a
/// node that is alive and holds a cut. A patch nobody watches is folded back
/// into the ground above it and drawn again when it is wanted, and the node
/// that was it is not the node that is.
fn the_beach(w: &World) -> Option<NodeIdx> {
    w.tree
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.alive)
        .find(|(_, n)| n.morphology.as_ref().is_some_and(|m| m.field.iter().any(|d| d.moved > 0.0)))
        .map(|(i, _)| NodeIdx(i as u32))
}

/// How high the beach is over its planet's mean sea, m.
fn height(w: &World, earth: NodeIdx, patch: NodeIdx) -> f64 {
    w.tree.offset_from(earth, patch, Vec3::ZERO).value.norm() - RADIUS
}

/// One frame of the film: the patch's own pieces and the water over it,
/// painted by how much of each is water — measured off the mixture, not
/// labelled (`film::fluidity`).
fn take(film: &mut Option<Film>, w: &World, patch: NodeIdx, half: f64) {
    let Some(film) = film.as_mut() else { return };
    let scene = film::of_node(w, patch);
    let values = film::fluidity(w, patch, &scene.bodies);
    let shot = Shot::framing(Vec3::ZERO, half * 1.4, 0.6, 0.35)
        .painted(Paint::Measured { values, range: Some((0.0, 1.0)) })
        .sized(320, 240);
    let (canvas, _) = film::shoot(w, patch, &shot);
    film.shoot(&canvas);
}

/// **The tidal half of the done-when, measured.** A wind raises a swell over
/// the sea in front of an Earth's beach on its equator. Once its sheet of
/// water is laid, a one-centimetre squiggle (which moved nothing) and a
/// five-centimetre channel (which took out sand, and so cannot be forgotten:
/// `docs/PLAY.md` §5.8) are cut into the sand where the tide reaches. Then a
/// lunar day of tide, and then a year.
///
/// **What is asserted is what holds**: the beach stays where it is over the
/// sea for a year; the channel's mass stays out; the squiggle is gone within
/// the year. **What is not met is reported, not asserted around**: the
/// squiggle is not gone by the next tide — the owner's decision is to report
/// it unmet, and `docs/PLAY.md` §7 records why and where the swash that would
/// erase it is scheduled. It slumps from 1 cm to 1.5 mm the moment it is cut
/// and then stays there under the swell.
#[test]
fn a_beach_under_a_tide_for_a_year() {
    let started = std::time::Instant::now();
    let (mut w, earth, patch, half, air) = a_tidal_beach(0xBEAD, WIND);
    let mut film = Film::open("tidal_beach");
    let density = match w.tree.nodes[patch.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
        Some(phys::recipe::Recipe::Tiled(t)) => t.density,
        _ => panic!("ground is tiled"),
    };
    // The sea over the beach first: a patch carved before its sheet is laid
    // is asked for one over ground already cut, and holds none.
    for _ in 0..3 {
        w.step_frame(1.0e6);
    }
    assert!(w.sheet_of(patch).is_some(), "the sea is not over the beach");
    let height0 = height(&w, earth, patch);
    take(&mut film, &w, patch, half);

    // A squiggle where the tide reaches, and a channel down the beach.
    let (sx, sy, squiggle) = (0.0, 0.3, 0.01);
    w.interact(Interaction::Mark { target: patch, deviation: Deviation::reshaped(sx, sy, 0.03, -squiggle) });
    let span = 0.05 / (2.0 * std::f64::consts::LN_2.sqrt());
    let mut y = -0.9;
    while y <= 0.9 {
        let mark = Deviation { span, ..Deviation::cut(0.4, y, span / half, -0.05 / 1.7725, 0.0) };
        let moved = density * mark.volume().abs();
        w.interact(Interaction::Mark { target: patch, deviation: Deviation { moved, ..mark } });
        y += span / half;
    }
    let (cut0, took) = channel_of(&w, patch);
    let cut_at = w.time;
    let line = [-0.6, -0.3, 0.0, 0.3, 0.6];
    let profile = |w: &World, patch: NodeIdx| line.iter().map(|&y| format!("{:.4}", depth_at(w, patch, half, y))).collect::<Vec<_>>().join(" ");
    println!("  a {squiggle} m squiggle and a channel 5 cm deep ({cut0:.3} m of marks, {took:.1} kg out of the beach) at t = {cut_at:.0} s");
    println!("  the channel's depth at y = {line:?}: {}", profile(&w, patch));
    let swell = |w: &World| {
        let o = w.ocean_of(earth).unwrap();
        (0..o.cells.len()).map(|k| o.wave_height(k)).fold(0.0f64, f64::max)
    };

    // The next tide: two of them, so the squiggle has had a whole one.
    let (mut wet, mut highest, mut fastest) = (0usize, 0.0f64, 0.0f64);
    let mut drift = 0.0f64;
    let frames = (2.0 * TIDE / FRAME) as usize;
    for f in 0..=frames {
        if let Some((_, flow)) = w.sea_over(patch, &Deviation::reshaped(sx, sy, 0.03, 0.0)) {
            wet += 1;
            fastest = fastest.max(flow);
        }
        highest = highest.max(swell(&w));
        drift = drift.max((height(&w, earth, patch) - height0).abs());
        if f % 12 == 0 {
            take(&mut film, &w, patch, half);
        }
        if f % 48 == 0 {
            println!(
                "  t {:6.0} s: the squiggle {:+.6} m, the channel {:.4} m of marks, the beach {:+.3} m over where it was",
                w.time,
                amplitude_of(&w, patch, sx, sy),
                channel_of(&w, patch).0,
                height(&w, earth, patch) - height0
            );
        }
        w.step_frame(1.0e6);
    }
    let (left, day) = (amplitude_of(&w, patch, sx, sy), channel_of(&w, patch));
    println!("  the channel's depth at y = {line:?}: {}", profile(&w, patch));
    println!(
        "  after {:.1} h of tide the squiggle is {left:+.6} m of {squiggle} ({:.0}% left), under the sea for {:.1} h of it at up to {fastest:.2} m/s; \
         the channel {:.4} m of marks and {:.1} kg out; the swell up to {highest:.2} m; the beach moved at most {drift:.3} m",
        (w.time - cut_at) / 3600.0,
        100.0 * left.abs() / squiggle,
        wet as f64 * FRAME / 3600.0,
        day.0,
        day.1
    );
    println!("  NOT MET: the squiggle is not gone by the next tide (the first clause of the tidal done-when)");
    assert!(drift < 1.0, "the beach wandered {drift} m over its sea in a day");

    // The wind drops, and the air with it: the scene states the air for the
    // tide and takes it away again (`Tree::take_away`). Left in a face of the
    // Earth whose solves are hours behind the world's clock, it is out of the
    // face within the month, and handing a node across a patch's edge in
    // anger is Phase 9's first scenario (`docs/BACKLOG.md`).
    let face = w.tree.nodes[air.get()].parent;
    assert!(w.tree.take_away(face, air), "the air could not be taken away");
    // The year, at a day a frame.
    w.pace_fixed(DAY);
    let mut months = 0;
    let mut wander = 0.0f64;
    let (mut swapped, mut lost) = (0, None);
    let mut patch = patch;
    for d in 0..365 {
        match the_beach(&w) {
            Some(now) => {
                if now != patch {
                    swapped += 1;
                    patch = now;
                }
                if d % 30 == 0 {
                    take(&mut film, &w, patch, half);
                    months += 1;
                }
                if d % 30 == 0 || w.tree.lca(earth, patch) != earth {
                    let mut chain = vec![patch.get()];
                    let mut n = patch;
                    while !w.tree.nodes[n.get()].parent.is_none() { n = w.tree.nodes[n.get()].parent; chain.push(n.get()); }
                    println!("  day {d}: the beach is node {}, chain up {:?}, earth is {}, alive {}; crossings {} sideways {} parted {} made {} rejoined {} folded {} splits {} merges {}", patch.get(), chain, earth.get(), w.tree.nodes[patch.get()].alive, w.stats.crossings, w.stats.crossings_sideways, w.stats.parted, w.stats.systems_made, w.stats.rejoined, w.stats.systems_folded, w.stats.splits, w.stats.merges);
                }
                wander = wander.max((height(&w, earth, patch) - height0).abs());
            }
            None => {
                lost = Some(d);
                break;
            }
        }
        w.step_frame(1.0e6);
    }
    assert!(lost.is_none(), "the channel was gone from the tree on day {lost:?}");
    take(&mut film, &w, patch, half);
    let (later, after) = (amplitude_of(&w, patch, sx, sy), channel_of(&w, patch));
    println!("  the node the beach was in changed {swapped} times in the year");
    println!("  the channel's depth at y = {line:?}: {}", profile(&w, patch));
    println!(
        "  a year on ({:.0} days): the squiggle {later:+.6} m; the channel {:.4} m of marks and {:.1} kg out, against {:.4} and {:.1} a day after; \
         the beach moved at most {wander:.3} m; {:.1} s computing",
        w.time / DAY,
        after.0,
        after.1,
        day.0,
        day.1,
        started.elapsed().as_secs_f64()
    );
    if let Some(f) = film.as_ref() {
        println!("  and {} frames of it are on disk ({months} of them one a month)", f.frames());
    }
    // **Not the bound it should be.** The Earth's ground is seven pieces, and
    // the moon's pull on them flexes the face under the beach by tens of
    // metres over a month, where a real solid tide is 0.3 m (a Love number
    // of 0.6 of the 0.54 m the sea's own equilibrium tide is here). What is
    // asserted is that it stays on the Earth; what it should be is reported
    // in the owner's words in `docs/PLAY.md` §7.
    println!("  NOT MET: the beach moved {wander:.1} m over its sea in the year, where the solid Earth's own tide is 0.3 m");
    assert!(wander < 1.0e3, "the beach wandered {wander} m over its sea in a year");
    assert!(later.abs() < 1e-9, "the squiggle is still there a year on: {later} m");
    assert!((after.1 - took).abs() < 1e-9 * took, "the channel's mass came back: {:.3} kg of {took:.3}", after.1);
    assert!(after.0 < 0.0 && depth_at(&w, patch, half, 0.3) > 0.0, "there is no channel a year on");
}
