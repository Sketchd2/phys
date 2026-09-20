//! A boundary crossing is the event — `docs/PLAY.md` D16.
//!
//! The tree could always re-express a node in a new frame: `reparent` has done
//! the frame change, the slot handling and the rekeying correctly since Phase
//! 1. What it could not do was *notice*. `reparent` was only ever called by
//! hand, through `Interaction::Rehome`, so a thing that left the region its
//! parent owns simply stayed its child. Measured before any of this existed: a
//! rocket at escape velocity leaves a 1 km forest node in 0.100 s, and a
//! hundred seconds later it is 1.12x10^6 m away, still a child of a node that
//! claims a kilometre, with the spatial index clamping it into a corner cell
//! where it can neither collide with nor exchange heat with anything.
//!
//! The done-when of the phase is the rocket: it leaves the forest, re-homes to
//! the planet and then to the star, and keeps correct gravity and correct
//! neighbours throughout.

use phys::engine::{default_spec, World};
use phys::ids::NodeIdx;
use phys::math::{v3, Vec3};
use phys::state::{Composition, Matter};
use phys::tree::Tree;
use phys::units::*;

const EARTH_RADIUS: f64 = 6.371e6;
const EARTH_MASS: f64 = 5.972e24;
const AU: f64 = 1.496e11;

/// A star to hang a system off. Its own radius is the photosphere, which is
/// why the planet below sits outside it — see `the_universe_has_no_outside`.
fn a_star() -> Tree {
    let mass = 1.989e30;
    let radius = 6.957e8;
    let mut m = Matter::neutral(mass, radius, 5772.0, Composition::primordial());
    m.internal_energy = 0.3 * G * mass * mass / radius;
    m.gravitational_binding = -0.6 * G * mass * mass / radius;
    m.luminosity = 3.828e26;
    let mut spec = default_spec(Tier::Stellar);
    spec.count = 8;
    Tree::new(0x5A11, m, Tier::Stellar, spec)
}

/// Star, planet at 1 AU, a 1 km patch of forest on its surface, and something
/// sitting in the forest.
fn a_system(climb: f64) -> (World, NodeIdx, NodeIdx, NodeIdx, NodeIdx) {
    let mut w = World::new(a_star(), 20.0);
    let star = w.tree.root;
    w.tree.refine(star);
    let planet = w.tree.promote(star, 0, default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.spec.count = 8;
        n.motion.offset = v3(AU, 0.0, 0.0);
        n.motion.velocity = v3(0.0, 29_780.0, 0.0);
    }
    w.tree.refine(planet);
    let forest = w.tree.promote(planet, 0, default_spec(Tier::Continuum));
    {
        let n = &mut w.tree.nodes[forest.get()];
        // **Cold, and that is not a workaround.** A 291 K parcel of primordial
        // gas is a puff of hydrogen: it disperses at its own sound speed, about
        // 1.2 km/s, so within a couple of hundred frames it is twelve
        // kilometres across, no longer one neighbourhood, and `resolve_extent`
        // splits it — correctly. A patch of ground is not that, and a test
        // about boundaries should not be a test about a gas cloud coming apart.
        n.matter = Matter::neutral(1.0e9, 1000.0, 3.0, Composition::primordial());
        n.matter.internal_energy = n.matter.thermal_energy();
        n.spec.count = 16;
        n.motion.offset = v3(0.0, 0.0, EARTH_RADIUS);
        // At rest on the surface. `promote` hands a child the velocity of the
        // body it came from, and a sampled planet's parcels carry kilometres
        // per second of turbulence — which is right for a parcel of a hot
        // interior and wrong for a patch of ground somebody is standing on.
        n.motion.velocity = Vec3::ZERO;
    }
    w.tree.refine(forest);
    let rocket = w.tree.promote(forest, 0, default_spec(Tier::Continuum));
    {
        let n = &mut w.tree.nodes[rocket.get()];
        n.matter = Matter::neutral(5.0e5, 10.0, 300.0, Composition::primordial());
        n.motion.offset = Vec3::ZERO;
        n.motion.velocity = v3(0.0, 0.0, climb);
    }
    (w, star, planet, forest, rocket)
}

fn escape_velocity() -> f64 {
    (2.0 * G * EARTH_MASS / EARTH_RADIUS).sqrt()
}

/// The measurement D16 was written from, now with something watching.
///
/// Before the crossing pass existed this test's world ran two thousand frames
/// and recorded zero re-parents.
#[test]
fn a_rocket_leaves_the_forest_and_the_engine_notices() {
    let (mut w, _, planet, forest, rocket) = a_system(escape_velocity());
    let mut left_at = None;
    for _ in 0..40 {
        w.step_frame(2000.0);
        if w.tree.nodes[rocket.get()].parent != forest && left_at.is_none() {
            let n = &w.tree.nodes[rocket.get()];
            left_at = Some((w.time, n.motion.offset.norm() - EARTH_RADIUS));
        }
    }
    let (t, altitude) = left_at.expect("the rocket never left a node claiming a kilometre");
    println!("  left the forest at t = {t:.3} s, {altitude:.0} m up, now a child of the planet");
    assert!(t <= 0.5, "it left the forest's volume at 0.1 s and was noticed at {t:.3} s");
    assert_eq!(
        w.tree.nodes[rocket.get()].parent,
        planet,
        "nothing the planet holds contains it, so the planet is where it belongs"
    );
    assert_eq!(w.stats.crossings, 1, "one boundary, one event");
    assert_eq!(w.stats.crossings_sideways, 0, "it landed in no sibling");
}

/// The phase's done-when. Forest, planet, star — with gravity and neighbours
/// correct at every stage.
#[test]
fn the_rocket_reaches_the_star_and_keeps_its_gravity() {
    let (mut w, star, planet, forest, rocket) = a_system(escape_velocity());

    // What its own ancestry says its gravity is, against the closed form. The
    // engine walks the chain and applies the shell theorem per step; at the
    // surface that is one term and it has to be g.
    let g_now = |w: &World| w.tree.gravity_at(rocket).norm();
    let from_planet = |w: &World| {
        w.tree.separation(planet, Vec3::ZERO, rocket, Vec3::ZERO).value.norm()
    };
    assert!(
        (g_now(&w) - G * EARTH_MASS / (EARTH_RADIUS * EARTH_RADIUS)).abs() < 1e-2,
        "on the ground it should feel one g, not {:.4}",
        g_now(&w)
    );

    let mut hops: Vec<(f64, NodeIdx, f64, f64)> = Vec::new();
    let mut parent = w.tree.nodes[rocket.get()].parent;
    // A tenth of a second at one second per second, for the first boundary...
    for _ in 0..40 {
        w.step_frame(2000.0);
        let p = w.tree.nodes[rocket.get()].parent;
        if p != parent {
            hops.push((w.time, p, from_planet(&w), g_now(&w)));
            parent = p;
        }
    }
    // ...and then a coarse pace for the month-long climb out of the planet's
    // grip. `step_frame` takes a wall-clock budget; how much world time a frame
    // covers is the pace.
    w.pace_fixed(2000.0);
    for _ in 0..4000 {
        w.step_frame(2000.0);
        let p = w.tree.nodes[rocket.get()].parent;
        if p != parent {
            hops.push((w.time, p, from_planet(&w), g_now(&w)));
            parent = p;
        }
        if parent == star {
            break;
        }
    }

    for (t, p, d, g) in &hops {
        let what = if *p == planet { "planet" } else if *p == star { "star" } else { "?" };
        println!("  t {t:12.3} s  -> {what:6}   {d:.4e} m from the planet   g {g:.6e} m/s2");
    }
    assert_eq!(hops.len(), 2, "two boundaries: out of the forest, out of the planet");
    assert_eq!(hops[0].1, planet);
    assert_eq!(hops[1].1, star, "it never reached the star");

    // **Gravity is correct throughout.** At the first hop it is on the surface
    // and the planet is the only term that matters; at the second it is at the
    // planet's Hill radius, which is *defined* as where the star's pull takes
    // over, so the star's term is what it must read.
    let surface_g = G * EARTH_MASS / (hops[0].2 * hops[0].2);
    assert!(
        (hops[0].3 - surface_g).abs() / surface_g < 1e-3,
        "crossing into the planet it read {:.6e} where the planet's own field is {surface_g:.6e}",
        hops[0].3
    );
    let star_g = G * 1.989e30 / (AU * AU);
    assert!(
        (hops[1].3 - star_g).abs() / star_g < 0.1,
        "crossing into the star it read {:.6e} where the star's field at 1 AU is {star_g:.6e}",
        hops[1].3
    );
    let planet_g = G * EARTH_MASS / (hops[1].2 * hops[1].2);
    assert!(
        star_g > 10.0 * planet_g,
        "the hand-over should happen where the star has taken over: star {star_g:.3e} \
         against planet {planet_g:.3e}"
    );

    // **And correct neighbours.** The index a node is found in is its parent's,
    // and it has to be in the one it actually lives in.
    let held_by = |w: &World, parent: NodeIdx| {
        w.tree.nodes[parent.get()]
            .children
            .iter()
            .any(|c| *c == rocket)
    };
    assert!(held_by(&w, star), "the star does not hold what it is the parent of");
    assert!(!held_by(&w, planet), "the planet still holds a slot for something that left");
    assert!(!held_by(&w, forest), "the forest still holds a slot for something that left");
}

/// Nothing is created or destroyed by a crossing.
#[test]
fn a_crossing_conserves_everything() {
    let (mut w, _, _, _, _) = a_system(escape_velocity());
    let before = w.conserved();
    for _ in 0..40 {
        w.step_frame(2000.0);
    }
    assert!(w.stats.crossings >= 1, "nothing crossed, so nothing was tested");
    let after = w.conserved();
    let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(b.abs()).max(1e-300);
    println!(
        "  {} crossing(s): energy {:.3e} -> {:.3e}, baryon {:.3e}",
        w.stats.crossings,
        before.energy,
        after.energy,
        rel(before.baryon, after.baryon)
    );
    assert!(rel(before.energy, after.energy) < 1e-9, "energy changed across a crossing");
    assert!(rel(before.baryon, after.baryon) < 1e-9, "baryon number changed across a crossing");
}

/// Standing on a planet is not leaving it.
///
/// **A surface is exactly where a sphere's boundary is**, so a thing resting on
/// the ground has its bulk across the boundary of everything under it. That is
/// why leaving asks whether it is *wholly* outside and arriving asks whether it
/// is *wholly* inside: in the band between, nothing happens. Without the band a
/// world's entire population is ejected into interplanetary space on the first
/// frame.
#[test]
fn standing_on_the_ground_is_not_leaving_it() {
    let (mut w, _, _, forest, rocket) = a_system(0.0);
    // Resting on the forest's own boundary rather than in the middle of it:
    // the hardest case, and the one a naive rule gets wrong.
    w.tree.nodes[rocket.get()].motion.offset = v3(0.0, 0.0, 1000.0);
    for _ in 0..200 {
        w.step_frame(2000.0);
    }
    println!(
        "  200 frames sitting on the boundary: {} crossings, {} splits, {} merges, forest radius {:.1}, reach {:.1}; still a child of the forest: {}",
        w.stats.crossings,
        w.stats.splits,
        w.stats.merges,
        w.tree.nodes[forest.get()].matter.radius,
        w.tree.contents_reach(forest, NodeIdx::NONE),
        w.tree.nodes[rocket.get()].parent == forest
    );
    assert_eq!(
        w.tree.nodes[rocket.get()].parent, forest,
        "something standing still was re-homed"
    );
    assert_eq!(w.stats.crossings, 0, "nothing moved and something crossed");
}

/// A node in the tail of its parent's own draw has not left anything.
///
/// The regression that matters most, because the first implementation of this
/// pass got it wrong and the suite said so. `matter.radius` is the *equivalent
/// uniform sphere*, so a centrally-concentrated profile legitimately puts
/// bodies well beyond it — `docs/BACKLOG.md` measures a Plummer tail at three
/// to four radii. Measured on the ladder `drill_to` builds, with the radius as
/// the bound: seven promoted parcels sat at 1.0 to 3.6 of their parent's radius
/// and at **0.28 to 0.93 of what the parent's other contents reach**, and every
/// one of them was re-homed. The ladder came apart in four frames.
#[test]
fn a_node_in_the_tail_of_its_parents_draw_has_not_left() {
    let mut w = World::new(phys::engine::galaxy(0x1234, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    // Promote the parcels that sit furthest out, which is where the tail is.
    let mut by_distance: Vec<(usize, f64)> = w.tree.nodes[root.get()]
        .bodies
        .iter()
        .enumerate()
        .map(|(i, b)| (i, b.pos.norm()))
        .collect();
    by_distance.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    let radius = w.tree.nodes[root.get()].matter.radius;
    let mut promoted = Vec::new();
    for &(slot, d) in by_distance.iter().take(6) {
        let c = w.tree.promote(root, slot, default_spec(tier.finer()));
        assert!(!c.is_none());
        promoted.push((c, d / radius));
    }
    let reach = w.tree.contents_reach(root, NodeIdx::NONE);
    for (c, ratio) in &promoted {
        println!(
            "  a parcel at {ratio:.2} of its parent's radius, {:.2} of what the parent holds",
            w.tree.nodes[c.get()].motion.offset.norm() / reach
        );
    }
    assert!(
        promoted.iter().any(|(_, r)| *r > 1.0),
        "this test needs parcels outside the radius to mean anything: {promoted:?}"
    );
    for _ in 0..20 {
        w.step_frame(2000.0);
    }
    for (c, _) in &promoted {
        assert_eq!(
            w.tree.nodes[c.get()].parent, root,
            "a parcel in the tail of its parent's own draw was re-homed out of it"
        );
    }
}

/// A crossing into a place that is not a node yet **generates** it.
///
/// D16's inward row: the detail about to be interacted with is produced by the
/// crossing rather than by somebody having visited the place first. Without it
/// a creature walking onto an unvisited patch becomes a direct child of the
/// planet — "rekeyed twice", which D16 says must not happen — and its
/// neighbours for that frame are other regions rather than the ground under
/// its feet.
#[test]
fn a_crossing_generates_the_place_it_arrives_at() {
    let (mut w, _, planet, forest, rocket) = a_system(0.0);

    // A second patch of ground beside the forest, still only one of the
    // planet's bodies. Stated by hand because a sampled planet has no tiling
    // and `docs/PLAY.md` D16 says so outright: spheres do not tile a surface,
    // which is why the *trigger* lands here and the tiling lands in Ground.
    // Overlapping the forest, because that is the only way two spheres can
    // share a border: D16 says so outright, and it is why the *tiling* is
    // Ground's and only the trigger is this phase's.
    let next_door = v3(0.0, 2000.0, EARTH_RADIUS);
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.bodies[1].pos = next_door;
        n.bodies[1].vel = Vec3::ZERO;
        n.bodies[1].radius = 2000.0;
        n.bodies[1].mass = 1.0e9;
    }
    w.tree.pin(planet);
    println!(
        "  the forest's own contents reach {:.0} m; the patch next door spans {:.0} to {:.0} m",
        w.tree.contents_reach(forest, NodeIdx::NONE),
        next_door.y - 2000.0,
        next_door.y + 2000.0
    );

    // Walk out of the forest towards it, slowly enough that the only thing
    // that happens is the crossing.
    // Walk out of the forest towards it.
    w.tree.nodes[rocket.get()].motion.velocity = v3(0.0, 900.0, 0.0);
    let before = w.tree.live_count();
    let mut crossed_at = None;
    for f in 0..200 {
        w.step_frame(2000.0);
        if w.tree.nodes[rocket.get()].parent != forest {
            let here = w.tree.separation(planet, Vec3::ZERO, rocket, Vec3::ZERO).value;
            crossed_at = Some((f, here, (here - next_door).norm()));
            break;
        }
    }
    let (f, here, to_next_door) = crossed_at.expect("it never left the forest");
    println!(
        "  crossed at frame {f}, {:.0} m along, {to_next_door:.0} m from the centre of a patch \
         claiming 2000",
        here.y
    );
    let arrived = w.tree.nodes[rocket.get()].parent;
    println!(
        "  it left the forest for slot 1, which had to be promoted to receive it: \
         {} -> {} live nodes, {} sideways, {} generated",
        before,
        w.tree.live_count(),
        w.stats.crossings_sideways,
        w.stats.crossings_generated
    );
    assert_ne!(arrived, forest, "it never left");
    assert_ne!(arrived, planet, "it became a direct child of the planet, which D16 forbids");
    assert_eq!(w.tree.nodes[arrived.get()].parent, planet, "it should be inside the planet's other patch");
    assert_eq!(w.stats.crossings_sideways, 1, "this is the sideways case");
    assert_eq!(w.stats.crossings_generated, 1, "the place had to be generated to be arrived at");
    assert!(
        w.tree.nodes[arrived.get()].is_materialised(),
        "the detail about to be interacted with should exist by the time it arrives"
    );
}

/// The universe has no outside, and says so rather than doing something.
#[test]
fn the_universe_has_no_outside() {
    // A planet at 1 AU is genuinely outside the volume its star claims — a
    // photosphere is 6.96e8 m and an orbit is 1.5e11 — so the measurement is
    // right and there is simply nowhere to send it.
    let (mut w, star, planet, _, _) = a_system(0.0);
    for _ in 0..20 {
        w.step_frame(2000.0);
    }
    println!(
        "  {} refusals in 20 frames; the planet is still a child of the star",
        w.stats.crossings_refused
    );
    assert!(w.stats.crossings_refused > 0, "the measurement did not fire at all");
    assert_eq!(w.tree.nodes[planet.get()].parent, star, "the root lost a child it cannot lose");
}
