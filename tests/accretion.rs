//! A planet, generated from origin.
//!
//! The owner's done-when for Phase 4: *starting with a loose ball of suitable
//! matter in orbit around a star and ending with a planet.* Nothing in it
//! states that a planet is there. A ball of rock is put in orbit, cold and
//! sub-virial, and left alone; it falls in under its own gravity, and the first
//! frame in which its own weight has beaten its strength and its matter is
//! packed densely enough to carry itself is the frame it becomes one — at
//! which point the engine writes down the surface it has, and from then on the
//! ground is generated from that rule rather than from the collapse.
//!
//! # Watching it
//!
//! ```sh
//! PHYS_FILM=/tmp/accretion cargo test --test accretion -- --nocapture
//! ```
//!
//! writes a numbered PNG per frame into that directory, one camera for the
//! whole sequence. With the variable unset nothing is drawn and the test is
//! identical. See `docs/VIEWING.md`.

use phys::chem::{Mixture, Phase};
use phys::engine::World;
use phys::film::NodeFilm;
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::observe::Observer;
use phys::recipe::Recipe;
use phys::render::Paint;
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter, Spread};
use phys::tree::Tree;
use phys::units::*;

/// A frame budget at which the scheduler defers nothing, microseconds.
///
/// **The owner's decision for Phase 5.** Phase 4 ran this collapse at 20 ms and
/// left the starved budget open because the test passed either way, while the
/// planet it printed was a property of the machine — frame 76, 86 or 87, the
/// radius moving 5% between runs. Under Phase 5's rule that no node is left
/// unsolved, only scheduled, it stopped passing either way: at 20 ms one task
/// fits a frame, the sun won it every frame, the ball was not integrated once
/// in fifty, and when the span it owed passed what it could integrate it was
/// redrawn from its equilibrium and never collapsed at all. At this budget
/// nothing is deferred, which both tests assert.
const GENEROUS_US: f64 = 1.0e8;

/// A star, and a loose ball of rock in orbit at one astronomical unit.
///
/// The ball is stated the way a scenario states anything: a mass, a size, a
/// temperature and what it is made of. What it is *not* told is that it is a
/// planet, that it is round, or that it has a surface.
fn a_ball_in_orbit(seed: u64) -> (World, NodeIdx) {
    let mass = 6.0e24;
    let radius = 2.0e7;
    let spec = SampleSpec::new(64, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let star = Matter::neutral(2.0e30, 7.0e8, 5800.0, Composition::primordial());
    let mut w = World::new(Tree::new(seed, star, Tier::Stellar, spec), 20.0);
    w.pace_fixed(500.0);
    let sun = w.tree.root;
    w.tree.refine(sun);
    let ball = w.tree.promote(sun, 0, spec);
    {
        let n = &mut w.tree.nodes[ball.get()];
        n.matter = Matter::neutral(mass, radius, 300.0, Composition::primordial());
        n.matter.gravitational_binding = -0.6 * G * mass * mass / radius;
        // Sub-virial: it holds a twentieth of what would hold it up, so it is
        // bound and falling in. Equilibrium would be half.
        n.matter.internal_energy = 0.05 * n.matter.gravitational_binding.abs();
        n.motion.offset = v3(1.496e11, 0.0, 0.0);
        n.motion.velocity = v3(0.0, 29_780.0, 0.0);
        n.spec.count = 256;
    }
    let silica = w
        .substances
        .intern(phys::material::substances::silica_arrangement())
        .expect("silica analyses");
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0);
    w.set_mixture(ball, mix);
    (w, ball)
}

/// How big the node is *now*, measured the way `summarise` measures it.
fn size_of(w: &World, idx: NodeIdx) -> f64 {
    let n = &w.tree.nodes[idx.get()];
    if n.bodies.is_empty() {
        return n.matter.radius;
    }
    Spread::of(n.bodies.iter().map(|b| (b.pos, b.mass, b.radius))).rms * phys::state::RMS_TO_RADIUS
}

/// A loose ball of rock in orbit becomes a planet, and nothing said it would.
#[test]
fn a_loose_ball_of_matter_in_orbit_becomes_a_planet() {
    let (mut w, ball) = a_ball_in_orbit(0xACC4);
    let mass0 = w.tree.nodes[ball.get()].matter.mass;
    let size0 = size_of(&w, ball);

    // Somebody is watching it happen, which is what makes it happen: a cloud
    // nobody is looking at is never materialised, and an unmaterialised cloud
    // has no contents to fall in.
    w.add_observer(Observer {
        anchor: ball,
        offset: v3(0.0, 0.0, 8.0e7),
        look: v3(0.0, 0.0, -1.0),
        angular_resolution: 1e-5,
        ..Default::default()
    });
    let mut film = NodeFilm::open("accretion", ball, 3.0e7, Paint::Speed { range: None });

    assert!(
        w.tree.nodes[ball.get()].morphology.is_none(),
        "it starts as a ball of rock and nothing else"
    );

    let mut became: Option<(usize, f64)> = None;
    let mut smallest = f64::INFINITY;
    for frame in 0..200 {
        w.step_frame(GENEROUS_US);
        film.take(&w);
        let size = size_of(&w, ball);
        smallest = smallest.min(size);
        if became.is_none() && w.tree.nodes[ball.get()].morphology.is_some() {
            became = Some((frame, size));
        }
        // Mass is conserved whatever else happens.
        let now = w.tree.nodes[ball.get()].matter.mass;
        assert!(
            (now - mass0).abs() <= 1e-9 * mass0,
            "the ball weighed {now:.6e} kg at frame {frame}, not {mass0:.6e}"
        );
    }

    assert_eq!(w.stats.tasks_deferred, 0, "the scheduler deferred work, so this measures the machine");
    let (frame, size) = became.expect("the ball never became a planet");
    let n = &w.tree.nodes[ball.get()];
    let density = n.matter.mass / (4.0 / 3.0 * std::f64::consts::PI * n.matter.radius.powi(3));
    let g = G * n.matter.mass / (n.matter.radius * n.matter.radius);
    println!(
        "  a ball of {:.3e} kg, {:.4e} m across, fell in for {:.3e} s and became a planet \
         at frame {frame}: {:.4e} m, {density:.0} kg/m^3, {g:.3} m/s^2 at its surface",
        mass0,
        size0,
        w.time,
        size
    );
    if film.recording() {
        println!("  and {} frames of it are on disk", film.frames());
    }

    assert!(size < 0.6 * size0, "it should have fallen a long way in: {size:.4e} m");
    assert!(
        matches!(
            n.morphology.as_ref().and_then(|m| m.recipe.as_ref()),
            Some(Recipe::Tiled(t)) if t.is_ball()
        ),
        "and what it became is a body with a surface"
    );
    // **Compacted past the packing its grains jam at, and no further than rock
    // goes.** Phase 4 stopped a planet at random loose packing, 1743 kg/m^3,
    // because what squeezes one further is an equation of state for a solid;
    // Phase 5 derived one (`eos::Condensed::solid`) and `compacted_radius`
    // balances it against the planet's own weight. Measured: the ball is
    // assessed at frame 9 with a binding of 0.2715 GM^2/R — less bound than a
    // uniform sphere's 0.6 — and compacts to 2814 kg/m^3, above silica's own
    // 2644 and short of Earth's 5514, which takes an iron core this ball does
    // not have.
    let rest = w.eos_of(ball).condensed().map(|c| c.rest_density).expect("a ball of rock is condensed");
    assert!(
        density > rest && density < 5514.0,
        "a planet of silica compacts past its rest density {rest:.0} and short of Earth's 5514 kg/m^3, \
         and this is {density:.0}"
    );
    assert!(g > 1.0, "and you could stand on it: {g:.3} m/s^2");
}

/// And the planet it became has ground, which the collapse did not draw.
///
/// "Generated from origin, and then shortcutted through the generated program."
/// The collapse's detail is gone; what is left is a rule, and the ground comes
/// out of the rule.
#[test]
fn the_planet_it_became_has_ground_under_the_recipe() {
    let (mut w, ball) = a_ball_in_orbit(0xACC5);
    w.add_observer(Observer {
        anchor: ball,
        offset: v3(0.0, 0.0, 8.0e7),
        look: v3(0.0, 0.0, -1.0),
        angular_resolution: 1e-5,
        ..Default::default()
    });
    for _ in 0..200 {
        w.step_frame(GENEROUS_US);
        if w.tree.nodes[ball.get()].morphology.is_some() {
            break;
        }
    }
    assert!(w.tree.nodes[ball.get()].morphology.is_some(), "it never became a planet");
    assert_eq!(w.stats.tasks_deferred, 0, "the scheduler deferred work, so this measures the machine");

    // The rule is a handful of numbers, and the collapse's two hundred and
    // fifty-six bodies are not in it.
    let state = w.tree.nodes[ball.get()].morphology.as_ref().unwrap().state_bytes();
    println!("  the planet it became is {state} B of rule");
    assert!(state < 512);

    // Descend onto it: the cells are the children, all the way down.
    let mut here = ball;
    let mut levels = 0;
    for _ in 0..12 {
        let Some((cell, side)) = w.tree.nodes[here.get()].morphology.as_ref().and_then(|m| {
            let Some(Recipe::Tiled(t)) = m.recipe.as_ref() else { return None };
            Some((t.cell_of_direction(v3(0.4, 1.0, 0.2))?, t.side))
        }) else {
            break;
        };
        if side <= 2.0 {
            break;
        }
        w.tree.refine(here);
        if cell >= w.tree.nodes[here.get()].bodies.len() {
            break;
        }
        let spec = w.tree.nodes[here.get()].spec;
        let child = w.tree.promote(here, cell, spec);
        if child.is_none() {
            break;
        }
        here = child;
        levels += 1;
    }
    let side = match w.tree.nodes[here.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
        Some(Recipe::Tiled(t)) => t.side,
        _ => f64::INFINITY,
    };
    println!("  and {levels} levels below it the ground is {side:.4} m across");
    assert!(levels >= 6, "the descent stopped after {levels} levels");
    assert!(side < 100.0, "and arrived at {side:.4} m");

    // The surface it presents is made of what it is made of.
    let surface = w.surface_of_node(here).clone();
    assert!(!surface.is_empty(), "the ground presents no surface");
    let p = &surface.pieces()[0];
    println!(
        "  it presents {} pieces of {} at {:.0} kg/m^3",
        surface.len(),
        p.material.name,
        p.material.density
    );
    assert!(p.material.density > 0.0 && p.material.strength() > 0.0);
}
