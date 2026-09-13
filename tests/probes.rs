//! Phase 0: the measurements `docs/PLAY.md` §6 says must exist before the plan
//! is designed further.
//!
//! `examples/probe.rs` prints these as a report; this pins them as assertions,
//! so a change to the engine has to walk past the number the plan was written
//! against. Several of them assert a **defect** — `docs/PLAY.md` §5.9's event
//! cap, D12's shape being blind to its conditions — because a characterisation
//! test is how you keep a known defect from moving while the suite stays green.
//! Each of those says, in the test, what it should become once it is fixed.

use phys::engine::{default_spec, galaxy, World, MAX_SUBSTEPS};
use phys::math::v3;
use phys::morph::{Environment, Event, EventKind, Morphology, Program};
use phys::rng::{Purpose, Stream};
use phys::solvers::dynamics::Dynamics;
use phys::solvers::frame::{Dof, Framework, Member};
use phys::topology::Material;
use phys::units::*;

const PROGRAMS: [Program; 6] = [
    Program::Tree,
    Program::Coral,
    Program::Tower,
    Program::Wall,
    Program::Terrain,
    Program::Settlement,
];

fn matured(program: Program, seed: u64) -> Morphology {
    let mut m = if program.is_planned() {
        Morphology::planned(program, 5.0e4, seed, 0x1234)
    } else {
        Morphology::new(program, seed, 0x1234, 0)
    };
    let mut env = Environment::default();
    env.labour = 0.02;
    env.reservoir_mass = 1.0e9;
    for _ in 0..400 {
        m.advance(YEAR * 0.25, &env);
    }
    m
}

// ---------------------------------------------------------------------------
// §3.4 — one clock puts the trajectory/ensemble boundary inside Continuum
// ---------------------------------------------------------------------------

/// At one second per second a `Continuum` node of a few tens of metres already
/// needs thousands of substeps, so it is crossed by its ensemble rather than
/// followed. This is the measurement §3.4 is written against, and it is what
/// makes the resolution floor real rather than arithmetic.
#[test]
fn the_trajectory_path_runs_out_inside_continuum() {
    let mut w = World::new(galaxy(0xAC70, 1e9), 20.0);
    w.pace_fixed(1.0 / 20.0);
    let span = w.frame_dt();
    assert!((span - 0.05).abs() < 1e-12, "frame span moved: {span}");

    let root = w.tree.root;
    let mut followed_coarsest: Option<f64> = None;
    let mut ensemble_finest: Option<f64> = None;
    for idx in w.drill_to(root, 1e-3, &default_spec) {
        let (tier, r) = {
            let n = &w.tree.nodes[idx.get()];
            (n.tier, n.matter.radius)
        };
        w.tree.refine(idx);
        let dt = w.node_dt(idx);
        let need = span / dt;
        if tier != Tier::Continuum {
            continue;
        }
        if need <= MAX_SUBSTEPS as f64 {
            followed_coarsest = Some(followed_coarsest.map_or(r, |x: f64| x.min(r)));
        } else {
            ensemble_finest = Some(ensemble_finest.map_or(r, |x: f64| x.max(r)));
        }
    }
    let followed = followed_coarsest.expect("some Continuum node is still followed");
    let ensemble = ensemble_finest.expect("some Continuum node falls to its ensemble");
    assert!(
        ensemble < followed,
        "the boundary should sit inside Continuum: followed down to {followed:.3e} m, \
         ensemble from {ensemble:.3e} m"
    );
}

/// The resolution floor §3.4 quotes, as arithmetic over the same constants the
/// engine uses. If `MAX_SUBSTEPS` or the frame rate move, these move with them.
#[test]
fn the_resolution_floor_is_what_the_plan_says() {
    let span = 0.05_f64;
    let floor = |c: f64| 4.0 * span * c / MAX_SUBSTEPS as f64;
    assert!((floor(340.0) - 0.2656).abs() < 1e-3, "air: {}", floor(340.0));
    assert!((floor(1500.0) - 1.1719).abs() < 1e-3, "water: {}", floor(1500.0));
    assert!((floor(5000.0) - 3.9062).abs() < 1e-3, "rock: {}", floor(5000.0));
}

// ---------------------------------------------------------------------------
// §5.9 — a structure regrows what was removed from it
// ---------------------------------------------------------------------------

/// **Characterisation of a defect.** `MAX_EVENTS` is 64 and compaction drops the
/// oldest half, so the sixty-fifth severance brings back every part the first
/// thirty-two named. When `docs/BACKLOG.md`'s "A structure regrows what was
/// removed from it" is fixed, this test should fail and be inverted: nothing
/// severed may ever return.
#[test]
fn the_sixty_fifth_severance_resurrects_the_first_thirty_two() {
    let mut m = matured(Program::Tree, 0xBEEF);
    let budget = 512;
    let sites: Vec<u32> = m.render(budget).site.iter().copied().take(65).collect();
    assert!(sites.len() == 65, "need 65 distinct sites, got {}", sites.len());

    for s in sites.iter().take(64) {
        m.record(Event { at: m.age, kind: EventKind::Severed, site: *s, magnitude: 0.001 }, 290.0);
    }
    let at64 = m.render(budget);
    let back64 = sites[..64].iter().filter(|s| at64.site.contains(s)).count();
    assert_eq!(back64, 0, "at 64 severances the log still suppresses every one");

    m.record(Event { at: m.age, kind: EventKind::Severed, site: sites[64], magnitude: 0.001 }, 290.0);
    let at65 = m.render(budget);
    let back65 = sites.iter().filter(|s| at65.site.contains(s)).count();
    assert!(
        back65 >= 32,
        "the cap should have resurrected the oldest half; {back65} of 65 came back"
    );
}

/// **Characterisation of a defect.** Only `render_branching` consults the event
/// log, so the four planned programs and terrain never suppress a severed site
/// at all — not after sixty-four edits, after one. When that is fixed this test
/// should fail and be inverted.
#[test]
fn four_programs_never_honour_a_severance() {
    for program in PROGRAMS {
        let mut m = matured(program, 0xBEEF);
        let sk = m.render(512);
        if sk.site.len() < 4 {
            continue;
        }
        let site = sk.site[sk.site.len() / 2];
        m.record(Event { at: m.age, kind: EventKind::Severed, site, magnitude: 0.001 }, 290.0);
        let after = m.render(512);
        let honoured = !after.site.contains(&site);
        let branching = matches!(program, Program::Tree | Program::Coral);
        assert_eq!(
            honoured, branching,
            "{}: honoured a severance = {honoured}, expected {branching} \
             (only the branching renderer consults `events`)",
            program.name()
        );
    }
}

// ---------------------------------------------------------------------------
// §5.8 — embodied energy
// ---------------------------------------------------------------------------

/// §5.8's sand-for-glass check rests entirely on manufactured things carrying
/// the energy that manufacturing put there. Every program that builds anything
/// does book it; terrain books zero, which is right, because rock was not made.
#[test]
fn every_building_program_books_its_embodied_energy() {
    for program in PROGRAMS {
        let m = matured(program, 0xB00C);
        let stored = m.stored_energy();
        if program == Program::Terrain {
            assert_eq!(stored, 0.0, "terrain should carry no embodied energy");
            assert!(m.built > 0.0, "terrain should still have mass");
            continue;
        }
        assert!(m.built > 0.0, "{}: built nothing", program.name());
        assert!(
            stored > 0.0,
            "{}: built {:.3e} kg and booked no embodied energy — §5.8's check \
             would be blind to a forgery of it",
            program.name(),
            m.built
        );
    }
}

/// The precision worry §5.8 raises, retired. Differencing a whole structure to
/// find its *smallest* member leaves about twelve significant digits, well clear
/// of the six the check needs.
#[test]
fn the_embodied_energy_check_has_the_digits() {
    for program in PROGRAMS {
        let m = matured(program, 0xD161);
        let whole = m.stored_energy();
        if whole <= 0.0 {
            continue;
        }
        let sk = m.render(512);
        let total: f64 = sk.mass.iter().sum();
        let smallest = sk
            .mass
            .iter()
            .map(|x| whole * x / total)
            .fold(f64::INFINITY, f64::min);
        let digits = 15.95 - (whole / smallest).log10();
        assert!(
            digits > 8.0,
            "{}: only {digits:.1} digits left after differencing for one member",
            program.name()
        );
    }
}

// ---------------------------------------------------------------------------
// D11 — an undesigned genome
// ---------------------------------------------------------------------------

/// D11's honesty test. A genome nobody designed must still render something
/// finite; if the space needs hand-tuning to work at all then coefficients have
/// merely replaced enum variants.
#[test]
fn an_undesigned_genome_renders_something_finite() {
    for program in PROGRAMS {
        let mut st = Stream::at(0x6E0E, 0x99, 0, Purpose::Structure);
        for i in 0..100 {
            let mut m = matured(program, 0x6E0E);
            for g in m.genome.iter_mut() {
                *g = st.uniform() as f32;
            }
            let sk = m.render(512);
            assert!(
                sk.pos.iter().all(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite()),
                "{} genome {i}: non-finite position",
                program.name()
            );
            assert!(
                sk.mass.iter().all(|x| x.is_finite() && *x >= 0.0)
                    && sk.radius.iter().all(|x| x.is_finite() && *x >= 0.0),
                "{} genome {i}: non-finite mass or radius",
                program.name()
            );
        }
    }
}

// ---------------------------------------------------------------------------
// D12 — shape against the conditions it grew in
// ---------------------------------------------------------------------------

/// **Characterisation of what D12 proposes to change.** `render` is pure in
/// `(genome, age, built, progress, events)` and never sees `Environment`, so a
/// tree grown in deep shade and one grown in full sun are the *same shape* once
/// they are the same size. Conditions reach form only through how much mass they
/// grew. When D12 lands this must fail.
#[test]
fn shape_is_blind_to_the_conditions_it_grew_in() {
    let mut grown = Vec::new();
    for light in [40.0, 400.0, 900.0] {
        let mut m = Morphology::new(Program::Tree, 0xF0BE, 0x1234, 0);
        let mut env = Environment::default();
        env.light_flux = light;
        env.reservoir_mass = 1.0e6;
        for _ in 0..400 {
            m.advance(YEAR * 0.25, &env);
        }
        grown.push(m);
    }
    // The conditions did change how much grew — otherwise this proves nothing.
    assert!(
        grown[2].built > grown[0].built * 100.0,
        "light should drive mass: {:.3e} against {:.3e}",
        grown[2].built,
        grown[0].built
    );

    let (age, built) = (grown[1].age, grown[1].built);
    let mut shapes = Vec::new();
    for m in &grown {
        let mut c = m.clone();
        c.age = age;
        c.built = built;
        shapes.push(c.render(256));
    }
    for (i, sk) in shapes.iter().enumerate().skip(1) {
        assert_eq!(sk.site.len(), shapes[0].site.len(), "part count differed at {i}");
        for (a, b) in sk.pos.iter().zip(shapes[0].pos.iter()) {
            assert!(
                a.x == b.x && a.y == b.y && a.z == b.z,
                "shape {i} differed from the first — D12 may have landed, in which \
                 case invert this test"
            );
        }
    }
}

/// **D12's stated risk, and it fires.** Two light histories with the same
/// integral and opposite ordering do not grow the same tree: the difference in
/// final mass is tens of percent. A per-species surface over *integrated*
/// conditions therefore cannot reproduce a life, and D12's shortcut has to be a
/// rate law over the current state instead.
#[test]
fn growth_is_not_markovian_in_the_light_integral() {
    let steps = 200;
    let run = |f: &dyn Fn(usize, usize) -> f64| {
        let mut m = Morphology::new(Program::Tree, 0x3A12, 0x1234, 0);
        let mut env = Environment::default();
        env.reservoir_mass = 1.0e6;
        for i in 0..steps {
            env.light_flux = f(i, steps);
            m.advance(YEAR * 0.5, &env);
        }
        m.built
    };
    let flat = run(&|_, _| 400.0);
    let rising = run(&|i, n| 80.0 + 640.0 * (i as f64) / (n as f64 - 1.0));
    let falling = run(&|i, n| 720.0 - 640.0 * (i as f64) / (n as f64 - 1.0));

    // The null control, because a test that cannot fail proves nothing: two
    // *identical* histories must come out identical, or the spread below is
    // measuring noise rather than ordering.
    let flat_again = run(&|_, _| 400.0);
    assert_eq!(flat, flat_again, "growth is not deterministic; the spread means nothing");

    let spread = (rising - falling).abs() / flat;
    assert!(
        spread > 0.2,
        "ordering barely mattered ({spread:.3}); if this has become small, growth \
         has become Markovian in the integral and D12's outcome surface is viable \
         after all"
    );
}

// ---------------------------------------------------------------------------
// D5 — can a limb reach a stride by bending?
// ---------------------------------------------------------------------------

/// **D5, confirmed and more strongly than it claimed.** A limb cannot be swung
/// through a stride by bending it: the member ruptures at about one percent of
/// tip travel, with the linear formulation still well inside its regime. So the
/// large rotation of locomotion has nowhere to live except a joint between
/// substructures, and corotational elements are not the answer — they would
/// permit a bend that the material does not.
#[test]
fn a_limb_ruptures_long_before_it_leaves_the_linear_regime() {
    let build = || Framework {
        joints: vec![v3(0.0, 0.0, 0.0), v3(0.0, 0.0, -0.4), v3(0.0, 0.0, -0.8)],
        members: vec![
            Member { a: 0, b: 1, radius: 0.03, truss: false, integrity: 1.0 },
            Member { a: 1, b: 2, radius: 0.03, truss: false, integrity: 1.0 },
        ],
        fixed: vec![true, false, false],
        material: Material::GREEN_WOOD,
        lumped: Vec::new(),
        mass_scale: 0.0,
        stiff_scale: 1.0,
    };

    // The largest load the limb survives, and how far it got.
    let mut best_travel = 0.0f64;
    let mut worst_ratio = 0.0f64;
    let mut broke_at = None;
    for force in [1.0, 5.0, 20.0, 60.0, 150.0, 400.0, 700.0, 1000.0] {
        let mut d = Dynamics::new(build());
        let mut load = vec![Dof::default(); 3];
        load[2].t = v3(force, 0.0, 0.0);
        let mut broken = 0usize;
        let mut diverged = false;
        for _ in 0..400 {
            let r = d.step(&load, 1.0e-3);
            broken += r.broken.len();
            if !r.displacement_ratio.is_finite() || r.displacement_ratio > 10.0 {
                diverged = true;
                break;
            }
        }
        if broken > 0 || diverged {
            broke_at = Some(force);
            break;
        }
        best_travel = best_travel.max(d.displacement[2].t.norm());
        worst_ratio = worst_ratio.max(d.displacement_ratio());
    }

    let broke_at = broke_at.expect("the limb should fail somewhere in this sweep");
    assert!(
        best_travel / 0.8 < 0.05,
        "the limb reached {:.3} of its own length before failing at {broke_at} N — \
         if this has grown, bending is a route to a stride after all",
        best_travel / 0.8
    );
    assert!(
        worst_ratio < 0.1,
        "chord rotation reached {worst_ratio:.3} before rupture; D5 assumed the \
         member breaks while still inside the linear regime"
    );
}

// ---------------------------------------------------------------------------
// D1/D10 — what a checkpoint costs
// ---------------------------------------------------------------------------

/// Only pinned nodes write their bodies, so a checkpoint of an *interactive*
/// subtree — where everything an actor touched is pinned — is the case that
/// matters. The marginal cost is about 181 bytes a body.
#[test]
fn a_pinned_body_costs_about_181_bytes_to_checkpoint() {
    let encode_with = |pin: bool| {
        let mut w = World::new(galaxy(0xC4EC, 1e9), 20.0);
        w.tree.nodes[0].spec.count = 4096;
        let root = w.tree.root;
        for idx in w.drill_to(root, 1.0e4, &default_spec) {
            w.tree.refine(idx);
            if pin {
                w.tree.nodes[idx.get()].pinned = true;
            }
        }
        let bodies: usize = w.tree.nodes.iter().filter(|n| n.alive).map(|n| n.bodies.len()).sum();
        (phys::persist::encode(w.view()).len(), bodies)
    };
    let (regen, bodies) = encode_with(false);
    let (pinned, _) = encode_with(true);
    assert!(bodies > 1000, "need a populated world; got {bodies} bodies");

    let per_body = (pinned - regen) as f64 / bodies as f64;
    assert!(
        (per_body - 181.0).abs() < 20.0,
        "a pinned body costs {per_body:.1} bytes, not the ~181 the plan is written \
         against"
    );
    assert!(
        regen < pinned / 10,
        "regenerable detail should cost far less than pinned: {regen} against {pinned}"
    );
}
