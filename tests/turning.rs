//! Phase 5 — what a turning body holds turns with it. The owner's decision:
//! a child held by its parent, in the fluid it describes or on its ground, is
//! carried between solves in the parent's turning frame, and anything nothing
//! holds keeps the straight line a thing no force acts on follows.

use phys::chem::{Arrangement, Bond, Element, Mixture, Order, Phase};
use phys::engine::World;
use phys::material::substances;
use phys::math::{v3, Vec3};
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::{turning_carry, Tree};
use phys::units::Tier;

const EARTH: f64 = 5.972e24;
const RADIUS: f64 = 6.371e6;
const DAY: f64 = 86_164.0;

/// **A held thing at rest on a turning body goes round with it, exactly.** A
/// quarter of a sidereal day turns a point on the equator a quarter of the way
/// round, at the ground's own speed; the straight line it used to be carried in
/// leaves it 5.49x10^6 m off the surface.
#[test]
fn a_held_thing_at_rest_goes_round_with_the_ground() {
    let w = v3(0.0, 0.0, std::f64::consts::TAU / DAY);
    let r = v3(RADIUS, 0.0, 0.0);
    let v = w.cross(r);
    let (at, vel) = turning_carry(w, r, v, 0.25 * DAY);
    let straight = r + v.scale(0.25 * DAY);
    println!(
        "  a quarter day: at {:?}, {:.6e} m from the centre; the straight line is {:.4e} m off the surface",
        at,
        at.norm(),
        straight.norm() - RADIUS
    );
    assert!((at - v3(0.0, RADIUS, 0.0)).norm() < 1e-6 * RADIUS);
    assert!((vel - w.cross(at)).norm() < 1e-9 * v.norm());
}

/// **What moves over the ground keeps its height.** Ten metres a second along
/// the surface for three hours stays on the sphere it started on; the straight
/// line in the turning frame would have climbed 915 m.
#[test]
fn a_held_thing_moving_over_the_ground_keeps_its_height() {
    let w = v3(0.0, 0.0, std::f64::consts::TAU / DAY);
    let r = v3(RADIUS + 10.0, 0.0, 0.0);
    let v = w.cross(r) + v3(0.0, 0.0, 10.0);
    let t = 3.0 * 3600.0;
    let (at, vel) = turning_carry(w, r, v, t);
    let climb = (10.0 * t).powi(2) / (2.0 * RADIUS);
    println!(
        "  three hours at 10 m/s: {:.6} m from where it started in height; a straight line would be {climb:.0} m up",
        at.norm() - r.norm()
    );
    assert!((at.norm() - r.norm()).abs() < 1e-6, "it did not keep its height");
    let rel = vel - w.cross(at);
    assert!((rel.norm() - 10.0).abs() < 1e-9, "its speed over the ground changed");
}

/// **A spinning planet's surface moves with its ground.** Each face of a tiled
/// Earth moves at the planet's own angular velocity crossed with where it is,
/// and what that leaves of the planet's angular momentum is the faces' own
/// turning. Measured before: points at the faces' centres were asked to hold
/// all of it, and went at 674.8 m/s over ground moving at 353.6.
#[test]
fn a_spinning_planets_faces_move_with_its_ground() {
    let spec = SampleSpec::new(65, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let mut w = World::new(
        Tree::new(0xA1B, Matter::neutral(EARTH, RADIUS, 290.0, Composition::crustal()), Tier::Planetary, spec),
        20.0,
    );
    let e = w.tree.root;
    let silica = w.substances.intern(substances::silica_arrangement()).unwrap();
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0);
    w.set_mixture(e, mix);
    {
        let n = &mut w.tree.nodes[e.get()];
        n.matter.spin = v3(0.0, 0.0, std::f64::consts::TAU / DAY * n.matter.moment_of_inertia());
        n.sync_spin_rate();
    }
    assert!(w.assess_surface(e));
    w.tree.refine(e);
    let n = &w.tree.nodes[e.get()];
    let omega = n.motion.spin_rate;
    let mut worst: f64 = 0.0;
    let mut l = Vec3::ZERO;
    for b in &n.bodies {
        worst = worst.max((b.vel - omega.cross(b.pos)).norm());
        l += b.pos.cross(b.momentum()) + b.spin;
    }
    let lost = (l - n.matter.spin).norm() / n.matter.spin.norm();
    println!("  worst slip of a face over its ground {worst:.3e} m/s; angular momentum {lost:.1e} of itself");
    assert!(worst < 1e-6, "a face slides over its own planet at {worst} m/s");
    assert!(lost < 1e-12, "the faces do not carry the planet's angular momentum");
}

/// Nitrogen, which is what the air is here.
fn nitrogen() -> Arrangement {
    Arrangement::molecule(vec![Element(7), Element(7)], vec![Bond::new(0, 1, Order::Triple)])
}

/// **An atmosphere is derived from the gas a planet holds.** The Earth's own
/// mass of nitrogen over an Earth with a sea: the pressure at the sea is its
/// weight over the area, the scale height `k T / m g`, and the density at the
/// sea their ratio — against the real 1.225 kg/m^3 and 8.5 km.
#[test]
fn an_atmosphere_is_the_gas_a_planet_holds() {
    let spec = SampleSpec::new(65, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let mut w = World::new(
        Tree::new(0xA1C, Matter::neutral(EARTH, RADIUS, 290.0, Composition::crustal()), Tier::Planetary, spec),
        20.0,
    );
    let e = w.tree.root;
    let silica = w.substances.intern(substances::silica_arrangement()).unwrap();
    let h2o = w.substances.intern(substances::water_arrangement()).unwrap();
    let n2 = w.substances.intern(nitrogen()).unwrap();
    let (water, air) = (1.4e21, 5.1e18);
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0 - (water + air) / EARTH);
    mix.add(h2o, Phase::Liquid, water / EARTH);
    mix.add(n2, Phase::Gas, air / EARTH);
    let composition = mix.composition(&w.substances).0;
    w.tree.nodes[e.get()].matter = Matter::neutral(EARTH, RADIUS, 290.0, composition);
    w.set_mixture(e, mix);
    assert!(w.atmosphere_of(e).is_none(), "no surface assessed, so nothing to stand the air on");
    assert!(w.assess_ocean(e));
    let (base, rho, height) = w.atmosphere_of(e).expect("an Earth with nitrogen and a sea has an atmosphere");
    let up = w.ambient_density(e, base + height);
    println!(
        "  at the sea {rho:.3} kg/m^3, scale height {:.2} km; one scale height up {:.3} of that",
        height / 1e3,
        up / rho
    );
    assert!((rho - 1.141).abs() < 0.01, "sea-level density {rho}");
    assert!((height - 8760.0).abs() < 50.0, "scale height {height}");
    assert!((up / rho - (-1.0f64).exp()).abs() < 1e-12);
    assert_eq!(w.ambient_density(e, base - 10.0), w.tree.nodes[e.get()].ocean.as_ref().unwrap().density);
}

/// A turning, tiled Earth made of silica, with the face under the equator
/// promoted into a node of its own and drawn.
fn an_earth_with_a_face(seed: u64) -> (World, phys::ids::NodeIdx, phys::ids::NodeIdx) {
    let spec = SampleSpec::new(65, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let mut w = World::new(
        Tree::new(seed, Matter::neutral(EARTH, RADIUS, 290.0, Composition::crustal()), Tier::Planetary, spec),
        20.0,
    );
    let e = w.tree.root;
    let silica = w.substances.intern(substances::silica_arrangement()).unwrap();
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0);
    w.set_mixture(e, mix);
    {
        let n = &mut w.tree.nodes[e.get()];
        n.matter.spin = v3(0.0, 0.0, std::f64::consts::TAU / DAY * n.matter.moment_of_inertia());
        n.sync_spin_rate();
    }
    assert!(w.assess_surface(e));
    w.tree.refine(e);
    let dir = v3(1.0, 0.0, 0.0);
    let f = match w.tree.nodes[e.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
        Some(phys::recipe::Recipe::Tiled(t)) => t.cell_of_direction(dir).unwrap(),
        _ => panic!("the Earth has no surface"),
    };
    let spec = w.tree.nodes[e.get()].spec;
    let face = w.tree.promote(e, f, spec);
    w.tree.refine(face);
    w.pace_fixed(60.0);
    (w, e, face)
}

/// **Ground holds together by its gravity against its own stiffness** — the
/// owner's decision for Phase 5. A face of a turning Earth, promoted into a
/// node of its own, stays where the planet draws it through a day of the
/// planet being solved, and goes round with the ground. Measured before: the
/// face met the planet's interior as a separate object, left at 76 m/s, and
/// was 1.46x10^6 m further out six hours later.
///
/// This is the face *as a node* — its centre, against the planet. What is
/// drawn inside it is `a_face_of_a_turning_planet_is_drawn_turning_with_it`.
#[test]
fn a_turning_planets_ground_holds_together() {
    let (mut w, e, face) = an_earth_with_a_face(0xA1D);
    let r0 = w.tree.offset_from(e, face, Vec3::ZERO).value.norm();
    let mut worst_r: f64 = 0.0;
    let mut worst_slip: f64 = 0.0;
    for _ in 0..(DAY / 60.0) as usize {
        w.step_frame(1.0e6);
        let at = w.tree.offset_from(e, face, Vec3::ZERO).value;
        let v = w.tree.velocity_from(e, face);
        let spin = w.tree.nodes[e.get()].motion.spin_rate;
        worst_r = worst_r.max((at.norm() - r0).abs());
        worst_slip = worst_slip.max((v - spin.cross(at)).norm());
    }
    println!(
        "  a day: the face's centre moved {worst_r:.3e} m at worst from {r0:.4e} m, and slid over its ground at {worst_slip:.3e} m/s; the Earth solved {} times",
        w.tree.nodes[e.get()].steps_taken
    );
    assert!(w.tree.nodes[e.get()].steps_taken > 0, "the Earth was never solved, so nothing here was tested");
    assert!(worst_r < 1.0, "the ground came apart: {worst_r} m");
    assert!(worst_slip < 0.01, "the face slid over its own planet at {worst_slip} m/s");
}

/// **A piece of a turning planet turns with it** — the owner's decision for
/// Phase 5, that a node's turning is an offset from its parent's. A face
/// promoted from a turning Earth and redrawn every frame is drawn where the
/// ground has turned to, going round at the ground's rate.
///
/// Measured before, with a node's turning read as its angular momentum over a
/// uniform sphere of its radius: the face turned at 23x the planet's rate, its
/// pieces were drawn at 0.63x, and its layout was redrawn where it started —
/// turned 0.0000 rad in six hours against the ground's 1.57.
#[test]
fn a_face_of_a_turning_planet_is_drawn_turning_with_it() {
    let (mut w, e, face) = an_earth_with_a_face(0xA1D);
    let spin = w.tree.nodes[e.get()].motion.spin_rate;
    let own = (w.tree.angular_velocity(face) - spin).norm();
    let across = |w: &World| {
        let fc = &w.tree.nodes[face.get()];
        (fc.bodies[63].pos - fc.bodies[0].pos).unit()
    };
    let before = across(&w);
    let t0 = w.tree.nodes[face.get()].time;
    let mut worst: f64 = 0.0;
    for _ in 0..(6 * 60) {
        w.step_frame(1.0e6);
        let fc = &w.tree.nodes[face.get()];
        let centre = w.tree.offset_at(e, face, fc.time);
        let fv = w.tree.velocity_at(face, fc.time);
        for b in fc.bodies.iter() {
            worst = worst.max((fv + b.vel - spin.cross(centre + b.pos)).norm());
        }
    }
    // Where the layout has turned to, against where the ground turning by the
    // same angle takes the same two pieces.
    let t = w.tree.nodes[face.get()].time - t0;
    let expected = phys::math::Quat::from_rate(spin, t).rotate(before);
    let off = across(&w).dot(expected).clamp(-1.0, 1.0).acos();
    println!(
        "  six hours: the face turns {own:.1e} rad/s apart from its planet; its layout is {off:.2e} rad \
         from where the ground turned it ({:.4} rad); its pieces slip over the ground at {worst:.1e} m/s",
        spin.norm() * t
    );
    assert!(own < 1e-12, "a piece of the planet turns on its own at {own} rad/s");
    assert!(off < 1e-6, "the layout is {off} rad from where the ground turned to");
    assert!(worst < 1e-2, "the face's pieces slip over their own ground at {worst} m/s");
}

/// **The world's books count what a child carries by moving** — the owner's
/// call for Phase 5. A promoted child speaks for the body that stood in for
/// it, and the books added what it holds in its own frame and never its going
/// round. A turning Earth with one face promoted then read a world momentum of
/// 1.9e26 kg m/s — exactly that face's — which moved by 2.2e26 in six hours as
/// the Earth was solved; a planet at rest has none.
#[test]
fn the_worlds_books_count_a_childs_motion() {
    let (mut w, e, face) = an_earth_with_a_face(0xA1D);
    let going_round = w.tree.nodes[face.get()].matter.mass * w.tree.velocity_from(e, face).norm();
    let p0 = w.conserved().momentum;
    let mut worst = p0.norm();
    for _ in 0..(6 * 60) {
        w.step_frame(1.0e6);
        worst = worst.max(w.conserved().momentum.norm());
    }
    println!(
        "  six hours: the world's momentum at most {worst:.2e} kg m/s, {:.1e} of one face's going round ({going_round:.2e})",
        worst / going_round
    );
    assert!(worst < 1e-9 * going_round, "the world has a momentum of {worst:.3e} kg m/s");
}
