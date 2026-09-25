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
/// What this measures is the face *as a node* — its centre, against the
/// planet. The pieces drawn inside it are a separate question: they are drawn
/// turning at the rate a uniform sphere of the face's radius gives its
/// angular momentum, 0.63 of the ground's, and the support that holds them is
/// derived for the ground's; measured, the face's pieces slip over the
/// turning ground at 5.6 to 8.2 km/s by the sixth hour of the air test. That
/// is with the owner, in `docs/PLAY.md` §7.
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
