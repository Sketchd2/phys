//! What may be asked of the world.
//!
//! `tests/view.rs` holds the outbound half of the seam: a client can look and
//! cannot step. This is the inbound half, and its gate is different. Looking is
//! safe because it changes nothing; asking is not, because the thing asking is
//! not trusted.
//!
//! The requirement is "measure only, never told — actors will lie". These tests
//! hold it three ways: an act cannot express an assertion, an act is clamped to
//! what the actor's body can do, and what the world knows about an actor is
//! what an instrument measured.

use phys::control::{decode, encode, Act, Actor, ActorId, Admin, Command, Outcome, Refusal, Roster};
use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::observe::{Instrument, Property, Quantity};
use phys::state::Composition;
use phys::units::*;

fn a_world() -> (World, NodeIdx) {
    let mut w = World::new(galaxy(0xAC70, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 512;
    w.time_rate = 0.05;
    let root = w.tree.root;
    let d = *w.drill_to(root, Tier::Continuum.max_radius(), &default_spec).last().unwrap();
    w.tree.refine(d);
    w.pace_to(d);
    (w, d)
}

/// An actor embodied in the node it is acting on, so reach is never the thing
/// under test unless a test says so.
fn embodied(w: &World, at: NodeIdx) -> Actor {
    let mut a = Actor::person(ActorId(1), at);
    a.reach = w.tree.nodes[at.get()].matter.radius * 10.0;
    a
}

// ---------------------------------------------------------------------------
// an act is not an assertion
// ---------------------------------------------------------------------------

/// The type-level guarantee, stated as a test so that adding an assertion
/// variant to `Act` has to walk past it.
///
/// There is no way to write `Act::SetTemperature`. Authoring exists — it is
/// `Admin::Author` — and an actor cannot construct a `Command::Admin` unless
/// the transport let it, which is the whole access-control model.
#[test]
fn an_actor_cannot_state_a_fact() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let a = embodied(&w, d);
    roster.admit(&mut w, a);

    let before = w.tree.nodes[d.get()].matter.temperature;

    // Everything an actor can send, sent at once. None of it sets anything.
    for act in [
        Act::Push { target: d, force: v3(1.0, 0.0, 0.0), seconds: 1.0 },
        Act::Heat { target: d, watts: 1.0, seconds: 1.0 },
        Act::Sense { target: d, instrument: Instrument::Imager, quantity: Quantity::Temperature },
    ] {
        roster.apply(&mut w, &Command::Act { actor: ActorId(1), act });
    }

    // The only way to set a temperature is the administrative path.
    let after_acting = w.tree.nodes[d.get()].matter.temperature;
    assert_eq!(
        after_acting, before,
        "no act may set a matter property; heat is delivered as energy through the mailbox"
    );

    roster.apply(
        &mut w,
        &Command::Admin(Admin::Author {
            target: d,
            property: Property::Temperature,
            value: before * 2.0,
        }),
    );
    let after_admin = w.tree.nodes[d.get()].matter.temperature;
    println!("  {before:.3e} K -> {after_acting:.3e} K by acting -> {after_admin:.3e} K by authoring");
    assert!(after_admin > after_acting, "authoring is the path that can set things");
    assert!(!w.audit.is_empty(), "and it is audited");
}

/// An actor cannot act as somebody else, and an unknown actor cannot act at all.
#[test]
fn an_unknown_actor_is_refused() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let a = embodied(&w, d);
    roster.admit(&mut w, a);

    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(999),
            act: Act::Push { target: d, force: v3(1.0, 0.0, 0.0), seconds: 1.0 },
        },
    );
    assert_eq!(out, Outcome::Refused(Refusal::NoSuchActor));
    assert_eq!(roster.refused, 1);
}

// ---------------------------------------------------------------------------
// an act is bounded by what the actor can do
// ---------------------------------------------------------------------------

#[test]
fn an_actor_cannot_push_harder_than_it_can_push() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let actor = embodied(&w, d);
    let (max_force, max_duration) = (actor.max_force, actor.max_duration);
    roster.admit(&mut w, actor);

    // Within its means.
    let ok = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Push { target: d, force: v3(max_force * 0.5, 0.0, 0.0), seconds: 1.0 },
        },
    );
    assert_eq!(ok, Outcome::Accepted);

    // A claim to push a planet.
    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Push {
                target: d,
                force: v3(1e18, 0.0, 0.0),
                seconds: max_duration * 1e6,
            },
        },
    );
    match out {
        Outcome::Clamped { fraction } => {
            println!("  asked for 1e18 N for {:.0} s; got {fraction:.3e} of it", max_duration * 1e6);
            // Force and duration are each clamped, so the impulse that reaches
            // the engine is bounded by the product of the two limits.
            let expected = (max_force / 1e18) * (1.0 / 1e6);
            assert!((fraction / expected - 1.0).abs() < 1e-9);
        }
        other => panic!("an impossible push should be clamped, got {other:?}"),
    }
    assert_eq!(roster.clamped, 1);
}

#[test]
fn power_and_flow_are_bounded_too() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let actor = embodied(&w, d);
    let max_power = actor.max_power;
    roster.admit(&mut w, actor);

    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Heat { target: d, watts: max_power * 1000.0, seconds: 1.0 },
        },
    );
    assert!(matches!(out, Outcome::Clamped { .. }), "got {out:?}");

    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Release {
                target: d,
                kg_per_second: 1e12,
                seconds: 1.0,
                composition: Composition::solar(),
                velocity: v3(0.0, 0.0, 0.0),
            },
        },
    );
    assert!(matches!(out, Outcome::Clamped { .. }), "got {out:?}");
}

/// Cooling is a negative rate, and must be clamped by magnitude like any other.
#[test]
fn cooling_is_bounded_by_the_same_limit_as_heating() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let a = embodied(&w, d);
    roster.admit(&mut w, a);

    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Heat { target: d, watts: -1e15, seconds: 1.0 },
        },
    );
    match out {
        Outcome::Clamped { fraction } => {
            println!("  asked to remove 1e15 W; got {fraction:.3e} of it");
            assert!(fraction > 0.0 && fraction < 1e-10);
        }
        other => panic!("expected a clamp, got {other:?}"),
    }
}

#[test]
fn an_actor_cannot_reach_what_is_far_away() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let mut actor = Actor::person(ActorId(1), d);
    actor.reach = 1e-3;
    roster.admit(&mut w, actor);

    let parent = w.tree.nodes[d.get()].parent;
    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Push { target: parent, force: v3(1.0, 0.0, 0.0), seconds: 1.0 },
        },
    );
    println!("  reaching for the parent node from 1 mm of reach: {out:?}");
    assert_eq!(out, Outcome::Refused(Refusal::OutOfReach));
}

#[test]
fn a_magnitude_that_is_not_a_number_is_refused() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let a = embodied(&w, d);
    roster.admit(&mut w, a);

    for act in [
        Act::Push { target: d, force: v3(f64::NAN, 0.0, 0.0), seconds: 1.0 },
        Act::Heat { target: d, watts: f64::INFINITY, seconds: 1.0 },
        Act::Heat { target: d, watts: 1.0, seconds: f64::NAN },
    ] {
        let out = roster.apply(&mut w, &Command::Act { actor: ActorId(1), act });
        assert_eq!(out, Outcome::Refused(Refusal::NotANumber));
    }
}

// ---------------------------------------------------------------------------
// what the world knows about an actor, it measured
// ---------------------------------------------------------------------------

#[test]
fn an_actor_cannot_sense_with_what_it_does_not_carry() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let mut actor = embodied(&w, d);
    actor.instruments = vec![Instrument::Imager];
    roster.admit(&mut w, actor);

    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Sense {
                target: d,
                instrument: Instrument::MassSpectrometer,
                quantity: Quantity::Temperature,
            },
        },
    );
    assert_eq!(out, Outcome::Refused(Refusal::NoSuchInstrument));

    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Sense {
                target: d,
                instrument: Instrument::Imager,
                quantity: Quantity::Temperature,
            },
        },
    );
    println!("  with the instrument it carries: {out:?}");
    assert!(out.accepted());
}

/// A reading keeps its own uncertainty. Flattening it to one number would put
/// the lie in the API.
#[test]
fn a_reading_arrives_whole() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let a = embodied(&w, d);
    roster.admit(&mut w, a);

    let out = roster.apply(
        &mut w,
        &Command::Act {
            actor: ActorId(1),
            act: Act::Sense {
                target: d,
                instrument: Instrument::Thermometer,
                quantity: Quantity::Temperature,
            },
        },
    );
    match out {
        Outcome::Sensed(r) => {
            println!("  {r:?}");
            match *r {
                phys::observe::Reading::Temperature { kelvin, uncertainty } => {
                    assert!(kelvin > 0.0);
                    assert!(uncertainty > 0.0, "a reading without an uncertainty is a lie");
                }
                other => panic!("a thermometer should report a temperature, got {other:?}"),
            }
        }
        other => panic!("expected a reading, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// administrative commands
// ---------------------------------------------------------------------------

#[test]
fn an_admin_command_needs_no_actor() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    assert!(roster.is_empty(), "nobody is connected");

    let out = roster.apply(&mut w, &Command::Admin(Admin::Dilate { target: d, rate: 100.0 }));
    assert_eq!(out, Outcome::Accepted);
    assert_eq!(w.tree.nodes[d.get()].bubble, 100.0);
    assert_eq!(w.bubbles().len(), 1);
}

#[test]
fn an_over_ambitious_bubble_reports_what_it_got() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let out = roster.apply(&mut w, &Command::Admin(Admin::Dilate { target: d, rate: 1e30 }));
    match out {
        Outcome::Clamped { fraction } => {
            println!("  asked 1e30x, got {}x — fraction {fraction:.3e}", w.tree.nodes[d.get()].bubble);
            assert_eq!(w.tree.nodes[d.get()].bubble, phys::dilation::MAX_BUBBLE);
        }
        other => panic!("expected a clamp, got {other:?}"),
    }

    let out = roster.apply(&mut w, &Command::Admin(Admin::Dilate { target: d, rate: -1.0 }));
    assert_eq!(out, Outcome::Refused(Refusal::ImpossibleRate));
}

#[test]
fn a_command_against_nothing_is_refused() {
    let (mut w, _d) = a_world();
    let mut roster = Roster::new();
    let out = roster.apply(&mut w, &Command::Admin(Admin::Pin { target: NodeIdx::NONE }));
    assert_eq!(out, Outcome::Refused(Refusal::NoSuchTarget));
}

// ---------------------------------------------------------------------------
// it crosses a process boundary
// ---------------------------------------------------------------------------

#[test]
fn every_command_survives_the_wire() {
    let cmds = [
        Command::Act {
            actor: ActorId(7),
            act: Act::Push { target: NodeIdx(3), force: v3(1.0, -2.0, 3.5), seconds: 0.25 },
        },
        Command::Act {
            actor: ActorId(u64::MAX),
            act: Act::Heat { target: NodeIdx(0), watts: -12.5, seconds: 4.0 },
        },
        Command::Act {
            actor: ActorId(1),
            act: Act::Release {
                target: NodeIdx(9),
                kg_per_second: 0.5,
                seconds: 2.0,
                composition: Composition::solar(),
                velocity: v3(0.0, 1.0, 0.0),
            },
        },
        Command::Act {
            actor: ActorId(2),
            act: Act::Sense {
                target: NodeIdx(4),
                instrument: Instrument::MassSpectrometer,
                quantity: Quantity::Spin,
            },
        },
        Command::Admin(Admin::Dilate { target: NodeIdx(5), rate: 250.0 }),
        Command::Admin(Admin::Author {
            target: NodeIdx(6),
            property: Property::TimeRate,
            value: 3.0,
        }),
        Command::Admin(Admin::Pin { target: NodeIdx(8) }),
    ];
    for c in &cmds {
        let bytes = encode(c);
        let back = decode(&bytes).unwrap_or_else(|e| panic!("{c:?} did not decode: {e}"));
        assert_eq!(back, *c, "{} bytes did not round-trip", bytes.len());
    }
    println!("  {} command shapes round-tripped", cmds.len());
}

#[test]
fn a_malformed_command_does_not_panic() {
    let good = encode(&Command::Admin(Admin::Dilate { target: NodeIdx(1), rate: 2.0 }));
    for cut in 0..good.len() {
        let _ = decode(&good[..cut]);
    }
    for i in 0..good.len() {
        for bit in 0..8 {
            let mut bad = good.clone();
            bad[i] ^= 1 << bit;
            let _ = decode(&bad);
        }
    }
    println!("  {} truncations and {} flips, no panics", good.len(), good.len() * 8);
}

/// A command decoded from bytes is applied through exactly the same path as one
/// built in process, which is what makes the transport a transport rather than
/// a second implementation.
#[test]
fn a_command_off_the_wire_is_applied_like_any_other() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let bytes = encode(&Command::Admin(Admin::Dilate { target: d, rate: 42.0 }));
    let cmd = decode(&bytes).expect("decode");
    assert_eq!(roster.apply(&mut w, &cmd), Outcome::Accepted);
    assert_eq!(w.tree.nodes[d.get()].bubble, 42.0);
    println!("  {} bytes of command set a {}x bubble", bytes.len(), 42.0);
}

// ---------------------------------------------------------------------------
// the loop closes
// ---------------------------------------------------------------------------

/// The whole seam, both directions, with only bytes crossing it.
///
/// An actor sends a command as bytes; the world applies it, steps, and hands
/// back a scene as bytes; the actor reads the consequence out of the scene. At
/// no point does either side hold the other's types. This is the shape an AI
/// driving an actor actually needs — act, then look, and learn from the
/// difference — and it is the reason the two halves were built to the same wire
/// format.
#[test]
fn an_actor_can_act_and_then_see_what_changed() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let a = embodied(&w, d);
    let max_force = a.max_force;
    roster.admit(&mut w, a);

    let mut client = phys::view::Client::new();
    let req = phys::view::ViewRequest::of(d);

    // Look first. Everything the actor knows, it got from here.
    let before = phys::view::decode(&phys::view::encode(&client.frame(&w, &req)))
        .expect("the first scene decoded");
    let mass_before = before.node().mass;
    let momentum_before = w.tree.nodes[d.get()].matter.momentum;

    // Act. Bytes in.
    let cmd = encode(&Command::Act {
        actor: ActorId(1),
        act: Act::Push { target: d, force: v3(0.0, 0.0, max_force), seconds: 1.0 },
    });
    let out = roster.apply(&mut w, &decode(&cmd).expect("the command decoded"));
    assert_eq!(out, Outcome::Accepted);

    // The push is posted to the mailbox and arrives at the speed of light, so
    // it is not felt until the world has run. That delay is the point: an
    // actor's act is an influence crossing space, not an assignment.
    let mut steps = 0;
    while w.tree.nodes[d.get()].matter.momentum == momentum_before && steps < 200 {
        w.step_frame(20_000.0);
        steps += 1;
    }
    let delivered = w.tree.nodes[d.get()].matter.momentum - momentum_before;
    println!(
        "  {} bytes of command; the impulse arrived after {steps} frame(s), \
         changing momentum by {:.3e} kg m/s",
        cmd.len(),
        delivered.norm()
    );
    assert!(
        delivered.norm() > 0.0,
        "the push never arrived — an act that cannot be felt is not an act"
    );
    assert!(
        delivered.norm() <= max_force * 1.0 * (1.0 + 1e-9),
        "more momentum arrived than the actor could possibly have delivered"
    );

    // Look again. Bytes out.
    let after = phys::view::decode(&phys::view::encode(&client.frame(&w, &req)))
        .expect("the second scene decoded");
    println!(
        "  scene: instant {:.3e} -> {:.3e} s, mass {:.3e} -> {:.3e} kg",
        before.instant,
        after.instant,
        mass_before,
        after.node().mass
    );
    assert!(after.instant > before.instant, "the world moved and the scene says so");
}

/// An actor learns its limits the same way it learns everything else: by trying
/// something and being told what happened.
#[test]
fn an_actor_discovers_its_own_limits_by_being_told() {
    let (mut w, d) = a_world();
    let mut roster = Roster::new();
    let a = embodied(&w, d);
    let truth = a.max_force;
    roster.admit(&mut w, a);

    // Binary search on the clamp fraction, which is all an actor is ever told.
    // Nothing here reads `max_force`; it is only used to check the answer.
    let (mut lo, mut hi) = (0.0f64, 1e9f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let out = roster.apply(
            &mut w,
            &Command::Act {
                actor: ActorId(1),
                act: Act::Push { target: d, force: v3(mid, 0.0, 0.0), seconds: 0.5 },
            },
        );
        match out {
            Outcome::Accepted => lo = mid,
            Outcome::Clamped { .. } => hi = mid,
            other => panic!("unexpected {other:?}"),
        }
    }
    println!("  actor inferred its limit as {lo:.6} N; it is really {truth} N");
    assert!(
        (lo - truth).abs() / truth < 1e-6,
        "an actor should be able to find its own limit from the answers alone"
    );
}
