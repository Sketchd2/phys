//! The other direction across the seam: what may be asked of the world.
//!
//! `view.rs` is what a client is allowed to *see*. This is what it is allowed
//! to *do*, and the two are deliberately asymmetric. A [`Scene`](crate::view::Scene)
//! is handed out freely because looking cannot change anything. A command is
//! checked, budgeted, clamped and answered, because the thing sending it is not
//! trusted.
//!
//! # Actors lie
//!
//! The requirement this module exists for is that an actor — a player, or a
//! model driving one — may say anything at all, and none of it may become true
//! by being said. Three separate mechanisms, because one would not be enough.
//!
//! **An act is not an assertion.** [`Act`] has no variant that states a fact.
//! There is no `SetHealth`, no `Declare`, no `IAmAt`. The engine's own
//! [`Interaction::Author`](crate::observe::Interaction::Author) can set a bulk
//! property directly and is exactly the kind of thing an actor must never
//! reach, so it is not in this enum — it is in [`Admin`], which an actor cannot
//! construct a [`Command`] for without authority. This is a type-level
//! guarantee rather than a validation rule: there is no code path to forget.
//!
//! **An act is bounded by what the actor can do.** A push is not "apply this
//! impulse", it is "push with this force, for this long", and it is clamped
//! against the [`Actor`]'s own limits before it reaches the engine. An actor
//! claiming to push a planet gets [`Outcome::Clamped`] and the number it
//! actually got. Nothing anywhere trusts the magnitude in the request.
//!
//! **What the world knows about an actor, it measured.** Sensing goes through
//! `World::measure`, which commits to the ledger and disturbs what it looks
//! at. An actor cannot sense with an instrument it does not carry, and it
//! learns its own limits the same way it learns everything else — by trying
//! something and being told what happened.
//!
//! # Why this is not just `Interaction`
//!
//! [`Interaction`](crate::observe::Interaction) is the engine's internal
//! vocabulary and the engine trusts it completely: `Impulse { dp }` applies
//! whatever momentum it is handed. That is correct for the engine's own use and
//! catastrophic as a public API. This module is the layer in between — every
//! command that leaves here has already been checked against an actor's means,
//! and only then becomes an `Interaction`.
//!
//! # Administrative commands are a separate vocabulary
//!
//! [`Admin`] holds the things that are not physics: time bubbles, direct
//! authoring, pinning. They are unphysical by construction, they apply
//! immediately rather than at the speed of light, and every one of them lands
//! in the audit trail. Keeping them in a different type from [`Act`] means the
//! question "can this sender do this?" is answered once, at the door, by which
//! constructor it was able to call.

use crate::ids::NodeIdx;
use crate::math::Vec3;
use crate::observe::{Instrument, Interaction, Property, Quantity, Reading};
use crate::state::Composition;
use crate::wire::{Reader, Result, Writer};

/// Identifies one actor. Assigned by the server; an actor does not choose it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ActorId(pub u64);

/// An embodied participant, and the limits of what it can do.
///
/// Every number here is a *capability*, not a state: the most this actor can
/// exert, deliver or carry. They are the server's record of the actor, never
/// the actor's claim about itself, which is the whole point.
#[derive(Debug, Clone)]
pub struct Actor {
    pub id: ActorId,
    /// The node this actor is embodied in — where it acts from, and the frame
    /// its senses are in.
    pub body: NodeIdx,
    /// Most force it can exert, newtons. A human is a few hundred.
    pub max_force: f64,
    /// Most power it can deliver or remove, watts.
    pub max_power: f64,
    /// Most mass it can release per second, kg/s.
    pub max_flow: f64,
    /// Longest a single act may run, seconds. Bounds the impulse a single
    /// command can deliver, so an actor cannot push gently for a geological
    /// age in one message.
    pub max_duration: f64,
    /// Instruments it carries. It cannot sense with what it does not have.
    pub instruments: Vec<Instrument>,
    /// How far it can reach, metres. Acts on anything further are refused.
    pub reach: f64,
    /// Where its senses sit, once [`Roster::admit`] has registered them.
    ///
    /// An actor that senses *is* an observer — it occupies a place, has a
    /// field of view and an angular resolution, and its looking is what makes
    /// the engine materialise what it is looking at. Set by `admit`, because
    /// an actor does not get to describe its own eyes.
    pub observer: Option<usize>,
}

impl Actor {
    /// A person-sized actor: a few hundred newtons, a few hundred watts, arms.
    pub fn person(id: ActorId, body: NodeIdx) -> Actor {
        Actor {
            id,
            body,
            max_force: 500.0,
            max_power: 400.0,
            max_flow: 1.0,
            max_duration: 10.0,
            instruments: vec![Instrument::Imager, Instrument::Thermometer],
            reach: 2.0,
            observer: None,
        }
    }

    /// The eyes this body implies. Person-scale by default: a wide field, a
    /// minute of arc, a tenth-second integration.
    pub fn senses(&self) -> crate::observe::Observer {
        crate::observe::Observer {
            anchor: self.body,
            angular_resolution: 3e-4,
            integration_time: 0.1,
            ..Default::default()
        }
    }
}

/// What an embodied actor may do.
///
/// Physical acts only, and every magnitude is a *rate* with a duration rather
/// than a total, so it can be checked against a capability. "Apply 10^12 N·s of
/// impulse" cannot be judged; "push with 400 N for 2 s" can.
#[derive(Debug, Clone, PartialEq)]
pub enum Act {
    /// Push on something.
    Push { target: NodeIdx, force: Vec3, seconds: f64 },
    /// Warm or cool it. Negative watts remove energy.
    Heat { target: NodeIdx, watts: f64, seconds: f64 },
    /// Release matter into it — exhaust, water, a reagent.
    Release {
        target: NodeIdx,
        kg_per_second: f64,
        seconds: f64,
        composition: Composition,
        velocity: Vec3,
    },
    /// Look at something. The only way anything ever learns anything, and an
    /// interaction in its own right because measurement disturbs.
    Sense { target: NodeIdx, instrument: Instrument, quantity: Quantity },
}

impl Act {
    pub fn target(&self) -> NodeIdx {
        match *self {
            Act::Push { target, .. }
            | Act::Heat { target, .. }
            | Act::Release { target, .. }
            | Act::Sense { target, .. } => target,
        }
    }
}

/// What an administrator may do. None of it is physics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Admin {
    /// Run a node and its subtree at `rate` seconds per second of world time.
    Dilate { target: NodeIdx, rate: f64 },
    /// Set a bulk property directly. The one path that can break conservation.
    Author { target: NodeIdx, property: Property, value: f64 },
    /// Keep a node's detail even when nobody is looking.
    Pin { target: NodeIdx },
}

impl Admin {
    pub fn target(&self) -> NodeIdx {
        match *self {
            Admin::Dilate { target, .. }
            | Admin::Author { target, .. }
            | Admin::Pin { target } => target,
        }
    }
}

/// One thing asked of the world.
///
/// The two variants are the whole access-control model: a sender that can only
/// build [`Command::Act`] cannot express an administrative act, whatever it
/// sends. A transport decides which constructor a connection is allowed to
/// reach; nothing downstream has to re-check.
#[derive(Debug, Clone, PartialEq)]
pub enum Command {
    Act { actor: ActorId, act: Act },
    Admin(Admin),
}

/// Why a command was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Refusal {
    /// No actor by that id. An actor cannot act as somebody else.
    NoSuchActor,
    /// The target does not exist, or is not alive.
    NoSuchTarget,
    /// Further away than this actor can reach.
    OutOfReach,
    /// The actor does not carry that instrument.
    NoSuchInstrument,
    /// A magnitude that is not a finite number. A caller bug.
    NotANumber,
    /// The measurement returned nothing.
    NothingToSense,
    /// A rate outside what the arithmetic can carry.
    ImpossibleRate,
}

impl std::fmt::Display for Refusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Refusal::NoSuchActor => "no such actor",
            Refusal::NoSuchTarget => "no such target",
            Refusal::OutOfReach => "out of reach",
            Refusal::NoSuchInstrument => "the actor does not carry that instrument",
            Refusal::NotANumber => "not a finite number",
            Refusal::NothingToSense => "nothing to sense",
            Refusal::ImpossibleRate => "impossible rate",
        };
        f.write_str(s)
    }
}

/// What happened.
///
/// An actor learns what it can do by being told, which is the only honest way:
/// publishing a capability table invites an actor to plan against it, and the
/// table would then be part of the API rather than a property of the body.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// Applied as asked.
    Accepted,
    /// Applied, but reduced to what the actor could actually manage. Carries
    /// the fraction that got through, so a model can learn its own limits by
    /// trying rather than by being handed a specification.
    Clamped { fraction: f64 },
    /// Not applied.
    Refused(Refusal),
    /// A measurement came back. The one command that returns information.
    ///
    /// Carries the whole [`Reading`], not a flattened number: every variant of
    /// it reports its own uncertainty, and an instrument that returns a value
    /// without one is lying. Reducing it here would put the lie in the API.
    Sensed(Box<Reading>),
}

impl Outcome {
    pub fn accepted(&self) -> bool {
        !matches!(self, Outcome::Refused(_))
    }
}

/// Clamp `asked` into `limit`, reporting the fraction that survived.
fn clamp_to(asked: f64, limit: f64) -> (f64, f64) {
    let limit = limit.max(0.0);
    if asked.abs() <= limit || asked == 0.0 {
        (asked, 1.0)
    } else {
        (asked.signum() * limit, limit / asked.abs())
    }
}

// ---------------------------------------------------------------------------
// applying
// ---------------------------------------------------------------------------

/// The actors a world knows about.
///
/// Kept beside the world rather than inside it: an actor is a property of the
/// *session*, not of the physics, and a world reloaded from disk has no opinion
/// about who is currently connected.
#[derive(Debug, Default)]
pub struct Roster {
    actors: Vec<Actor>,
    pub accepted: u64,
    pub clamped: u64,
    pub refused: u64,
}

impl Roster {
    pub fn new() -> Roster {
        Roster::default()
    }

    /// Admit an actor, registering the senses its body implies.
    ///
    /// Takes the world because admitting an actor changes it: something is now
    /// looking, and what is looked at has to be resolved. An actor cannot
    /// describe its own eyes — `Actor::senses` does, from its body.
    pub fn admit(&mut self, world: &mut crate::engine::World, mut actor: Actor) -> ActorId {
        let id = actor.id;
        actor.observer = Some(world.add_observer(actor.senses()));
        self.actors.retain(|a| a.id != id);
        self.actors.push(actor);
        id
    }

    pub fn get(&self, id: ActorId) -> Option<&Actor> {
        self.actors.iter().find(|a| a.id == id)
    }

    pub fn len(&self) -> usize {
        self.actors.len()
    }

    pub fn is_empty(&self) -> bool {
        self.actors.is_empty()
    }

    /// Apply one command to a world.
    ///
    /// Administrative commands are applied as given — whoever could construct
    /// an [`Admin`] has already been trusted by the transport. Actor commands
    /// are checked against the actor's body and its means first, and what
    /// reaches the engine is never the magnitude that was sent.
    pub fn apply(&mut self, world: &mut crate::engine::World, cmd: &Command) -> Outcome {
        let out = match cmd {
            Command::Admin(a) => self.apply_admin(world, *a),
            Command::Act { actor, act } => self.apply_act(world, *actor, act),
        };
        match out {
            Outcome::Refused(_) => self.refused += 1,
            Outcome::Clamped { .. } => self.clamped += 1,
            _ => self.accepted += 1,
        }
        out
    }

    fn alive(world: &crate::engine::World, idx: NodeIdx) -> bool {
        !idx.is_none() && idx.get() < world.tree.nodes.len() && world.tree.nodes[idx.get()].alive
    }

    fn apply_admin(&mut self, world: &mut crate::engine::World, a: Admin) -> Outcome {
        if !Self::alive(world, a.target()) {
            return Outcome::Refused(Refusal::NoSuchTarget);
        }
        match a {
            Admin::Dilate { target, rate } => match world.dilate(target, rate) {
                Some(applied) => {
                    if applied == rate {
                        Outcome::Accepted
                    } else {
                        Outcome::Clamped { fraction: applied / rate }
                    }
                }
                None => Outcome::Refused(Refusal::ImpossibleRate),
            },
            Admin::Author { target, property, value } => {
                if !value.is_finite() {
                    return Outcome::Refused(Refusal::NotANumber);
                }
                world.interact(Interaction::Author { target, property, value });
                Outcome::Accepted
            }
            Admin::Pin { target } => {
                world.interact(Interaction::Pin { target });
                Outcome::Accepted
            }
        }
    }

    fn apply_act(
        &mut self,
        world: &mut crate::engine::World,
        id: ActorId,
        act: &Act,
    ) -> Outcome {
        let Some(actor) = self.get(id).cloned() else {
            return Outcome::Refused(Refusal::NoSuchActor);
        };
        let target = act.target();
        if !Self::alive(world, target) || !Self::alive(world, actor.body) {
            return Outcome::Refused(Refusal::NoSuchTarget);
        }
        // Reach is checked against the separation the engine computes, not
        // against anything in the request.
        let sep = world
            .tree
            .separation(actor.body, Vec3::ZERO, target, Vec3::ZERO)
            .value
            .norm();
        if sep > actor.reach {
            return Outcome::Refused(Refusal::OutOfReach);
        }

        match act {
            Act::Push { target, force, seconds } => {
                if !force.is_finite() || !seconds.is_finite() {
                    return Outcome::Refused(Refusal::NotANumber);
                }
                let magnitude = force.norm();
                let (f, ff) = clamp_to(magnitude, actor.max_force);
                let (s, fs) = clamp_to(seconds.max(0.0), actor.max_duration);
                // Impulse is force times time, and both halves were checked, so
                // the product cannot exceed what this body could deliver.
                let dp = if magnitude > 0.0 {
                    force.scale(f / magnitude * s)
                } else {
                    Vec3::ZERO
                };
                world.interact(Interaction::Impulse { target: *target, dp });
                let fraction = ff * fs;
                if fraction >= 1.0 {
                    Outcome::Accepted
                } else {
                    Outcome::Clamped { fraction }
                }
            }
            Act::Heat { target, watts, seconds } => {
                if !watts.is_finite() || !seconds.is_finite() {
                    return Outcome::Refused(Refusal::NotANumber);
                }
                let (p, fp) = clamp_to(*watts, actor.max_power);
                let (s, fs) = clamp_to(seconds.max(0.0), actor.max_duration);
                let joules = p * s;
                if joules >= 0.0 {
                    world.interact(Interaction::Deposit { target: *target, joules, radius: 0.0 });
                } else {
                    world.interact(Interaction::Extract { target: *target, joules: -joules });
                }
                let fraction = fp * fs;
                if fraction >= 1.0 {
                    Outcome::Accepted
                } else {
                    Outcome::Clamped { fraction }
                }
            }
            Act::Release { target, kg_per_second, seconds, composition, velocity } => {
                if !kg_per_second.is_finite() || !seconds.is_finite() || !velocity.is_finite() {
                    return Outcome::Refused(Refusal::NotANumber);
                }
                let (r, fr) = clamp_to(kg_per_second.max(0.0), actor.max_flow);
                let (s, fs) = clamp_to(seconds.max(0.0), actor.max_duration);
                world.interact(Interaction::Inject {
                    target: *target,
                    mass: r * s,
                    composition: *composition,
                    velocity: *velocity,
                });
                let fraction = fr * fs;
                if fraction >= 1.0 {
                    Outcome::Accepted
                } else {
                    Outcome::Clamped { fraction }
                }
            }
            Act::Sense { target, instrument, quantity } => {
                if !actor.instruments.contains(instrument) {
                    return Outcome::Refused(Refusal::NoSuchInstrument);
                }
                if actor.observer.is_none() {
                    // Only reachable for an `Actor` assembled by hand rather
                    // than admitted. Sensing without senses is not a physics
                    // failure, it is a bookkeeping one.
                    return Outcome::Refused(Refusal::NoSuchInstrument);
                }
                match world.measure(*target, *instrument, *quantity) {
                    Some(r) => Outcome::Sensed(Box::new(r)),
                    None => Outcome::Refused(Refusal::NothingToSense),
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// the wire
// ---------------------------------------------------------------------------

const INSTRUMENTS: [Instrument; 6] = [
    Instrument::Imager,
    Instrument::Spectrometer,
    Instrument::Thermometer,
    Instrument::ParticleDetector,
    Instrument::Interferometer,
    Instrument::MassSpectrometer,
];
const QUANTITIES: [Quantity; 7] = [
    Quantity::DecayTime,
    Quantity::QuantumLevel,
    Quantity::Position,
    Quantity::Momentum,
    Quantity::Spin,
    Quantity::PhotonCount,
    Quantity::Temperature,
];
const PROPERTIES: [Property; 6] = [
    Property::Mass,
    Property::Temperature,
    Property::Radius,
    Property::Charge,
    Property::Luminosity,
    Property::TimeRate,
];

pub fn encode(c: &Command) -> Vec<u8> {
    let mut w = Writer::new();
    w.header();
    match c {
        Command::Act { actor, act } => {
            w.u8(0);
            w.u64(actor.0);
            match act {
                Act::Push { target, force, seconds } => {
                    w.u8(0);
                    w.u32(target.0);
                    w.vec3(*force);
                    w.f64(*seconds);
                }
                Act::Heat { target, watts, seconds } => {
                    w.u8(1);
                    w.u32(target.0);
                    w.f64(*watts);
                    w.f64(*seconds);
                }
                Act::Release { target, kg_per_second, seconds, composition, velocity } => {
                    w.u8(2);
                    w.u32(target.0);
                    w.f64(*kg_per_second);
                    w.f64(*seconds);
                    for v in composition.0 {
                        w.f64(v);
                    }
                    w.vec3(*velocity);
                }
                Act::Sense { target, instrument, quantity } => {
                    w.u8(3);
                    w.u32(target.0);
                    w.u8(INSTRUMENTS.iter().position(|x| x == instrument).unwrap_or(0) as u8);
                    w.u8(QUANTITIES.iter().position(|x| x == quantity).unwrap_or(0) as u8);
                }
            }
        }
        Command::Admin(a) => {
            w.u8(1);
            match a {
                Admin::Dilate { target, rate } => {
                    w.u8(0);
                    w.u32(target.0);
                    w.f64(*rate);
                }
                Admin::Author { target, property, value } => {
                    w.u8(1);
                    w.u32(target.0);
                    w.u8(PROPERTIES.iter().position(|x| x == property).unwrap_or(0) as u8);
                    w.f64(*value);
                }
                Admin::Pin { target } => {
                    w.u8(2);
                    w.u32(target.0);
                }
            }
        }
    }
    w.finish()
}

pub fn decode(bytes: &[u8]) -> Result<Command> {
    let mut r = Reader::new(bytes);
    r.header()?;
    let out = match r.tag("command", 2)? {
        0 => {
            let actor = ActorId(r.u64()?);
            let act = match r.tag("act", 4)? {
                0 => Act::Push {
                    target: NodeIdx(r.u32()?),
                    force: r.vec3()?,
                    seconds: r.f64()?,
                },
                1 => Act::Heat {
                    target: NodeIdx(r.u32()?),
                    watts: r.f64()?,
                    seconds: r.f64()?,
                },
                2 => {
                    let target = NodeIdx(r.u32()?);
                    let kg_per_second = r.f64()?;
                    let seconds = r.f64()?;
                    let mut c = [0.0f64; crate::units::NSPECIES];
                    for slot in c.iter_mut() {
                        *slot = r.f64()?;
                    }
                    Act::Release {
                        target,
                        kg_per_second,
                        seconds,
                        composition: Composition(c),
                        velocity: r.vec3()?,
                    }
                }
                _ => Act::Sense {
                    target: NodeIdx(r.u32()?),
                    instrument: INSTRUMENTS[r.tag("instrument", 6)? as usize],
                    quantity: QUANTITIES[r.tag("quantity", 7)? as usize],
                },
            };
            Command::Act { actor, act }
        }
        _ => Command::Admin(match r.tag("admin", 3)? {
            0 => Admin::Dilate { target: NodeIdx(r.u32()?), rate: r.f64()? },
            1 => Admin::Author {
                target: NodeIdx(r.u32()?),
                property: PROPERTIES[r.tag("property", 6)? as usize],
                value: r.f64()?,
            },
            _ => Admin::Pin { target: NodeIdx(r.u32()?) },
        }),
    };
    r.finish()?;
    Ok(out)
}
