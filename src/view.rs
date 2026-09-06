//! What a client is allowed to see.
//!
//! # The seam
//!
//! Until now the viewer *was* the engine: it held a `World`, called `refine` to
//! make detail, read `tree.nodes[i].bodies` directly, and reached into
//! `stats` for its readouts. That works exactly as long as everything lives in
//! one process, and stops working the moment the authoritative loop is
//! somewhere else — which is the whole direction of the architecture.
//!
//! So this module defines the boundary. A client sends a [`ViewRequest`] and
//! receives a [`Scene`], and a `Scene` is the *only* thing it gets. It contains
//! no solvers, no `Tree`, no `Body`, no aggregate — nothing that could be
//! stepped. A renderer holding one can draw a picture and cannot advance
//! physics, which is the property worth having.
//!
//! # Everything is node-relative
//!
//! Positions come back in units of the node's own radius. This is not a
//! convenience: it is what lets one renderer draw a galaxy and a nucleus with
//! the same code, and it means the numbers crossing the boundary are of order
//! one whether the node is fifteen kiloparsecs across or four femtometres. The
//! scale itself travels once, as `node.radius`, in full precision.
//!
//! # Already interpolated
//!
//! A node the last frame could not bring all the way to the instant is behind
//! by its lateness, and `render` carries its bodies forward at their own
//! velocities before handing them over — the same closed-form step
//! `World::render_lag` describes. The client is told the instant it is looking
//! at and does not have to know that some of the world was solved earlier.
//!
//! # It crosses a process boundary
//!
//! `Scene` encodes with the same wire format a world file uses, so the split is
//! real rather than notional: `phys-headless` writes scenes from one process and
//! reads them in another with no engine involved. That is also how the cost of
//! a client gets measured, which is the number Phase 3 needs.

use crate::ids::{NodeIdx, PathKey};
use crate::solvers::SolverKind;
use crate::units::Tier;
use crate::wire::{Reader, Result, Writer};

/// What a client is asking to look at.
#[derive(Debug, Clone, Copy)]
pub struct ViewRequest {
    /// The node to draw.
    pub node: NodeIdx,
    /// Largest number of bodies to send. The *client's* budget — a phone and a
    /// workstation ask the same world for different amounts. Zero means all.
    pub max_bodies: usize,
    /// How many ancestors of `node` to describe, for navigation context.
    pub trail: usize,
    /// Send bodies only if the node has been re-solved since this instant.
    ///
    /// This is what stops bandwidth scaling with how much world is on screen.
    /// A node that was coasted has bodies the client can carry forward itself —
    /// they each travel with their own velocity — so there is nothing new to
    /// send, and the answer is a few hundred bytes of facts instead of a body
    /// list. A client asks with the instant of its last update; a client that
    /// wants everything asks with negative infinity.
    ///
    /// It is the same idea as `WorldStore::save_since`, pointed the other way:
    /// traffic proportional to *events* rather than to *entities*.
    pub since: f64,
}

impl Default for ViewRequest {
    fn default() -> ViewRequest {
        ViewRequest { node: NodeIdx::NONE, max_bodies: 0, trail: 16, since: f64::NEG_INFINITY }
    }
}

impl ViewRequest {
    pub fn of(node: NodeIdx) -> ViewRequest {
        ViewRequest { node, ..Default::default() }
    }
}

/// One drawable thing.
///
/// `f32` throughout, deliberately. Position and radius are in units of the
/// node's radius so they are of order one; the three channel scalars are for
/// colouring and their true ranges travel separately in [`NodeFacts`]. The one
/// bound worth stating: a single body's mass must stay under `f32::MAX`
/// (3.4×10^38 kg). The largest the engine produces is a galaxy's super-particle
/// at around 10^35 kg, three orders of magnitude clear.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Speck {
    /// In units of the node's radius.
    pub pos: [f32; 3],
    /// Node radii per second. Carried so the client can advance a coasting
    /// node itself, which is what makes `ViewRequest::since` worth having.
    pub vel: [f32; 3],
    pub radius: f32,
    pub mass: f32,
    pub temperature: f32,
    pub kind: u8,
}

impl Speck {
    /// Node radii per second. The channel a renderer usually colours by.
    pub fn speed(&self) -> f32 {
        (self.vel[0] * self.vel[0] + self.vel[1] * self.vel[1] + self.vel[2] * self.vel[2]).sqrt()
    }

    /// Where this body will be `dt` seconds after the scene's instant.
    ///
    /// The same closed-form carry the engine performs internally, available to
    /// a client that is holding a scene older than the frame it is drawing.
    pub fn at(&self, dt: f32) -> [f32; 3] {
        [
            self.pos[0] + self.vel[0] * dt,
            self.pos[1] + self.vel[1] * dt,
            self.pos[2] + self.vel[2] * dt,
        ]
    }
}

/// The node being looked at, in full precision.
///
/// Full precision because these span the ladder: a galaxy weighs 10^39 kg and a
/// nucleon 10^-27, and an `f32` reports the first as infinity.
#[derive(Debug, Clone, Copy)]
pub struct NodeFacts {
    pub key: PathKey,
    pub tier: u8,
    pub solver: u8,
    pub depth: u32,
    pub mass: f64,
    pub radius: f64,
    pub temperature: f64,
    pub internal_energy: f64,
    pub binding_energy: f64,
    pub luminosity: f64,
    pub charge: f64,
    /// How often this node is re-solved, seconds.
    pub cadence: f64,
    /// The sub-step its physics needs when it is solved, seconds.
    pub timestep: f64,
    /// How long before its detail stops meaning anything. Infinite for
    /// structured matter, which never forgets.
    pub mixing_time: f64,
    /// How far behind the instant its last solve left it, seconds.
    pub lag: f64,
    /// Bodies the node actually holds, which may exceed the number sent.
    pub body_count: u32,
    /// How many of those belong to a structure rather than to loose matter.
    pub structural_parts: u32,
    pub materialised: bool,
    pub pinned: bool,
}

/// The world the node is in.
#[derive(Debug, Clone, Copy, Default)]
pub struct WorldFacts {
    pub time: f64,
    /// Simulated seconds the last frame covered.
    pub frame_span: f64,
    /// Fraction of the asked-for pace actually being sustained.
    pub time_throttle: f64,
    /// Worst lateness anywhere, in units of that node's own characteristic
    /// time. Under one means every node was re-solved before it had changed.
    pub worst_lateness: f64,
    pub live_nodes: u32,
    /// Nodes carried to the instant in closed form rather than solved.
    pub coasted: u32,
    pub materialised_bodies: u64,
    pub detail_bytes: u64,
    pub thermalised: u64,
    pub detail_debt: f64,
    pub frames: u64,
}

impl Default for NodeFacts {
    fn default() -> NodeFacts {
        NodeFacts {
            key: PathKey::ROOT,
            tier: 0,
            solver: 0,
            depth: 0,
            mass: 0.0,
            radius: 0.0,
            temperature: 0.0,
            internal_energy: 0.0,
            binding_energy: 0.0,
            luminosity: 0.0,
            charge: 0.0,
            cadence: f64::INFINITY,
            timestep: 0.0,
            mixing_time: f64::INFINITY,
            lag: 0.0,
            body_count: 0,
            structural_parts: 0,
            materialised: false,
            pinned: false,
        }
    }
}

/// One rung of the ladder above the node being watched.
#[derive(Debug, Clone, Copy, Default)]
pub struct TrailStep {
    pub tier: u8,
    pub radius: f64,
}

/// Everything a client gets.
#[derive(Debug, Clone, Default)]
pub struct Scene {
    /// The instant every position in here has been carried to.
    pub instant: f64,
    pub node: NodeFacts,
    pub world: WorldFacts,
    /// Nearest ancestor first.
    pub trail: Vec<TrailStep>,
    /// False when the node had not been re-solved since the instant the client
    /// asked about, so nothing was sent and it should carry what it has.
    pub bodies_included: bool,
    pub bodies: Vec<Speck>,
}

impl Scene {
    pub fn tier(&self) -> Tier {
        Tier::ALL[(self.node.tier as usize).min(Tier::ALL.len() - 1)]
    }

    /// Range of a channel across the bodies present, for labelling a ramp.
    /// Returns `(lo, hi)`; equal when there is nothing to range over.
    pub fn channel_range(&self, pick: fn(&Speck) -> f32) -> (f32, f32) {
        let mut lo = f32::INFINITY;
        let mut hi = f32::NEG_INFINITY;
        for b in &self.bodies {
            let v = pick(b);
            if v.is_finite() {
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        if lo > hi {
            (0.0, 0.0)
        } else {
            (lo, hi)
        }
    }
}

// ---------------------------------------------------------------------------
// the wire
// ---------------------------------------------------------------------------

/// A body on the wire: two bytes per position and velocity component, two for
/// the radius, one each for mass, temperature and kind.
const SPECK_BYTES: usize = 2 * 3 + 2 * 3 + 2 + 1 + 1 + 1;
const TRAIL_BYTES: usize = 1 + 8;

/// Smallest half-width the position range is allowed to collapse to, in node
/// radii. Guards against a divide-by-zero on a node whose bodies are all at the
/// origin.
const MIN_SPAN: f32 = 1e-6;

/// # Why quantising is safe here, and how far it goes
///
/// Positions are stored as `i16` across a span the scene *measures* rather than
/// assumes. A fixed range looked safe and was not: bodies do not sit inside
/// their node's nominal radius, because `prolong` samples profiles with tails —
/// a Plummer sphere puts outliers many radii out — and every one of them
/// clamped to the edge, collapsing the outer half of a galaxy onto a cube.
///
/// With the span measured, the step is `2 * span / 65536`. The question is
/// whether that is small against *the node's own resolution*: the engine never
/// claims to know where anything is more precisely than one resolution element,
/// which for `n` bodies is `radius / n^(1/3)`. The step is under one percent of
/// an element while
///
/// ```text
///     2 * span * n^(1/3) / 65536 <= 0.01
/// ```
///
/// which for a span of ten radii holds to `n` of about 35,000 — comfortably
/// past anything `prolong` builds. So the wire is two orders of magnitude finer
/// than the physics it describes, and quantisation is provably not the limiting
/// error. `quantisation_is_finer_than_the_physics` measures it rather than
/// trusting this arithmetic.
///
/// Past that count the honest response is to send fewer bodies rather than more
/// precision: a client that cannot resolve 35,000 specks is not helped by
/// knowing exactly where they are.
#[inline]
fn q_pos(v: f32, span: f32) -> i16 {
    let c = (v / span).clamp(-1.0, 1.0);
    (c * i16::MAX as f32) as i16
}
#[inline]
fn dq_pos(v: i16, span: f32) -> f32 {
    v as f32 / i16::MAX as f32 * span
}

/// Quantise into a range as an unsigned integer of `bits` bits.
#[inline]
fn q_range(v: f32, lo: f32, hi: f32, max: f32) -> u16 {
    if !(hi > lo) || !v.is_finite() {
        return 0;
    }
    let t = ((v - lo) / (hi - lo)).clamp(0.0, 1.0);
    (t * max) as u16
}
#[inline]
fn dq_range(q: u16, lo: f32, hi: f32, max: f32) -> f32 {
    if !(hi > lo) {
        return lo;
    }
    lo + (q as f32 / max) * (hi - lo)
}

/// Channel ranges travel once per scene so the per-body bytes can be small.
/// Logs where the quantity spans orders of magnitude, because a linear ramp
/// over a mass range of 10^6 puts every body in the first bucket.
struct Ranges {
    pos: f32,
    vel: f32,
    radius: (f32, f32),
    mass: (f32, f32),
    temperature: (f32, f32),
}

fn log1p_safe(v: f32) -> f32 {
    if v > 0.0 {
        v.ln()
    } else {
        0.0
    }
}

impl Ranges {
    fn of(bodies: &[Speck]) -> Ranges {
        let mut pos: f32 = 0.0;
        let mut vel: f32 = 0.0;
        let mut r = (f32::INFINITY, f32::NEG_INFINITY);
        let mut m = (f32::INFINITY, f32::NEG_INFINITY);
        let mut t = (f32::INFINITY, f32::NEG_INFINITY);
        for b in bodies {
            for c in b.pos {
                if c.is_finite() {
                    pos = pos.max(c.abs());
                }
            }
            for c in b.vel {
                if c.is_finite() {
                    vel = vel.max(c.abs());
                }
            }
            if b.radius.is_finite() {
                let lr = log1p_safe(b.radius);
                r = (r.0.min(lr), r.1.max(lr));
            }
            if b.mass.is_finite() {
                let lm = log1p_safe(b.mass);
                m = (m.0.min(lm), m.1.max(lm));
            }
            if b.temperature.is_finite() {
                t = (t.0.min(b.temperature), t.1.max(b.temperature));
            }
        }
        let fix = |p: (f32, f32)| if p.0 > p.1 { (0.0, 0.0) } else { p };
        Ranges {
            pos: pos.max(MIN_SPAN),
            vel: if vel > 0.0 { vel } else { 1.0 },
            radius: fix(r),
            mass: fix(m),
            temperature: fix(t),
        }
    }
}

pub fn encode(s: &Scene) -> Vec<u8> {
    let mut w = Writer::new();
    w.header();
    w.f64(s.instant);

    let n = &s.node;
    w.u128(n.key.0);
    w.u8(n.tier);
    w.u8(n.solver);
    w.u32(n.depth);
    for v in [
        n.mass,
        n.radius,
        n.temperature,
        n.internal_energy,
        n.binding_energy,
        n.luminosity,
        n.charge,
        n.cadence,
        n.timestep,
        n.mixing_time,
        n.lag,
    ] {
        w.f64(v);
    }
    w.u32(n.body_count);
    w.u32(n.structural_parts);
    w.bool(n.materialised);
    w.bool(n.pinned);

    let d = &s.world;
    for v in [d.time, d.frame_span, d.time_throttle, d.worst_lateness, d.detail_debt] {
        w.f64(v);
    }
    w.u32(d.live_nodes);
    w.u32(d.coasted);
    w.u64(d.materialised_bodies);
    w.u64(d.detail_bytes);
    w.u64(d.thermalised);
    w.u64(d.frames);

    w.seq(s.trail.len());
    for t in &s.trail {
        w.u8(t.tier);
        w.f64(t.radius);
    }

    w.bool(s.bodies_included);
    let r = Ranges::of(&s.bodies);
    w.u32(r.pos.to_bits());
    w.u32(r.vel.to_bits());
    for (lo, hi) in [r.radius, r.mass, r.temperature] {
        w.u32(lo.to_bits());
        w.u32(hi.to_bits());
    }
    w.seq(s.bodies.len());
    for b in &s.bodies {
        for c in b.pos {
            w.u16(q_pos(c, r.pos) as u16);
        }
        for c in b.vel {
            w.u16(q_pos(c, r.vel) as u16);
        }
        w.u16(q_range(log1p_safe(b.radius), r.radius.0, r.radius.1, u16::MAX as f32));
        w.u8(q_range(log1p_safe(b.mass), r.mass.0, r.mass.1, 255.0) as u8);
        w.u8(q_range(b.temperature, r.temperature.0, r.temperature.1, 255.0) as u8);
        w.u8(b.kind);
    }
    w.finish()
}

pub fn decode(bytes: &[u8]) -> Result<Scene> {
    let mut r = Reader::new(bytes);
    r.header()?;
    let instant = r.f64()?;

    let mut node = NodeFacts { key: PathKey(r.u128()?), ..Default::default() };
    node.tier = r.tag("tier", Tier::ALL.len() as u8)?;
    node.solver = r.u8()?;
    node.depth = r.u32()?;
    node.mass = r.f64()?;
    node.radius = r.f64()?;
    node.temperature = r.f64()?;
    node.internal_energy = r.f64()?;
    node.binding_energy = r.f64()?;
    node.luminosity = r.f64()?;
    node.charge = r.f64()?;
    node.cadence = r.f64()?;
    node.timestep = r.f64()?;
    node.mixing_time = r.f64()?;
    node.lag = r.f64()?;
    node.body_count = r.u32()?;
    node.structural_parts = r.u32()?;
    node.materialised = r.bool()?;
    node.pinned = r.bool()?;

    let world = WorldFacts {
        time: r.f64()?,
        frame_span: r.f64()?,
        time_throttle: r.f64()?,
        worst_lateness: r.f64()?,
        detail_debt: r.f64()?,
        live_nodes: r.u32()?,
        coasted: r.u32()?,
        materialised_bodies: r.u64()?,
        detail_bytes: r.u64()?,
        thermalised: r.u64()?,
        frames: r.u64()?,
    };

    let n = r.seq("trail", TRAIL_BYTES)?;
    let mut trail = Vec::with_capacity(n);
    for _ in 0..n {
        trail.push(TrailStep {
            tier: r.tag("trail tier", Tier::ALL.len() as u8)?,
            radius: r.f64()?,
        });
    }

    let bodies_included = r.bool()?;
    let pos_span = f32::from_bits(r.u32()?);
    let vel_scale = f32::from_bits(r.u32()?);
    let mut ranges = [(0.0f32, 0.0f32); 3];
    for slot in ranges.iter_mut() {
        *slot = (f32::from_bits(r.u32()?), f32::from_bits(r.u32()?));
    }
    let n = r.seq("specks", SPECK_BYTES)?;
    let mut bodies = Vec::with_capacity(n);
    for _ in 0..n {
        let pos = [
            dq_pos(r.u16()? as i16, pos_span),
            dq_pos(r.u16()? as i16, pos_span),
            dq_pos(r.u16()? as i16, pos_span),
        ];
        let vel = [
            dq_pos(r.u16()? as i16, vel_scale),
            dq_pos(r.u16()? as i16, vel_scale),
            dq_pos(r.u16()? as i16, vel_scale),
        ];
        let radius = dq_range(r.u16()?, ranges[0].0, ranges[0].1, u16::MAX as f32).exp();
        let mass = dq_range(r.u8()? as u16, ranges[1].0, ranges[1].1, 255.0).exp();
        let temperature = dq_range(r.u8()? as u16, ranges[2].0, ranges[2].1, 255.0);
        bodies.push(Speck { pos, vel, radius, mass, temperature, kind: r.u8()? });
    }

    r.finish()?;
    Ok(Scene { instant, node, world, trail, bodies_included, bodies })
}

// ---------------------------------------------------------------------------
// producing one
// ---------------------------------------------------------------------------

impl crate::engine::World {
    /// Draw a scene.
    ///
    /// Takes `&self`. That is the whole point: rendering cannot materialise,
    /// cannot step, and cannot change the world in any way. A client that wants
    /// finer detail asks for it as a *command* — which is the other direction
    /// across this boundary, and stays explicit.
    pub fn render(&self, req: &ViewRequest) -> Scene {
        let mut scene = Scene { instant: self.time, ..Default::default() };

        scene.world = WorldFacts {
            time: self.time,
            frame_span: self.frame_dt(),
            time_throttle: self.time_throttle,
            worst_lateness: self.stats.worst_lateness,
            detail_debt: self.stats.detail_debt,
            live_nodes: self.tree.live_count() as u32,
            coasted: self.stats.coasted as u32,
            materialised_bodies: self.tree.materialised_bodies() as u64,
            detail_bytes: self.tree.detail_bytes() as u64,
            thermalised: self.stats.thermalised,
            frames: self.stats.frames,
        };

        let idx = req.node;
        if idx.is_none() || idx.get() >= self.tree.nodes.len() || !self.tree.nodes[idx.get()].alive
        {
            return scene;
        }
        let n = &self.tree.nodes[idx.get()];

        scene.node = NodeFacts {
            key: n.key,
            tier: n.tier.index() as u8,
            solver: crate::solvers::for_tier(n.tier) as u8,
            depth: n.depth,
            mass: n.agg.mass,
            radius: n.agg.radius,
            temperature: n.agg.temperature,
            internal_energy: n.agg.internal_energy,
            binding_energy: n.agg.binding_energy,
            luminosity: n.agg.luminosity,
            charge: n.agg.charge,
            cadence: self.node_cadence(idx),
            timestep: self.node_dt(idx),
            mixing_time: self.mixing_time(idx),
            lag: self.render_lag(idx),
            body_count: n.bodies.len() as u32,
            structural_parts: n.last_report.structural_parts as u32,
            materialised: n.is_materialised(),
            pinned: n.pinned,
        };

        // The ladder above, nearest first.
        let mut cur = n.parent;
        while !cur.is_none() && scene.trail.len() < req.trail {
            let a = &self.tree.nodes[cur.get()];
            scene.trail.push(TrailStep { tier: a.tier.index() as u8, radius: a.agg.radius });
            cur = a.parent;
        }

        // Bodies are sent only when there is something new to say about them.
        // A node that has not been re-solved since the client last asked has
        // bodies the client can carry forward itself, using the velocities it
        // already holds — so the answer is the facts above and nothing else.
        scene.bodies_included = n.last_solved > req.since;
        if !scene.bodies_included {
            return scene;
        }

        // Carried to the instant and scaled to the node.
        let lag = self.render_lag(idx);
        let inv = 1.0 / n.agg.radius.max(1e-300);
        let want = if req.max_bodies == 0 { n.bodies.len() } else { req.max_bodies };
        // Take every k-th rather than the first k: a prefix of a sampled
        // population is not a sample of it, and the first thousand parcels of a
        // disc are all in one place.
        let stride = if want == 0 { 1 } else { n.bodies.len().div_ceil(want.max(1)).max(1) };
        scene.bodies.reserve(n.bodies.len().div_ceil(stride));
        for b in n.bodies.iter().step_by(stride) {
            let p = b.pos + b.vel.scale(lag);
            scene.bodies.push(Speck {
                pos: [(p.x * inv) as f32, (p.y * inv) as f32, (p.z * inv) as f32],
                vel: [
                    (b.vel.x * inv) as f32,
                    (b.vel.y * inv) as f32,
                    (b.vel.z * inv) as f32,
                ],
                radius: (b.radius * inv) as f32,
                mass: b.mass as f32,
                temperature: b.temperature as f32,
                kind: b.kind as u8,
            });
        }
        scene
    }
}

/// Solver names, for a client that has no engine to ask.
pub const SOLVER_NAMES: [&str; 5] =
    ["Barnes-Hut", "gravity + SPH", "SPH", "molecular dynamics", "statistical"];

pub fn solver_name(tag: u8) -> &'static str {
    SOLVER_NAMES.get(tag as usize).copied().unwrap_or("?")
}

/// So a client can name a tier without linking the engine's units.
pub const TIER_NAMES: [&str; 7] = [
    "galactic",
    "stellar",
    "planetary",
    "continuum",
    "molecular",
    "atomic",
    "nuclear",
];

pub fn tier_name(tag: u8) -> &'static str {
    TIER_NAMES.get(tag as usize).copied().unwrap_or("?")
}

/// Body-kind names, likewise.
pub const KIND_NAMES: [&str; 12] = [
    "super-particle",
    "star",
    "compact object",
    "planet",
    "gas parcel",
    "grain",
    "molecule",
    "atom",
    "nucleus",
    "nucleon",
    "electron",
    "photon",
];

pub fn kind_name(tag: u8) -> &'static str {
    KIND_NAMES.get(tag as usize).copied().unwrap_or("?")
}

impl SolverKind {
    pub fn tag(self) -> u8 {
        self as u8
    }
}
