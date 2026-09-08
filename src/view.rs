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
//! receives a [`Scene`], and a `Scene` is the *only* thing it gets. It carries
//! no solver and no clock: a renderer holding one can draw a picture and cannot
//! advance physics, which is the property worth having. Note what that
//! guarantee is *not*. A [`Recipe`] does hand the client a node's matter, because
//! the client needs one to generate its own scenery — the invariant is
//! "cannot step time", not "cannot see state", and it survives intact, since
//! [`crate::sampler`] samples an instant and has no time argument to give it.
//!
//! # Everything is node-relative
//!
//! Positions come back in units of the node's own radius. This is not a
//! convenience: it is what lets one renderer draw a galaxy and a nucleus with
//! the same code, and it means the numbers crossing the boundary are of order
//! one whether the node is fifteen kiloparsecs across or four femtometres. The
//! scale itself travels once, as `node.radius`, in full precision.
//!
//! When a scene holds several nodes, each one's `offset` and `scale` are in
//! radii of `nodes[0]` — the *frame* node the request was rooted at — so a
//! renderer places them all in one space and then scales each node's specks by
//! its own `scale`.
//!
//! # Already interpolated
//!
//! A node the last frame could not bring all the way to the instant is behind
//! by its lateness, and `render` carries its bodies forward at their own
//! velocities before handing them over — the same closed-form step
//! `World::render_lag` describes. The client is told the instant it is looking
//! at and does not have to know that some of the world was solved earlier.
//!
//! # Three things that stop bandwidth tracking scene complexity
//!
//! A city is the hard case: ten thousand things on screen, almost none of them
//! doing anything. Sending every body every frame is quadratic nonsense, so
//! three mechanisms between them make traffic track *events* and *screen*
//! rather than *entities* and *world*.
//!
//! 1. **Recipes.** Detail the server does not itself hold is not sent as
//!    bodies; the ~300 bytes that would *generate* those bodies is sent
//!    instead, and the client runs the same deterministic sampler. See
//!    [`Recipe`].
//! 2. **Volume queries with distance LOD.** [`Volume`] asks for everything
//!    near an eye and gives each node detail in proportion to the angle it
//!    subtends. A building a pixel wide comes back as bulk matter. Cost is
//!    bounded by the screen, not by the world.
//! 3. **Per-client deltas.** A [`Client`] remembers what it sent and sends only
//!    what changed, plus the keys that left the query. See [`Detail::Unchanged`]
//!    and [`Scene::dropped`].
//!
//! # It crosses a process boundary
//!
//! `Scene` encodes with the same wire format a world file uses, so the split is
//! real rather than notional: `phys-headless` writes scenes from one process and
//! reads them in another with no engine involved. That is also how the cost of
//! a client gets measured, which is the number Phase 3 needs.

use std::collections::HashMap;

use crate::ids::{NodeIdx, PathKey};
use crate::math::{Quat, Vec3};
use crate::solvers::SolverKind;
use crate::units::Tier;
use crate::wire::{Reader, Result, Writer};

/// A sphere to look inside, and how finely.
///
/// Coordinates are in radii of the request's frame node, the same units
/// everything else in a [`Scene`] uses.
#[derive(Debug, Clone, Copy)]
pub struct Volume {
    /// Where the client is looking from.
    pub eye: [f64; 3],
    /// How far out to gather nodes.
    pub reach: f64,
    /// Half-angle, radians, below which a node is described but not detailed.
    ///
    /// This is the whole LOD rule. A node's half-angle is `radius / distance`,
    /// so the threshold is a statement about pixels: at a 60° field of view
    /// across 1080 pixels, one pixel is about 1e-3 radians, and anything under
    /// that cannot be drawn as more than a dot however many bodies it has.
    pub detail_angle: f64,
    /// Most nodes to return. The query stops here, largest angle first, so the
    /// cap costs the things furthest away.
    pub max_nodes: usize,
}

impl Default for Volume {
    fn default() -> Volume {
        Volume { eye: [0.0; 3], reach: f64::INFINITY, detail_angle: 1e-3, max_nodes: 256 }
    }
}

/// What a client is asking to look at.
#[derive(Debug, Clone)]
pub struct ViewRequest {
    /// The node to draw, and the frame everything else is expressed in.
    pub node: NodeIdx,
    /// Largest number of bodies to send. The *client's* budget — a phone and a
    /// workstation ask the same world for different amounts. Zero means all.
    ///
    /// In a volume query this is the budget for the whole scene, shared out by
    /// angular size, so a client's cost per frame is a number it chooses rather
    /// than a number the world hands it.
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
    /// traffic proportional to *events* rather than to *entities*. A [`Client`]
    /// tracks this per node instead, which is strictly better; this field is
    /// the floor under it for a client that keeps no state.
    pub since: f64,
    /// Look inside a sphere rather than at one node.
    pub volume: Option<Volume>,
    /// Whether the client can run [`Recipe::build`].
    ///
    /// A client that cannot — a thin renderer with no sampler linked — sets
    /// this false and gets bodies or nothing. Nothing is the honest answer for
    /// an unmaterialised node: the server does not have its bodies either, and
    /// making them to answer a *view* would be the view mutating the world.
    pub allow_recipes: bool,
}

impl Default for ViewRequest {
    fn default() -> ViewRequest {
        ViewRequest {
            node: NodeIdx::NONE,
            max_bodies: 0,
            trail: 16,
            since: f64::NEG_INFINITY,
            volume: None,
            allow_recipes: true,
        }
    }
}

impl ViewRequest {
    pub fn of(node: NodeIdx) -> ViewRequest {
        ViewRequest { node, ..Default::default() }
    }

    /// Look inside a sphere around `eye`, in frame-node radii.
    pub fn within(node: NodeIdx, eye: [f64; 3], reach: f64) -> ViewRequest {
        ViewRequest {
            node,
            volume: Some(Volume { eye, reach, ..Default::default() }),
            ..Default::default()
        }
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
    /// node itself, which is what makes [`ViewRequest::since`] worth having.
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

/// Turn engine bodies into specks, in node-relative units, carried by `lag`.
fn specks_of(bodies: &[crate::state::Body], radius: f64, lag: f64, stride: usize) -> Vec<Speck> {
    let inv = 1.0 / radius.max(1e-300);
    let stride = stride.max(1);
    let mut out = Vec::with_capacity(bodies.len().div_ceil(stride));
    for b in bodies.iter().step_by(stride) {
        let p = b.pos + b.vel.scale(lag);
        out.push(Speck {
            pos: [(p.x * inv) as f32, (p.y * inv) as f32, (p.z * inv) as f32],
            vel: [(b.vel.x * inv) as f32, (b.vel.y * inv) as f32, (b.vel.z * inv) as f32],
            radius: (b.radius * inv) as f32,
            mass: b.mass as f32,
            temperature: b.temperature as f32,
            kind: b.kind as u8,
        });
    }
    out
}

/// Take every `stride`-th of `want` from `have`.
fn stride_for(have: usize, want: usize) -> usize {
    if want == 0 || have == 0 {
        1
    } else {
        have.div_ceil(want.max(1)).max(1)
    }
}

// ---------------------------------------------------------------------------
// recipes — detail as instructions rather than as bodies
// ---------------------------------------------------------------------------

/// How to build a node's detail, instead of the detail itself.
///
/// # Why this is a large win and not a cheat
///
/// A node the server has not materialised has no bodies *on the server either*.
/// Today such a node renders as an empty scene: the client is told a building
/// exists and given nothing to draw, because the only way to get bodies would
/// be for `render` to call `refine` — a view mutating the world, which the
/// `&self` on [`crate::engine::World::render`] exists to forbid.
///
/// A recipe is the way out. [`crate::sampler`] is deterministic in
/// `(matter, spec, world_seed, path_key, epoch)` and nothing else, so those
/// ~300 bytes *are* the bodies, losslessly, whatever the count. Ten thousand
/// procedural buildings cost three megabytes of recipes once instead of a
/// hundred and twenty megabytes of bodies every time one comes into view, and
/// the server pays nothing at all: it never built them.
///
/// # What it is never used for
///
/// Only for detail the server does not hold. A node that is materialised has
/// bodies that have since been *stepped*, and stepped bodies are not what the
/// sampler would draw from the current matter — so those are sent
/// explicitly. A node that is pinned or has stored detail was altered by an
/// interaction and by definition cannot be regenerated; those are explicit too.
///
/// That rule is also the answer to the obvious worry about float determinism
/// across machines. A recipe is used exactly where there is no server-side
/// truth for the client's version to disagree with, so two clients drawing the
/// same untouched scenery a few ulps apart is a difference nobody can observe.
/// The moment anything *happens* there, the node materialises and explicit
/// bodies take over.
///
/// # What is checked
///
/// [`Recipe::build`] verifies that the mass it generated sums to the matter
/// mass it was given. That catches the failure that matters — a blob from a
/// different build of the engine, decoded into plausible-looking nonsense —
/// and costs no extra bytes on the wire. `checksum` catches the cruder case of
/// a mangled blob before it is decoded at all.
#[derive(Debug, Clone, PartialEq)]
pub struct Recipe {
    pub key: PathKey,
    pub epoch: u32,
    pub seed: u64,
    /// Bodies this will produce, before `stride`.
    pub count: u32,
    /// Take every `stride`-th, to honour the client's body budget.
    pub stride: u32,
    /// Seconds to carry the generated bodies forward by, so a recipe lands at
    /// the scene's instant like everything else.
    pub lag: f64,
    /// FNV-1a over `blob`.
    pub checksum: u64,
    /// Matter, spec and morphology, in the crate's own encoding. Opaque
    /// here on purpose: the client opens it through [`Recipe::build`] and the
    /// wire never has to know what is in it.
    blob: Vec<u8>,
}

/// Largest recipe blob accepted from the wire. A node's matter and a spec are
/// about 290 bytes; a morphology adds its event log. Well clear of both, and
/// far under anything that could be used to make a client allocate.
const MAX_RECIPE_BYTES: usize = 1 << 16;

fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        h ^= *b as u64;
        h = h.wrapping_mul(0x1000_0000_01b3);
    }
    h
}

impl Recipe {
    /// Bytes this recipe occupies, excluding its fixed header.
    pub fn bytes(&self) -> usize {
        self.blob.len()
    }

    /// Run the sampler and produce the specks the server would have sent.
    ///
    /// Node-relative, carried to the scene's instant, strided to the client's
    /// budget — indistinguishable from [`Detail::Explicit`] once built.
    pub fn build(&self) -> Result<Vec<Speck>> {
        if fnv1a(&self.blob) != self.checksum {
            return Err(crate::wire::WireError::BadRecipe { what: "checksum" });
        }
        let mut r = Reader::new(&self.blob);
        let matter = crate::persist::get_matter(&mut r)?;
        let spec = crate::persist::get_spec(&mut r)?;
        let morph = if r.bool()? { Some(crate::persist::get_morphology(&mut r)?) } else { None };
        r.finish()?;

        let bodies = match &morph {
            Some(m) => {
                crate::sampler::sample_structured(&matter, m, spec.count, self.seed, self.key.0, self.epoch).0
            }
            None => crate::sampler::sample(&matter, spec, self.seed, self.key.0, self.epoch).0,
        };

        // The check that is worth making. If this build of the engine samples
        // differently from the one that wrote the recipe, the conserved
        // quantity it was built against is the first thing to disagree.
        let got: f64 = bodies.iter().map(|b| b.mass).sum();
        if matter.mass > 0.0 && ((got - matter.mass) / matter.mass).abs() > 1e-9 {
            return Err(crate::wire::WireError::BadRecipe { what: "mass" });
        }

        Ok(specks_of(&bodies, matter.radius, self.lag, self.stride as usize))
    }
}

/// What a node's detail arrived as.
#[derive(Debug, Clone, PartialEq)]
pub enum Detail {
    /// Nothing new. Either the node has not been re-solved since the client
    /// last heard about it, or it is too small on screen to be worth detailing.
    /// A client holding bodies for it should carry them forward; a client
    /// holding none should draw it as bulk matter.
    Unchanged,
    /// Bodies, quantised.
    Explicit(Vec<Speck>),
    /// Instructions to make them. See [`Recipe`].
    Recipe(Recipe),
}

impl Detail {
    pub fn is_unchanged(&self) -> bool {
        matches!(self, Detail::Unchanged)
    }
    /// Specks, if they are already here. A recipe returns nothing until built.
    pub fn specks(&self) -> &[Speck] {
        match self {
            Detail::Explicit(b) => b,
            _ => &[],
        }
    }
}

// ---------------------------------------------------------------------------
// facts
// ---------------------------------------------------------------------------

/// A node, in full precision.
///
/// Full precision because these span the ladder: a galaxy weighs 10^39 kg and a
/// nucleon 10^-27, and an `f32` reports the first as infinity.
#[derive(Debug, Clone, Copy)]
pub struct NodeFacts {
    pub key: PathKey,
    /// Bumped whenever a recorded interaction changed what is inside this
    /// node. Detail generated at an older epoch is *gone*, not merely stale,
    /// which is why a client tracks it separately from the instant.
    pub epoch: u32,
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
    /// Nodes whose own physics needs a finer step than the frame can afford,
    /// and which cannot be crossed by ensemble because something is watching
    /// them. Distinct from lateness, which recovers: this does not. Non-zero
    /// means the world is running at a pace something in it cannot be
    /// integrated at, and it will still be non-zero next frame.
    pub unreachable: u64,
    pub detail_debt: f64,
    pub frames: u64,
}

impl Default for NodeFacts {
    fn default() -> NodeFacts {
        NodeFacts {
            key: PathKey::ROOT,
            epoch: 0,
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

/// One node in a scene, placed in the frame node's coordinates.
#[derive(Debug, Clone)]
pub struct NodeView {
    pub facts: NodeFacts,
    /// False when the server knew the client already had these facts and sent
    /// only the key. Everything but `facts.key` is then default-valued, and a
    /// client that keeps state should merge in the copy it holds.
    ///
    /// This matters more than it sounds. A node's facts are 141 bytes and its
    /// placement is 21; for a neighbourhood of two hundred things that nobody
    /// is touching, re-sending the facts every frame *is* the steady-state
    /// bill, and [`Detail::Unchanged`] alone does nothing about it.
    ///
    /// One field is excluded from the comparison: `facts.lag` is a per-frame
    /// derivative rather than a property, so a held copy's `lag` is as of the
    /// last frame that carried facts.
    pub facts_included: bool,
    /// Centre, in radii of the frame node. Zero for the frame node itself.
    pub offset: [f32; 3],
    /// This node's radius, in radii of the frame node. One for the frame node.
    /// Multiply a speck's `pos` and `radius` by this to place it in the frame.
    pub scale: f32,
    /// Half-angle subtended at the eye, radians. Zero when the request named
    /// no volume, because then there is no eye to subtend at.
    pub angular_size: f32,
    pub detail: Detail,
}

impl Default for NodeView {
    fn default() -> NodeView {
        NodeView {
            facts: NodeFacts::default(),
            facts_included: true,
            offset: [0.0; 3],
            scale: 1.0,
            angular_size: 0.0,
            detail: Detail::Unchanged,
        }
    }
}

/// Everything a client gets.
#[derive(Debug, Clone, Default)]
pub struct Scene {
    /// The instant every position in here has been carried to.
    pub instant: f64,
    pub world: WorldFacts,
    /// Nearest ancestor of the frame node first.
    pub trail: Vec<TrailStep>,
    /// The frame node first, then whatever else the volume query gathered,
    /// largest on screen first.
    pub nodes: Vec<NodeView>,
    /// Nodes the client was holding that are no longer in the query — it has
    /// walked away from them, or they were culled. The client should forget
    /// them. Empty unless the request came through a [`Client`].
    pub dropped: Vec<PathKey>,
}

impl Scene {
    /// The node the request was rooted at, and the frame everything is in.
    pub fn frame(&self) -> Option<&NodeView> {
        self.nodes.first()
    }

    /// Facts about the frame node. A default-valued `NodeFacts` when the
    /// request named nothing that exists, which is what an empty scene is.
    pub fn node(&self) -> NodeFacts {
        self.nodes.first().map(|n| n.facts).unwrap_or_default()
    }

    /// The frame node's specks. Empty until [`Scene::materialise`] has run, if
    /// they arrived as a recipe.
    pub fn bodies(&self) -> &[Speck] {
        self.nodes.first().map(|n| n.detail.specks()).unwrap_or(&[])
    }

    /// False when the frame node had nothing new to say and the client should
    /// carry forward what it already holds.
    pub fn bodies_included(&self) -> bool {
        self.nodes.first().is_some_and(|n| !n.detail.is_unchanged())
    }

    pub fn tier(&self) -> Tier {
        Tier::ALL[(self.node().tier as usize).min(Tier::ALL.len() - 1)]
    }

    /// Turn every recipe in the scene into the specks it describes.
    ///
    /// The client's half of the first bandwidth mechanism, and the only place
    /// that runs the sampler outside the engine. Returns how many specks were
    /// built. On failure the scene is left exactly as it was — a partially
    /// expanded scene would be worse than an unexpanded one — and the caller
    /// should re-ask with `allow_recipes: false`.
    pub fn materialise(&mut self) -> Result<usize> {
        let mut built = Vec::with_capacity(self.nodes.len());
        for (i, n) in self.nodes.iter().enumerate() {
            if let Detail::Recipe(r) = &n.detail {
                built.push((i, r.build()?));
            }
        }
        let total = built.iter().map(|(_, b)| b.len()).sum();
        for (i, b) in built {
            self.nodes[i].detail = Detail::Explicit(b);
        }
        Ok(total)
    }

    /// Fill in facts the server left out, and remember the ones it sent.
    ///
    /// The client's half of the facts elision, and the mirror of
    /// [`Scene::materialise`]. A client that keeps state between frames must
    /// call this before reading `facts` on anything, because a node the server
    /// judged unchanged arrives as a bare key. Returns how many were filled in.
    ///
    /// `known` is the client's own table and is updated in place; a client that
    /// drops a node should drop it from here too — [`Scene::dropped`] says
    /// which, and `forget_dropped` does it.
    pub fn merge_facts(&mut self, known: &mut HashMap<PathKey, NodeFacts>) -> usize {
        let mut filled = 0;
        for n in self.nodes.iter_mut() {
            if n.facts_included {
                known.insert(n.facts.key, n.facts);
            } else if let Some(f) = known.get(&n.facts.key) {
                // `lag` was never in the comparison, so the held copy's is as
                // of the last frame that carried facts. Say so by zeroing it
                // rather than reporting a stale number as current.
                n.facts = NodeFacts { lag: 0.0, ..*f };
                filled += 1;
            }
        }
        filled
    }

    /// Drop from a client-side table everything this scene says to forget.
    pub fn forget_dropped<T>(&self, known: &mut HashMap<PathKey, T>) {
        for k in &self.dropped {
            known.remove(k);
        }
    }

    /// Specks in the scene, whichever node they belong to.
    pub fn speck_count(&self) -> usize {
        self.nodes.iter().map(|n| n.detail.specks().len()).sum()
    }

    /// Range of a channel across the frame node's bodies, for labelling a ramp.
    /// Returns `(lo, hi)`; equal when there is nothing to range over.
    pub fn channel_range(&self, pick: fn(&Speck) -> f32) -> (f32, f32) {
        let mut lo = f32::INFINITY;
        let mut hi = f32::NEG_INFINITY;
        for b in self.bodies() {
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
/// A node view with `Detail::Unchanged` and nothing else: facts, placement,
/// and a one-byte tag. Everything else only makes it larger.
const NODEVIEW_MIN_BYTES: usize = 1 + 16 + 12 + 4 + 4 + 1;
const KEY_BYTES: usize = 16;

/// Smallest half-width the position range is allowed to collapse to, in node
/// radii. Guards against a divide-by-zero on a node whose bodies are all at the
/// origin.
const MIN_SPAN: f32 = 1e-6;

/// # Why quantising is safe here, and how far it goes
///
/// Positions are stored as `i16` across a span the scene *measures* rather than
/// assumes. A fixed range looked safe and was not: bodies do not sit inside
/// their node's nominal radius, because `sample` samples profiles with tails —
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
/// past anything `sample` builds. So the wire is two orders of magnitude finer
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

/// Channel ranges travel once per node so the per-body bytes can be small.
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

/// Identity of a node's facts, for deciding whether to send them again.
///
/// `lag` is left out on purpose: it moves every frame on a node that is
/// otherwise doing nothing, and hashing it would mean the facts were never
/// once judged unchanged.
fn facts_hash(n: &NodeFacts) -> u64 {
    let mut w = Writer::new();
    put_facts(&mut w, &NodeFacts { lag: 0.0, ..*n });
    fnv1a(&w.finish())
}

fn put_facts(w: &mut Writer, n: &NodeFacts) {
    w.u128(n.key.0);
    w.u32(n.epoch);
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
}

fn get_facts(r: &mut Reader) -> Result<NodeFacts> {
    let mut n = NodeFacts { key: PathKey(r.u128()?), ..Default::default() };
    n.epoch = r.u32()?;
    n.tier = r.tag("tier", Tier::ALL.len() as u8)?;
    n.solver = r.u8()?;
    n.depth = r.u32()?;
    n.mass = r.f64()?;
    n.radius = r.f64()?;
    n.temperature = r.f64()?;
    n.internal_energy = r.f64()?;
    n.binding_energy = r.f64()?;
    n.luminosity = r.f64()?;
    n.charge = r.f64()?;
    n.cadence = r.f64()?;
    n.timestep = r.f64()?;
    n.mixing_time = r.f64()?;
    n.lag = r.f64()?;
    n.body_count = r.u32()?;
    n.structural_parts = r.u32()?;
    n.materialised = r.bool()?;
    n.pinned = r.bool()?;
    Ok(n)
}

fn put_specks(w: &mut Writer, bodies: &[Speck]) {
    let r = Ranges::of(bodies);
    w.u32(r.pos.to_bits());
    w.u32(r.vel.to_bits());
    for (lo, hi) in [r.radius, r.mass, r.temperature] {
        w.u32(lo.to_bits());
        w.u32(hi.to_bits());
    }
    w.seq(bodies.len());
    for b in bodies {
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
}

fn get_specks(r: &mut Reader) -> Result<Vec<Speck>> {
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
    Ok(bodies)
}

pub fn encode(s: &Scene) -> Vec<u8> {
    let mut w = Writer::new();
    w.header();
    w.f64(s.instant);

    let d = &s.world;
    for v in [d.time, d.frame_span, d.time_throttle, d.worst_lateness, d.detail_debt] {
        w.f64(v);
    }
    w.u32(d.live_nodes);
    w.u32(d.coasted);
    w.u64(d.materialised_bodies);
    w.u64(d.detail_bytes);
    w.u64(d.thermalised);
    w.u64(d.unreachable);
    w.u64(d.frames);

    w.seq(s.trail.len());
    for t in &s.trail {
        w.u8(t.tier);
        w.f64(t.radius);
    }

    w.seq(s.nodes.len());
    for v in &s.nodes {
        w.bool(v.facts_included);
        if v.facts_included {
            put_facts(&mut w, &v.facts);
        } else {
            w.u128(v.facts.key.0);
        }
        for c in v.offset {
            w.u32(c.to_bits());
        }
        w.u32(v.scale.to_bits());
        w.u32(v.angular_size.to_bits());
        match &v.detail {
            Detail::Unchanged => w.u8(0),
            Detail::Explicit(b) => {
                w.u8(1);
                put_specks(&mut w, b);
            }
            Detail::Recipe(r) => {
                w.u8(2);
                w.u128(r.key.0);
                w.u32(r.epoch);
                w.u64(r.seed);
                w.u32(r.count);
                w.u32(r.stride);
                w.f64(r.lag);
                w.u64(r.checksum);
                w.bytes(&r.blob);
            }
        }
    }

    w.seq(s.dropped.len());
    for k in &s.dropped {
        w.u128(k.0);
    }
    w.finish()
}

pub fn decode(bytes: &[u8]) -> Result<Scene> {
    let mut r = Reader::new(bytes);
    r.header()?;
    let instant = r.f64()?;

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
        unreachable: r.u64()?,
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

    let n = r.seq("nodes", NODEVIEW_MIN_BYTES)?;
    let mut nodes = Vec::with_capacity(n);
    for _ in 0..n {
        let facts_included = r.bool()?;
        let facts = if facts_included {
            get_facts(&mut r)?
        } else {
            NodeFacts { key: PathKey(r.u128()?), ..Default::default() }
        };
        let offset = [
            f32::from_bits(r.u32()?),
            f32::from_bits(r.u32()?),
            f32::from_bits(r.u32()?),
        ];
        let scale = f32::from_bits(r.u32()?);
        let angular_size = f32::from_bits(r.u32()?);
        let detail = match r.tag("detail", 3)? {
            0 => Detail::Unchanged,
            1 => Detail::Explicit(get_specks(&mut r)?),
            _ => {
                let key = PathKey(r.u128()?);
                let epoch = r.u32()?;
                let seed = r.u64()?;
                let count = r.u32()?;
                let stride = r.u32()?;
                let lag = r.f64()?;
                let checksum = r.u64()?;
                let blob = r.bytes()?;
                if blob.len() > MAX_RECIPE_BYTES {
                    return Err(crate::wire::WireError::TooLong {
                        what: "recipe",
                        len: blob.len() as u64,
                        remaining: MAX_RECIPE_BYTES,
                    });
                }
                Detail::Recipe(Recipe { key, epoch, seed, count, stride, lag, checksum, blob })
            }
        };
        nodes.push(NodeView { facts, facts_included, offset, scale, angular_size, detail });
    }

    let n = r.seq("dropped", KEY_BYTES)?;
    let mut dropped = Vec::with_capacity(n);
    for _ in 0..n {
        dropped.push(PathKey(r.u128()?));
    }

    r.finish()?;
    Ok(Scene { instant, world, trail, nodes, dropped })
}

// ---------------------------------------------------------------------------
// per-client state
// ---------------------------------------------------------------------------

/// What one client is known to be holding for one node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Held {
    pub epoch: u32,
    /// Instant the detail it holds was solved at.
    pub at: f64,
    /// Hash of the facts it was last sent, so unchanged facts are not resent.
    pub facts: u64,
}

/// One connection's memory of what it has already been told.
///
/// # Why the server keeps this and not the client
///
/// [`ViewRequest::since`] is one instant for the whole request, so a client
/// watching a hundred nodes must ask with the oldest of them and re-receive the
/// ninety-nine that had not changed. Tracking per node fixes that, and it has
/// to live on the server because the *decision* — send or not — is made there.
///
/// It is a few dozen bytes per node the client can see, which is the price of
/// making traffic proportional to what happened rather than to what is
/// on screen. A client that disconnects takes its `Client` with it.
#[derive(Debug, Default)]
pub struct Client {
    held: HashMap<PathKey, Held>,
    /// Bytes this client has been sent, encoded.
    pub sent_bytes: u64,
    pub frames: u64,
    /// Nodes described but not detailed, because nothing had changed or they
    /// were too small on screen.
    pub unchanged: u64,
    pub explicit: u64,
    pub recipes: u64,
    pub dropped: u64,
}

impl Client {
    pub fn new() -> Client {
        Client::default()
    }

    /// Nodes this client is believed to hold detail for.
    pub fn holding(&self) -> usize {
        self.held.len()
    }

    pub fn holds(&self, key: PathKey) -> Option<Held> {
        self.held.get(&key).copied()
    }

    /// Throw away everything. The client reconnected, or fell far enough
    /// behind that starting again is cheaper than catching up.
    pub fn reset(&mut self) {
        self.held.clear();
    }

    /// Render one frame for this client, sending only what it does not have.
    ///
    /// The returned scene's [`Scene::dropped`] names the keys that left the
    /// query since the last frame, so the client can release them.
    pub fn frame(&mut self, world: &crate::engine::World, req: &ViewRequest) -> Scene {
        let mut scene = world.render_tracked(req, Some(&self.held));

        // Anything we believed it held that is not in this scene has gone out
        // of the query — walked away from, or culled by the node cap.
        let mut present: Vec<PathKey> = Vec::with_capacity(scene.nodes.len());
        for n in &scene.nodes {
            present.push(n.facts.key);
        }
        scene.dropped = self
            .held
            .keys()
            .copied()
            .filter(|k| !present.contains(k))
            .collect();
        scene.dropped.sort_by_key(|k| k.0);
        for k in &scene.dropped {
            self.held.remove(k);
        }
        self.dropped += scene.dropped.len() as u64;

        for n in &scene.nodes {
            match &n.detail {
                Detail::Unchanged => {
                    self.unchanged += 1;
                    match self.held.get_mut(&n.facts.key) {
                        // Keep holding what it holds; only the instant moves
                        // on, and only for a node it has detail for.
                        Some(h) => {
                            h.at = h.at.max(scene.instant - n.facts.lag);
                            if n.facts_included {
                                h.facts = facts_hash(&n.facts);
                                h.epoch = n.facts.epoch;
                            }
                        }
                        // A node too small on screen to detail is still a node
                        // the client now knows about, and remembering that it
                        // has the facts is what stops them being re-sent every
                        // frame for the rest of the neighbourhood's life.
                        None => {
                            self.held.insert(
                                n.facts.key,
                                Held {
                                    epoch: n.facts.epoch,
                                    at: f64::NEG_INFINITY,
                                    facts: facts_hash(&n.facts),
                                },
                            );
                        }
                    }
                }
                Detail::Explicit(_) | Detail::Recipe(_) => {
                    if matches!(n.detail, Detail::Recipe(_)) {
                        self.recipes += 1;
                    } else {
                        self.explicit += 1;
                    }
                    // A node whose detail is being sent always has its facts
                    // sent with it, so the stored hash is the one just written.
                    self.held.insert(
                        n.facts.key,
                        Held {
                            epoch: n.facts.epoch,
                            at: scene.instant,
                            facts: facts_hash(&n.facts),
                        },
                    );
                }
            }
        }

        self.frames += 1;
        self.sent_bytes += encode(&scene).len() as u64;
        scene
    }
}

// ---------------------------------------------------------------------------
// producing a scene
// ---------------------------------------------------------------------------

/// A node the gather step found, before it was decided what to send about it.
struct Candidate {
    idx: NodeIdx,
    /// Centre relative to the frame node's centre, metres, in the frame node's
    /// orientation.
    offset: Vec3,
    /// Half-angle at the eye, radians. Infinite when there is no volume, which
    /// is how "always detail this one" falls out without a special case.
    angle: f64,
}

impl crate::engine::World {
    /// Draw a scene.
    ///
    /// Takes `&self`. That is the whole point: rendering cannot materialise,
    /// cannot step, and cannot change the world in any way. A client that wants
    /// finer detail asks for it as a *command* — which is the other direction
    /// across this boundary, and stays explicit.
    pub fn render(&self, req: &ViewRequest) -> Scene {
        self.render_tracked(req, None)
    }

    /// As [`World::render`], but told what the client already holds.
    ///
    /// `held` is `None` for a stateless client, which falls back to the single
    /// [`ViewRequest::since`] instant for every node. Use [`Client`] rather
    /// than calling this directly; it is public so a transport can keep its own
    /// table if it wants to.
    pub fn render_tracked(
        &self,
        req: &ViewRequest,
        held: Option<&HashMap<PathKey, Held>>,
    ) -> Scene {
        let mut scene = Scene {
            instant: self.time,
            world: WorldFacts {
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
                unreachable: self.stats.unreachable,
                frames: self.stats.frames,
            },
            ..Default::default()
        };

        let idx = req.node;
        if idx.is_none() || idx.get() >= self.tree.nodes.len() || !self.tree.nodes[idx.get()].alive
        {
            return scene;
        }
        let frame = &self.tree.nodes[idx.get()];
        let frame_radius = frame.matter.radius.max(1e-300);

        // The ladder above, nearest first.
        let mut cur = frame.parent;
        while !cur.is_none() && scene.trail.len() < req.trail {
            let a = &self.tree.nodes[cur.get()];
            scene.trail.push(TrailStep { tier: a.tier.index() as u8, radius: a.matter.radius });
            cur = a.parent;
        }

        let found = self.gather(idx, frame_radius, req.volume.as_ref());
        let budgets = share_bodies(self, &found, req.max_bodies);

        scene.nodes.reserve(found.len());
        for (c, want) in found.iter().zip(budgets) {
            let n = &self.tree.nodes[c.idx.get()];
            let facts = self.facts_of(c.idx);
            let detail = self.detail_of(c.idx, req, held, want, c.angle);
            // Facts travel only when they are news. For a neighbourhood of
            // things nobody is touching this is the difference between a few
            // hundred bytes a frame and thirty kilobytes.
            let facts_included = match held {
                // Detail always travels with its facts. A client that is being
                // handed bodies is being handed the frame they belong in, and
                // the bytes are lost in the noise of the bodies themselves.
                _ if !detail.is_unchanged() => true,
                Some(h) => h.get(&n.key).is_none_or(|h| h.facts != facts_hash(&facts)),
                None => true,
            };
            scene.nodes.push(NodeView {
                facts,
                facts_included,
                offset: [
                    (c.offset.x / frame_radius) as f32,
                    (c.offset.y / frame_radius) as f32,
                    (c.offset.z / frame_radius) as f32,
                ],
                scale: (n.matter.radius / frame_radius) as f32,
                angular_size: if c.angle.is_finite() { c.angle as f32 } else { 0.0 },
                detail,
            });
        }
        scene
    }

    fn facts_of(&self, idx: NodeIdx) -> NodeFacts {
        let n = &self.tree.nodes[idx.get()];
        NodeFacts {
            key: n.key,
            epoch: n.epoch,
            tier: n.tier.index() as u8,
            solver: crate::solvers::for_tier(n.tier) as u8,
            depth: n.depth,
            mass: n.matter.mass,
            radius: n.matter.radius,
            temperature: n.matter.temperature,
            internal_energy: n.matter.internal_energy,
            binding_energy: n.matter.binding_energy,
            luminosity: n.matter.luminosity,
            charge: n.matter.charge,
            cadence: self.node_cadence(idx),
            timestep: self.node_dt(idx),
            mixing_time: self.mixing_time(idx),
            lag: self.render_lag(idx),
            body_count: n.bodies.len() as u32,
            structural_parts: n.last_report.structural_parts as u32,
            materialised: n.is_materialised(),
            pinned: n.pinned,
        }
    }

    /// Everything within the volume, largest on screen first, capped.
    ///
    /// Without a volume this is just the node asked for, which is the whole of
    /// the old behaviour and costs one push.
    fn gather(&self, root: NodeIdx, frame_radius: f64, vol: Option<&Volume>) -> Vec<Candidate> {
        let Some(vol) = vol else {
            return vec![Candidate { idx: root, offset: Vec3::ZERO, angle: f64::INFINITY }];
        };
        let eye = Vec3 { x: vol.eye[0], y: vol.eye[1], z: vol.eye[2] }.scale(frame_radius);
        let reach = vol.reach * frame_radius;

        let mut out = vec![Candidate { idx: root, offset: Vec3::ZERO, angle: f64::INFINITY }];
        // Descend through *promoted* children only. An unmaterialised node's
        // contents are a recipe, not a subtree — which is exactly why a city
        // of untouched buildings costs nothing to walk.
        let mut stack = vec![(root, Vec3::ZERO, Quat::IDENTITY)];
        while let Some((i, pos, rot)) = stack.pop() {
            let n = &self.tree.nodes[i.get()];
            for &c in &n.children {
                if c.is_none() || c.get() >= self.tree.nodes.len() {
                    continue;
                }
                let child = &self.tree.nodes[c.get()];
                if !child.alive {
                    continue;
                }
                let cpos = pos + rot.rotate(child.motion.offset);
                let crot = rot.then(child.motion.orientation);
                let d = (cpos - eye).norm();
                // Children are contained in their parent, so a node whose own
                // sphere misses the reach takes its whole subtree with it.
                if d - child.matter.radius > reach {
                    continue;
                }
                let angle = child.matter.radius / d.max(1e-300);
                out.push(Candidate { idx: c, offset: cpos, angle });
                // Below the detail angle there is nothing to draw and nothing
                // inside it can be bigger, so the walk stops. This is what
                // bounds the cost of a query by the screen rather than by the
                // world.
                if angle >= vol.detail_angle {
                    stack.push((c, cpos, crot));
                }
            }
        }

        // Largest first, so the cap costs the things furthest away. The frame
        // node keeps its place at the head: it is what was asked for.
        out[1..].sort_by(|a, b| b.angle.total_cmp(&a.angle));
        out.truncate(vol.max_nodes.max(1));
        out
    }

    /// What to say about one node.
    fn detail_of(
        &self,
        idx: NodeIdx,
        req: &ViewRequest,
        held: Option<&HashMap<PathKey, Held>>,
        want: usize,
        angle: f64,
    ) -> Detail {
        let n = &self.tree.nodes[idx.get()];

        // Too small on screen to draw as anything but a dot. The facts still
        // travel; the bodies do not.
        if let Some(v) = &req.volume {
            if angle < v.detail_angle {
                return Detail::Unchanged;
            }
        }

        // Does the client already have this, at this epoch, solved no earlier
        // than it was told? Then it can carry its own bodies forward.
        let fresh = match held {
            Some(h) => h
                .get(&n.key)
                .is_some_and(|h| h.epoch == n.epoch && n.last_solved <= h.at),
            None => n.last_solved <= req.since,
        };
        if fresh {
            return Detail::Unchanged;
        }

        let lag = self.render_lag(idx);
        if n.is_materialised() {
            let stride = stride_for(n.bodies.len(), want);
            return Detail::Explicit(specks_of(&n.bodies, n.matter.radius, lag, stride));
        }

        // Nothing materialised. The server has no bodies to send — but it has
        // the ~300 bytes that make them, and the client has the same sampler.
        if !req.allow_recipes || n.pinned || self.tree.persisted.contains_key(&n.key) {
            return Detail::Unchanged;
        }
        match Recipe::of(n, self.tree.world_seed, lag, want) {
            Some(r) => Detail::Recipe(r),
            None => Detail::Unchanged,
        }
    }
}

impl Recipe {
    /// The instructions that would produce this node's detail, if there are
    /// any. `None` when the node's detail is not derivable — which is checked
    /// by the caller, but is cheap to state twice for something this easy to
    /// get wrong.
    fn of(n: &crate::tree::Node, seed: u64, lag: f64, want: usize) -> Option<Recipe> {
        if n.pinned || n.spec.count == 0 {
            return None;
        }
        let mut w = Writer::new();
        crate::persist::put_matter(&mut w, &n.matter);
        crate::persist::put_spec(&mut w, &n.spec);
        match &n.morphology {
            Some(m) => {
                w.bool(true);
                crate::persist::put_morphology(&mut w, m);
            }
            None => w.bool(false),
        }
        let blob = w.finish();
        if blob.len() > MAX_RECIPE_BYTES {
            return None;
        }
        Some(Recipe {
            key: n.key,
            epoch: n.epoch,
            seed,
            count: n.spec.count as u32,
            stride: stride_for(n.spec.count, want) as u32,
            lag,
            checksum: fnv1a(&blob),
            blob,
        })
    }
}

/// Share one body budget across the nodes a query found.
///
/// Weighted by solid angle — the square of the half-angle — because that is
/// what a body is competing for: screen area. A node twice as wide gets four
/// times the specks, which is the right answer if the goal is even density on
/// the display rather than even density in the world.
///
/// A budget of zero means unlimited, and every node gets everything it has.
fn share_bodies(w: &crate::engine::World, found: &[Candidate], budget: usize) -> Vec<usize> {
    if budget == 0 {
        return vec![0; found.len()];
    }
    if found.len() == 1 {
        return vec![budget];
    }
    let weight = |c: &Candidate| {
        if c.angle.is_finite() {
            (c.angle * c.angle).max(0.0)
        } else {
            1.0
        }
    };
    // The frame node was asked for by name, so it is never squeezed out by a
    // crowd of small things around it: it takes a fixed half before the rest
    // compete for what is left.
    let head = (budget / 2).max(1);
    let total: f64 = found[1..].iter().map(weight).sum();
    let mut out = Vec::with_capacity(found.len());
    out.push(head.min(node_bodies(w, found[0].idx).max(1)));
    let rest = budget.saturating_sub(out[0]);
    for c in &found[1..] {
        let share = if total > 0.0 { rest as f64 * weight(c) / total } else { 0.0 };
        // A node worth drawing at all is worth more than one speck; below a
        // handful the shape is gone and the bytes were wasted either way.
        out.push((share as usize).max(8));
    }
    out
}

fn node_bodies(w: &crate::engine::World, idx: NodeIdx) -> usize {
    let n = &w.tree.nodes[idx.get()];
    if n.is_materialised() {
        n.bodies.len()
    } else {
        n.spec.count
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
