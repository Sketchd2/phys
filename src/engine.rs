//! The orchestrator: what actually happens in 50 milliseconds.
//!
//! Each frame the engine does five things, in this order:
//!
//! 1. **Survey.** Walk the live tree; for every node work out how far it is
//!    from each observer, what tier that observer needs there, and what it
//!    would cost to provide it. Emit a `Task` for every discrepancy.
//! 2. **Plan.** Hand the tasks to the frame budget, which fits what it can into
//!    the wall clock and reports the rest as detail debt.
//! 3. **Execute.** Materialise, coarsen, promote and step, in planned order.
//! 4. **Deliver.** Hand over influences whose light has finally arrived.
//! 5. **Record.** Push each active node's state into its history ring, so that
//!    future observers have a past light cone to read from.
//!
//! The ordering matters. Surveying before planning means the budget sees the
//! whole demand and can choose; executing before delivering means an influence
//! never lands on a node that has not yet been stepped to the right time.

use crate::budget::{cost, FrameBudget, Plan, Task, TaskKind};
use crate::causal::{CausalGate, Clock, History, Influence, InfluenceKind, Mailbox, Moment};
use crate::neighbourhood::{Occupant, Reservoir};
use crate::ids::{EntityId, NodeIdx, PathKey};
use crate::math::Vec3;
use crate::observe::*;
use crate::sampler::SampleSpec;
use crate::rng::{Purpose, Stream};
use crate::solvers::{self, SolverKind};
use crate::state::{Matter, Body};
use crate::tree::{Residency, Tree};
use crate::units::*;
use std::collections::HashMap;

/// How many individual joint failures a structure remembers by name. Beyond
/// this the damage is real but anonymous — it shows in the mass, not in the
/// record of which twig went.
/// Most sub-steps one node may take to reach the world instant in a frame.
///
/// Not a stability limit — the stability limit is [`World::node_dt`]. This is
/// the point at which following the trajectory stops being the cheaper way to
/// answer the question, and the node is crossed by its ensemble instead. Set
/// high enough that anything watchable is integrated properly, and low enough
/// that a nucleus resolved inside a galaxy cannot stall the frame.
pub const MAX_SUBSTEPS: u32 = 256;

/// Substeps a molecular node will take inside one `advance_node` call.
///
/// Separate from [`MAX_SUBSTEPS`], and smaller, because it is bounding a
/// different thing: that one bounds passes over a *span*, this one bounds work
/// inside a single pass. A force field can ask for a million substeps once it
/// has looked at the configuration, and this is the point past which the node
/// stops covering the whole step and covers the part it can integrate stably
/// instead.
pub const MD_MAX_SUBSTEPS: u32 = 64;

/// Refinement error above which a node resolves itself, with or without an
/// audience.
///
/// `refinement_error` is about 0.05 for a node the matter describes
/// perfectly and climbs past one when the matter has started lying: an
/// unresolved Jeans length, or a dynamical time shorter than the frame. Set
/// just above the quiet value, so "something is happening here" is what
/// triggers detail and idleness is what does not.
pub const REFINE_THRESHOLD: f64 = 1.0;

/// Slowest the world may run relative to the pace it was asked for.
///
/// A millionfold slowdown has said everything slowing down can say. Past it
/// the honest answer is not a slower clock but a staler world, and that is
/// what `stats.worst_lateness` is for.
pub const MIN_TIME_THROTTLE: f64 = 1e-6;

pub const NOTABLE_BREAKS: usize = 48;

/// Substeps a single `shake` call may take. A structure whose period is far
/// shorter than the frame it is asked to cover would otherwise spend the whole
/// frame budget resolving motion nobody can see; past this it is integrated
/// coarsely and the report says so through its own convergence flag.
pub const MAX_SHAKE_STEPS: f64 = 240.0;

/// How many structures may be integrated through time at once. Past this the
/// least recently started is dropped: it is a stand in front of an observer,
/// not a forest.
pub const MAX_SHAKEN: usize = 64;

/// Below this much standing mass a structure is rubble, not a structure, and
/// there is nothing meaningful left to analyse.
pub const COLLAPSE_MASS: f64 = 1e-6;

/// What one step of matter evolution did.
#[derive(Debug, Clone, Copy, Default)]
pub struct MatterReport {
    /// Nodes whose thermal account moved.
    pub nodes: usize,
    /// Energy radiated away, J. Leaves the world: the sky is not a node.
    pub radiated: f64,
    /// Energy absorbed from incident light, J.
    pub absorbed: f64,
    /// Energy a node was asked to radiate and did not have, J.
    ///
    /// Non-zero means something is being made to shine out of reserves it does
    /// not hold — usually a node whose luminosity was authored rather than
    /// derived. Zero for a world nobody has adjusted, which is the assertion
    /// worth making in a test.
    pub radiation_deficit: f64,
}

/// What one step of falling debris did.
#[derive(Debug, Clone, Copy, Default)]
pub struct FallReport {
    /// Contacts detected this step.
    pub contacts: usize,
    /// Members of standing structures that took an impulse.
    pub struck_members: usize,
    /// Joints those impulses broke.
    pub secondary_breaks: usize,
    /// Mass those breaks brought down, kg.
    pub secondary_mass: f64,
    /// Members that failed inside a piece while it was falling.
    pub broken_while_falling: usize,
    /// Pieces that came to rest this step.
    pub settled: usize,
    /// Pieces still in the air.
    pub still_falling: usize,
    /// Highest utilisation any struck structure reached, so an impact that did
    /// nothing can be told from one that never arrived.
    pub peak_utilisation: f64,
    /// Largest single impulse delivered, N s.
    pub largest_impulse: f64,
}

/// Where the ground is, in a structure's own frame.
///
/// Not zero. A generated structure is recentred on its own centre of mass, so
/// its foundations sit at whatever negative height that put them; assuming zero
/// makes debris fall through the floor or land in mid-air depending on the
/// structure. What anchors a structure is what it is standing on.
fn ground_of(topo: &crate::topology::Topology) -> f64 {
    let mut lowest = f64::INFINITY;
    for i in 0..topo.support.len() {
        if topo.support[i] == crate::morph::NO_SUPPORT && topo.joints[i].radius > 0.0 {
            lowest = lowest.min(topo.base[i].z.min(topo.tip[i].z));
        }
    }
    if lowest.is_finite() {
        lowest
    } else {
        topo.base.iter().map(|p| p.z).fold(f64::INFINITY, f64::min)
    }
}

/// How long a piece may fall before it is written off as litter. A limb that
/// has been in the air for this long has either landed somewhere the collision
/// search cannot see or is falling forever, and neither is worth a frame.
pub const MAX_FALL_SECONDS: f64 = 12.0;

/// How many pieces may be falling at once. A crown fire breaks thousands of
/// joints; simulating every twig's descent would spend the whole budget on
/// debris nobody is looking at, and the mass is already accounted for whether
/// its fall is drawn or not.
pub const MAX_FALLING: usize = 24;

/// What a run of dynamics did.
#[derive(Debug, Clone, Default)]
pub struct ShakeOutcome {
    /// Substeps taken.
    pub steps: u32,
    /// Conjugate-gradient iterations across all of them.
    pub iterations: u32,
    /// Joints that failed while the structure was moving.
    pub broken_joints: usize,
    /// Pieces that came away as their own falling objects.
    pub detached_pieces: usize,
    /// Structural mass that fell off.
    pub detached_mass: f64,
    /// Kinetic energy the structure is carrying, J.
    pub kinetic: f64,
    /// Elastic energy stored in it, J.
    pub strain: f64,
    /// Strain energy released by members that failed, J.
    pub released: f64,
    /// Energy removed by damping, J.
    pub dissipated: f64,
    /// Largest nodal displacement, m.
    pub displacement: f64,
    /// Largest chord rotation of any member, radians. Above about 0.1 the
    /// small-displacement assumption is spent and the restoring force is being
    /// overestimated.
    pub displacement_ratio: f64,
    /// A step failed to converge and the run stopped early.
    pub diverged: bool,
}

/// What an insult did to a structure.
#[derive(Debug, Default, Clone, Copy)]
pub struct DamageOutcome {
    /// Joints that failed.
    pub broken_joints: usize,
    /// Structural mass that fell off and is now litter in the same node.
    pub detached_mass: f64,
    /// Structural mass destroyed outright — burned or vaporised. The atoms stay
    /// in the node as combustion products.
    pub consumed_mass: f64,
    /// Free energy liberated by that destruction, J.
    pub energy_released: f64,
    /// Energy the insult itself delivered, J.
    pub energy_delivered: f64,
    /// Highest joint utilisation reached, whether or not anything broke. Below
    /// 1 the structure rode it out.
    pub peak_utilisation: f64,
    /// Whether the structure was statically indeterminate and needed a solve.
    pub indeterminate: bool,
    /// Conjugate-gradient iterations used; zero on the exact path.
    pub solver_iterations: u32,
    /// The structure no longer exists.
    pub collapsed: bool,
    /// Pieces that came away as their own falling objects.
    pub detached_pieces: usize,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct EngineStats {
    pub frames: u64,
    pub sim_time: f64,
    pub tasks_run: u64,
    pub tasks_deferred: u64,
    pub bodies_stepped: u64,
    pub detail_debt: f64,
    pub worst_causality_violation: f64,
    pub last_frame_us: f64,
    pub materialised_bodies: usize,
    pub live_nodes: usize,
    /// Worst lateness anywhere in the world, in units of that node's own
    /// characteristic time. Under one, every node was re-solved before it had
    /// changed appreciably. This is the engine's honest "is the world keeping
    /// up" number, and unlike a frame time it is scale-free.
    pub worst_lateness: f64,
    /// Nodes that were due this frame and did not fit in the budget.
    pub overdue: usize,
    /// Nodes carried across the frame by their ensemble rather than their
    /// trajectory, because the span was too long to integrate.
    pub thermalised: u64,
    /// Nodes whose own physics needs a smaller step than the frame can afford,
    /// and which cannot be thermalised out of the difficulty because something
    /// is watching them — pinned, bubbled, or built upon.
    ///
    /// Unlike [`Self::overdue`], which is a node that lost a race for this
    /// frame's budget and will win a later one, this is a standing condition:
    /// the shortfall recurs every frame and the node's lateness grows without
    /// bound. Non-zero here means the world is being asked to run at a pace
    /// that something in it cannot be integrated at, and the answer is to slow
    /// the pace, coarsen the node, or stop watching it — not to wait.
    pub unreachable: u64,
    /// Nodes carried forward in closed form without being re-solved. The
    /// overwhelming majority, every frame, and the reason one instant is
    /// affordable across thirty-eight orders of magnitude.
    pub coasted: usize,
    /// Seconds of interior evolution handed out beyond what the world clock
    /// paid for, summed over every bubbled node.
    ///
    /// The honest headline for a world somebody has been adjusting. Energy is
    /// still conserved inside a bubble — a tree aged a century really does burn
    /// a century of its own reserves — but what it *received* across its
    /// boundary was one year of sunlight, and this is the size of that gap in
    /// the only unit that covers every process at once. Zero for a world nobody
    /// has touched, which is the assertion worth making in a test.
    pub bubble_seconds: f64,
    /// Nodes currently carrying a bubble factor other than one.
    pub bubbled: usize,
    /// Nodes that are made of something, and so run chemistry each frame.
    pub reacting_nodes: usize,
    /// Energy radiated out of the world, J, summed over nodes and frames. The
    /// sky is not a node, so this genuinely leaves rather than moving.
    pub radiated: f64,
    /// Energy absorbed from incident light, J, likewise summed.
    pub absorbed: f64,
    /// Energy nodes were asked to radiate and did not hold. Non-zero means
    /// something is shining out of reserves it does not have.
    pub radiation_deficit: f64,
    /// Mass fractions moved into and out of solution, summed over nodes and
    /// frames. Not a mass — a node's fraction — so it is a measure of how much
    /// chemistry is happening rather than of how much matter there is.
    pub dissolved: f64,
    pub precipitated: f64,
    /// Melting, freezing, boiling and condensing, likewise.
    pub phase_changed: f64,
    /// The worst ratio of a node's actual extent to the radius it claims, over
    /// every node the frame advanced. One means the contents exactly fill what
    /// the node says it is; above one they have outgrown it, and every length
    /// derived from that radius — SPH's smoothing length, the gravity
    /// softening, the LOD's angular size, the neighbour grid's spacing — is
    /// wrong by this factor. See `state::Spread`.
    pub worst_occupancy: f64,
    /// The node it was measured on, so a number worth chasing says where to
    /// look. `PathKey` and not `EntityId`: this is a measurement of a place.
    /// `None` until a frame has advanced something, which is not the same as
    /// zero and should not be spelled like it.
    pub worst_occupancy_at: Option<PathKey>,
    /// Boundary crossings, summed over frames — `docs/PLAY.md` D16. A node
    /// that left the region its parent owns and was re-homed to whatever owns
    /// it now.
    pub crossings: u64,
    /// Of those, the ones that landed *sideways*, in something the arbiter
    /// already held, rather than outward into the arbiter itself. The creature
    /// walking from the forest into the desert, which D16 says must never
    /// become a direct child of the planet on its way.
    pub crossings_sideways: u64,
    /// Of those, the ones whose destination had to be **generated** to be
    /// arrived at: the inward row of D16's table, where the place a thing
    /// crossed into was still only one of its parent's bodies and is promoted
    /// to meet it.
    pub crossings_generated: u64,
    /// Nodes that stopped being one neighbourhood and became two, summed over
    /// frames. `docs/BACKLOG.md`'s "a node cannot split when its contents
    /// spread out", which Phase 1 measured and connected to nothing.
    pub splits: u64,
    /// Sibling nodes that became one neighbourhood again and were folded back
    /// into one.
    pub merges: u64,
    /// How far outside its parent's contents the worst crosser was, as a
    /// multiple of what those contents reach.
    ///
    /// **The crossing pass is the detector now.** A node that steps over a
    /// boundary crosses at a ratio a hair above one; a node *flung* over it
    /// arrives with a number that says so. `worst_occupancy` used to report
    /// this class of fault — `docs/BACKLOG.md` measures the ladder at 10^6 and
    /// a nucleus at 10^22 — and it did so only because nothing re-homed the
    /// victim, so it sat inside a node claiming a metre and the ratio kept
    /// climbing. Re-homing it fixes the tree and would have made the fault
    /// invisible, which is why the measurement moved here rather than going
    /// away.
    pub worst_crossing: f64,
    /// The node it was measured on. `PathKey` and not `EntityId`: a
    /// measurement of a place, and the place it names is the one that flung it.
    pub worst_crossing_at: Option<PathKey>,
    /// Nodes measured as having left the region their parent owns, with nowhere
    /// to go: a child of the root, which has no outside. Not an error and not
    /// silent — it is how the sampler faults `docs/BACKLOG.md` records at
    /// 10^5 and 10^11 radii announce themselves once the measurement runs every
    /// frame.
    pub crossings_refused: u64,
    /// Nodes crossed by their **ensemble** rather than followed, summed over
    /// frames. `docs/PLAY.md` §3.7: the resolution floor is real and derivable,
    /// and the engine should *report* reaching it rather than silently dropping
    /// a node to its equilibrium — the discipline `displacement_ratio` already
    /// applies to the small-displacement regime.
    ///
    /// Distinct from [`Self::unreachable`], which counts the nodes that could
    /// *not* be dropped: pinned, bubbled, or with something built on them, so
    /// somebody is deliberately watching them run. Those fall behind instead.
    /// Between them the two account for every node the floor caught.
    pub ensembled: u64,
    /// The most substeps any node asked for, **uncapped**, over frames.
    ///
    /// The number to compare against `MAX_SUBSTEPS`: at or below it the node is
    /// followed, above it the trajectory is not merely expensive but the wrong
    /// answer. Reported rather than clamped, because a node wanting 289 million
    /// substeps and a node wanting 257 are both "capped" and are not the same
    /// situation.
    pub worst_substeps: f64,
    /// The node that asked. A `PathKey`, because this is a measurement of a
    /// place; `None` until some node has been surveyed.
    pub worst_substeps_at: Option<PathKey>,
    /// Overlaps resolved into an impulse pair, summed over nodes and frames.
    /// Counts contacts, not newton-seconds, for the same reason
    /// `exchange_crossings` counts crossings: a coupling that silently stops
    /// happening looks identical to one that has nothing to do.
    pub contacts_resolved: u64,
    /// Surfaces baked. `PLAY.md` D13 stores a boundary and regenerates it when
    /// `epoch` moves, so this counts arrangements that changed, not frames.
    pub surfaces_baked: u64,
    /// Cells promoted because an observer came close enough to want them. See
    /// [`World::approach`].
    pub approaches: u64,
    /// Surfaces whose materials did not reconcile with the node's own solid
    /// pools. `PLAY.md` D18's invariant; non-zero means a node is presenting
    /// something its bulk contradicts.
    pub surface_mismatches: u64,
    /// The last such disagreement, so a number worth chasing says where to look.
    pub worst_surface_mismatch: Option<crate::shape::Mismatch>,
    /// Node-steps where the fluid solver priced matter the gas law does not
    /// describe. `PLAY.md` §7's ninth Phase 2 item; **Water** is what fixes it.
    pub eos_outside_validity: u64,
    /// Where, so a number worth chasing says where to look.
    pub eos_outside_validity_at: Option<crate::ids::PathKey>,
    /// Boundaries heat actually crossed, summed over nodes and frames. Counts
    /// transfers, not joules: it answers "is anything talking to its
    /// neighbours at all", which is the question a coupling that silently does
    /// nothing would otherwise pass every test on.
    pub exchange_crossings: u64,
}

/// Where the world clock's span per frame comes from.
///
/// `refresh_pace` runs at the top of every frame, which is what makes "zooming
/// in slows time" an arithmetic consequence rather than a policy — and which
/// also meant an assignment to `pace` was silently discarded on the next frame
/// unless the caller knew to blank `paced_to` first. That incantation worked by
/// accident of an early return, read as a bug wherever it appeared, and nothing
/// stopped a refactor removing the return it depended on. This is the choice
/// made explicit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PaceMode {
    /// Follow `paced_to`'s own timescale, coupling resolution and time rate.
    ///
    /// Was the default. `docs/PLAY.md` D1 retired it as one: a shared world runs
    /// at one second per second, and a clock dragged slower by whoever is
    /// looking most closely is a single-observer answer to a shared question.
    /// Reached now only through an explicit [`World::pace_to`], which is what
    /// it is for — and [`World::time_throttle`] belongs to this mode with it.
    Follow,
    /// A fixed span per frame, set by the caller and left alone. What a world
    /// is, at `1.0`.
    #[default]
    Fixed,
}

/// Cells along a side of each face of a planet's ocean.
///
/// Sixteen puts a cell at 630 km on an Earth, which resolves a tide — the
/// forcing is the planet's own size — and the shallow-water waves that carry
/// it, which are a few thousand kilometres long. Not a physical number and not
/// pretending to be one: it is the resolution, as `SampleSpec::count` is.
pub const OCEAN_CELLS: usize = 16;

/// How fast a liquid's contents move, or will once they are let go: the
/// greater of `sqrt(2 g H)` for a column of depth `H` in a field `g`, the
/// liquid's bulk motion, and the fastest thing driving it — `drivers`, the
/// walls that move with a body through it. Weakly-compressible SPH runs a
/// liquid at ten times this, and the scheduler and the solver have to agree on
/// it — see `World::signal_speed_of`.
///
/// **Bulk motion is the rms speed, not the fastest parcel.** The sound speed
/// is the stiffness, so taking it from the fastest parcel is a feedback: one
/// parcel squirting from between a sinking ball and a floor stiffened the whole
/// bucket, which threw it faster, which stiffened it further — measured, from
/// 4.8 m/s to 147 m/s in two frames. WCSPH's own rule is that the sound speed
/// comes from the scene's velocity scale, not from its noise.
pub fn flow_scale<'a>(
    bodies: impl Iterator<Item = &'a crate::state::Body>,
    field: crate::math::Vec3,
    spacing: f64,
    drivers: f64,
) -> f64 {
    let g = field.norm();
    let (mut lo, mut hi, mut sum, mut mass) = (f64::INFINITY, f64::NEG_INFINITY, 0.0, 0.0);
    for b in bodies {
        let h = if g > 0.0 { -b.pos.dot(field) / g } else { 0.0 };
        lo = lo.min(h);
        hi = hi.max(h);
        sum += b.mass * b.vel.norm2();
        mass += b.mass;
    }
    let depth = if hi >= lo { (hi - lo).max(spacing) } else { spacing };
    let bulk = if mass > 0.0 { (sum / mass).sqrt() } else { 0.0 };
    bulk.max(drivers).max((2.0 * g * depth).sqrt())
}

/// The world.
pub struct World {
    pub tree: Tree,
    pub mailbox: Mailbox,
    pub ledger: Ledger,
    pub observers: Vec<Observer>,
    pub budget: FrameBudget,
    pub gate: CausalGate,
    /// The world instant, seconds since the scenario epoch.
    ///
    /// There is one, and everything in the world is at it. What differs
    /// between a galaxy arm and a nucleus is not what time it is for them, it
    /// is how often each is re-solved and how it was carried here — the arm in
    /// one closed-form step, the nucleus through ten thousand.
    pub time: f64,
    /// Simulated seconds the world advances per frame at a time rate of one.
    ///
    /// **One, unless something asked otherwise.** A world runs at one second per
    /// second — `docs/PLAY.md` D1 — and that is what a `World` is built with.
    ///
    /// It can still be taken from whatever is being watched, through
    /// [`World::pace_to`], and a frame then covers about one characteristic time
    /// of that node: a galaxy advances millennia, a carbon atom femtoseconds.
    /// That is a fine answer to a single observer's question and the wrong one
    /// to a shared question, so it survives as a tool for single-player
    /// exploration and offline study rather than as what a world does.
    pub pace: f64,
    /// Multiplier on `pace`. The user's time control.
    pub time_rate: f64,
    /// Fraction of the asked-for pace the engine is actually sustaining.
    ///
    /// One when everything due fits in the frame. Less when it does not — and
    /// the response to not fitting is to advance *less simulated time*, not to
    /// do less physics. That is the difference between a world that slows down
    /// under load and a world that stops: at the pace of a galaxy nothing's
    /// trajectory can be integrated at all, and an engine without this number
    /// answers by deferring every solve and standing still.
    ///
    /// **It belongs to [`PaceMode::Follow`] and applies only there.**
    /// `docs/PLAY.md` D1 says overload must show up as a *staler* world rather
    /// than a slower one, and this is the second mechanism that slows one —
    /// the clock would still drag under load however the pace was set. So it
    /// goes where its own justification goes: the paragraph above is an
    /// argument about a galactic pace, and a galactic pace is now something a
    /// single observer asks for. A world at one second per second gets staler
    /// under load instead, in `worst_lateness` and detail debt.
    ///
    /// D1 also notes what that leaves open, and it is not built: the knapsack
    /// wants a lateness ceiling for anything an actor is interacting with, with
    /// resolution surrendered to hold it. Until then a fixed-pace world under
    /// sustained overload gets arbitrarily stale rather than arbitrarily slow,
    /// which is the trade D1 chooses.
    pub time_throttle: f64,
    /// The node the pace is taken from, re-read at the start of every frame.
    ///
    /// Held as a node rather than a number because a node's cadence changes
    /// under it: materialising a galaxy into twenty thousand stars shortens its
    /// characteristic time by two orders of magnitude, and a pace fixed when it
    /// was still unmaterialised would then be asking for a span its own
    /// stars could not be integrated across.
    pub paced_to: NodeIdx,
    /// Retained history, only for nodes something might observe from a
    /// distance. Keeping it keyed by path rather than in the node itself means
    /// a node can be coarsened and rebuilt without losing its past.
    /// Retained history, only for nodes something might observe from a
    /// distance.
    ///
    /// Keyed by *address*, not by name, and deliberately. A history is
    /// bookkeeping the scheduler creates — which nodes get one is decided by
    /// the frame budget from a wall-clock allowance — so naming a node for it
    /// would make identity depend on how fast the machine is, and
    /// `next_entity` is persisted. Measured: a slower machine saved a
    /// different world. See `docs/BACKLOG.md`.
    pub histories: HashMap<PathKey, History>,
    /// The same, for the same reason: a clock is something the budget made,
    /// not something that happened.
    pub clocks: HashMap<PathKey, Clock>,
    /// Address to identity. The one index a move has to migrate.
    ///
    /// Every other side table is keyed by [`EntityId`], so `reparent` leaves
    /// them alone. This one cannot be: a node discarded and rebuilt has to
    /// recover its name, and its address is the only thing it comes back with.
    /// See [`EntityId`] for why that is one line rather than one per table.
    ///
    /// Entries are issued lazily — a node that nothing has ever recorded
    /// against has no identity and needs none.
    pub identities: HashMap<PathKey, EntityId>,
    /// Next identity to issue. Monotonic, never reused, persisted.
    pub next_entity: u64,
    pub stats: EngineStats,
    /// Whether the clock follows a node or is driven by hand.
    pub pace_mode: PaceMode,
    /// Actions that deliberately broke conservation, with what they cost.
    pub audit: Vec<AuthorEvent>,
    /// Whether cost estimates assume the GPU path.
    pub gpu: bool,
    /// Construction rate available to planned programs, fraction of a design
    /// per second. Zero means no crews are working.
    pub labour_rate: f64,
    /// Growth steps refused because their transaction did not balance.
    pub rejected_growth_steps: u64,
    /// Per-node environment overrides, keyed by path so they survive the node
    /// being coarsened and rebuilt.
    pub environments: HashMap<EntityId, crate::morph::Environment>,
    /// Every substance this world has ever analysed.
    ///
    /// World state, not scenery: a node's mixture names substances by id, so a
    /// world reloaded without its catalogue would be pointing at nothing.
    pub substances: crate::chem::Registry,
    /// The structures currently being integrated through time.
    ///
    /// A bounded set, deliberately. Dynamics is expensive and it is only worth
    /// anything to somebody watching — a tree swaying in a forest nobody is
    /// looking at is indistinguishable from one standing still, and the engine
    /// already declines to materialise what nobody can see. This is the same
    /// rule applied to motion rather than to detail: the stand in front of the
    /// observer moves, and the forest behind them does not.
    shaking: Vec<(NodeIdx, crate::solvers::structure::DynamicStructure)>,
    /// Pieces that have come away and are still falling, with the node they
    /// fell from.
    falling: Vec<(NodeIdx, crate::solvers::structure::Fragment)>,
    history_depth: usize,
}

impl World {
    pub fn new(tree: Tree, ups: f64) -> World {
        let mut w = World {
            tree,
            mailbox: Mailbox::new(),
            ledger: Ledger::new(),
            observers: Vec::new(),
            budget: FrameBudget::ups(ups),
            gate: CausalGate::new(1e3 * YEAR),
            time: 0.0,
            pace: 1.0,
            time_rate: 1.0,
            time_throttle: 1.0,
            paced_to: NodeIdx::NONE,
            histories: HashMap::new(),
            clocks: HashMap::new(),
            identities: HashMap::new(),
            next_entity: 1,
            stats: EngineStats::default(),
            pace_mode: PaceMode::Fixed,
            audit: Vec::new(),
            gpu: false,
            labour_rate: 0.0,
            rejected_growth_steps: 0,
            environments: HashMap::new(),
            substances: crate::chem::Registry::new(),
            shaking: Vec::new(),
            falling: Vec::new(),
            history_depth: 64,
        };
        // A world runs at one second per second. `docs/PLAY.md` D1: that is
        // what a shared world *is*, and the clock must not be dragged slower by
        // whoever is looking most closely — a player inspecting a rifle bolt
        // must not slow down the war.
        //
        // This used to call `pace_to(root)`, so a world was watchable the moment
        // it was built without anybody having to know what timescale it lived
        // on. That is still worth having and is still one call away; it is a
        // single-player convenience rather than what a world is. `paced_to` is
        // pointed at the root regardless, so `pace_to` has a subject to return
        // to and a viewer has something sensible to offer.
        w.paced_to = w.tree.root;
        w.pace_realtime();
        w
    }

    /// Everything durable about this world, ready to hand to a store.
    ///
    /// Cheap in the sense that matters: it moves no fine detail that could be
    /// regenerated, because `persist` declines to write it.
    ///
    /// **`&mut self`, and that is the point.** A save writes each node's matter
    /// and throws away the bodies of every unpinned one, while a materialised
    /// node's authority is its *bodies* — so a world checkpointed mid-solve
    /// used to lose whatever the solver had done since the node was
    /// materialised. Measured at 2.35x10^-8 of the root's energy on the
    /// reference world, and zero on one at rest, which is why it went unnoticed
    /// for so long.
    ///
    /// The fix is [`Tree::settle`], and taking the world by `&mut` is what
    /// makes it impossible to route around: there is no way to obtain a
    /// `WorldView` of a world whose matter has not been brought into step with
    /// its own detail. That mattered more than the borrow was worth, because
    /// `docs/PLAY.md` D16 turns the mid-flight checkpoint from a rarity into
    /// the ordinary case — a crossing is a save point in all but name.
    pub fn view(&mut self) -> crate::persist::WorldView<'_> {
        self.settle_all();
        crate::persist::WorldView {
            tree: &self.tree,
            ledger: &self.ledger,
            time: self.time,
            pace: self.pace,
            time_rate: self.time_rate,
            time_throttle: self.time_throttle,
            paced_to: self.paced_to,
            pace_mode: self.pace_mode,
            labour_rate: self.labour_rate,
            rejected_growth_steps: self.rejected_growth_steps,
            environments: &self.environments,
            identities: &self.identities,
            next_entity: self.next_entity,
            substances: &self.substances,
            audit: &self.audit,
            mailbox: &self.mailbox,
        }
    }

    /// Bring every materialised node's matter into step with its own detail.
    ///
    /// The pre-pass [`World::view`] exists for. Deepest node first, because a
    /// parent's stand-in body is written from its promoted child's matter and
    /// that matter has to be current before it is read — settling a parent
    /// before its child would summarise last frame's child into this frame's
    /// parent.
    ///
    /// Costs one `summarise` per materialised node per save and nothing at all
    /// per frame, which is the trade `docs/BACKLOG.md` set out: the alternative
    /// was keeping every solved node's matter in step as it went, which is
    /// precisely the work the materialised-detail design exists to avoid.
    ///
    /// A node whose detail says nothing new keeps its matter untouched — see
    /// the idempotence rule in [`Tree::settle`] — so a world nobody has
    /// disturbed still saves bit-for-bit as the world it was.
    pub fn settle_all(&mut self) {
        let mut order: Vec<(u32, NodeIdx)> = self
            .tree
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.alive && n.is_materialised())
            .map(|(i, n)| (n.depth, NodeIdx(i as u32)))
            .collect();
        order.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.0.cmp(&b.1.0)));
        for (_, idx) in order {
            self.tree.settle(idx);
        }
    }

    /// Rebuild a world from a snapshot.
    ///
    /// `ups` comes from the caller rather than from the file: the frame budget
    /// describes the machine the world is being run on, not the world. A save
    /// made on a workstation opens on a laptop at the laptop's budget.
    pub fn from_snapshot(s: crate::persist::Snapshot, ups: f64) -> World {
        let mut w = World::new(s.tree, ups);
        w.ledger = s.ledger;
        w.time = s.time;
        w.pace = s.pace;
        w.time_rate = s.time_rate;
        w.time_throttle = s.time_throttle;
        w.paced_to = s.paced_to;
        w.pace_mode = s.pace_mode;
        w.labour_rate = s.labour_rate;
        w.rejected_growth_steps = s.rejected_growth_steps;
        w.environments = s.environments;
        w.identities = s.identities;
        w.next_entity = s.next_entity.max(1);
        w.substances = s.substances;
        w.audit = s.audit;
        w.mailbox = crate::causal::Mailbox::restore(s.in_flight, s.delivered, s.in_flight_peak);
        w.stats.sim_time = s.time;
        w
    }

    /// Place a structure and give it conditions to grow in.
    /// Seed a structure on a node.
    ///
    /// `env` is `None` for the ordinary case: the environment is then *derived*
    /// from physics every frame — light from what is shining on the node, water
    /// from the liquid phase of what it is made of, temperature from its own
    /// matter — so a tree planted in a place that later freezes stops growing
    /// without anyone arranging it.
    ///
    /// `Some(env)` **pins** an authored environment for that node forever, and
    /// is for scenarios that are placing a situation rather than simulating one
    /// — a lit planetary surface with no star in the tree to light it. It was
    /// once the only option, and it quietly defeated the derivation: every
    /// planted node had an override and none of them ever felt the weather.
    pub fn plant(
        &mut self,
        idx: NodeIdx,
        program: crate::morph::Program,
        env: Option<crate::morph::Environment>,
    ) {
        let id = self.identify(idx);
        self.tree.plant(idx, program);
        match env {
            Some(e) => {
                self.environments.insert(id, e);
            }
            None => {
                self.environments.remove(&id);
            }
        }
        self.rewrite_recipe(idx);
    }

    /// Ask of every node that might have become one whether it is a body now.
    ///
    /// **Once a frame, and cheap because of the order the questions are in.** A
    /// node that already has a recipe is skipped on a pointer test, and one
    /// nobody has said what it is made of is skipped on an emptiness test — so
    /// the cost on an ordinary world is a handful of comparisons. What is left
    /// is a node made of something, which is asked whether its matter is packed
    /// and then, only if it is, whether its own weight has beaten its strength.
    ///
    /// This is what turns a cloud into a planet while the world runs, which is
    /// the whole of "generated from origin": nothing states that a planet is
    /// there, and the first frame in which one is, is the frame its own
    /// collapse made it one.
    fn derive_layouts(&mut self) {
        for i in 0..self.tree.nodes.len() {
            let n = &self.tree.nodes[i];
            if !n.alive || n.morphology.is_some() || n.matter.mixture.is_empty() {
                continue;
            }
            if self.assess_surface(NodeIdx(i as u32)) {
                self.tree.stats.layouts_derived += 1;
            }
        }
    }

    /// Follow every observer down the surface they are over.
    ///
    /// **"Nothing is generated until approached."** `docs/PLAY.md` D6: "A patch
    /// refines into sub-patches as an observer descends and coarsens behind
    /// them, and the terrain regenerates bit-identically. A planet nobody has
    /// visited costs its `Matter` and nothing else."
    ///
    /// Materialising is the scheduler's business and already follows angular
    /// size; what the scheduler has never done is *promote*, and a patch's
    /// cells are bodies until one of them is promoted into a node of its own.
    /// So this walks from each observer's anchor down towards the observer,
    /// promoting the cell it is over for as long as the patch is coarser than
    /// the observer can resolve. Twenty-four levels take an Earth from six
    /// thousand kilometres to a square metre, and each level is one refine and
    /// one promote.
    ///
    /// It issues no identity, which is the rule: which nodes a frame promotes
    /// depends on a wall-clock allowance, so anything a promotion named would
    /// make identity depend on how fast the machine is.
    ///
    /// Coarsening behind the observer needs nothing here. A patch is a
    /// structure with a recipe, so `World::collapsible` already says its detail
    /// may be released, and the scheduler releases it as soon as the observer's
    /// acuity on it drops.
    fn approach(&mut self) {
        if self.observers.is_empty() {
            return;
        }
        for i in 0..self.observers.len() {
            let obs = self.observers[i];
            let mut here = obs.anchor;
            // Where the observer is, as a direction from the centre of the
            // body it is over. A patch's address is in the planet's own axes,
            // so that is the frame the question has to be asked in.
            let from_centre = {
                let mut ball = obs.anchor;
                let mut found = NodeIdx::NONE;
                while !ball.is_none() {
                    let is_ball = matches!(
                        self.tree.nodes[ball.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()),
                        Some(crate::recipe::Recipe::Tiled(t)) if t.is_ball()
                    );
                    if is_ball {
                        found = ball;
                        break;
                    }
                    ball = self.tree.nodes[ball.get()].parent;
                }
                if found.is_none() {
                    continue;
                }
                self.tree.separation(found, Vec3::ZERO, obs.anchor, obs.offset).value
            };
            // A bound on the walk rather than the thing that ends it: the walk
            // ends because the patch is fine enough or because the observer is
            // not over one of its cells.
            for _ in 0..64 {
                if here.is_none() || !self.tree.nodes[here.get()].alive {
                    break;
                }
                let is_tiled = matches!(
                    self.tree.nodes[here.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()),
                    Some(crate::recipe::Recipe::Tiled(_))
                );
                if !is_tiled {
                    break;
                }
                // How far the observer is from this patch's *surface*, and how
                // finely it wants to see it.
                //
                // Not the distance to the node's origin, which for a patch is
                // its centre of mass and for a whole face of a planet is 1.9
                // million metres underground. A descent measured against that
                // thinks an observer standing on the ground is halfway to the
                // core and stops eleven levels early.
                let Some(surface) = self.tree.nodes[here.get()].morphology.as_ref().and_then(|m| {
                    let Some(crate::recipe::Recipe::Tiled(t)) = m.recipe.as_ref() else {
                        return None;
                    };
                    Some((t.clamped_direction(from_centre).scale(t.sphere), t.side))
                }) else {
                    break;
                };
                let d = (from_centre - surface.0).norm().max(1e-30);
                let wanted = obs.linear_resolution(d);
                if !(surface.1 > wanted) {
                    break;
                }
                // The cell the observer is over, from the parameterisation
                // rather than from a search. See `Tiled::cell_of_direction`.
                self.tree.refine(here);
                let Some(best) = self.tree.nodes[here.get()].morphology.as_ref().and_then(|m| {
                    let Some(crate::recipe::Recipe::Tiled(t)) = m.recipe.as_ref() else {
                        return None;
                    };
                    t.cell_of_direction(from_centre)
                }) else {
                    break;
                };
                if best >= self.tree.nodes[here.get()].bodies.len() {
                    break;
                }
                let spec = self.tree.nodes[here.get()].spec;
                let child = self.tree.promote(here, best, spec);
                if child.is_none() {
                    break;
                }
                self.tree.nodes[child.get()].residency = Residency::Observed;
                self.stats.approaches += 1;
                here = child;
            }
        }
    }

    /// Fraction of a volume a poured pile of grains fills before it begins to
    /// carry load rather than flow.
    ///
    /// Random loose packing. A universal of sphere packing rather than a
    /// property of any material — the same kind of number as Turnbull's 0.45 —
    /// and what [`World::assess_surface`] uses to tell a planet from a cloud
    /// that is still falling in.
    const RANDOM_LOOSE_PACKING: f64 = crate::sampler::RANDOM_LOOSE_PACKING;

    /// How fast this node was cooling when it last passed a temperature, K/s.
    ///
    /// **Read off its own recorded past.** `History` keeps a temperature per
    /// node per frame for as long as the causal window holds, which is exactly
    /// the historical datum a formation needs: a melt that froze did so at
    /// whatever rate it happened to be losing heat at, and that rate is not
    /// recoverable from the node's state afterwards — a cold rock looks like a
    /// cold rock however it got there.
    ///
    /// `None` when the node has no recorded past that crosses the temperature
    /// downwards, which is the honest answer for something that was simply
    /// always cold.
    pub fn cooling_rate_through(&self, idx: NodeIdx, through: f64) -> Option<f64> {
        let key = self.tree.nodes.get(idx.get())?.key;
        let history = self.histories.get(&key)?;
        let mut previous: Option<crate::causal::Moment> = None;
        let mut found = None;
        for m in history.moments() {
            if let Some(p) = previous {
                if p.temperature > through && m.temperature <= through {
                    let dt = m.t - p.t;
                    if dt > 0.0 {
                        found = Some((p.temperature - m.temperature) / dt);
                    }
                }
            }
            previous = Some(*m);
        }
        found
    }

    /// What a node's own state and its own past say about how it is laid out.
    ///
    /// **The program is derived, not chosen.** `docs/PLAY.md` D11 asks for a
    /// genome rather than a species table, and the owner's call on this phase
    /// went further: a program is "a function that takes in the current and
    /// historical data, and decides on how it is laid out", derived from the
    /// axioms rather than engineered. This is the freezing half of that
    /// function; [`World::assess_surface`] is the rounding half.
    ///
    /// **Freezing is an event, and the event is the gate.** Nothing here asks
    /// whether a node *looks* like it froze, because it cannot be told from
    /// looking: a cold rock is a cold rock however it got there, and a gate on
    /// the state alone would hand a recipe to every solid node in the world.
    /// So this is called from `react_all`, on the node that just moved mass
    /// from liquid to solid — the one moment at which the fact is available.
    /// `react_all` runs for every node carrying a mixture whatever its
    /// residency, which is what keeps the answer from depending on who was
    /// watching.
    ///
    /// **The rate is measured, and drawn from equilibrium when there is nothing
    /// to measure.** A node with a recorded past has its real cooling rate
    /// through its own melting point, and that is used. A node with none gets
    /// the rate its own radiative balance implies at that temperature, which is
    /// D14's rule for matter nobody made — "the sampler draws formation
    /// conditions from the equilibrium the node is in" — and is the same
    /// arithmetic [`crate::material::Formation::of_matter`] does. What comes
    /// out is a grain, and grains are what a melt lays down.
    fn derive_frozen_layout(&mut self, idx: NodeIdx) -> bool {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return false;
        }
        if self.tree.nodes[idx.get()].morphology.is_some() {
            return false;
        }
        let (mixture, mass, radius, temperature) = {
            let n = &self.tree.nodes[idx.get()];
            (n.matter.mixture, n.matter.mass, n.matter.radius, n.matter.temperature)
        };
        if mixture.is_empty() || !(mass > 0.0) || !(radius > 0.0) {
            return false;
        }
        // More of it is solid than not: a melt that has only begun to freeze is
        // still a melt, and what it is laid out as is not settled until the
        // solid is the thing that is there.
        let solid = mixture
            .entries()
            .iter()
            .find(|p| p.phase == crate::chem::Phase::Solid && p.fraction > 0.5)
            .copied();
        let Some(solid) = solid else { return false };
        let Some(sub) = self.substances.get(solid.substance) else { return false };
        let props = sub.props;
        if !(props.melting_point > 0.0) || temperature >= props.melting_point {
            return false;
        }
        // Both rates a freezing needs. The front speed is set by how fast
        // latent heat can leave through the surface either way; only the
        // cooling rate has a measured alternative.
        let area = 4.0 * std::f64::consts::PI * radius * radius;
        let power = crate::units::SIGMA_SB * area * props.melting_point.powi(4);
        let atoms = props.atoms_per_unit.max(1) as f64;
        let heat_of_fusion = crate::units::K_B
            * crate::units::N_AVOGADRO
            * props.melting_point
            * atoms
            / (props.unit_mass * crate::units::N_AVOGADRO).max(1e-30);
        let specific_heat = 3.0 * crate::units::K_B * atoms / props.unit_mass.max(1e-30);
        let front = power / (area * (heat_of_fusion * props.density.max(1e-6)).max(1e-30));
        let cooling = self
            .cooling_rate_through(idx, props.melting_point)
            .filter(|r| *r > 0.0)
            .unwrap_or(power / (mass * specific_heat).max(1e-30));
        if !(cooling > 0.0) || !(front > 0.0) {
            return false;
        }
        let formation = crate::material::Formation::cooled(cooling, front);
        let grain = crate::material::grain_scale(&props, formation);
        if !(grain > 0.0) || !grain.is_finite() {
            return false;
        }
        let density = self.tree.nodes[idx.get()].matter.density().max(1e-9);
        let side = (mass / density).cbrt();
        if !(side > 0.0) || !side.is_finite() {
            return false;
        }
        let key = self.tree.nodes[idx.get()].key;
        let seed = self.tree.world_seed;
        let mut m = crate::morph::Morphology::new(crate::morph::Program::Terrain, seed, key.0, 0);
        m.built = mass;
        m.design_mass = mass;
        m.progress = 1.0;
        m.recipe = Some(crate::recipe::Recipe::Granular(crate::recipe::Granular {
            grain,
            density,
            side,
        }));
        let extent = m.extent();
        // Anything promoted out of the body list goes back into the node before
        // the list is discarded, or it is a live node nothing can reach.
        self.tree.shed_children(idx);
        let n = &mut self.tree.nodes[idx.get()];
        n.matter.radius = extent.max(1e-30);
        n.morphology = Some(m);
        n.bodies.clear();
        n.topology = None;
        self.tree.stats.structures += 1;
        self.tree.stats.layouts_derived += 1;
        true
    }

    /// Ask whether a node is a body its own gravity has rounded, and if it is,
    /// write down the surface it has.
    ///
    /// **Nothing here knows what a planet is.** `docs/PLAY.md` D6: "A planetary
    /// surface is not a new kind of thing." What makes a body round is that its
    /// own weight has overcome the strength of what it is made of, and both
    /// sides of that are measured: the central pressure of a self-gravitating
    /// sphere is `3 G M^2 / (8 pi R^4)`, and the strength is D14's, read off
    /// the node's own mixture at the node's own size. A 500 m asteroid comes
    /// out at 2.2e-4 Pa against rock's 10^8 and stays a potato; an Earth comes
    /// out at 1.7e11 and is a ball. The transition is around two hundred
    /// kilometres of rock, which is where it is.
    ///
    /// A node with no mixture is not refused on a guess: it has not been said
    /// what it is made of, so there is no fact of the matter about its surface,
    /// and it gets none. That is the second axiom rather than a gap — a
    /// scenario that wants ground says what its planet is made of, exactly as
    /// it already says what it is composed of.
    ///
    /// Solidity is deliberately *not* part of the gate. A magma ocean is still
    /// a surface and still tiles; whether anything can stand on it is the
    /// mixture's business, and `bake_surface` already answers it by presenting
    /// no solid for a node that holds none.
    pub fn assess_surface(&mut self, idx: NodeIdx) -> bool {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return false;
        }
        if self.tree.nodes[idx.get()].morphology.is_some() {
            return false;
        }
        // **How big it is now, not how big it was when it last coarsened.**
        // While a node is materialised its bodies are the authority and its
        // matter is the summary made when it last folded them back, so a
        // collapsing cloud's `matter.radius` is stale by as long as it has been
        // collapsing. Asking the stale one let a ball fall to a twelfth of the
        // density it should have been caught at before anything noticed.
        //
        // Measured the way `summarise` measures it, so the two cannot disagree.
        let (mass, radius, spec_count) = {
            let n = &self.tree.nodes[idx.get()];
            let radius = if n.bodies.is_empty() {
                n.matter.radius
            } else {
                let spread =
                    crate::state::Spread::of(n.bodies.iter().map(|b| (b.pos, b.mass, b.radius)));
                (spread.rms * crate::state::RMS_TO_RADIUS).max(1e-30)
            };
            (n.matter.mass, radius, n.spec.count)
        };
        if !(mass > 0.0) || !(radius > 0.0) {
            return false;
        }
        // **Its matter has to be packed densely enough to carry itself**, and
        // this is tested first because it is the half that refuses most nodes
        // and it costs a division, where deriving a material runs a nucleation
        // march.
        //
        // The pressure test alone is satisfied by a *cloud*, and trivially:
        // strength goes as the square of how much of the volume is filled, so a
        // ball of rock grains at half a kilogram a cubic metre has no strength
        // to overcome and is declared round while it is still falling in.
        // Measured: an Earth's worth of rock spread over 1.5e8 m passed the
        // pressure test at its first frame and froze there, thirty times too
        // large, with its collapse never run.
        //
        // What tells a planet from a cloud is that a planet is *condensed*: its
        // bulk density is the density of what it is made of. The threshold is
        // not a choice — it is where a pile of grains stops flowing and starts
        // carrying load, which is random loose packing, a universal of sphere
        // packing in the same family as Turnbull's 0.45.
        //
        // Read off the formation rather than off the material, because a
        // measured material's density already *is* the node's bulk density —
        // the packing is folded into it — so dividing one by the other is
        // identically one and answers nothing.
        let packed = self.packing_at(idx, radius).unwrap_or(0.0);
        if !(packed >= Self::RANDOM_LOOSE_PACKING) {
            return false;
        }
        // It is a body now, so its matter is brought into step with the detail
        // that made it one before anything is written down. `Tree::settle` is
        // `coarsen` without the destruction.
        self.tree.settle(idx);
        let Some(material) = self.material_of(idx) else { return false };
        let pressure =
            3.0 * crate::units::G * mass * mass / (8.0 * std::f64::consts::PI * radius.powi(4));
        if !(pressure > material.strength_of(radius, self.tree.nodes[idx.get()].matter.temperature))
        {
            return false;
        }

        // **And then it compacts.** What passed the test above is a pile whose
        // own weight exceeds what its grains can carry, so the grains give: the
        // pores close, and past full density the solid itself is squeezed until
        // its equation of state carries the weight. Without this a planet
        // stopped where its grains jammed — measured at 1743 kg/m^3 for an
        // Earth's worth of silicate, against the 2644 its own substance is.
        let radius = self.compacted_radius(idx, radius);

        let key = self.tree.nodes[idx.get()].key;
        let seed = self.tree.world_seed;
        let mut m = crate::morph::Morphology::new(crate::morph::Program::Terrain, seed, key.0, 0);
        m.built = mass;
        m.design_mass = mass;
        m.progress = 1.0;
        let mut relief = [0.0f32; 8];
        relief.copy_from_slice(&m.genome);
        // How many cells a patch divides into, from the resolution policy the
        // node already carries. Fixed in the recipe rather than read from
        // whoever is looking, because the address has to mean the same thing at
        // every level of detail.
        let cells = ((spec_count as f64).sqrt().floor() as usize).clamp(2, 32) as u8;
        m.recipe = Some(crate::recipe::Recipe::Tiled(crate::recipe::Tiled {
            density: (mass / (4.0 / 3.0 * std::f64::consts::PI * radius.powi(3))).max(1e-9),
            sphere: radius,
            face: crate::recipe::WHOLE_BALL,
            level: 0,
            u: 0,
            v: 0,
            cells,
            side: radius * crate::recipe::FACE_SIDE,
            depth: radius * crate::recipe::FACE_SIDE * crate::recipe::SLAB_ASPECT,
            relief,
        }));
        let n = &mut self.tree.nodes[idx.get()];
        // A ball's radius is its own, and the recipe says so. Stated rather
        // than taken from `extent()` so that acquiring a surface cannot move a
        // planet's radius by a hair and disturb everything standing on it.
        n.matter.radius = radius;
        n.morphology = Some(m);
        n.bodies.clear();
        n.topology = None;
        self.tree.stats.structures += 1;
        true
    }

    /// Squeeze a self-gravitating ball of condensed matter to the radius where
    /// its equation of state carries its own weight, and return that radius.
    ///
    /// **The balance is the virial theorem's**, which needs no density profile:
    /// for a body in hydrostatic equilibrium the volume-averaged pressure is
    /// `-W / 3V`, and the mean density is what carries it. The condensed
    /// phase's `p(rho)` is Murnaghan's (`eos.rs`), so the answer is where
    /// `p(M / V(r)) = -W(r) / 3V(r)`, and the left side rises as `r^-12`
    /// against the right's `r^-4`, so there is exactly one such radius below
    /// the one at rest density.
    ///
    /// **The binding follows homologously**: a contraction by `R/r` scales
    /// every separation by the same factor, so `W(r) = W(R) R / r` whatever
    /// the profile, and the energy released goes into the node's internal
    /// account, which is where the heat of a planet's formation belongs. The
    /// total is unchanged to rounding.
    ///
    /// Returns `radius` unchanged for anything whose matter is not condensed,
    /// and for anything the pressure does not reach full density.
    pub fn compacted_radius(&mut self, idx: NodeIdx, radius: f64) -> f64 {
        let Some(c) = self.eos_of(idx).condensed() else { return radius };
        let n = &self.tree.nodes[idx.get()];
        let mass = n.matter.mass;
        let binding = n.matter.gravitational_binding;
        if !(mass > 0.0) || !(binding < 0.0) || !(radius > 0.0) {
            return radius;
        }
        let volume = |r: f64| 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
        let excess = |r: f64| {
            let weight = -binding * radius / r / (3.0 * volume(r));
            c.pressure(mass / volume(r)) - weight
        };
        // At rest density the solid carries nothing, so the weight wins there;
        // bisect in the logarithm from there down.
        let at_rest = (mass / c.rest_density / (4.0 / 3.0 * std::f64::consts::PI)).cbrt();
        if at_rest >= radius {
            // Still porous when its own weight is carried: nothing to do here
            // that the packing test above did not already decide.
            return radius;
        }
        let (mut hi, mut lo) = (at_rest, at_rest * 1e-3);
        if !(excess(lo) > 0.0) {
            return radius;
        }
        for _ in 0..200 {
            let mid = (hi * lo).sqrt();
            if excess(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let r = (hi * lo).sqrt();
        let n = &mut self.tree.nodes[idx.get()];
        let squeezed = binding * radius / r;
        n.matter.internal_energy += binding - squeezed;
        n.matter.gravitational_binding = squeezed;
        n.matter.radius = r;
        r
    }

    /// Write the node's generated program against what is actually measured
    /// here.
    ///
    /// `Tree` has no registry and no environments, so it writes a recipe
    /// against the conditions it can see. This is where the measured ones
    /// arrive: the material read off the node's own mixture, and the fluid,
    /// light and flow the node is actually in. D11's whole claim is that those
    /// are measurements rather than columns — "a coral is in water because its
    /// node's mixture is water; nobody tells it" — so this is the line that
    /// makes it true.
    pub fn rewrite_recipe(&mut self, idx: NodeIdx) {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return;
        }
        let env = self.environment_at(idx);
        // **The material the thing is made of, and not the one the node
        // measures.** `Material::measured` reads a node's mixture and gets its
        // *packing* from the node's bulk density against the substance's — a
        // measurement that is right for a node which is the material (a plank,
        // a block of stone) and wrong for one that is an arrangement of it,
        // because a structure's volume is mostly the space between its parts.
        // Measured on the six-panel box: a crate of solid planks reads as 4.8%
        // packed, so its wood derives at 77 kg/m^3.
        //
        // And using it here would be *circular*: the recipe turns a mass into
        // a size using the density, the size becomes the node's radius, and
        // the radius is what the packing is measured from. The loop has every
        // size as a fixed point, which is another way of saying it decides
        // nothing.
        //
        // So a generator uses the material's own stated formation. D11's
        // density column is therefore still on `Program` and this is the
        // measurement that says why: **the porosity of a deposited solid is
        // not derivable from a node's bulk density once the node contains
        // void**, and the engine has no other way to ask.
        let material = match self.tree.nodes[idx.get()].morphology.as_ref() {
            Some(m) => m.program.material(),
            None => return,
        };
        let Some(m) = self.tree.nodes[idx.get()].morphology.as_mut() else {
            return;
        };
        // A parts list is not regenerated: somebody placed those parts, and
        // rewriting the rule would be the engine deciding it knew better.
        if m.is_assembled() {
            return;
        }
        m.regenerate(&env, &material);
        let extent = m.extent();
        let stored = m.stored_energy();
        let n = &mut self.tree.nodes[idx.get()];
        n.matter.radius = extent.max(1e-30);
        n.matter.chemical_energy = stored;
        self.tree.retier(idx);
    }

    /// Give a node a structure that is already there, at a stated mass. The
    /// counterpart to [`Self::plant`] for terrain and for anything that was
    /// standing before the world was looked at. See [`Tree::emplace`].
    pub fn emplace(
        &mut self,
        idx: NodeIdx,
        program: crate::morph::Program,
        built: f64,
        env: Option<crate::morph::Environment>,
    ) {
        let id = self.identify(idx);
        self.tree.emplace(idx, program, built);
        match env {
            Some(e) => {
                self.environments.insert(id, e);
            }
            None => {
                self.environments.remove(&id);
            }
        }
        self.rewrite_recipe(idx);
    }

    /// State a composite: one node whose recipe is the parts it is made of.
    ///
    /// `docs/PLAY.md` D15, the forward direction. See [`Tree::assemble`] for
    /// what the recipe is and why `program` is provenance rather than species.
    pub fn assemble(
        &mut self,
        idx: NodeIdx,
        program: crate::morph::Program,
        parts: crate::assembly::Assembly,
        env: Option<crate::morph::Environment>,
    ) {
        let id = self.identify(idx);
        self.tree.assemble(idx, program, parts);
        match env {
            Some(e) => {
                self.environments.insert(id, e);
            }
            None => {
                self.environments.remove(&id);
            }
        }
    }

    /// Attach a promoted child back into its parent's recipe. **D15's join.**
    ///
    /// Weld, glue, nail and grown-together are not four features: they differ
    /// only in what the join is made of, and therefore in its strength under
    /// D14's Griffith law. So this takes a substance and a contact area, and
    /// everything else about how hard it is to undo follows from the same
    /// measurement every other strength in the engine comes from.
    ///
    /// **Explicit, never automatic.** A part that has come off and is merely
    /// lying against the thing it came from stays its own node until something
    /// puts it back — a decision taken deliberately, because "touching and at
    /// rest" would re-absorb a fence panel the wind had just torn off, and a
    /// world that quietly reassembles itself is worse than one that leaves a
    /// pile of parts on the ground.
    ///
    /// The child must actually be a child: joining across the tree would put a
    /// load path over a node boundary, which is the gap `BACKLOG.md` calls
    /// "Structures cannot span promoted children" and D5's substructuring
    /// closes in **Bodies**.
    ///
    /// Returns the site name the part now has in the recipe.
    pub fn join(
        &mut self,
        composite: NodeIdx,
        child: NodeIdx,
        join: crate::chem::SubstanceId,
        area: f64,
    ) -> Option<u32> {
        if composite.is_none() || child.is_none() {
            return None;
        }
        if self.tree.nodes[child.get()].parent != composite {
            return None;
        }
        let slot = self.tree.nodes[child.get()].slot as usize;
        let (at, facing, mass, radius, substance, child_parts) = {
            let c = &self.tree.nodes[child.get()];
            let substance = c
                .matter
                .mixture
                .entries()
                .iter()
                .filter(|p| p.phase == crate::chem::Phase::Solid && p.fraction > 0.0)
                .max_by(|a, b| a.fraction.total_cmp(&b.fraction))
                .map(|p| p.substance)
                .unwrap_or(crate::chem::SubstanceId::UNSPECIATED);
            (
                c.motion.offset,
                // A part brings the way it is facing with it, which is the
                // whole of what an orientation is for: a plank joined on
                // sideways is a plank lying sideways.
                c.motion.orientation,
                c.matter.mass,
                c.matter.radius,
                substance,
                c.morphology.as_ref().and_then(|m| m.assembly().cloned()),
            )
        };

        // The part's shape is its own if it has one, and a sphere of its
        // bounding radius if it does not. A node with no recipe has no shape to
        // bring, and inventing a plate for it would be the engine deciding what
        // an undescribed thing looks like — which is the first axiom's line.
        // `Vec3::ZERO` half-extents *are* the sphere, so this states nothing.
        let half = match child_parts.as_ref().and_then(|a| a.parts.first()) {
            Some(p) => p.half,
            None => crate::math::Vec3::ZERO,
        };

        let Some(m) = self.tree.nodes[composite.get()].morphology.as_mut() else {
            return None;
        };
        if m.assembly().is_none() {
            m.set_assembly(Default::default());
        }
        let assembly = m.assembly_mut().expect("just set");
        let site = assembly.next_site();
        // Joined to whichever part it is actually touching — the nearest one,
        // measured, rather than to whatever happened to be added first.
        let to = assembly
            .parts
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| (a.pos - at).norm2().total_cmp(&(b.pos - at).norm2()))
            .map(|(i, _)| i as u16)
            .unwrap_or(crate::assembly::UNJOINED);
        let mut part = crate::state::Body::solid(at, half, mass, substance, site);
        part.orientation = facing;
        if !part.is_boxed() {
            part.radius = radius;
        }
        assembly.parts.push(part);
        assembly.joins.resize(assembly.parts.len(), crate::assembly::Join::NONE);
        if to != crate::assembly::UNJOINED {
            let last = assembly.parts.len() - 1;
            assembly.joins[last] = crate::assembly::Join::new(to, join, area);
        }
        m.built = assembly.mass();
        m.design_mass = m.built;

        // The child's matter is already counted in the composite's — a promoted
        // child is one of its parent's bodies seen at finer resolution, not an
        // addition to it — so absorbing it moves nothing. What goes is the
        // node, and the slot it was standing in.
        self.tree.release_subtree(child);
        if let Some(c) = self.tree.nodes[composite.get()].children.get_mut(slot) {
            *c = NodeIdx::NONE;
        }
        let radius = self.tree.nodes[composite.get()]
            .morphology
            .as_ref()
            .map(|m| m.extent())
            .unwrap_or(0.0);
        {
            // Any *other* part that had been promoted out of this composite is
            // folded back first: the recipe's slots are about to be renumbered
            // and a child pointing at one of them would be unreachable. See
            // `Tree::shed_children`.
            self.tree.shed_children(composite);
            let n = &mut self.tree.nodes[composite.get()];
            n.matter.radius = radius.max(1e-30);
            // The parts moved, so the body list they generate has, and the
            // joins are what the solver reads next frame.
            n.bodies.clear();
            n.topology = None;
        }
        self.tree.retier(composite);
        self.tree.bump_epoch(composite);
        self.tree.record_edit(composite);
        self.disturb(composite);
        self.tree.stats.joins += 1;
        Some(site)
    }

    /// One part comes away and becomes a node of its own, at that moment.
    ///
    /// **D15's break, which is the same transform run backwards.** The part
    /// keeps its slot in the recipe and gains a node — a promoted child, which
    /// is what this engine means by "a thing of its own" — carrying its own
    /// single-part recipe, so a wall that comes off a box is still a wall and
    /// not an anonymous lump of mass. See [`crate::assembly::Assembly::break_join`]
    /// for why the recipe keeps six parts and five joins rather than shrinking
    /// to five.
    ///
    /// Returns the new node, or `NONE` if there was no such part or nothing
    /// holding it on.
    pub fn detach(&mut self, composite: NodeIdx, site: u32) -> NodeIdx {
        if composite.is_none() || !self.tree.nodes[composite.get()].alive {
            return NodeIdx::NONE;
        }
        let Some(slot) = self.tree.nodes[composite.get()]
            .morphology
            .as_ref()
            .and_then(|m| m.assembly())
            .and_then(|a| a.index_of_site(site))
        else {
            return NodeIdx::NONE;
        };
        // The bodies have to exist for there to be a slot to promote: the parts
        // *are* the body list, in recipe order, which is what makes the site
        // and the slot the same thing.
        self.tree.refine(composite);
        if slot >= self.tree.nodes[composite.get()].bodies.len() {
            return NodeIdx::NONE;
        }

        let broke = {
            let Some(m) = self.tree.nodes[composite.get()].morphology.as_mut() else {
                return NodeIdx::NONE;
            };
            let Some(a) = m.assembly_mut() else { return NodeIdx::NONE };
            let broke = a.break_join(site);
            if broke {
                let taken = a.taken(&[site]);
                let at = m.age;
                m.events.push(crate::morph::Event {
                    at,
                    kind: crate::morph::EventKind::Severed,
                    site,
                    magnitude: 0.0,
                });
                Some(taken)
            } else {
                None
            }
        };
        let Some(taken) = broke else { return NodeIdx::NONE };

        let program = self.tree.nodes[composite.get()]
            .morphology
            .as_ref()
            .map(|m| m.program)
            .unwrap_or(crate::morph::Program::Wall);
        let spec = self.tree.nodes[composite.get()].spec;
        let child = self.tree.promote(composite, slot, spec);
        if child.is_none() {
            return NodeIdx::NONE;
        }
        // The piece is what it was: the same solid, the same material, now
        // about its own centre.
        self.tree.assemble(child, program, taken);
        self.identify(child);
        self.tree.record_edit(composite);
        self.disturb(composite);
        self.tree.stats.detachments += 1;
        child
    }

    /// The name of whatever lives at this address, if it has been given one.
    ///
    /// Lookup only. A node nothing has recorded against has no identity, and
    /// asking does not create one — which is what keeps the index the size of
    /// the world's *history* rather than the size of its tree.
    pub fn identity_of(&self, key: PathKey) -> Option<EntityId> {
        self.identities.get(&key).copied()
    }

    /// The name of whatever lives at this address, issuing one if it has none.
    ///
    /// Called on the write path of anything that keys a side table: giving a
    /// node chemistry, an environment, a clock or a history is the moment it
    /// becomes a thing worth naming. Reading those tables uses
    /// [`Self::identity_of`] instead, so a lookup never grows the index.
    pub fn issue_identity(&mut self, key: PathKey) -> EntityId {
        if let Some(id) = self.identities.get(&key) {
            return *id;
        }
        let id = EntityId(self.next_entity);
        self.next_entity += 1;
        self.identities.insert(key, id);
        id
    }

    /// The identity of a live node, issuing one if it has none.
    pub fn identify(&mut self, idx: NodeIdx) -> EntityId {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return EntityId::NONE;
        }
        let key = self.tree.nodes[idx.get()].key;
        self.issue_identity(key)
    }

    /// The identity of a live node, without issuing one.
    pub fn identity(&self, idx: NodeIdx) -> Option<EntityId> {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return None;
        }
        self.identity_of(self.tree.nodes[idx.get()].key)
    }

    pub fn add_observer(&mut self, o: Observer) -> usize {
        self.observers.push(o);
        self.observers.len() - 1
    }

    // -----------------------------------------------------------------
    // the frame
    // -----------------------------------------------------------------

    /// Advance the world by one frame.
    ///
    /// `wall_us` is what the caller is prepared to spend. The engine will spend
    /// up to that and no more; if the work does not fit, detail is dropped and
    /// `stats.detail_debt` records how much value went unserved.
    ///
    /// The frame has three passes, and the order is the design:
    ///
    /// 1. **Survey and solve.** Every node's *lateness* — how long since it was
    ///    last re-solved, in units of its own characteristic time — is computed
    ///    against the new instant. Nodes at or past one are due, and the budget
    ///    takes the most overdue that fit. Nothing is skipped for being off
    ///    screen; a node passed over grows more overdue and wins next frame.
    /// 2. **Coast.** Every node that was not solved is carried to the same
    ///    instant in closed form. Position under constant velocity and
    ///    orientation under constant spin are exact solutions, so this costs
    ///    one add per node and introduces no error at all.
    /// 3. **Commit.** The world instant moves, and every live node is at it.
    ///
    /// The pass that used to be here — decide what the observer can see, and
    /// only simulate that — is gone. It was the wrong organising principle: a
    /// tree falls whether or not anyone is pointing at it. What the observer
    /// still controls is *resolution*, which feeds back into the cadence on its
    /// own, because a node represented more finely comes due more often.
    pub fn step_frame(&mut self, wall_us: f64) -> Plan {
        let t0 = std::time::Instant::now();
        self.budget.target_us = wall_us;
        self.refresh_pace();

        let horizon = self.time + self.frame_dt();
        // A body that has become one gets its surface, and a surface is
        // generated by being approached and by nothing else.
        self.derive_layouts();
        self.approach();
        self.hold();
        let tasks = self.survey(horizon);
        let bytes = self.tree.detail_bytes();
        let plan = self.budget.plan(tasks, bytes);
        let achieved = self.execute(&plan, horizon);
        let coasted = self.coast_to(horizon);
        // Everything has moved by now — solved nodes in `execute`, everything
        // else in `coast_to` — so this is the first moment at which "where is
        // it" has one answer for the whole world. `docs/PLAY.md` D16.
        self.cross_boundaries();
        // And the same measurement one level down: a node whose contents have
        // stopped being one neighbourhood is describing two places at once.
        self.resolve_extents(&plan);
        // Both of those can make a node out of a body, and a node made from a
        // body is where that body was at its parent's contents instant — which
        // may be behind the world's (`Node::carried`). Carried here, so that
        // every node the frame ends with is at the frame's instant and not
        // only the ones that existed when it started. Measured without it: an
        // Atomic node made in the last frame's crossing or split passes stood
        // at 2.9e-17 s in a world at 1 s.
        for i in 0..self.tree.nodes.len() {
            if self.tree.nodes[i].alive {
                self.tree.carry(NodeIdx(i as u32), horizon);
            }
        }

        self.deliver_influences(horizon);
        // Chemistry runs on the span the frame actually covered, after the
        // influences that changed the temperatures it reads. It costs at most
        // eight comparisons per node that is made of something, and nothing at
        // all for every node that is not.
        // Matter evolution before chemistry, and the order is not arbitrary:
        // the thermal step decides what temperature the phase step is deciding
        // at. Warming a node and then asking whether its ice has melted is the
        // right way round; asking first and warming after delays every thaw by
        // a frame.
        let span = horizon - self.time;
        let evolution = self.evolve_matter(span);
        self.stats.radiated += evolution.radiated;
        self.stats.absorbed += evolution.absorbed;
        self.stats.radiation_deficit += evolution.radiation_deficit;
        let chemistry = self.react_all(span);
        // What has happened to a surface feels its own physics, on the world
        // clock and not on whether anybody is looking. `docs/PLAY.md` §5.
        self.weather(span);
        // And a planet's ocean answers what is outside it. `ocean.rs`.
        self.tides(span);
        self.time = horizon;
        self.record_histories();

        self.retime(achieved);
        let (worst_lateness, overdue) = self.lateness_report();
        let actual = t0.elapsed().as_secs_f64() * 1e6;
        self.budget.observe_frame(plan.planned_us, actual);
        self.stats.frames += 1;
        self.stats.sim_time = self.time;
        self.stats.tasks_run += plan.accepted.len() as u64;
        self.stats.tasks_deferred += plan.deferred as u64;
        self.stats.detail_debt = plan.unmet_value;
        self.stats.last_frame_us = actual;
        self.stats.materialised_bodies = self.tree.materialised_bodies();
        self.stats.live_nodes = self.tree.live_count();
        self.stats.worst_lateness = worst_lateness;
        self.stats.overdue = overdue;
        self.stats.coasted = coasted;
        self.stats.dissolved += chemistry.dissolved;
        self.stats.precipitated += chemistry.precipitated;
        self.stats.phase_changed += chemistry.melted
            + chemistry.frozen
            + chemistry.boiled
            + chemistry.condensed;
        self.stats.reacting_nodes = self
            .tree
            .nodes
            .iter()
            .filter(|n| n.alive && !n.matter.mixture.is_empty())
            .count();
        self.stats.bubbled = self.bubbles().len();
        plan
    }

    /// How much world time one frame covers.
    ///
    /// This used to be the smallest timestep any node in the world wanted, and
    /// that one line is why the engine could not run every scale at once: it
    /// made a resolved nucleus set the pace for the galaxy containing it. It is
    /// now a *rate* — how fast the user has asked simulated time to run —
    /// because a node no longer has to be stepped at the world's cadence. It
    /// only has to arrive at the world's instant, and most nodes get there in
    /// closed form.
    ///
    /// The one hard cap left is causality: no node may be advanced past the
    /// arrival of an influence, or the influence would land in its past.
    pub fn frame_dt(&self) -> f64 {
        let mut dt = self.pace * self.time_rate * self.applied_throttle();
        if !(dt > 0.0) || !dt.is_finite() {
            dt = 1e-30;
        }
        if let Some(next) = self.mailbox.next_arrival() {
            dt = dt.min((next - self.time).max(0.0).max(1e-30));
        }
        dt
    }

    /// Adjust how much simulated time the next frame may cover.
    ///
    /// The frame rate is the invariant, so when the work does not fit something
    /// has to give, and there are two different things that can: how much
    /// simulated time passes, and how much of the world is resolved. This
    /// controls the first, and it has to be careful not to answer a question
    /// that belongs to the second.
    ///
    /// The signal is *shortfall on work that actually ran*: a node the frame
    /// accepted, and integrated, and which still did not reach the instant. That
    /// says the span was too long, and shortening it fixes it. Work that was
    /// deferred entirely says something else — that more of the world is
    /// resolved than can be simulated — and slowing time does not help, it just
    /// stops the clock while the debt stays exactly where it was. That is the
    /// coarsener's problem, and it shows up as detail debt and lateness.
    ///
    /// Bounded below, because a world that has slowed by a factor of a million
    /// has said everything it can say by slowing further, and "frozen" is a
    /// worse answer than "slow, with some of it stale".
    /// The throttle as it actually bears on the clock.
    ///
    /// One in [`PaceMode::Fixed`], whatever the stored value is: a world runs at
    /// one second per second and overload makes it staler rather than slower.
    /// The stored value is left alone rather than reset, so a session that
    /// switches to `pace_to` and back does not lose what it had learned.
    fn applied_throttle(&self) -> f64 {
        match self.pace_mode {
            PaceMode::Fixed => 1.0,
            PaceMode::Follow => self.time_throttle,
        }
    }

    fn retime(&mut self, achieved: f64) {
        // Nothing to learn in `Fixed`: the number would decay against a clock it
        // is not allowed to move, and then bite the moment somebody paced to a
        // node.
        if self.pace_mode == PaceMode::Fixed {
            return;
        }
        if achieved < 0.95 {
            self.time_throttle =
                (self.time_throttle * achieved.clamp(0.05, 1.0)).max(MIN_TIME_THROTTLE);
        } else if self.time_throttle < 1.0 {
            // Recover slowly. Being late is visible and being early is not.
            self.time_throttle = (self.time_throttle * 1.5).min(1.0);
        }
    }

    /// Run the world at the pace of a particular node: one frame covers about
    /// one of its characteristic times.
    ///
    /// This is what makes the same engine watchable at both ends of the ladder.
    /// Paced to a galaxy, a frame is a few thousand years and the spiral turns;
    /// paced to a carbon atom, a frame is a femtosecond and the bonds vibrate.
    /// Neither choice changes the physics — only how much of it happens between
    /// two pictures.
    pub fn pace_to(&mut self, idx: NodeIdx) {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return;
        }
        self.pace_mode = PaceMode::Follow;
        self.paced_to = idx;
        self.refresh_pace();
    }

    /// One second of world time per second of wall time, which is what a world
    /// is. `docs/PLAY.md` D1.
    ///
    /// The span a frame covers is therefore the frame's own length: at twenty
    /// updates a second, fifty milliseconds — which is the arithmetic §3.4
    /// states and the number its resolution-floor table is computed at.
    ///
    /// This was got wrong when D1 was first built: the pace was left at the
    /// `1.0` the field is initialised with, which is one second *per frame* and
    /// so twenty times real time at twenty updates a second. It read as "one
    /// second per second" and was not. `World::resolution_floor` is what caught
    /// it — the floor came out twenty times coarser than §3.4's table, and the
    /// formula was right.
    ///
    /// It reads the frame rate once. A caller that changes it — `step_frame`
    /// takes a wall budget every call — asks again, and the budget being spent
    /// on one frame is deliberately *not* what sets the clock: a frame that is
    /// given less time to work in advances the same span and gets staler, which
    /// is the trade D1 makes.
    pub fn pace_realtime(&mut self) {
        self.pace_mode = PaceMode::Fixed;
        let span = self.budget.target_us * 1e-6;
        self.pace = if span > 0.0 && span.is_finite() { span } else { 1.0 };
    }

    /// Drive the clock by hand: `span` seconds of world time per frame,
    /// regardless of what anything is doing.
    ///
    /// For growth demos, scripted scenarios, replay and tests — anything that
    /// needs the clock to be an input rather than a consequence. It is the
    /// honest counterpart to [`World::pace_to`], and the two are mutually
    /// exclusive by construction rather than by convention.
    pub fn pace_fixed(&mut self, span: f64) {
        if !(span > 0.0) || !span.is_finite() {
            return;
        }
        self.pace_mode = PaceMode::Fixed;
        self.pace = span;
    }

    /// Re-read the pace from the node it is taken from.
    ///
    /// Does nothing in [`PaceMode::Fixed`], which is how an assignment to
    /// `pace` sticks. It used to be that `paced_to = NodeIdx::NONE` was the
    /// incantation for the same thing, which worked by accident of an early
    /// return and read as a bug at every call site.
    ///
    /// Called at the top of every frame, because the subject's cadence moves:
    /// the moment a galaxy is materialised into its stars, the span a frame may
    /// cover has to come down with it or the stars cannot be integrated across
    /// one. This is the feedback that makes "zooming in slows time" an
    /// arithmetic consequence rather than a policy.
    pub fn refresh_pace(&mut self) {
        if self.pace_mode == PaceMode::Fixed {
            return;
        }
        let idx = self.paced_to;
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return;
        }
        let tau = self.node_pace(idx);
        if tau.is_finite() && tau > 0.0 {
            self.pace = tau;
        }
    }

    /// How fast this node's interior runs, per second of world time.
    ///
    /// The product, up the chain of frames, of every level's relativistic
    /// dilation and every level's administrative bubble — which is the chain
    /// rule, not an approximation, because each node's frame velocity is
    /// measured in its parent's frame. See `dilation.rs` for why the three
    /// factors live in one product and are reported separately.
    ///
    /// Costs one multiply per level of depth, and depth is about twelve at the
    /// bottom of the ladder.
    pub fn time_rate_of(&self, idx: NodeIdx) -> crate::dilation::TimeRate {
        let mut rate = crate::dilation::TimeRate::default();
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return rate;
        }
        let mut cur = idx;
        let mut guard = 0;
        while !cur.is_none() && guard < 64 {
            let n = &self.tree.nodes[cur.get()];
            let mut here = crate::dilation::physical_rate(
                n.motion.velocity,
                n.matter.external_potential,
                n.matter.mass,
            );
            // Through the gate, not straight from the field. `Node::bubble` is
            // public and a value written directly would otherwise be reported
            // here and not applied by `total()`.
            here.bubble = crate::dilation::accept_bubble(n.bubble).unwrap_or(1.0);
            rate = rate.compose(here);
            cur = n.parent;
            guard += 1;
        }
        rate
    }

    /// Local seconds per coordinate second, as one number.
    pub fn local_rate(&self, idx: NodeIdx) -> f64 {
        self.time_rate_of(idx).total()
    }

    /// How much world time a frame should cover, for someone watching `idx`.
    ///
    /// # Why this is not the cadence
    ///
    /// It was, and that was a bug worth writing down, because the two questions
    /// look identical and are not.
    ///
    /// [`World::node_cadence`] asks *may this representation go stale* — and for
    /// a node held as matter alone the honest answer can be "not for a
    /// very long time". A ball of ten-thousand-kelvin hydrogen has no bulk
    /// motion in its own rest frame (`promote` sets the momentum to zero, which
    /// is what a rest frame means), barely spins, and is not being stirred, so
    /// nothing a node's matter reports about it — mass, radius, temperature —
    /// changes at all. Its measured cadence was 5.2x10^18 seconds. That is not
    /// wrong: a hundred and sixty billion years is genuinely how long that
    /// *description* stays accurate.
    ///
    /// The pace asks something else: *how fast should the clock run for someone
    /// looking at this*. Answering it with the cadence let a client watching an
    /// unmaterialised node set the world clock to 2.6x10^11 seconds a frame
    /// instead of 6.2 — a factor of four times ten to the tenth, decided by
    /// whether the thing being watched happened to have been materialised yet.
    /// Recipes made that reachable from outside: a client can now be looking at
    /// scenery the engine never built.
    ///
    /// So the pace is additionally bounded by how long the node's *interior*
    /// takes to rearrange, which is one resolution element at the internal
    /// random speed. [`Matter::velocity_dispersion`] is that speed and is
    /// floored at the thermal speed, so it is never zero for warm matter — it
    /// reported 8.2 km/s for the node above, giving four days rather than a
    /// hundred and sixty billion years.
    ///
    /// For a materialised node the cadence is already the shorter of the two
    /// and the bound does nothing, which is the property that keeps "zooming in
    /// slows time" an arithmetic consequence: time slows when detail is
    /// *resolved*, not merely when something small is pointed at.
    pub fn node_pace(&self, idx: NodeIdx) -> f64 {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() || !self.tree.nodes[idx.get()].alive
        {
            return f64::INFINITY;
        }
        let cadence = self.node_cadence(idx);
        let n = &self.tree.nodes[idx.get()];
        let churn = n.matter.velocity_dispersion();
        if !(churn > 0.0) || !churn.is_finite() {
            return cadence;
        }
        cadence.min(self.node_resolution(idx) / churn)
    }

    /// The length scale a node is currently represented at, metres.
    ///
    /// Not the node's radius — the size of the smallest thing it is currently
    /// showing. A planet held as matter alone is represented at its own
    /// radius; the same planet split into four thousand parcels is represented
    /// at four hundred kilometres, and has to be re-solved sixteen times as
    /// often for the difference to mean anything.
    pub fn node_resolution(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        let parts = n.bodies.len();
        if parts > 1 {
            n.matter.radius / (parts as f64).cbrt()
        } else {
            n.matter.radius
        }
    }

    /// How long this node may be left alone before its state is visibly stale.
    ///
    /// For a node held as matter alone this is [`Matter::characteristic_time`]
    /// — how long before the thing moves, turns, or rearranges by its own size.
    /// For a materialised node it is the bodies that are represented, so it is
    /// the bodies that set the cadence: the time for the fastest of them to
    /// cross one resolution element. Thermal motion counts in that case and not
    /// in the first, and the difference is not a fudge — at matter resolution
    /// thermal motion is a temperature, and at parcel resolution it is
    /// something you can watch happen.
    pub fn node_cadence(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        if n.is_materialised() {
            let h = self.node_resolution(idx);
            // The bodies' own speeds, and nothing else. Seeding this with the
            // matter's sound speed looked harmless and was not: a galaxy's
            // "sound speed" is a gas-pressure formula applied to a collisionless
            // stellar system, and it saturated at 0.577c, so a materialised
            // galaxy claimed to need re-solving four orders of magnitude more
            // often than its own stars could justify. Sampling samples
            // thermal motion into the bodies already, so where a sound speed is
            // meaningful it is in here anyway.
            let v = n.bodies.iter().map(|b| b.vel.norm()).fold(0.0f64, f64::max);
            let moving = if v > 0.0 { h / v } else { f64::INFINITY };
            // **And anything out of balance, which speeds alone cannot see.**
            // Contents at rest under a net force cover a resolution element in
            // `sqrt(2 h / a)` from standing still. Measured without this: a
            // ball put in a bucket of water was coasted for four seconds and
            // the bucket never solved once, because nothing in it was moving
            // yet. Not yet measured is infinite, which is due now.
            let pushed = if n.unrest > 0.0 { (2.0 * h / n.unrest).sqrt() } else { f64::INFINITY };
            // **A thing standing in a liquid is coupled to it at the liquid's
            // signal speed**, through a wall whose spring rings in `gap / c`
            // (`hydro::Wall`). The child moves on its own clock between the
            // parent's solves, so a liquid coasted while a child in it moves is
            // a child sinking into water nobody is solving. Measured: a ball
            // coasted two frames into a coasted bucket, and the next solve
            // found parcels a spacing inside it and threw it out at 14.5 m/s.
            let coupled = if n.matter.mixture.in_phase(crate::chem::Phase::Liquid) > 0.0
                && n.children.iter().any(|c| !c.is_none())
            {
                let c = self.signal_speed_of(idx);
                if c > 0.0 { h / c } else { f64::INFINITY }
            } else {
                f64::INFINITY
            };
            return moving.min(pushed).min(coupled);
        }
        n.matter.characteristic_time(n.matter.radius)
    }

    /// How many of its own characteristic times a node has gone unsolved.
    ///
    /// One means it has just come due. Ten means the world has moved on nine
    /// cadences without asking it what it is doing, which is what "stale" means
    /// in a way that is the same number for a nucleus and a galaxy.
    pub fn lateness(&self, idx: NodeIdx, horizon: f64) -> f64 {
        let cadence = self.node_cadence(idx);
        if !(cadence > 0.0) {
            return 1e9;
        }
        if !cadence.is_finite() {
            return 0.0;
        }
        // In the node's own time, not the world's. The cadence is a local
        // quantity — a resolution element at the speeds measured inside the
        // node — so the elapsed span has to be local too. A node in a hundred-
        // fold bubble has lived a hundred times as long since it was last
        // solved, is a hundred times as late, and is ranked accordingly.
        let elapsed = (horizon - self.tree.nodes[idx.get()].last_solved) * self.local_rate(idx);
        (elapsed / cadence).clamp(0.0, 1e9)
    }

    /// How many sub-steps this node needs to cross one frame.
    ///
    /// Capped at [`MAX_SUBSTEPS`], beyond which following the trajectory is not
    /// merely expensive but the wrong answer, and the node is crossed by its
    /// ensemble instead.
    pub fn substeps(&self, idx: NodeIdx) -> u32 {
        let h = self.node_dt(idx);
        if !(h > 0.0) || !h.is_finite() {
            return 1;
        }
        let span = self.frame_dt();
        ((span / h).ceil().clamp(1.0, MAX_SUBSTEPS as f64)) as u32
    }

    /// The finest a node of this signal speed can be resolved and still be
    /// *followed*, at the pace the world is currently running. Metres.
    ///
    /// `docs/PLAY.md` §3.4 derives it by rearranging `node_dt`: a node is
    /// followed while `frame_span <= node_dt * MAX_SUBSTEPS`, and the term that
    /// binds inside `Continuum` is `0.25 * h / c_signal`, so
    ///
    /// ```text
    ///     h >= 4 * frame_span * c_signal / MAX_SUBSTEPS
    /// ```
    ///
    /// At the 50 ms frame a world runs at, that is `c_signal / 1280`:
    ///
    /// ```text
    ///     air    340 m/s     0.27 m
    ///     water 1500 m/s     1.2 m
    ///     rock  5000 m/s     3.9 m
    /// ```
    ///
    /// The floor is a property of the *material* rather than of the tier, which
    /// is why it is asked this way round. It also only became visible when the
    /// clock stopped following the observer: zooming in used to slow time until
    /// the substeps fit, which is D1 doing what it was chosen to do — landing
    /// the cost somewhere honest instead of in a silently slower world.
    pub fn resolution_floor(&self, signal_speed: f64) -> f64 {
        if !(signal_speed > 0.0) || !signal_speed.is_finite() {
            return 0.0;
        }
        4.0 * self.frame_dt() * signal_speed / MAX_SUBSTEPS as f64
    }

    /// The floor for what a particular node is made of, measured from its own
    /// contents rather than assumed.
    pub fn resolution_floor_of(&self, idx: NodeIdx) -> f64 {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return 0.0;
        }
        let n = &self.tree.nodes[idx.get()];
        let flow = n.bodies.iter().map(|b| b.vel.norm()).fold(0.0f64, f64::max);
        self.resolution_floor(self.signal_speed_of(idx).max(flow))
    }

    /// The equation of state a node's matter answers to, from what it is made
    /// of. See `eos.rs`.
    pub fn eos_of(&self, idx: NodeIdx) -> crate::eos::Eos {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return crate::eos::Eos::Gas;
        }
        crate::eos::Eos::of_matter(&self.tree.nodes[idx.get()].matter, &self.substances)
    }

    /// Speed at which a disturbance crosses a node's matter, m/s.
    ///
    /// `Matter::sound_speed` for a gas and for anything nobody has described,
    /// and the condensed phase's own for a liquid or a solid. §4.2 measured
    /// what the gas law says about water — 820 m/s, capped by the velocity
    /// dispersion — where the engine's own water, priced by its vaporisation
    /// energy, carries sound at 1569 m/s against a real 1500.
    pub fn sound_speed_of(&self, idx: NodeIdx) -> f64 {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return 0.0;
        }
        self.eos_of(idx).sound_speed(&self.tree.nodes[idx.get()].matter)
    }

    /// Speed at which a disturbance crosses a node's contents *as they are
    /// being solved*, m/s: [`World::sound_speed_of`], except for a liquid,
    /// which weakly-compressible SPH runs at ten times its flow scale
    /// (`flow_scale`), capped at its own sound speed.
    ///
    /// **The scheduler has to ask this and not the physical speed**, because it
    /// is what sets the step the solver can take. Measured with the physical
    /// one: a bucket of 1200 parcels was priced at 1569 m/s, `advance_to` cut
    /// each frame into 6765 calls of a few microseconds, and a frame took 21 s
    /// of which the solve needed 0.4. It is also §3.7's second response to the
    /// resolution floor, which is the reason WCSPH was chosen: the floor for a
    /// liquid moves with the speed the liquid is actually integrated at.
    pub fn signal_speed_of(&self, idx: NodeIdx) -> f64 {
        let physical = self.sound_speed_of(idx);
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return physical;
        }
        let n = &self.tree.nodes[idx.get()];
        if !(n.matter.mixture.in_phase(crate::chem::Phase::Liquid) > 0.0) {
            // A solid carries no signal a fluid solver has to follow: packed, it
            // is rigid and not solved as a fluid at all (`advance_node`);
            // below its rest density it presses with nothing and its parcels
            // move ballistically. §3.5: solids escape the floor.
            if !n.matter.gas_law_applies() && self.eos_of(idx).condensed().is_some() {
                return 0.0;
            }
            return physical;
        }
        if n.bodies.is_empty() {
            return physical;
        }
        let mask = n.structural_mask();
        let eos = if mask.is_some() {
            crate::eos::Eos::of_loose(&n.matter.mixture, &self.substances)
        } else {
            crate::eos::Eos::of_matter(&n.matter, &self.substances)
        };
        let Some(c) = eos.condensed() else { return physical };
        let field = if mask.is_some() { n.gravity } else { crate::math::Vec3::ZERO };
        let loose = n.bodies.iter().enumerate().filter(|(i, _)| {
            !mask.as_ref().map(|m| m.get(*i).copied().unwrap_or(false)).unwrap_or(false)
                && n.children.get(*i).map(|c| c.is_none()).unwrap_or(true)
        });
        let first = loose.clone().next().map(|(_, b)| b.mass).unwrap_or(0.0);
        let spacing = if c.rest_density > 0.0 { (first / c.rest_density).cbrt() } else { 0.0 };
        let drivers = n
            .children
            .iter()
            .filter(|c| !c.is_none())
            .map(|c| self.tree.nodes[c.get()].motion.velocity.norm())
            .fold(0.0f64, f64::max);
        let scale = flow_scale(loose.map(|(_, b)| b), field, spacing, drivers);
        (10.0 * scale).min(c.sound_speed(c.rest_density))
    }

    /// May this node be resolved at the pace the world is currently running?
    ///
    /// Resolving something the frame cannot integrate is not detail, it is a
    /// picture redrawn from scratch every frame — the node would thermalise
    /// immediately and lose whatever the resolution was for. Watching molecules
    /// vibrate requires slowing time down, and this is where the engine says so
    /// rather than pretending otherwise.
    pub fn can_resolve(&self, idx: NodeIdx) -> bool {
        let h = self.node_dt(idx);
        h.is_finite() && h > 0.0 && self.frame_dt() <= h * MAX_SUBSTEPS as f64
    }

    /// How long before this node's detail stops being *this* state and becomes
    /// merely *a* state of the same matter.
    ///
    /// The persistence rule, and the replacement for "discard it when nobody is
    /// looking". Detail is kept because something happened in it, and released
    /// once the node has had time to forget — at which point remembering it
    /// buys nothing and regenerating it costs nothing, because a fresh
    /// maximum-entropy draw and the stored sample are the same distribution.
    ///
    /// Three cases, and the differences between them are physical:
    ///
    /// * **Structured matter never forgets.** A tree that lost a branch has
    ///   lost it; no amount of waiting rearranges wood into a state the sampler
    ///   could have drawn. Anything carrying a morphology or a topology, and
    ///   anything somebody has touched, mixes in infinite time.
    /// * **Ordered motion does not mix.** A rotating disc keeps its stars in
    ///   their lanes however fast they are going, so the mixing speed is the
    ///   random part of the internal motion with the rotation taken out.
    /// * **Everything else mixes on a crossing time** — the time for a
    ///   constituent moving at that random speed to cross one resolution
    ///   element. A gas parcel forgets in milliseconds; a galaxy never.
    pub fn mixing_time(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        if n.pinned || n.morphology.is_some() || n.topology.is_some() {
            return f64::INFINITY;
        }
        // Anything the world has built structure on top of is remembered too.
        // A node's promoted children are not a sample of it — they are specific
        // objects with their own histories, and coarsening the parent releases
        // the whole subtree beneath it. One galaxy arm quietly forgetting would
        // take the star, the planet and the nucleus somebody was looking at
        // with it.
        if n.children.iter().any(|c| !c.is_none()) {
            return f64::INFINITY;
        }
        let dispersion = n.matter.velocity_dispersion();
        let ordered = n.matter.angular_velocity().norm() * n.matter.radius;
        let random2 = dispersion * dispersion - ordered * ordered;
        if !(random2 > 0.0) {
            return f64::INFINITY;
        }
        let h = self.node_resolution(idx);
        let v = random2.sqrt();
        if v > 0.0 {
            h / v
        } else {
            f64::INFINITY
        }
    }

    /// Record that something happened here, so the node's detail is kept for a
    /// mixing time rather than released at the next opportunity.
    pub fn disturb(&mut self, idx: NodeIdx) {
        if idx.is_none() {
            return;
        }
        let t = self.time;
        let mut cur = idx;
        // Upwards as well: a changed child means the parent's sample no longer
        // describes what is inside it either.
        while !cur.is_none() {
            let n = &mut self.tree.nodes[cur.get()];
            n.last_disturbed = t;
            cur = n.parent;
        }
    }

    fn lateness_report(&self) -> (f64, usize) {
        let mut worst = 0.0f64;
        let mut overdue = 0;
        for i in 0..self.tree.nodes.len() {
            let idx = NodeIdx(i as u32);
            if !self.tree.nodes[i].alive {
                continue;
            }
            let l = self.lateness(idx, self.time);
            if l >= 1.0 {
                overdue += 1;
            }
            worst = worst.max(l);
        }
        (worst, overdue)
    }

    /// Build the frame's candidate work list.
    fn survey(&mut self, horizon: f64) -> Vec<Task> {
        let mut tasks = Vec::new();
        let live: Vec<NodeIdx> = (0..self.tree.nodes.len())
            .map(|i| NodeIdx(i as u32))
            .filter(|i| self.tree.nodes[i.get()].alive)
            .collect();

        for idx in live {
            let (tier, materialised, count, radius) = {
                let n = &self.tree.nodes[idx.get()];
                (n.tier, n.is_materialised(), n.bodies.len(), n.matter.radius)
            };

            // What the observers want is *resolution*, not whether the physics
            // runs. This loop decides how finely the node is drawn; the cadence
            // below then follows from that on its own.
            let mut acuity = 0.0f64;
            let mut urgency = 0.0f64;
            let mut wanted_tier = tier;
            for obs in &self.observers {
                let sep = self
                    .tree
                    .separation(obs.anchor, obs.offset, idx, Vec3::ZERO);
                let d = sep.value.norm().max(1e-30);
                if !self.gate.reaches(d) {
                    // Outside the light cone: nothing that happens here can
                    // reach the observer within the horizon, so there is no
                    // resolution to be gained by drawing it finer.
                    continue;
                }
                let theta = crate::coords::angular_size(radius, d);
                let s = (theta / obs.angular_resolution).min(1e6) * obs.priority;
                if s > acuity {
                    acuity = s;
                    wanted_tier = obs.required_tier(d);
                }
                urgency = urgency.max(self.gate.urgency(d));
            }
            let residency = if acuity > 1.0 {
                Residency::Observed
            } else if urgency > 0.0 {
                Residency::Causal
            } else {
                Residency::Speculative
            };
            {
                let n = &mut self.tree.nodes[idx.get()];
                if n.residency != Residency::Pinned {
                    n.residency = residency;
                }
            }

            let error = self.refinement_error(idx);
            let lateness = self.lateness(idx, horizon);

            // Materialise when the matter cannot express what the node is
            // doing — and when the world's pace leaves room to actually
            // integrate it.
            //
            // Two independent reasons, and the second is the one that was
            // missing. `acuity` is somebody wanting to see it; `error` is the
            // node's own state saying the matter is no longer an adequate
            // description of it — an unresolved Jeans length, a dynamical time
            // shorter than the frame. A cloud collapsing in the dark refines
            // because it is collapsing, not because anyone turned to look.
            let wants_detail = (acuity > 0.5 && wanted_tier > tier) || error > REFINE_THRESHOLD;
            if !materialised && wants_detail && self.can_resolve(idx) {
                let n_children = self.tree.nodes[idx.get()].spec.count;
                tasks.push(Task {
                    node: idx,
                    kind: TaskKind::Materialise,
                    cost_us: cost::materialise_us(n_children),
                    lateness: acuity,
                    urgency: urgency.max(0.01),
                    error,
                    bytes: (n_children * std::mem::size_of::<Body>()) as i64,
                });
            }

            // Release detail the node has forgotten. Negative bytes: this task
            // gives resources back, so the planner always accepts it.
            //
            // The old rule was "coarsen when nobody is looking", and it is the
            // single line that made this an observer-driven engine: an
            // unobserved world discarded all its detail on the first frame and
            // then did no physics at all, because there was nothing left to
            // step. Detail now goes when something can rebuild it — see
            // [`World::collapsible`], which is a mixing time for matter and a
            // recipe for a structure — and not before, whoever is or is not
            // watching.
            if materialised && acuity < 0.25 && self.collapsible(idx) {
                tasks.push(Task {
                    node: idx,
                    kind: TaskKind::Coarsen,
                    cost_us: cost::COARSEN_US * count as f64,
                    lateness: 1.0,
                    urgency: 1.0,
                    error: 1.0,
                    bytes: -((count * std::mem::size_of::<Body>()) as i64),
                });
            }

            // Growth advances whether or not anything is materialised — in
            // fact especially when nothing is. This is the payoff of the
            // matter representation: a forest of 10^9 trees held as 10^4
            // nodes costs 10^4 ODE steps, so growth can run on the entire world
            // every frame while the fine structure stays unbuilt.
            if self.tree.nodes[idx.get()].morphology.is_some() {
                tasks.push(Task {
                    node: idx,
                    kind: TaskKind::Grow,
                    cost_us: cost::GROW_US,
                    lateness: 1.0,
                    urgency: 1.0,
                    error: 1.0,
                    bytes: 0,
                });
            }

            // Re-solve what is materialised, once it has come due.
            if materialised && lateness >= 1.0 {
                let kind = solvers::for_tier(tier);
                tasks.push(Task {
                    node: idx,
                    kind: TaskKind::Step,
                    // One solver pass. How many passes the node gets is
                    // decided at execution time, out of whatever the frame has
                    // left — see `World::execute`. Quoting the whole crossing
                    // here was worse than quoting one pass: it priced every
                    // node out of every frame, and a planner that accepts
                    // nothing leaves a world that is perfectly on time and
                    // completely still.
                    cost_us: cost::step_us(kind, count, self.gpu),
                    lateness,
                    urgency: urgency.max(0.05),
                    error,
                    bytes: 0,
                });
            }
        }
        tasks
    }

    /// How wrong is it to leave this node coarse?
    ///
    /// Two physical criteria, both of which are about structure the matter
    /// cannot represent: an unresolved Jeans length (the node is about to
    /// fragment) and a short dynamical time relative to the frame step (the
    /// node is evolving faster than we are looking at it).
    fn refinement_error(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        let a = &n.matter;
        let mut e = 0.05;
        let jeans = a.jeans_length();
        if jeans.is_finite() && jeans < 2.0 * a.radius {
            e += (2.0 * a.radius / jeans.max(1e-30)).min(1e3);
        }
        let dyn_t = a.dynamical_time();
        if dyn_t.is_finite() && dyn_t > 0.0 {
            e += (n.tier.dt() / dyn_t).min(1e3);
        }
        e
    }

    /// Run the plan, and report the worst fraction of the frame any solved node
    /// actually got across. One means everything that ran arrived.
    fn execute(&mut self, plan: &Plan, horizon: f64) -> f64 {
        // How many sub-steps each solved node may take.
        //
        // The plan decided *which* nodes run; this decides *how far* each one
        // gets, out of what the frame has left after the one pass each of them
        // was costed at. Shared equally rather than by lateness, because it has
        // to be computed without a clock: a scheduler that depended on how fast
        // the machine happened to be running that frame would not replay.
        let per_pass: f64 = plan
            .accepted
            .iter()
            .filter(|t| t.kind == TaskKind::Step)
            .map(|t| self.budget.adjust(t.cost_us))
            .sum();
        let allowance = if per_pass > 0.0 {
            ((self.budget.sim_budget_us() / per_pass).floor()).clamp(1.0, MAX_SUBSTEPS as f64) as u32
        } else {
            1
        };
        let mut achieved = 1.0f64;
        for task in &plan.accepted {
            if !self.tree.nodes[task.node.get()].alive {
                continue;
            }
            let started_at = self.tree.nodes[task.node.get()].time;
            match task.kind {
                TaskKind::Materialise => {
                    self.tree.refine(task.node);
                    self.disturb(task.node);
                }
                TaskKind::Coarsen => {
                    self.tree.coarsen(task.node);
                }
                TaskKind::Promote => {}
                TaskKind::Step => {
                    self.advance_to(task.node, horizon, allowance);
                    let asked = horizon - started_at;
                    if asked > 0.0 {
                        let got = self.tree.nodes[task.node.get()].time - started_at;
                        achieved = achieved.min((got / asked).clamp(0.0, 1.0));
                    }
                }
                TaskKind::Grow => {
                    let n = &self.tree.nodes[task.node.get()];
                    let dt = horizon - n.last_grown;
                    if dt > 0.0 {
                        // Growth is the archetypal thing a bubble is for — a
                        // century of a tree in an afternoon — so it runs on the
                        // node's own clock like the rest of its interior.
                        let rate = self.local_rate(task.node);
                        let physical = self.time_rate_of(task.node).physical();
                        self.stats.bubble_seconds += dt * (rate - physical).abs();
                        let local = dt * rate;
                        self.grow_node(task.node, local);
                        self.tree.nodes[task.node.get()].last_grown = horizon;
                    }
                }
                TaskKind::Observe => {}
            }
        }
        achieved
    }

    /// Bring a node toward the world instant by integrating it.
    ///
    /// Sub-steps at whatever its physics needs, which is a completely separate
    /// question from how often it is *scheduled*: a node may be visited once a
    /// frame and take four hundred steps to get across, or be visited once in a
    /// thousand frames and take one.
    ///
    /// Two ways it can fall short, and they mean different things:
    ///
    /// * **The span is unreachable in principle** — a resolved nucleus asked to
    ///   cover a millisecond would need 10^19 steps. Following the trajectory
    ///   is then not merely expensive, it is the wrong answer, and the node is
    ///   crossed by its ensemble instead. See [`World::thermalise`].
    /// * **The frame ran out of allowance.** The node integrates as far as it
    ///   can and stops, and `last_solved` records where it got to, so its
    ///   lateness says exactly how far behind it is. `coast_to` then carries
    ///   its frame the rest of the way in closed form: the node moves, its
    ///   insides do not, and nobody is told a story about either.
    pub fn advance_to(&mut self, idx: NodeIdx, horizon: f64, allowance: u32) {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return;
        }
        let span = horizon - self.tree.nodes[idx.get()].time;
        if !(span > 0.0) {
            return;
        }
        // `node_dt` is what the node's own physics needs, measured on the
        // node's own clock. The loop below walks *coordinate* time toward the
        // horizon, so the stable coordinate step is that divided by the rate:
        // a node running a hundred times faster inside a bubble covers a
        // hundredth as much world time per pass, and needs a hundred times as
        // many passes to cross the same frame. That is the honest price of a
        // bubble, and it is paid in the scheduler rather than hidden.
        let rate = self.local_rate(idx);
        let h0 = self.node_dt(idx) / rate;
        // What this node would need to be *followed*, before any cap. See
        // `Stats::worst_substeps`: reported rather than clamped, because the
        // engine dropping a node to its ensemble is a decision and should read
        // as one.
        if h0 > 0.0 && h0.is_finite() {
            let wanted = span / h0;
            if wanted > self.stats.worst_substeps {
                self.stats.worst_substeps = wanted;
                self.stats.worst_substeps_at = Some(self.tree.nodes[idx.get()].key);
            }
        }
        if h0 > 0.0 && h0.is_finite() && span / h0 > MAX_SUBSTEPS as f64 && self.forgettable(idx) {
            self.stats.ensembled += 1;
            self.thermalise(idx, horizon);
            return;
        }
        let started_at = self.tree.nodes[idx.get()].time;
        let mut steps = 0u32;
        while self.tree.nodes[idx.get()].time < horizon && steps < allowance {
            let remaining = horizon - self.tree.nodes[idx.get()].time;
            let h = (self.node_dt(idx) / rate).min(remaining);
            if !(h > 0.0) {
                break;
            }
            self.advance_node(idx, h);
            steps += 1;
        }
        // The gate above is a *prediction*, made from `node_dt` — the estimate
        // cheap enough to pay for every node every frame. A force field can
        // demand far less than that estimate once it has looked at the actual
        // configuration, and only `advance_node` pays the force evaluation that
        // finds out. So ask the question again on the way out, using what the
        // loop achieved rather than what it planned: the same test, against a
        // step that has been measured instead of guessed.
        //
        // This is what stops a node falling behind *forever*. Short of budget
        // and unreachable in principle look identical for one frame — the node
        // did not reach the horizon either way — and they are not the same
        // thing. A node that ran out of allowance catches up when the frame is
        // cheaper. A node whose own physics needs more than `MAX_SUBSTEPS`
        // passes per span never catches up, because the deficit is per frame
        // and reappears in the next one. That node is crossed by its ensemble,
        // which is exact for the only thing still observable at that cadence,
        // and is at the horizon afterwards rather than permanently behind it.
        if steps > 0 && self.tree.nodes[idx.get()].time < horizon {
            let achieved = (self.tree.nodes[idx.get()].time - started_at) / steps as f64;
            if !(achieved > 0.0) || span / achieved > MAX_SUBSTEPS as f64 {
                if self.forgettable(idx) {
                    self.stats.ensembled += 1;
                    self.thermalise(idx, horizon);
                    return;
                }
                // Not forgettable — pinned, edited, bubbled, or with something
                // built on it. Somebody is deliberately watching this node run, so the
                // one thing that must not happen is replacing the trajectory
                // they asked for with a draw from its equilibrium. It falls
                // behind, and goes on falling behind for as long as the world
                // runs at a pace its own physics cannot afford. That is a real
                // and permanent condition rather than a transient one, so it is
                // counted rather than left to be inferred from a lateness that
                // never comes down.
                self.stats.unreachable += 1;
            }
        }
        let n = &mut self.tree.nodes[idx.get()];
        let slack = horizon.abs() * 1e-12;
        if n.time + slack >= horizon {
            n.time = horizon;
        }
        n.last_solved = n.time;
    }

    /// How far a node's contents must be carried forward to be *drawn* at the
    /// world instant.
    ///
    /// A node the frame did not solve all the way has bodies that are behind
    /// where the node itself is. Drawing them there shows a world that stutters
    /// at exactly the rate the scheduler skips things, which is the one
    /// artefact the whole closed-form design exists to avoid. Carrying each
    /// body at its own velocity is the same exact solution applied one level
    /// down, and it costs one multiply-add per body at draw time.
    ///
    /// It is interpolation, not simulation: the forces are not re-evaluated, so
    /// this is right for as long as the node is not badly overdue, and its
    /// lateness is the number that says whether it is.
    pub fn render_lag(&self, idx: NodeIdx) -> f64 {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return 0.0;
        }
        let raw = (self.time - self.tree.nodes[idx.get()].last_solved).max(0.0);
        // Never carry further than one of the node's own characteristic times.
        //
        // A cadence is *defined* as how long the fastest body takes to cross one
        // resolution element, so this caps the extrapolation at one element —
        // which is exactly as far as a straight line is a refinement of the
        // truth rather than a different picture.
        //
        // The bound is not theoretical tidiness. `last_solved` can legitimately
        // be far in the past while the bodies are valid *now*: materialising a
        // node samples it from the matter as it currently is, and `refine`
        // has no clock to say so. Without this cap, rendering a freshly
        // materialised node in an old world flung every body to 10^8 node radii
        // — bodies that were sitting, correctly, at three and a half.
        let cadence = self.node_cadence(idx);
        if cadence.is_finite() && cadence > 0.0 {
            raw.min(cadence)
        } else {
            raw
        }
    }

    /// May this node's detail be thrown away and drawn again?
    ///
    /// No, if somebody has touched it — a tree you broke is not a
    /// representative sample of anything, and neither is one whose *recipe*
    /// carries a break ([`Node::contains_edit`]), because the draw would mend
    /// it. No, if something finer has been built on it, because releasing it
    /// would take the whole subtree with it.
    ///
    /// A node that fails this test may still be *collapsible* — see
    /// [`World::collapsible`], which is the question the scheduler asks before
    /// releasing detail, and which a broken box passes.
    pub fn forgettable(&self, idx: NodeIdx) -> bool {
        let n = &self.tree.nodes[idx.get()];
        // A bubble is somebody deliberately watching this node run. Crossing it
        // by ensemble is the right answer for matter nobody is following, and
        // exactly the wrong one here: the whole point of speeding a region up
        // is to see what it *does*, and thermalising it would replace that with
        // a fresh draw from its equilibrium. A bubbled node falls behind
        // honestly instead, and its lateness says by how much.
        !n.pinned
            && !n.contains_edit
            && n.bubble == 1.0
            && !n.children.iter().any(|c| !c.is_none())
    }

    /// May this node's detail be released, because something can rebuild it?
    ///
    /// **Not the same question as [`World::forgettable`], and `docs/PLAY.md`
    /// D19 is the decision that separates them.** Forgetting is a fresh draw
    /// from the ensemble, and it mends whatever happened to the node — so a
    /// broken thing must never be forgotten. Collapsing is throwing the detail
    /// away and running its *description* again, and a broken thing survives
    /// that unchanged, because the break is in the description.
    ///
    /// The scheduler used to gate its `Coarsen` task on mixing alone, which
    /// made the two the same question and cost both answers: a node carrying a
    /// morphology mixes in infinite time (correctly — no amount of waiting
    /// rearranges wood back into a sample), so **no grown or built thing ever
    /// coarsened, ever**. Measured on a thirty-year oak: 2 232 bodies resident
    /// for the life of the world, against 288 bytes of genome and event log
    /// that regenerate them exactly.
    ///
    /// Three refusals, and each is a different thing that cannot be rebuilt:
    ///
    /// * **Pinned** — this node's *own* detail was altered body by body, and
    ///   only the store has it.
    /// * **Promoted children** — releasing the subtree would take specific
    ///   objects with their own histories, not a sample of anything.
    /// * **A topology with no morphology to regenerate it** — structure that
    ///   nothing can draw again.
    ///
    /// Otherwise: a morphology regenerates its bodies from `(matter, genome,
    /// age, events)` and may collapse whenever nothing is watching, and
    /// everything else may collapse once it has had a mixing time to forget
    /// what put it there.
    pub fn collapsible(&self, idx: NodeIdx) -> bool {
        let n = &self.tree.nodes[idx.get()];
        if n.pinned || n.children.iter().any(|c| !c.is_none()) {
            return false;
        }
        if n.morphology.is_some() {
            return true;
        }
        if n.topology.is_some() {
            return false;
        }
        self.time - n.last_disturbed > self.mixing_time(idx)
    }

    /// Cross a span too long to integrate, by ensemble instead of trajectory.
    ///
    /// This is the step that makes one shared instant affordable across
    /// thirty-eight orders of magnitude. A resolved nucleus asked to cross
    /// fifty milliseconds would need 10^21 steps; it also does not need them,
    /// because over that span it has sampled its accessible states 10^21 times
    /// and where it ends up is a draw from its equilibrium ensemble, not the
    /// endpoint of a trajectory. So the detail is summarised back to matter,
    /// that matter is carried across in closed form, and the detail
    /// is drawn again at the far end.
    ///
    /// Both halves are things the engine already guarantees: summarising is
    /// conservative to within `IDEMPOTENT_TOLERANCE` and sampling is a
    /// maximum-entropy sample of the same conserved tuple, which is exactly
    /// what "a fresh draw from the ensemble" means.
    ///
    /// The node is left coarse rather than immediately re-drawn. Re-drawing it
    /// here would hide the cost from the frame budget and, worse, would do it
    /// again next frame and every frame after: if the world is running at a
    /// pace this node cannot be integrated at, the honest thing is to stop
    /// pretending to resolve it. The survey will materialise it again when
    /// something wants it *and* the pace allows it — see the gate in
    /// [`World::can_resolve`].
    ///
    /// Detail that has been *touched*, or that something finer has been built
    /// on, is exempt. A tree somebody broke is not a representative sample of
    /// anything, and re-drawing it would silently mend it; releasing a node
    /// with promoted children would take the whole subtree under it. Those
    /// nodes fall behind honestly instead — and their lateness says so.
    fn thermalise(&mut self, idx: NodeIdx, horizon: f64) {
        if !self.forgettable(idx) {
            return;
        }
        let was_materialised = self.tree.nodes[idx.get()].is_materialised();
        if was_materialised {
            self.tree.coarsen(idx);
        }
        self.tree.carry(idx, horizon);
        let n = &mut self.tree.nodes[idx.get()];
        n.time = horizon;
        n.last_solved = horizon;
        n.epoch = n.epoch.wrapping_add(1);
        self.stats.thermalised += 1;
    }

    /// Carry every node that was not solved to the world instant.
    ///
    /// Free, and exact. A node's offset under constant velocity and its
    /// orientation under constant spin are both closed-form solutions, so there
    /// is no approximation here at all — only the assumption that nothing
    /// changed the velocity or the spin, which is precisely what "it was not
    /// due to be re-solved" asserts.
    ///
    /// Returns how many nodes were carried rather than solved. It is nearly all
    /// of them, nearly every frame, and that is the point.
    ///
    /// **Where a node is, not what it holds.** Every node's motion is carried to
    /// the world instant here. Its contents' clock follows only where its
    /// contents are closed-form too, which is a node holding only matter; a
    /// node with bodies is never left unsolved — the owner's rule for Phase 5 —
    /// so it keeps the time its bodies were solved to and is scheduled, and
    /// when it runs it covers all of it. See `Node::carried`.
    fn coast_to(&mut self, horizon: f64) -> usize {
        let mut coasted = 0;
        for i in 0..self.tree.nodes.len() {
            if !self.tree.nodes[i].alive {
                continue;
            }
            let idx = NodeIdx(i as u32);
            self.tree.carry(idx, horizon);
            let owed = horizon - self.tree.nodes[i].time;
            // **What it holds catches up wherever it is allowed to.** The
            // owner's rule for Phase 5 measures causality on what a node holds,
            // and a node whose contents are further behind than its own physics
            // could ever integrate is the gate `advance_to` applies — except
            // that one only runs for a node the budget chose. A node it never
            // chose went on owing: measured before this, three unpinned
            // children of a galaxy paced at 1e15 s a frame against their
            // 3.2e12 s step, on a budget that accepts one task, left two of
            // them at t = 0. So the same test, against the same step, for
            // everything the frame did not solve. A node that may not be
            // redrawn — pinned, edited, bubbled, holding children — stays
            // behind and is what `check_causality` reports; so is one that
            // could integrate its span and was simply not chosen, which is
            // scheduled and ranked by its lateness.
            if owed > 0.0 && !self.tree.nodes[i].bodies.is_empty() {
                let h0 = self.node_dt(idx) / self.local_rate(idx);
                if h0 > 0.0 && h0.is_finite() && owed / h0 > MAX_SUBSTEPS as f64 && self.forgettable(idx) {
                    self.stats.ensembled += 1;
                    self.thermalise(idx, horizon);
                }
            }
            let n = &mut self.tree.nodes[i];
            let dt = horizon - n.time;
            if !(dt > 0.0) || !n.bodies.is_empty() {
                continue;
            }
            n.time = horizon;
            coasted += 1;
            let key = n.key;
            // Same two clocks as `advance_node`: the frame took its own
            // kinematic share inside `advance`, and the node's own clock takes
            // the whole chain.
            let physical = self.time_rate_of(NodeIdx(i as u32)).physical();
            if let Some(c) = self.clocks.get_mut(&key) {
                c.time = horizon;
                c.proper_time += dt * physical;
            }
        }
        coasted
    }

    /// The sub-step a node's own physics needs, in seconds.
    ///
    /// A stability limit and nothing else. It used to double as the node's
    /// scheduling interval, which conflated two questions that have different
    /// answers by many orders of magnitude: how *finely* a node must be
    /// integrated when it is integrated, and how *often* it is worth
    /// integrating at all. The second is now [`World::node_cadence`].
    ///
    /// Every characteristic time the node has, not just the gravitational one.
    /// A node's parcels have to be stepped faster than a sound wave crosses the
    /// gap between them, or the pressure force acts across a distance the
    /// information could not have travelled — which is not a small error, it is
    /// a solver that heats its own contents. It showed up as a continuum node
    /// three metres across whose gas reached two thirds of light speed while
    /// the conservation check reported no drift at all, because every
    /// individual step was conserving the energy the previous one had invented.
    pub fn node_dt(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        let parts = if n.bodies.is_empty() { n.spec.count } else { n.bodies.len() };
        let flow = n.bodies.iter().map(|b| b.vel.norm()).fold(0.0f64, f64::max);
        // The signal crosses at the speed the node's *own* equation of state
        // gives it, not the gas law's: see `World::sound_speed_of`.
        let signal = self.signal_speed_of(idx).max(flow);
        let spacing = n.matter.radius / (parts.max(1) as f64).cbrt();
        let crossing = if signal > 0.0 && signal.is_finite() { spacing / signal } else { f64::INFINITY };
        let mut natural = n
            .tier
            .dt()
            .min(n.matter.dynamical_time() / 50.0)
            .min(0.25 * crossing);
        // Where a force field decides the timestep, ask the force field. The
        // engine already has a function whose entire job is "what step does
        // this system need"; the scheduler was not calling it.
        if matches!(
            solvers::for_tier(n.tier),
            SolverKind::MolecularDynamics
        ) && !n.bodies.is_empty()
        {
            // The cheap bound here; `advance_node` substeps to the
            // configuration-aware one, which costs a force evaluation and must
            // not be paid on every scheduling decision.
            natural = natural.min(solvers::md::stable_dt(&n.bodies) * 8.0);
        }
        natural
    }

    /// Run the tier's solver over a node's materialised bodies.
    pub fn advance_node(&mut self, idx: NodeIdx, dt: f64) -> solvers::SolveReport {
        let (tier, key, epoch, radius, count, tick) = {
            let n = &self.tree.nodes[idx.get()];
            (n.tier, n.key, n.epoch, n.matter.radius, n.bodies.len(), n.steps_taken)
        };
        // `!(dt > 0.0)` rather than `dt <= 0.0` so a NaN span is refused rather
        // than passed through: NaN fails every comparison, so the old spelling
        // let it past, and a node's clock never recovers from one.
        if count == 0 || !(dt > 0.0) || !dt.is_finite() {
            return solvers::SolveReport::default();
        }
        // `dt` is coordinate time — the span the world clock moved. What the
        // node's *interior* experiences is that span on the node's own clock,
        // which is where relativity and any bubble enter. Everything below the
        // solver call therefore runs on `local`; everything about where the
        // node *is* stays on `dt`. See `dilation.rs` for why the split falls
        // this way round.
        let rate = self.local_rate(idx);
        let dt = dt * rate;
        let seed = self.tree.world_seed;

        // The child is the real thing and its body is a stand-in, so the solver
        // has to see where the child actually is before it computes anything.
        // See `Tree::sync_children`.
        let promoted = self.tree.sync_children(idx);
        let before = self.tree.stand_in_velocities(idx, &promoted);

        // What is ordered is not the tier solver's business. `docs/PLAY.md`
        // §3.3: **the tier says which regime the disordered contents are in,
        // and the node's own state says which contents are ordered.**
        //
        // A node holds both. `sample_structured` lays a structure's members out
        // first and then "the unstructured remainder: litter, air, rubble" —
        // one node, two kinds of thing — and until this, `for_tier` was handed
        // the lot. A building, a wolf and a boulder are all `Continuum`, and
        // `for_tier(Continuum)` is `Hydro`.
        //
        // Measured, on a forty-year-old tree standing on a planet, advanced for
        // one twentieth of a second: its members reached 1.9x10^8 m/s — 64% of
        // the speed of light — and travelled 9.6x10^6 m. The tree is 6.3 m
        // across. SPH reads `Matter` through a gas equation of state, so a
        // solid handed to it bursts from its own pressure before anything else
        // happens; and the same shape of error is waiting in every other tier's
        // solver, which would integrate a joined member as a free particle.
        //
        // The members are left where they are, which is what a standing
        // structure does. What *moves* them is the structural path — `damage`
        // asks whether it stands up under a load, `shake` asks what it does
        // while the load is on it — and those are driven by a caller with
        // mechanisms in hand rather than by the scheduler. Putting them on the
        // frame loop is a scheduling question and is not this.
        let ordered = self.tree.nodes[idx.get()].structural_mask();

        // Read before the bodies are borrowed, because the report below cannot
        // reach the node once they are. See the `SolverKind::Hydro` arm.
        let gas_law = self.tree.nodes[idx.get()].matter.gas_law_applies();
        let stand_ins = self.tree.nodes[idx.get()]
            .children
            .iter()
            .filter(|c| !c.is_none())
            .count();
        let mut eos_suspect = false;
        // What each body answers to when squeezed, in the order of the node's
        // bodies: a stand-in for a promoted child answers with the child's own
        // matter, and everything else with the node's. `eos.rs`.
        //
        // This is what the stand-in detonation in `BACKLOG.md` wanted, seen
        // from the parent's side: a 48-tonne box and a 10 kg ball in a room
        // were a hot dense gas to their parent's SPH and gained a factor of
        // three a frame with nothing touching. Priced as what they are — two
        // solids, each far below its own rest density at the parent's
        // smoothing length — they carry no pressure at all.
        let body_eos: Vec<crate::eos::Eos> = if matches!(solvers::for_tier(tier), SolverKind::Hydro) {
            let n = &self.tree.nodes[idx.get()];
            // A structure is its solids, so what is loose around it is the
            // rest of the mixture.
            let own = if ordered.is_some() {
                crate::eos::Eos::of_loose(&n.matter.mixture, &self.substances)
            } else {
                crate::eos::Eos::of_matter(&n.matter, &self.substances)
            };
            (0..n.bodies.len())
                .map(|i| match n.children.get(i) {
                    Some(c) if !c.is_none() => {
                        crate::eos::Eos::of_matter(&self.tree.nodes[c.get()].matter, &self.substances)
                    }
                    _ => own,
                })
                .collect()
        } else {
            Vec::new()
        };

        // What the loose contents can rest on, and whether they rest at all.
        // The node's own ordered members are walls to them (`hydro::Wall`),
        // and **contents with something inside their own node to rest on carry
        // their own weight** — the node's stored field. A fluid region with
        // nothing ordered in it is held up by the fluid around it, which is
        // outside the node, and gets no uniform field: that keeps every gas
        // node the engine already had exactly as it was.
        //
        // A member is a slab where the generator stated one and a beam from its
        // base to its tip otherwise — a body alone is a sphere or a box, and a
        // branch is neither.
        let (mut walls, field) = match &ordered {
            Some(mask) => {
                let n = &self.tree.nodes[idx.get()];
                let topo = n.topology.as_ref();
                let walls: Vec<solvers::hydro::Wall> = n
                    .bodies
                    .iter()
                    .enumerate()
                    .zip(mask.iter())
                    .filter(|(_, o)| **o)
                    .map(|((i, b), _)| {
                        let beam = topo.and_then(|t| Some((*t.base.get(i)?, *t.tip.get(i)?)));
                        match beam {
                            Some((base, tip)) if !b.is_boxed() && (tip - base).norm() > 0.0 => {
                                solvers::hydro::Wall::capsule(base, tip, b.radius)
                            }
                            _ => solvers::hydro::Wall::of(b),
                        }
                    })
                    .collect();
                (walls, n.gravity)
            }
            None => (Vec::new(), crate::math::Vec3::ZERO),
        };
        // Whether the node holds a liquid, and which of its bodies stand in for
        // promoted children — those answer with their own law and are not the
        // node's liquid.
        let liquid_node = self.tree.nodes[idx.get()].matter.mixture.in_phase(crate::chem::Phase::Liquid) > 0.0;
        let stand_in: Vec<bool> =
            self.tree.nodes[idx.get()].children.iter().map(|c| !c.is_none()).collect();
        // **A packed solid is not a fluid, whether or not anything built it.**
        // `docs/PLAY.md` §3.5: solids escape the resolution floor, and it binds
        // free fluid only. A structure's members were already kept from the
        // tier solver by §3.3's dispatch; a sampled lump of solid — a rock, a
        // wooden ball — had no topology to say so and went through SPH at its
        // own elastic sound speed, which since it has an equation of state is
        // the real one: 11 km/s for wood, cutting a frame into 23,000 solves of
        // four bodies that did not move. Its motion is the node's, and inside
        // it there is nothing for the tier solver to do. Anything with a
        // promoted child in it is left to the solver, which is how the child
        // feels its parent.
        let rigid = {
            let n = &self.tree.nodes[idx.get()];
            !liquid_node
                && ordered.is_none()
                && !stand_in.iter().any(|s| *s)
                && matches!(solvers::for_tier(tier), SolverKind::Hydro)
                && crate::eos::Eos::of_matter(&n.matter, &self.substances).condensed().is_some()
                && crate::sampler::is_packed(&n.matter, n.rest_density)
        };

        let bodies = &mut self.tree.nodes[idx.get()].bodies;
        // Solve the disordered contents in place where there are no ordered
        // ones, which is every node that is not a structure and costs nothing.
        let mut loose: Vec<crate::state::Body> = Vec::new();
        let mut loose_of: Vec<usize> = Vec::new();
        if let Some(mask) = &ordered {
            for (i, b) in bodies.iter().enumerate() {
                if !mask.get(i).copied().unwrap_or(false) {
                    loose.push(*b);
                    loose_of.push(i);
                }
            }
        }
        let partitioned = ordered.is_some();
        let bodies: &mut Vec<crate::state::Body> = if partitioned {
            &mut loose
        } else {
            &mut self.tree.nodes[idx.get()].bodies
        };
        let count = bodies.len();
        // A structure with no loose contents at all — a solid object, whose
        // every body is a member. There is nothing for the tier solver to do,
        // and its members are not its to move.
        //
        // This used to `return` here, and that was wrong in a way nothing
        // caught, because everything below runs on the *node* rather than on
        // its contents: the node's own clock, its motion, the contacts it is in
        // and the heat crossing its boundaries. Skipping the solver skipped
        // those too, so a node made entirely of ordered matter never moved,
        // never aged, never collided and never exchanged — a solid rock adrift
        // in space, frozen, while a rock with a single speck of dust in it
        // behaved perfectly.
        //
        // Measured, on a 48-tonne box given 1 m/s for one second of frames:
        // it moved 0.000000 m and its clock read 0.0 s after a hundred steps.
        // The same box with one loose body added moved 1.000000 m. A structure
        // whose members are not the solver's to move is not the same statement
        // as a node that is not there.
        let report = if count == 0 {
            solvers::SolveReport::default()
        } else if rigid {
            // Nothing inside it to integrate, and the whole span covered.
            solvers::SolveReport { dt_used: dt, ..Default::default() }
        } else {
            match solvers::for_tier(tier) {
            SolverKind::Gravity | SolverKind::GravityHydro => {
                let params = solvers::gravity::GravityParams {
                    theta: 0.5,
                    softening: radius / (count as f64).cbrt() * 0.3,
                    retarded: true,
                    post_newtonian: tier == Tier::Planetary,
                    // The quadrupole is worth its cache line only where an
                    // observer could resolve the difference.
                    quadrupole: tier >= Tier::Planetary,
                };
                solvers::gravity::step_leapfrog(bodies, dt, params)
            }
            SolverKind::Hydro => {
                // **The equation of state, applied outside its validity, says
                // so.** `PLAY.md` §7's ninth Phase 2 item, on §3.7's precedent
                // that a node crossed by its ensemble reports it rather than
                // doing it quietly.
                //
                // SPH prices a node's `Matter` through `pressure`, which is an
                // ideal gas plus radiation and the only equation of state the
                // engine has. There are two ways to be outside it and both are
                // now measurable:
                //
                // - **The matter is not a gas.** D17 put a mixture on every
                //   `Matter`, so a node that has been described knows its own
                //   phase. A bucket of water prices at 4x10^8 Pa.
                // - **The node is mostly vacuum with solids in it.** A
                //   `Continuum` node whose every body is the stand-in for a
                //   promoted child has no contents of its own, and pricing that
                //   as a hot dense gas detonates it: measured, on a 12 m root
                //   holding a 48-tonne box and a 10 kg ball with **zero
                //   collisions**, the ball's speed goes 5.36 -> 7.96 -> 20.9 ->
                //   63.3 -> 194 -> 566 m/s, roughly threefold a frame, against
                //   a control of exactly 5.3572 every frame with the root not
                //   advanced.
                //
                // Phase 5 gave both an equation of state (`eos.rs`). What is
                // left to report is matter the engine still cannot
                // price: a stand-in or a condensed node whose substance has no
                // condensed equation of state, which an undescribed child is.
                let mut eos: Vec<crate::eos::Eos> = if partitioned {
                    loose_of.iter().map(|&i| body_eos[i]).collect()
                } else {
                    body_eos.clone()
                };
                let priced = |e: &crate::eos::Eos| matches!(e, crate::eos::Eos::Condensed(_));
                eos_suspect = (!gas_law && !eos.iter().any(priced))
                    || (count > 0 && stand_ins >= count && !eos.iter().all(priced));

                // **Weakly-compressible SPH**, `docs/PLAY.md` §4.2 and §3.7's
                // second response to the resolution floor: a liquid's physical
                // sound speed is replaced by one ten times the fastest thing in
                // it — its fastest parcel, or the `sqrt(2 g H)` a column of
                // depth `H` reaches when it falls — which keeps its density
                // within a per cent of rest (Mach 0.1 squared) and is the
                // standard treatment for a free surface. For a 2 m/s wave it
                // is 20 m/s instead of 1500, and §4.2's scene goes from 30,000
                // substeps to 80. It is an approximation with a known error
                // bound: the liquid is a per cent more compressible than it
                // is, and nothing else changes.
                //
                // Only for a liquid, and only for the node's own contents. A
                // solid is not a free-surface flow and keeps the speed it has;
                // a stand-in answers with its child's own law.
                //
                // And the smoothing length is the liquid's own spacing, not the
                // node's radius over its count: a liquid parcel sits where its
                // rest density puts it (`sampler::packed_positions`), and the
                // kernel has to span its neighbours rather than the node.
                // **A solid child standing in the fluid is a wall to it**, and
                // takes back every push it gives: see `hydro::Wall`. Only where
                // the node holds a liquid — a stand-in in a gas keeps the
                // pairwise coupling it has always had.
                if liquid_node {
                    for (k, e) in eos.iter().enumerate() {
                        let slot = if partitioned { loose_of[k] } else { k };
                        if stand_in.get(slot).copied().unwrap_or(false) && priced(e) {
                            let mut w = solvers::hydro::Wall::of(&bodies[k]);
                            w.owner = Some(k as u32);
                            walls.push(w);
                        }
                    }
                }
                let spacing = {
                    let (mut m, mut k) = (0.0, 0usize);
                    let mut rest = 0.0;
                    for (q, (b, e)) in bodies.iter().zip(eos.iter()).enumerate() {
                        let slot = if partitioned { loose_of[q] } else { q };
                        if stand_in.get(slot).copied().unwrap_or(false) {
                            continue;
                        }
                        if let crate::eos::Eos::Condensed(c) = e {
                            m += b.mass;
                            rest += c.rest_density * b.mass;
                            k += 1;
                        }
                    }
                    (k > 0 && m > 0.0).then(|| (m / k as f64 / (rest / m)).cbrt())
                };
                if liquid_node {
                    let drivers = walls
                        .iter()
                        .filter_map(|w| w.owner)
                        .map(|o| bodies[o as usize].vel.norm())
                        .fold(0.0f64, f64::max);
                    let scale = flow_scale(
                        bodies
                            .iter()
                            .enumerate()
                            .filter(|(k, _)| !walls.iter().any(|w| w.owner == Some(*k as u32)))
                            .zip(eos.iter())
                            .filter(|(_, e)| priced(e))
                            .map(|((_, b), _)| b),
                        field,
                        spacing.unwrap_or(0.0),
                        drivers,
                    );
                    for (k, e) in eos.iter_mut().enumerate() {
                        let slot = if partitioned { loose_of[k] } else { k };
                        if stand_in.get(slot).copied().unwrap_or(false) {
                            continue;
                        }
                        if let crate::eos::Eos::Condensed(c) = *e {
                            let artificial = (10.0 * scale).min(c.sound_speed(c.rest_density));
                            *e = crate::eos::Eos::Condensed(crate::eos::Condensed {
                                bulk_modulus: c.rest_density * artificial * artificial,
                                stiffening: crate::eos::LIQUID_STIFFENING,
                                ..c
                            });
                        }
                    }
                }
                let params = solvers::hydro::HydroParams {
                    h: match spacing {
                        Some(s) => 1.3 * s,
                        None => radius / (count as f64).cbrt() * 1.2,
                    },
                    gravity: field,
                    ..Default::default()
                };
                // Substep to what the Courant condition allows, for the same
                // reason molecular dynamics does below and with the same words:
                // the scheduler's timestep answers to causality and to the
                // tier; a pressure wave answers to neither, and a fluid handed
                // a step longer than its own signal crossing time does not
                // integrate inaccurately, it detonates.
                //
                // `hydro::courant_dt` has been here the whole time and nothing
                // called it. Measured, on the loose contents of a tree standing
                // on a planet: the stable step is 1.7x10^-5 s and the frame
                // asked for 0.05 — **three thousand times over** — and the
                // parcels left at 4x10^7 m/s.
                //
                // The cap is a budget and not a licence, again as below: a node
                // that cannot afford the whole span covers the part it can
                // integrate stably, reports it in `dt_used`, and lets the
                // shortfall become lateness the scheduler can see.
                let stable = solvers::hydro::courant_dt_with(bodies, params, 0.3, &eos, &walls).max(1e-30);
                let wanted = (dt / stable).ceil();
                let substeps = (wanted.clamp(1.0, MAX_SUBSTEPS as f64) as u32).max(1);
                let h = (dt / substeps as f64).min(stable);
                let mut total = solvers::SolveReport::default();
                for k in 0..substeps {
                    let r = solvers::hydro::step_with(bodies, h, params, &eos, &walls);
                    if k == 0 {
                        total = r;
                    } else {
                        total.after = r.after;
                        total.steps += r.steps;
                        total.interactions += r.interactions;
                        total.non_mechanical_energy += r.non_mechanical_energy;
                        total.unrest = r.unrest;
                    }
                }
                total.dt_used = h * substeps as f64;
                total
            }
            SolverKind::MolecularDynamics => {
                let params = solvers::md::MdParams::default();
                // Substep to whatever the force field needs. The scheduler's
                // timestep answers to causality and to the tier; the Lennard-
                // Jones potential answers to neither, and a molecular system
                // handed a step longer than its own vibrational period does not
                // integrate inaccurately, it detonates.
                let stable = solvers::md::configuration_dt(bodies, params).max(1e-24);
                // The cap on that substepping is a budget, not a licence. When
                // the span needs more substeps than one pass will pay for, the
                // node covers the part it can integrate *stably* and stops
                // there; `dt_used` reports how far it got and the node's clock
                // follows it, so the shortfall becomes lateness the scheduler
                // can see. Stretching the step to fit the cap instead — which
                // is what `dt / substeps` did unconditionally — is precisely
                // the detonation described above, delivered by the guard that
                // exists to prevent it.
                let wanted = (dt / stable).ceil();
                // `.max(1)` after the cast, not just the clamp before it: a
                // float clamp propagates NaN, and `NaN as u32` is zero.
                let substeps = (wanted.clamp(1.0, MD_MAX_SUBSTEPS as f64) as u32).max(1);
                let h = (dt / substeps as f64).min(stable);
                let mut total = solvers::SolveReport::default();
                for k in 0..substeps {
                    let r = solvers::md::step(bodies, h, params, seed, key.0, epoch, tick + k as u64);
                    if k == 0 {
                        total = r;
                    } else {
                        total.after = r.after;
                        total.steps += r.steps;
                        total.interactions += r.interactions;
                        total.non_mechanical_energy += r.non_mechanical_energy;
                    }
                }
                total.dt_used = h * substeps as f64;
                total
            }
            SolverKind::Statistical => self.advance_statistical(idx, dt),
            }
        };

        // What the solve left out of balance, for contents held in a field —
        // the only case where it is not zero at rest. See `Node::unrest`.
        self.tree.nodes[idx.get()].unrest =
            if field != crate::math::Vec3::ZERO && report.unrest.is_finite() { report.unrest } else { 0.0 };

        // Reported once the solver has let go of the bodies. See the
        // `SolverKind::Hydro` arm for what this means and why it is reported
        // rather than corrected.
        if eos_suspect {
            self.stats.eos_outside_validity += 1;
            self.stats.eos_outside_validity_at = Some(self.tree.nodes[idx.get()].key);
        }

        // The disordered contents were solved in a buffer; put them back.
        if partitioned {
            let n = &mut self.tree.nodes[idx.get()];
            for (k, &i) in loose_of.iter().enumerate() {
                if let (Some(dst), Some(src)) = (n.bodies.get_mut(i), loose.get(k)) {
                    *dst = *src;
                }
            }
        }

        // ... and the force it computed on each stand-in is handed to the
        // child it stands for. Without this the force lands on the body and is
        // discarded by the next sync, which is why two promoted things could
        // not affect each other at all.
        if !before.is_empty() {
            // The pull each child received, before it is handed on: what the
            // fluid around it pushes back against. See `buoy_children`.
            let pulls: Vec<(NodeIdx, Vec3)> = {
                let n = &self.tree.nodes[idx.get()];
                before
                    .iter()
                    .filter_map(|(slot, was)| {
                        let c = n.children.get(*slot).copied()?;
                        let now = n.bodies.get(*slot)?.vel;
                        (!c.is_none()).then_some((c, now - *was))
                    })
                    .collect()
            };
            self.tree.apply_body_forces(idx, &before);
            self.buoy_children(idx, &pulls);
        }

        self.stats.bodies_stepped += count as u64;
        // A solver that could not cover the whole span says so in `dt_used`.
        // The node's clock has to agree with it: advancing by the span that was
        // *asked for* rather than the one that was *integrated* would make the
        // node's lateness a fiction, and lateness is the one number the
        // scheduler uses to decide what to do about a node that is struggling.
        let dt = if report.dt_used.is_finite() && report.dt_used > 0.0 {
            report.dt_used.min(dt)
        } else {
            dt
        };
        // What the node holds may no longer be what the node says it is. Cheap
        // here and nowhere else: this pass has just touched every one of those
        // bodies, so the measurement rides a cache line that is already warm,
        // and it is on the node's own cadence — which is what
        // `BACKLOG.md` asks for, since a node nobody is advancing is not
        // spreading either.
        let spread = self.tree.spread(idx);
        let occupancy = spread.occupancy(self.tree.nodes[idx.get()].matter.radius);
        if occupancy > self.stats.worst_occupancy {
            self.stats.worst_occupancy = occupancy;
            self.stats.worst_occupancy_at = Some(key);
        }

        // Whatever the solver has just driven into whatever else. Before the
        // exchange, so the heat a collision makes is there to be conducted
        // away in the same pass rather than a frame later.
        self.contact_within(idx);

        // Heat crosses the boundaries between the things this node holds, on
        // the same span the solver just integrated. After the solve rather than
        // before it: the solver moves them, and what is next to what is a
        // question about where they ended up.
        self.exchange_within(idx, dt);

        // Back to coordinate time for everything the *parent* observes.
        let coordinate = dt / rate;
        let physical_rate = self.time_rate_of(idx).physical();
        self.stats.bubble_seconds += coordinate * (rate - physical_rate).abs();
        let n = &mut self.tree.nodes[idx.get()];
        n.time += coordinate;
        n.steps_taken += 1;
        // The solve may have moved the angular momentum, the mass or the
        // radius, and the angular velocity is derived from all three. Re-derive
        // before carrying the frame forward, so the orientation this span
        // advances by is the one the node's contents actually imply.
        n.sync_spin_rate();
        // Its motion is carried to the world instant by `Tree::carry`, on its
        // own clock; this advances what it holds. See `Node::carried`.
        let node_time = n.time;
        let clock = self
            .clocks
            .entry(key)
            .or_insert_with(|| Clock::new(node_time, tier.dt()));
        clock.time = node_time;
        // Two proper times, deliberately. `Motion::proper_time` is the frame's
        // own kinematic share against its immediate parent, which is what makes
        // a `Motion` meaningful on its own. `Clock::proper_time` is what a clock
        // actually sitting on the node reads: the whole chain, gravity
        // included. A bubble is excluded from both, because proper time is a
        // physical reading and an administrator speeding a region up does not
        // change what its clocks say.
        let physical = self.time_rate_of(idx).physical();
        if let Some(c) = self.clocks.get_mut(&key) {
            c.proper_time += coordinate * physical;
        }
        report
    }

    /// The material a node presents to something that runs into it, if it
    /// presents one at all.
    ///
    /// **Measured first, and only then inherited.** `docs/PLAY.md` D13: a
    /// node's material follows from what it *is*, not from what generated it.
    /// §2A measured what the old rule cost — this read `topology.material` or
    /// `morphology.material()` and returned `None` for everything else, so a
    /// rock, a boulder and a ball of wood had no surface and could not collide.
    /// That is a *provenance test standing in for a state measurement*, and
    /// §3.3 had already made the same call correctly one layer down:
    /// `structural_mask` asks what a node's joints measure, not what generated
    /// it.
    ///
    /// So the order is:
    ///
    /// 1. **What the matter is made of**, through `Material::measured` over the
    ///    node's own mixture, which D17 put on `Matter` for exactly this. Only
    ///    the *solid* pools count, and the formation conditions come from the
    ///    node's own radiative balance — see `Formation::of_matter`. A node
    ///    that has been described answers from its description, whatever made
    ///    it.
    /// 2. **What a structure was built out of**, for a node whose chemistry
    ///    nobody has stated. `Morphology::material` is itself a per-species
    ///    table and one of the columns D11 exists to move off `Program`; it
    ///    stays as the fallback until every generator seeds a mixture.
    /// 3. **Nothing.** A gas parcel or a star cluster has no surface, and that
    ///    is now measured rather than assumed: an undescribed node, or one
    ///    whose mixture holds no solid pool at all, gets `None`.
    fn surface_of(&self, idx: NodeIdx) -> Option<crate::neighbourhood::Resilience> {
        Some(crate::neighbourhood::Resilience::of(&self.contact_material(idx)?))
    }

    /// The material a contact with this node reads, **without measuring it
    /// again**.
    ///
    /// The baked surface already carries the material of every piece, and
    /// measuring one is not cheap: `grain_scale` marches a thousand steps of a
    /// nucleation integral, and `Material::measured` blends over every solid
    /// pool. Doing that per contact per frame took the steady-state frame from
    /// 25 ms to 121 ms — measured, on `frames_stay_within_budget`.
    ///
    /// So the surface is the cache, which is what `PLAY.md` D13 means by
    /// "derived once and stored": a caller that has baked a surface reads the
    /// answer off it, and one that has not falls back to what built the node.
    fn contact_material(&self, idx: NodeIdx) -> Option<crate::material::Material> {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return None;
        }
        let n = &self.tree.nodes[idx.get()];
        if let Some(s) = &n.surface {
            if n.surface_epoch == n.epoch {
                return s.pieces().first().map(|p| p.material);
            }
        }
        match &n.topology {
            Some(t) => Some(t.material),
            None => n.morphology.as_ref().map(|m| m.material()),
        }
    }

    /// The boundary this node presents, baked if it is stale. `PLAY.md` D18.
    ///
    /// **The generator emits the pieces; nothing infers a decomposition.** That
    /// is the clause that keeps convex decomposition — normally the hard,
    /// unsolved half of this problem — from arising at all. There are two
    /// generators and no third:
    ///
    /// - **A structure states its members.** `Topology` carries `base`, `tip`
    ///   and a cross-section radius per member, which *are* a capsule, and the
    ///   recipe is what put them there. A wall with a doorway emits the pieces
    ///   around the opening because the recipe put the opening there.
    /// - **Unstructured solid matter states one piece.** A rock is one filled
    ///   solid with no cavity in it, so it presents one hull over its bodies —
    ///   or, unmaterialised, the sphere of its own radius. Emitting one piece
    ///   per body would be emitting four thousand of them, and
    ///   `PERFORMANCE.md`'s narrow-phase row prices that at 1.2 ms for a single
    ///   contact.
    ///
    /// **A node that is not solid gets an empty surface**, which is D13's own
    /// line rather than a fallback: a liquid's surface is a property of its
    /// container and the field it is in, and a cloud of gas has none for the
    /// same reason it has no shape.
    ///
    /// The result is stored on the node and kept until its `epoch` moves. See
    /// `Node::surface`.
    pub fn surface_of_node(&mut self, idx: NodeIdx) -> &crate::shape::Surface {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return Self::nothing();
        }
        let stale = {
            let n = &self.tree.nodes[idx.get()];
            n.surface.is_none() || n.surface_epoch != n.epoch
        };
        if stale {
            self.stats.surfaces_baked += 1;
            let baked = self.bake_surface(idx);
            // **D18's invariant, checked at bake time.** A node's surface
            // materials are a partition of its solid pools, and nothing in D17
            // or D18 alone stops the two drifting: a node whose mixture is all
            // water by mass while its primitives present steel would be
            // *telling* a contact something its own bulk contradicts, which is
            // the second axiom failing by way of having two answers to one
            // question.
            //
            // It is cheap here and unpleasant to retrofit once both
            // representations exist and have drifted, which is why D18 asks for
            // it asserted rather than documented. Reported rather than
            // panicking, on §3.7's precedent: a node that cannot describe
            // itself says so.
            let n = &self.tree.nodes[idx.get()];
            if let Err(why) = baked.reconcile(&n.matter.mixture, n.matter.mass, 1e-6) {
                self.stats.surface_mismatches += 1;
                self.stats.worst_surface_mismatch = Some(why);
            }
            let epoch = self.tree.nodes[idx.get()].epoch;
            let n = &mut self.tree.nodes[idx.get()];
            n.surface_epoch = epoch;
            n.surface = Some(baked);
        }
        self.tree.nodes[idx.get()]
            .surface
            .as_ref()
            .unwrap_or_else(|| Self::nothing())
    }

    /// The empty surface, for a node that presents none.
    fn nothing() -> &'static crate::shape::Surface {
        static EMPTY: std::sync::OnceLock<crate::shape::Surface> = std::sync::OnceLock::new();
        EMPTY.get_or_init(crate::shape::Surface::default)
    }

    /// Emit the pieces. See [`World::surface_of_node`] for what decides which.
    fn bake_surface(&self, idx: NodeIdx) -> crate::shape::Surface {
        use crate::chem::{Phase, SubstanceId};
        use crate::shape::{Hull, Piece, Surface};

        let n = &self.tree.nodes[idx.get()];
        // What it is made of, and which substance that came from, so the bake
        // can be checked against the node's own pools.
        let measured = self.material_of(idx);
        let substance = n
            .matter
            .mixture
            .entries()
            .iter()
            .filter(|p| p.phase == Phase::Solid && p.fraction > 0.0)
            .max_by(|a, b| a.fraction.total_cmp(&b.fraction))
            .map(|p| p.substance)
            .unwrap_or(SubstanceId::UNSPECIATED);
        let material = match measured {
            Some(m) => m,
            // A node nobody has described falls back to what built it, exactly
            // as `surface_of` does — and to nothing at all if nothing built it,
            // because an undescribed cloud has no boundary to present.
            None => match (&n.topology, &n.morphology) {
                (Some(t), _) => t.material,
                (None, Some(m)) => m.material(),
                (None, None) => return Surface::default(),
            },
        };
        let solid_mass = if n.matter.is_described() {
            n.matter.mass * n.matter.solid_fraction()
        } else {
            n.matter.mass
        };

        // An assembly states its parts, each with its own material.
        //
        // First, and not merely as an optimisation: D18 attaches a material per
        // primitive, and this is the only thing in the engine that has ever had
        // more than one. Everything below — a tree, a sampled rock — is one
        // substance throughout and takes the node's single measured material,
        // which is why that path was enough until a box of oak panels on a
        // steel frame existed to break it. The parts are solids of stated size
        // rather than members inferred from joints, so nothing here needs the
        // topology at all.
        if let Some(a) = n.morphology.as_ref().and_then(|m| m.assembly()) {
            if !a.is_empty() {
                let pieces = a.pieces(&n.matter, &self.substances, material, substance);
                if !pieces.is_empty() {
                    return Surface::new(pieces);
                }
            }
        }

        // A structure states its members.
        if let (Some(mask), Some(t)) = (n.structural_mask(), n.topology.as_ref()) {
            let members = mask.iter().filter(|o| **o).count();
            if members > 0 {
                let share = solid_mass / members as f64;
                let mut pieces = Vec::with_capacity(members);
                for (i, ordered) in mask.iter().enumerate() {
                    if !ordered {
                        continue;
                    }
                    let radius = t.joints.get(i).map(|j| j.radius).unwrap_or(0.0);
                    let base = t.base.get(i).copied().unwrap_or(Vec3::ZERO);
                    let tip = t.tip.get(i).copied().unwrap_or(Vec3::ZERO);
                    // **A part that states a box presents one.** A coursed
                    // wall's blocks and a patch of ground's columns are filled
                    // solids of stated extent, and a capsule over the same two
                    // endpoints is a bead: measured at a 0.169 m scallop
                    // between blocks on a 1.2 m panel, and at a gap on every
                    // corner of a terrain grid, which is what something
                    // standing on the ground falls into.
                    //
                    // `Body::hull` is the one place that turns a body into
                    // geometry, so this never has to ask what kind of body it
                    // is holding.
                    let boxed = n.bodies.get(i).is_some_and(|b| b.is_boxed());
                    let hull = if boxed {
                        n.bodies[i].hull()
                    } else if (tip - base).norm2() > 0.0 {
                        Hull::capsule(base, tip, radius)
                    } else if let Some(b) = n.bodies.get(i) {
                        Hull::sphere(b.pos, b.radius.max(radius))
                    } else {
                        continue;
                    };
                    let mass = n.bodies.get(i).map(|b| b.mass).unwrap_or(share);
                    pieces.push(Piece { hull, material, substance, mass });
                }
                if !pieces.is_empty() {
                    return Surface::new(pieces);
                }
            }
        }

        // Unstructured matter: solid or nothing.
        if n.matter.is_described() && n.matter.solid_fraction() <= 0.0 {
            return Surface::default();
        }
        // **One piece, and it is a sphere.** A rock is one filled solid with no
        // cavity, so a single primitive is the honest description of it — and
        // the alternative is expensive twice over: hulling four thousand
        // sampled bodies costs an O(n^2) bake, and the hull it produces then
        // costs a four-thousand-sphere support query on every contact.
        // `PERFORMANCE.md`'s narrow-phase row prices that at 1.2 ms for one
        // pair.
        //
        // Materialised, the radius is *measured* rather than assumed: `Spread`
        // is Phase 1's own measurement of where a node's contents actually are,
        // and its furthest reach is the smallest sphere that contains them.
        // Unmaterialised, the node's own radius is all there is, and it is
        // exactly what the matter claims.
        let hull = if n.bodies.is_empty() {
            Hull::sphere(Vec3::ZERO, n.matter.radius)
        } else {
            let spread = crate::state::Spread::of(
                n.bodies.iter().map(|b| (b.pos, b.mass, b.radius)),
            );
            Hull::sphere(spread.centre, spread.furthest.max(1e-30))
        };
        Surface::new(vec![Piece { hull, material, substance, mass: solid_mass }])
    }

    /// What a node is made of, measured from its own matter. `PLAY.md` D13.
    ///
    /// `None` for matter nobody has described and for matter that is described
    /// and holds no solid — a cup of water has a ceramic material and the water
    /// in it has none, which is D13's own worked case.
    pub fn material_of(&self, idx: NodeIdx) -> Option<crate::material::Material> {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return None;
        }
        let n = &self.tree.nodes[idx.get()];
        let mixture = n.matter.mixture;
        if mixture.is_empty() {
            return None;
        }
        // The formation conditions of the first solid pool stand for the node's
        // — they are a property of the *node's* history rather than of one
        // substance in it, and `Material::measured` applies them to each pool
        // in turn.
        let solid = mixture
            .entries()
            .iter()
            .find(|p| p.phase == crate::chem::Phase::Solid && p.fraction > 0.0)?;
        let props = self.substances.get(solid.substance)?.props;
        // The same correction the packing gets: a node that states the volume
        // it fills is measured against that rather than against the sphere it
        // claims. See `World::packing_at`.
        let mut matter = n.matter;
        if let Some(v) = n
            .morphology
            .as_ref()
            .and_then(|m| m.recipe.as_ref().map(|r| r.solid_volume(m.growth())))
        {
            if v > 0.0 {
                matter.radius = (0.75 * v / std::f64::consts::PI).cbrt();
            }
        }
        let formation = crate::material::Formation::of_matter(&matter, &props);
        crate::material::Material::measured(&mixture, &self.substances, formation)
    }

    /// How much of this node's volume its own solid actually fills, 0 to 1.
    ///
    /// The packing `Formation::of_matter` measures, on its own rather than
    /// folded into a material's density — which is where it usually is, and
    /// which makes it useless for asking *whether* a node is condensed: a
    /// measured material's density is the node's own bulk density, so dividing
    /// one by the other is identically one.
    pub fn packing_of(&self, idx: NodeIdx) -> Option<f64> {
        let r = self.tree.nodes.get(idx.get())?.matter.radius;
        self.packing_at(idx, r)
    }

    /// As [`World::packing_of`], at a radius the caller has measured for
    /// itself.
    pub fn packing_at(&self, idx: NodeIdx, radius: f64) -> Option<f64> {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return None;
        }
        let n = &self.tree.nodes[idx.get()];
        let mut matter = n.matter;
        matter.radius = radius.max(1e-30);
        // **Against the volume the thing actually fills**, where it states one.
        // `Formation::of_matter` reads the packing off `matter.density()`,
        // which is the mass over the sphere the node claims — right for a node
        // that is the material and wrong for one that is an arrangement of it.
        // A recipe knows the volume it laid down, so the radius handed to the
        // formation is the sphere of *that* volume rather than the node's own.
        if let Some(v) = n
            .morphology
            .as_ref()
            .and_then(|m| m.recipe.as_ref().map(|r| r.solid_volume(m.growth())))
        {
            if v > 0.0 {
                matter.radius = (0.75 * v / std::f64::consts::PI).cbrt();
            }
        }
        let mixture = matter.mixture;
        let solid = mixture
            .entries()
            .iter()
            .find(|p| p.phase == crate::chem::Phase::Solid && p.fraction > 0.0)?;
        let props = self.substances.get(solid.substance)?.props;
        Some(crate::material::Formation::of_matter(&matter, &props).packing())
    }

    /// Fill in what this structure's joins are made of.
    ///
    /// `docs/PLAY.md` D15: "weld, glue and grown-together are not three
    /// features — they differ only in what the join is made of and therefore in
    /// its strength under D14". The substance travels with the joint from the
    /// recipe; the *material* is derived from it here, because this is where
    /// the substance registry is and because a material is a derivation rather
    /// than a lookup. There is no table of glues and nothing anywhere names
    /// one.
    ///
    /// Derived once per distinct substance rather than once per joint — the
    /// third axiom's "derived once for the kind and then run for each
    /// individual of it". A structure has thousands of joints and one or two
    /// kinds of join.
    ///
    /// Formation conditions are the node's own, exactly as they are for its
    /// bulk material and for each of its parts: glue set in a cold damp shed is
    /// not glue set in a kiln, and that is a property of where the thing is
    /// rather than of the glue.
    fn derive_bonds(&self, idx: NodeIdx, topo: &mut crate::topology::Topology) {
        use crate::chem::SubstanceId;
        topo.bonds.clear();
        // The node's own formation conditions, measured against the volume it
        // states rather than the sphere it claims — the same correction
        // `World::material_of` makes, because a joint's history is the node's
        // history and two answers to one question is the thing to avoid.
        let mut matter = self.tree.nodes[idx.get()].matter;
        if let Some(v) = self.tree.nodes[idx.get()]
            .morphology
            .as_ref()
            .and_then(|m| m.recipe.as_ref().map(|r| r.solid_volume(m.growth())))
        {
            if v > 0.0 {
                matter.radius = (0.75 * v / std::f64::consts::PI).cbrt();
            }
        }
        for j in &topo.joints {
            if j.bond == SubstanceId::UNSPECIATED || topo.bonds.iter().any(|(s, _)| *s == j.bond) {
                continue;
            }
            let Some(sub) = self.substances.get(j.bond) else { continue };
            let formation = crate::material::Formation::of_matter(&matter, &sub.props);
            topo.bonds.push((j.bond, crate::material::Material::of(&sub.props, formation)));
        }
    }

    /// Resolve the overlaps inside one node. `docs/PLAY.md` D3, the impulsive
    /// half.
    ///
    /// Overlap is `pairs(0.0)` — the gap a pair is separated by, at or below
    /// zero — so the same index that says what is *near* what says what is
    /// *inside* what, with no second traversal and no second notion of
    /// adjacency. That was the point of building one primitive.
    ///
    /// # Where velocity is written, and why not through the mailbox
    ///
    /// The impulse goes straight onto `Motion::velocity` for a promoted child,
    /// which is the same path D4 established for the force its parent's solver
    /// computes, and the only field that means "how fast this node is moving in
    /// its parent's frame". Posting it as momentum through `causal::Influence`
    /// would land on `matter.momentum` — the node's *internal* momentum, the
    /// motion of its contents about their own centre — which is a different
    /// quantity that happens to have a similar name.
    ///
    /// The heat does go through the mailbox, as an `Exchange`, because heat is
    /// energy arriving in another node's books and that is exactly what the
    /// mailbox is for. Velocity is not: the parent already owns where its
    /// children are going, and everything in this pass is inside one node,
    /// which is one causal cell by construction.
    ///
    /// # What it does not do
    ///
    /// It does not push interpenetrating things apart. A contact that is
    /// already separating is left alone — resolving it again is how two
    /// overlapping things end up vibrating against each other forever — but
    /// nothing here removes an existing overlap either. Positional correction
    /// is a solver concern and `PLAY.md` does not call for one.
    fn contact_within(&mut self, idx: NodeIdx) -> u64 {
        use crate::neighbourhood::{contact, Side};
        // Bake every promoted child's surface *before* the index is built.
        // Baking writes to the node and the index borrows the tree, and the
        // alternative — rebuilding a proxy inside the loop — is the per-frame
        // derivation `PLAY.md` D13 retires.
        let children: Vec<NodeIdx> = self.tree.nodes[idx.get()]
            .children
            .iter()
            .copied()
            .filter(|c| !c.is_none())
            .collect();
        for c in children {
            self.surface_of_node(c);
        }
        // And this node's own, which is what its loose bodies present.
        self.surface_of_node(idx);
        let nb = self.tree.neighbourhood(idx);
        if nb.len() < 2 {
            return 0;
        }
        let Some(overlaps) = nb.pairs(0.0) else { return 0 };
        if overlaps.is_empty() {
            return 0;
        }
        // The material a plain body presents is its own node's: a structure's
        // members are made of what the structure is made of.
        let mine = self.surface_of(idx);

        let side_of = |w: &World, i: usize| -> Option<(Occupant, Side)> {
            let (occ, pos, radius) = nb.at(i)?;
            let side = match occ {
                Occupant::Body(k) => {
                    let n = &w.tree.nodes[idx.get()];
                    let b = n.bodies.get(k as usize)?;
                    // A member of this node is a capsule, for the same reason a
                    // promoted child's members are: `Topology` carries `base`,
                    // `tip` and a cross-section radius, and the sphere at the
                    // midpoint is 164x shorter than the beam the renderer draws.
                    // This is the side a limb landing on a tree strikes.
                    //
                    // **Unless the recipe stated a box**, which is what a
                    // panel, a slab of ground or a floor is. The beam through a
                    // bucket's floor has the floor's seam for a radius, and a
                    // ball sinking through the water found it 4 cm under the
                    // floor's top — and went through. A stated solid is its own
                    // shape, and `Body::hull` is what turns one into geometry.
                    let member = n.topology.as_ref().and_then(|t| {
                        let i = k as usize;
                        let r = t.joints.get(i)?.radius;
                        if r > 0.0 && b.is_boxed() {
                            return Some(b.hull());
                        }
                        let (base, tip) = (*t.base.get(i)?, *t.tip.get(i)?);
                        (r > 0.0 && (tip - base).norm2() > 0.0)
                            .then(|| crate::shape::Hull::capsule(base, tip, r))
                    });
                    // A member of a structure whose contents carry the node's
                    // weight is held by the same support that weight assumes,
                    // so a contact meets it as immovable: see `held` below.
                    let mass = if n.gravity != crate::math::Vec3::ZERO { f64::INFINITY } else { b.mass };
                    match member {
                        Some(h) => Side::shaped(
                            pos,
                            b.vel,
                            mass,
                            b.heat_capacity(),
                            mine?,
                            vec![h],
                        ),
                        None => Side::sphere(pos, b.vel, b.mass, radius, b.heat_capacity(), mine?),
                    }
                }
                Occupant::Child(c) => {
                    let n = &w.tree.nodes[c.get()];
                    // A promoted child is the engine's rigid body: one
                    // velocity, one spin, and `apply_contact` has always put an
                    // impulse straight onto them. What it lacked was a shape.
                    //
                    // **The stored one**, baked above. It is in the child's own
                    // frame, so it is turned by the child's orientation and
                    // *then* moved to where the child is — the same composition
                    // `Motion::body_to_parent` uses for a body-fixed point.
                    // Translating alone left a spinning box's walls in the axes
                    // they were built in.
                    let baked = match &n.surface {
                        Some(s) if !s.is_empty() => s,
                        // A node that presents no surface is not collidable,
                        // which is D13's own line for a liquid or a gas.
                        _ => return None,
                    };
                    let shape = baked.placed(n.motion.orientation, pos);
                    let _ = radius;
                    let mut side = Side::shaped(
                        pos,
                        n.motion.velocity,
                        n.matter.mass,
                        n.matter.heat_capacity(),
                        w.surface_of(c)?,
                        shape,
                    );
                    // D18 puts the material on the piece, so the contact reads
                    // the one it struck rather than an average of the house.
                    side.per_piece = baked
                        .pieces()
                        .iter()
                        .map(|p| crate::neighbourhood::Resilience::of(&p.material))
                        .collect();
                    side
                }
            };
            Some((occ, side))
        };

        // Which of this node's bodies are liquid it holds loose — the parcels
        // weakly-compressible SPH couples to a solid child through the child's
        // own wall (`hydro::Wall`).
        let liquid_loose: Option<Vec<bool>> = {
            let n = &self.tree.nodes[idx.get()];
            (n.matter.mixture.in_phase(crate::chem::Phase::Liquid) > 0.0).then(|| {
                let mask = n.structural_mask();
                (0..n.bodies.len())
                    .map(|k| !mask.as_ref().map(|m| m[k]).unwrap_or(false))
                    .collect()
            })
        };

        let mut resolved = 0u64;
        let now = self.time;
        for (i, j) in overlaps {
            let (Some((oa, a)), Some((ob, b))) = (side_of(self, i), side_of(self, j)) else {
                continue;
            };
            // **And a solid child in a liquid is the solver's already.** The
            // parcels meet it as a wall and give every push back, which is its
            // buoyancy and its drag; an impulse here as well would couple it
            // twice. Measured before this line: 371 contacts resolved between
            // a floating ball and the water holding it up, on top of the
            // pressure that was holding it up.
            if let Some(loose) = &liquid_loose {
                let pair = match (oa, ob) {
                    (Occupant::Child(c), Occupant::Body(k)) | (Occupant::Body(k), Occupant::Child(c)) => {
                        Some((c, k))
                    }
                    _ => None,
                };
                if let Some((c, k)) = pair {
                    let solid = crate::eos::Eos::of_matter(&self.tree.nodes[c.get()].matter, &self.substances)
                        .condensed()
                        .is_some();
                    if solid && loose.get(k as usize).copied().unwrap_or(false) {
                        continue;
                    }
                }
            }
            // Two bodies of the same node are not in contact, they are in it
            // *together*, and whatever couples them is already running: the
            // structure solver if the node is a structure, the tier's own
            // solver if it is a continuum. Contact is for what that coupling
            // does not reach — a promoted child against another, or against the
            // bodies of the node it is sitting in. That is the same line
            // `drop_fragments` already draws for a limb landing on a tree,
            // arrived at from the other direction.
            //
            // Without this a structure would take its own joints as collisions.
            // Measured, on materialised programs: a `Wall` has 527 overlapping
            // pairs among 55 members, because a wall is courses of blocks
            // packed against each other, and it would come apart on the frame
            // it was materialised. A `Tower`, a `Tree` and a `Settlement` have
            // none — their members are long and thin and sit a member-length
            // apart — so the hazard is real without being universal, which is
            // exactly the kind that gets shipped.
            if matches!((oa, ob), (Occupant::Body(_), Occupant::Body(_))) {
                continue;
            }
            // **A member held by its structure takes no velocity from a
            // contact.** In a node whose contents carry its weight, the frame
            // is held up by definition — that is what makes the weight a
            // weight — and a member is the structure doing the holding. Giving
            // it the impulse instead handed a bucket's floor panel a velocity
            // nothing ever integrates: every contact split its impulse with a
            // floor that appeared to recede, and a ball resting on it gained
            // downward speed without end. The momentum goes to the support, as
            // the momentum gravity adds comes from it.
            let held = |o: Occupant| match o {
                Occupant::Body(k) => {
                    self.tree.nodes[idx.get()].gravity != crate::math::Vec3::ZERO
                        && self.tree.nodes[idx.get()]
                            .topology
                            .as_ref()
                            .and_then(|t| t.joints.get(k as usize))
                            .map(|j| j.radius > 0.0)
                            .unwrap_or(false)
                }
                _ => false,
            };
            let (held_a, held_b) = (held(oa), held(ob));
            // **And against what holds it, a child is moved out of the
            // overlap.** An impulse removes the approach and leaves the overlap
            // where it was; for a thing resting in a field, each solve adds
            // `g dt` of approach back and the overlap deepens by that every
            // time — measured, a ball twice as dense as water came to the floor
            // of a bucket, took 378 contacts and went through it. A child owns
            // where it is, so the child moves. Only against a held member,
            // because that is where momentum already goes to the support:
            // between two free things a shift of position would move `sum r x
            // p`, and a contact between them conserves that exactly.
            if held_a != held_b {
                if let Some((depth, axis)) = crate::neighbourhood::overlap(&a, &b) {
                    let (child, sign) = if held_a { (ob, 1.0) } else { (oa, -1.0) };
                    if let Occupant::Child(ch) = child {
                        if !ch.is_none() {
                            let n = &mut self.tree.nodes[ch.get()];
                            n.motion.offset = n.motion.offset + axis.scale(sign * depth);
                        }
                    }
                }
            }
            let Some(c) = contact(&a, &b) else { continue };
            let total = c.normal + c.friction;
            if !total.is_finite() {
                continue;
            }
            // `a` takes the negative of every impulse `b` takes. Written this
            // way round so momentum conservation is a property of the code
            // rather than something a test has to keep watch on — and a held
            // member's share goes to its support, with its heat still arriving:
            // what the contact did not give back is made at the interface,
            // whoever is holding what.
            let zero = crate::math::Vec3::ZERO;
            if held_a {
                apply_contact(self, idx, oa, zero, zero, c.heat_a, now);
            } else {
                apply_contact(self, idx, oa, total.scale(-1.0), c.spin_a, c.heat_a, now);
            }
            if held_b {
                apply_contact(self, idx, ob, zero, zero, c.heat_b, now);
            } else {
                apply_contact(self, idx, ob, total, c.spin_b, c.heat_b, now);
            }
            resolved += 1;
        }
        self.stats.contacts_resolved += resolved;
        resolved
    }

    /// Move heat between the things inside one node that are next to each
    /// other. `docs/PLAY.md` D3, the transport half.
    ///
    /// This is what "a hot node beside a cold one equilibrates without either
    /// being told the other exists" is made of, and the striking thing about it
    /// is how little of it is here: the adjacency comes from
    /// [`Neighbourhood`](crate::neighbourhood::Neighbourhood), the transport
    /// from [`exchange`](crate::neighbourhood::exchange), the coefficient from
    /// [`radiative_conductance`](crate::neighbourhood::radiative_conductance),
    /// and what is left is bookkeeping. Conduction and diffusion will be two
    /// more coefficients and no more code than that, which is the whole claim
    /// D3 was making.
    ///
    /// # Bodies and children are not two cases
    ///
    /// An occupant is a materialised body or a promoted child, and the physics
    /// does not distinguish them — but where the answer has to be *written*
    /// does. A body is inside this node and inside its clock, so it is written
    /// directly. A child is a node of its own with its own clock and its own
    /// rate, and writing into another node's matter behind its clock is exactly
    /// the bug `causal::Mailbox` exists to prevent, so a child's share is
    /// posted as an [`InfluenceKind::Exchange`] and arrives when *its* clock
    /// reaches it.
    ///
    /// Both spellings conserve, and they conserve for the same reason rather
    /// than by separate arrangement. A promoted child's stand-in body carries a
    /// copy of the child's `internal_energy`, refreshed by `sync_children` at
    /// the top of every solve, so energy handed to the child reappears in this
    /// node's own body list next frame — and energy taken out of a plain body
    /// here and given to a child leaves the body list by one route and comes
    /// back by the other. The sum over this node's bodies is unchanged either
    /// way.
    ///
    /// # Why the potentials are updated as the walk goes
    ///
    /// One occupant can have twenty neighbours, and reading every temperature
    /// once and applying every transfer against those stale readings — Jacobi —
    /// lets a cold thing with twenty hot neighbours be given twenty separate
    /// shares of "the whole way to equilibrium" and end up hotter than any of
    /// them. Updating in place as each pair is settled costs nothing, cannot
    /// overshoot because each transfer is individually bounded by that pair's
    /// own equilibrium, and is deterministic because `pairs` is. Energy is
    /// conserved under either, but only one of them is *right*.
    fn exchange_within(&mut self, idx: NodeIdx, dt: f64) -> u64 {
        if !(dt > 0.0) || !dt.is_finite() {
            return 0;
        }
        // Only where a blackbody is what the node's contents actually are.
        // Identical gate, and identical reason, to `evolve_matter`: a galactic
        // node's `temperature` is a velocity dispersion, and two star clusters
        // do not radiate at each other as blackbodies the size of star
        // clusters. The tier is the statement of which physics applies.
        if self.tree.nodes[idx.get()].tier < crate::units::Tier::Planetary {
            return 0;
        }
        let nb = self.tree.neighbourhood(idx);
        if nb.len() < 2 {
            return 0;
        }
        // The node's own resolution, not the index's reach: see
        // `Neighbourhood::resolution`. Above this length the structure belongs
        // to the parent, and the parent's luminosity field — `environment_at`
        // — is already carrying it.
        let Some(pairs) = nb.pairs(nb.resolution()) else {
            return 0;
        };
        if pairs.is_empty() {
            return 0;
        }

        // Read every potential once, then let the walk move them.
        let mut side: Vec<Reservoir> = Vec::with_capacity(nb.len());
        // And how well each conducts, W/m/K — `transport.rs`. A body is made of
        // its own substance where it names one and of its node's mixture where
        // it does not; a child is made of its own matter. Matter nobody has
        // described conducts nothing, rather than something invented.
        let mut conducts: Vec<f64> = Vec::with_capacity(nb.len());
        let own = {
            let m = &self.tree.nodes[idx.get()].matter;
            crate::transport::conductivity_of(&m.mixture, &self.substances, m.temperature)
        };
        for i in 0..nb.len() {
            let (occ, _, _) = nb.at(i).expect("index is in range");
            match occ {
                Occupant::Body(k) => {
                    let b = &self.tree.nodes[idx.get()].bodies[k as usize];
                    side.push(Reservoir::new(b.temperature, b.heat_capacity()));
                    let named = self.substances.get(b.substance).map(|s| {
                        let phase = if b.temperature >= s.props.boiling_point {
                            crate::chem::Phase::Gas
                        } else if b.temperature >= s.props.melting_point {
                            crate::chem::Phase::Liquid
                        } else {
                            crate::chem::Phase::Solid
                        };
                        crate::transport::conductivity(&s.props, phase, b.temperature)
                    });
                    conducts.push(named.or(own).unwrap_or(0.0));
                }
                Occupant::Child(c) => {
                    let m = &self.tree.nodes[c.get()].matter;
                    side.push(Reservoir::new(m.temperature, m.heat_capacity()));
                    conducts.push(
                        crate::transport::conductivity_of(&m.mixture, &self.substances, m.temperature)
                            .unwrap_or(0.0),
                    );
                }
            }
        }
        let mut moved = vec![0.0f64; nb.len()];

        let mut crossings = 0u64;
        for (i, j) in pairs {
            let (_, pi, ri) = nb.at(i).expect("pair index is in range");
            let (_, pj, rj) = nb.at(j).expect("pair index is in range");
            let d = (pj - pi).norm();
            let area = crate::neighbourhood::radiative_area(ri, rj, d);
            // Radiation across whatever gap there is, and conduction across
            // the face they share when there is none: D3's "three calls to one
            // function with different coefficients", two of them now.
            let g = crate::neighbourhood::radiative_conductance(
                side[i].potential,
                side[j].potential,
                area,
            ) + match (nb.at(i).map(|o| o.0), nb.at(j).map(|o| o.0)) {
                // Two parcels of this node's own medium are neighbouring cells
                // of it; anything involving a separate object touches only
                // where the two surfaces meet. `transport.rs`.
                (Some(Occupant::Body(_)), Some(Occupant::Body(_))) => {
                    crate::transport::continuum_conductance(conducts[i], conducts[j], d)
                }
                _ => crate::transport::conductive_conductance(conducts[i], ri, conducts[j], rj, d),
            };
            let q = crate::neighbourhood::exchange(side[i], side[j], g, dt);
            if !q.is_finite() || q == 0.0 {
                continue;
            }
            // Move the potentials with the heat, so the next pair sees the
            // state this one left behind.
            if side[i].capacity.is_finite() && side[i].capacity > 0.0 {
                side[i].potential -= q / side[i].capacity;
            }
            if side[j].capacity.is_finite() && side[j].capacity > 0.0 {
                side[j].potential += q / side[j].capacity;
            }
            moved[i] -= q;
            moved[j] += q;
            crossings += 1;
        }

        let now = self.time;
        for i in 0..nb.len() {
            if moved[i] == 0.0 || !moved[i].is_finite() {
                continue;
            }
            let (occ, pos, _) = nb.at(i).expect("index is in range");
            match occ {
                Occupant::Body(k) => {
                    if let Some(b) = self.tree.nodes[idx.get()].bodies.get_mut(k as usize) {
                        b.add_heat(moved[i]);
                    }
                }
                Occupant::Child(c) => {
                    // The separation from the node's centre is the light delay
                    // that applies; within a node it is far below a frame, and
                    // it stops being so exactly when the node is large enough
                    // that it should be.
                    self.mailbox.post(
                        c,
                        now,
                        pos.norm(),
                        InfluenceKind::Exchange,
                        moved[i],
                        Vec3::ZERO,
                    );
                }
            }
        }
        self.stats.exchange_crossings += crossings;
        crossings
    }

    /// The nuclear tier does not integrate trajectories; it samples events.
    fn advance_statistical(&mut self, idx: NodeIdx, dt: f64) -> solvers::SolveReport {
        let (key, epoch, count) = {
            let n = &self.tree.nodes[idx.get()];
            (n.key, n.epoch, n.bodies.len())
        };
        let seed = self.tree.world_seed;
        let mut stream = Stream::at(seed, key.0, epoch, Purpose::Decay);
        let before = solvers::measure(&self.tree.nodes[idx.get()].bodies, 0.0);
        let mut released = 0.0;
        let bodies = &mut self.tree.nodes[idx.get()].bodies;
        for b in bodies.iter_mut() {
            // Free neutrons decay; this is the one process fast enough to
            // matter on the timescales a user watching a nucleus experiences.
            if b.kind == crate::state::BodyKind::Nucleon && b.charge == 0.0 {
                let iso = solvers::nuclear::Isotope::Neutron;
                let n_nuclei = (b.mass / M_NEUTRON).max(0.0);
                let decays = iso.sample_decays(n_nuclei, dt, &mut stream);
                if decays > 0.0 {
                    released += decays * iso.q_value();
                    b.charge += decays * E_CHARGE;
                    b.internal_energy += decays * iso.q_value();
                }
            }
            b.pos += b.vel.scale(dt);
        }
        let after = solvers::measure(&self.tree.nodes[idx.get()].bodies, 0.0);
        solvers::SolveReport {
            steps: 1,
            interactions: count as u64,
            dt_used: dt,
            before,
            after,
            non_mechanical_energy: released,
            unrest: 0.0,
        }
    }

    /// Advance one structure's developmental state.
    ///
    /// The environment is read off the node's own matter, so a structure in
    /// a cold or crowded node grows slowly without anyone having to arrange it.
    /// The transaction is validated before it is applied: a growth program
    /// cannot mint free energy or order, it can only trade for them.
    pub fn grow_node(&mut self, idx: NodeIdx, dt: f64) -> Option<crate::morph::GrowthStep> {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive || dt <= 0.0 {
            return None;
        }
        let env = self.environment_at(idx);
        let node = &mut self.tree.nodes[idx.get()];
        let morph = node.morphology.as_mut()?;
        let was = (morph.built, morph.progress);
        let txn = morph.advance(dt, &env);
        if txn.validate().is_err() {
            // A program that cannot balance its books does not get to run. This
            // is a bug in the program, not a condition to be smoothed over.
            self.rejected_growth_steps += 1;
            return None;
        }
        let extent = morph.extent().max(1e-30);
        let stored = morph.stored_energy();
        let (morph_built, morph_progress) = (morph.built, morph.progress);

        // Apply the growth step to the node's matter. Mass moves *within* the node
        // — carbon from its air into its wood — so mass, composition and baryon
        // number are all unchanged, and only the energy and entropy accounts
        // move. What crosses the boundary is energy, and it is booked.
        node.matter.chemical_energy = stored;
        // Only the thermalised share stays. What was re-radiated has left the
        // node, and adding it here would cook a forest in a season.
        node.matter.internal_energy += txn.heat_released;
        node.matter.entropy += txn.entropy_local;
        node.matter.entropy_exported += txn.entropy_exported;
        node.matter.radius = extent;
        node.matter.luminosity = crate::state::stefan_boltzmann(extent, node.matter.temperature);

        // **Only if the structure it would generate has actually changed.** A
        // recipe is a pure function of itself and `(built, progress)`, so when
        // neither moved the drawing is the same drawing and there is nothing
        // stale to throw away.
        //
        // Discarding unconditionally was cheap for a tree, which grows every
        // frame anyway, and ruinous for ground, which grows never: a planet's
        // whole surface tree was folded back on every frame it was looked at,
        // so an observer descending it re-promoted seven patches a frame and
        // arrived nowhere. Measured: 42 promotions over six frames with two
        // live nodes at the end of them.
        let changed = (morph_built, morph_progress) != was;
        if changed {
            let node = &mut self.tree.nodes[idx.get()];
            node.bodies.clear();
            // The children it had promoted out of it are folded back before
            // the fresh body list replaces their slots.
            self.tree.shed_children(idx);
            // It grew, so it may not be the size of thing it was. See
            // `Tree::retier`.
            self.tree.retier(idx);
        }
        self.tree.stats.growth_steps += 1;
        self.tree.stats.external_energy_absorbed += txn.net_boundary_flux();
        Some(txn)
    }

    /// Load a structure and see what survives.
    ///
    /// The whole point of topology: nothing here says "lightning destroys a
    /// tree" or "wet snow breaks branches". The insult produces forces,
    /// temperatures and energy deposition, and then the ordinary stress
    /// calculation decides what fails. A limb comes down because the moment at
    /// its base exceeded what its cross-section could carry, which is also why
    /// real limbs come down.
    ///
    /// Damage persists: broken joints become events in the developmental state,
    /// so the structure regenerates broken for the rest of its life.
    pub fn damage(
        &mut self,
        idx: NodeIdx,
        mechanisms: &[crate::solvers::structure::Mechanism],
    ) -> DamageOutcome {
        use crate::solvers::structure as st;
        let mut out = DamageOutcome::default();
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return out;
        }
        if self.tree.nodes[idx.get()].morphology.is_none() {
            return out;
        }
        self.disturb(idx);
        self.tree.refine(idx);

        let ambient = self.tree.nodes[idx.get()].matter.temperature;
        let bodies = self.tree.nodes[idx.get()].bodies.clone();
        let mut topo = match self.tree.nodes[idx.get()].topology.clone() {
            Some(t) => t,
            None => return out,
        };
        self.derive_bonds(idx, &mut topo);
        let structural_mass: f64 = self
            .tree
            .nodes[idx.get()]
            .morphology
            .as_ref()
            .map(|m| m.built)
            .unwrap_or(0.0);
        // Nothing left to load. Without this the density correction divides by
        // a vanishing volume, member radii go to zero, and the stresses come
        // back as 10^16 — a collapsed structure reported as an infinitely
        // overloaded one.
        if structural_mass <= COLLAPSE_MASS {
            out.collapsed = true;
            return out;
        }

        // Gravity is always present; the caller supplies whatever else is
        // happening. Mechanisms compose, so a structure can be burning, iced
        // and in a gale at once and the solver never learns those words.
        let mut field = st::LoadField::new(bodies.len(), ambient);
        for m in mechanisms {
            field.apply(m, &bodies, &topo);
        }
        // Gravity last, and once. It acts on the accreted mass as well as the
        // structure's own, so it has to follow anything that adds mass — and
        // applying it on both sides of that would weigh the structure twice.
        //
        // In the field the node is actually in, derived from what it is inside.
        // A structure on a moon carries a sixth of the weight it would on Earth
        // and fails under a sixth of the snow.
        field.apply(&st::weather::gravity(self.tree.gravity_at(idx)), &bodies, &topo);

        let (loads, indeterminate, iters) = st::analyse_with(&bodies, &topo, &field);
        let failures = st::apply_failures(&bodies, &mut topo, &loads, &field);

        out.peak_utilisation = failures.peak_utilisation;
        out.broken_joints = failures.broken_sites.len();
        out.detached_mass = failures.detached_mass;
        out.consumed_mass = failures.consumed_mass;
        out.energy_delivered = failures.energy_delivered;
        out.indeterminate = indeterminate;
        out.solver_iterations = iters;
        let insult_report = &failures;

        // Fold the damage into the developmental state, so it survives the
        // structure being discarded and rebuilt.
        // The event log records what is *observable*, not every twig. A crown
        // fire breaks thousands of joints; a bounded log that kept the oldest
        // of them would fill with the first few hundred twigs and drop the
        // major limb that went later. So breaks are ranked by the mass they
        // were carrying and only the significant ones are named — the rest are
        // captured by the structure's mass loss, which is all an observer could
        // detect anyway.
        //
        // This is the same resolution-scoped rule the measurement ledger uses:
        // commit what someone could tell apart, regenerate the rest.
        let mut ranked: Vec<(f64, u32)> = insult_report
            .broken_sites
            .iter()
            .chain(failures.broken_sites.iter())
            .map(|&site| {
                let carried = topo
                    .site
                    .iter()
                    .position(|&s| s == site)
                    .and_then(|i| loads.get(i))
                    .map(|l| l.carried)
                    .unwrap_or(0.0);
                (carried, site)
            })
            .collect();
        ranked.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut sites: Vec<u32> = ranked.into_iter().map(|(_, s)| s).collect();
        sites.dedup();
        sites.truncate(NOTABLE_BREAKS);

        // Everything that lost its support is now its own falling object.
        if !failures.broken_members.is_empty() {
            let cut = st::detach(&topo, &failures.broken_members);
            for piece in &cut.pieces {
                if self.falling.len() >= MAX_FALLING {
                    break;
                }
                if let Some(frag) = st::Fragment::new(&bodies, &topo, piece) {
                    out.detached_pieces += 1;
                    self.falling.push((idx, frag));
                }
            }
        }

        // **An assembled thing loses a part, not a fraction.** D15's break: the
        // panel whose seam failed becomes a node of its own at that moment, and
        // the recipe keeps six parts with five joins and a break. `sever_many`
        // is the grown case — a tree has no parts list to take a member out of,
        // so it books the loss as a fraction of its structural mass and
        // regenerates a smaller tree — and running both would take the mass
        // twice.
        if self.tree.nodes[idx.get()].morphology.as_ref().is_some_and(|m| m.is_assembled()) {
            let mut lost = 0.0;
            for &site in &sites {
                let part = self.detach(idx, site);
                if !part.is_none() {
                    lost += self.tree.nodes[part.get()].matter.mass;
                    out.detached_pieces += 1;
                }
            }
            out.detached_mass = lost;
            self.tree.stats.damage_events += 1;
            // Deliberately not clearing the body list: the parts that just
            // became nodes are promoted children of those very slots, and
            // clearing `children` would orphan them. The recipe changed in its
            // joins and not in its geometry, so the bodies it generates are
            // still the bodies it generated.
            return out;
        }

        let node = &mut self.tree.nodes[idx.get()];
        let temperature = node.matter.temperature;
        let made_of = node.matter.composition;
        if let Some(m) = node.morphology.as_mut() {
            if !sites.is_empty() && structural_mass > 0.0 {
                let fraction = (failures.detached_mass / structural_mass).clamp(0.0, 1.0);
                let txn = m.sever_many(&sites, fraction, made_of);
                out.detached_mass = txn.mass_detached;
            }
            if insult_report.consumed_mass > 0.0 {
                let burn = m.consume(insult_report.consumed_mass, temperature, made_of);
                if burn.validate().is_ok() {
                    // Burning releases the free energy the wood was holding.
                    // The atoms stay in the node as combustion products, so mass
                    // and baryon number are untouched; only the energy moves.
                    node.matter.chemical_energy -= burn.energy_released;
                    node.matter.internal_energy += burn.heat_released;
                    node.matter.entropy += burn.entropy_local;
                    node.matter.entropy_exported += burn.entropy_exported;
                    out.energy_released = burn.energy_released;
                } else {
                    self.rejected_growth_steps += 1;
                }
            }
            node.matter.chemical_energy = m.stored_energy();
            node.matter.radius = m.extent().max(1e-30);
        }
        // The structure it would generate has changed.
        node.bodies.clear();
        node.topology = None;
        // The borrow of `node` ends here; the children it had promoted out of
        // it are folded back before the fresh body list replaces their slots.
        self.tree.shed_children(idx);
        // What is left of it may be a different size of thing.
        self.tree.retier(idx);
        // `docs/PLAY.md` D19: the break is now in `Morphology::events`, which
        // regenerates it, so this is an edit and *not* a pin. The node keeps no
        // body list, collapses to its recipe when nobody is watching, and comes
        // back with the same members missing. What it must not do is be
        // redrawn from the ensemble, which would quietly mend it — that is the
        // one thing the flag stops.
        self.tree.record_edit(idx);
        self.tree.stats.damage_events += 1;
        out
    }

    /// Integrate a structure through real time under a set of mechanisms.
    ///
    /// [`World::damage`] asks whether a structure stands up under a load.
    /// This asks what it *does* while that load is on it, which is a different
    /// question with a different answer: a gust a quasi-static check passes at
    /// 60% utilisation can break the same member outright, because a load
    /// arriving suddenly deflects a structure about twice as far as the same
    /// load standing still.
    ///
    /// The structure's dynamic state persists between calls — it has to, since
    /// what it does next depends on how it is already moving — and is
    /// discarded the moment the observer looks at something else.
    pub fn shake(
        &mut self,
        idx: NodeIdx,
        mechanisms: &[crate::solvers::structure::Mechanism],
        seconds: f64,
    ) -> ShakeOutcome {
        use crate::solvers::structure as st;
        let mut out = ShakeOutcome::default();
        if idx.is_none() || seconds <= 0.0 || !self.tree.nodes[idx.get()].alive {
            return out;
        }
        if self.tree.nodes[idx.get()].morphology.is_none() {
            return out;
        }
        self.disturb(idx);
        self.tree.refine(idx);
        let ambient = self.tree.nodes[idx.get()].matter.temperature;
        let bodies = self.tree.nodes[idx.get()].bodies.clone();
        let mut topo = match self.tree.nodes[idx.get()].topology.clone() {
            Some(t) => t,
            None => return out,
        };
        self.derive_bonds(idx, &mut topo);
        let topo = topo;

        // Rebuild whenever the structure itself has changed. Carrying a stale
        // dynamic state across a regeneration would let a tree that has lost a
        // limb keep swinging it.
        let existing = self
            .shaking
            .iter()
            .position(|(n, ds)| *n == idx && ds.tip_node.len() == bodies.len());
        let slot = match existing {
            Some(k) => k,
            None => {
                self.shaking.retain(|(n, _)| *n != idx);
                let Some(ds) = st::dynamic_structure(&bodies, &topo) else {
                    return out;
                };
                if self.shaking.len() >= MAX_SHAKEN {
                    // The oldest goes: whatever the observer stopped watching
                    // first is what they are least likely to look back at.
                    self.shaking.remove(0);
                }
                self.shaking.push((idx, ds));
                self.shaking.len() - 1
            }
        };
        let ds = &mut self.shaking[slot].1;

        let mut field = st::LoadField::new(bodies.len(), ambient);
        for m in mechanisms {
            field.apply(m, &bodies, &topo);
        }
        field.apply(&st::weather::gravity(self.tree.gravity_at(idx)), &bodies, &topo);

        // Substep to whatever the structure's own period demands. A sway with
        // a two-second period is resolved by a twentieth of a second; a steel
        // frame's is a hundred times shorter, and integrating it at the frame
        // rate would report a structure that does not move.
        let period = ds.dynamics.natural_period();
        let target = if period > 0.0 { period / 20.0 } else { seconds };
        let steps = (seconds / target).ceil().max(1.0).min(MAX_SHAKE_STEPS) as usize;
        let h = seconds / steps as f64;
        let mut broken: Vec<usize> = Vec::new();
        for _ in 0..steps {
            let rep = ds.advance(&field, h);
            out.iterations += rep.iterations;
            out.released += rep.released;
            out.dissipated += rep.dissipated;
            broken.extend(rep.broken);
            if !rep.converged {
                out.diverged = true;
                break;
            }
        }
        out.steps = steps as u32;
        out.kinetic = ds.dynamics.kinetic_energy();
        out.strain = ds.dynamics.strain_energy();
        out.displacement = ds
            .dynamics
            .displacement
            .iter()
            .map(|d| d.t.norm())
            .fold(0.0f64, f64::max);
        out.displacement_ratio = ds.dynamics.displacement_ratio();

        if broken.is_empty() {
            return out;
        }

        // Fold the failures into the developmental state, exactly as the static
        // path does, so a limb lost in a gust survives the structure being
        // discarded and rebuilt.
        let members = ds.failed_members(&broken);
        out.broken_joints = members.len();

        // What came away is its own object now, with its own roots, and it is
        // falling. Re-rooting is not cosmetic: the static analysis walks the
        // support forest towards its roots, and a branch still carrying its old
        // support index would be analysed as though the trunk were holding it
        // up — which is exactly what stopped being true.
        let cut = st::detach(&topo, &members);
        for piece in &cut.pieces {
            if self.falling.len() >= MAX_FALLING {
                break;
            }
            if let Some(frag) = st::Fragment::new(&bodies, &topo, piece) {
                out.detached_pieces += 1;
                self.falling.push((idx, frag));
            }
        }
        let mut ranked: Vec<(f64, u32)> = members
            .iter()
            .filter_map(|&m| {
                let i = m as usize;
                topo.site.get(i).map(|&site| (bodies[i].mass, site))
            })
            .collect();
        ranked.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut sites: Vec<u32> = ranked.iter().map(|&(_, s)| s).collect();
        sites.dedup();
        sites.truncate(NOTABLE_BREAKS);
        let detached: f64 = ranked.iter().map(|&(m, _)| m).sum();

        let structural_mass = self.tree.nodes[idx.get()]
            .morphology
            .as_ref()
            .map(|m| m.built)
            .unwrap_or(0.0);
        let node = &mut self.tree.nodes[idx.get()];
        let made_of = node.matter.composition;
        if let Some(m) = node.morphology.as_mut() {
            if !sites.is_empty() && structural_mass > 0.0 {
                let fraction = (detached / structural_mass).clamp(0.0, 1.0);
                let txn = m.sever_many(&sites, fraction, made_of);
                out.detached_mass = txn.mass_detached;
            }
            node.matter.chemical_energy = m.stored_energy();
            node.matter.radius = m.extent().max(1e-30);
        }
        node.bodies.clear();
        node.topology = None;
        // The borrow of `node` ends here; the children it had promoted out of
        // it are folded back before the fresh body list replaces their slots.
        self.tree.shed_children(idx);
        self.shaking.retain(|(n, _)| *n != idx);
        // Likewise: a structure that has shed members is smaller than it was.
        self.tree.retier(idx);
        // D19 again, and for the same reason as the strike path above.
        self.tree.record_edit(idx);
        self.tree.stats.damage_events += 1;
        out
    }

    /// Pieces currently falling from a node.
    pub fn falling(&self) -> &[(NodeIdx, crate::solvers::structure::Fragment)] {
        &self.falling
    }

    /// Advance every falling piece, and hand what they hit to the structures
    /// they hit it with.
    ///
    /// This is what makes a break more than a bookkeeping entry. A limb that
    /// comes away is re-rooted into its own object — its old support index
    /// stopped being true the moment it broke — and then it falls, and what it
    /// lands on gets the impulse. The impulse goes in through the ordinary
    /// mechanism vocabulary, so what happens next is the same stress
    /// calculation that decides everything else, and a limb heavy enough to
    /// break what it lands on produces another falling limb.
    pub fn drop_fragments(&mut self, dt: f64) -> FallReport {
        use crate::solvers::frame::Dof;
        use crate::solvers::structure as st;
        let mut report = FallReport::default();
        if self.falling.is_empty() || dt <= 0.0 {
            return report;
        }

        // What the debris is falling onto has to exist. `damage` clears a
        // node's detail when it breaks something, so without this the
        // collision test runs against nothing and every piece falls through
        // the tree it came off.
        let nodes: Vec<NodeIdx> = {
            let mut v: Vec<NodeIdx> = self.falling.iter().map(|(n, _)| *n).collect();
            v.sort_unstable_by_key(|n| n.get());
            v.dedup();
            v
        };
        let mut struck: HashMap<usize, (Vec<crate::state::Body>, crate::topology::Topology)> =
            HashMap::new();
        // The field each of them is actually falling in, derived from what they
        // are inside rather than assumed. `PLAY.md` D6: the constant goes, and g
        // comes from the mass and radius of the thing underneath. Gathered here
        // because the loop below holds `self.falling` mutably.
        let mut field: HashMap<usize, crate::math::Vec3> = HashMap::new();
        for node in nodes {
            let bodies = self.tree.refine(node).to_vec();
            if let Some(topo) = self.tree.nodes[node.get()].topology.clone() {
                struck.insert(node.get(), (bodies, topo));
            }
            field.insert(node.get(), self.tree.gravity_at(node));
        }

        // What else a piece falling out of each node could reach. `docs/PLAY.md`
        // Phase 1 asks for "a branch lands on the next tree", and until this a
        // piece could strike only the structure it came off: `struck` was keyed
        // by the falling piece's own node and nothing else was tested. Measured,
        // on two trees 8 m apart whose crowns overlap by two metres — radius
        // 10.4 and 10.5 — nine limbs came off one in a gale, made 225 contacts
        // and struck 160 members, every one of them in the tree they fell from.
        //
        // The candidates are the node's siblings, which is where D3's adjacency
        // index earns its place: the parent already knows what is next to what,
        // and a sibling with a topology is a thing that can be landed on. Each
        // carries the offset into the falling piece's frame, so the geometry is
        // compared in one frame rather than two.
        //
        // Index 0 is always the node itself, which is what makes the ground
        // test — a property of one structure — stay attached to the right one.
        let mut targets: HashMap<usize, Vec<(NodeIdx, crate::math::Vec3)>> = HashMap::new();
        for node in self.falling.iter().map(|(n, _)| *n).collect::<Vec<_>>() {
            if targets.contains_key(&node.get()) {
                continue;
            }
            let mut list = vec![(node, crate::math::Vec3::ZERO)];
            let parent = self.tree.nodes[node.get()].parent;
            if !parent.is_none() {
                let nb = self.tree.neighbourhood(parent);
                let reach = nb.reach();
                if let Some(near) = nb.near(self.tree.nodes[node.get()].motion.offset, reach) {
                    for occ in near {
                        let crate::neighbourhood::Occupant::Child(c) = occ else { continue };
                        if c == node || c.is_none() || !self.tree.nodes[c.get()].alive {
                            continue;
                        }
                        if self.tree.nodes[c.get()].topology.is_none()
                            && self.tree.nodes[c.get()].morphology.is_none()
                        {
                            continue;
                        }
                        // Where the neighbour sits, seen from the falling
                        // piece's frame.
                        let offset = self.tree.separation(node, crate::math::Vec3::ZERO, c, crate::math::Vec3::ZERO);
                        list.push((c, offset.value));
                    }
                }
            }
            for (c, _) in list.iter().skip(1).map(|(c, o)| (*c, *o)).collect::<Vec<_>>() {
                let bodies = self.tree.refine(c).to_vec();
                if let Some(topo) = self.tree.nodes[c.get()].topology.clone() {
                    struck.insert(c.get(), (bodies, topo));
                }
            }
            targets.insert(node.get(), list);
        }

        let mut strikes: Vec<(NodeIdx, u32, crate::math::Vec3)> = Vec::new();
        for (node, frag) in self.falling.iter_mut() {
            frag.age += dt;
            let n = frag.dynamics.dynamics.frame.joints.len();
            let mut load = vec![Dof::default(); n];
            // A piece falls in the field of whatever it is inside: Earth's
            // surface on Earth, a sixth of it on a moon, and essentially nothing
            // adrift between stars. It used to be 9.80665 m/s^2 down the
            // structure's own negative z, everywhere.
            let g = field.get(&node.get()).copied().unwrap_or(crate::math::Vec3::ZERO);
            for i in 0..n {
                let m = frag.dynamics.dynamics.frame.lumped[i].t.z;
                load[i].t = g.scale(m);
            }
            let rep = frag.dynamics.dynamics.step(&load, dt);
            report.broken_while_falling += rep.broken.len();

            let candidates = targets.get(&node.get()).cloned().unwrap_or_default();
            let mut contacts = Vec::new();
            for (k, (target, offset)) in candidates.iter().enumerate() {
                match struck.get(&target.get()) {
                    Some((bodies, topo)) => contacts.extend(frag.contacts_on(
                        bodies,
                        topo,
                        ground_of(topo),
                        k as u32,
                        *offset,
                    )),
                    None if k == 0 => contacts.extend(frag.contacts_on(
                        &[],
                        &crate::topology::Topology::default(),
                        0.0,
                        0,
                        crate::math::Vec3::ZERO,
                    )),
                    None => {}
                }
            }
            report.contacts += contacts.len();
            // Derived, not chosen. This was `0.15` with the comment "Wood on
            // wood: it does not bounce", which is true at the speed it was
            // tuned at and less true at every other one — the same two pieces
            // of timber return a third of the approach at a walking pace and a
            // seventh of it at twenty metres a second. `docs/PLAY.md` D3 is
            // partly about this constant.
            //
            // The speed is the fastest of the contacts rather than the mean:
            // `resolve` takes one restitution for the whole set, and it is the
            // hardest contact that decides whether the piece bounces or stays.
            let closing = contacts
                .iter()
                .map(|c| c.closing.abs())
                .fold(0.0f64, f64::max);
            let mine = crate::neighbourhood::Resilience::of(&frag.topo.material);
            let theirs = struck
                .get(&node.get())
                .map(|(_, t)| crate::neighbourhood::Resilience::of(&t.material))
                .unwrap_or(mine);
            let e = crate::neighbourhood::restitution(&mine, &theirs, closing);
            for (target, member, impulse) in frag.resolve(&contacts, e) {
                let hit = candidates
                    .get(target as usize)
                    .map(|(n, _)| *n)
                    .unwrap_or(*node);
                strikes.push((hit, member, impulse));
            }
        }

        // One analysis per node per step, not one per contact. Every impulse
        // that arrived in the same step arrived together, and asking the
        // structure to answer for them one at a time is both wrong — each
        // answer regenerates the structure and renumbers the members the next
        // impulse was aimed at — and unaffordable, at ninety contacts a step.
        //
        // An impulse is also not a force. The mechanism takes a force, and what
        // a struck member feels is the impulse spread over the step it arrived
        // in: at fifty steps a second that is fifty times what handing the
        // impulse over unconverted would suggest, and it is the difference
        // between a limb that breaks what it lands on and one that settles into
        // it without a sound.
        let mut by_node: HashMap<usize, (NodeIdx, Vec<st::Mechanism>)> = HashMap::new();
        for (node, member, impulse) in strikes {
            let entry = by_node.entry(node.get()).or_insert_with(|| (node, Vec::new()));
            entry.1.push(st::Mechanism::PointImpulse {
                at: member,
                impulse: impulse.scale(1.0 / dt),
            });
            report.largest_impulse = report.largest_impulse.max(impulse.norm());
            report.struck_members += 1;
        }
        let mut ordered: Vec<(usize, (NodeIdx, Vec<st::Mechanism>))> = by_node.into_iter().collect();
        ordered.sort_by_key(|(k, _)| *k);
        for (_, (node, mechanisms)) in ordered {
            let out = self.damage(node, &mechanisms);
            report.peak_utilisation = report.peak_utilisation.max(out.peak_utilisation);
            report.secondary_breaks += out.broken_joints;
            report.secondary_mass += out.detached_mass;
        }

        // Pieces that have landed and stopped stop being simulated: a branch
        // lying on the ground is litter, and the node already accounts for its
        // mass. Keeping it would spend a frame budget on something that has
        // finished happening.
        let before = self.falling.len();
        self.falling
            .retain(|(_, f)| !f.at_rest() && f.age < MAX_FALL_SECONDS);
        report.settled = before - self.falling.len();
        report.still_falling = self.falling.len();
        report
    }

    /// Stop integrating everything, releasing the dynamic state.
    pub fn stop_dynamics(&mut self) {
        self.shaking.clear();
        self.falling.clear();
    }

    /// The dynamic state of one node, if it is being integrated.
    pub fn shaken(&self, idx: NodeIdx) -> Option<&crate::solvers::structure::DynamicStructure> {
        self.shaking.iter().find(|(n, _)| *n == idx).map(|(_, ds)| ds)
    }

    /// Read the growth environment off a node and its surroundings.
    ///
    /// A scenario can override this per node — placing a lit planetary surface
    /// is authoring, not physics, and deriving insolation from the galaxy's
    /// own luminosity gives a correct answer (about 10^-4 W/m^2) that is
    /// correct precisely because a tree in interstellar space does not grow.
    pub fn environment_at(&self, idx: NodeIdx) -> crate::morph::Environment {
        let n = &self.tree.nodes[idx.get()];
        if let Some(env) = self.identities.get(&n.key).and_then(|id| self.environments.get(id)) {
            let mut env = *env;
            // A stated environment overrides what it states and no more. Saying
            // nothing about the feedstock is not saying "nothing", so the
            // measurement stands — which is what keeps a scenario that plants a
            // tree with `Environment::default()` from building it out of
            // hydrogen and helium.
            if env.feedstock.is_none() {
                env.feedstock = n.matter.composition;
            }
            return env;
        }
        // Illumination from the parent's luminosity at this node's distance —
        // so a structure in the shade of its own node's parent really is in the
        // shade, without a separate lighting system.
        let light = if !n.parent.is_none() {
            let p = &self.tree.nodes[n.parent.get()];
            let d = n.motion.offset.norm().max(p.matter.radius * 0.01).max(1e-6);
            (p.matter.luminosity / (4.0 * std::f64::consts::PI * d * d)).min(1400.0)
        } else {
            crate::morph::Environment::default().light_flux
        };
        // Water availability, measured rather than declared — and measured
        // without the engine ever being told which substance is water, which it
        // deliberately does not know. What a growing thing needs is a *mobile
        // solvent*, and that is the liquid phase, whatever it happens to be.
        //
        // This is what makes a desert a desert. Nothing anywhere tests for a
        // biome: a patch whose mixture holds no liquid grows nothing, because
        // `thermal_factor` and this between them leave no growth to be had. And
        // it is why a patch freezes into one: `react_all` moves mass from
        // liquid to solid as the temperature falls, so the same ground that
        // supported a forest in summer supports nothing in winter, through the
        // phase machinery that was written for salt dissolving in a beaker.
        //
        // A node with no mixture at all has not had its chemistry described, so
        // there is nothing to measure and the fallback is "unlimited". That is
        // the honest answer to no information, and it is what everything built
        // before mixtures existed relies on.
        let water = if n.matter.mixture.is_empty() {
            1.0
        } else {
            n.matter.mixture.in_phase(crate::chem::Phase::Liquid)
        };
        // Competition for the same ground. A node already mostly structure has
        // little room left, which is what stops a forest growing without bound
        // and what makes a clearing fill in faster than a thicket.
        let structural = n.morphology.as_ref().map(|m| m.built).unwrap_or(0.0);
        let crowding = if n.matter.mass > 0.0 {
            (structural / n.matter.mass).clamp(0.0, 1.0)
        } else {
            0.0
        };
        crate::morph::Environment {
            light_flux: light,
            temperature: n.matter.temperature,
            water,
            crowding,
            // A structure can only be built out of matter that is actually
            // available: the node's *unstructured* remainder, not its total.
            // Using the total lets a node grow a structure many times its own
            // mass, because the limit then applies per step rather than
            // cumulatively — a one-kilogram node grew a three-tonne tree.
            reservoir_mass: (n.matter.mass
                - n.morphology.as_ref().map(|m| m.built).unwrap_or(0.0))
            .max(0.0),
            labour: self.labour_rate,
            // What the node is standing in, measured: the mass of everything
            // that is not solid, over the volume it occupies. A node nobody has
            // described is in air, which is the honest answer to no
            // information and what everything built before mixtures existed
            // relies on.
            fluid_density: {
                let fluid = if n.matter.mixture.is_empty() {
                    0.0
                } else {
                    n.matter.mixture.in_phase(crate::chem::Phase::Liquid)
                        + n.matter.mixture.in_phase(crate::chem::Phase::Gas)
                };
                let volume = n.matter.volume();
                if fluid > 0.0 && volume > 0.0 {
                    (n.matter.mass * fluid / volume).max(1.225)
                } else {
                    1.225
                }
            },
            // The flow this node has met. Its own bulk motion relative to what
            // contains it is the wind it feels, and a node that has never moved
            // feels the one everything is proportioned against by default.
            flow_speed: {
                let moving = n.motion.velocity.norm();
                if moving > 0.0 {
                    moving
                } else {
                    crate::morph::Environment::default().flow_speed
                }
            },
            // What there is here to build out of. D11's `substrate` column,
            // measured: a structure is made of what its node holds.
            feedstock: n.matter.composition,
        }
    }

    fn deliver_influences(&mut self, until: f64) {
        let arrivals = self.mailbox.drain_until(until);
        for inf in arrivals {
            self.apply_influence(inf);
        }
    }

    fn apply_influence(&mut self, inf: Influence) {
        if inf.target.is_none() || !self.tree.nodes[inf.target.get()].alive {
            return;
        }
        let n = &mut self.tree.nodes[inf.target.get()];
        match inf.kind {
            InfluenceKind::Radiation | InfluenceKind::Blast => {
                n.matter.add_heat(inf.energy);
                n.matter.momentum += inf.momentum;
            }
            InfluenceKind::Impact | InfluenceKind::UserImpulse => {
                n.matter.momentum += inf.momentum;
                n.matter.add_heat(inf.energy);
            }
            InfluenceKind::Probe => {}
            InfluenceKind::Exchange => {
                n.matter.add_heat(inf.energy);
                n.matter.momentum += inf.momentum;
            }
        }
        // An exchange with a neighbour is ordinary physics between two things
        // the engine already knows, replayable from the same seeds, so it does
        // not pin. Everything else here is information from outside that no
        // amount of re-sampling would reproduce, and pinning is how the node
        // says so. Pinning on an exchange would pin every node with a warm
        // neighbour, and its whole ancestry, permanently — `Tree::pin` is
        // one-way — which is axiom four exactly inverted.
        let idx = inf.target;
        if inf.kind != InfluenceKind::Exchange {
            self.tree.pin(idx);
            self.disturb(idx);
        }
        if !self.tree.nodes[idx.get()].bodies.is_empty() {
            // Distribute the impulse over the existing bodies rather than
            // discarding them — throwing away detail a user is looking at, in
            // response to that user poking it, is the worst possible moment.
            let n = &mut self.tree.nodes[idx.get()];
            let total = n.matter.mass.max(1e-300);
            for b in n.bodies.iter_mut() {
                let f = b.mass / total;
                if b.mass > 0.0 {
                    b.vel += inf.momentum.scale(f / b.mass);
                }
                // `add_heat` rather than a bare `internal_energy +=`: the
                // matter's temperature moved on the line above, and a body list
                // whose temperatures did not follow is a materialisation that
                // no longer summarises to its own matter.
                b.add_heat(inf.energy * f);
            }
        }
    }

    fn record_histories(&mut self) {
        let depth = self.history_depth;
        let entries: Vec<(PathKey, Moment)> = self
            .tree
            .nodes
            .iter()
            .filter(|n| n.alive && n.residency.rank() >= Residency::Causal.rank())
            .map(|n| {
                (
                    n.key,
                    Moment {
                        t: self.time,
                        offset: n.motion.offset,
                        velocity: n.motion.velocity,
                        mass: n.matter.mass,
                        luminosity: n.matter.luminosity,
                        temperature: n.matter.temperature,
                    },
                )
            })
            .collect();
        for (key, snap) in entries {
            self.histories
                .entry(key)
                .or_insert_with(|| History::new(depth))
                .push(snap);
        }
    }

    // -----------------------------------------------------------------
    // interaction
    // -----------------------------------------------------------------

    /// Apply a user interaction. Everything the user can do goes through here,
    /// so there is exactly one place where the world can change for
    /// non-physical reasons — and it is audited.
    pub fn interact(&mut self, action: Interaction) {
        match action {
            Interaction::Impulse { target, dp } => {
                let d = self.observer_distance(target);
                self.mailbox.post(
                    target,
                    self.time,
                    d,
                    InfluenceKind::UserImpulse,
                    0.0,
                    dp,
                );
            }
            Interaction::Deposit {
                target,
                joules,
                radius: _,
            } => {
                let d = self.observer_distance(target);
                self.mailbox
                    .post(target, self.time, d, InfluenceKind::Radiation, joules, Vec3::ZERO);
            }
            Interaction::Extract { target, joules } => {
                let d = self.observer_distance(target);
                self.mailbox.post(
                    target,
                    self.time,
                    d,
                    InfluenceKind::Radiation,
                    -joules,
                    Vec3::ZERO,
                );
            }
            Interaction::Inject {
                target,
                mass,
                composition,
                velocity,
            } => {
                if target.is_none() || !self.tree.nodes[target.get()].alive {
                    return;
                }
                let n = &mut self.tree.nodes[target.get()];
                let old = n.matter.mass;
                n.matter.composition =
                    crate::state::Composition::blend(n.matter.composition, old, composition, mass);
                n.matter.mass += mass;
                n.matter.momentum += velocity.scale(mass);
                n.matter.baryon_number = n.matter.mass * n.matter.composition.nucleons_per_kg();
                self.tree.pin(target);
                self.tree.bump_epoch(target);
                self.disturb(target);
            }
            Interaction::Measure {
                target,
                instrument,
                quantity,
            } => {
                self.measure(target, instrument, quantity);
            }
            Interaction::Pin { target } => {
                self.tree.pin(target);
                self.disturb(target);
            }
            Interaction::Author {
                target,
                property,
                value,
            } => self.author(target, property, value),
            // Immediate, unlike `Impulse` and `Deposit` above, which are posted
            // to the mailbox and arrive at the speed of light. A bubble is not
            // an influence travelling through the world to reach the node — it
            // is a change to how fast the engine agrees to run it, made by
            // somebody standing outside the simulation. There is no distance
            // for it to cross. Same reason `Pin` and `Author` apply here.
            Interaction::Dilate { target, rate } => {
                self.dilate(target, rate);
            }
            Interaction::Mark { target, deviation } => {
                self.mark(target, deviation);
            }
        }
    }

    /// Put a deviation into a node's field.
    ///
    /// `docs/PLAY.md` §5.1: a mark is not a pin. The old rule was that anything
    /// touched is persisted outright and for ever, which at play resolution
    /// makes a footprint as permanent as a felled trunk. A mark is remembered
    /// exactly while it is there and decays at the rate its own material, flux
    /// and geometry set — see [`World::weather`].
    ///
    /// The mass it carried away leaves the node here, at the moment the mark is
    /// made, because that is when it left: a cut bank is lighter afterwards.
    /// Nothing debits it later and nothing has to remember to, which is the
    /// property §5.8 asks for.
    pub fn mark(&mut self, target: NodeIdx, deviation: crate::erode::Deviation) -> bool {
        if target.is_none() || !self.tree.nodes[target.get()].alive {
            return false;
        }
        if self.tree.nodes[target.get()].morphology.is_none() {
            return false;
        }
        let moved = deviation.moved.abs();
        {
            let n = &mut self.tree.nodes[target.get()];
            if let Some(m) = n.morphology.as_mut() {
                m.field.push(deviation);
                m.built = (m.built - moved).max(0.0);
            }
            n.matter.mass = (n.matter.mass - moved).max(0.0);
        }
        self.disturb(target);
        true
    }

    /// Let every stored deviation feel its own physics for `dt` seconds.
    ///
    /// `docs/PLAY.md` §5, once a frame. Three things happen to a deviation and
    /// the node decides none of them:
    ///
    /// * its amplitude falls at [`crate::erode::relaxation_rate`], derived from
    ///   the material's cohesion, the flux the node is actually in, and the
    ///   feature's own half-width;
    /// * it is **dropped** when the amplitude is below what the coarse level can
    ///   represent — which is §5.1's "forgetting becomes the deviation reaching
    ///   zero" — but only if dropping it leaves the conserved tuple unchanged;
    /// * it is **kept for ever** otherwise, because §5.8 forbids forgetting to
    ///   be a source: a deviation that took mass away would put it back.
    ///
    /// **Decay runs on the world clock, not on the absence of an audience**
    /// (§5.5). What observation changes is only whether the detail has to be
    /// stored, and a node currently resolved for somebody keeps its field
    /// whatever the amplitude has fallen to, because dropping it in front of
    /// them is the pop §5.4 forbids.
    fn weather(&mut self, dt: f64) {
        if !(dt > 0.0) {
            return;
        }
        for i in 0..self.tree.nodes.len() {
            let idx = NodeIdx(i as u32);
            {
                let n = &self.tree.nodes[i];
                if !n.alive {
                    continue;
                }
                match n.morphology.as_ref() {
                    Some(m) if !m.field.is_empty() => {}
                    _ => continue,
                }
            }
            let env = self.environment_at(idx);
            let (temperature, mass, watched) = {
                let n = &self.tree.nodes[i];
                (
                    n.matter.temperature,
                    n.matter.mass,
                    n.pinned || n.residency.rank() >= Residency::Observed.rank(),
                )
            };
            // What it takes to detach one grain of what this is made of, which
            // is the material's strength at the size of its own worst flaw
            // rather than the whole piece's. Measured off the node's mixture
            // where it has one; the program's material otherwise.
            let material = match self.material_of(idx) {
                Some(m) => m,
                None => match self.tree.nodes[i].morphology.as_ref() {
                    Some(m) => m.program.material(),
                    None => continue,
                },
            };
            let cohesion = match self.loose_grain_strength(idx, &material, env.fluid_density) {
                Some(loose) => loose,
                None => material.strength_of(material.flaw_size, temperature),
            };
            let local = dt * self.local_rate(idx);
            let mut dropped = 0u64;
            let n = &mut self.tree.nodes[i];
            let Some(m) = n.morphology.as_mut() else { continue };
            for d in m.field.iter_mut() {
                let rate =
                    crate::erode::relaxation_rate(cohesion, env.fluid_density, env.flow_speed, d.span);
                if rate > 0.0 && d.conservative(mass) {
                    d.amplitude *= (-rate * local).exp();
                }
            }
            if !watched {
                let floor = crate::tree::IDEMPOTENT_TOLERANCE * n.matter.radius.abs().max(1e-300);
                let before = m.field.len();
                m.field.retain(|d| !(d.conservative(mass) && d.amplitude.abs() <= floor));
                dropped = (before - m.field.len()) as u64;
            }
            if dropped > 0 {
                self.tree.stats.deviations_forgotten += dropped;
                self.disturb(idx);
            }
        }
    }

    /// Give a node the ocean its own matter describes, if it describes one and
    /// has none yet. Returns whether it has one afterwards.
    ///
    /// **Measured, not stated.** A planetary node whose mixture holds a liquid
    /// on a solid majority has water on rock, and that water spreads over the
    /// node's surface to the depth its own volume gives — one depth everywhere
    /// while ground on a sphere has no relief (`ocean.rs`). Its surface gravity
    /// is its own, and its seabed drag is the log-law against the ground's own
    /// grain, the deposited flaw scale `Material::measured` derives.
    pub fn assess_ocean(&mut self, idx: NodeIdx) -> bool {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return false;
        }
        if self.tree.nodes[idx.get()].ocean.is_some() {
            return true;
        }
        let (liquid, solid, mass, radius, tier) = {
            let n = &self.tree.nodes[idx.get()];
            let m = &n.matter.mixture;
            (m.in_phase(crate::chem::Phase::Liquid), m.in_phase(crate::chem::Phase::Solid), n.matter.mass, n.matter.radius, n.tier)
        };
        if tier != crate::units::Tier::Planetary || !(liquid > 0.0) || solid < 0.5 || !(radius > 0.0) {
            return false;
        }
        // The liquid's volume, pool by pool at its own rest density, and the
        // surface tension it has, by mass.
        let mut volume = 0.0;
        let (mut liquid_mass, mut tension) = (0.0, 0.0);
        for p in self.tree.nodes[idx.get()].matter.mixture.entries() {
            if p.phase != crate::chem::Phase::Liquid {
                continue;
            }
            let Some(sub) = self.substances.get(p.substance) else { continue };
            let Some(c) = crate::eos::Condensed::liquid(&sub.props) else { continue };
            let m = p.fraction * mass;
            volume += m / c.rest_density;
            liquid_mass += m;
            tension += m * crate::erode::surface_tension(&sub.props);
        }
        let density = if volume > 0.0 { liquid_mass / volume } else { 0.0 };
        let tension = if liquid_mass > 0.0 { tension / liquid_mass } else { 0.0 };
        let area = 4.0 * std::f64::consts::PI * radius * radius;
        let depth = volume / area;
        if !(depth > 0.0) {
            return false;
        }
        let g = crate::units::G * mass / (radius * radius);
        let grain = self.seabed_grain(idx);
        let drag = crate::ocean::log_law_drag(depth, grain);
        let mut ocean = crate::ocean::Ocean::new(OCEAN_CELLS, radius, depth, g, drag, density, tension);
        ocean.time = self.tree.nodes[idx.get()].time;
        self.tree.nodes[idx.get()].ocean = Some(Box::new(ocean));
        true
    }

    /// The grain an ocean's bed is rough with, metres.
    ///
    /// **The grain the ground's own freezing laid down, where it froze** — the
    /// owner's decision for Phase 5. A seabed is rough with what it is made
    /// of, and what a melt is made of once it has frozen is the grain
    /// `derive_frozen_layout` wrote into its `Recipe::Granular`, found on the
    /// node or anywhere up its ancestry the way `is_loose` finds it. Where no
    /// freezing event ever happened there is no grain to read, and the flaw
    /// scale `Material::measured` derives stands in: for a planet nobody saw
    /// freeze that is the deposited increment, 3.9 km on an Earth, which
    /// describes how the planet was laid down rather than its bed. Measured on
    /// that Earth: a drag of 2.57e-2 from the flaw against 7.56e-4 from the
    /// 1.66 cm grain a silicate melt froze to in `tests/accretion.rs`.
    pub fn seabed_grain(&self, idx: NodeIdx) -> f64 {
        let mut at = idx;
        while !at.is_none() {
            if let Some(crate::recipe::Recipe::Granular(g)) =
                self.tree.nodes[at.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref())
            {
                if g.grain > 0.0 && g.grain.is_finite() {
                    return g.grain;
                }
            }
            at = self.tree.nodes[at.get()].parent;
        }
        self.material_of(idx).map(|m| m.flaw_size).unwrap_or(0.0)
    }

    /// The equilibrium elevation of each cell of a node's ocean — the tidal
    /// potential over `-g` — from everything else the node's parent holds, at
    /// where it is now, in the node's own axes.
    ///
    /// Its parent's bodies where the parent is materialised, the planet's own
    /// stand-in left out; the parent's mass at its own centre where it is not.
    /// Exact, not the quadrupole: `ocean::tidal_potential` subtracts only what
    /// moves the whole planet.
    pub fn tidal_equilibrium(&self, idx: NodeIdx) -> Vec<f64> {
        let n = &self.tree.nodes[idx.get()];
        let Some(ocean) = n.ocean.as_ref() else { return Vec::new() };
        let parent = n.parent;
        let mut sources: Vec<(Vec3, f64)> = Vec::new();
        // Offsets are in root-aligned axes, and what takes a vector derived
        // from them into the planet's frame is the whole composition to the
        // root (`Tree::gravity_at`), not the planet's own facing alone.
        let into = self.tree.axes_from(self.tree.root, idx).conjugate();
        if !parent.is_none() {
            let p = &self.tree.nodes[parent.get()];
            if p.bodies.is_empty() {
                let m = (p.matter.mass - n.matter.mass).max(0.0);
                sources.push((into.rotate(Vec3::ZERO - n.motion.offset), m));
            } else {
                // Where each body is at the planet's instant. The parent's
                // bodies are at the time its contents were solved to, and the
                // planet is carried to the world's (`Node::carried`); between
                // the two, a body goes at its own velocity. Reading it where
                // the last solve left it held a moon still — an Earth-moon
                // pair's cadence is 3.7 days, and the tide it drove came out at
                // half a sidereal day rather than half a lunar one.
                let since = n.carried - p.time;
                for (k, b) in p.bodies.iter().enumerate() {
                    if k == n.slot as usize {
                        continue;
                    }
                    let at = b.pos + b.vel.scale(since);
                    sources.push((into.rotate(at - n.motion.offset), b.mass));
                }
            }
        }
        ocean
            .cells
            .iter()
            .map(|c| {
                let x = c.up.scale(ocean.radius);
                let phi: f64 = sources.iter().map(|(r, m)| crate::ocean::tidal_potential(x, *r, *m)).sum();
                -phi / ocean.g
            })
            .collect()
    }

    /// Carry every ocean in the world forward by a span of world time, on its
    /// own node's clock and at its own stable step.
    ///
    /// Cheap by construction: an Earth's ocean at sixteen cells a face side is
    /// 1536 cells stepped every thousand seconds or so, which at one second per
    /// second is once in twenty thousand frames. It runs whether or not anyone
    /// is watching, like weather, because a tide nobody watched still came in.
    fn tides(&mut self, span: f64) {
        if !(span > 0.0) {
            return;
        }
        for i in 0..self.tree.nodes.len() {
            let idx = NodeIdx(i as u32);
            if !self.tree.nodes[i].alive || self.tree.nodes[i].tier != crate::units::Tier::Planetary {
                continue;
            }
            if self.tree.nodes[i].ocean.is_none() && !self.assess_ocean(idx) {
                continue;
            }
            let local = span * self.local_rate(idx);
            let spin = self.tree.nodes[i].motion.spin_rate;
            let stable = self.tree.nodes[i].ocean.as_ref().map(|o| o.stable_step()).unwrap_or(0.0);
            if !(stable > 0.0) {
                continue;
            }
            // The forcing moves with the moon and the planet's turning, both
            // of which the frame has already advanced; within a span it is
            // held, which is exact for the forcing's own period being far
            // longer than a stable step.
            let eq = self.tidal_equilibrium(idx);
            let Some(ocean) = self.tree.nodes[i].ocean.as_mut() else { continue };
            // Whole stable steps, and the remainder carried: a frame is far
            // shorter than a step, so most frames take none.
            let owed = local + ocean.owed;
            let steps = (owed / stable).floor() as u64;
            for _ in 0..steps {
                ocean.step(stable, &eq, spin);
            }
            ocean.owed = owed - steps as f64 * stable;
            // The wind works every frame rather than every stable step: its
            // stress is a rate, and a sea at the beach grows over seconds.
            self.blow(idx, local);
        }
    }

    /// The winds over a planet's ocean: every gas node it holds whose air
    /// reaches the sea, as a wind over each cell its footprint covers, in the
    /// planet's own axes, with the node each came from.
    ///
    /// **The air's own motion is the wind**, relative to the water under it:
    /// its velocity less the surface's turning and the current there. A node's
    /// velocity is its centre's, so that is the height the wind is measured
    /// at — the ten metres a scene states is the meteorologist's `U10` — and a
    /// node whose centre is not above the sea is not air over it. What it
    /// blows over is the disc its sphere cuts from the sea surface. Its
    /// density is its own.
    ///
    /// The first rule written was the mean height of the column over the
    /// footprint's centre, `(h + r) / 2`, which reads a thousand-kilometre
    /// mass of air's wind at five hundred kilometres — by the log-law over a
    /// calm sea, a drag coefficient of about 1e-4 for a wind the scene had
    /// stated at ten metres.
    ///
    /// **Where the air comes from is the scene's for now** — the owner's
    /// decision for Phase 5. `docs/PLAY.md` derives prevailing wind from a
    /// planet's rotation and its equator-to-pole gradient, and no phase has
    /// specified that yet; until the root `Program` carries a climate, a scene
    /// that wants a wind states a node of air and pushes it.
    pub fn winds_over(&self, idx: NodeIdx) -> Vec<(crate::ocean::Wind, NodeIdx)> {
        let mut out = Vec::new();
        let n = &self.tree.nodes[idx.get()];
        let Some(ocean) = n.ocean.as_ref() else { return out };
        let into = self.tree.axes_from(self.tree.root, idx).conjugate();
        let spin = n.motion.spin_rate;
        // Air anywhere under the planet: over open sea it is the planet's own
        // child, and near a shore it is held by the patch of ground it is over.
        let mut under: Vec<NodeIdx> = Vec::new();
        let mut stack: Vec<NodeIdx> = n.children.iter().copied().filter(|c| !c.is_none()).collect();
        while let Some(c) = stack.pop() {
            if !self.tree.nodes[c.get()].alive {
                continue;
            }
            under.push(c);
            stack.extend(self.tree.nodes[c.get()].children.iter().copied().filter(|c| !c.is_none()));
        }
        for c in under {
            let air = &self.tree.nodes[c.get()];
            if air.matter.mixture.is_empty() || air.matter.mixture.in_phase(crate::chem::Phase::Gas) < 0.5 {
                continue;
            }
            let r = air.matter.radius;
            let volume = air.matter.volume();
            if !(r > 0.0) || !(volume > 0.0) {
                continue;
            }
            let at = into.rotate(self.tree.offset_from(idx, c, Vec3::ZERO).value);
            let h = at.norm() - ocean.radius;
            if !(h > 0.0) || !(h < r) {
                continue;
            }
            let footprint = (r * r - h * h).sqrt();
            let height = h;
            let velocity = into.rotate(self.tree.velocity_from(idx, c));
            let density = air.matter.mass / volume;
            let reach = (footprint / ocean.radius).min(std::f64::consts::PI).cos();
            let dir = at.unit();
            let mut covered: Vec<(usize, f64)> = ocean
                .cells
                .iter()
                .enumerate()
                .filter(|(_, cell)| cell.up.dot(dir) >= reach)
                .map(|(k, cell)| (k, cell.area))
                .collect();
            if covered.is_empty() {
                // Smaller than a cell: the one it stands over, over its own
                // footprint.
                let k = ocean.cell_of(dir);
                let disc = std::f64::consts::PI * footprint * footprint;
                covered.push((k, disc.min(ocean.cells[k].area)));
            }
            for (k, area) in covered {
                let cell = &ocean.cells[k];
                let surface = spin.cross(cell.up.scale(ocean.radius)) + cell.velocity;
                let rel = velocity - surface;
                let rel = rel - cell.up.scale(rel.dot(cell.up));
                out.push((
                    crate::ocean::Wind { cell: k, velocity: rel, height, air_density: density, area },
                    c,
                ));
            }
        }
        out
    }

    /// Let the winds over a planet's ocean work on it for `span` s of the
    /// planet's time, and carry its sea.
    ///
    /// **The books close the way a contact's do** (`apply_contact`): the air
    /// loses the momentum its stress carries, the planet's own contents take
    /// it — its bodies where it is materialised, its matter where it is not —
    /// and the kinetic energy the air gave up, less what those contents took
    /// up as motion, is heat in them. The ocean's currents and sea are a
    /// description of part of that energy, not a second account of it.
    fn blow(&mut self, idx: NodeIdx, span: f64) {
        if !(span > 0.0) {
            return;
        }
        let winds = self.winds_over(idx);
        if !winds.is_empty() {
            let list: Vec<crate::ocean::Wind> = winds.iter().map(|w| w.0).collect();
            let blown = match self.tree.nodes[idx.get()].ocean.as_mut() {
                Some(o) => o.raise(span, &list),
                None => return,
            };
            let to_root = self.tree.axes_from(self.tree.root, idx);
            let radius = self.tree.nodes[idx.get()].ocean.as_ref().map(|o| o.radius).unwrap_or(0.0);
            for ((w, air), b) in winds.iter().zip(blown.iter()) {
                let dp = to_root.rotate(b.momentum);
                if !(dp.norm() > 0.0) {
                    continue;
                }
                let where_ = to_root.rotate(self.tree.nodes[idx.get()].ocean.as_ref().unwrap().cells[w.cell].up.scale(radius));
                // The air.
                let a = &mut self.tree.nodes[air.get()];
                let m = a.matter.mass.max(1e-300);
                let v0 = a.motion.velocity;
                let v1 = v0 - dp.scale(1.0 / m);
                a.motion.velocity = v1;
                let lost = 0.5 * m * (v0.norm2() - v1.norm2());
                // The planet's contents.
                let taken = self.take_reaction(idx, dp, where_);
                self.add_heat_to(idx, lost - taken);
            }
        }
        // The sea travels, in steps it is stable at.
        let Some(ocean) = self.tree.nodes[idx.get()].ocean.as_mut() else { return };
        let mut left = span;
        while left > 0.0 {
            let step = ocean.wave_step().min(left);
            if !(step > 0.0) || !step.is_finite() {
                break;
            }
            ocean.propagate(step);
            left -= step;
        }
    }

    /// The fluid a body holds around itself, as it is described rather than
    /// drawn: `(base, surface density, scale height)` of its atmosphere, m,
    /// kg/m^3 and m, or `None` for a node that has none.
    ///
    /// **Derived from its own gas pool**, the way its ocean is from its liquid
    /// one. The gas's weight over the surface is the pressure at its base,
    /// `p0 = M_gas g / 4 pi R^2`; an isothermal column at the node's own
    /// temperature falls off over `H = k T / m g` with `m` the gas's mean
    /// molecular mass; and the base density is then `p0 / g H`. The base is
    /// the sea's surface where there is one and the ground's where there is
    /// not.
    ///
    /// Only for a body whose surface has been assessed — an ocean, or ground a
    /// rounded body carries. A node whose bodies *are* its fluid, a bucket of
    /// water or a cloud, already pushes on what is in it through its solver,
    /// and a second account of the same pressure would double it.
    pub fn atmosphere_of(&self, idx: NodeIdx) -> Option<(f64, f64, f64)> {
        let n = &self.tree.nodes[idx.get()];
        let base = match (n.ocean.as_ref(), n.morphology.as_ref().and_then(|m| m.recipe.as_ref())) {
            (Some(o), _) => o.radius,
            (None, Some(crate::recipe::Recipe::Tiled(t))) if t.is_ball() => t.sphere,
            _ => return None,
        };
        let (mut gas, mut per_molecule) = (0.0, 0.0);
        for p in n.matter.mixture.entries() {
            if p.phase != crate::chem::Phase::Gas {
                continue;
            }
            let Some(sub) = self.substances.get(p.substance) else { continue };
            if !(sub.props.unit_mass > 0.0) {
                continue;
            }
            gas += p.fraction;
            per_molecule += p.fraction / sub.props.unit_mass;
        }
        if !(gas > 0.0) || !(per_molecule > 0.0) || !(base > 0.0) {
            return None;
        }
        let molecule = gas / per_molecule;
        let mass = gas * n.matter.mass;
        let g = crate::units::G * n.matter.mass / (base * base);
        let height = crate::units::K_B * n.matter.temperature.max(1.0) / (molecule * g);
        let density = mass / (4.0 * std::f64::consts::PI * base * base * height);
        (height > 0.0 && density > 0.0).then_some((base, density, height))
    }

    /// The density of the fluid a node describes around itself at `r` from
    /// its centre, kg/m^3: its sea below the sea's surface, its atmosphere
    /// above it, and nothing inside the ground or where it describes neither.
    pub fn ambient_density(&self, idx: NodeIdx, r: f64) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        if let Some(o) = n.ocean.as_ref() {
            if r < o.radius {
                return if r > o.radius - o.depth { o.density } else { 0.0 };
            }
        }
        match self.atmosphere_of(idx) {
            Some((base, rho, height)) if r >= base => rho * (-(r - base) / height).exp(),
            _ => 0.0,
        }
    }

    /// The radius of the surface a node describes — its sea's where it has
    /// one, its ground's where it is a rounded body — or `None`.
    pub fn surface_radius(&self, idx: NodeIdx) -> Option<f64> {
        let n = &self.tree.nodes[idx.get()];
        if let Some(o) = n.ocean.as_ref() {
            return Some(o.radius);
        }
        match n.morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
            Some(crate::recipe::Recipe::Tiled(t)) if t.is_ball() => Some(t.sphere),
            _ => None,
        }
    }

    /// The nearest of a node and its ancestors that describes a surface, and
    /// so the fluid around it: a sea, or a rounded body's ground and air.
    pub fn describing(&self, mut idx: NodeIdx) -> Option<NodeIdx> {
        while !idx.is_none() {
            if self.surface_radius(idx).is_some() {
                return Some(idx);
            }
            idx = self.tree.nodes[idx.get()].parent;
        }
        None
    }

    /// The density of the fluid around a child of `parent`, kg/m^3, from the
    /// nearest ancestor that describes one, at the child's distance from that
    /// ancestor's centre.
    ///
    /// **Not where the parent's own bodies are that fluid.** A patch of shore
    /// whose water is drawn as parcels already pushes on what is in it through
    /// its solver; the sea it describes is left out of the account below the
    /// sea's surface there, and the air above it is not.
    pub fn ambient_around(&self, parent: NodeIdx, child: NodeIdx) -> f64 {
        let Some(anc) = self.describing(parent) else { return 0.0 };
        let r = self.tree.offset_from(anc, child, Vec3::ZERO).value.norm();
        if anc != parent {
            let p = &self.tree.nodes[parent.get()];
            let resolves = p.is_materialised() && p.matter.mixture.in_phase(crate::chem::Phase::Liquid) > 0.0;
            let below = self.tree.nodes[anc.get()].ocean.as_ref().map(|o| r < o.radius).unwrap_or(false);
            if resolves && below {
                return 0.0;
            }
        }
        self.ambient_density(anc, r)
    }

    /// Whether a node holds a child, so that the child is carried in the
    /// turning frame its parent is in (`Node::turning`).
    ///
    /// **Held is support, measured**: a child is held where it touches or is
    /// inside what holds it, or where the fluid around it is at least as dense
    /// as it is, so the fluid floats it. What holds is a body that describes a
    /// surface — the ground, a patch of it, the sea and what is in it are all
    /// inside that — or anything such a body holds in turn, which is how a
    /// thing standing on a patch of ground goes round with the planet the
    /// patch is part of. A parcel of air at the density around it is held; a
    /// stone resting on the ground is held; a moon, or anything denser than
    /// the thin air it is passing through, is not, and is carried in the
    /// straight line a thing no force holds follows.
    pub fn holds(&self, parent: NodeIdx, child: NodeIdx) -> bool {
        let c = &self.tree.nodes[child.get()];
        let r = c.motion.offset.norm();
        let reach = match self.surface_radius(parent) {
            Some(surface) => surface,
            None if self.tree.nodes[parent.get()].turning != Vec3::ZERO => {
                self.tree.nodes[parent.get()].matter.radius
            }
            None => return false,
        };
        if r - c.matter.radius <= reach {
            return true;
        }
        let volume = c.matter.volume();
        volume > 0.0 && self.ambient_around(parent, child) >= c.matter.mass / volume
    }

    /// A node's own spin, in root-aligned axes — the axes offsets are in. A
    /// node's spin rate is kept in its parent's.
    fn spin_in_root(&self, idx: NodeIdx) -> Vec3 {
        let n = &self.tree.nodes[idx.get()];
        if n.parent.is_none() {
            n.motion.spin_rate
        } else {
            self.tree.axes_from(self.tree.root, n.parent).rotate(n.motion.spin_rate)
        }
    }

    /// Derive, for every node, the frame it is carried in until its parent
    /// next solves it, from the root down: where its parent holds it, the
    /// frame the parent is itself carried in if something holds the parent,
    /// and the parent's own spin if nothing does — a patch of ground goes round
    /// with its planet rather than by any spin of its own, and so does what
    /// stands on it. Measured with the other order: a face patch of an Earth
    /// carried a stray spin of 1e-37 rad/s of its own, and the air it held
    /// was carried at that. Zero where it is not held. Once a frame, before
    /// anything is solved.
    fn hold(&mut self) {
        let mut stack = vec![self.tree.root];
        while let Some(p) = stack.pop() {
            if p.is_none() || !self.tree.nodes[p.get()].alive {
                continue;
            }
            // A held parent goes round in the frame it is carried in, whatever
            // it is doing itself; only a parent nothing holds — a planet, a
            // moon — turns by its own spin.
            let carried_in = self.tree.nodes[p.get()].turning;
            let frame = if carried_in != Vec3::ZERO { carried_in } else { self.spin_in_root(p) };
            let children: Vec<NodeIdx> =
                self.tree.nodes[p.get()].children.iter().copied().filter(|c| !c.is_none()).collect();
            for c in children {
                if !self.tree.nodes[c.get()].alive {
                    continue;
                }
                let turning = if frame != Vec3::ZERO && self.holds(p, c) { frame } else { Vec3::ZERO };
                self.tree.nodes[c.get()].turning = turning;
                stack.push(c);
            }
        }
    }

    /// Archimedes, on the promoted children of `idx`: the weight of the fluid
    /// each displaces, pushing it up, with the reaction on the parent's own
    /// contents and the work booked as their heat, so the totals do not move.
    ///
    /// **Against the pull each child actually received**, `pulls` — the change
    /// the parent's solve made to its stand-in's velocity — rather than
    /// against the field `Tree::gravity_at` would give. The two differ wherever
    /// the parent's mass is drawn coarsely: a planet of nine parcels pulls a
    /// thing at its surface by whatever those nine happen to add up to, not by
    /// what a smooth sphere would. The fluid's weight is
    /// `rho V` times the same pull, so a thing at the density around it feels
    /// exactly nothing however lumpy the pull is. The gate in
    /// [`World::atmosphere_of`] is what keeps the pull a gravitational one:
    /// where the parent's own bodies are its fluid, the solver already pushes.
    ///
    /// **The owner's decision for Phase 5**, and it is what holds air up: a
    /// mass of air stated over a sea had nothing under it, and measured, a day
    /// later it had fallen through the planet and was moving at 882 m/s. At the
    /// density around it a parcel of air weighs nothing, a stone sinks and a
    /// ball of wood floats — one law, not a rule for air.
    ///
    /// The fluid goes round with the body that describes it, and so does what
    /// it holds: that is the turning carry's account (`Node::turning`), which
    /// is the exact solution for a thing whose forces supply its centripetal,
    /// so the fluid's own acceleration is not counted a second time here.
    fn buoy_children(&mut self, idx: NodeIdx, pulls: &[(NodeIdx, Vec3)]) {
        if self.describing(idx).is_none() {
            return;
        }
        for &(c, pull) in pulls {
            if !self.tree.nodes[c.get()].alive || !pull.is_finite() {
                continue;
            }
            let at = self.tree.position_at(c, self.tree.nodes[idx.get()].time);
            let rho = self.ambient_around(idx, c);
            let volume = self.tree.nodes[c.get()].matter.volume();
            if !(rho > 0.0) || !(volume > 0.0) {
                continue;
            }
            let dp = pull.scale(-rho * volume);
            let n = &mut self.tree.nodes[c.get()];
            let m = n.matter.mass.max(1e-300);
            let v0 = n.motion.velocity;
            n.motion.velocity = v0 + dp.scale(1.0 / m);
            let gained = 0.5 * m * (n.motion.velocity.norm2() - v0.norm2());
            let taken = self.take_reaction(idx, Vec3::ZERO - dp, at);
            self.add_heat_to(idx, -(gained + taken));
        }
    }

    /// Give a node's own contents an impulse at a place in root-aligned axes
    /// relative to its centre, and return the kinetic energy they took up, J.
    fn take_reaction(&mut self, idx: NodeIdx, dp: Vec3, at: Vec3) -> f64 {
        let n = &mut self.tree.nodes[idx.get()];
        let free: Vec<usize> = (0..n.bodies.len())
            .filter(|&s| {
                let c = n.child_of(s);
                c.is_none()
            })
            .collect();
        let mass: f64 = free.iter().map(|&s| n.bodies[s].mass).sum();
        if mass > 0.0 {
            let dv = dp.scale(1.0 / mass);
            let mut gained = 0.0;
            for &s in &free {
                let b = &mut n.bodies[s];
                let before = 0.5 * b.mass * b.vel.norm2();
                b.vel += dv;
                gained += 0.5 * b.mass * b.vel.norm2() - before;
            }
            return gained;
        }
        let m = n.matter.mass.max(1e-300);
        let before = 0.5 * n.matter.momentum.norm2() / m;
        n.matter.momentum += dp;
        n.matter.spin += at.cross(dp);
        0.5 * n.matter.momentum.norm2() / m - before
    }

    /// Heat into a node's own contents, J: its free bodies by mass where it is
    /// materialised, its matter where it is not.
    fn add_heat_to(&mut self, idx: NodeIdx, joules: f64) {
        if joules == 0.0 {
            return;
        }
        let n = &mut self.tree.nodes[idx.get()];
        let free: Vec<usize> = (0..n.bodies.len()).filter(|&s| n.child_of(s).is_none()).collect();
        let mass: f64 = free.iter().map(|&s| n.bodies[s].mass).sum();
        if mass > 0.0 {
            for &s in &free {
                let share = joules * n.bodies[s].mass / mass;
                n.bodies[s].add_heat(share);
            }
        } else {
            n.matter.internal_energy += joules;
        }
    }

    /// Whether a node's ground is an uncemented aggregate: laid down, and never
    /// seen to freeze.
    ///
    /// **History decides**, which is Phase 5's choice for what tells a poured
    /// pile from a cemented one when nothing in the engine states the neck area
    /// at a grain contact. A layout written by a freezing event
    /// (`Recipe::Granular`, see `World::assess_surface`) is a rock whose grains
    /// grew into each other, and it keeps the strength Griffith gives it. Ground
    /// — a `Recipe::Tiled` surface — with no such event anywhere up its
    /// ancestry was laid down and never bonded, and what holds a grain of it is
    /// its own weight. Grown and built things are not ground: their joints are
    /// explicit and carry their own bonds.
    pub fn is_loose(&self, idx: NodeIdx) -> bool {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return false;
        }
        let tiled = matches!(
            self.tree.nodes[idx.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()),
            Some(crate::recipe::Recipe::Tiled(_))
        );
        if !tiled {
            return false;
        }
        let mut at = idx;
        while !at.is_none() {
            if matches!(
                self.tree.nodes[at.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref()),
                Some(crate::recipe::Recipe::Granular(_))
            ) {
                return false;
            }
            at = self.tree.nodes[at.get()].parent;
        }
        true
    }

    /// The stress it takes to lift a grain of a loose node's ground, Pa, or
    /// `None` for ground that is not loose (`World::is_loose`).
    ///
    /// Its weight in the fluid over it, tipped out of its pocket
    /// (`erode::loose_grain_threshold`), and — where its pores hold liquid and
    /// air both — the menisci between its grains (`erode::capillary_cohesion`).
    /// Which of dry, damp and drowned it is, is read off its own state: the
    /// liquid it holds against the pore space its solid leaves at its packing.
    /// Less liquid than pore space is damp; more is water standing over it,
    /// and no meniscus survives that.
    pub fn loose_grain_strength(
        &self,
        idx: NodeIdx,
        material: &crate::material::Material,
        fluid_density: f64,
    ) -> Option<f64> {
        if !self.is_loose(idx) {
            return None;
        }
        let n = &self.tree.nodes[idx.get()];
        let mix = &n.matter.mixture;
        let (mut solid_mass, mut solid_volume) = (0.0, 0.0);
        let (mut liquid_mass, mut liquid_volume, mut tension) = (0.0, 0.0, 0.0);
        for p in mix.entries() {
            let Some(sub) = self.substances.get(p.substance) else { continue };
            let Some(c) = crate::eos::Condensed::of(&sub.props, p.phase) else { continue };
            let m = p.fraction * n.matter.mass;
            match p.phase {
                crate::chem::Phase::Solid => {
                    solid_mass += m;
                    solid_volume += m / c.rest_density;
                }
                crate::chem::Phase::Liquid => {
                    liquid_mass += m;
                    liquid_volume += m / c.rest_density;
                    tension += m * crate::erode::surface_tension(&sub.props);
                }
                _ => {}
            }
        }
        if !(solid_volume > 0.0) {
            return None;
        }
        let grain_density = solid_mass / solid_volume;
        let packing = (material.density / grain_density).clamp(1e-3, crate::erode::CLOSE_PACKING);
        let g = n.gravity.norm();
        let pores = solid_volume * (1.0 / packing - 1.0);
        let damp = liquid_volume > 0.0 && liquid_volume < pores;
        let buoyancy = if damp { 0.0 } else { fluid_density };
        let weight = crate::erode::loose_grain_threshold(grain_density, buoyancy, g, material.flaw_size, packing);
        let menisci = if damp && liquid_mass > 0.0 {
            crate::erode::capillary_cohesion(tension / liquid_mass, material.flaw_size, packing)
        } else {
            0.0
        };
        Some(weight + menisci)
    }

    /// Distance from the primary observer to a node — the light-travel distance
    /// an interaction has to cross before it takes effect.
    fn observer_distance(&self, target: NodeIdx) -> f64 {
        match self.observers.first() {
            Some(o) => self
                .tree
                .separation(o.anchor, o.offset, target, Vec3::ZERO)
                .value
                .norm(),
            None => 0.0,
        }
    }

    /// Perform a measurement, committing the outcome to the ledger.
    pub fn measure(
        &mut self,
        target: NodeIdx,
        instrument: Instrument,
        quantity: Quantity,
    ) -> Option<Reading> {
        if target.is_none() || !self.tree.nodes[target.get()].alive {
            return None;
        }
        let (key, matter, epoch) = {
            let n = &self.tree.nodes[target.get()];
            (n.key, n.matter, n.epoch)
        };
        let obs = *self.observers.first()?;
        let sep = self.tree.separation(obs.anchor, obs.offset, target, Vec3::ZERO);
        let d = sep.value.norm().max(1e-30);

        let view = match self.histories.get(&key) {
            Some(h) if !h.is_empty() => h.retarded(obs.offset, self.time),
            _ => crate::causal::RetardedView {
                snapshot: Moment {
                    t: self.time,
                    offset: sep.value + obs.offset,
                    velocity: self.tree.velocity_from(self.tree.root, target),
                    mass: matter.mass,
                    luminosity: matter.luminosity,
                    temperature: matter.temperature,
                },
                t_retarded: self.time - d / C,
                delay: d / C,
                distance: d,
                within_history: false,
            },
        };

        let seed = self.tree.world_seed;
        let time = self.time;
        // The reading itself is drawn from a stream addressed by the ledger
        // sequence, so repeating a measurement is a *new* measurement (with a
        // new disturbance), while re-querying a committed fact is not.
        let fact = self.ledger.get_or_sample(
            key,
            quantity,
            time,
            seed,
            epoch,
            |s| match quantity {
                Quantity::DecayTime => {
                    solvers::nuclear::Isotope::Neutron.sample_lifetime(s)
                }
                Quantity::Temperature => matter.temperature * (1.0 + 1e-3 * s.normal()),
                _ => s.uniform(),
            },
        );
        let _ = fact;

        let mut stream = Stream::at(seed, key.0, epoch, Purpose::PhotonEmission);
        let reading = read(instrument, &obs, &view, &matter, &mut stream);

        // Measurement disturbs. An interferometer deposits real energy, and the
        // engine applies it rather than reporting a free lunch.
        if let Reading::Position { disturbance, .. } = reading {
            if disturbance > 0.0 {
                let n = &mut self.tree.nodes[target.get()];
                n.matter.add_heat(disturbance);
                self.tree.pin(target);
                self.disturb(target);
            }
        }
        Some(reading)
    }

    /// Move a node under a different parent, carrying everything keyed by its
    /// path across with it.
    ///
    /// [`Tree::reparent`] does the frame change, the slot handling and the
    /// rekeying, and migrates the pinned detail it owns. Everything else
    /// addressed by `PathKey` lives on `World`, and this is where it follows.
    ///
    /// **The list below is the enumeration point.** A node's identity is spread
    /// across several side tables — each justified individually, by memory and
    /// by surviving a coarsen — and nothing else in the engine enumerates them
    /// together. A sixth table added without a line here would silently be left
    /// behind by every move, and the symptom would be a thrown object arriving
    /// without its chemistry. Add the line with the table.
    pub fn reparent(&mut self, node: NodeIdx, new_parent: NodeIdx) -> bool {
        let Some(moved) = self.tree.reparent(node, new_parent) else {
            return false;
        };
        // What a node *is* — its chemistry, its environment — is keyed by
        // `EntityId`, which a move does not change, so those tables are not
        // touched here. That is the point of issuing identity rather than
        // deriving it: a new side table of that kind is safe without anyone
        // remembering to add a line to this function.
        //
        // Three things are keyed by address and do move. The identity index,
        // because a node discarded and rebuilt recovers its name from its
        // address and nothing else. And the clock and the history, which are
        // keyed by address on purpose — they are bookkeeping the frame budget
        // created, and naming a node for one would make identity depend on how
        // fast the machine is. Their contents still have to travel: a clock
        // that changed rooms did not un-tick.
        //
        // Collected then reinserted, in two passes: within one move an old
        // address and a new one can name different nodes, so mutating in place
        // could overwrite an entry that had not been read yet.
        macro_rules! migrate {
            ($table:expr) => {{
                let mut taken = Vec::new();
                for (old, new) in &moved.keys {
                    if let Some(v) = $table.remove(old) {
                        taken.push((*new, v));
                    }
                }
                for (new, v) in taken {
                    $table.insert(new, v);
                }
            }};
        }
        migrate!(self.identities);
        migrate!(self.clocks);
        migrate!(self.histories);

        // Both ends changed by hand, so both need their detail kept rather than
        // regenerated, and their neighbours told.
        self.disturb(moved.from);
        self.disturb(moved.to);
        self.disturb(moved.moved);
        true
    }

    /// Reconcile every node the frame advanced with what it is actually
    /// holding — `docs/BACKLOG.md`'s split, and its inverse.
    ///
    /// **On the node's own cadence, which is what the backlog asks for.** Only
    /// a node whose contents have *moved* can have stopped being one
    /// neighbourhood, and the nodes that moved are exactly the ones the plan
    /// accepted a step for. A node nobody is advancing is not spreading either,
    /// so paying a grid build for it every frame would buy nothing.
    ///
    /// The merge is checked against the split's own siblings rather than over
    /// every pair in the parent: two nodes that are one neighbourhood are
    /// overlapping, and the cheap overlap test in [`Tree::would_merge`] rejects
    /// the rest before any grid is built.
    fn resolve_extents(&mut self, plan: &Plan) -> usize {
        let advanced: Vec<NodeIdx> = plan
            .accepted
            .iter()
            .filter(|t| t.kind == TaskKind::Step)
            .map(|t| t.node)
            .collect();
        let mut changed = 0;
        for node in advanced {
            if node.is_none() || !self.tree.nodes[node.get()].alive {
                continue;
            }
            // A structure's contents are members of a recipe, and whether they
            // are one thing is a question about its joins rather than about
            // their spacing. D15's break is what separates those, and it makes
            // nodes of the pieces itself.
            if self.tree.nodes[node.get()].morphology.is_some() {
                continue;
            }
            if let Some(new) = self.tree.resolve_extent(node) {
                self.stats.splits += 1;
                self.disturb(node);
                self.disturb(new);
                changed += 1;
                continue;
            }
            // Nothing split, so ask the other question: has this node come back
            // together with one of its siblings?
            let parent = self.tree.nodes[node.get()].parent;
            if parent.is_none() {
                continue;
            }
            let siblings: Vec<NodeIdx> = self.tree.nodes[parent.get()]
                .children
                .iter()
                .copied()
                .filter(|c| !c.is_none() && *c != node)
                .collect();
            for other in siblings {
                if self.tree.merge(node, other) {
                    self.stats.merges += 1;
                    self.disturb(node);
                    changed += 1;
                    break;
                }
            }
        }
        changed
    }

    /// Re-home everything that has left the region its parent owns.
    ///
    /// `docs/PLAY.md` D16's trigger, and the thing that was missing: `reparent`
    /// did the work correctly from the day it was written and was only ever
    /// called by hand, through `Interaction::Rehome`. Measured before this
    /// existed, a rocket at escape velocity left a 1 km forest node in 0.100 s
    /// and was still its child a hundred seconds and 1.12x10^6 m later.
    ///
    /// One pass, once a frame, after everything has moved. A node crosses at
    /// most one boundary per frame — the rocket reaches the planet on one frame
    /// and the star on a later one — which is both cheaper and more honest than
    /// iterating to a fixed point: a crossing is an event, and two of them are
    /// two events.
    ///
    /// The cost is a subtraction and a cube root per live node per frame, on a
    /// list `coast_to` has just walked anyway. Nodes appended during the pass —
    /// a body promoted to be arrived at — are deliberately not examined until
    /// the next frame; they have not moved yet.
    pub fn cross_boundaries(&mut self) -> usize {
        let mut crossed = 0;
        for i in 0..self.tree.nodes.len() {
            if self.cross(NodeIdx(i as u32)).is_some() {
                crossed += 1;
            }
        }
        crossed
    }

    /// Re-home one node, if it has left.
    ///
    /// **The parent arbitrates but is not necessarily the destination**, and
    /// that distinction is what keeps this one mechanism rather than two. A
    /// creature walking from forest to desert lands inside a sibling and should
    /// never become a direct child of the planet — it would be rekeyed twice,
    /// and for one frame its neighbours would be *other regions* rather than
    /// the ground under its feet. A rocket climbing out of the forest lands
    /// inside no sibling and genuinely belongs to the planet.
    ///
    /// So the question the arbiter answers is "which of the things I hold
    /// contains this now", over its bodies and its promoted children together —
    /// and where the answer is a body, that body is **promoted to meet it**,
    /// which is the inward row of D16's table. The detail about to be
    /// interacted with is generated by the crossing rather than by somebody
    /// having visited the place first.
    ///
    /// Where several occupants contain it, the *smallest* wins: "whatever now
    /// contains you" means the finest thing that does, or a creature would land
    /// in a continent when it is standing in a field.
    pub fn cross(&mut self, node: NodeIdx) -> Option<crate::crossing::Crossed> {
        use crate::crossing::{standing, Direction, Standing};
        if node.is_none() || !self.tree.nodes[node.get()].alive {
            return None;
        }
        let (parent, offset, radius) = {
            let n = &self.tree.nodes[node.get()];
            (n.parent, n.motion.offset, n.matter.radius)
        };
        if parent.is_none() {
            return None;
        }
        // Two bounds, cheap one first. `Tree::domain` is a subtraction and a
        // cube root; `Tree::contents_reach` is a pass over everything the
        // parent holds, and in a healthy world it is paid for a handful of
        // nodes rather than for all of them.
        if standing(offset, radius, self.tree.domain(parent)) != Standing::Outside {
            return None;
        }
        let arbiter = self.tree.nodes[parent.get()].parent;
        if arbiter.is_none() {
            // Checked *before* the reach, which costs a sweep of the parent's
            // contents: a child of the root has nowhere to go whatever the
            // measurement says, and every scenario on the shelf has some.
            // The universe has no outside. Counted rather than ignored: a
            // child of the root measuring as escaped every frame is a real
            // signal, and it is the one `docs/BACKLOG.md`'s sampler faults
            // produce — a granite block's promoted contents sit 2.7x10^3 of
            // their parent's radius out, a nucleus's 8x10^22.
            self.stats.crossings_refused += 1;
            return None;
        }
        // Beyond the sphere the parent claims is not the same as beyond what
        // the parent holds. See `Tree::contents_reach` for the ladder this was
        // measured on, where seven parcels in the tail of their parent's own
        // draw were re-homed and the ladder came apart.
        let reach = self.tree.contents_reach(parent, node);
        if standing(offset, radius, reach) != Standing::Outside {
            return None;
        }
        // How far over it went, recorded before anything is moved. See
        // `EngineStats::worst_crossing`.
        let over = (offset.norm() - radius) / reach.max(1e-300);
        if over > self.stats.worst_crossing {
            self.stats.worst_crossing = over;
            self.stats.worst_crossing_at = Some(self.tree.nodes[node.get()].key);
        }

        // Where it is, in the frame of the node that is about to decide.
        let here = self.tree.offset_from(arbiter, node, Vec3::ZERO).value;
        let slots = self.tree.nodes[arbiter.get()].bodies.len();
        let mut best: Option<(usize, f64)> = None;
        for slot in 0..slots {
            let Some((pos, claim)) = self.tree.occupant_claim(arbiter, slot) else {
                continue;
            };
            // The place it has just left cannot be the place it arrives at, and
            // it does not need excluding by name: it was measured as wholly
            // outside that volume a few lines ago, so it fails this test too.
            if standing(here - pos, radius, claim) != Standing::Inside {
                continue;
            }
            if best.is_none_or(|(_, c)| claim < c) {
                best = Some((slot, claim));
            }
        }

        let (destination, direction, generated) = match best {
            Some((slot, _)) => {
                let existing = self.tree.nodes[arbiter.get()].child_of(slot);
                if existing.is_none() {
                    let spec = self.tree.nodes[arbiter.get()].spec;
                    let made = self.tree.promote(arbiter, slot, spec);
                    if made.is_none() {
                        (arbiter, Direction::Outward, false)
                    } else {
                        // Deliberately **not** named. Only an event may name a
                        // node, and a crossing is not one however much it looks
                        // like one: which nodes a frame advances depends on a
                        // wall-clock allowance, so where everything is when the
                        // pass runs depends on how fast the machine is. Issuing
                        // an `EntityId` here made `next_entity` — which is
                        // persisted — a function of the budget, and
                        // `identity_does_not_depend_on_how_fast_the_machine_is`
                        // caught it within a frame of the pass existing. The
                        // place gets a name when something names it.
                        (made, Direction::Sideways, true)
                    }
                } else {
                    (existing, Direction::Sideways, false)
                }
            }
            None => (arbiter, Direction::Outward, false),
        };
        if destination == node || self.tree.is_ancestor(node, destination) {
            return None;
        }
        // The detail about to be met, generated before the arrival rather than
        // after it. `reparent` needs the destination materialised anyway — it
        // has to take a slot in a body list — so this is where the inward case
        // is paid for rather than an extra pass.
        self.tree.refine(destination);
        if !self.reparent(node, destination) {
            return None;
        }
        self.stats.crossings += 1;
        if direction == Direction::Sideways {
            self.stats.crossings_sideways += 1;
        }
        if generated {
            self.stats.crossings_generated += 1;
        }
        Some(crate::crossing::Crossed {
            node,
            from: parent,
            to: destination,
            direction,
            generated,
        })
    }

    /// Direct authoring. The one path that can violate conservation — so it
    /// records exactly how much it violated it by.
    fn author(&mut self, target: NodeIdx, property: Property, value: f64) {
        if target.is_none() || !self.tree.nodes[target.get()].alive {
            return;
        }
        let before = self.tree.nodes[target.get()].matter.total_energy();
        let key = {
            let n = &mut self.tree.nodes[target.get()];
            match property {
                Property::Mass => n.matter.mass = value.max(0.0),
                Property::Temperature => n.matter.set_temperature(value),
                Property::Radius => n.matter.radius = value.max(1e-30),
                Property::Charge => n.matter.charge = value,
                Property::Luminosity => n.matter.luminosity = value.max(0.0),
                // Reachable only if someone calls `author` directly; `dilate`
                // is the way in. It used to clamp, which meant this path and
                // that one disagreed: `dilate(-5.0)` was refused while
                // `Author { TimeRate, -5.0 }` silently became a millionfold
                // *slowdown*. Two ways to write one field must not have two
                // opinions about what is legal.
                Property::TimeRate => {
                    if let Some(r) = crate::dilation::accept_bubble(value) {
                        n.bubble = r;
                    }
                }
            }
            n.key
        };
        let after = self.tree.nodes[target.get()].matter.total_energy();
        self.audit.push(AuthorEvent {
            key,
            property,
            delta_energy: after - before,
            time: self.time,
        });
        if property == Property::Radius {
            // An authored size is still a size. `Interaction::Author` is the one
            // path that sets a value by hand and it is audited for exactly that
            // reason, but the tier is derived from the radius either way.
            self.tree.retier(target);
        }
        self.tree.pin(target);
        self.tree.bump_epoch(target);
        self.disturb(target);
    }

    /// Put a node and its subtree in a time bubble.
    ///
    /// Returns the rate actually applied, which differs from the one asked for
    /// when it had to be clamped. Returns `None` for a target that does not
    /// exist or a rate that is not a positive finite number — a caller bug, not
    /// an over-ambitious administrator, and worth telling apart.
    ///
    /// # Why this is not `author`
    ///
    /// It goes in the same audit trail, for the same reason: somebody reached
    /// in and changed something the physics did not. But it does three things
    /// `author` must not.
    ///
    /// It does not **pin**. Pinning says "this detail was altered and can no
    /// longer be regenerated", and a bubble alters nothing about the detail —
    /// the node's contents are exactly what the sampler would draw. Pinning
    /// every bubbled node would make a balancing pass permanently expensive.
    ///
    /// It does not **bump the epoch**. An epoch bump means the old detail is
    /// gone for good; a bubble does not invalidate anything, it only changes
    /// how fast what is there proceeds.
    ///
    /// It does not **disturb**. Disturbance resets the mixing-time clock that
    /// decides when detail may be released, and a bubble is not an event in the
    /// node's history — it is a statement about the observer's patience.
    pub fn dilate(&mut self, target: NodeIdx, rate: f64) -> Option<f64> {
        if target.is_none()
            || target.get() >= self.tree.nodes.len()
            || !self.tree.nodes[target.get()].alive
        {
            return None;
        }
        let accepted = crate::dilation::accept_bubble(rate)?;
        let key = {
            let n = &mut self.tree.nodes[target.get()];
            n.bubble = accepted;
            n.key
        };
        // Setting a rate injects no energy — the divergence a bubble causes is
        // a flow, not a step, and `EngineStats::bubble_seconds` is where it
        // accumulates. Recording zero here is the truthful entry, not a
        // placeholder.
        self.audit.push(AuthorEvent {
            key,
            property: Property::TimeRate,
            delta_energy: 0.0,
            time: self.time,
        });
        Some(accepted)
    }

    /// What a node is made of, by substance. Empty for anything nobody has
    /// given a composition to and nothing has handed one down.
    ///
    /// A read of the node's own matter since `docs/PLAY.md` D17. It used to
    /// consult a side table keyed by `EntityId`, which is why it needed an
    /// identity to answer at all — and why almost every node in the world
    /// answered "nothing".
    pub fn mixture_of(&self, idx: NodeIdx) -> crate::chem::Mixture {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return crate::chem::Mixture::new();
        }
        self.tree.nodes[idx.get()].matter.mixture
    }

    /// Say what a node is made of.
    ///
    /// The mixture is *speciation*: it says which substances account for the
    /// node's mass, and the elemental account it implies has to agree with the
    /// matter's own composition. This does not check that — nothing can,
    /// cheaply, for a partially speciated node — but `Mixture::composition` is
    /// how a caller finds out, and `tests/chem.rs` asserts it for the cases
    /// the engine builds itself.
    pub fn set_mixture(&mut self, idx: NodeIdx, mix: crate::chem::Mixture) {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return;
        }
        self.tree.nodes[idx.get()].matter.mixture = mix;
        self.refresh_rest_density(idx);
    }

    /// Derive and store the density a node's condensed matter rests at: its
    /// liquid's if it holds any, since a liquid is what lies and flows, and its
    /// solid's otherwise. Zero for matter the gas law describes. See
    /// `Node::rest_density`, and `eos.rs` for the law.
    ///
    /// Stored because a regeneration reads it and a client has no registry;
    /// refreshed wherever a node's mixture changes.
    pub fn refresh_rest_density(&mut self, idx: NodeIdx) {
        if idx.is_none() || idx.get() >= self.tree.nodes.len() {
            return;
        }
        let density = {
            let m = &self.tree.nodes[idx.get()].matter;
            if m.gas_law_applies() {
                0.0
            } else {
                let of = |want: crate::chem::Phase| {
                    let (mut mass, mut volume) = (0.0, 0.0);
                    for p in m.mixture.entries().iter().filter(|p| p.phase == want) {
                        let Some(s) = self.substances.get(p.substance) else { continue };
                        let Some(c) = crate::eos::Condensed::of(&s.props, want) else { continue };
                        mass += p.fraction;
                        volume += p.fraction / c.rest_density;
                    }
                    if volume > 0.0 { mass / volume } else { 0.0 }
                };
                let liquid = of(crate::chem::Phase::Liquid);
                if liquid > 0.0 { liquid } else { of(crate::chem::Phase::Solid) }
            }
        };
        self.tree.nodes[idx.get()].rest_density = density;
    }

    /// Run one pass of chemistry over every node that is made of something.
    ///
    /// O(substances) per node, which is at most eight, so this costs about what
    /// a single conditional per node costs and is affordable every frame for
    /// every node. Nodes with no mixture — every galaxy, star and planet the
    /// engine builds — are not visited at all.
    ///
    /// Runs on the node's *local* clock, so a region in a time bubble reacts
    /// faster along with everything else about it, and on the node's mixing
    /// time, because how fast a node stirs itself is exactly what sets how fast
    /// a solute reaches the far side of it.
    /// Advance every node's matter by the laws that need no detail.
    ///
    /// The counterpart to `advance_node`: that one runs a solver over bodies,
    /// this one runs physics over [`Matter`], and like growth it does so
    /// whether or not anything is materialised — in fact especially when
    /// nothing is. Without it a node that nobody is looking at is frozen but
    /// moving: `coast_to` advances its position and nothing advances its state,
    /// so a planet left alone for a century comes back at the same temperature
    /// it was, and every star in the world radiates into every scene while
    /// spending nothing.
    ///
    /// **One law here, and it is deliberately one.** A node radiates by
    /// Stefan-Boltzmann and absorbs what falls on it, and the difference goes
    /// into its thermal account. That is the whole of it, and it is enough to
    /// produce a great deal that is not written anywhere:
    ///
    /// * A hot thing cools, and a lit thing warms. Both directions, one
    ///   subtraction, no branch deciding which.
    /// * A node in shadow cools below freezing, `react_all` moves its liquid
    ///   into the solid phase with the latent heat booked, and
    ///   `environment_at` then measures no mobile solvent — so growth stops.
    ///   Snow accumulating and a winter are the same event seen twice, and
    ///   neither is written down.
    /// * Light returns, the solid melts, the liquid is measurable again and
    ///   growth resumes. Snowmelt, likewise unwritten.
    ///
    /// What it deliberately does *not* do is move anything between nodes. A
    /// node warms its own contents and radiates into the void, not onto its
    /// neighbour, because the engine has no notion of which nodes are adjacent
    /// — see the audit entry in the backlog. Until it does, meltwater has
    /// nowhere to run and a fire cannot spread to the next tree.
    fn evolve_matter(&mut self, dt: f64) -> MatterReport {
        let mut report = MatterReport::default();
        if !(dt > 0.0) {
            return report;
        }
        for i in 0..self.tree.nodes.len() {
            let idx = NodeIdx(i as u32);
            if !self.tree.nodes[i].alive {
                continue;
            }
            let local = dt * self.local_rate(idx);
            if !(local > 0.0) {
                continue;
            }
            // Only where a blackbody is what the node actually is.
            //
            // A node's `temperature` means two different things depending on
            // what it holds. For a planet or a rock it is a thermodynamic
            // temperature and Stefan-Boltzmann applies. For a galaxy or a star
            // cluster it is a *velocity dispersion* wearing the same field — a
            // way of saying how fast the members move relative to each other —
            // and a galaxy does not radiate as a 10^6 K blackbody the size of a
            // galaxy. Applied to those, this law cools the whole world to
            // nothing in a few frames, which is how the mistake was found.
            //
            // The tier is the honest place to draw it, because the tier is
            // exactly the statement of which physics describes the node.
            if self.tree.nodes[i].tier < crate::units::Tier::Planetary {
                continue;
            }
            let env = self.environment_at(idx);
            let n = &mut self.tree.nodes[i];
            let r = n.matter.radius;
            if !(r > 0.0) || !n.matter.is_finite() {
                continue;
            }
            // Radiated from the whole surface; absorbed over the cross-section
            // the node presents to what is illuminating it. The factor of four
            // between the two is why a body's equilibrium temperature is what
            // it is, and it falls out rather than being put in.
            let area = std::f64::consts::PI * r * r;
            let radiated = crate::state::stefan_boltzmann(r, n.matter.temperature) * local;
            let absorbed = env.light_flux * area * local;
            let net = absorbed - radiated;
            if net == 0.0 {
                continue;
            }
            // A node cannot radiate away more than it has. Floored rather than
            // clamped silently: the shortfall is reported, because a world
            // where this is large is one whose nodes are being asked to shine
            // brighter than their own reserves.
            let available = n.matter.internal_energy;
            let applied = if net < 0.0 && -net > available {
                report.radiation_deficit += -net - available;
                -available
            } else {
                net
            };
            n.matter.add_heat(applied);
            n.matter.luminosity = crate::state::stefan_boltzmann(r, n.matter.temperature);
            report.radiated += radiated;
            report.absorbed += absorbed;
            report.nodes += 1;
        }
        report
    }

    fn react_all(&mut self, dt: f64) -> crate::chem::ReactionReport {
        let mut total = crate::chem::ReactionReport::default();
        if !(dt > 0.0) {
            return total;
        }
        // Gated on a mixture that is actually there, which is what `PLAY.md`
        // D17 asks for: `react` costs 0.069 µs per node for one substance and
        // 0.934 for eight, and running it over `UNSPECIATED` would spend
        // between 1.4% and 18.7% of a frame at 10^4 nodes saying nothing. The
        // gate used to be "does a side table hold an entry for this node's
        // identity"; it is now a field on the node's own matter.
        let live: Vec<NodeIdx> = (0..self.tree.nodes.len())
            .map(|i| NodeIdx(i as u32))
            .filter(|i| {
                let n = &self.tree.nodes[i.get()];
                n.alive && !n.matter.mixture.is_empty()
            })
            .collect();
        if live.is_empty() {
            return total;
        }
        for idx in live {
            let (temperature, mass, mut mix) = {
                let n = &self.tree.nodes[idx.get()];
                (n.matter.temperature, n.matter.mass, n.matter.mixture)
            };
            let local = dt * self.local_rate(idx);
            let tau = self.mixing_time(idx);
            let r = crate::chem::react(&mut mix, &self.substances, temperature, local, tau);
            self.tree.nodes[idx.get()].matter.mixture = mix;
            if r.quiet() && r.heat == 0.0 {
                continue;
            }
            // A phase moved, so what the condensed matter rests at may have.
            self.refresh_rest_density(idx);
            // Latent heat is real energy and comes out of the node's own
            // internal account. Positive `heat` was absorbed by the matter, so
            // it leaves the thermal store.
            let n = &mut self.tree.nodes[idx.get()];
            n.matter.internal_energy = (n.matter.internal_energy - r.heat * mass).max(0.0);
            // **A melt that froze lays down grains.** The one moment at which
            // "this froze" is a fact rather than a guess is the moment the mass
            // moves, so the layout is derived here rather than in the frame's
            // sweep — see `World::derive_frozen_layout`.
            if r.frozen > 0.0 && self.tree.nodes[idx.get()].morphology.is_none() {
                self.derive_frozen_layout(idx);
            }
            total.melted += r.melted;
            total.frozen += r.frozen;
            total.boiled += r.boiled;
            total.condensed += r.condensed;
            total.dissolved += r.dissolved;
            total.precipitated += r.precipitated;
            total.heat += r.heat;
            total.unresolved += r.unresolved;
        }
        total
    }

    /// Nodes currently carrying a bubble, with the factor each was given.
    pub fn bubbles(&self) -> Vec<(NodeIdx, f64)> {
        (0..self.tree.nodes.len())
            .map(|i| NodeIdx(i as u32))
            .filter(|i| {
                let n = &self.tree.nodes[i.get()];
                n.alive && n.bubble != 1.0
            })
            .map(|i| (i, self.tree.nodes[i.get()].bubble))
            .collect()
    }

    /// Everything the given observer can currently see, nearest first.
    pub fn look(&mut self, observer: usize, instrument: Instrument) -> Vec<Sighting> {
        let obs = match self.observers.get(observer) {
            Some(o) => *o,
            None => return Vec::new(),
        };
        let mut out = Vec::new();
        let live: Vec<NodeIdx> = (0..self.tree.nodes.len())
            .map(|i| NodeIdx(i as u32))
            .filter(|i| self.tree.nodes[i.get()].alive && *i != obs.anchor)
            .collect();

        for idx in live {
            let (key, matter) = {
                let n = &self.tree.nodes[idx.get()];
                (n.key, n.matter)
            };
            let sep = self.tree.separation(obs.anchor, obs.offset, idx, Vec3::ZERO);
            let d = sep.value.norm();
            if d <= 0.0 || !obs.sees(sep.value) {
                continue;
            }
            let view = match self.histories.get(&key) {
                Some(h) if !h.is_empty() => h.retarded(obs.offset, self.time),
                _ => crate::causal::RetardedView {
                    snapshot: Moment {
                        t: self.time,
                        offset: sep.value + obs.offset,
                        velocity: self.tree.velocity_from(self.tree.root, idx),
                        mass: matter.mass,
                        luminosity: matter.luminosity,
                        temperature: matter.temperature,
                    },
                    t_retarded: self.time - d / C,
                    delay: d / C,
                    distance: d,
                    within_history: false,
                },
            };
            let theta = crate::coords::angular_size(matter.radius, d);
            let dir = (view.snapshot.offset - obs.offset).unit();
            let dop = crate::coords::doppler(view.snapshot.velocity - obs.velocity, dir);
            let mut stream = Stream::at(self.tree.world_seed, key.0, 0, Purpose::PhotonEmission);
            let reading = read(instrument, &obs, &view, &matter, &mut stream);
            out.push(Sighting {
                node: idx,
                key,
                view,
                angular_size: theta,
                resolved: theta >= obs.angular_resolution,
                required_tier: obs.required_tier(d),
                flux: flux(matter.luminosity, d, dop),
                doppler: dop,
                reading,
            });
        }
        out.sort_by(|a, b| {
            a.view
                .distance
                .partial_cmp(&b.view.distance)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        out
    }

    /// Drill from a node down to a target *size* along the most massive child
    /// at each step, materialising and promoting as it goes.
    ///
    /// This is the "zoom in" primitive: the path from a galaxy to a nucleus is
    /// one call, and the engine materialises exactly the chain of nodes along
    /// the way and nothing else. That chain is a few thousand bodies, not 10^66.
    ///
    /// `target` is a physical scale — the smallest thing that has to exist.
    /// From a viewpoint it is [`crate::observe::Observer::linear_resolution`]
    /// at the distance to the node, the smallest separation that observer
    /// could tell apart. From a physical interaction it is whatever scale the
    /// interaction acts at: a contact patch, a wavelength, a blast radius.
    ///
    /// **The target is a length rather than a tier deliberately.** A tier is a
    /// label `promote` derives from the promoted *body's* radius, and a label
    /// can disagree with what the node it names actually holds — a node full
    /// of nucleons can be marked `Atomic` and then never satisfy
    /// `tier >= Tier::Nuclear` however far the descent goes. A radius cannot
    /// disagree with itself, and it shrinks at every step, so the loop ends
    /// because it arrived rather than because it ran out of patience.
    pub fn drill_to(
        &mut self,
        from: NodeIdx,
        target: f64,
        specs: &dyn Fn(Tier) -> SampleSpec,
    ) -> Vec<NodeIdx> {
        let mut path = vec![from];
        let mut cur = from;
        // Still capped. The cap is now a guard against a tree that misbehaves
        // — a promoted body no larger than its parent, say — rather than the
        // thing that ends an ordinary descent.
        for _ in 0..64 {
            let (tier, radius) = {
                let n = &self.tree.nodes[cur.get()];
                (n.tier, n.matter.radius)
            };
            // Written so that a non-finite radius or target stops rather than
            // spinning: `!(a >= b)` is false only when the comparison holds.
            if !(radius >= target) {
                break;
            }
            self.tree.refine(cur);
            let best = {
                let n = &self.tree.nodes[cur.get()];
                let mut bi = 0usize;
                let mut bm = -1.0;
                for (i, b) in n.bodies.iter().enumerate() {
                    if b.mass > bm {
                        bm = b.mass;
                        bi = i;
                    }
                }
                bi
            };
            if self.tree.nodes[cur.get()].bodies.is_empty() {
                break;
            }
            let spec = specs(tier.finer());
            let child = self.tree.promote(cur, best, spec);
            if child.is_none() {
                break;
            }
            self.tree.nodes[child.get()].residency = Residency::Observed;
            path.push(child);
            cur = child;
        }
        path
    }

    /// Total conserved quantities over the whole live world.
    pub fn conserved(&self) -> crate::state::Conserved {
        self.tree.total_conserved()
    }

    /// Largest causality violation between any two materialised nodes: the
    /// invariant the scheduler exists to protect.
    ///
    /// **Measured on what the nodes hold** — `Node::time`, not
    /// `Node::carried` — which is the owner's decision for Phase 5. Where a
    /// node is, is carried to the world instant every frame and cannot be
    /// skewed; what it holds is behind for as long as it waits to be solved,
    /// and catches up by its ensemble wherever it may be redrawn (`coast_to`).
    /// A node that may not — pinned, edited, bubbled, holding children — and
    /// that its own physics cannot carry at the world's pace stays behind,
    /// and this is where it shows. Measured in `phys-demo`: a pinned Atomic
    /// node of 56 bodies at 2.4e-17 s in a world paced at 8.9e7 s a frame,
    /// 5.5e12 s behind a neighbour, where the coasting this replaced reported
    /// zero by relabelling frozen contents as current.
    pub fn check_causality(&self) -> f64 {
        let live: Vec<NodeIdx> = (0..self.tree.nodes.len())
            .map(|i| NodeIdx(i as u32))
            .filter(|i| self.tree.nodes[i.get()].alive)
            .collect();
        let mut worst: f64 = 0.0;
        for (i, a) in live.iter().enumerate() {
            for b in live.iter().skip(i + 1) {
                // Only disjoint regions can violate causality with respect to
                // each other; a node and its own ancestor are the same matter
                // described at two resolutions. See `Tree::sibling_separations`.
                if self.tree.is_ancestor(*a, *b) || self.tree.is_ancestor(*b, *a) {
                    continue;
                }
                let na = &self.tree.nodes[a.get()];
                let nb = &self.tree.nodes[b.get()];
                let gap = (self
                    .tree
                    .separation(*a, Vec3::ZERO, *b, Vec3::ZERO)
                    .value
                    .norm()
                    - na.matter.radius
                    - nb.matter.radius)
                    .max(0.0);
                worst = worst.max(crate::causal::causality_violation(na.time, nb.time, gap));
            }
        }
        worst
    }

    /// A one-line summary for the debug overlay.
    pub fn summary(&self) -> String {
        format!(
            "t={} nodes={} bodies={} detail={:.1} MB ledger={} facts ({:.1} kB) frame={:.1} ms debt={:.2e} cons_err={:.2e}",
            fmt_time(self.time),
            self.stats.live_nodes,
            self.stats.materialised_bodies,
            self.tree.detail_bytes() as f64 / 1e6,
            self.ledger.len(),
            self.ledger.bytes() as f64 / 1e3,
            self.stats.last_frame_us / 1e3,
            self.stats.detail_debt,
            self.tree.stats.worst_conservation_error,
        )
    }
}

/// Build the standard scenario: a disc galaxy inside a dark halo.
///
/// The dark matter is modelled as a static NFW-like potential rather than as
/// mixed into the baryonic composition. That is both the standard approach in
/// galaxy simulation and the only one that keeps the books straight here: dark
/// matter carries no baryon number, and folding it into the composition would
/// have the engine report ~10^67 baryons for a galaxy that contains 10^66.
pub fn galaxy(world_seed: u64, stars: f64) -> Tree {
    let mass_stars = stars * 0.8 * M_SUN;
    let mass_gas = mass_stars * 0.15;
    let mass_dark = (mass_stars + mass_gas) * 8.0;
    let total = mass_stars + mass_gas;
    let radius = 15.0 * KPC;

    let mut matter = Matter::neutral(total, radius, 1e4, crate::state::Composition::primordial());
    // A galaxy's internal energy is dominated by orbital motion, not heat:
    // the virial theorem sets it, so derive it rather than inventing a number.
    // The halo dominates the potential, so the velocity dispersion is set by
    // the *total* enclosed mass even though only the baryons are represented.
    let enclosed = total + mass_dark;
    let sigma = (G * enclosed / (2.0 * radius)).sqrt();
    matter.internal_energy = 0.5 * total * sigma * sigma;
    // The baryons' own binding, which refinement can and must reproduce...
    matter.gravitational_binding = -0.6 * G * total * total / radius;
    // ...and the halo's grip on them, which it cannot, because the halo is not
    // made of the thing being refined.
    matter.external_potential = -G * total * mass_dark / radius;
    // Angular momentum of a rotationally supported disc.
    matter.spin = crate::math::v3(0.0, 0.0, 0.7 * total * sigma * radius);
    matter.luminosity = stars * 3.828e26 * 0.3;

    let spec = SampleSpec {
        count: 20_000,
        profile: crate::sampler::Profile::Disk {
            scale_height_ratio: 0.12,
        },
        spectrum: crate::sampler::MassSpectrum::Equal,
        kind: crate::state::BodyKind::Super,
        composition_scatter: 0.15,
        turbulent_fraction: 0.3,
    };
    Tree::new(world_seed, matter, Tier::Galactic, spec)
}

/// The default refinement policy, and the budgeted form of it.
///
/// Both live in [`crate::sampler`] because [`crate::tree::Tree::promote`] needs
/// them: a child's tier follows from its *size*, so a caller that guesses the
/// tier from the parent can be wrong, and the node then materialises under a
/// policy meant for a different scale. The tier has to be able to reach for its
/// own policy, and the tier is decided below this module.
pub use crate::sampler::{budgeted_spec, default_spec};

/// Give one side of a contact its impulse, its couple and its share of the heat.
///
/// A free function rather than a method because `contact_within` is holding a
/// borrow of the node's neighbourhood while it works, and the two halves of a
/// contact have to be applied to different places — a body in this node, or a
/// node of its own one level down.
fn apply_contact(
    w: &mut World,
    parent: NodeIdx,
    who: crate::neighbourhood::Occupant,
    impulse: Vec3,
    spin: Vec3,
    heat: f64,
    now: f64,
) {
    match who {
        crate::neighbourhood::Occupant::Body(k) => {
            if let Some(b) = w.tree.nodes[parent.get()].bodies.get_mut(k as usize) {
                if b.mass > 0.0 {
                    b.vel += impulse.scale(1.0 / b.mass);
                }
                b.spin += spin;
                if heat != 0.0 {
                    b.add_heat(heat);
                }
            }
        }
        crate::neighbourhood::Occupant::Child(c) => {
            if c.is_none() || !w.tree.nodes[c.get()].alive {
                return;
            }
            // The impulse lands at the parent's instant and the child's motion
            // is carried to a later one, so its position there moves by the
            // change for the difference. See `Node::carried`.
            let since = (w.tree.nodes[c.get()].carried - w.tree.nodes[parent.get()].time).max(0.0);
            let n = &mut w.tree.nodes[c.get()];
            if n.matter.mass > 0.0 {
                let dv = impulse.scale(1.0 / n.matter.mass);
                n.motion.velocity = n.motion.velocity + dv;
                n.motion.offset = n.motion.offset + dv.scale(since);
            }
            n.matter.spin += spin;
            // Angular momentum arrived, so the angular velocity it implies has
            // changed. Without this line the node banks the momentum and never
            // turns, which is what `PLAY.md` §2A measured.
            n.sync_spin_rate();
            if heat != 0.0 {
                // Through the mailbox, like every other joule crossing into
                // another node's books. See `World::exchange_within`.
                let d = n.motion.offset.norm();
                w.mailbox
                    .post(c, now, d, InfluenceKind::Exchange, heat, Vec3::ZERO);
            }
        }
    }
}
