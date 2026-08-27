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
use crate::causal::{CausalGate, Clock, History, Influence, InfluenceKind, Mailbox, Snapshot};
use crate::ids::{NodeIdx, PathKey};
use crate::math::Vec3;
use crate::observe::*;
use crate::prolong::ProlongSpec;
use crate::rng::{Purpose, Stream};
use crate::solvers::{self, SolverKind};
use crate::state::{Aggregate, Body};
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

/// Refinement error above which a node resolves itself, with or without an
/// audience.
///
/// `refinement_error` is about 0.05 for a node the aggregate describes
/// perfectly and climbs past one when the bulk state has started lying: an
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
        if topo.support[i] == crate::morph::NO_SUPPORT && topo.bonds[i].radius > 0.0 {
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
    /// Nodes carried forward in closed form without being re-solved. The
    /// overwhelming majority, every frame, and the reason one instant is
    /// affordable across thirty-eight orders of magnitude.
    pub coasted: usize,
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
    /// Set from what is being watched rather than from any solver's stability
    /// limit. One frame covers about one characteristic time of the node the
    /// view is paced to, so a galaxy advances millennia per frame and a carbon
    /// atom femtoseconds, and both are watchable. See [`World::pace_to`].
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
    pub time_throttle: f64,
    /// The node the pace is taken from, re-read at the start of every frame.
    ///
    /// Held as a node rather than a number because a node's cadence changes
    /// under it: materialising a galaxy into twenty thousand stars shortens its
    /// characteristic time by two orders of magnitude, and a pace fixed when it
    /// was still a single aggregate would then be asking for a span its own
    /// stars could not be integrated across.
    pub paced_to: NodeIdx,
    /// Retained history, only for nodes something might observe from a
    /// distance. Keeping it keyed by path rather than in the node itself means
    /// a node can be coarsened and rebuilt without losing its past.
    pub histories: HashMap<PathKey, History>,
    pub clocks: HashMap<PathKey, Clock>,
    pub stats: EngineStats,
    /// Actions that deliberately broke conservation, with what they cost.
    pub audit: Vec<AuthorEvent>,
    /// Whether cost estimates assume the GPU path.
    pub gpu: bool,
    /// Construction rate available to planned programs, fraction of a design
    /// per second. Zero means no crews are working.
    pub labour_rate: f64,
    /// Growth steps refused because their transaction did not balance.
    pub rejected_transactions: u64,
    /// Per-node environment overrides, keyed by path so they survive the node
    /// being coarsened and rebuilt.
    pub environments: HashMap<PathKey, crate::morph::Environment>,
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
            stats: EngineStats::default(),
            audit: Vec::new(),
            gpu: false,
            labour_rate: 0.0,
            rejected_transactions: 0,
            environments: HashMap::new(),
            shaking: Vec::new(),
            falling: Vec::new(),
            history_depth: 64,
        };
        // Start paced to the root, so a world is watchable the moment it is
        // built without anybody having to know what timescale it lives on.
        let root = w.tree.root;
        w.pace_to(root);
        w
    }

    /// Place a structure and give it conditions to grow in.
    pub fn plant(
        &mut self,
        idx: NodeIdx,
        program: crate::morph::Program,
        env: crate::morph::Environment,
    ) {
        let key = self.tree.nodes[idx.get()].key;
        self.tree.plant(idx, program);
        self.environments.insert(key, env);
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
        let tasks = self.survey(horizon);
        let bytes = self.tree.detail_bytes();
        let plan = self.budget.plan(tasks, bytes);
        let achieved = self.execute(&plan, horizon);
        let coasted = self.coast_to(horizon);

        self.deliver_influences(horizon);
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
        let mut dt = self.pace * self.time_rate * self.time_throttle;
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
    fn retime(&mut self, achieved: f64) {
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
        self.paced_to = idx;
        self.refresh_pace();
    }

    /// Re-read the pace from the node it is taken from.
    ///
    /// Called at the top of every frame, because the subject's cadence moves:
    /// the moment a galaxy is materialised into its stars, the span a frame may
    /// cover has to come down with it or the stars cannot be integrated across
    /// one. This is the feedback that makes "zooming in slows time" an
    /// arithmetic consequence rather than a policy.
    pub fn refresh_pace(&mut self) {
        let idx = self.paced_to;
        if idx.is_none() || !self.tree.nodes[idx.get()].alive {
            return;
        }
        let tau = self.node_cadence(idx);
        if tau.is_finite() && tau > 0.0 {
            self.pace = tau;
        }
    }

    /// The length scale a node is currently represented at, metres.
    ///
    /// Not the node's radius — the size of the smallest thing it is currently
    /// showing. A planet held as a single aggregate is represented at its own
    /// radius; the same planet split into four thousand parcels is represented
    /// at four hundred kilometres, and has to be re-solved sixteen times as
    /// often for the difference to mean anything.
    pub fn node_resolution(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        let parts = n.bodies.len();
        if parts > 1 {
            n.agg.radius / (parts as f64).cbrt()
        } else {
            n.agg.radius
        }
    }

    /// How long this node may be left alone before its state is visibly stale.
    ///
    /// For a node held as bulk state this is [`Aggregate::characteristic_time`]
    /// — how long before the thing moves, turns, or rearranges by its own size.
    /// For a materialised node it is the bodies that are represented, so it is
    /// the bodies that set the cadence: the time for the fastest of them to
    /// cross one resolution element. Thermal motion counts in that case and not
    /// in the first, and the difference is not a fudge — at bulk resolution
    /// thermal motion is a temperature, and at parcel resolution it is
    /// something you can watch happen.
    pub fn node_cadence(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        if n.is_materialised() {
            let h = self.node_resolution(idx);
            // The bodies' own speeds, and nothing else. Seeding this with the
            // aggregate's sound speed looked harmless and was not: a galaxy's
            // "sound speed" is a gas-pressure formula applied to a collisionless
            // stellar system, and it saturated at 0.577c, so a materialised
            // galaxy claimed to need re-solving four orders of magnitude more
            // often than its own stars could justify. Prolongation samples
            // thermal motion into the bodies already, so where a sound speed is
            // meaningful it is in here anyway.
            let v = n.bodies.iter().map(|b| b.vel.norm()).fold(0.0f64, f64::max);
            return if v > 0.0 { h / v } else { f64::INFINITY };
        }
        n.agg.characteristic_time(n.agg.radius)
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
        ((horizon - self.tree.nodes[idx.get()].last_solved) / cadence).clamp(0.0, 1e9)
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
    /// merely *a* state of the same bulk.
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
        let dispersion = n.agg.velocity_dispersion();
        let ordered = n.agg.angular_velocity().norm() * n.agg.radius;
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
                (n.tier, n.is_materialised(), n.bodies.len(), n.agg.radius)
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

            // Materialise when the bulk state cannot express what the node is
            // doing — and when the world's pace leaves room to actually
            // integrate it.
            //
            // Two independent reasons, and the second is the one that was
            // missing. `acuity` is somebody wanting to see it; `error` is the
            // node's own state saying the aggregate is no longer an adequate
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
            // step. Detail now goes when the node has had a mixing time to
            // forget what put it there — see [`World::mixing_time`] — and not
            // before, whoever is or is not watching.
            let forgotten =
                self.time - self.tree.nodes[idx.get()].last_disturbed > self.mixing_time(idx);
            if materialised && forgotten && acuity < 0.25 && !self.tree.nodes[idx.get()].pinned {
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
            // aggregate representation: a forest of 10^9 trees held as 10^4
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
    /// Two physical criteria, both of which are about structure the bulk state
    /// cannot represent: an unresolved Jeans length (the node is about to
    /// fragment) and a short dynamical time relative to the frame step (the
    /// node is evolving faster than we are looking at it).
    fn refinement_error(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        let a = &n.agg;
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
                        self.grow_node(task.node, dt);
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
        let h0 = self.node_dt(idx);
        if h0 > 0.0 && h0.is_finite() && span / h0 > MAX_SUBSTEPS as f64 && self.forgettable(idx) {
            self.thermalise(idx, horizon);
            return;
        }
        let mut steps = 0u32;
        while self.tree.nodes[idx.get()].time < horizon && steps < allowance {
            let remaining = horizon - self.tree.nodes[idx.get()].time;
            let h = self.node_dt(idx).min(remaining);
            if !(h > 0.0) {
                break;
            }
            self.advance_node(idx, h);
            steps += 1;
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
        (self.time - self.tree.nodes[idx.get()].last_solved).max(0.0)
    }

    /// May this node's detail be thrown away and drawn again?
    ///
    /// No, if somebody has touched it — a tree you broke is not a
    /// representative sample of anything. No, if something finer has been built
    /// on it, because releasing it would take the whole subtree with it.
    pub fn forgettable(&self, idx: NodeIdx) -> bool {
        let n = &self.tree.nodes[idx.get()];
        !n.pinned && !n.children.iter().any(|c| !c.is_none())
    }

    /// Cross a span too long to integrate, by ensemble instead of trajectory.
    ///
    /// This is the step that makes one shared instant affordable across
    /// thirty-eight orders of magnitude. A resolved nucleus asked to cross
    /// fifty milliseconds would need 10^21 steps; it also does not need them,
    /// because over that span it has sampled its accessible states 10^21 times
    /// and where it ends up is a draw from its equilibrium ensemble, not the
    /// endpoint of a trajectory. So the detail is restricted back to the bulk
    /// state, the bulk state is carried across in closed form, and the detail
    /// is drawn again at the far end.
    ///
    /// Both halves are things the engine already guarantees: restriction is
    /// conservative to within `IDEMPOTENT_TOLERANCE` and prolongation is a
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
        let n = &mut self.tree.nodes[idx.get()];
        let dt = horizon - n.time;
        if dt > 0.0 {
            n.frame.advance(dt);
        }
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
    fn coast_to(&mut self, horizon: f64) -> usize {
        let mut coasted = 0;
        for i in 0..self.tree.nodes.len() {
            let n = &mut self.tree.nodes[i];
            if !n.alive {
                continue;
            }
            let dt = horizon - n.time;
            if !(dt > 0.0) {
                continue;
            }
            n.frame.advance(dt);
            n.time = horizon;
            coasted += 1;
            let key = n.key;
            let velocity = n.frame.velocity;
            if let Some(c) = self.clocks.get_mut(&key) {
                c.time = horizon;
                c.proper_time += crate::coords::proper_time_step(dt, velocity);
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
        let mut natural = n
            .tier
            .dt()
            .min(n.agg.dynamical_time() / 50.0)
            .min(0.25 * n.agg.signal_crossing(parts, flow));
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
            (n.tier, n.key, n.epoch, n.agg.radius, n.bodies.len(), n.steps_taken)
        };
        if count == 0 || dt <= 0.0 {
            return solvers::SolveReport::default();
        }
        let seed = self.tree.world_seed;
        let bodies = &mut self.tree.nodes[idx.get()].bodies;

        let report = match solvers::for_tier(tier) {
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
                let params = solvers::hydro::HydroParams {
                    h: radius / (count as f64).cbrt() * 1.2,
                    ..Default::default()
                };
                solvers::hydro::step(bodies, dt, params)
            }
            SolverKind::MolecularDynamics => {
                let params = solvers::md::MdParams::default();
                // Substep to whatever the force field needs. The scheduler's
                // timestep answers to causality and to the tier; the Lennard-
                // Jones potential answers to neither, and a molecular system
                // handed a step longer than its own vibrational period does not
                // integrate inaccurately, it detonates.
                let stable = solvers::md::configuration_dt(bodies, params).max(1e-24);
                let substeps = ((dt / stable).ceil() as u32).clamp(1, 64);
                let h = dt / substeps as f64;
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
                total.dt_used = dt;
                total
            }
            SolverKind::Statistical => self.advance_statistical(idx, dt),
        };

        self.stats.bodies_stepped += count as u64;
        let n = &mut self.tree.nodes[idx.get()];
        n.time += dt;
        n.steps_taken += 1;
        n.frame.advance(dt);
        let clock = self
            .clocks
            .entry(key)
            .or_insert_with(|| Clock::new(n.time, tier.dt()));
        clock.time = n.time;
        clock.proper_time += crate::coords::proper_time_step(dt, n.frame.velocity);
        report
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
        }
    }

    /// Advance one structure's developmental state.
    ///
    /// The environment is read off the node's own aggregate, so a structure in
    /// a cold or crowded node grows slowly without anyone having to arrange it.
    /// The transaction is validated before it is applied: a growth program
    /// cannot mint free energy or order, it can only trade for them.
    pub fn grow_node(&mut self, idx: NodeIdx, dt: f64) -> Option<crate::morph::Transaction> {
        if idx.is_none() || !self.tree.nodes[idx.get()].alive || dt <= 0.0 {
            return None;
        }
        let env = self.environment_at(idx);
        let node = &mut self.tree.nodes[idx.get()];
        let morph = node.morphology.as_mut()?;
        let txn = morph.advance(dt, &env);
        if txn.validate().is_err() {
            // A program that cannot balance its books does not get to run. This
            // is a bug in the program, not a condition to be smoothed over.
            self.rejected_transactions += 1;
            return None;
        }
        let extent = morph.extent().max(1e-30);
        let stored = morph.stored_energy();

        // Apply the transaction to the aggregate. Mass moves *within* the node
        // — carbon from its air into its wood — so mass, composition and baryon
        // number are all unchanged, and only the energy and entropy accounts
        // move. What crosses the boundary is energy, and it is booked.
        node.agg.chemical_energy = stored;
        // Only the thermalised share stays. What was re-radiated has left the
        // node, and adding it here would cook a forest in a season.
        node.agg.internal_energy += txn.heat_released;
        node.agg.entropy += txn.entropy_local;
        node.agg.entropy_exported += txn.entropy_exported;
        node.agg.radius = extent;
        node.agg.luminosity = crate::state::stefan_boltzmann(extent, node.agg.temperature);

        // The structure it would generate has changed, so any materialised copy
        // is stale. Discarding it is correct and cheap — it is regenerable.
        node.bodies.clear();
        node.children.clear();

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

        let ambient = self.tree.nodes[idx.get()].agg.temperature;
        let bodies = self.tree.nodes[idx.get()].bodies.clone();
        let mut topo = match self.tree.nodes[idx.get()].topology.clone() {
            Some(t) => t,
            None => return out,
        };
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
        field.apply(&st::weather::gravity(), &bodies, &topo);

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

        let node = &mut self.tree.nodes[idx.get()];
        let temperature = node.agg.temperature;
        if let Some(m) = node.morphology.as_mut() {
            if !sites.is_empty() && structural_mass > 0.0 {
                let fraction = (failures.detached_mass / structural_mass).clamp(0.0, 1.0);
                let txn = m.sever_many(&sites, fraction);
                out.detached_mass = txn.mass_detached;
            }
            if insult_report.consumed_mass > 0.0 {
                let burn = m.consume(insult_report.consumed_mass, temperature);
                if burn.validate().is_ok() {
                    // Burning releases the free energy the wood was holding.
                    // The atoms stay in the node as combustion products, so mass
                    // and baryon number are untouched; only the energy moves.
                    node.agg.chemical_energy -= burn.energy_released;
                    node.agg.internal_energy += burn.heat_released;
                    node.agg.entropy += burn.entropy_local;
                    node.agg.entropy_exported += burn.entropy_exported;
                    out.energy_released = burn.energy_released;
                } else {
                    self.rejected_transactions += 1;
                }
            }
            node.agg.chemical_energy = m.stored_energy();
            node.agg.radius = m.extent().max(1e-30);
        }
        // The structure it would generate has changed.
        node.bodies.clear();
        node.topology = None;
        node.children.clear();
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
        let ambient = self.tree.nodes[idx.get()].agg.temperature;
        let bodies = self.tree.nodes[idx.get()].bodies.clone();
        let topo = match self.tree.nodes[idx.get()].topology.clone() {
            Some(t) => t,
            None => return out,
        };

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
        field.apply(&st::weather::gravity(), &bodies, &topo);

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
        if let Some(m) = node.morphology.as_mut() {
            if !sites.is_empty() && structural_mass > 0.0 {
                let fraction = (detached / structural_mass).clamp(0.0, 1.0);
                let txn = m.sever_many(&sites, fraction);
                out.detached_mass = txn.mass_detached;
            }
            node.agg.chemical_energy = m.stored_energy();
            node.agg.radius = m.extent().max(1e-30);
        }
        node.bodies.clear();
        node.topology = None;
        node.children.clear();
        self.shaking.retain(|(n, _)| *n != idx);
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
        for node in nodes {
            let bodies = self.tree.refine(node).to_vec();
            if let Some(topo) = self.tree.nodes[node.get()].topology.clone() {
                struck.insert(node.get(), (bodies, topo));
            }
        }

        let mut strikes: Vec<(NodeIdx, u32, crate::math::Vec3)> = Vec::new();
        for (node, frag) in self.falling.iter_mut() {
            frag.age += dt;
            let n = frag.dynamics.dynamics.frame.nodes.len();
            let mut load = vec![Dof::default(); n];
            for i in 0..n {
                let m = frag.dynamics.dynamics.frame.lumped[i].t.z;
                load[i].t = st::G_EARTH.scale(m);
            }
            let rep = frag.dynamics.dynamics.step(&load, dt);
            report.broken_while_falling += rep.broken.len();

            let contacts = match struck.get(&node.get()) {
                Some((bodies, topo)) => frag.contacts(bodies, topo, ground_of(topo)),
                None => frag.contacts(&[], &crate::topology::Topology::default(), 0.0),
            };
            report.contacts += contacts.len();
            // Wood on wood: it does not bounce.
            for (member, impulse) in frag.resolve(&contacts, 0.15) {
                strikes.push((*node, member, impulse));
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
    pub fn settle(&mut self) {
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
    /// bulk luminosity gives a correct answer (about 10^-4 W/m^2) that is
    /// correct precisely because a tree in interstellar space does not grow.
    pub fn environment_at(&self, idx: NodeIdx) -> crate::morph::Environment {
        let n = &self.tree.nodes[idx.get()];
        if let Some(env) = self.environments.get(&n.key) {
            return *env;
        }
        // Illumination from the parent's luminosity at this node's distance —
        // so a structure in the shade of its own node's parent really is in the
        // shade, without a separate lighting system.
        let light = if !n.parent.is_none() {
            let p = &self.tree.nodes[n.parent.get()];
            let d = n.frame.offset.norm().max(p.agg.radius * 0.01).max(1e-6);
            (p.agg.luminosity / (4.0 * std::f64::consts::PI * d * d)).min(1400.0)
        } else {
            crate::morph::Environment::default().light_flux
        };
        crate::morph::Environment {
            light_flux: light,
            temperature: n.agg.temperature,
            water: 1.0,
            crowding: 0.0,
            // A structure can only be built out of matter that is actually
            // available: the node's *unstructured* remainder, not its total.
            // Using the total lets a node grow a structure many times its own
            // mass, because the limit then applies per step rather than
            // cumulatively — a one-kilogram node grew a three-tonne tree.
            reservoir_mass: (n.agg.mass
                - n.morphology.as_ref().map(|m| m.built).unwrap_or(0.0))
            .max(0.0),
            labour: self.labour_rate,
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
                n.agg.add_heat(inf.energy);
                n.agg.momentum += inf.momentum;
            }
            InfluenceKind::Impact | InfluenceKind::UserImpulse => {
                n.agg.momentum += inf.momentum;
                n.agg.add_heat(inf.energy);
            }
            InfluenceKind::Probe => {}
        }
        // The node's procedural detail no longer represents its bulk state.
        let idx = inf.target;
        self.tree.pin(idx);
        self.disturb(idx);
        if !self.tree.nodes[idx.get()].bodies.is_empty() {
            // Distribute the impulse over the existing bodies rather than
            // discarding them — throwing away detail a user is looking at, in
            // response to that user poking it, is the worst possible moment.
            let n = &mut self.tree.nodes[idx.get()];
            let total = n.agg.mass.max(1e-300);
            for b in n.bodies.iter_mut() {
                let f = b.mass / total;
                if b.mass > 0.0 {
                    b.vel += inf.momentum.scale(f / b.mass);
                }
                b.internal_energy += inf.energy * f;
            }
        }
    }

    fn record_histories(&mut self) {
        let depth = self.history_depth;
        let entries: Vec<(PathKey, Snapshot)> = self
            .tree
            .nodes
            .iter()
            .filter(|n| n.alive && n.residency.rank() >= Residency::Causal.rank())
            .map(|n| {
                (
                    n.key,
                    Snapshot {
                        t: self.time,
                        offset: n.frame.offset,
                        velocity: n.frame.velocity,
                        mass: n.agg.mass,
                        luminosity: n.agg.luminosity,
                        temperature: n.agg.temperature,
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
                let old = n.agg.mass;
                n.agg.composition =
                    crate::state::Composition::blend(n.agg.composition, old, composition, mass);
                n.agg.mass += mass;
                n.agg.momentum += velocity.scale(mass);
                n.agg.baryon_number = n.agg.mass * n.agg.composition.nucleons_per_kg();
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
        }
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
        let (key, agg, epoch) = {
            let n = &self.tree.nodes[target.get()];
            (n.key, n.agg, n.epoch)
        };
        let obs = *self.observers.first()?;
        let sep = self.tree.separation(obs.anchor, obs.offset, target, Vec3::ZERO);
        let d = sep.value.norm().max(1e-30);

        let view = match self.histories.get(&key) {
            Some(h) if !h.is_empty() => h.retarded(obs.offset, self.time),
            _ => crate::causal::RetardedView {
                snapshot: Snapshot {
                    t: self.time,
                    offset: sep.value + obs.offset,
                    velocity: self.tree.velocity_from(self.tree.root, target),
                    mass: agg.mass,
                    luminosity: agg.luminosity,
                    temperature: agg.temperature,
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
                Quantity::Temperature => agg.temperature * (1.0 + 1e-3 * s.normal()),
                _ => s.uniform(),
            },
        );
        let _ = fact;

        let mut stream = Stream::at(seed, key.0, epoch, Purpose::PhotonEmission);
        let reading = read(instrument, &obs, &view, &agg, &mut stream);

        // Measurement disturbs. An interferometer deposits real energy, and the
        // engine applies it rather than reporting a free lunch.
        if let Reading::Position { disturbance, .. } = reading {
            if disturbance > 0.0 {
                let n = &mut self.tree.nodes[target.get()];
                n.agg.add_heat(disturbance);
                self.tree.pin(target);
                self.disturb(target);
            }
        }
        Some(reading)
    }

    /// Direct authoring. The one path that can violate conservation — so it
    /// records exactly how much it violated it by.
    fn author(&mut self, target: NodeIdx, property: Property, value: f64) {
        if target.is_none() || !self.tree.nodes[target.get()].alive {
            return;
        }
        let before = self.tree.nodes[target.get()].agg.total_energy();
        let key = {
            let n = &mut self.tree.nodes[target.get()];
            match property {
                Property::Mass => n.agg.mass = value.max(0.0),
                Property::Temperature => n.agg.set_temperature(value),
                Property::Radius => n.agg.radius = value.max(1e-30),
                Property::Charge => n.agg.charge = value,
                Property::Luminosity => n.agg.luminosity = value.max(0.0),
            }
            n.key
        };
        let after = self.tree.nodes[target.get()].agg.total_energy();
        self.audit.push(AuthorEvent {
            key,
            property,
            delta_energy: after - before,
            time: self.time,
        });
        self.tree.pin(target);
        self.tree.bump_epoch(target);
        self.disturb(target);
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
            let (key, agg) = {
                let n = &self.tree.nodes[idx.get()];
                (n.key, n.agg)
            };
            let sep = self.tree.separation(obs.anchor, obs.offset, idx, Vec3::ZERO);
            let d = sep.value.norm();
            if d <= 0.0 || !obs.sees(sep.value) {
                continue;
            }
            let view = match self.histories.get(&key) {
                Some(h) if !h.is_empty() => h.retarded(obs.offset, self.time),
                _ => crate::causal::RetardedView {
                    snapshot: Snapshot {
                        t: self.time,
                        offset: sep.value + obs.offset,
                        velocity: self.tree.velocity_from(self.tree.root, idx),
                        mass: agg.mass,
                        luminosity: agg.luminosity,
                        temperature: agg.temperature,
                    },
                    t_retarded: self.time - d / C,
                    delay: d / C,
                    distance: d,
                    within_history: false,
                },
            };
            let theta = crate::coords::angular_size(agg.radius, d);
            let dir = (view.snapshot.offset - obs.offset).unit();
            let dop = crate::coords::doppler(view.snapshot.velocity - obs.velocity, dir);
            let mut stream = Stream::at(self.tree.world_seed, key.0, 0, Purpose::PhotonEmission);
            let reading = read(instrument, &obs, &view, &agg, &mut stream);
            out.push(Sighting {
                node: idx,
                key,
                view,
                angular_size: theta,
                resolved: theta >= obs.angular_resolution,
                required_tier: obs.required_tier(d),
                flux: flux(agg.luminosity, d, dop),
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

    /// Drill from a node down to the requested tier along the most massive
    /// child at each step, materialising and promoting as it goes.
    ///
    /// This is the "zoom in" primitive: the path from a galaxy to a nucleus is
    /// one call, and the engine materialises exactly the chain of nodes along
    /// the way and nothing else. That chain is a few thousand bodies, not 10^66.
    pub fn drill(&mut self, from: NodeIdx, to_tier: Tier, specs: &dyn Fn(Tier) -> ProlongSpec) -> Vec<NodeIdx> {
        let mut path = vec![from];
        let mut cur = from;
        for _ in 0..64 {
            let tier = self.tree.nodes[cur.get()].tier;
            if tier >= to_tier {
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
                    - na.agg.radius
                    - nb.agg.radius)
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

    let mut agg = Aggregate::neutral(total, radius, 1e4, crate::state::Composition::primordial());
    // A galaxy's internal energy is dominated by orbital motion, not heat:
    // the virial theorem sets it, so derive it rather than inventing a number.
    // The halo dominates the potential, so the velocity dispersion is set by
    // the *total* enclosed mass even though only the baryons are represented.
    let enclosed = total + mass_dark;
    let sigma = (G * enclosed / (2.0 * radius)).sqrt();
    agg.internal_energy = 0.5 * total * sigma * sigma;
    // The baryons' own binding, which refinement can and must reproduce...
    agg.binding_energy = -0.6 * G * total * total / radius;
    // ...and the halo's grip on them, which it cannot, because the halo is not
    // made of the thing being refined.
    agg.external_potential = -G * total * mass_dark / radius;
    // Angular momentum of a rotationally supported disc.
    agg.spin = crate::math::v3(0.0, 0.0, 0.7 * total * sigma * radius);
    agg.luminosity = stars * 3.828e26 * 0.3;

    let spec = ProlongSpec {
        count: 20_000,
        profile: crate::prolong::Profile::Disk {
            scale_height_ratio: 0.12,
        },
        spectrum: crate::prolong::MassSpectrum::Equal,
        kind: crate::state::BodyKind::Super,
        composition_scatter: 0.15,
        turbulent_fraction: 0.3,
    };
    Tree::new(world_seed, agg, Tier::Galactic, spec)
}

/// The default refinement policy, and the budgeted form of it.
///
/// Both live in [`crate::prolong`] because [`crate::tree::Tree::promote`] needs
/// them: a child's tier follows from its *size*, so a caller that guesses the
/// tier from the parent can be wrong, and the node then materialises under a
/// policy meant for a different scale. The tier has to be able to reach for its
/// own policy, and the tier is decided below this module.
pub use crate::prolong::{budgeted_spec, default_spec};
