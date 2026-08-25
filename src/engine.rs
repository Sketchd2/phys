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
}

/// The world.
pub struct World {
    pub tree: Tree,
    pub mailbox: Mailbox,
    pub ledger: Ledger,
    pub observers: Vec<Observer>,
    pub budget: FrameBudget,
    pub gate: CausalGate,
    /// Coordinate time at the root, seconds since the scenario epoch.
    pub time: f64,
    /// Simulated seconds per wall-clock second. The user's time control.
    pub time_rate: f64,
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
        World {
            tree,
            mailbox: Mailbox::new(),
            ledger: Ledger::new(),
            observers: Vec::new(),
            budget: FrameBudget::ups(ups),
            gate: CausalGate::new(1e3 * YEAR),
            time: 0.0,
            time_rate: 1.0,
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
        }
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

    /// Advance the world by one frame's worth of simulated time.
    ///
    /// `wall_us` is what the caller is prepared to spend. The engine will spend
    /// up to that and no more; if the work does not fit, detail is dropped and
    /// `stats.detail_debt` records how much value went unserved.
    pub fn step_frame(&mut self, wall_us: f64) -> Plan {
        let t0 = std::time::Instant::now();
        self.budget.target_us = wall_us;

        let tasks = self.survey();
        let bytes = self.tree.detail_bytes();
        let plan = self.budget.plan(tasks, bytes);
        self.execute(&plan);

        let dt = self.frame_dt();
        self.deliver_influences(self.time + dt);
        self.time += dt;
        self.record_histories();

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
        plan
    }

    /// How much simulated time one frame covers.
    ///
    /// Not a free choice: the timestep is fixed by accuracy (`docs/PHYSICS.md`
    /// shows the dt^2 convergence), so what the frame budget really controls is
    /// how many steps are affordable, and therefore how fast simulated time
    /// advances. Zooming into a nucleus does not slow the frame rate — it slows
    /// *time*, which is both physically evocative and the honest thing to do.
    pub fn frame_dt(&self) -> f64 {
        let mut dt = f64::INFINITY;
        for node in self.tree.nodes.iter().filter(|n| n.alive) {
            let natural = node.tier.dt();
            let dyn_t = node.agg.dynamical_time();
            let by_physics = natural.min(dyn_t / 50.0);
            dt = dt.min(by_physics);
        }
        if !dt.is_finite() {
            dt = 1.0;
        }
        // The causal ceiling: no node may be advanced past the next influence
        // arrival, or the influence would land in its past.
        if let Some(next) = self.mailbox.next_arrival() {
            dt = dt.min((next - self.time).max(0.0).max(1e-30));
        }
        dt * self.time_rate
    }

    /// Build the frame's candidate work list.
    fn survey(&mut self) -> Vec<Task> {
        let mut tasks = Vec::new();
        let live: Vec<NodeIdx> = (0..self.tree.nodes.len())
            .map(|i| NodeIdx(i as u32))
            .filter(|i| self.tree.nodes[i.get()].alive)
            .collect();

        for idx in live {
            let (tier, materialised, count, _key, radius) = {
                let n = &self.tree.nodes[idx.get()];
                (n.tier, n.is_materialised(), n.bodies.len(), n.key, n.agg.radius)
            };

            // How much does anyone care about this node?
            let mut salience = 0.0f64;
            let mut urgency = 0.0f64;
            let mut wanted_tier = tier;
            for obs in &self.observers {
                let sep = self
                    .tree
                    .separation(obs.anchor, obs.offset, idx, Vec3::ZERO);
                let d = sep.value.norm().max(1e-30);
                if !self.gate.reaches(d) {
                    // Outside the light cone: whatever happens here cannot
                    // reach the observer within the horizon, so it does not
                    // need resolving no matter how interesting it is.
                    continue;
                }
                let theta = crate::coords::angular_size(radius, d);
                let s = (theta / obs.angular_resolution).min(1e6) * obs.priority;
                if s > salience {
                    salience = s;
                    wanted_tier = obs.required_tier(d);
                }
                urgency = urgency.max(self.gate.urgency(d));
            }
            let residency = if salience > 1.0 {
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

            // Materialise when an observer needs finer structure than the bulk
            // state can express.
            if !materialised && wanted_tier > tier && salience > 0.5 {
                let n_children = self.tree.nodes[idx.get()].spec.count;
                tasks.push(Task {
                    node: idx,
                    kind: TaskKind::Materialise,
                    cost_us: cost::materialise_us(n_children),
                    salience,
                    urgency: urgency.max(0.01),
                    error,
                    novelty: 1.0,
                    bytes: (n_children * std::mem::size_of::<Body>()) as i64,
                });
            }

            // Coarsen when nobody needs the detail. Negative bytes: this task
            // gives resources back, so the planner always accepts it.
            if materialised && salience < 0.25 && !self.tree.nodes[idx.get()].pinned {
                tasks.push(Task {
                    node: idx,
                    kind: TaskKind::Coarsen,
                    cost_us: cost::COARSEN_US * count as f64,
                    salience: 1.0,
                    urgency: 1.0,
                    error: 1.0,
                    novelty: 0.0,
                    bytes: -((count * std::mem::size_of::<Body>()) as i64),
                });
            }

            // Record the causal lookahead this node is entitled to, from its
            // nearest sibling. Nodes with no siblings are causally isolated
            // within their parent and may run as fast as their physics allows.
            {
                let parent = self.tree.nodes[idx.get()].parent;
                if !parent.is_none() {
                    for (child, gap) in self.tree.sibling_separations(parent) {
                        let key = self.tree.nodes[child.get()].key;
                        let tier = self.tree.nodes[child.get()].tier;
                        let t = self.tree.nodes[child.get()].time;
                        let c = self
                            .clocks
                            .entry(key)
                            .or_insert_with(|| Clock::new(t, tier.dt()));
                        c.nearest_neighbour = gap;
                    }
                }
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
                    salience: salience.max(1.0),
                    urgency: 1.0,
                    error: 1.0,
                    novelty: 0.0,
                    bytes: 0,
                });
            }

            // Step whatever is materialised.
            if materialised {
                let kind = solvers::for_tier(tier);
                tasks.push(Task {
                    node: idx,
                    kind: TaskKind::Step,
                    cost_us: cost::step_us(kind, count, self.gpu),
                    salience: salience.max(0.05),
                    urgency: urgency.max(0.05),
                    error,
                    novelty: 0.1,
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

    fn execute(&mut self, plan: &Plan) {
        for task in &plan.accepted {
            if !self.tree.nodes[task.node.get()].alive {
                continue;
            }
            match task.kind {
                TaskKind::Materialise => {
                    self.tree.refine(task.node);
                }
                TaskKind::Coarsen => {
                    self.tree.coarsen(task.node);
                }
                TaskKind::Promote => {}
                TaskKind::Step => {
                    let dt = self.node_dt(task.node);
                    self.advance_node(task.node, dt);
                }
                TaskKind::Grow => {
                    let dt = self.frame_dt();
                    self.grow_node(task.node, dt);
                }
                TaskKind::Observe => {}
            }
        }
    }

    /// The step a node may take: the smaller of what its physics wants and what
    /// causality permits.
    pub fn node_dt(&self, idx: NodeIdx) -> f64 {
        let n = &self.tree.nodes[idx.get()];
        let natural = n.tier.dt().min(n.agg.dynamical_time() / 50.0);
        let clock = self
            .clocks
            .get(&n.key)
            .copied()
            .unwrap_or_else(|| Clock::new(n.time, natural));
        let next = self.mailbox.next_arrival().unwrap_or(f64::INFINITY);
        clock.safe_step(1e-30, next).min(natural) * self.time_rate
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
                solvers::md::step(bodies, dt, params, seed, key.0, epoch, tick)
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
            }
            Interaction::Measure {
                target,
                instrument,
                quantity,
            } => {
                self.measure(target, instrument, quantity);
            }
            Interaction::Pin { target } => self.tree.pin(target),
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

/// The default refinement policy: how each tier splits into the next.
///
/// This is the engine's "what is the universe made of" table. Each entry says
/// how many children a node of that tier has, how they are arranged, and how
/// mass is divided among them. Changing a line here changes the character of
/// the world at that scale and nothing else — the conservation machinery,
/// scheduler and observation path are all indifferent to it.
pub fn default_spec(tier: Tier) -> ProlongSpec {
    use crate::prolong::{MassSpectrum, Profile};
    use crate::state::BodyKind;
    match tier {
        // A galaxy resolves into star-forming complexes and dark matter.
        Tier::Galactic => ProlongSpec {
            count: 20_000,
            profile: Profile::Disk { scale_height_ratio: 0.12 },
            spectrum: MassSpectrum::Equal,
            kind: BodyKind::Super,
            composition_scatter: 0.15,
            turbulent_fraction: 0.3,
        },
        // A complex resolves into individual stars, drawn from the IMF — which
        // is where "a statistical stand-in" becomes "a particular star".
        Tier::Stellar => ProlongSpec {
            count: 4_000,
            profile: Profile::Plummer,
            spectrum: MassSpectrum::Kroupa { min_msun: 0.08, max_msun: 60.0 },
            kind: BodyKind::Star,
            composition_scatter: 0.05,
            turbulent_fraction: 0.45,
        },
        // A star resolves into its structural shells; a planet into its layers.
        Tier::Planetary => ProlongSpec {
            count: 2_000,
            profile: Profile::Plummer,
            spectrum: MassSpectrum::PowerLaw { alpha: -1.5, ratio: 30.0 },
            kind: BodyKind::GasParcel,
            composition_scatter: 0.02,
            turbulent_fraction: 0.6,
        },
        // Bulk matter resolves into fluid parcels or grains.
        Tier::Continuum => ProlongSpec {
            count: 4_000,
            profile: Profile::Uniform,
            spectrum: MassSpectrum::Equal,
            kind: BodyKind::Grain,
            composition_scatter: 0.01,
            turbulent_fraction: 0.2,
        },
        // A parcel resolves into molecules.
        Tier::Molecular => ProlongSpec {
            count: 8_000,
            profile: Profile::Uniform,
            spectrum: MassSpectrum::Species,
            kind: BodyKind::Molecule,
            composition_scatter: 0.0,
            turbulent_fraction: 0.0,
        },
        // A molecule resolves into atoms.
        Tier::Atomic => ProlongSpec {
            count: 64,
            profile: Profile::Lattice,
            spectrum: MassSpectrum::Species,
            kind: BodyKind::Atom,
            composition_scatter: 0.0,
            turbulent_fraction: 0.0,
        },
        // An atom resolves into a nucleus of nucleons. Below this the engine
        // stops producing trajectories and switches to the statistical
        // description in `solvers::quantum` — not because it runs out of
        // compute, but because there is nothing else there to describe.
        Tier::Nuclear => ProlongSpec {
            count: 56,
            profile: Profile::WoodsSaxon,
            spectrum: MassSpectrum::Equal,
            kind: BodyKind::Nucleon,
            composition_scatter: 0.0,
            turbulent_fraction: 0.0,
        },
    }
}
