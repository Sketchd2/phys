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
            history_depth: 64,
        }
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
