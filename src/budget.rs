//! The frame budget: how 20 updates per second becomes a guarantee instead of
//! a hope.
//!
//! # The inversion
//!
//! A conventional simulation decides what to compute and then takes however
//! long it takes. This one is given 50 milliseconds and decides what fits.
//! Every candidate piece of work — materialise this cloud, step that cluster,
//! resolve those molecules — is a `Task` with an estimated cost and an
//! estimated *value*, and each frame the scheduler solves a knapsack against
//! the wall clock.
//!
//! The consequence is that the frame rate never degrades. What degrades first
//! is the *rate of simulated time*, and then detail, in the order that leaves
//! the world least stale — the most overdue work first. A user who zooms too
//! far does not get a slideshow; they get a world running in slow motion with
//! a visible "detail debt" indicator, which is a far better failure mode.
//!
//! # Value
//!
//! ```text
//!   value = lateness * urgency * error
//! ```
//!
//! * **lateness** — how long since the node was last solved, measured in its
//!   own characteristic times. One means the node has just come due; ten means
//!   it is nine cadences stale. This is the term that decides nearly every
//!   frame, and it is self-correcting: a node passed over grows more overdue
//!   and outranks its neighbours next frame, so nothing can be starved.
//! * **urgency** — how close the work is to the causal horizon. Something whose
//!   influence arrives next frame outranks something whose influence arrives in
//!   an hour, even if the latter is closer.
//! * **error** — the estimated inaccuracy of *not* doing the work. A node in
//!   free fall with a resolved Jeans length is cheap to skip; one about to
//!   fragment is not.
//!
//! Value used to be led by observer salience — how much of the screen the work
//! subtended. That put the camera in charge of the physics, which is the wrong
//! way round: a tree falls whether or not anyone is pointing at it, and a
//! reactor that only runs while it is on screen is not a simulation. What the
//! camera legitimately controls is *resolution*, and resolution feeds back into
//! lateness on its own, because a node represented more finely has a shorter
//! characteristic time and therefore comes due more often.

use crate::ids::NodeIdx;

/// The unit of schedulable work.
#[derive(Debug, Clone, Copy)]
pub struct Task {
    pub node: NodeIdx,
    pub kind: TaskKind,
    /// Estimated cost in microseconds, from the calibrated model.
    pub cost_us: f64,
    /// How many of this node's own characteristic times have passed since it
    /// was last solved. Below one it is not yet due.
    pub lateness: f64,
    pub urgency: f64,
    pub error: f64,
    /// Bytes of working set this task will add. Checked against the VRAM cap
    /// separately from time — on a 6 GB card memory runs out before time does.
    pub bytes: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TaskKind {
    Materialise,
    Step,
    Coarsen,
    Promote,
    Observe,
    /// Advance a structure's developmental state. Runs on the aggregate, never
    /// on the structure, so it costs O(1) per node no matter how elaborate the
    /// thing being grown is.
    Grow,
}

impl Task {
    #[inline]
    pub fn value(&self) -> f64 {
        self.lateness * self.urgency * self.error
    }

    /// Value per microsecond — the quantity a greedy knapsack sorts on.
    #[inline]
    pub fn density(&self) -> f64 {
        let c = self.cost_us.max(1e-6);
        self.value() / c
    }
}

/// Outcome of planning one frame.
#[derive(Debug, Clone, Default)]
pub struct Plan {
    pub accepted: Vec<Task>,
    pub deferred: usize,
    pub planned_us: f64,
    pub planned_bytes: i64,
    /// Microseconds the deferred work would have cost. With `planned_us` this
    /// is the frame's total *demand*, which is what the engine needs in order
    /// to decide how much simulated time it can afford next frame.
    pub deferred_us: f64,
    /// Total value the frame had to leave on the table. This is the engine's
    /// honest "how far behind am I" number, and it is what the UI shows as
    /// detail debt.
    pub unmet_value: f64,
    /// Set when the frame accepted work it could not afford, because the
    /// cheapest indivisible piece of work cost more than the whole budget.
    pub overrun: bool,
}

/// A frame budget with a self-calibrating cost model.
///
/// The cost model starts from the constants in `docs/PERFORMANCE.md` and then
/// corrects itself from measurements, because the true cost per body depends on
/// cache behaviour, on the particular machine, and on what else is resident.
/// A fixed cost model is wrong within minutes of leaving the developer's
/// machine.
pub struct FrameBudget {
    /// Wall-clock target per frame, microseconds. 50,000 for 20 UPS.
    pub target_us: f64,
    /// Fraction of the frame the simulation may use; the rest is rendering,
    /// input and slack.
    pub sim_fraction: f64,
    /// Hard cap on materialised detail, bytes.
    ///
    /// `u64` rather than `usize` because the cap describes a *device* budget,
    /// not the host's address space: a 32-bit target (wasm32) cannot even name
    /// four gigabytes, and the constant overflowed at compile time there.
    pub memory_cap: u64,
    /// Multiplicative correction learned from measured frames.
    calibration: f64,
    /// Exponential moving average of the achieved frame time.
    pub measured_us: f64,
    pub frames: u64,
    pub overruns: u64,
}

impl FrameBudget {
    /// 20 updates per second, with 70% of each frame available to physics.
    pub fn ups(ups: f64) -> FrameBudget {
        FrameBudget {
            target_us: 1e6 / ups,
            sim_fraction: 0.7,
            memory_cap: 4 * 1024 * 1024 * 1024,
            calibration: 1.0,
            measured_us: 0.0,
            frames: 0,
            overruns: 0,
        }
    }

    #[inline]
    pub fn sim_budget_us(&self) -> f64 {
        self.target_us * self.sim_fraction
    }

    /// Correct a raw cost estimate by the learned calibration.
    #[inline]
    pub fn adjust(&self, raw_us: f64) -> f64 {
        raw_us * self.calibration
    }

    /// Feed back the measured frame time and update the calibration.
    ///
    /// Deliberately asymmetric: it reacts fast to overruns and slowly to
    /// underruns. Being late is visible to the user and being early is not, so
    /// the controller should be eager to back off and reluctant to push.
    pub fn observe_frame(&mut self, planned_us: f64, actual_us: f64) {
        self.frames += 1;
        self.measured_us = if self.frames == 1 {
            actual_us
        } else {
            0.9 * self.measured_us + 0.1 * actual_us
        };
        if actual_us > self.target_us {
            self.overruns += 1;
        }
        if planned_us > 1.0 && actual_us > 0.0 {
            let ratio = actual_us / planned_us;
            let rate = if ratio > 1.0 { 0.35 } else { 0.05 };
            self.calibration = (self.calibration * (1.0 - rate) + ratio * self.calibration * rate)
                .clamp(0.05, 20.0);
        }
    }

    pub fn calibration(&self) -> f64 {
        self.calibration
    }

    /// Greedy knapsack by value density.
    ///
    /// Greedy rather than exact: the exact 0/1 knapsack is NP-hard and the
    /// greedy solution is within a factor of 2 of optimal, which is far inside
    /// the error of the cost estimates themselves. Spending frame time to plan
    /// the frame better than the estimates justify would be self-defeating.
    ///
    /// Tasks are sorted with a deterministic tie-break on node index so that
    /// two runs of the same scenario schedule identically — a scheduler that
    /// depends on hash order would silently break replay.
    pub fn plan(&self, mut tasks: Vec<Task>, current_bytes: usize) -> Plan {
        tasks.sort_by(|a, b| {
            b.density()
                .partial_cmp(&a.density())
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.node.0.cmp(&b.node.0))
                .then((a.kind as u8).cmp(&(b.kind as u8)))
        });
        let budget = self.sim_budget_us();
        let mut plan = Plan::default();
        let mut bytes = current_bytes as i64;
        let cap = self.memory_cap as i64;
        // Highest density first, so this is the one to fall back on if nothing
        // fits at all.
        let first_choice = tasks.first().copied();

        for t in tasks {
            let cost = self.adjust(t.cost_us);
            // Freeing work (coarsening) is always accepted: it costs little and
            // buys back the resource the next frame will need. Growth is too —
            // deferring it would make simulated time run at different rates for
            // different structures depending on what the camera happened to be
            // pointing at, which is a far worse artefact than a late frame.
            let frees = t.bytes < 0 || t.kind == TaskKind::Grow;
            if !frees && plan.planned_us + cost > budget {
                plan.deferred += 1;
                plan.deferred_us += cost;
                plan.unmet_value += t.value();
                continue;
            }
            if !frees && bytes + t.bytes > cap {
                plan.deferred += 1;
                plan.deferred_us += cost;
                plan.unmet_value += t.value();
                continue;
            }
            plan.planned_us += cost;
            bytes += t.bytes;
            plan.planned_bytes += t.bytes;
            plan.accepted.push(t);
        }
        // A frame that accepts nothing is worse than a frame that runs late.
        //
        // Work is indivisible at this level — a solver pass over a node is one
        // pass or none — so a world whose cheapest node costs more than a whole
        // frame would defer every task, every frame, forever, and stand
        // perfectly still at a perfect twenty updates a second. Take the best
        // one and overrun, and say so.
        if plan.accepted.is_empty() {
            if let Some(best) = first_choice {
                plan.planned_us += self.adjust(best.cost_us);
                plan.planned_bytes += best.bytes;
                plan.deferred = plan.deferred.saturating_sub(1);
                plan.deferred_us = (plan.deferred_us - self.adjust(best.cost_us)).max(0.0);
                plan.unmet_value = (plan.unmet_value - best.value()).max(0.0);
                plan.overrun = true;
                plan.accepted.push(best);
            }
        }
        plan
    }
}

/// Cost model constants, in microseconds.
///
/// Anchored to the single-core measurements in `docs/PERFORMANCE.md`: a
/// Barnes-Hut force evaluation costs about 11.5 us per body per step on one
/// Ryzen 5 3600 core at n = 4000. The GPU factors are the ratios assumed for
/// the RTX 2060 path and are stated as assumptions, not measurements.
pub mod cost {
    /// Per body, per step, CPU single core.
    pub const GRAVITY_STEP_US: f64 = 11.5;
    pub const HYDRO_STEP_US: f64 = 3.2;
    pub const MD_STEP_US: f64 = 1.8;
    pub const STATISTICAL_STEP_US: f64 = 0.05;
    /// Per node, per growth step. One ODE evaluation on an aggregate.
    pub const GROW_US: f64 = 0.12;
    /// Per body produced.
    pub const MATERIALISE_US: f64 = 0.35;
    pub const COARSEN_US: f64 = 0.08;
    /// Speedup assumed for the RTX 2060 path over one CPU core.
    pub const GPU_SPEEDUP: f64 = 60.0;

    pub fn step_us(kind: crate::solvers::SolverKind, n: usize, gpu: bool) -> f64 {
        let per = match kind {
            crate::solvers::SolverKind::Gravity => GRAVITY_STEP_US,
            crate::solvers::SolverKind::GravityHydro => GRAVITY_STEP_US + HYDRO_STEP_US,
            crate::solvers::SolverKind::Hydro => HYDRO_STEP_US,
            crate::solvers::SolverKind::MolecularDynamics => MD_STEP_US,
            crate::solvers::SolverKind::Statistical => STATISTICAL_STEP_US,
        };
        let raw = per * n as f64;
        if gpu {
            raw / GPU_SPEEDUP
        } else {
            raw
        }
    }

    pub fn materialise_us(n: usize) -> f64 {
        MATERIALISE_US * n as f64
    }
}
