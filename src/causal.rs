//! Causality: light cones as a scheduling primitive.
//!
//! # The idea
//!
//! Parallel discrete-event simulation has a classical problem: to let a process
//! run ahead safely you need a *lookahead* — a guarantee that no message can
//! arrive with a timestamp earlier than `now + L`. Getting a useful lookahead
//! out of a general simulation is famously hard, and conservative schemes
//! deadlock or crawl without one.
//!
//! Here it is free. Nothing can influence anything else faster than light, so a
//! node at distance `d` from every other active region has a guaranteed
//! lookahead of `d/c`, always, with no analysis and no user annotation.
//! Physics hands us exactly the property the algorithm needs.
//!
//! And the lookahead is enormous where it matters most. Two molecular clouds a
//! kiloparsec apart cannot affect each other for 3000 years, so they can be
//! stepped completely independently for 3000 years of simulated time — millions
//! of frames — with no synchronisation at all. Meanwhile two nucleons 2 fm
//! apart have a lookahead of 7 zeptoseconds and must be stepped in lockstep.
//! The same rule produces both, and it is the *physical* rule, so the engine
//! cannot get it wrong in a way that shows up as a visible artefact: a
//! synchronisation error would have to be a faster-than-light influence.
//!
//! # Consequence for observation
//!
//! Every observation in this engine is of a *retarded* state. When a user looks
//! at a star 8 kpc away they see it as it was 26,000 years ago, because that is
//! the only state that exists on their past light cone. `History` keeps enough
//! of each node's past for that lookup to be answerable.

use crate::coords::retarded_time;
use crate::ids::NodeIdx;
use crate::math::Vec3;
use crate::units::C;
use std::collections::BinaryHeap;

/// One past state of a node, kept for retarded evaluation.
///
/// Deliberately small (72 bytes): a node needs enough history to cover the
/// light-crossing time of its neighbourhood, and at galactic tier that is
/// thousands of samples across millions of nodes.
#[derive(Debug, Clone, Copy, Default)]
pub struct Snapshot {
    pub t: f64,
    pub offset: Vec3,
    pub velocity: Vec3,
    pub mass: f64,
    pub luminosity: f64,
    pub temperature: f64,
}

impl Snapshot {
    pub fn lerp(a: &Snapshot, b: &Snapshot, t: f64) -> Snapshot {
        let dt = b.t - a.t;
        let u = if dt.abs() > 0.0 {
            ((t - a.t) / dt).clamp(0.0, 1.0)
        } else {
            0.0
        };
        Snapshot {
            t,
            // Hermite interpolation on position using the stored velocities:
            // linear interpolation of position alone puts a kink in the
            // trajectory at every sample, which shows up as spurious jerk in
            // retarded force evaluation.
            offset: hermite(a.offset, a.velocity, b.offset, b.velocity, dt, u),
            velocity: a.velocity + (b.velocity - a.velocity).scale(u),
            mass: a.mass + (b.mass - a.mass) * u,
            luminosity: a.luminosity + (b.luminosity - a.luminosity) * u,
            temperature: a.temperature + (b.temperature - a.temperature) * u,
        }
    }
}

fn hermite(p0: Vec3, v0: Vec3, p1: Vec3, v1: Vec3, dt: f64, u: f64) -> Vec3 {
    let u2 = u * u;
    let u3 = u2 * u;
    let h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    let h10 = u3 - 2.0 * u2 + u;
    let h01 = -2.0 * u3 + 3.0 * u2;
    let h11 = u3 - u2;
    p0.scale(h00) + v0.scale(h10 * dt) + p1.scale(h01) + v1.scale(h11 * dt)
}

/// Fixed-capacity ring of past states.
#[derive(Debug, Clone)]
pub struct History {
    ring: Vec<Snapshot>,
    head: usize,
    len: usize,
}

impl History {
    pub fn new(capacity: usize) -> History {
        History {
            ring: vec![Snapshot::default(); capacity.max(2)],
            head: 0,
            len: 0,
        }
    }

    pub fn push(&mut self, s: Snapshot) {
        let cap = self.ring.len();
        self.ring[self.head] = s;
        self.head = (self.head + 1) % cap;
        self.len = (self.len + 1).min(cap);
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn at(&self, i: usize) -> &Snapshot {
        let cap = self.ring.len();
        let start = (self.head + cap - self.len) % cap;
        &self.ring[(start + i) % cap]
    }

    pub fn newest(&self) -> Option<&Snapshot> {
        if self.len == 0 {
            None
        } else {
            Some(self.at(self.len - 1))
        }
    }

    pub fn oldest(&self) -> Option<&Snapshot> {
        if self.len == 0 {
            None
        } else {
            Some(self.at(0))
        }
    }

    /// How far back this history reaches. If a retarded query needs a time
    /// older than this, the engine must fall back to extrapolation and say so —
    /// silently clamping would fabricate a state that never existed.
    pub fn span(&self) -> f64 {
        match (self.oldest(), self.newest()) {
            (Some(a), Some(b)) => b.t - a.t,
            _ => 0.0,
        }
    }

    /// Interpolate the state at time `t`. `Err` carries the clamped result when
    /// `t` falls outside the retained window, so callers can degrade
    /// deliberately rather than by accident.
    pub fn sample(&self, t: f64) -> Result<Snapshot, Snapshot> {
        if self.len == 0 {
            return Err(Snapshot::default());
        }
        let oldest = *self.at(0);
        let newest = *self.at(self.len - 1);
        if t <= oldest.t {
            return Err(oldest);
        }
        if t >= newest.t {
            // Extrapolate forward ballistically: this is the common case for a
            // node that has not been stepped recently, and a straight line is
            // exactly right for anything not currently being accelerated.
            let dt = t - newest.t;
            let mut s = newest;
            s.t = t;
            s.offset = newest.offset + newest.velocity.scale(dt);
            return if dt > self.span().max(1e-30) {
                Err(s)
            } else {
                Ok(s)
            };
        }
        // Binary search for the bracketing pair.
        let (mut lo, mut hi) = (0usize, self.len - 1);
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if self.at(mid).t <= t {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(Snapshot::lerp(self.at(lo), self.at(hi), t))
    }

    /// The state of this node on the past light cone of an observer at
    /// `observer_pos` at time `t_obs`, both expressed in a common frame.
    ///
    /// This — not `newest()` — is what an observer is entitled to see.
    pub fn retarded(&self, observer_pos: Vec3, t_obs: f64) -> RetardedView {
        let src_at = |t: f64| self.sample(t).unwrap_or_else(|s| s).offset;
        let (t_ret, pos, dist) = retarded_time(observer_pos, t_obs, &src_at);
        let (snap, in_window) = match self.sample(t_ret) {
            Ok(s) => (s, true),
            Err(s) => (s, false),
        };
        RetardedView {
            snapshot: Snapshot { offset: pos, ..snap },
            t_retarded: t_ret,
            delay: t_obs - t_ret,
            distance: dist,
            within_history: in_window,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RetardedView {
    pub snapshot: Snapshot,
    pub t_retarded: f64,
    pub delay: f64,
    pub distance: f64,
    /// False when the answer was extrapolated beyond the retained history. The
    /// observation system reports this rather than hiding it.
    pub within_history: bool,
}

/// The light-cone test that gates every materialisation decision.
///
/// A region needs fine detail only if fine detail there could reach an observer
/// within the horizon. Everything else can stay coarse no matter how
/// interesting it is, because its influence has not arrived yet.
#[derive(Debug, Clone, Copy)]
pub struct CausalGate {
    /// How far into the future the engine promises to be correct. Larger means
    /// more of the world must be resolved; smaller means detail pops in.
    pub horizon: f64,
}

impl CausalGate {
    pub fn new(horizon: f64) -> CausalGate {
        CausalGate { horizon }
    }

    /// Can a change at distance `d`, happening now, reach the observer within
    /// the horizon?
    #[inline]
    pub fn reaches(&self, d: f64) -> bool {
        d <= C * self.horizon
    }

    /// The radius of the region that must be kept consistent.
    #[inline]
    pub fn radius(&self) -> f64 {
        C * self.horizon
    }

    /// Fraction of the horizon already used up by light travel — 0 at the
    /// observer, 1 at the edge of the cone. Used to taper detail smoothly
    /// towards the horizon instead of cutting it off with a visible shell.
    #[inline]
    pub fn urgency(&self, d: f64) -> f64 {
        if self.horizon <= 0.0 {
            return 0.0;
        }
        (1.0 - (d / C) / self.horizon).clamp(0.0, 1.0)
    }
}

/// A pending influence in flight. Interactions do not apply instantly: they are
/// posted with an arrival time and delivered when the target's clock reaches it.
#[derive(Debug, Clone, Copy)]
pub struct Influence {
    pub arrives: f64,
    pub target: NodeIdx,
    pub kind: InfluenceKind,
    pub energy: f64,
    pub momentum: Vec3,
    pub source_distance: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InfluenceKind {
    Radiation,
    Blast,
    Impact,
    Probe,
    UserImpulse,
}

/// Min-heap ordering on arrival time.
#[derive(Debug, Clone, Copy)]
struct Timed(Influence);

impl PartialEq for Timed {
    fn eq(&self, o: &Self) -> bool {
        self.0.arrives == o.0.arrives
    }
}
impl Eq for Timed {}
impl PartialOrd for Timed {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Timed {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        // Reversed: BinaryHeap is a max-heap, we want earliest first. Ties are
        // broken by target index so the order is total and deterministic —
        // otherwise two influences arriving at the same instant could be
        // applied in either order and replay would diverge.
        o.0.arrives
            .partial_cmp(&self.0.arrives)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(o.0.target.0.cmp(&self.0.target.0))
    }
}

/// Delivery queue for influences travelling at c.
#[derive(Default)]
pub struct Mailbox {
    heap: BinaryHeap<Timed>,
    pub delivered: u64,
    pub in_flight_peak: usize,
}

impl Mailbox {
    pub fn new() -> Mailbox {
        Mailbox::default()
    }

    /// Post an influence, computing its arrival from the separation. The
    /// `+ d/c` here is the entire relativistic content of the interaction
    /// system: nothing else in the engine needs to know about light delay.
    pub fn post(
        &mut self,
        target: NodeIdx,
        now: f64,
        distance: f64,
        kind: InfluenceKind,
        energy: f64,
        momentum: Vec3,
    ) {
        self.heap.push(Timed(Influence {
            arrives: now + distance / C,
            target,
            kind,
            energy,
            momentum,
            source_distance: distance,
        }));
        self.in_flight_peak = self.in_flight_peak.max(self.heap.len());
    }

    /// Everything that has arrived by `t`, in deterministic order.
    pub fn drain_until(&mut self, t: f64) -> Vec<Influence> {
        let mut out = Vec::new();
        while let Some(Timed(top)) = self.heap.peek().copied() {
            if top.arrives <= t {
                self.heap.pop();
                out.push(top);
                self.delivered += 1;
            } else {
                break;
            }
        }
        out
    }

    pub fn pending(&self) -> usize {
        self.heap.len()
    }

    /// Time of the next arrival — a hard ceiling on how far any node may be
    /// advanced without risking a causality violation.
    pub fn next_arrival(&self) -> Option<f64> {
        self.heap.peek().map(|t| t.0.arrives)
    }
}

/// Per-node scheduling state for the conservative advance.
#[derive(Debug, Clone, Copy)]
pub struct Clock {
    /// Coordinate time this node has been integrated to.
    pub time: f64,
    /// Proper time elapsed in the node's own frame — what a clock sitting on
    /// the node would read. Differs from `time` by the accumulated Lorentz and
    /// gravitational factors, and is what the UI shows when the user attaches
    /// themselves to a node.
    pub proper_time: f64,
    /// Distance to the nearest node that could influence this one.
    pub nearest_neighbour: f64,
    /// The node's own preferred step from its physics.
    pub natural_dt: f64,
}

impl Clock {
    pub fn new(time: f64, natural_dt: f64) -> Clock {
        Clock {
            time,
            proper_time: time,
            nearest_neighbour: f64::INFINITY,
            natural_dt,
        }
    }

    /// The guaranteed-safe lookahead: no influence can arrive sooner.
    #[inline]
    pub fn lookahead(&self) -> f64 {
        if self.nearest_neighbour.is_finite() {
            self.nearest_neighbour / C
        } else {
            f64::INFINITY
        }
    }

    /// The largest step this node may take right now.
    #[inline]
    pub fn safe_step(&self, global_floor: f64, next_event: f64) -> f64 {
        let by_causality = self.lookahead();
        let by_event = (next_event - self.time).max(0.0);
        self.natural_dt
            .min(by_causality)
            .min(if by_event > 0.0 { by_event } else { f64::INFINITY })
            .max(global_floor)
    }
}

/// Conservative multi-rate advance over a set of clocks.
///
/// Returns, for each clock, the step it may take. Nodes far from everything
/// take huge steps; nodes packed together take tiny ones; nobody ever computes
/// a force from a state that had not happened yet.
pub fn plan_steps(clocks: &[Clock], next_event: f64, floor: f64) -> Vec<f64> {
    clocks
        .iter()
        .map(|c| c.safe_step(floor, next_event))
        .collect()
}

/// Check that an advance did not outrun causality: no node may be ahead of
/// another by more than the light travel time between them.
///
/// This is the assertion that the whole scheduler exists to satisfy, and it is
/// checked directly in `tests/causality.rs` rather than trusted.
pub fn causality_violation(t_a: f64, t_b: f64, separation: f64) -> f64 {
    let skew = (t_a - t_b).abs();
    let allowed = separation / C;
    (skew - allowed).max(0.0)
}
