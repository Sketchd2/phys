//! How fast a node's own clock runs.
//!
//! # What was wrong
//!
//! The engine computed proper time and then threw it away. `Frame::advance`
//! accumulated `dt/γ` into `proper_time`, `Clock` carried the number, the UI
//! could display it — and nothing ever *read* it. Every node's chemistry,
//! growth, decay and dynamics were stepped with coordinate `dt` regardless of
//! how fast the node was moving or how deep a well it sat in.
//! `coords::gravitational_dilation` was written, documented, and called from
//! nowhere at all.
//!
//! So a muon crossing the atmosphere at 0.999c decayed on the atmosphere's
//! clock instead of its own, which is the textbook case relativity exists to
//! get right, and the engine had every ingredient for it and used none of them.
//!
//! # The split that fixes it
//!
//! A node has two motions and they answer to different clocks:
//!
//! * **Its trajectory through its parent's frame** — where it is, which way it
//!   points — advances on *coordinate* time. That is what "the parent watched
//!   it move" means, and it is what `Frame::advance` has always done.
//! * **Its interior** — what it is doing, burning, growing, becoming —
//!   advances on *local* time, which is coordinate time multiplied by the rate
//!   below.
//!
//! Get that split wrong in either direction and something visible breaks: tie
//! the interior to coordinate time and fast things age too quickly; tie the
//! trajectory to local time and a fast object appears to crawl.
//!
//! # One rule, three sources
//!
//! ```text
//!     rate(node) = product over the node and every ancestor of
//!                      (1/γ) · √(1 + 2Φ/c²) · bubble
//! ```
//!
//! The product is not an approximation, it is the chain rule. Each node's frame
//! velocity is measured in its *parent's* frame, so `dτ_child/dt_parent` is a
//! per-level quantity and the ratio to the root's clock is what you get by
//! multiplying along the chain — which is also why nesting is free and a
//! nucleon twelve frames deep costs twelve multiplies.
//!
//! Two of the three factors are physics. The third is not:
//!
//! * `1/γ` — special relativity, from the node's velocity in its parent frame.
//! * `√(1 + 2Φ/c²)` — gravitational dilation, to first post-Newtonian order,
//!   from the potential the node sits in.
//! * `bubble` — an administrator saying "run this faster than the universe".
//!
//! Keeping the third in the same product as the first two is deliberate. There
//! is exactly one place that decides how fast a node's interior runs, so a
//! solver cannot accidentally respect one kind of dilation and ignore another;
//! and because [`TimeRate`] reports the three separately, nothing can lose
//! track of which part was physics and which part was somebody's thumb on the
//! scale.
//!
//! # What a bubble deliberately does not do
//!
//! It does not move the object faster. A bubble scales the interior only, and
//! the reason is not caution — it is that *external* speed-up is already
//! available and is called "going faster". An object that traverses the world
//! at a hundred times its velocity is an object with a hundred times the
//! velocity, which the engine models, and which relativity then dilates. What
//! velocity cannot buy is a tree that ages a century while the world watches
//! for a year, and that is precisely what the bubble is for.
//!
//! It also keeps the causal structure intact. `causal.rs` gates influence
//! delivery on light travel time in coordinate seconds; a node whose *position*
//! advanced at a hundred times the rate would outrun the influences it had
//! already emitted, and the mailbox would deliver its past to its own future.

use crate::math::Vec3;
use crate::units::C;

/// A bubble factor outside this range is refused.
///
/// Not a physics limit — a numerics one. The interior's sub-step shrinks in
/// proportion to the rate, so a rate of 10^12 asks for sub-steps that underflow
/// against the world clock, and the node would stall rather than run fast. The
/// bound is wide enough for anything a balancing pass wants (a century in an
/// afternoon is about 10^4) and narrow enough that the arithmetic stays sound.
pub const MAX_BUBBLE: f64 = 1e6;
pub const MIN_BUBBLE: f64 = 1e-6;

/// How fast one node's interior runs, and why.
///
/// Reported as three factors rather than one number so that "this reactor ran
/// hot" and "somebody sped this reactor up" never look the same in a log.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimeRate {
    /// Special-relativistic, `1/γ`. At most 1, and 1 for anything slow.
    pub kinematic: f64,
    /// Gravitational, `√(1 + 2Φ/c²)`. At most 1, and 1 far from mass.
    pub gravitational: f64,
    /// Administrative. Exactly 1 unless somebody set it.
    pub bubble: f64,
}

impl Default for TimeRate {
    fn default() -> TimeRate {
        TimeRate { kinematic: 1.0, gravitational: 1.0, bubble: 1.0 }
    }
}

impl TimeRate {
    /// Local seconds per coordinate second.
    pub fn total(&self) -> f64 {
        let r = self.kinematic * self.gravitational * self.bubble;
        if r.is_finite() && r > 0.0 {
            r
        } else {
            // A rate of zero or NaN would freeze or poison a node rather than
            // slow it. Physics cannot produce either — `gamma` saturates and
            // `gravitational_dilation` floors — so this only catches arithmetic
            // that has already gone wrong, and it fails to "normal".
            1.0
        }
    }

    /// The part relativity is responsible for.
    pub fn physical(&self) -> f64 {
        self.kinematic * self.gravitational
    }

    /// True when nobody has interfered, so a caller can assert a world is
    /// still describing physics.
    pub fn is_physical(&self) -> bool {
        self.bubble == 1.0
    }

    /// Fold another frame's contribution in, walking up the tree.
    pub fn compose(self, o: TimeRate) -> TimeRate {
        TimeRate {
            kinematic: self.kinematic * o.kinematic,
            gravitational: self.gravitational * o.gravitational,
            bubble: self.bubble * o.bubble,
        }
    }
}

/// One frame's physical contribution: its velocity in its parent's frame, and
/// the potential energy it holds by sitting where it does.
///
/// `potential_energy` is what `Aggregate::external_potential` stores — an
/// energy in joules, not a potential — so it is divided by the mass here. A
/// massless or empty node contributes nothing, which is correct: there is
/// nothing there to have a clock.
pub fn physical_rate(velocity: Vec3, potential_energy: f64, mass: f64) -> TimeRate {
    let phi = if mass > 0.0 && potential_energy.is_finite() {
        potential_energy / mass
    } else {
        0.0
    };
    TimeRate {
        kinematic: 1.0 / crate::coords::gamma(velocity),
        gravitational: crate::coords::gravitational_dilation(phi),
        bubble: 1.0,
    }
}

/// Clamp a requested bubble factor to what the arithmetic can carry.
///
/// Returns the accepted value; `None` for a request that is not a positive
/// finite number at all, which is a caller bug rather than an over-ambitious
/// admin and should be reported as one.
pub fn accept_bubble(rate: f64) -> Option<f64> {
    if !rate.is_finite() || rate <= 0.0 {
        return None;
    }
    Some(rate.clamp(MIN_BUBBLE, MAX_BUBBLE))
}

/// The speed at which the kinematic factor becomes visible at a given
/// tolerance — a convenience for tests and for deciding whether it is worth
/// reporting. `1/γ = 1 - tol` at `v/c = sqrt(1-(1-tol)^2)`.
pub fn speed_for_dilation(tolerance: f64) -> f64 {
    let f = (1.0 - tolerance).clamp(0.0, 1.0);
    C * (1.0 - f * f).sqrt()
}
