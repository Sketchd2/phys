//! # phys — a multiscale galaxy engine
//!
//! A physics engine that presents a galaxy consistent down to the subatomic
//! level, in real time, on a consumer GPU.
//!
//! ## The problem, stated honestly
//!
//! A galaxy of 10^9 stars contains about 1.2 x 10^66 baryons. An RTX 2060 can
//! touch roughly 10^7 particles per frame at 20 updates per second. The gap is
//! a factor of 10^59. No optimisation closes a gap like that, and no amount of
//! hardware will: it is 45 orders of magnitude beyond Moore's law's remaining
//! runway.
//!
//! So this engine does not simulate a galaxy. It is *indistinguishable from*
//! one, to any observer inside it, at any resolution they can actually achieve
//! — and it makes that a precise, testable claim rather than a slogan.
//!
//! ## The four ideas
//!
//! 1. **Lazy materialisation with exact conservation** (`prolong`, `state`).
//!    Detail is generated on demand from a seeded distribution and destroyed
//!    when nobody is looking. The generator is constrained so that coarsening
//!    the generated detail returns the original bulk state *exactly* — energy,
//!    momentum, angular momentum, charge, baryon and lepton number. No
//!    experiment performed at the coarse scale can detect the deception.
//!
//! 2. **Causality as a scheduling primitive** (`causal`). Nothing propagates
//!    faster than light, so a region at distance d has a guaranteed lookahead
//!    of d/c. Conservative parallel discrete-event simulation normally struggles
//!    to find a lookahead; relativity hands us one for free, and it is enormous
//!    exactly where the distances are large.
//!
//! 3. **Observation as commitment** (`observe`). An unobserved quantity has no
//!    value; a measured one is recorded permanently in a ledger. At the
//!    subatomic tier this is not an approximation of quantum mechanics, it *is*
//!    quantum mechanics, which is why the trick survives all the way down.
//!
//! 4. **A frame budget that spends detail, not time** (`budget`). Each frame
//!    gets 50 ms and decides what fits. Frame rate is the invariant; fidelity
//!    is the free variable.
//!
//! See `docs/DESIGN.md` for the full argument and `docs/PERFORMANCE.md` for the
//! measured numbers behind the budget.

pub mod budget;
pub mod causal;
pub mod coords;
pub mod engine;
pub mod ids;
pub mod math;
pub mod morph;
pub mod observe;
pub mod prolong;
pub mod rng;
pub mod solvers;
pub mod state;
pub mod tree;
pub mod units;

pub use engine::World;
