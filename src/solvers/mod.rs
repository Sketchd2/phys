//! Per-tier solvers.
//!
//! Each tier of the ladder gets the solver that is *right* for it rather than
//! one universal method run at different resolutions. A single method cannot
//! span this range: symplectic N-body integration is correct for a galaxy and
//! meaningless for a nucleus; the Schrödinger equation is correct for an atom
//! and computationally absurd for a molecular cloud.
//!
//! What holds them together is the interface, not the method: every solver
//! takes a set of bodies plus a timestep, and every solver is required to
//! report the conserved tuple it started with and ended with. A solver that
//! cannot state its energy budget cannot be part of the ladder, because the
//! scale-transition guarantee in `prolong.rs` would have nothing to stand on.

pub mod gravity;
pub mod hydro;
pub mod md;
pub mod nuclear;
pub mod quantum;

use crate::state::{Body, Conserved};
use crate::units::Tier;

/// What a solver did, and what it cost. The scheduler uses `cost_estimate` to
/// budget the next frame; the tests use `before`/`after` to assert conservation.
#[derive(Debug, Clone, Copy, Default)]
pub struct SolveReport {
    pub steps: u32,
    pub interactions: u64,
    pub dt_used: f64,
    pub before: Conserved,
    pub after: Conserved,
    /// Energy deliberately injected or removed by physics that is not
    /// mechanical — radiative cooling, fusion, decay. Subtracted before the
    /// conservation check, and reported so nothing hides in it.
    pub non_mechanical_energy: f64,
}

impl SolveReport {
    /// Relative energy drift after accounting for declared sources and sinks.
    pub fn drift(&self) -> f64 {
        let expected = self.before.energy + self.non_mechanical_energy;
        let scale = expected.abs().max(self.after.energy.abs());
        if scale > 0.0 {
            (self.after.energy - expected).abs() / scale
        } else {
            0.0
        }
    }

    pub fn momentum_drift(&self) -> f64 {
        let d = (self.after.momentum - self.before.momentum).norm();
        let s = self
            .before
            .momentum
            .norm()
            .max(self.after.momentum.norm())
            .max(self.before.energy.abs() / crate::units::C);
        if s > 0.0 {
            d / s
        } else {
            0.0
        }
    }
}

/// Pick the solver for a tier. The `Continuum` tier is deliberately served by
/// the hydro solver rather than a separate FEM path: at the resolutions the
/// budget allows, a meshless particle method and a mesh method are
/// indistinguishable to an observer, and the particle method composes with the
/// tiers on either side of it without a remeshing step.
pub fn for_tier(tier: Tier) -> SolverKind {
    match tier {
        Tier::Galactic | Tier::Stellar => SolverKind::Gravity,
        Tier::Planetary => SolverKind::GravityHydro,
        Tier::Continuum => SolverKind::Hydro,
        Tier::Molecular | Tier::Atomic => SolverKind::MolecularDynamics,
        Tier::Nuclear => SolverKind::Statistical,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverKind {
    Gravity,
    GravityHydro,
    Hydro,
    MolecularDynamics,
    Statistical,
}

impl SolverKind {
    /// Estimated cost in "interaction units" for `n` bodies over one step.
    /// Feeds the frame budget in `budget.rs`; the constants come from the
    /// measurements in `docs/PERFORMANCE.md`.
    pub fn cost(self, n: usize) -> f64 {
        let n = n as f64;
        match self {
            SolverKind::Gravity => n * n.max(2.0).log2() * 1.4,
            SolverKind::GravityHydro => n * n.max(2.0).log2() * 2.2,
            SolverKind::Hydro => n * 48.0,
            SolverKind::MolecularDynamics => n * 96.0,
            SolverKind::Statistical => n * 4.0,
        }
    }
}

/// Total conserved quantities of a materialised set, about the origin.
pub fn measure(bodies: &[Body], potential: f64) -> Conserved {
    let mut c = crate::state::restrict(bodies, potential).conserved();
    // `restrict` measures spin about the centre of mass; for a solver check we
    // want angular momentum about a fixed origin, which is what actually has to
    // be conserved under internal forces.
    c.angular_momentum = crate::state::total_spin(bodies, crate::math::Vec3::ZERO);
    c
}
