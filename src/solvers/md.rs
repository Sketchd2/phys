//! Molecular dynamics for the molecular and atomic tiers.
//!
//! Velocity Verlet with a Lennard-Jones + shifted-force Coulomb potential and a
//! cell-list neighbour search. Short-ranged, so cost is linear in particle
//! count — which is why the molecular tier can afford more particles per frame
//! than the galactic tier can, despite each one being far smaller.
//!
//! The thermostat is deliberately a *deterministic* Langevin variant: the
//! random force is drawn from the node's own seeded stream rather than from
//! system entropy, so a molecular simulation replays exactly. A thermostat that
//! breaks replay would break the engine's core guarantee for the sake of a
//! convenience.

use crate::math::Vec3;
use crate::rng::{Purpose, Stream};
use crate::solvers::hydro::NeighbourGrid;
use crate::solvers::SolveReport;
use crate::state::Body;
use crate::units::*;

/// Lennard-Jones parameters per species: `(sigma [m], epsilon [J])`.
/// Values are the standard UFF-like set, adequate for the qualitative
/// behaviour an observer can actually see at this tier.
pub fn lj_params(s: Species) -> (f64, f64) {
    match s {
        Species::Hydrogen => (2.571e-10, 0.0184 * EV),
        Species::Helium => (2.640e-10, 0.0024 * EV),
        Species::Carbon => (3.431e-10, 0.0046 * EV),
        Species::Nitrogen => (3.261e-10, 0.0031 * EV),
        Species::Oxygen => (3.118e-10, 0.0026 * EV),
        Species::Silicon => (3.826e-10, 0.0175 * EV),
        Species::Iron => (2.912e-10, 0.0056 * EV),
        Species::Other => (3.500e-10, 0.0050 * EV),
    }
}

/// Dominant species of a body, for force-field lookup.
fn dominant(b: &Body) -> Species {
    let mut best = Species::Hydrogen;
    let mut bv = -1.0;
    for s in Species::ALL {
        let v = b.composition.get(s);
        if v > bv {
            bv = v;
            best = s;
        }
    }
    best
}

#[derive(Debug, Clone, Copy)]
pub struct MdParams {
    /// Interaction cutoff. 2.5 sigma is the convention; beyond it the LJ tail
    /// contributes less than 1% of the energy.
    pub cutoff: f64,
    /// Coulomb screening length. Real MD uses Ewald summation for the
    /// long-range part; here the engine leans on the tier above instead — a
    /// charge imbalance large enough to matter at long range is, by
    /// construction, visible in the parent node's aggregate charge.
    pub debye: f64,
    pub thermostat: Option<f64>,
    /// Langevin friction, 1/s.
    pub friction: f64,
}

impl Default for MdParams {
    fn default() -> Self {
        MdParams {
            cutoff: 1.0e-9,
            debye: 1.0e-9,
            thermostat: None,
            friction: 1e12,
        }
    }
}

/// Lennard-Jones force magnitude along the separation (positive = repulsive).
#[inline]
pub fn lj_force(r: f64, sigma: f64, epsilon: f64) -> f64 {
    if r <= 0.0 {
        return 0.0;
    }
    let sr6 = (sigma / r).powi(6);
    24.0 * epsilon * (2.0 * sr6 * sr6 - sr6) / r
}

#[inline]
pub fn lj_potential(r: f64, sigma: f64, epsilon: f64) -> f64 {
    if r <= 0.0 {
        return 0.0;
    }
    let sr6 = (sigma / r).powi(6);
    4.0 * epsilon * (sr6 * sr6 - sr6)
}

/// One velocity-Verlet step.
/// One velocity-Verlet step.
///
/// `tick` is the node's step counter. It must advance, or the thermostat's
/// noise repeats — see `Stream::split`.
pub fn step(
    bodies: &mut [Body],
    dt: f64,
    params: MdParams,
    world_seed: u64,
    path_key: u128,
    epoch: u32,
    tick: u64,
) -> SolveReport {
    let before = crate::solvers::measure(bodies, 0.0);
    let n = bodies.len();
    if n == 0 || dt == 0.0 {
        return SolveReport {
            before,
            after: before,
            dt_used: dt,
            ..Default::default()
        };
    }

    let acc = forces(bodies, params);
    for (b, a) in bodies.iter_mut().zip(&acc) {
        b.vel += a.scale(0.5 * dt);
        b.pos += b.vel.scale(dt);
    }
    let acc2 = forces(bodies, params);
    let mut interactions = 0u64;
    for (b, a) in bodies.iter_mut().zip(&acc2) {
        b.vel += a.scale(0.5 * dt);
        interactions += 1;
    }

    let mut thermostat_energy = 0.0;
    if let Some(target) = params.thermostat {
        let base = Stream::at(world_seed, path_key, epoch, Purpose::ThermalNoise);
        let before_ke = crate::state::kinetic_energy_of(bodies);
        for (i, b) in bodies.iter_mut().enumerate() {
            if b.mass <= 0.0 {
                continue;
            }
            // Address includes both the step and the body, so no two draws in
            // the whole simulation share a stream position.
            let mut stream = base.split(tick.wrapping_mul(0x9E37_79B9).wrapping_add(i as u64));
            // Langevin: friction plus a fluctuation of exactly the size the
            // fluctuation-dissipation theorem requires. Getting the ratio
            // wrong is the standard way to build a thermostat that silently
            // pumps energy.
            let gamma = params.friction;
            let sigma = (2.0 * gamma * K_B * target / b.mass).sqrt() * dt.sqrt();
            let noise = stream.normal3().scale(sigma);
            b.vel = b.vel.scale((1.0 - gamma * dt).max(0.0)) + noise;
        }
        thermostat_energy = crate::state::kinetic_energy_of(bodies) - before_ke;
    }

    let after = crate::solvers::measure(bodies, 0.0);
    SolveReport {
        steps: 1,
        interactions,
        dt_used: dt,
        before,
        after,
        non_mechanical_energy: thermostat_energy,
    }
}

/// Pairwise forces via cell lists.
pub fn forces(bodies: &[Body], params: MdParams) -> Vec<Vec3> {
    let n = bodies.len();
    let mut acc = vec![Vec3::ZERO; n];
    if n == 0 {
        return acc;
    }
    let grid = NeighbourGrid::build(bodies, params.cutoff);
    let mut nb = Vec::with_capacity(128);
    for i in 0..n {
        grid.neighbours(bodies[i].pos, &mut nb);
        let bi = bodies[i];
        let (si, ei) = lj_params(dominant(&bi));
        let mut a = Vec3::ZERO;
        for &jj in nb.iter() {
            let j = jj as usize;
            if j == i {
                continue;
            }
            let bj = bodies[j];
            let d = bi.pos - bj.pos;
            let r = d.norm();
            if r <= 0.0 || r > params.cutoff {
                continue;
            }
            let (sj, ej) = lj_params(dominant(&bj));
            // Lorentz-Berthelot mixing.
            let sigma = 0.5 * (si + sj);
            let epsilon = (ei * ej).sqrt();
            let mut f = lj_force(r, sigma, epsilon);
            if bi.charge != 0.0 && bj.charge != 0.0 {
                // Screened Coulomb, shifted so the force goes smoothly to zero
                // at the cutoff — an unshifted cutoff puts an impulse into
                // every pair that crosses it and heats the system.
                let scr = (-r / params.debye).exp();
                let fc = K_COULOMB * bi.charge * bj.charge / (r * r) * scr;
                let fc_cut = K_COULOMB * bi.charge * bj.charge
                    / (params.cutoff * params.cutoff)
                    * (-params.cutoff / params.debye).exp();
                f += fc - fc_cut;
            }
            if bi.mass > 0.0 {
                a += d.scale(f / (r * bi.mass));
            }
        }
        acc[i] = a;
    }
    acc
}

/// Potential energy of the configuration.
pub fn potential_energy(bodies: &[Body], params: MdParams) -> f64 {
    let n = bodies.len();
    let grid = NeighbourGrid::build(bodies, params.cutoff);
    let mut nb = Vec::with_capacity(128);
    let mut total = 0.0;
    for i in 0..n {
        grid.neighbours(bodies[i].pos, &mut nb);
        let bi = bodies[i];
        let (si, ei) = lj_params(dominant(&bi));
        for &jj in nb.iter() {
            let j = jj as usize;
            if j <= i {
                continue;
            }
            let bj = bodies[j];
            let r = (bi.pos - bj.pos).norm();
            if r <= 0.0 || r > params.cutoff {
                continue;
            }
            let (sj, ej) = lj_params(dominant(&bj));
            total += lj_potential(r, 0.5 * (si + sj), (ei * ej).sqrt());
            if bi.charge != 0.0 && bj.charge != 0.0 {
                total += K_COULOMB * bi.charge * bj.charge / r * (-r / params.debye).exp();
            }
        }
    }
    total
}

/// Stable timestep: a fraction of the fastest vibrational period. Getting this
/// wrong is immediately obvious — the system explodes within a few hundred
/// steps.
pub fn stable_dt(bodies: &[Body]) -> f64 {
    let mut dt: f64 = 1e-14;
    for b in bodies {
        let (sigma, epsilon) = lj_params(dominant(b));
        if b.mass > 0.0 && epsilon > 0.0 {
            let omega = (72.0 * epsilon / (b.mass * sigma * sigma)).sqrt();
            dt = dt.min(0.02 * std::f64::consts::TAU / omega);
        }
    }
    dt
}

/// Instantaneous kinetic temperature of a materialised set.
pub fn temperature_of(bodies: &[Body]) -> f64 {
    if bodies.is_empty() {
        return 0.0;
    }
    let ke = crate::state::kinetic_energy_of(bodies);
    let dof = 3.0 * bodies.len() as f64;
    2.0 * ke / (dof * K_B)
}
