//! Smoothed-particle hydrodynamics for the gas and continuum tiers.
//!
//! SPH rather than a grid, for one architectural reason: it is *meshless*, so
//! it composes with the tiers on either side of it without a remeshing step. A
//! gas parcel that condenses into a protostar is promoted to a `Planetary` node
//! by changing its interpretation, not by interpolating it onto a new grid —
//! and interpolation between grids is exactly where conservation goes to die.
//!
//! The formulation is the standard density-entropy one with Monaghan
//! artificial viscosity, which conserves momentum and angular momentum exactly
//! (forces are pairwise and antisymmetric) and energy to integrator accuracy.

use crate::math::Vec3;
use crate::solvers::SolveReport;
use crate::state::Body;
use crate::units::*;
use std::collections::HashMap;

/// Cubic spline kernel (Monaghan & Lattanzio 1985), normalised in 3D.
#[inline]
pub fn kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        sigma * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}

/// Radial derivative of the kernel, `dW/dr` (negative inside the support).
#[inline]
pub fn kernel_grad(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        sigma * -0.75 * (2.0 - q).powi(2) / h
    } else {
        0.0
    }
}

#[derive(Debug, Clone, Copy)]
pub struct HydroParams {
    /// Smoothing length. Set from the mean interparticle spacing so that each
    /// particle has ~50 neighbours, the standard compromise between noise and
    /// resolution.
    pub h: f64,
    pub gamma: f64,
    /// Monaghan viscosity coefficients. Without these, shocks are not captured
    /// and supernova blast waves simply pass through each other.
    pub alpha: f64,
    pub beta: f64,
    /// Enable optically-thin radiative cooling. This is what lets gas collapse:
    /// without a cooling channel, compression heats gas until pressure stops
    /// it, and no star ever forms.
    pub cooling: bool,
}

impl Default for HydroParams {
    fn default() -> Self {
        HydroParams {
            h: 1.0,
            gamma: 5.0 / 3.0,
            alpha: 1.0,
            beta: 2.0,
            cooling: true,
        }
    }
}

/// Uniform spatial hash for neighbour finding — O(n) build, O(1) query.
/// The cell size is `2h`, the kernel support radius, so a query touches 27
/// cells and nothing outside them can be a neighbour.
pub struct NeighbourGrid {
    cells: HashMap<(i64, i64, i64), Vec<u32>>,
    spacing: f64,
}

impl NeighbourGrid {
    pub fn build(bodies: &[Body], spacing: f64) -> NeighbourGrid {
        let mut cells: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
        let s = spacing.max(1e-30);
        for (i, b) in bodies.iter().enumerate() {
            cells.entry(key_of(b.pos, s)).or_default().push(i as u32);
        }
        NeighbourGrid { cells, spacing: s }
    }

    pub fn neighbours(&self, pos: Vec3, out: &mut Vec<u32>) {
        out.clear();
        let (cx, cy, cz) = key_of(pos, self.spacing);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(v) = self.cells.get(&(cx + dx, cy + dy, cz + dz)) {
                        out.extend_from_slice(v);
                    }
                }
            }
        }
        // Deterministic order regardless of hash iteration order — otherwise
        // the floating-point sum below depends on the hasher's seed and replay
        // diverges between runs.
        out.sort_unstable();
    }

    pub fn occupied_cells(&self) -> usize {
        self.cells.len()
    }
}

#[inline]
fn key_of(p: Vec3, s: f64) -> (i64, i64, i64) {
    (
        (p.x / s).floor() as i64,
        (p.y / s).floor() as i64,
        (p.z / s).floor() as i64,
    )
}

/// Densities by kernel summation.
pub fn densities(bodies: &[Body], params: HydroParams) -> Vec<f64> {
    let grid = NeighbourGrid::build(bodies, 2.0 * params.h);
    let mut out = vec![0.0; bodies.len()];
    let mut nb = Vec::with_capacity(128);
    for (i, b) in bodies.iter().enumerate() {
        grid.neighbours(b.pos, &mut nb);
        let mut rho = 0.0;
        for &j in nb.iter() {
            let o = &bodies[j as usize];
            let r = (b.pos - o.pos).norm();
            rho += o.mass * kernel(r, params.h);
        }
        out[i] = rho;
    }
    out
}

/// One SPH step: densities, pressures, forces, then a velocity-Verlet update.
///
/// Returns the energy explicitly radiated away, so the caller can subtract it
/// before checking conservation instead of quietly tolerating a drift.
pub fn step(bodies: &mut [Body], dt: f64, params: HydroParams) -> SolveReport {
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

    let rho = densities(bodies, params);
    let mut pressure = vec![0.0; n];
    let mut cs = vec![0.0; n];
    for i in 0..n {
        let b = &bodies[i];
        let mu = b.composition.mean_molecular_mass(b.temperature);
        let p = if mu > 0.0 {
            rho[i] * K_B * b.temperature / mu
        } else {
            0.0
        } + A_RAD * b.temperature.powi(4) / 3.0;
        pressure[i] = p;
        cs[i] = if rho[i] > 0.0 {
            (params.gamma * p / rho[i]).sqrt()
        } else {
            0.0
        };
    }

    let grid = NeighbourGrid::build(bodies, 2.0 * params.h);
    let mut acc = vec![Vec3::ZERO; n];
    let mut du = vec![0.0; n];
    let mut nb = Vec::with_capacity(128);
    let mut interactions = 0u64;

    for i in 0..n {
        grid.neighbours(bodies[i].pos, &mut nb);
        let bi = bodies[i];
        if rho[i] <= 0.0 {
            continue;
        }
        for &jj in nb.iter() {
            let j = jj as usize;
            if j == i || rho[j] <= 0.0 {
                continue;
            }
            let bj = bodies[j];
            let d = bi.pos - bj.pos;
            let r = d.norm();
            if r <= 0.0 || r >= 2.0 * params.h {
                continue;
            }
            let grad = kernel_grad(r, params.h);
            let dir = d.scale(1.0 / r);

            // Symmetric pressure form: the force on i from j is exactly minus
            // the force on j from i, so momentum and angular momentum are
            // conserved to machine precision rather than to truncation.
            let term = pressure[i] / (rho[i] * rho[i]) + pressure[j] / (rho[j] * rho[j]);

            // Monaghan artificial viscosity, active only in compression.
            let v_ij = bi.vel - bj.vel;
            let vr = v_ij.dot(d);
            let visc = if vr < 0.0 {
                let h = params.h;
                let mu_ij = h * vr / (r * r + 0.01 * h * h);
                let c_bar = 0.5 * (cs[i] + cs[j]);
                let rho_bar = 0.5 * (rho[i] + rho[j]);
                (-params.alpha * c_bar * mu_ij + params.beta * mu_ij * mu_ij) / rho_bar
            } else {
                0.0
            };

            let f = bj.mass * (term + visc) * grad;
            acc[i] += dir.scale(-f);
            // Viscous heating: the energy the viscosity removes from bulk
            // motion reappears as heat. Dropping this term is the classic way
            // to lose 10% of a shock's energy.
            du[i] += 0.5 * bj.mass * visc * grad * v_ij.dot(dir);
            interactions += 1;
        }
    }

    let mut radiated = 0.0;
    for i in 0..n {
        let b = &mut bodies[i];
        b.vel += acc[i].scale(dt);
        b.pos += b.vel.scale(dt);
        let heat = du[i] * dt * b.mass;
        b.internal_energy += heat;
        let mu = b.composition.mean_molecular_mass(b.temperature);
        let particles = if mu > 0.0 { b.mass / mu } else { 0.0 };
        if particles > 0.0 {
            b.temperature = (b.temperature + heat / (1.5 * particles * K_B)).max(2.725);
        }
        if params.cooling && rho[i] > 0.0 {
            let loss = cooling_rate(b.temperature, rho[i], b.composition.metallicity())
                * (b.mass / rho[i])
                * dt;
            let capped = loss.min(b.internal_energy.max(0.0) * 0.5);
            b.internal_energy -= capped;
            radiated += capped;
            if particles > 0.0 {
                b.temperature = (b.temperature - capped / (1.5 * particles * K_B)).max(2.725);
            }
        }
    }

    let after = crate::solvers::measure(bodies, 0.0);
    SolveReport {
        steps: 1,
        interactions,
        dt_used: dt,
        before,
        after,
        non_mechanical_energy: -radiated,
    }
}

/// Optically-thin cooling rate in W/m^3.
///
/// A three-regime fit to the standard collisional-ionisation-equilibrium curve:
/// molecular/atomic line cooling below 10^4 K, the Lyman-alpha peak around
/// 10^4-10^5 K, and bremsstrahlung above 10^7 K. Metallicity scales the line
/// cooling, which is why the first generation of stars formed differently from
/// later ones — a difference this engine reproduces for free.
pub fn cooling_rate(t: f64, rho: f64, metallicity: f64) -> f64 {
    if t <= 10.0 {
        return 0.0;
    }
    let n = rho / M_PROTON; // number density, m^-3
    let lambda = if t < 1e4 {
        // Molecular/fine-structure cooling, strongly metallicity-dependent.
        1e-40 * (t / 100.0).powf(2.0) * (0.01 + metallicity * 30.0)
    } else if t < 1e7 {
        // Line cooling: peak near 10^5 K.
        let x = (t / 1e5).ln();
        1e-35 * (-x * x * 0.5).exp() * (0.1 + metallicity * 30.0)
    } else {
        // Free-free.
        2.3e-40 * t.sqrt()
    };
    lambda * n * n
}

/// Courant condition, including the viscous signal speed.
pub fn courant_dt(bodies: &[Body], params: HydroParams, cfl: f64) -> f64 {
    let mut dt = f64::INFINITY;
    for b in bodies {
        let mu = b.composition.mean_molecular_mass(b.temperature);
        let rho = b.mass / (4.0 / 3.0 * std::f64::consts::PI * params.h.powi(3));
        let p = if mu > 0.0 { rho * K_B * b.temperature / mu } else { 0.0 };
        let c = if rho > 0.0 { (params.gamma * p / rho).sqrt() } else { 0.0 };
        let v = b.vel.norm();
        let signal = c + v + 1e-30;
        dt = dt.min(cfl * params.h / signal);
    }
    dt
}

/// Jeans criterion: does this parcel have to be refined, or may it stay coarse?
///
/// A node whose Jeans length is unresolved is on the verge of collapsing into
/// structure the engine cannot see, and refusing to refine it is how a
/// simulation quietly produces a galaxy with no stars in it.
pub fn needs_refinement(rho: f64, temperature: f64, mu: f64, h: f64) -> bool {
    if rho <= 0.0 {
        return false;
    }
    let cs = (1.6667 * K_B * temperature / mu).sqrt();
    let jeans = cs * (std::f64::consts::PI / (G * rho)).sqrt();
    jeans < 4.0 * h
}
