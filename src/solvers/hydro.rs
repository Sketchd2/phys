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

use crate::eos::Eos;
use crate::math::Vec3;
use crate::neighbourhood::NeighbourGrid;
use crate::solvers::SolveReport;
use crate::state::Body;
use crate::units::*;

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
    step_with(bodies, dt, params, &[])
}

/// The equation of state body `i` answers to: `eos[i]`, or the gas law where
/// the slice does not reach — so an empty slice is every body a gas, which is
/// what [`step`] has always been.
#[inline]
fn eos_at(eos: &[Eos], i: usize) -> Eos {
    eos.get(i).copied().unwrap_or(Eos::Gas)
}

/// Pressure and sound speed of one body at a density.
fn state_of(b: &Body, rho: f64, eos: Eos, gamma: f64) -> (f64, f64) {
    match eos {
        Eos::Condensed(c) => (c.pressure(rho), c.sound_speed(rho)),
        Eos::Gas => {
            let mu = b.composition.mean_molecular_mass(b.temperature);
            let p = if mu > 0.0 { rho * K_B * b.temperature / mu } else { 0.0 }
                + A_RAD * b.temperature.powi(4) / 3.0;
            let cs = if rho > 0.0 { (gamma * p / rho).sqrt() } else { 0.0 };
            (p, cs)
        }
    }
}

/// [`step`], with an equation of state per body.
///
/// `docs/PLAY.md` Phase 5's second piece. A body whose entry is
/// [`Eos::Condensed`] is priced by Tait or Murnaghan rather than by the gas
/// law, and two more things change with it, both because they are gas laws too:
///
/// - **The compression work is carried.** A barotropic pressure does work that
///   the gas branch's temperature-based pressure never booked, and the books
///   close only if it goes somewhere: the standard SPH energy equation,
///   `du_i = (p_i / rho_i^2) sum_j m_j v_ij . grad W_ij`, puts it in the body's
///   internal energy. It is elastic, not heat, so it does not move the
///   temperature.
/// - **No optically thin cooling.** `cooling_rate` is a law for a gas thin
///   enough that every photon escapes. Asked about water it answers 10^20
///   W/m^3, and the half-the-energy cap below turned that into a liquid losing
///   half its internal energy every substep. Condensed matter is optically
///   thick and radiates from its surface, which `evolve_matter` already does.
pub fn step_with(bodies: &mut [Body], dt: f64, params: HydroParams, eos: &[Eos]) -> SolveReport {
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
        let (p, c) = state_of(&bodies[i], rho[i], eos_at(eos, i), params.gamma);
        pressure[i] = p;
        cs[i] = c;
    }

    let grid = NeighbourGrid::build(bodies, 2.0 * params.h);
    let mut acc = vec![Vec3::ZERO; n];
    let mut du = vec![0.0; n];
    // Compression work on condensed bodies, J/kg/s. See `step_with`.
    let mut dw = vec![0.0; n];
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
            if matches!(eos_at(eos, i), Eos::Condensed(_)) {
                dw[i] += pressure[i] / (rho[i] * rho[i]) * bj.mass * grad * v_ij.dot(dir);
            }
            interactions += 1;
        }
    }

    let mut radiated = 0.0;
    for i in 0..n {
        let b = &mut bodies[i];
        b.vel += acc[i].scale(dt);
        b.pos += b.vel.scale(dt);
        let heat = du[i] * dt * b.mass;
        b.internal_energy += heat + dw[i] * dt * b.mass;
        let mu = b.composition.mean_molecular_mass(b.temperature);
        let particles = if mu > 0.0 { b.mass / mu } else { 0.0 };
        if particles > 0.0 {
            b.temperature = (b.temperature + heat / (1.5 * particles * K_B)).max(2.725);
        }
        if params.cooling && rho[i] > 0.0 && !matches!(eos_at(eos, i), Eos::Condensed(_)) {
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
    courant_dt_with(bodies, params, cfl, &[])
}

/// [`courant_dt`], with an equation of state per body — as [`step_with`].
pub fn courant_dt_with(bodies: &[Body], params: HydroParams, cfl: f64, eos: &[Eos]) -> f64 {
    let mut dt = f64::INFINITY;
    for (i, b) in bodies.iter().enumerate() {
        let rho = b.mass / (4.0 / 3.0 * std::f64::consts::PI * params.h.powi(3));
        let c = match eos_at(eos, i) {
            Eos::Condensed(c) => c.sound_speed(rho),
            Eos::Gas => {
                let mu = b.composition.mean_molecular_mass(b.temperature);
                let p = if mu > 0.0 { rho * K_B * b.temperature / mu } else { 0.0 };
                if rho > 0.0 { (params.gamma * p / rho).sqrt() } else { 0.0 }
            }
        };
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
