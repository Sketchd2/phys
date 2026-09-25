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
    /// A uniform field the contents are held in, m/s^2, in the node's axes.
    ///
    /// Zero for a fluid region that nothing inside its own node holds up —
    /// its weight is carried by the fluid around it, which is outside the
    /// node — and the node's own `gravity` for contents resting on something
    /// the node holds. See `World::advance_node`.
    pub gravity: Vec3,
}

impl Default for HydroParams {
    fn default() -> Self {
        HydroParams {
            h: 1.0,
            gamma: 5.0 / 3.0,
            alpha: 1.0,
            beta: 2.0,
            cooling: true,
            gravity: Vec3::ZERO,
        }
    }
}

/// Something solid the fluid cannot enter: a sphere of `radius`, or a box of
/// `half`-extents turned by `orientation`, in the node's frame. What a node's
/// own ordered members present to its loose contents.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wall {
    pub centre: Vec3,
    pub radius: f64,
    pub half: Vec3,
    pub orientation: crate::math::Quat,
}

impl Wall {
    pub fn of(b: &Body) -> Wall {
        Wall { centre: b.pos, radius: b.radius, half: b.half, orientation: b.orientation }
    }

    /// Signed distance from a point to the surface, and the outward normal.
    pub fn distance(&self, p: Vec3) -> (f64, Vec3) {
        let d = p - self.centre;
        if self.half == Vec3::ZERO {
            let r = d.norm();
            let n = if r > 0.0 { d.scale(1.0 / r) } else { crate::math::v3(0.0, 0.0, 1.0) };
            return (r - self.radius, n);
        }
        let q = self.orientation.conjugate().rotate(d);
        let h = self.half;
        let e = crate::math::v3(q.x.abs() - h.x, q.y.abs() - h.y, q.z.abs() - h.z);
        let out = crate::math::v3(e.x.max(0.0), e.y.max(0.0), e.z.max(0.0));
        let outside = out.norm();
        let (dist, local) = if outside > 0.0 {
            (outside, crate::math::v3(out.x * q.x.signum(), out.y * q.y.signum(), out.z * q.z.signum()).scale(1.0 / outside))
        } else {
            // Inside: the nearest face.
            let m = e.x.max(e.y).max(e.z);
            let n = if m == e.x {
                crate::math::v3(q.x.signum(), 0.0, 0.0)
            } else if m == e.y {
                crate::math::v3(0.0, q.y.signum(), 0.0)
            } else {
                crate::math::v3(0.0, 0.0, q.z.signum())
            };
            (m, n)
        };
        (dist, self.orientation.rotate(local))
    }
}

/// The fraction of the kernel's mass that lies beyond a plane `d` from its
/// centre: `integral_{z > d} W dV`, for the cubic spline, as a function of
/// `d / h` alone.
///
/// **What a wall is to a density sum.** A parcel beside a wall has no
/// neighbours on the wall's side, so the sum over its neighbours comes out
/// short by exactly this much of a full kernel, and Tait's floor turns that
/// deficit into zero pressure. Measured, in a bucket two layers deep: the
/// bottom layer read below rest, pressed with nothing, and the top layer sank
/// into it until all 120 parcels lay in one sheet on the floor. A wall stands
/// where fluid would otherwise be, and counting it as the rest density it
/// displaces is the semi-analytic wall of Kulasegaram and Ferrand — derived
/// from the kernel, which is the only number in it. Half at contact, zero past
/// two smoothing lengths. The free surface is deliberately *not* given one,
/// because a deficit there is what makes it free.
pub fn beyond_plane(d: f64, h: f64) -> f64 {
    static TABLE: std::sync::OnceLock<Vec<f64>> = std::sync::OnceLock::new();
    const N: usize = 256;
    let table = TABLE.get_or_init(|| {
        // A(z) = 2 pi integral_{|z|}^{2} W(r, 1) r dr, then the tail of A.
        let steps = 2000;
        let dz = 2.0 / steps as f64;
        let area: Vec<f64> = (0..=steps)
            .map(|k| {
                let z = k as f64 * dz;
                let m = 400;
                let dr = (2.0 - z) / m as f64;
                let mut a = 0.0;
                for j in 0..m {
                    let r = z + (j as f64 + 0.5) * dr;
                    a += kernel(r, 1.0) * r * dr;
                }
                2.0 * std::f64::consts::PI * a
            })
            .collect();
        // tail[k] = integral_{z_k}^{2} A(z) dz, trapezoid from the far end.
        let mut tail = vec![0.0; steps + 1];
        for k in (0..steps).rev() {
            tail[k] = tail[k + 1] + 0.5 * (area[k] + area[k + 1]) * dz;
        }
        (0..=N).map(|i| tail[(i * steps) / N]).collect()
    });
    if !(h > 0.0) {
        return 0.0;
    }
    let q = (d / h).max(0.0);
    if q >= 2.0 {
        return 0.0;
    }
    let x = q / 2.0 * N as f64;
    let i = (x.floor() as usize).min(N - 1);
    let f = x - i as f64;
    table[i] * (1.0 - f) + table[i + 1] * f
}

/// How far the discrete kernel sum on a cubic lattice of `spacing` departs from
/// the density it is summing: `spacing^3 sum_n W(|n| spacing, h)`.
///
/// Derived, not tuned. A kernel integrates to one and a lattice sum of it does
/// not quite, and the difference is a pure function of `h / spacing`; for a
/// gas it is noise in a density that is free to be anything, and for a liquid
/// at rest it is a pressure, seven times over.
///
/// Measured, for the cubic spline, against `h / spacing`:
///
/// ```text
///   1.0  0.999972     1.5  1.001795     4.0  1.000008
///   1.2  1.000810     2.0  1.000079     6.0  1.000000
///   1.3  0.997262     3.0  0.999979
/// ```
///
/// Past four spacings it is the integral to eight parts in a million, and it
/// is taken as one there rather than summed: a stand-in for a 10 kg ball in a
/// 12 m room has a smoothing length of 110 of its own spacings, and summing
/// that lattice cost 221^3 kernel evaluations per body per substep.
pub fn lattice_sum(h: f64, spacing: f64) -> f64 {
    if !(h > 0.0) || !(spacing > 0.0) || h >= 4.0 * spacing {
        return 1.0;
    }
    let k = (2.0 * h / spacing).ceil() as i64;
    let mut sum = 0.0;
    for ix in -k..=k {
        for iy in -k..=k {
            for iz in -k..=k {
                let r = spacing * ((ix * ix + iy * iy + iz * iz) as f64).sqrt();
                sum += kernel(r, h);
            }
        }
    }
    sum * spacing.powi(3)
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
    step_with(bodies, dt, params, &[], &[])
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
///
/// And three things a liquid needs that a gas does not, which together are
/// weakly-compressible SPH (`docs/PLAY.md` §4.2) — the stiffness itself is the
/// caller's, since it is an equation of state:
///
/// - **A condensed body's density is summed over condensed neighbours**, and
///   divided by [`lattice_sum`] so that a liquid at rest reads its rest
///   density. Counting the air above it would put the free surface at the
///   wrong density, and counting the water below a gas parcel would give the
///   air a thousand times the pressure it has.
/// - **`params.gravity`** acts on every body.
/// - **`walls`** push a body out along their normal once it is within half a
///   parcel spacing (`h / 1.3`) of them — where a cell of liquid resting on a
///   floor has its centre, and where `sampler::packed_positions` lays it — at
///   the stiffness its own sound speed gives: `c^2 (g - d) / g^2` for a gap
///   `g`. Reaching out a whole spacing instead, which was the first version,
///   put every bottom parcel of a fresh draw half a spacing into a spring and
///   launched a lone one 2.8 m into the air.
///
/// Both of the last two do work on the contents from outside them, and that is
/// booked in `non_mechanical_energy` rather than left as drift.
pub fn step_with(
    bodies: &mut [Body],
    dt: f64,
    params: HydroParams,
    eos: &[Eos],
    walls: &[Wall],
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

    let condensed: Vec<bool> = (0..n).map(|i| matches!(eos_at(eos, i), Eos::Condensed(_))).collect();
    let mut rho = if condensed.iter().any(|c| *c) {
        let grid = NeighbourGrid::build(bodies, 2.0 * params.h);
        let mut out = vec![0.0; n];
        let mut nb = Vec::with_capacity(128);
        for (i, b) in bodies.iter().enumerate() {
            grid.neighbours(b.pos, &mut nb);
            out[i] = nb
                .iter()
                .map(|&j| j as usize)
                .filter(|&j| condensed[j] == condensed[i])
                .map(|j| bodies[j].mass * kernel((b.pos - bodies[j].pos).norm(), params.h))
                .sum();
        }
        out
    } else {
        densities(bodies, params)
    };
    // The kernel's own bias on a lattice at the spacing the rest density
    // gives, per condensed body, and then the walls' share of its sum: see
    // `beyond_plane`.
    for i in 0..n {
        if let Eos::Condensed(c) = eos_at(eos, i) {
            let spacing = (bodies[i].mass / c.rest_density).cbrt();
            rho[i] /= lattice_sum(params.h, spacing);
            for w in walls {
                if (bodies[i].pos - w.centre).norm() > w.radius.max(w.half.norm()) + 2.0 * params.h {
                    continue;
                }
                let (d, _) = w.distance(bodies[i].pos);
                rho[i] += c.rest_density * beyond_plane(d, params.h);
            }
        }
    }
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

    // The field, and the walls, as accelerations from outside.
    let spacing = 0.5 * params.h / 1.3;
    for i in 0..n {
        acc[i] += params.gravity;
        if walls.is_empty() {
            continue;
        }
        let p = bodies[i].pos;
        for w in walls {
            // Cheap reject before the exact distance.
            let reach = w.radius.max(w.half.norm()) + spacing;
            if (p - w.centre).norm() > reach {
                continue;
            }
            let (d, normal) = w.distance(p);
            if d < spacing {
                let c = cs[i].max(1e-30);
                acc[i] += normal.scale(c * c * (spacing - d) / (spacing * spacing));
                // And the wall's share of the artificial viscosity. A spring
                // with nothing to damp it rings for ever: measured, a bucket's
                // bottom layer bounced at +-0.5 m/s for the whole of a second
                // while the parcels above it, which the pairwise viscosity
                // does reach, were still. A dashpot at the spring's own
                // frequency `c / g` — damping ratio a half — against the
                // normal velocity, and what it takes out of the motion goes
                // into the parcel as heat, exactly as the pairwise term's does.
                let vn = bodies[i].vel.dot(normal);
                let damp = -(c / spacing) * vn;
                acc[i] += normal.scale(damp);
                du[i] += -damp * vn;
            }
        }
    }

    let mut radiated = 0.0;
    let mut external = 0.0;
    for i in 0..n {
        let b = &mut bodies[i];
        // Work done from outside: the field and the walls, less the pairwise
        // part, which is internal and conserves. Measured on the velocity the
        // step actually moves the body with.
        // The wall's dashpot is not in this: its work is booked as heat in
        // `du`, which the conservation check already sees.
        let outside = params.gravity + {
            let mut a = Vec3::ZERO;
            for w in walls {
                let (d, normal) = w.distance(b.pos);
                if d < spacing {
                    let c = cs[i].max(1e-30);
                    a += normal.scale(c * c * (spacing - d) / (spacing * spacing));
                }
            }
            a
        };
        b.vel += acc[i].scale(dt);
        external += b.mass * outside.dot(b.vel) * dt;
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
        non_mechanical_energy: external - radiated,
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
    // A body in a field accelerates across its own smoothing length in
    // `sqrt(h / g)`, and a step longer than that lets it fall through
    // whatever holds it up.
    let g = params.gravity.norm();
    let mut dt = if g > 0.0 { cfl * (params.h / g).sqrt() } else { f64::INFINITY };
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
