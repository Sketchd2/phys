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

/// Pairs excluded from the nonbonded sum.
///
/// Two atoms held at a covalent bond length are far inside each other's
/// Lennard-Jones radius — a hydrogen pair sits at 0.74 angstroms with a sigma
/// of 2.57 — so the repulsive term between them is enormous and entirely
/// spurious. The bond already describes that interaction; leaving the van der
/// Waals term in as well double-counts it and then some, and a molecule
/// assembled with both blows itself apart in a few hundred femtoseconds.
///
/// Every force field excludes bonded neighbours for this reason. Both the
/// directly bonded pairs and the pairs across an angle are excluded, which is
/// the standard 1-2 and 1-3 treatment.
#[derive(Debug, Clone, Default)]
pub struct Exclusions {
    /// Sorted neighbour list per atom.
    by_atom: Vec<Vec<u32>>,
}

impl Exclusions {
    pub fn none() -> Exclusions {
        Exclusions::default()
    }

    #[inline]
    pub fn excluded(&self, i: usize, j: u32) -> bool {
        match self.by_atom.get(i) {
            Some(list) => list.binary_search(&j).is_ok(),
            None => false,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.by_atom.is_empty()
    }
}

impl Bonded {
    /// The 1-2 and 1-3 pairs this bond set implies.
    pub fn exclusions(&self, n: usize) -> Exclusions {
        let mut by_atom = vec![Vec::new(); n];
        let mut add = |a: u32, b: u32| {
            if (a as usize) < n && (b as usize) < n && a != b {
                by_atom[a as usize].push(b);
                by_atom[b as usize].push(a);
            }
        };
        for b in &self.bonds {
            add(b.a, b.b);
        }
        for a in &self.angles {
            add(a.a, a.b);
            add(a.b, a.c);
            add(a.a, a.c);
        }
        for list in by_atom.iter_mut() {
            list.sort_unstable();
            list.dedup();
        }
        Exclusions { by_atom }
    }
}

/// Pairwise forces via cell lists.
pub fn forces(bodies: &[Body], params: MdParams) -> Vec<Vec3> {
    forces_excluding(bodies, params, &Exclusions::none())
}

/// As [`forces`], skipping pairs the bond set already describes.
pub fn forces_excluding(bodies: &[Body], params: MdParams, skip: &Exclusions) -> Vec<Vec3> {
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
            if j == i || skip.excluded(i, jj) {
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
    potential_energy_excluding(bodies, params, &Exclusions::none())
}

/// As [`potential_energy`], skipping pairs the bond set already describes.
pub fn potential_energy_excluding(bodies: &[Body], params: MdParams, skip: &Exclusions) -> f64 {
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
            if j <= i || skip.excluded(i, jj) {
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

/// Stable timestep for a bonded system.
///
/// A covalent bond is two orders of magnitude stiffer than the van der Waals
/// interaction beside it, so the unbonded timestep integrates the bonds
/// unstably — visibly, within a few hundred steps. A fiftieth of the shortest
/// vibrational period keeps Verlet's energy error under a part in ten thousand.
pub fn stable_dt_bonded(bodies: &[Body], bonded: &Bonded) -> f64 {
    let period = bonded.shortest_period(bodies);
    if period.is_finite() && period > 0.0 {
        stable_dt(bodies).min(period / 50.0)
    } else {
        stable_dt(bodies)
    }
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

// ---------------------------------------------------------------------------
// Bonded interactions
// ---------------------------------------------------------------------------

/// Covalent bonds and the angles between them.
///
/// # Why the particle tier needed this
///
/// Lennard-Jones plus screened Coulomb is a fine description of a gas and a
/// reasonable one of a liquid. It is not a description of a *molecule*. Nothing
/// in it distinguishes two hydrogens that are bonded from two that happen to be
/// near each other, so a water molecule handed to the molecular tier came apart
/// the moment anything warmed it: the bond that held it was a van der Waals
/// well two hundred times too shallow.
///
/// This is the particle-tier half of the same gap the structural tier had. A
/// bond decided what broke and never what moved.
///
/// # Morse, not harmonic
///
/// ```text
///     V(r) = D_e [1 - e^{-a(r - r0)}]^2,     a = sqrt(k / 2 D_e)
/// ```
///
/// A harmonic bond has the right stiffness near equilibrium and is infinitely
/// strong: it can be stretched across the simulation box and will still pull
/// back. Dissociation then has to be bolted on as a rule — break at some
/// extension, or above some energy — and the rule is a number nobody can
/// defend.
///
/// The Morse potential has the same curvature at the bottom of the well, so it
/// reproduces the vibrational frequency exactly, and it flattens out at `D_e`.
/// A molecule given more than its dissociation energy comes apart because the
/// potential runs out, not because a threshold fired. That is the same
/// principle as the rest of the engine: the failure is in the physics, not in a
/// branch above it.
#[derive(Debug, Clone, Default)]
pub struct Bonded {
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
}

/// A Morse bond between two particles.
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    pub a: u32,
    pub b: u32,
    /// Equilibrium separation, m.
    pub r0: f64,
    /// Well depth, J. The energy needed to pull the bond apart from rest.
    pub well: f64,
    /// Range parameter, 1/m. `sqrt(k / 2 D_e)` for a force constant `k`.
    pub alpha: f64,
}

impl Bond {
    /// A bond with a given harmonic force constant `k` (N/m) and well depth.
    pub fn new(a: u32, b: u32, r0: f64, well: f64, k: f64) -> Bond {
        let alpha = if well > 0.0 { (k / (2.0 * well)).sqrt() } else { 0.0 };
        Bond { a, b, r0, well, alpha }
    }

    /// Harmonic force constant at the bottom of the well, N/m.
    #[inline]
    pub fn force_constant(&self) -> f64 {
        2.0 * self.well * self.alpha * self.alpha
    }

    /// Potential energy at separation `r`, measured from the well bottom.
    #[inline]
    pub fn energy(&self, r: f64) -> f64 {
        let x = 1.0 - (-self.alpha * (r - self.r0)).exp();
        self.well * x * x
    }

    /// Attractive force magnitude at `r`. Positive pulls the pair together.
    #[inline]
    pub fn tension(&self, r: f64) -> f64 {
        let e = (-self.alpha * (r - self.r0)).exp();
        2.0 * self.alpha * self.well * (1.0 - e) * e
    }

    /// Separation at which the restoring force peaks, `r0 + ln2/alpha`. Past
    /// this the bond is on its way apart: pulling harder makes it weaker, which
    /// is what dissociation *is*.
    #[inline]
    pub fn inflection(&self) -> f64 {
        if self.alpha > 0.0 {
            self.r0 + std::f64::consts::LN_2 / self.alpha
        } else {
            f64::INFINITY
        }
    }
}

/// A harmonic bend at particle `b`, between bonds `b-a` and `b-c`.
///
/// Angles are harmonic rather than Morse because they do not dissociate: a
/// molecule loses its shape by breaking a bond, not by opening an angle to
/// infinity.
#[derive(Debug, Clone, Copy)]
pub struct Angle {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    /// Rest angle at `b`, radians.
    pub rest: f64,
    /// Bending constant, J/rad^2.
    pub stiffness: f64,
}

/// Spectroscopic constants for a covalent bond: `(r0 [m], D_e [J], k [N/m])`.
///
/// Measured values for the pairs the engine actually materialises, and a
/// generic single bond for the rest. The point of using real numbers is that
/// the vibrational frequency, the dissociation energy and the bond length are
/// then not three independent knobs — fixing any two fixes the third, and the
/// tests check that the solver reproduces all of them.
pub fn covalent(a: Species, b: Species) -> (f64, f64, f64) {
    use Species::*;
    let (lo, hi) = if (a as u8) <= (b as u8) { (a, b) } else { (b, a) };
    let (r0_ang, de_ev, k) = match (lo, hi) {
        (Hydrogen, Hydrogen) => (0.741, 4.75, 575.0),
        (Hydrogen, Carbon) => (1.090, 4.28, 490.0),
        (Hydrogen, Nitrogen) => (1.010, 4.05, 630.0),
        (Hydrogen, Oxygen) => (0.958, 4.81, 845.0),
        (Carbon, Carbon) => (1.540, 3.60, 450.0),
        (Carbon, Nitrogen) => (1.470, 3.17, 500.0),
        (Carbon, Oxygen) => (1.430, 3.70, 500.0),
        (Nitrogen, Nitrogen) => (1.098, 9.79, 2295.0),
        (Oxygen, Oxygen) => (1.208, 5.16, 1140.0),
        (Silicon, Oxygen) => (1.630, 8.30, 600.0),
        (Silicon, Silicon) => (2.330, 3.21, 200.0),
        (Iron, Iron) => (2.480, 1.15, 140.0),
        _ => (1.500, 3.50, 400.0),
    };
    (r0_ang * 1.0e-10, de_ev * EV, k)
}

impl Bonded {
    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty() && self.angles.is_empty()
    }

    /// Bond two particles using the constants for their dominant species.
    pub fn bond(&mut self, bodies: &[Body], a: u32, b: u32) {
        let (r0, well, k) = covalent(dominant(&bodies[a as usize]), dominant(&bodies[b as usize]));
        self.bonds.push(Bond::new(a, b, r0, well, k));
    }

    /// Constrain the angle at `b` to whatever it currently is.
    ///
    /// Taking the rest angle from the configuration rather than a table is the
    /// honest default: the geometry the generator produced is the geometry it
    /// meant, and inventing a tetrahedral angle for it would silently deform
    /// every molecule at the first step.
    pub fn bend(&mut self, bodies: &[Body], a: u32, b: u32, c: u32, stiffness: f64) {
        let u = bodies[a as usize].pos - bodies[b as usize].pos;
        let v = bodies[c as usize].pos - bodies[b as usize].pos;
        let rest = angle_between(u, v);
        self.angles.push(Angle { a, b, c, rest, stiffness });
    }

    /// Accelerations from the bonded terms alone.
    pub fn accelerations(&self, bodies: &[Body]) -> Vec<Vec3> {
        let mut acc = vec![Vec3::ZERO; bodies.len()];
        for f in self.forces(bodies).iter().enumerate() {
            let (i, force) = f;
            if bodies[i].mass > 0.0 {
                acc[i] = force.scale(1.0 / bodies[i].mass);
            }
        }
        acc
    }

    /// Forces from the bonded terms, in newtons.
    pub fn forces(&self, bodies: &[Body]) -> Vec<Vec3> {
        let n = bodies.len();
        let mut f = vec![Vec3::ZERO; n];
        for b in &self.bonds {
            let (i, j) = (b.a as usize, b.b as usize);
            if i >= n || j >= n {
                continue;
            }
            let d = bodies[j].pos - bodies[i].pos;
            let r = d.norm();
            if r <= 0.0 {
                continue;
            }
            // Equal and opposite along the line of centres, so the bonded
            // terms cannot move the centre of mass or add angular momentum.
            let pull = d.scale(b.tension(r) / r);
            f[i] += pull;
            f[j] -= pull;
        }
        for a in &self.angles {
            let (i, j, k) = (a.a as usize, a.b as usize, a.c as usize);
            if i >= n || j >= n || k >= n {
                continue;
            }
            let u = bodies[i].pos - bodies[j].pos;
            let v = bodies[k].pos - bodies[j].pos;
            let (lu, lv) = (u.norm(), v.norm());
            if lu <= 0.0 || lv <= 0.0 {
                continue;
            }
            let (uh, vh) = (u.scale(1.0 / lu), v.scale(1.0 / lv));
            let cos = uh.dot(vh).clamp(-1.0, 1.0);
            let sin = (1.0 - cos * cos).sqrt();
            // At exactly straight or exactly folded the bend direction is
            // undefined. It is also a stationary point of the potential, so
            // there is no force to miss by skipping it.
            if sin < 1.0e-7 {
                continue;
            }
            let theta = cos.acos();
            // F = -dV/dtheta * grad(theta), with V = k(theta - rest)^2 / 2.
            let dv = a.stiffness * (theta - a.rest);
            let grad_i = (vh - uh.scale(cos)).scale(-1.0 / (lu * sin));
            let grad_k = (uh - vh.scale(cos)).scale(-1.0 / (lv * sin));
            let fi = grad_i.scale(-dv);
            let fk = grad_k.scale(-dv);
            f[i] += fi;
            f[k] += fk;
            // The centre takes the reaction, which is what keeps the bend an
            // internal force.
            f[j] -= fi + fk;
        }
        f
    }

    /// Potential energy stored in the bonded terms, J.
    pub fn energy(&self, bodies: &[Body]) -> f64 {
        let n = bodies.len();
        let mut total = 0.0;
        for b in &self.bonds {
            let (i, j) = (b.a as usize, b.b as usize);
            if i >= n || j >= n {
                continue;
            }
            total += b.energy((bodies[j].pos - bodies[i].pos).norm());
        }
        for a in &self.angles {
            let (i, j, k) = (a.a as usize, a.b as usize, a.c as usize);
            if i >= n || j >= n || k >= n {
                continue;
            }
            let theta = angle_between(bodies[i].pos - bodies[j].pos, bodies[k].pos - bodies[j].pos);
            let d = theta - a.rest;
            total += 0.5 * a.stiffness * d * d;
        }
        total
    }

    /// Bonds stretched past the peak of their restoring force — on their way
    /// apart rather than merely stretched.
    pub fn dissociating(&self, bodies: &[Body]) -> Vec<usize> {
        self.bonds
            .iter()
            .enumerate()
            .filter(|(_, b)| {
                let (i, j) = (b.a as usize, b.b as usize);
                i < bodies.len()
                    && j < bodies.len()
                    && (bodies[j].pos - bodies[i].pos).norm() > b.inflection()
            })
            .map(|(i, _)| i)
            .collect()
    }

    /// Shortest vibrational period in the set, which is what bounds the
    /// timestep. A covalent bond is two orders of magnitude stiffer than the
    /// van der Waals interaction beside it, so a bonded system that reuses the
    /// unbonded timestep integrates its own bonds unstably.
    pub fn shortest_period(&self, bodies: &[Body]) -> f64 {
        let mut shortest = f64::INFINITY;
        for b in &self.bonds {
            let (i, j) = (b.a as usize, b.b as usize);
            if i >= bodies.len() || j >= bodies.len() {
                continue;
            }
            let (mi, mj) = (bodies[i].mass, bodies[j].mass);
            if mi <= 0.0 || mj <= 0.0 {
                continue;
            }
            let reduced = mi * mj / (mi + mj);
            let k = b.force_constant();
            if k > 0.0 {
                shortest = shortest.min(std::f64::consts::TAU * (reduced / k).sqrt());
            }
        }
        shortest
    }
}

fn angle_between(u: Vec3, v: Vec3) -> f64 {
    let (lu, lv) = (u.norm(), v.norm());
    if lu <= 0.0 || lv <= 0.0 {
        return 0.0;
    }
    (u.dot(v) / (lu * lv)).clamp(-1.0, 1.0).acos()
}

/// As [`step`], with covalent bonds.
///
/// Same velocity Verlet, with the bonded accelerations added to the pairwise
/// ones. Verlet is symplectic, so a bonded molecule's vibrational energy stays
/// bounded rather than drifting, which matters far more here than at the
/// unbonded tier: a bond oscillates a hundred times faster than anything else
/// in the system and would be the first thing to accumulate error.
pub fn step_bonded(
    bodies: &mut [Body],
    bonded: &Bonded,
    dt: f64,
    params: MdParams,
    world_seed: u64,
    path_key: u128,
    epoch: u32,
    tick: u64,
) -> SolveReport {
    if bonded.is_empty() {
        return step(bodies, dt, params, world_seed, path_key, epoch, tick);
    }
    let skip = bonded.exclusions(bodies.len());
    let potential = |b: &[Body]| {
        bonded.energy(b) + potential_energy_excluding(b, params, &skip)
    };
    let before = crate::solvers::measure(bodies, potential(bodies));
    let n = bodies.len();
    if n == 0 || dt == 0.0 {
        return SolveReport { before, after: before, dt_used: dt, ..Default::default() };
    }

    let acc = total_accelerations(bodies, bonded, params, &skip);
    for (b, a) in bodies.iter_mut().zip(&acc) {
        b.vel += a.scale(0.5 * dt);
        b.pos += b.vel.scale(dt);
    }
    let acc2 = total_accelerations(bodies, bonded, params, &skip);
    for (b, a) in bodies.iter_mut().zip(&acc2) {
        b.vel += a.scale(0.5 * dt);
    }

    let after = crate::solvers::measure(bodies, potential(bodies));
    SolveReport {
        steps: 1,
        interactions: n as u64,
        dt_used: dt,
        before,
        after,
        non_mechanical_energy: 0.0,
    }
}

fn total_accelerations(
    bodies: &[Body],
    bonded: &Bonded,
    params: MdParams,
    skip: &Exclusions,
) -> Vec<Vec3> {
    let mut acc = forces_excluding(bodies, params, skip);
    for (i, f) in bonded.forces(bodies).iter().enumerate() {
        if bodies[i].mass > 0.0 {
            acc[i] += f.scale(1.0 / bodies[i].mass);
        }
    }
    acc
}
