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

/// How much van der Waals interaction survives between each pair.
///
/// Two atoms held at a covalent bond length are far inside each other's
/// Lennard-Jones radius — a hydrogen pair sits at 0.74 angstroms with a sigma
/// of 2.57 — so the repulsive term between them is enormous and entirely
/// spurious. The bond already describes that interaction; leaving the van der
/// Waals term in as well double-counts it and then some, and a molecule
/// assembled with both blows itself apart in a few hundred femtoseconds.
///
/// Every force field removes bonded neighbours for this reason. A *reactive*
/// force field cannot do it with a switch that flips, because the same wall
/// that would tear a molecule apart also stops two atoms ever reaching the
/// distance at which they could bond: they meet hundreds of electron volts of
/// repulsion and bounce, and no chemistry ever happens. So a bonded pair's
/// dispersion is faded out over the same range the bond is faded in — zero
/// where the bond is doing the work, one where it is not. Pairs across an angle
/// are removed outright, as usual, since nothing brings them together or apart
/// except the bonds either side of them.
#[derive(Debug, Clone, Default)]
pub struct Exclusions {
    /// Per atom, sorted by partner: `(partner, shield_inner, shield_outer)`.
    /// A pair removed outright carries an infinite range, which weighs zero
    /// everywhere.
    by_atom: Vec<Vec<(u32, f64, f64)>>,
}

impl Exclusions {
    pub fn none() -> Exclusions {
        Exclusions::default()
    }

    /// Weight on the nonbonded interaction between `i` and `j` at separation
    /// `r`: one for an ordinary pair, zero for one the bonds already describe.
    #[inline]
    pub fn weight(&self, i: usize, j: u32, r: f64) -> f64 {
        self.weight_slope(i, j, r).0
    }

    /// The weight and its derivative with respect to separation.
    ///
    /// The derivative is not optional. A distance-dependent weight makes the
    /// pair potential `w(r) V(r)`, so the force is `w f - w' V`, and dropping
    /// the second term leaves a force that is not the gradient of anything —
    /// which shows up as steady, unexplained heating exactly where atoms are
    /// meeting each other.
    #[inline]
    pub fn weight_slope(&self, i: usize, j: u32, r: f64) -> (f64, f64) {
        let Some(list) = self.by_atom.get(i) else {
            return (1.0, 0.0);
        };
        match list.binary_search_by_key(&j, |e| e.0) {
            Ok(k) => {
                let (_, inner, outer) = list[k];
                let (s, ds) = switch(r, inner, outer);
                (1.0 - s, -ds)
            }
            Err(_) => (1.0, 0.0),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.by_atom.is_empty()
    }
}

impl Bonded {
    /// The 1-2 and 1-3 pairs this bond set implies, with their shielding.
    pub fn exclusions(&self, n: usize) -> Exclusions {
        let mut by_atom: Vec<Vec<(u32, f64, f64)>> = vec![Vec::new(); n];
        let mut add = |a: u32, b: u32, inner: f64, outer: f64| {
            if (a as usize) < n && (b as usize) < n && a != b {
                by_atom[a as usize].push((b, inner, outer));
                by_atom[b as usize].push((a, inner, outer));
            }
        };
        for b in &self.bonds {
            add(b.a, b.b, b.shield_inner, b.shield_outer);
        }
        for a in &self.angles {
            // The far pair of an angle is held by the two bonds either side and
            // never approaches on its own, so it comes out entirely.
            add(a.a, a.c, f64::INFINITY, f64::INFINITY);
        }
        for list in by_atom.iter_mut() {
            // Sort by partner, and where a pair appears twice keep the harder
            // exclusion — an angle's outright removal beats a bond's fade.
            list.sort_by(|x, y| x.0.cmp(&y.0).then(y.2.total_cmp(&x.2)));
            list.dedup_by_key(|e| e.0);
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
            if j == i {
                continue;
            }
            let bj = bodies[j];
            let d = bi.pos - bj.pos;
            let r = d.norm();
            if r <= 0.0 || r > params.cutoff {
                continue;
            }
            let (w, dw) = skip.weight_slope(i, jj, r);
            if w <= 0.0 && dw == 0.0 {
                continue;
            }
            let (sj, ej) = lj_params(dominant(&bj));
            // Lorentz-Berthelot mixing.
            let sigma = 0.5 * (si + sj);
            let epsilon = (ei * ej).sqrt();
            let mut f = lj_force(r, sigma, epsilon);
            let mut v = lj_potential(r, sigma, epsilon);
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
                v += K_COULOMB * bi.charge * bj.charge / r * scr;
            }
            if bi.mass > 0.0 {
                // -d(w V)/dr = w f - w' V.
                a += d.scale((w * f - dw * v) / (r * bi.mass));
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
            if j <= i {
                continue;
            }
            let bj = bodies[j];
            let r = (bi.pos - bj.pos).norm();
            if r <= 0.0 || r > params.cutoff {
                continue;
            }
            let w = skip.weight(i, jj, r);
            if w <= 0.0 {
                continue;
            }
            let (sj, ej) = lj_params(dominant(&bj));
            total += w * lj_potential(r, 0.5 * (si + sj), (ei * ej).sqrt());
            if bi.charge != 0.0 && bj.charge != 0.0 {
                total += w * K_COULOMB * bi.charge * bj.charge / r * (-r / params.debye).exp();
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

/// As [`stable_dt`], but bounded by the configuration the particles are
/// actually in rather than by the well they would sit in at equilibrium.
///
/// The vibrational period at the bottom of the Lennard-Jones well is the right
/// bound for a system that is *near* the bottom of it. A system that is not —
/// atoms packed closer than their own radii, which is what refining a solid
/// too far produces — sits on a wall where the force is five orders of
/// magnitude larger, and the equilibrium period is meaningless there.
///
/// The acceleration criterion is the standard answer: no particle may be
/// allowed to move more than a small fraction of its neighbour spacing in one
/// step, so `dt <= sqrt(2 eta d / a)` for the acceleration each one is actually
/// under. Without it, a compressed configuration does not integrate
/// inaccurately, it detonates — 10^25 m/s within two hundred steps, with the
/// conservation check reporting a drift of exactly 1.0 and nobody watching it.
pub fn configuration_dt(bodies: &[Body], params: MdParams) -> f64 {
    let base = stable_dt(bodies);
    if bodies.len() < 2 {
        return base;
    }
    let acc = forces(bodies, params);
    let mut spacing = f64::INFINITY;
    let grid = NeighbourGrid::build(bodies, params.cutoff);
    let mut nb = Vec::with_capacity(64);
    for i in 0..bodies.len() {
        grid.neighbours(bodies[i].pos, &mut nb);
        for &jj in nb.iter() {
            let j = jj as usize;
            if j == i {
                continue;
            }
            let r = (bodies[j].pos - bodies[i].pos).norm();
            if r > 0.0 {
                spacing = spacing.min(r);
            }
        }
    }
    if !spacing.is_finite() || spacing <= 0.0 {
        return base;
    }
    let mut limit = base;
    for a in &acc {
        let mag = a.norm();
        if mag > 0.0 {
            limit = limit.min((2.0 * 0.02 * spacing / mag).sqrt());
        }
    }
    limit.max(1e-24)
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
// Reactive chemistry
// ---------------------------------------------------------------------------

/// Covalent bonds that form, hold and break during the simulation.
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
/// # Morse, not harmonic
///
/// ```text
///     V(r) = f(r) D_e { [1 - e^{-a(r - r0)}]^2 - 1 },     a = sqrt(k / 2 D_e)
/// ```
///
/// A harmonic bond has the right stiffness near equilibrium and is infinitely
/// strong: it can be stretched across the simulation box and will still pull
/// back. Dissociation then has to be bolted on as a rule — break at some
/// extension, or above some energy — and the rule is a number nobody can
/// defend. The Morse potential has the same curvature at the bottom of the
/// well, so it reproduces the vibrational frequency exactly, and it rises to
/// zero at infinity from a depth of exactly `D_e`. A molecule given more than
/// its dissociation energy comes apart because the potential ran out, not
/// because a threshold fired.
///
/// # How a bond can appear without energy appearing with it
///
/// The hard part of reactive chemistry is not deciding when a bond exists. It
/// is that adding a term to the potential ordinarily *changes the potential*,
/// and a simulation that gains energy every time two atoms meet is worthless
/// however plausible its chemistry.
///
/// The factor `f(r)` above is what makes this safe. It is one out to `inner`,
/// zero beyond `outer`, and a smooth half-cosine between — so a bond at the
/// edge of its range contributes exactly zero energy and exactly zero force.
/// Bonds are created only inside that range and destroyed only outside it,
/// which means **the moment of creation and the moment of destruction are both
/// energetically invisible**. Conservation across a topology change is not
/// patched up afterwards; there is nothing to patch. The same switch multiplies
/// the angle terms, so a bend fades out with the bond that defined it.
///
/// What the bond then releases is real: forming an O–H bond drops the potential
/// by 4.8 eV and that energy arrives as kinetic energy of the pair. An
/// exothermic reaction warms the gas, and the engine's own conserved tuple
/// shows exactly where the heat came from.
#[derive(Debug, Clone, Default)]
pub struct Bonded {
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
}

/// A Morse bond between two particles, switched off smoothly at range.
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
    /// Separation below which the bond acts at full strength, m.
    pub inner: f64,
    /// Separation beyond which it contributes nothing at all, m.
    pub outer: f64,
    /// Below this the pair's van der Waals interaction is fully shielded, m.
    pub shield_inner: f64,
    /// Above this it acts at full strength, m. This is also the range at which
    /// the pair is considered a bond at all.
    pub shield_outer: f64,
}

/// Where a bond's two switches sit: `(inner, outer, shield_outer)`.
///
/// Both ranges are derived from where the potentials are actually negligible
/// rather than from a round multiple of the bond length, because the two scales
/// are unrelated and the gap between them is where reactive chemistry either
/// works or does not.
///
/// * The Morse switch begins where the potential is within 2% of zero — at
///   `r0 + ln(100)/alpha` — so what the switch reshapes is a tail worth a
///   fiftieth of the well and not the well itself.
/// * It ends, and the dispersion shield begins, past the Lennard-Jones minimum.
///   The repulsive wall of a pair at covalent distance is hundreds of electron
///   volts; the two atoms would bounce off it long before they reached anything
///   the bond could hold them at, and no chemistry would ever happen. Handing
///   over at the minimum means neither term is doing anything much at the
///   crossover, so the handover leaves no barrier of its own.
///
/// An earlier version handed over at 2.4 bond lengths, which for hydrogen is
/// deep inside the van der Waals wall. It left a 0.026 eV activation barrier
/// that was nobody's chemistry — the atoms could not reach each other, and the
/// only sign of it was that a test which expected a bond never got one.
fn ranges(r0: f64, alpha: f64, sigma: f64) -> (f64, f64, f64) {
    let negligible = if alpha > 0.0 { r0 + 100f64.ln() / alpha } else { 2.4 * r0 };
    let minimum = 2f64.powf(1.0 / 6.0) * sigma;
    let outer = (1.1 * negligible).max(1.25 * minimum);
    let inner = negligible.min(0.85 * outer).max(1.6 * r0);
    (inner, outer, 1.35 * outer)
}

/// The switch, and its derivative.
///
/// A raised cosine: one below `inner`, zero above `outer`, with zero slope at
/// both ends so the *force* is continuous as well as the energy. A linear ramp
/// would be continuous in energy and leave a step in the force, which shows up
/// as a slow and entirely artificial heating.
#[inline]
fn switch(r: f64, inner: f64, outer: f64) -> (f64, f64) {
    if r <= inner {
        (1.0, 0.0)
    } else if r >= outer {
        (0.0, 0.0)
    } else {
        let span = outer - inner;
        let x = std::f64::consts::PI * (r - inner) / span;
        (0.5 * (1.0 + x.cos()), -0.5 * std::f64::consts::PI / span * x.sin())
    }
}

impl Bond {
    /// A bond with a given harmonic force constant `k` (N/m) and well depth.
    ///
    /// The switch runs from 1.6 to 2.4 times the equilibrium length: far enough
    /// out that the well depth and the vibrational frequency are untouched —
    /// `f = 1` well past where the potential has any curvature left — and short
    /// enough that the candidate search stays local.
    pub fn new(a: u32, b: u32, r0: f64, well: f64, k: f64, sigma: f64) -> Bond {
        let alpha = if well > 0.0 { (k / (2.0 * well)).sqrt() } else { 0.0 };
        let (inner, outer, shield_outer) = ranges(r0, alpha, sigma);
        Bond { a, b, r0, well, alpha, inner, outer, shield_inner: outer, shield_outer }
    }

    /// How much of the pair's van der Waals interaction survives at `r`.
    ///
    /// Zero where the bond is doing the work, one where it is not, and a smooth
    /// ramp between — so the two descriptions hand over continuously instead of
    /// double-counting inside the bond and leaving a cliff outside it.
    #[inline]
    pub fn dispersion_weight(&self, r: f64) -> f64 {
        1.0 - switch(r, self.shield_inner, self.shield_outer).0
    }

    /// Harmonic force constant at the bottom of the well, N/m.
    #[inline]
    pub fn force_constant(&self) -> f64 {
        2.0 * self.well * self.alpha * self.alpha
    }

    /// Morse potential referenced to the dissociation limit: `-D_e` at rest,
    /// zero at infinity. Unswitched; [`Bond::energy`] applies the switch.
    #[inline]
    pub fn morse(&self, r: f64) -> f64 {
        let x = 1.0 - (-self.alpha * (r - self.r0)).exp();
        self.well * (x * x - 1.0)
    }

    /// Potential energy at separation `r`, switched off at range.
    #[inline]
    pub fn energy(&self, r: f64) -> f64 {
        let (s, _) = switch(r, self.inner, self.outer);
        if s == 0.0 {
            0.0
        } else {
            s * self.morse(r)
        }
    }

    /// Attractive force magnitude at `r`. Positive pulls the pair together.
    #[inline]
    pub fn tension(&self, r: f64) -> f64 {
        let (s, ds) = switch(r, self.inner, self.outer);
        if s == 0.0 && ds == 0.0 {
            return 0.0;
        }
        let e = (-self.alpha * (r - self.r0)).exp();
        let dmorse = 2.0 * self.alpha * self.well * (1.0 - e) * e;
        // Tension is dV/dr, and V = s(r) morse(r), so the product rule gives
        // both terms. The switch term is not a correction to be subtracted:
        // where the switch is closing, it is *lifting* a negative potential
        // towards zero, which is extra restoring force, and getting its sign
        // wrong makes the bond weakest exactly where it is being asked to hold.
        s * dmorse + ds * self.morse(r)
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
/// infinity. They carry the two bonds' switch ranges so that when either bond
/// leaves range, the bend leaves with it, continuously.
#[derive(Debug, Clone, Copy)]
pub struct Angle {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    /// Rest angle at `b`, radians.
    pub rest: f64,
    /// Bending constant, J/rad^2.
    pub stiffness: f64,
    /// Switch range of the `b-a` bond.
    pub inner_a: f64,
    pub outer_a: f64,
    /// Switch range of the `b-c` bond.
    pub inner_c: f64,
    pub outer_c: f64,
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

/// How many covalent bonds a species will hold.
///
/// This is what stops a hydrogen atom acquiring five neighbours. The
/// alternative — a many-body bond order that weakens every bond as an atom
/// becomes over-coordinated — is what a Tersoff or Brenner potential does, and
/// it is the right answer for a solver that has to get the energetics of
/// intermediate coordination right. The question being asked here is only
/// whether a molecule holds together and can react, and a valence count answers
/// it with a rule nobody has to calibrate.
pub fn valence(s: Species) -> usize {
    match s {
        Species::Hydrogen => 1,
        Species::Helium => 0,
        Species::Carbon => 4,
        Species::Nitrogen => 3,
        Species::Oxygen => 2,
        Species::Silicon => 4,
        Species::Iron => 6,
        Species::Other => 2,
    }
}

/// Rest angle at an atom holding `bonds` bonds, radians.
///
/// Electron-pair repulsion, with the two cases where a lone pair closes the
/// angle down from the ideal. That is the difference between water at 104.5
/// degrees and a linear triatomic, and it is visible in every property water
/// has.
pub fn bond_angle(s: Species, bonds: usize) -> f64 {
    let degrees: f64 = match (s, bonds) {
        (Species::Oxygen, 2) => 104.5,
        (Species::Nitrogen, 3) => 107.0,
        (_, 0 | 1 | 2) => 180.0,
        (_, 3) => 120.0,
        _ => 109.4712206,
    };
    degrees.to_radians()
}

/// Bending force constant at an atom, J/rad^2. Measured values; water's bend at
/// 0.70 aJ/rad^2 is the one most people would recognise.
pub fn bend_constant(s: Species) -> f64 {
    match s {
        Species::Oxygen => 4.37 * EV,
        Species::Nitrogen => 4.00 * EV,
        Species::Carbon => 3.90 * EV,
        Species::Silicon => 2.20 * EV,
        _ => 3.00 * EV,
    }
}

/// What one pass of chemistry did.
#[derive(Debug, Clone, Copy, Default)]
pub struct Reaction {
    pub formed: usize,
    pub broken: usize,
    /// Potential energy the pass itself changed, J. It must be zero — bonds are
    /// only ever created and destroyed where their switch is — and it is
    /// reported so that "must be" is something a test can check.
    pub energy_change: f64,
}

impl Bonded {
    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty() && self.angles.is_empty()
    }

    /// Bond two particles using the constants for their dominant species.
    pub fn bond(&mut self, bodies: &[Body], a: u32, b: u32) {
        let (sa, sb) = (dominant(&bodies[a as usize]), dominant(&bodies[b as usize]));
        let (r0, well, k) = covalent(sa, sb);
        let sigma = 0.5 * (lj_params(sa).0 + lj_params(sb).0);
        self.bonds.push(Bond::new(a, b, r0, well, k, sigma));
    }

    /// The range at which a pair of species starts counting as bonded.
    pub fn capture_radius(a: Species, b: Species) -> f64 {
        let (r0, well, k) = covalent(a, b);
        let alpha = if well > 0.0 { (k / (2.0 * well)).sqrt() } else { 0.0 };
        let sigma = 0.5 * (lj_params(a).0 + lj_params(b).0);
        ranges(r0, alpha, sigma).2
    }

    /// Constrain the angle at `b` to whatever it currently is.
    ///
    /// Taking the rest angle from the configuration rather than a table is the
    /// honest default for a bend imposed by hand: the geometry the caller
    /// produced is the geometry it meant, and inventing a tetrahedral angle for
    /// it would silently deform the molecule on the first step. Angles built by
    /// [`Bonded::react`] use the species' own rest angle instead, because there
    /// the geometry is an accident of how the atoms happened to meet.
    pub fn bend(&mut self, bodies: &[Body], a: u32, b: u32, c: u32, stiffness: f64) {
        let u = bodies[a as usize].pos - bodies[b as usize].pos;
        let v = bodies[c as usize].pos - bodies[b as usize].pos;
        let rest = angle_between(u, v);
        let (ia, oa) = self.range_of(a, b);
        let (ic, oc) = self.range_of(b, c);
        self.angles.push(Angle {
            a,
            b,
            c,
            rest,
            stiffness,
            inner_a: ia,
            outer_a: oa,
            inner_c: ic,
            outer_c: oc,
        });
    }

    fn range_of(&self, a: u32, b: u32) -> (f64, f64) {
        for bond in &self.bonds {
            if (bond.a == a && bond.b == b) || (bond.a == b && bond.b == a) {
                return (bond.inner, bond.outer);
            }
        }
        // No such bond: a switch that never switches, so a hand-placed bend
        // behaves exactly as it did before ranges existed.
        (f64::INFINITY, f64::INFINITY)
    }

    /// Number of bonds currently held by each atom.
    pub fn coordination(&self, n: usize) -> Vec<usize> {
        let mut z = vec![0usize; n];
        for b in &self.bonds {
            if (b.a as usize) < n && (b.b as usize) < n {
                z[b.a as usize] += 1;
                z[b.b as usize] += 1;
            }
        }
        z
    }

    /// Let the chemistry change.
    ///
    /// Bonds past their outer range are removed, and new ones form between
    /// atoms that have come within range and still have valence free. Both
    /// happen where the switch is zero, so neither changes the energy — which
    /// is reported rather than assumed, and asserted in `tests/bonded.rs`.
    ///
    /// The closest pairs bond first, so the outcome does not depend on the
    /// order atoms happen to sit in the array. Ties break on index, so it
    /// replays exactly.
    pub fn react(&mut self, bodies: &[Body]) -> Reaction {
        let n = bodies.len();
        let before = self.energy(bodies);
        let mut report = Reaction::default();

        // Anything past its range is holding nothing.
        let kept: Vec<Bond> = self
            .bonds
            .iter()
            .copied()
            .filter(|b| {
                let (i, j) = (b.a as usize, b.b as usize);
                i < n && j < n && (bodies[j].pos - bodies[i].pos).norm() < b.shield_outer
            })
            .collect();
        report.broken = self.bonds.len() - kept.len();
        self.bonds = kept;

        // Candidate pairs, from the same cell lists the nonbonded sum uses.
        let mut reach: f64 = 0.0;
        for b in bodies.iter() {
            let si = dominant(b);
            for s in Species::ALL {
                reach = reach.max(Bonded::capture_radius(si, s));
            }
        }
        if reach <= 0.0 || n == 0 {
            report.energy_change = self.energy(bodies) - before;
            return report;
        }
        let grid = NeighbourGrid::build(bodies, reach);
        let mut nb = Vec::with_capacity(64);
        let mut z = self.coordination(n);
        let mut bonded: std::collections::HashSet<(u32, u32)> = self
            .bonds
            .iter()
            .map(|b| (b.a.min(b.b), b.a.max(b.b)))
            .collect();

        let mut candidates: Vec<(f64, u32, u32)> = Vec::new();
        for i in 0..n {
            if valence(dominant(&bodies[i])) == 0 {
                continue;
            }
            grid.neighbours(bodies[i].pos, &mut nb);
            for &jj in nb.iter() {
                let j = jj as usize;
                if j <= i {
                    continue;
                }
                if bonded.contains(&(i as u32, jj)) {
                    continue;
                }
                let capture =
                    Bonded::capture_radius(dominant(&bodies[i]), dominant(&bodies[j]));
                let r = (bodies[j].pos - bodies[i].pos).norm();
                if r < capture {
                    candidates.push((r, i as u32, jj));
                }
            }
        }
        candidates.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));

        for (_, i, j) in candidates {
            let (a, b) = (i as usize, j as usize);
            if z[a] >= valence(dominant(&bodies[a])) || z[b] >= valence(dominant(&bodies[b])) {
                continue;
            }
            self.bond(bodies, i, j);
            bonded.insert((i, j));
            z[a] += 1;
            z[b] += 1;
            report.formed += 1;
        }

        if report.formed > 0 || report.broken > 0 {
            self.rebuild_angles(bodies);
        }
        report.energy_change = self.energy(bodies) - before;
        report
    }

    /// Rebuild every bend from the current bond list.
    ///
    /// Rest angles come from the species and its coordination rather than from
    /// the geometry: two hydrogens that have just found an oxygen are wherever
    /// they happened to arrive, and the molecule's job is to pull them to 104.5
    /// degrees, not to memorise the accident.
    pub fn rebuild_angles(&mut self, bodies: &[Body]) {
        let n = bodies.len();
        let mut neighbours: Vec<Vec<(u32, f64, f64)>> = vec![Vec::new(); n];
        for b in &self.bonds {
            let (i, j) = (b.a as usize, b.b as usize);
            if i >= n || j >= n {
                continue;
            }
            neighbours[i].push((b.b, b.inner, b.outer));
            neighbours[j].push((b.a, b.inner, b.outer));
        }
        self.angles.clear();
        for centre in 0..n {
            let list = &neighbours[centre];
            if list.len() < 2 {
                continue;
            }
            let species = dominant(&bodies[centre]);
            let rest = bond_angle(species, list.len());
            let stiffness = bend_constant(species);
            for x in 0..list.len() {
                for y in x + 1..list.len() {
                    self.angles.push(Angle {
                        a: list[x].0,
                        b: centre as u32,
                        c: list[y].0,
                        rest,
                        stiffness,
                        inner_a: list[x].1,
                        outer_a: list[x].2,
                        inner_c: list[y].1,
                        outer_c: list[y].2,
                    });
                }
            }
        }
    }

    /// Accelerations from the bonded terms alone.
    pub fn accelerations(&self, bodies: &[Body]) -> Vec<Vec3> {
        let mut acc = vec![Vec3::ZERO; bodies.len()];
        for (i, force) in self.forces(bodies).iter().enumerate() {
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
            // Equal and opposite along the line of centres, so the bonded terms
            // cannot move the centre of mass or add angular momentum.
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
            let (su, dsu) = switch(lu, a.inner_a, a.outer_a);
            let (sv, dsv) = switch(lv, a.inner_c, a.outer_c);
            let s = su * sv;
            if s == 0.0 && dsu == 0.0 && dsv == 0.0 {
                continue;
            }
            let (uh, vh) = (u.scale(1.0 / lu), v.scale(1.0 / lv));
            let cos = uh.dot(vh).clamp(-1.0, 1.0);
            let sin = (1.0 - cos * cos).sqrt();
            let delta = cos.acos() - a.rest;
            let bend = 0.5 * a.stiffness * delta * delta;

            // Radial part: the switch fading the bend in and out. It acts along
            // each bond, so it stays internal like everything else here.
            let radial_i = uh.scale(-dsu * sv * bend);
            let radial_k = vh.scale(-dsv * su * bend);
            f[i] += radial_i;
            f[k] += radial_k;
            f[j] -= radial_i + radial_k;

            // Angular part. At exactly straight or exactly folded the bend
            // direction is undefined; it is also a stationary point of the
            // potential, so there is no force to miss by skipping it.
            if sin < 1.0e-7 || s == 0.0 {
                continue;
            }
            // F = -dV/dtheta grad(theta), with V = s k (theta - rest)^2 / 2.
            let dv = s * a.stiffness * delta;
            let grad_i = (vh - uh.scale(cos)).scale(-1.0 / (lu * sin));
            let grad_k = (uh - vh.scale(cos)).scale(-1.0 / (lv * sin));
            let fi = grad_i.scale(-dv);
            let fk = grad_k.scale(-dv);
            f[i] += fi;
            f[k] += fk;
            // The centre takes the reaction, which keeps the bend internal.
            f[j] -= fi + fk;
        }
        f
    }

    /// Potential energy stored in the bonded terms, J.
    ///
    /// Referenced to infinite separation, so a bound molecule has *negative*
    /// bonded energy and forming a bond releases exactly its well depth into
    /// the particles' motion.
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
            let u = bodies[i].pos - bodies[j].pos;
            let v = bodies[k].pos - bodies[j].pos;
            let (su, _) = switch(u.norm(), a.inner_a, a.outer_a);
            let (sv, _) = switch(v.norm(), a.inner_c, a.outer_c);
            if su * sv == 0.0 {
                continue;
            }
            let d = angle_between(u, v) - a.rest;
            total += su * sv * 0.5 * a.stiffness * d * d;
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

    /// Shortest period the chemistry *could* produce, whether or not any bond
    /// exists yet. A reactive run has to be integrated at a timestep the bonds
    /// it is about to form will survive, not the one its current bonds need.
    pub fn reachable_period(bodies: &[Body]) -> f64 {
        let mut shortest = f64::INFINITY;
        let mut seen: Vec<Species> = Vec::new();
        for b in bodies.iter() {
            let s = dominant(b);
            if !seen.contains(&s) {
                seen.push(s);
            }
        }
        for b in bodies.iter() {
            if b.mass <= 0.0 {
                continue;
            }
            for other in bodies.iter() {
                if other.mass <= 0.0 {
                    continue;
                }
                let (_, well, k) = covalent(dominant(b), dominant(other));
                if k <= 0.0 || well <= 0.0 {
                    continue;
                }
                let reduced = b.mass * other.mass / (b.mass + other.mass);
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

/// Stable timestep for a *reactive* system, bounded by the stiffest bond the
/// atoms present could form rather than by the ones they already have.
pub fn stable_dt_reactive(bodies: &[Body]) -> f64 {
    let period = Bonded::reachable_period(bodies);
    if period.is_finite() && period > 0.0 {
        stable_dt(bodies).min(period / 50.0)
    } else {
        stable_dt(bodies)
    }
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
    let potential = |b: &[Body]| bonded.energy(b) + potential_energy_excluding(b, params, &skip);
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

/// A step of *reactive* molecular dynamics: integrate, then let the chemistry
/// change.
///
/// The reaction pass runs after the integration so that the forces acting over
/// the step are the ones the report's conserved tuple was measured against.
/// Because bonds are created and destroyed only where their switch is zero, the
/// pass cannot change the potential energy, which is what makes it safe to run
/// inside the step at all.
pub fn step_reactive(
    bodies: &mut [Body],
    bonded: &mut Bonded,
    dt: f64,
    params: MdParams,
    world_seed: u64,
    path_key: u128,
    epoch: u32,
    tick: u64,
) -> (SolveReport, Reaction) {
    let report = if bonded.is_empty() {
        bonded.react(bodies);
        step_bonded(bodies, bonded, dt, params, world_seed, path_key, epoch, tick)
    } else {
        step_bonded(bodies, bonded, dt, params, world_seed, path_key, epoch, tick)
    };
    let reaction = bonded.react(bodies);
    (report, reaction)
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
