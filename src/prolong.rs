//! Prolongation: manufacturing fine detail that is *provably* consistent with
//! the coarse state it came from.
//!
//! # The problem
//!
//! A node says "I am 3x10^4 solar masses of 20 K molecular gas, 12 pc across,
//! spinning this way, moving that way". The camera zooms in. We now need a
//! million gas parcels. They must:
//!
//! 1. reproduce the parent's conserved tuple **exactly** — not statistically,
//!    exactly, or energy and momentum drift every time the user pans;
//! 2. look right — the correct density profile, the correct velocity
//!    distribution, the correct correlations, or the deception is visible;
//! 3. be **reproducible** — pan away, come back, get the same million parcels,
//!    or the world is not a world.
//!
//! # The method
//!
//! Sample from the maximum-entropy distribution consistent with the coarse
//! state (that is requirement 2), then apply an exact linear projection onto
//! the constraint surface (requirement 1). Requirement 3 comes free from
//! `rng.rs`: the sample is a pure function of the node's path key.
//!
//! The projection is the interesting part. Done naively, fixing the momentum
//! breaks the angular momentum, fixing the angular momentum breaks the energy,
//! and so on around the loop forever. The fix is to build the corrections in an
//! order where each one lives in the null space of the previous constraints:
//!
//! ```text
//!   centre positions        =>  sum m r = 0
//!   subtract mean velocity  =>  sum m v = 0
//!   subtract rigid rotation =>  sum m r x v = 0        (leaves sum m v = 0)
//!   scale residual by s     =>  both still 0           (scaling is linear)
//!   add rigid rotation w_t  =>  L = L_target exactly   (adds no momentum)
//!   add bulk drift v_b      =>  P = P_target exactly   (adds no L about com)
//! ```
//!
//! Each step is in the kernel of the constraints already satisfied, so nothing
//! is ever undone. And because the residual field δ has zero angular momentum
//! by construction, it is *energetically orthogonal* to the rigid rotation we
//! add — `Σ m δ·(ω×r) = ω·Σ m r×δ = ω·L_res = 0` — which is what lets us solve
//! for the energy scale `s` in closed form instead of iterating.

use crate::math::{det_sum_by, det_sum_v3_by, Vec3};
use crate::rng::{Purpose, Stream};
use crate::state::{mutual_gravitational_energy, Aggregate, Body, BodyKind, Composition};
use crate::units::*;

/// Spatial arrangement of the children. Chosen by tier and by what the node is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Profile {
    /// Uniform ball. Gas parcels, unstructured debris.
    Uniform,
    /// Plummer sphere — the standard equilibrium model for a star cluster or a
    /// dark matter subhalo. Finite central density, so no cusp singularity.
    Plummer,
    /// Exponential disk with a sech^2 vertical profile. Galactic tier.
    Disk { scale_height_ratio: f64 },
    /// Thin shell. Supernova remnants, electron shells.
    Shell,
    /// Woods-Saxon: the empirical nuclear density profile.
    WoodsSaxon,
    /// Regular lattice with thermal displacement. Solids.
    Lattice,
}

/// How mass is divided among the children.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MassSpectrum {
    /// Equal shares. Statistical super-particles.
    Equal,
    /// Kroupa (2001) initial mass function. Turns a cloud into a cluster.
    Kroupa { min_msun: f64, max_msun: f64 },
    /// `dN/dm ∝ m^alpha`. Cloud fragmentation, debris, dust grains.
    PowerLaw { alpha: f64, ratio: f64 },
    /// Masses set by the composition — one child per nucleus, correct
    /// per-species masses. Used at molecular and finer tiers.
    Species,
}

/// Everything needed to turn one aggregate into many bodies.
#[derive(Debug, Clone, Copy)]
pub struct ProlongSpec {
    pub count: usize,
    pub profile: Profile,
    pub spectrum: MassSpectrum,
    pub kind: BodyKind,
    /// Fractional scatter applied to per-child composition before exact
    /// renormalisation. Zero gives chemically identical children.
    pub composition_scatter: f64,
    /// Fraction of the velocity dispersion that is coherent turbulence rather
    /// than thermal noise. Drives the ISM's Larson relations.
    pub turbulent_fraction: f64,
}

impl ProlongSpec {
    pub fn new(count: usize, profile: Profile, spectrum: MassSpectrum, kind: BodyKind) -> Self {
        ProlongSpec {
            count,
            profile,
            spectrum,
            kind,
            composition_scatter: 0.0,
            turbulent_fraction: 0.0,
        }
    }
}

/// What actually happened during a prolongation. The engine keeps these; the
/// consistency tests assert on them; the debug UI shows them.
#[derive(Debug, Clone, Copy)]
pub struct ProlongReport {
    pub count: usize,
    /// Worst relative error across the conserved tuple. Must stay at round-off.
    pub conservation_error: f64,
    /// True if the requested radius was overridden because the requested
    /// angular momentum could not fit inside it at the requested energy.
    pub radius_overridden: bool,
    pub realised_radius: f64,
    /// Fraction of internal energy that ended up as coherent rotation.
    pub rotational_fraction: f64,
    /// Set when the target state was thermodynamically impossible and had to be
    /// relaxed. A non-empty reason is a bug in whatever produced the aggregate.
    pub relaxations: u32,
    /// How much of the energy budget had to be absorbed into the children's
    /// internal account rather than their motion, as a fraction of the
    /// parent's internal energy. Near zero when the velocity projection
    /// converged cleanly; a large value means the requested state was awkward
    /// (extreme mass ratios, near-degenerate geometry) but still exact.
    pub internal_energy_residual: f64,
    /// Factor applied to member radii to make the enclosed volume agree with
    /// the structural mass at the material's density. Far from 1 means the
    /// program's nominal proportions disagree with its own density.
    pub radius_correction: f64,
    /// What the on-creation design pass achieved.
    pub design: crate::solvers::structure::DesignReport,
    /// How many of the emitted parts belong to the structure itself, as
    /// opposed to the unstructured matter surrounding it in the same node.
    pub structural_parts: usize,
    /// Natural magnitudes the error above is measured against.
    pub scales: crate::state::Scales,
    /// The self-potential the sampler settled on. The engine must use this same
    /// value when restricting these bodies, or the two directions are measuring
    /// different quantities and the round trip is not a round trip.
    pub potential: f64,
}

impl Default for ProlongReport {
    fn default() -> Self {
        ProlongReport {
            count: 0,
            conservation_error: 0.0,
            radius_overridden: false,
            realised_radius: 0.0,
            rotational_fraction: 0.0,
            relaxations: 0,
            radius_correction: 1.0,
            design: crate::solvers::structure::DesignReport::default(),
            structural_parts: 0,
            internal_energy_residual: 0.0,
            scales: crate::state::Scales::unit(),
            potential: 0.0,
        }
    }
}

/// Inertia tensor straight from slices, with no intermediate allocation.
fn inertia_of_slices(pos: &[Vec3], masses: &[f64]) -> crate::math::Mat3 {
    let mut m = crate::math::Mat3::zero();
    for (p, &mass) in pos.iter().zip(masses) {
        let r2 = p.norm2();
        let o = p.outer(*p);
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                m.0[i][j] += mass * (r2 * delta - o.0[i][j]);
            }
        }
    }
    m
}

/// The prolongation operator P.
///
/// Deterministic in `(world_seed, path_key, epoch)` — call it a thousand times,
/// on a thousand machines, get the same bodies.
pub fn prolong(
    agg: &Aggregate,
    spec: ProlongSpec,
    world_seed: u64,
    path_key: u128,
    epoch: u32,
) -> (Vec<Body>, ProlongReport) {
    let n = spec.count.max(1);
    let mut report = ProlongReport {
        count: n,
        ..Default::default()
    };

    if agg.mass <= 0.0 || !agg.is_finite() {
        return (Vec::new(), report);
    }

    // ---- 1. masses ------------------------------------------------------
    let mut masses = sample_masses(agg, spec, n, world_seed, path_key, epoch);
    let m_sum = det_sum_by(n, &|i| masses[i]);
    let k = agg.mass / m_sum;
    for m in masses.iter_mut() {
        *m *= k; // exact by construction: sum is now agg.mass to round-off
    }

    // ---- 2. positions ---------------------------------------------------
    let mut pos = sample_positions(spec, n, world_seed, path_key, epoch);

    // Centre so that sum m r = 0, then scale to hit the requested radius.
    recentre(&mut pos, &masses, agg.mass);
    let mut scale = radius_scale(&pos, &masses, agg.mass, agg.radius);
    for p in pos.iter_mut() {
        *p = p.scale(scale);
    }

    // ---- 3. energy budget, with the potential re-derived ----------------
    //
    // The parent's *total* energy is the invariant. How it splits between
    // thermal, rotational and gravitational is not — that split is re-derived
    // from the geometry we just sampled, and whatever the potential turns out
    // to be, the kinetic budget absorbs the difference. This is why zooming in
    // and out does not leak energy even though the children's positions were
    // never stored: we are not required to reproduce the old potential, only
    // the old total.
    let softening = agg.radius / (n as f64).cbrt() * 0.1;
    let mut phi = potential_estimate(&pos, &masses, softening, spec.profile, agg);
    let mut random_ke_target = agg.internal_energy + agg.binding_energy - phi;

    // A configuration can be too tightly bound to hold the energy it claims.
    // Physically the answer is that it must be bigger, so make it bigger.
    let mut guard = 0;
    while random_ke_target <= 0.0 && guard < 32 {
        for p in pos.iter_mut() {
            *p = p.scale(1.5);
        }
        scale *= 1.5;
        phi = potential_estimate(&pos, &masses, softening * scale, spec.profile, agg);
        random_ke_target = agg.internal_energy + agg.binding_energy - phi;
        report.radius_overridden = true;
        report.relaxations += 1;
        guard += 1;
    }
    if random_ke_target <= 0.0 {
        random_ke_target = agg.internal_energy.abs().max(1e-30);
        report.relaxations += 1;
    }

    // ---- 4. rotation feasibility ---------------------------------------
    //
    // L^2/(2I) is the energy locked up in rigid rotation. If that exceeds the
    // whole random-energy budget, the body cannot be this compact while
    // carrying this much spin — so, again, it gets bigger. (A real object in
    // this situation would shed mass; that belongs to the tier solver, not to
    // the sampler, so here we simply refuse to build an impossible state.)
    let mut inertia = inertia_tensor_of(&pos, &masses);
    let mut omega = inertia.solve(agg.spin).unwrap_or(Vec3::ZERO);
    let mut ke_rot = 0.5 * omega.dot(agg.spin);
    guard = 0;
    while ke_rot > 0.95 * random_ke_target && ke_rot > 0.0 && guard < 32 {
        let need = (ke_rot / (0.5 * random_ke_target)).sqrt().max(1.2);
        for p in pos.iter_mut() {
            *p = p.scale(need);
        }
        scale *= need;
        phi = potential_estimate(&pos, &masses, softening * scale, spec.profile, agg);
        random_ke_target =
            (agg.internal_energy + agg.binding_energy - phi).max(agg.internal_energy.abs().max(1e-30));
        inertia = inertia_tensor_of(&pos, &masses);
        omega = inertia.solve(agg.spin).unwrap_or(Vec3::ZERO);
        ke_rot = 0.5 * omega.dot(agg.spin);
        report.radius_overridden = true;
        report.relaxations += 1;
        guard += 1;
    }

    // ---- 5. velocities: sample, project, then polish --------------------
    //
    // The projection below is exact in Newtonian mechanics. The engine's
    // functionals are relativistic (`p = gamma m v`, `K = (gamma-1) m c^2`),
    // so the closed-form answer is a first guess that is wrong at O(v^2/c^2) —
    // one part in 10^6 for a galactic disk, which is far too large to accept
    // in a quantity that gets round-tripped thousands of times. So we finish
    // with a fixed-point polish against the *exact* functionals, the same ones
    // `restrict` uses. It converges in three passes because each correction
    // lives in the null space of the others.
    let masses_for_comp = masses.clone();
    let resid = sample_velocities(agg, spec, n, &pos, world_seed, path_key, epoch);

    let parts = Parts {
        pos,
        masses,
        comps: sample_compositions(agg, spec, &masses_for_comp, world_seed, path_key, epoch),
        radii: Vec::new(),
        kind: spec.kind,
    };
    let bodies = close_books(agg, parts, resid, phi, random_ke_target, omega, ke_rot, scale, &mut report);

    (bodies, report)
}

/// Inverse of `prolong` for the purposes of the round trip: see
/// `state::restrict`. Kept here as documentation of the pairing.
pub use crate::state::restrict;

// ---------------------------------------------------------------------------
// samplers
// ---------------------------------------------------------------------------

fn sample_masses(
    agg: &Aggregate,
    spec: ProlongSpec,
    n: usize,
    seed: u64,
    key: u128,
    epoch: u32,
) -> Vec<f64> {
    let mut st = Stream::at(seed, key, epoch, Purpose::Masses);
    let mut m = Vec::with_capacity(n);
    match spec.spectrum {
        MassSpectrum::Equal => m.resize(n, agg.mass / n as f64),
        MassSpectrum::Kroupa { min_msun, max_msun } => {
            for _ in 0..n {
                m.push(kroupa_sample(&mut st, min_msun, max_msun) * M_SUN);
            }
        }
        MassSpectrum::PowerLaw { alpha, ratio } => {
            let lo = agg.mass / n as f64 / ratio.max(1.0);
            let hi = lo * ratio.max(1.0);
            for _ in 0..n {
                m.push(st.power_law(lo, hi, alpha));
            }
        }
        MassSpectrum::Species => {
            // One child per nucleus, drawn from the composition. The masses are
            // then physical, not statistical.
            let mut weights = [0.0; NSPECIES];
            for (i, s) in Species::ALL.iter().enumerate() {
                weights[i] = agg.composition.get(*s) / s.mass_kg();
            }
            for _ in 0..n {
                let i = st.weighted(&weights);
                m.push(Species::ALL[i].mass_kg());
            }
        }
    }
    m
}

/// Kroupa (2001) three-part power law, sampled by inverse CDF.
fn kroupa_sample(st: &mut Stream, lo: f64, hi: f64) -> f64 {
    let breaks = [0.08f64, 0.5f64];
    let slopes = [-0.3f64, -1.3f64, -2.3f64];
    let seg = [
        (lo.max(0.01), breaks[0].min(hi), slopes[0]),
        (breaks[0].max(lo), breaks[1].min(hi), slopes[1]),
        (breaks[1].max(lo), hi.max(breaks[1]), slopes[2]),
    ];
    let mut weights = [0.0f64; 3];
    for (i, (a, b, al)) in seg.iter().enumerate() {
        if b > a {
            let p = al + 1.0;
            weights[i] = if p.abs() < 1e-9 {
                (b / a).ln()
            } else {
                (b.powf(p) - a.powf(p)) / p
            };
        }
    }
    let i = st.weighted(&weights);
    let (a, b, al) = seg[i];
    if b <= a {
        return a;
    }
    st.power_law(a, b, al)
}

fn sample_positions(spec: ProlongSpec, n: usize, seed: u64, key: u128, epoch: u32) -> Vec<Vec3> {
    let mut st = Stream::at(seed, key, epoch, Purpose::Positions);
    let mut out = Vec::with_capacity(n);
    match spec.profile {
        Profile::Uniform => {
            for _ in 0..n {
                out.push(st.in_ball());
            }
        }
        Profile::Plummer => {
            for _ in 0..n {
                // Inverse CDF of the Plummer sphere: r = a / sqrt(u^{-2/3} - 1)
                let u = st.range(1e-9, 1.0 - 1e-9);
                let r = 1.0 / (u.powf(-2.0 / 3.0) - 1.0).max(1e-12).sqrt();
                out.push(st.direction().scale(r.min(12.0)));
            }
        }
        Profile::Disk { scale_height_ratio } => {
            for _ in 0..n {
                // Exponential in R (inverse CDF via a two-term approximation),
                // sech^2 in z — the observed structure of a stellar disk.
                let u = st.range(1e-9, 1.0 - 1e-9);
                let r = -(1.0 - u).ln();
                let phi = st.range(0.0, std::f64::consts::TAU);
                let uz = st.range(1e-9, 1.0 - 1e-9);
                let z = scale_height_ratio * 0.5 * ((uz / (1.0 - uz)).ln());
                out.push(crate::math::v3(r * phi.cos(), r * phi.sin(), z));
            }
        }
        Profile::Shell => {
            for _ in 0..n {
                let r = 1.0 + 0.05 * st.normal();
                out.push(st.direction().scale(r));
            }
        }
        Profile::WoodsSaxon => {
            // rho(r) = rho0 / (1 + exp((r - R)/a)); rejection-sample it.
            for _ in 0..n {
                let mut r = 1.0;
                for _ in 0..32 {
                    let cand = st.uniform().cbrt() * 1.4;
                    let p = 1.0 / (1.0 + ((cand - 1.0) / 0.11).exp());
                    if st.uniform() < p {
                        r = cand;
                        break;
                    }
                }
                out.push(st.direction().scale(r));
            }
        }
        Profile::Lattice => {
            // Simple cubic with thermal displacement — enough to be a solid.
            let side = (n as f64).cbrt().ceil() as usize;
            let step = 2.0 / side as f64;
            for i in 0..n {
                let ix = i % side;
                let iy = (i / side) % side;
                let iz = i / (side * side);
                let base = crate::math::v3(
                    ix as f64 * step - 1.0,
                    iy as f64 * step - 1.0,
                    iz as f64 * step - 1.0,
                );
                out.push(base + st.normal3().scale(step * 0.05));
            }
        }
    }
    out
}

fn sample_velocities(
    agg: &Aggregate,
    spec: ProlongSpec,
    n: usize,
    pos: &[Vec3],
    seed: u64,
    key: u128,
    epoch: u32,
) -> Vec<Vec3> {
    let mut st = Stream::at(seed, key, epoch, Purpose::Velocities);
    let mut out = Vec::with_capacity(n);
    let turb = spec.turbulent_fraction.clamp(0.0, 1.0);
    let extent = {
        let mut m: f64 = 0.0;
        for p in pos.iter() {
            m = m.max(p.max_abs());
        }
        if m > 0.0 { m } else { 1.0 }
    };
    for i in 0..n {
        // Thermal part: isotropic Gaussian, i.e. Maxwell-Boltzmann.
        let thermal = st.normal3();
        if turb > 0.0 {
            // Turbulent part: a large-scale solenoidal field, giving neighbours
            // correlated velocities. Without this, "gas" looks like a gas of
            // independent points rather than a fluid, and the Larson
            // size-linewidth relation comes out flat.
            //
            // The positions are normalised before entering the trigonometry.
            // Feeding raw metres to `sin` at galactic scales evaluates the
            // sine of ~10^20 radians, where f64 has no fractional bits left and
            // the "field" is pure noise uncorrelated with position — which
            // silently defeats the entire point of the term.
            let p = pos[i].scale(1.0 / extent);
            let k1 = crate::math::v3(1.7, -0.9, 1.1);
            let k2 = crate::math::v3(-0.6, 1.5, 0.8);
            let s1 = (p.dot(k1)).sin();
            let s2 = (p.dot(k2)).cos();
            let swirl = crate::math::v3(s1 * k2.z - s2 * k1.z, s2 * k1.x - s1 * k2.x, s1 * k2.y - s2 * k1.y);
            out.push(thermal.scale(1.0 - turb) + swirl.scale(turb));
        } else {
            out.push(thermal);
        }
    }
    let _ = agg;
    out
}

fn sample_compositions(
    agg: &Aggregate,
    spec: ProlongSpec,
    masses: &[f64],
    seed: u64,
    key: u128,
    epoch: u32,
) -> Vec<Composition> {
    let n = masses.len();
    if spec.composition_scatter <= 0.0 {
        return vec![agg.composition; n];
    }
    let mut st = Stream::at(seed, key, epoch, Purpose::Composition);
    let mut comps: Vec<Composition> = Vec::with_capacity(n);
    for _ in 0..n {
        let mut c = agg.composition.0;
        for v in c.iter_mut() {
            if *v > 0.0 {
                *v *= (1.0 + spec.composition_scatter * st.normal()).max(1e-6);
            }
        }
        comps.push(Composition(c).normalised());
    }
    // Exact rebalancing: the mass-weighted mean must equal the parent's
    // composition to the last bit, because baryon number, lepton number and
    // charge are all derived from it.
    //
    // Rescaling each species and then renormalising each body are two
    // constraints that fight: the renormalisation undoes part of the rescaling.
    // Alternating them converges geometrically (this is Sinkhorn scaling), and
    // eight rounds takes the residual to round-off. Doing it once — the obvious
    // implementation — leaves a 10^-3 error in the composition, which
    // propagates straight into baryon number and shows up as the world quietly
    // gaining nucleons every time a user zooms in.
    let total_mass = det_sum_by(n, &|i| masses[i]);
    for _round in 0..8 {
        for s in 0..NSPECIES {
            let have = det_sum_by(n, &|i| masses[i] * comps[i].0[s]) / total_mass;
            let want = agg.composition.0[s];
            if have > 1e-300 && want > 0.0 {
                let f = want / have;
                for c in comps.iter_mut() {
                    c.0[s] *= f;
                }
            } else if want <= 0.0 {
                for c in comps.iter_mut() {
                    c.0[s] = 0.0;
                }
            }
        }
        for c in comps.iter_mut() {
            *c = c.normalised();
        }
    }
    comps
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

/// Shift positions so the centre of mass is at the origin, returning the shift
/// that was applied — callers with parallel geometry (joint locations, member
/// endpoints) must apply the same one or their geometry silently detaches from
/// the bodies it describes.
fn recentre(pos: &mut [Vec3], masses: &[f64], total: f64) -> Vec3 {
    let n = pos.len();
    let com = det_sum_v3_by(n, &|i| pos[i].scale(masses[i])).scale(1.0 / total);
    for p in pos.iter_mut() {
        *p -= com;
    }
    com
}

/// Choose the scale factor that makes `restrict` report exactly `target`.
/// `restrict` uses `1.291 * rms` (the RMS-to-uniform-sphere conversion), so we
/// invert that same expression rather than guessing.
fn radius_scale(pos: &[Vec3], masses: &[f64], total: f64, target: f64) -> f64 {
    let n = pos.len();
    let r2 = det_sum_by(n, &|i| masses[i] * pos[i].norm2()) / total;
    let rms = r2.max(0.0).sqrt();
    if rms > 0.0 && target > 0.0 {
        target / (1.291 * rms)
    } else {
        1.0
    }
}

fn inertia_tensor_of(pos: &[Vec3], masses: &[f64]) -> crate::math::Mat3 {
    inertia_of_slices(pos, masses)
}

/// Gravitational binding energy of the sampled configuration. Exact for small
/// sets; for large sets we use the analytic form for the profile, which is what
/// the sampled set converges to and avoids an O(n^2) blow-up during zoom.
fn potential_estimate(
    pos: &[Vec3],
    masses: &[f64],
    softening: f64,
    profile: Profile,
    agg: &Aggregate,
) -> f64 {
    let n = pos.len();
    if n <= 512 {
        let bodies: Vec<Body> = pos
            .iter()
            .zip(masses)
            .map(|(p, m)| Body {
                pos: *p,
                mass: *m,
                ..Default::default()
            })
            .collect();
        return mutual_gravitational_energy(&bodies, softening);
    }
    // Analytic coefficients: U = -k G M^2 / R.
    let k = match profile {
        Profile::Uniform => 0.6,
        Profile::Plummer => 0.5 * std::f64::consts::FRAC_PI_2 / 2.0,
        Profile::Disk { .. } => 0.4,
        Profile::Shell => 0.5,
        Profile::WoodsSaxon => 0.72,
        Profile::Lattice => 0.6,
    };
    let r = {
        let total = det_sum_by(n, &|i| masses[i]);
        let r2 = det_sum_by(n, &|i| masses[i] * pos[i].norm2()) / total;
        (r2.max(0.0).sqrt() * 1.291).max(1e-30)
    };
    -k * G * agg.mass * agg.mass / r
}

fn child_radius(agg: &Aggregate, spec: ProlongSpec, mass: f64, n: usize) -> f64 {
    match spec.kind {
        BodyKind::Star => {
            // Main-sequence mass-radius relation, both branches.
            let m = mass / M_SUN;
            let r = if m < 1.0 {
                m.powf(0.8)
            } else {
                m.powf(0.57)
            };
            r * R_SUN
        }
        BodyKind::Nucleus | BodyKind::Nucleon => {
            let a = (mass / AMU).max(1.0);
            1.2e-15 * a.cbrt()
        }
        BodyKind::Atom => BOHR,
        BodyKind::Electron => LAMBDA_COMPTON_E,
        BodyKind::Photon => 0.0,
        _ => agg.radius / (n as f64).cbrt() * 0.5,
    }
}

/// Largest share of a materialisation budget that unstructured matter may take,
/// however much of the mass it represents.
pub const LITTER_BUDGET_SHARE: f64 = 0.15;

/// Geometry handed to [`close_books`]: where the parts are, how heavy they are,
/// and what they are made of. Produced either by the max-entropy samplers above
/// or by a developmental program in `morph.rs`.
pub(crate) struct Parts {
    pub pos: Vec<Vec3>,
    pub masses: Vec<f64>,
    pub comps: Vec<Composition>,
    /// Per-part radii. Empty means "derive them from the spec".
    pub radii: Vec<f64>,
    pub kind: BodyKind,
}

/// Project a set of parts onto the conserved-quantity constraint surface and
/// emit bodies.
///
/// This is the single implementation of the engine's central guarantee, shared
/// by every way of producing geometry. A structured object generated by a
/// growth program goes through exactly the same projection as a sampled gas
/// cloud, which is the reason adding morphology did not require a second
/// conservation story — only a second way of choosing where the parts go.
#[allow(clippy::too_many_arguments)]
pub(crate) fn close_books(
    agg: &Aggregate,
    parts: Parts,
    resid: Vec<Vec3>,
    phi: f64,
    random_ke_target: f64,
    omega: Vec3,
    ke_rot: f64,
    scale: f64,
    report: &mut ProlongReport,
) -> Vec<Body> {
    let Parts { pos, masses, comps, radii, kind } = parts;
    let n = pos.len();
    let spec = ProlongSpec::new(n.max(1), Profile::Uniform, MassSpectrum::Equal, kind);
    let inertia = inertia_of_slices(&pos, &masses);
    let mut resid = resid;

    // (a) remove net momentum
    let p_res = det_sum_v3_by(n, &|i| resid[i].scale(masses[i]));
    let v_mean = p_res.scale(1.0 / agg.mass);
    for v in resid.iter_mut() {
        *v -= v_mean;
    }
    // (b) remove net angular momentum
    let l_res = det_sum_v3_by(n, &|i| pos[i].cross(resid[i].scale(masses[i])));
    if let Some(w) = inertia.solve(l_res) {
        for i in 0..n {
            resid[i] -= w.cross(pos[i]);
        }
        let p2 = det_sum_v3_by(n, &|i| resid[i].scale(masses[i])).scale(1.0 / agg.mass);
        for v in resid.iter_mut() {
            *v -= p2;
        }
    }
    // `resid` now carries zero momentum and zero angular momentum, so it can be
    // scaled freely without disturbing either constraint.

    let ke_res = det_sum_by(n, &|i| 0.5 * masses[i] * resid[i].norm2());
    let mut s = if ke_res > 0.0 {
        ((random_ke_target - ke_rot).max(0.0) / ke_res).sqrt()
    } else {
        0.0
    };
    let mut om = omega;
    let mut vb = agg.momentum.scale(1.0 / agg.mass);

    // Targets, stated in exactly the form `restrict` will measure them.
    let p_target = agg.momentum;
    let l_target = agg.spin;
    let k_target = random_ke_target + crate::state::bulk_kinetic(agg.mass, agg.momentum);

    let mut vel = vec![Vec3::ZERO; n];
    let build = |vel: &mut Vec<Vec3>, s: f64, om: Vec3, vb: Vec3| {
        for i in 0..n {
            vel[i] = resid[i].scale(s) + om.cross(pos[i]) + vb;
        }
    };
    let com_now = det_sum_v3_by(n, &|i| pos[i].scale(masses[i])).scale(1.0 / agg.mass);
    let mut gmass = vec![0.0f64; n];
    for _pass in 0..10 {
        build(&mut vel, s, om, vb);
        // The functionals are evaluated straight from the slices. Materialising
        // a `Vec<Body>` per pass (three per pass, at 184 bytes each) dominated
        // the cost of prolongation — it was more expensive than every physical
        // computation in the sampler put together.
        for i in 0..n {
            gmass[i] = crate::coords::gamma(vel[i]) * masses[i];
        }
        let p_now = det_sum_v3_by(n, &|i| vel[i].scale(gmass[i]));
        let l_now = det_sum_v3_by(n, &|i| (pos[i] - com_now).cross(vel[i].scale(gmass[i])));
        let k_now = det_sum_by(n, &|i| (crate::coords::gamma(vel[i]) - 1.0) * masses[i] * crate::units::C2);

        // Energy: evaluate at s = 0 to isolate the part that does not scale.
        let k_fixed = det_sum_by(n, &|i| {
            let v = om.cross(pos[i]) + vb;
            (crate::coords::gamma(v) - 1.0) * masses[i] * crate::units::C2
        });

        let gm = det_sum_by(n, &|i| gmass[i]);
        if gm > 0.0 {
            vb += (p_target - p_now).scale(1.0 / gm);
        }
        let inertia_rel = inertia_of_slices(&pos, &gmass);
        if let Some(dw) = inertia_rel.solve(l_target - l_now) {
            om += dw;
        }
        let denom = k_now - k_fixed;
        if denom > 0.0 && (k_target - k_fixed) > 0.0 {
            s *= ((k_target - k_fixed) / denom).sqrt();
        }
    }
    build(&mut vel, s, om, vb);

    report.realised_radius = agg.radius * scale;
    report.rotational_fraction = if random_ke_target > 0.0 {
        ke_rot / random_ke_target
    } else {
        0.0
    };
    report.potential = phi;

    // ---- 6. composition, charge, internal energy ------------------------
    let mut bodies: Vec<Body> = Vec::with_capacity(n);
    for i in 0..n {
        let frac = masses[i] / agg.mass;
        bodies.push(Body {
            pos: pos[i] + agg.com,
            vel: vel[i],
            mass: masses[i],
            radius: if radii.is_empty() {
                child_radius(agg, spec, masses[i], n)
            } else {
                radii[i]
            },
            charge: agg.charge * frac,
            internal_energy: 0.0,
            temperature: agg.temperature,
            composition: comps[i],
            spin: Vec3::ZERO,
            slot: i as u32,
            kind: spec.kind,
        });
    }

    // ---- 6.4 close the energy books exactly -----------------------------
    //
    // The velocity scaling above hits the energy target by iteration, and
    // iteration can fall short: an extreme mass spectrum makes the inertia
    // tensor ill-conditioned, and the rotational and residual components stop
    // separating cleanly. Rather than iterate harder and hope, close the
    // remaining gap *algebraically*.
    //
    // `restrict` computes the parent's internal energy as
    // `E_kinetic + sum(U_i) - K_bulk`, so `sum(U_i)` is a free parameter that
    // appears linearly and affects neither momentum nor angular momentum.
    // Solving for it makes total energy exact by construction for any
    // configuration whatsoever, and turns the velocity iteration from a
    // correctness requirement into a quality one: the better it converges, the
    // more of the energy sits where it physically belongs (bulk motion) rather
    // than in the internal account.
    {
        let e_kin = crate::state::kinetic_energy_of(&bodies);
        let target_internal_sum = crate::state::bulk_kinetic(agg.mass, agg.momentum)
            + agg.internal_energy
            + agg.binding_energy
            - phi
            - e_kin;
        report.internal_energy_residual = if agg.internal_energy.abs() > 0.0 {
            target_internal_sum / agg.internal_energy.abs()
        } else {
            0.0
        };
        for b in bodies.iter_mut() {
            b.internal_energy = target_internal_sum * (b.mass / agg.mass);
        }
    }

    // ---- 6.5 residual angular momentum becomes intrinsic spin -----------
    //
    // Some configurations simply cannot carry the requested angular momentum
    // as orbital motion: two point masses have no moment of inertia about their
    // own axis, and a collinear set has none about the line. The inertia solve
    // correctly reports this as singular. The angular momentum still has to go
    // somewhere, and the honest place is the children's intrinsic spin — which
    // is exactly what the parent's spin *was* before we refined it, one level
    // down. Distributing it by mass fraction makes L exact to round-off in
    // every configuration, degenerate or not, and disturbs neither momentum nor
    // energy (intrinsic spin energy stays in the parent's internal budget).
    {
        let com_now = det_sum_v3_by(n, &|i| bodies[i].pos.scale(bodies[i].mass)).scale(1.0 / agg.mass);
        let l_now = crate::state::total_spin(&bodies, com_now);
        let residual = agg.spin - l_now;
        for b in bodies.iter_mut() {
            b.spin += residual.scale(b.mass / agg.mass);
        }
    }

    // ---- 7. verify, do not assume ---------------------------------------
    //
    // Measured with the same potential estimator the engine will use when it
    // coarsens these bodies back down, because "conserved" is only meaningful
    // relative to a fixed definition of the terms.
    let mut back = crate::state::restrict(&bodies, phi);
    // These are properties the children cannot report on — the node's
    // surroundings, the free energy locked in its structure, and what it has
    // already dumped into the environment — so they pass through both
    // directions unchanged.
    back.external_potential = agg.external_potential;
    back.chemical_energy = agg.chemical_energy;
    back.entropy_exported = agg.entropy_exported;
    let scales = crate::state::Scales::of(&bodies);
    report.conservation_error = back.conserved().error_against(&agg.conserved(), &scales);
    report.scales = scales;

    bodies
}

/// Materialise a structure from its developmental state instead of from a
/// distribution.
///
/// The only thing that changes is where the parts come from. Everything after
/// that — the momentum, angular-momentum and energy projection, the algebraic
/// energy close, the intrinsic-spin residual, the verification — is the same
/// [`close_books`] the sampled path uses. A generated oak is held to exactly
/// the conservation standard a gas cloud is, because it goes through exactly
/// the same code.
///
/// Two things differ from `prolong`, and both follow from the same principle:
/// **the developmental state is the authority on the structure's geometry.**
///
/// * The positions are *not* rescaled to hit some independently-stored radius.
///   The morphology decides how big the thing is, and the aggregate's radius is
///   kept equal to `Morphology::extent()` by the growth step, so the two agree
///   by construction rather than by correction.
/// * There is no gravitational relaxation loop. A tree is held together by
///   chemistry, not self-gravity, and expanding it until it is gravitationally
///   comfortable would be nonsense — its size is set by how much it has grown.
pub fn prolong_structured(
    agg: &Aggregate,
    morph: &crate::morph::Morphology,
    budget: usize,
    world_seed: u64,
    path_key: u128,
    epoch: u32,
) -> (Vec<Body>, crate::topology::Topology, ProlongReport) {
    let mut report = ProlongReport {
        count: budget,
        ..Default::default()
    };
    if agg.mass <= 0.0 || !agg.is_finite() {
        return (Vec::new(), crate::topology::Topology::default(), report);
    }

    // ---- 1. geometry from the program -----------------------------------
    //
    // A node holding a structure is not made *only* of that structure: a forest
    // node has air, soil and leaf litter in it too. So the mass splits, and the
    // unstructured remainder is sampled the ordinary way and appended. Assuming
    // the node is nothing but structure is fine right up until a limb is
    // severed, at which point the structure's mass drops while the node's does
    // not, and the missing mass gets silently redistributed into the surviving
    // branches — a tree that grows heavier every time you prune it.
    let structural = morph.built.clamp(0.0, agg.mass);
    let residual = (agg.mass - structural).max(0.0);
    let residual_frac = residual / agg.mass;

    // Budget is allocated by *salience*, not by mass. Unstructured matter is
    // interchangeable — any sample of it is as good as any other — so a handful
    // of parcels represents it as well as thousands would. Structure is the
    // opposite: its specific arrangement is the whole information content, and
    // it is what an observer is looking at. Splitting the budget by mass gave a
    // 60-tonne forest floor 94% of it and left a three-tonne tree with 500
    // members.
    let litter_share = residual_frac.min(LITTER_BUDGET_SHARE);
    let skel = morph.render(((budget as f64) * (1.0 - litter_share)).round().max(1.0) as usize);
    let n_struct = skel.len();
    if n_struct == 0 {
        return (Vec::new(), crate::topology::Topology::default(), report);
    }

    let skel_geom = skel.clone();
    let mut masses = skel.mass;
    let m_sum = det_sum_by(n_struct, &|i| masses[i]);
    let k = if m_sum > 0.0 { structural / m_sum } else { 0.0 };
    for m in masses.iter_mut() {
        *m *= k;
    }

    let mut pos_all = skel.pos;
    let mut radii_all = skel.radius;
    let mut comps: Vec<Composition> = vec![morph.program.substrate(); n_struct];

    // The unstructured remainder: litter, air, rubble. Sampled from the same
    // max-entropy machinery every other node uses.
    if residual > 0.0 && residual_frac > 1e-12 {
        let n_res = ((budget as f64) * litter_share).round().max(1.0) as usize;
        let mut st = Stream::at(world_seed, path_key, epoch, Purpose::Positions);
        let each = residual / n_res as f64;
        for _ in 0..n_res {
            pos_all.push(st.in_ball());
            masses.push(each);
            radii_all.push(0.0);
            comps.push(agg.composition);
        }
    }
    let n = pos_all.len();
    report.count = n;
    report.structural_parts = n_struct;

    // The skeleton is generated in units of the structure's extent. Scale it so
    // that `restrict` reports exactly `agg.radius` — which the growth step has
    // already set to `morph.extent()`, so this is a unit conversion rather than
    // a correction to the shape.
    let mut pos = pos_all;
    let com_shift = recentre(&mut pos, &masses, agg.mass);
    let scale = radius_scale(&pos, &masses, agg.mass, agg.radius);
    for p in pos.iter_mut() {
        *p = p.scale(scale);
    }

    // Member radii are set by the *density*, not by the position scale.
    //
    // Scaling radii geometrically alongside positions looks harmless and is
    // not: the resulting members enclose a volume with no relation to the mass
    // they carry. It produced a 13 m tree with a half-metre trunk radius —
    // three times the whole tree's volume in the trunk alone. Section modulus
    // goes as r^3, so a factor of two in radius is a factor of eight in
    // apparent strength, and every structural conclusion drawn from it would
    // have been wrong.
    let mut radii: Vec<f64> = radii_all.iter().map(|r| r * scale).collect();
    {
        let target_volume = structural / morph.program.density();
        let mut current = 0.0;
        for i in 0..n_struct.min(radii.len()) {
            let len = (skel_geom.tip[i] - skel_geom.base[i]).norm() * scale;
            current += std::f64::consts::PI * radii[i] * radii[i] * len;
        }
        if current > 0.0 && target_volume > 0.0 {
            let f = (target_volume / current).sqrt();
            for r in radii.iter_mut().take(n_struct) {
                *r *= f;
            }
            report.radius_correction = f;
        }
    }

    // ---- 2. energy budget ------------------------------------------------
    let softening = agg.radius / (n as f64).cbrt() * 0.1;
    let phi = potential_estimate(&pos, &masses, softening, Profile::Uniform, agg);
    // Same expression the sampled path uses. Self-gravity is negligible for a
    // structure, so this is essentially the thermal budget, but it is written
    // the same way so that the two paths cannot drift apart.
    let random_ke_target = (agg.internal_energy + agg.binding_energy - phi)
        .max(agg.internal_energy.abs().max(1e-30));

    // ---- 3. thermal jitter, then the shared projection --------------------
    //
    // A structure is not motionless: its parts vibrate at the ambient
    // temperature. Sampling that jitter and then projecting it is what lets a
    // tree carry a temperature, radiate, and respond to being hit.
    let mut st = Stream::at(world_seed, path_key, epoch, Purpose::ThermalNoise);
    let resid: Vec<Vec3> = (0..n).map(|_| st.normal3()).collect();

    let inertia = inertia_of_slices(&pos, &masses);
    let omega = inertia.solve(agg.spin).unwrap_or(Vec3::ZERO);
    let ke_rot = 0.5 * omega.dot(agg.spin);

    let member_radii = radii.clone();
    let parts = Parts {
        pos,
        masses,
        comps,
        radii,
        kind: morph.body_kind(),
    };
    let bodies = close_books(
        agg,
        parts,
        resid,
        phi,
        random_ke_target,
        omega,
        ke_rot,
        scale,
        &mut report,
    );

    // The joints, in the same index space and the same units as the bodies.
    // Regenerated from the program exactly as the geometry is, so cohesion
    // costs no more to store than shape does.
    let mut topo = crate::topology::Topology::from_skeleton(
        &skel_geom,
        morph.material(),
        com_shift,
        scale,
        &member_radii,
        bodies.len(),
    );

    // Analysed, then proportioned, on creation.
    //
    // The generator decides where members go. It has no way to know what any of
    // them will carry, so the radii it produces are a shape scaled as a group
    // to match the structural mass — which leaves a few members at the point of
    // failure and most of the material in members doing nothing. One analysis
    // and a few passes of fully-stressed sizing fixes both, at the mass the
    // structure already had.
    //
    // This runs every time the structure is materialised, and it must: the
    // radii are regenerated from the program like everything else, so the
    // optimisation has to be part of the regeneration or it would be lost the
    // first time anybody looked away.
    let mut bodies = bodies;
    let (cases, passes) = if topo.is_determinate() {
        (design_cases(morph, &bodies, &topo, 3), crate::solvers::structure::DESIGN_PASSES)
    } else {
        // A redundant structure costs a conjugate-gradient solve per case per
        // pass, and the frame budget has to survive the observer looking at a
        // city. Two directions and two passes is what fits.
        (design_cases(morph, &bodies, &topo, 2), 2)
    };
    report.design = crate::solvers::structure::optimise(&mut bodies, &mut topo, &cases, passes);
    (bodies, topo, report)
}

/// The loads a structure is proportioned against when it is created.
///
/// Its own weight; a vertical overload standing in for anything that settles on
/// it — snow, fruit, a floor's contents; and the flow its program says it lives
/// in, from several directions. Sizing against a single case produces a
/// structure that is optimal for that case and brittle in every other, which is
/// how a tree proportioned only for a westerly ends up falling over in an
/// easterly.
fn design_cases(
    morph: &crate::morph::Morphology,
    bodies: &[Body],
    topo: &crate::topology::Topology,
    directions: usize,
) -> Vec<crate::solvers::structure::LoadField> {
    use crate::solvers::structure as st;
    let ambient = 290.0;
    let mut cases = Vec::with_capacity(directions + 1);

    // Vertical overload. Not a weather preset — a body-force field at two and a
    // half gravities, which is what any load that settles on a structure looks
    // like to the members carrying it: snow on a crown, a fruit crop, a floor's
    // contents. Four gravities was tried and is worse: the vertical case then
    // dominates the envelope and the structure is starved of what it needs to
    // stand up in a wind.
    let mut heavy = st::LoadField::new(bodies.len(), ambient);
    heavy.apply(
        &st::Mechanism::BodyAcceleration(st::G_EARTH.scale(2.5)),
        bodies,
        topo,
    );
    cases.push(heavy);

    let (speed, density) = morph.program.design_flow();
    for k in 0..directions.max(1) {
        let angle = std::f64::consts::TAU * k as f64 / directions.max(1) as f64;
        let mut field = st::LoadField::new(bodies.len(), ambient);
        if speed > 0.0 {
            field.apply(
                &st::Mechanism::FlowDrag {
                    velocity: crate::math::v3(speed * angle.cos(), speed * angle.sin(), 0.0),
                    fluid_density: density,
                    drag_coefficient: 1.2,
                },
                bodies,
                topo,
            );
        }
        field.apply(&st::weather::gravity(), bodies, topo);
        cases.push(field);
    }
    cases
}
