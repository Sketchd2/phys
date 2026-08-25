//! Structural loading and failure, for arbitrary structures.
//!
//! # What is general and what is a preset
//!
//! The solver's vocabulary is *mechanisms*: a body-force field, drag in a
//! moving fluid, mass accreting on upward-facing surfaces, energy conducted
//! along a path, a thermal field, a point impulse. Snow, wind, lightning and
//! bushfire are not concepts this module knows about — they are constructors in
//! [`weather`] that produce mechanisms, and a user can write their own without
//! touching anything here. An earlier version had the weather *in* the solver,
//! which meant every new load case was a new enum arm inside the physics.
//!
//! # Two solvers, chosen by the structure
//!
//! If the load path is a forest — every part supported by at most one other,
//! no alternative routes — the internal forces are exactly determined by
//! statics alone. Accumulate force and moment from the leaves inward and every
//! joint's load is known in O(n) with no solve at all.
//!
//! Add a single brace and that stops being true. The structure becomes
//! statically indeterminate: how load divides between the routes depends on
//! their relative stiffness, and no amount of force balance will tell you. So
//! [`Topology::is_determinate`] decides, and the redundant case solves the
//! spring network for displacements by matrix-free conjugate gradient, converts
//! the tie forces into relieving loads, and *then* runs the same exact
//! accumulation. The fast path is not an approximation of the general one; it
//! is the general one with an empty tie list.

use crate::math::Vec3;
use crate::morph::NO_SUPPORT;
use crate::state::Body;
use crate::topology::Topology;

/// Standard gravity, pointing down the z axis.
pub const G_EARTH: Vec3 = Vec3 { x: 0.0, y: 0.0, z: -9.80665 };

/// Everything acting on a structure, per part.
///
/// Mechanisms write into this; the solver reads it. Keeping them separate is
/// what lets an arbitrary number of simultaneous loads compose — a structure
/// can be burning, loaded with ice and in a gale at the same time, and the
/// solver never learns that any of those words exist.
#[derive(Debug, Clone)]
pub struct LoadField {
    /// External force on each part, N.
    pub force: Vec<Vec3>,
    /// Mass accreted onto each part, kg. Carried as load *and* as inertia.
    pub added_mass: Vec<f64>,
    /// Temperature of each part, K.
    pub temperature: Vec<f64>,
    /// Energy deposited into each part this step, J.
    pub deposited: Vec<f64>,
    /// Parts destroyed outright by a mechanism, before any stress is computed.
    pub destroyed: Vec<bool>,
    /// Total energy the mechanisms delivered, J.
    pub energy_delivered: f64,
    /// Ambient temperature the field was built against, K.
    pub ambient: f64,
}

impl LoadField {
    pub fn new(parts: usize, ambient: f64) -> LoadField {
        LoadField {
            force: vec![Vec3::ZERO; parts],
            added_mass: vec![0.0; parts],
            temperature: vec![ambient; parts],
            deposited: vec![0.0; parts],
            destroyed: vec![false; parts],
            energy_delivered: 0.0,
            ambient,
        }
    }

    pub fn len(&self) -> usize {
        self.force.len()
    }
    pub fn is_empty(&self) -> bool {
        self.force.is_empty()
    }

    /// Apply a mechanism into this field. Mechanisms accumulate.
    pub fn apply(&mut self, m: &Mechanism, bodies: &[Body], topo: &Topology) {
        m.apply(bodies, topo, self);
    }
}

/// A physical way of loading a structure.
///
/// Deliberately mechanisms rather than scenarios. Each one is a generic
/// physical process with parameters; the named weather in [`weather`] is built
/// out of these.
#[derive(Debug, Clone, Copy)]
pub enum Mechanism {
    /// A uniform acceleration field. Gravity, but also a manoeuvring vehicle or
    /// a centrifuge.
    BodyAcceleration(Vec3),

    /// Drag from a fluid moving at `velocity` relative to the structure.
    /// Wind, water, ash-laden air — the fluid's density and the drag
    /// coefficient are parameters, so nothing here is specific to air.
    FlowDrag {
        velocity: Vec3,
        fluid_density: f64,
        drag_coefficient: f64,
    },

    /// Mass settling on upward-facing surfaces. Snow, ice, ash, dust, spray.
    ///
    /// `areal_mass` is what would accumulate per square metre of *silhouette*
    /// if it all stuck. `capacity` bounds what actually adheres — below it the
    /// surface holds everything, above it the excess sheds. `footprint` is the
    /// silhouette; pass zero to estimate it from the geometry.
    SurfaceAccretion {
        areal_mass: f64,
        capacity: f64,
        footprint: f64,
        material_density: f64,
    },

    /// Energy conducted through the structure from an entry point to ground,
    /// distributed by each member's resistance. Lightning, fault current, a
    /// beam boring through.
    ConductedEnergy { joules: f64, entry: u32 },

    /// A thermal environment that parts equilibrate towards, with a time
    /// constant set by their own thermal mass. Fire, radiant heating, cryogenic
    /// immersion — the sign of the difference is not assumed anywhere.
    ThermalField {
        temperature: f64,
        /// Only parts below this height above the structure's base are exposed.
        /// Use `f64::INFINITY` to expose everything.
        ceiling: f64,
        duration: f64,
        /// Convective coefficient, W/m^2/K. Higher couples faster.
        coupling: f64,
    },

    /// A localised impulse. Impact, blast, a falling neighbour.
    PointImpulse { at: u32, impulse: Vec3 },
}

impl Mechanism {
    fn apply(&self, bodies: &[Body], topo: &Topology, out: &mut LoadField) {
        let n = bodies.len().min(topo.support.len()).min(out.len());
        match *self {
            Mechanism::BodyAcceleration(a) => {
                for i in 0..n {
                    out.force[i] += a.scale(bodies[i].mass + out.added_mass[i]);
                }
            }

            Mechanism::FlowDrag {
                velocity,
                fluid_density,
                drag_coefficient,
            } => {
                let speed = velocity.norm();
                if speed <= 0.0 {
                    return;
                }
                let d = velocity.scale(1.0 / speed);
                for i in 0..n {
                    let (axis, len) = member_axis(topo, i);
                    if len <= 0.0 {
                        continue;
                    }
                    // Only the part of the member across the flow presents area.
                    let across = (1.0 - axis.dot(d).abs()).max(0.0);
                    let area = 2.0 * topo.bonds[i].radius * len * across;
                    out.force[i] +=
                        d.scale(0.5 * fluid_density * speed * speed * drag_coefficient * area);
                }
            }

            Mechanism::SurfaceAccretion {
                areal_mass,
                capacity,
                footprint,
                material_density,
            } => {
                // Distribute over upward-facing projected area, but bound the
                // total by the structure's silhouette. Summing per-member areas
                // over-counts by the crown's area index — branches shade one
                // another, and what falls between them reaches the ground.
                let mut per_member = vec![0.0f64; n];
                let mut total = 0.0;
                let mut r2max: f64 = 0.0;
                for i in 0..n {
                    let (axis, len) = member_axis(topo, i);
                    if len <= 0.0 {
                        continue;
                    }
                    let upward = (1.0 - axis.z.abs()).max(0.0);
                    let a = 2.0 * topo.bonds[i].radius * len * upward;
                    per_member[i] = a;
                    total += a;
                    r2max = r2max
                        .max(bodies[i].pos.x * bodies[i].pos.x + bodies[i].pos.y * bodies[i].pos.y);
                }
                let silhouette = if footprint > 0.0 {
                    footprint
                } else {
                    std::f64::consts::PI * r2max * 0.45
                };
                // Adhesion saturates: a surface holds only so much before the
                // rest sheds.
                let held = if capacity > 0.0 {
                    capacity * (1.0 - (-areal_mass / capacity).exp())
                } else {
                    areal_mass
                };
                let intercepted = silhouette * held;
                let share = if total > 0.0 { intercepted / total } else { 0.0 };
                let _ = material_density;
                for i in 0..n {
                    out.added_mass[i] += per_member[i] * share;
                }
            }

            Mechanism::ConductedEnergy { joules, entry } => {
                // The channel follows the support chain to ground. Energy goes
                // where the resistance is, and resistance goes as length over
                // area, so a thin member takes far more energy per kilogram
                // than a thick one — which is why a discharge destroys twigs
                // and leaves the trunk standing.
                let mut path = Vec::new();
                let mut cur = entry as usize;
                let mut guard = 0;
                while cur < n && guard <= n && topo.bonds[cur].radius > 0.0 {
                    path.push(cur);
                    let s = topo.support[cur];
                    if s == NO_SUPPORT {
                        break;
                    }
                    cur = s as usize;
                    guard += 1;
                }
                let resistance: f64 = path.iter().map(|&i| member_resistance(topo, i)).sum();
                out.energy_delivered += joules;
                for &i in &path {
                    let share = if resistance > 0.0 {
                        member_resistance(topo, i) / resistance
                    } else {
                        0.0
                    };
                    let deposited = joules * share;
                    out.deposited[i] += deposited;
                    let needed = bodies[i].mass * topo.material.destruction_enthalpy;
                    if needed > 0.0 && deposited >= needed {
                        out.destroyed[i] = true;
                        out.temperature[i] = topo.material.thermal_gone.max(1500.0);
                    } else if needed > 0.0 {
                        let f = deposited / needed;
                        out.temperature[i] += f * (1500.0 - out.ambient);
                    }
                }
            }

            Mechanism::ThermalField {
                temperature,
                ceiling,
                duration,
                coupling,
            } => {
                let ground = structure_base(topo);
                for i in 0..n {
                    if topo.bonds[i].radius <= 0.0 {
                        continue;
                    }
                    if bodies[i].pos.z - ground > ceiling {
                        continue;
                    }
                    // Newtonian heating with a time constant from the member's
                    // own thermal mass: tau = rho c r / h. Thin members reach
                    // the environment and thick ones barely notice, which is
                    // the whole reason a ground fire takes the understory.
                    let r = topo.bonds[i].radius.max(1e-6);
                    // Lumped thermal time constant, tau = rho c r / h. The char
                    // layer and the moisture a live member carries are what
                    // make `h` an *effective* coefficient well below the raw
                    // convective-plus-radiative value.
                    let tau = topo.material.density * topo.material.specific_heat * r
                        / coupling.max(1e-9);
                    let reached = 1.0 - (-duration / tau.max(1e-9)).exp();
                    let before = out.temperature[i];
                    out.temperature[i] = before + reached * (temperature - before);
                    let dq = bodies[i].mass
                        * topo.material.specific_heat
                        * (out.temperature[i] - before);
                    out.deposited[i] += dq;
                    out.energy_delivered += dq.max(0.0);
                    if topo.material.combustible && out.temperature[i] >= topo.material.thermal_gone
                    {
                        out.destroyed[i] = true;
                    }
                }
            }

            Mechanism::PointImpulse { at, impulse } => {
                if (at as usize) < n {
                    out.force[at as usize] += impulse;
                    out.energy_delivered +=
                        impulse.norm2() / (2.0 * bodies[at as usize].mass.max(1e-12));
                }
            }
        }
    }
}

#[inline]
fn member_axis(topo: &Topology, i: usize) -> (Vec3, f64) {
    let d = topo.tip[i] - topo.base[i];
    let len = d.norm();
    if len > 0.0 {
        (d.scale(1.0 / len), len)
    } else {
        (Vec3::ZERO, 0.0)
    }
}

#[inline]
fn member_resistance(topo: &Topology, i: usize) -> f64 {
    let (_, len) = member_axis(topo, i);
    let a = topo.bonds[i].area().max(1e-12);
    topo.material.resistivity * len.max(1e-9) / a
}

fn structure_base(topo: &Topology) -> f64 {
    topo.base
        .iter()
        .zip(topo.bonds.iter())
        .filter(|(_, b)| b.radius > 0.0)
        .map(|(p, _)| p.z)
        .fold(f64::INFINITY, f64::min)
}

/// Per-joint state after a loading analysis.
#[derive(Debug, Clone, Copy)]
pub struct JointLoad {
    /// Peak fibre stress at the joint, Pa.
    pub stress: f64,
    /// Stress over the strength the joint still has. At or above 1 it fails.
    pub utilisation: f64,
    pub force: Vec3,
    pub moment: Vec3,
    /// Mass supported through this joint, kg.
    pub carried: f64,
}

/// Outcome of loading a structure.
#[derive(Debug, Clone, Default)]
pub struct FailureReport {
    pub broken_sites: Vec<u32>,
    pub detached: Vec<u32>,
    pub detached_mass: f64,
    pub peak_utilisation: f64,
    pub peak_at: u32,
    pub consumed_mass: f64,
    pub energy_delivered: f64,
    /// True when the redundant solver ran.
    pub indeterminate: bool,
    /// Conjugate-gradient iterations used. Zero on the exact path.
    pub solver_iterations: u32,
}

/// Compute the internal loads throughout a structure.
///
/// Dispatches on whether the structure is statically determinate. Both paths
/// end in the same exact accumulation; the redundant one first works out how
/// much load the alternative routes take away from the primary path.
pub fn analyse(bodies: &[Body], topo: &Topology, loads: &LoadField) -> Vec<JointLoad> {
    analyse_with(bodies, topo, loads).0
}

/// As [`analyse`], also reporting how the answer was obtained.
pub fn analyse_with(
    bodies: &[Body],
    topo: &Topology,
    loads: &LoadField,
) -> (Vec<JointLoad>, bool, u32) {
    let n = bodies.len().min(topo.support.len()).min(loads.len());
    let mut external = loads.force.clone();
    external.truncate(n);

    let mut indeterminate = false;
    let mut iterations = 0;
    if !topo.is_determinate() {
        indeterminate = true;
        let (tie_forces, iters) = solve_tie_forces(bodies, topo, &external);
        iterations = iters;
        for (t, f) in topo.ties.iter().zip(&tie_forces) {
            if (t.a as usize) < n {
                external[t.a as usize] += *f;
            }
            if (t.b as usize) < n {
                external[t.b as usize] -= *f;
            }
        }
    }

    (accumulate(bodies, topo, &external, loads, n), indeterminate, iterations)
}

/// The exact O(n) pass over a support forest.
fn accumulate(
    bodies: &[Body],
    topo: &Topology,
    external: &[Vec3],
    loads: &LoadField,
    n: usize,
) -> Vec<JointLoad> {
    let mut force = vec![Vec3::ZERO; n];
    let mut torque = vec![Vec3::ZERO; n];
    let mut carried = vec![0.0f64; n];

    for i in 0..n {
        let f = external.get(i).copied().unwrap_or(Vec3::ZERO);
        force[i] = f;
        torque[i] = bodies[i].pos.cross(f);
        carried[i] = bodies[i].mass + loads.added_mass.get(i).copied().unwrap_or(0.0);
    }
    // Parts are emitted parents-first, so one reverse pass visits every child
    // before its parent.
    for i in (0..n).rev() {
        let p = topo.support[i];
        if p != NO_SUPPORT && (p as usize) < n {
            let pi = p as usize;
            let (f, t, c) = (force[i], torque[i], carried[i]);
            force[pi] += f;
            torque[pi] += t;
            carried[pi] += c;
        }
    }

    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let bond = topo.bonds.get(i);
        let (at, radius) = match bond {
            Some(b) => (b.at, b.radius),
            None => (bodies[i].pos, bodies[i].radius),
        };
        let moment = torque[i] - at.cross(force[i]);
        let (axis, _) = member_axis(topo, i);
        let bending = if axis.norm2() > 0.0 {
            (moment - axis.scale(moment.dot(axis))).norm()
        } else {
            moment.norm()
        };

        let section = std::f64::consts::PI * radius.powi(3) / 4.0;
        let area = std::f64::consts::PI * radius * radius;
        let sigma_bend = if section > 0.0 { bending / section } else { 0.0 };
        let axial = force[i].dot(axis);
        let sigma_axial = if area > 0.0 { axial.abs() / area } else { 0.0 };
        let stress = sigma_bend + sigma_axial;

        let t = loads.temperature.get(i).copied().unwrap_or(loads.ambient);
        let integrity = bond.map(|b| b.integrity).unwrap_or(1.0);
        // Tension is the weak direction for brittle materials, and it is what
        // decides whether masonry topples or merely settles.
        let tensile = axial > 0.0;
        let ratio = if tensile { topo.material.tensile_ratio } else { 1.0 };
        let strength = topo.material.rupture * topo.material.strength_at(t) * integrity * ratio;
        let utilisation = if strength > 0.0 {
            stress / strength
        } else if stress > 0.0 {
            f64::INFINITY
        } else {
            0.0
        };

        out.push(JointLoad {
            stress,
            utilisation,
            force: force[i],
            moment,
            carried: carried[i],
        });
    }
    out
}

/// Solve a statically indeterminate structure for the force each tie carries.
///
/// Matrix-free conjugate gradient on the axial spring network. Every
/// connection — support bonds and ties alike — contributes `k = EA/L` along its
/// own axis; anchored parts are held fixed. The result is the stiffness-weighted
/// load distribution that statics alone cannot supply.
///
/// Matrix-free because assembling a sparse `3n x 3n` matrix for a structure that
/// may be rebuilt every frame costs more than the solve does.
pub fn solve_tie_forces(bodies: &[Body], topo: &Topology, external: &[Vec3]) -> (Vec<Vec3>, u32) {
    let n = bodies.len().min(topo.support.len());
    if topo.ties.is_empty() || n == 0 {
        return (Vec::new(), 0);
    }
    let e = topo.material.stiffness;

    // Each connection is an anisotropic spring: stiff along its own axis
    // (`EA/L`) and much softer across it (`3EI/L^3`, the cantilever stiffness).
    // The ratio is `3r^2/4L^2`, which for a slender member is four orders of
    // magnitude — so slender structures behave as pin-jointed trusses without
    // that having to be assumed, and stubby ones do not.
    //
    // The transverse term is not a detail. Without it the system is singular in
    // every direction no member happens to point along, and conjugate gradient
    // wanders off into the null space.
    let stiffness = |area: f64, radius: f64, len: f64| -> (f64, f64) {
        let inertia = std::f64::consts::PI * radius.powi(4) / 4.0;
        let axial = e * area / len.max(1e-9);
        let transverse = 3.0 * e * inertia / len.max(1e-9).powi(3);
        (axial, transverse)
    };

    // (a, b, axial k, transverse k, axis). `b == usize::MAX` means "ground".
    let mut springs: Vec<(usize, usize, f64, f64, Vec3)> = Vec::new();
    const GROUND: usize = usize::MAX;

    for i in 0..n {
        if topo.bonds[i].radius <= 0.0 {
            continue;
        }
        let (axis, len) = member_axis(topo, i);
        if len <= 0.0 {
            continue;
        }
        let (ka, kt) = stiffness(topo.bonds[i].area(), topo.bonds[i].radius, len);
        let p = topo.support[i];
        if p == NO_SUPPORT {
            // Anchored: the member's base is a fixed point in the ground, and
            // the member itself is the spring between that point and the part.
            // Treating an anchored part as simply immovable is what made the
            // truss unable to deflect at all, and therefore carry nothing in
            // its redundant members.
            springs.push((i, GROUND, ka, kt, axis));
        } else if (p as usize) < n {
            springs.push((i, p as usize, ka, kt, axis));
        }
    }
    for t in &topo.ties {
        let (a, b) = (t.a as usize, t.b as usize);
        if a >= n || b >= n || t.integrity <= 0.0 {
            continue;
        }
        let d = bodies[b].pos - bodies[a].pos;
        let len = d.norm();
        if len <= 0.0 {
            continue;
        }
        let radius = (t.area / std::f64::consts::PI).max(0.0).sqrt();
        let (ka, kt) = stiffness(t.area, radius, len);
        springs.push((a, b, ka, kt, d.scale(1.0 / len)));
    }

    // A part touched by no spring has an empty row in the stiffness matrix but
    // a non-zero load — an inconsistent system that conjugate gradient can
    // never satisfy, so it runs to its iteration cap and the whole solve is
    // discarded. Loose litter in the same node as a structure is exactly such a
    // part, so this is the common case, not a corner one.
    let mut participates = vec![false; n];
    for &(a, b, _, _, _) in &springs {
        participates[a] = true;
        if b != GROUND {
            participates[b] = true;
        }
    }

    let apply = |u: &[Vec3], out: &mut Vec<Vec3>| {
        for v in out.iter_mut() {
            *v = Vec3::ZERO;
        }
        for &(a, b, ka, kt, axis) in &springs {
            let rel = if b == GROUND { u[a] } else { u[a] - u[b] };
            let along = axis.scale(rel.dot(axis));
            let across = rel - along;
            let f = along.scale(ka) + across.scale(kt);
            out[a] += f;
            if b != GROUND {
                out[b] -= f;
            }
        }
        for i in 0..n {
            if !participates[i] {
                out[i] = Vec3::ZERO;
            }
        }
    };

    let b_vec: Vec<Vec3> = (0..n)
        .map(|i| {
            if participates[i] {
                external.get(i).copied().unwrap_or(Vec3::ZERO)
            } else {
                Vec3::ZERO
            }
        })
        .collect();

    let mut u = vec![Vec3::ZERO; n];
    let mut r = b_vec.clone();
    let mut p = r.clone();
    let mut rs = dot(&r, &r);
    let tol2 = rs * 1e-16;
    let mut ap = vec![Vec3::ZERO; n];
    let mut iters = 0u32;
    let mut converged = rs <= 0.0;
    if rs > 0.0 {
        let max_iter = (n * 3).clamp(64, 4000);
        for _ in 0..max_iter {
            iters += 1;
            apply(&p, &mut ap);
            let denom = dot(&p, &ap);
            if !(denom.abs() > 1e-300) {
                break;
            }
            let alpha = rs / denom;
            if !alpha.is_finite() {
                break;
            }
            for i in 0..n {
                u[i] += p[i].scale(alpha);
                r[i] -= ap[i].scale(alpha);
            }
            let rs_new = dot(&r, &r);
            if !rs_new.is_finite() {
                // The system was too ill-conditioned to solve; fall back to the
                // determinate answer rather than returning nonsense.
                return (vec![Vec3::ZERO; topo.ties.len()], iters);
            }
            if rs_new <= tol2 {
                converged = true;
                break;
            }
            let beta = rs_new / rs;
            for i in 0..n {
                p[i] = r[i] + p[i].scale(beta);
            }
            rs = rs_new;
        }
    }
    // An unconverged solve is not a solution. Using one anyway produced tie
    // forces off by ten orders of magnitude, which the accumulation then turned
    // into stresses of 10^14 and demolished the structure on the first frame.
    // Falling back to the determinate answer — ties carry nothing — is
    // conservative, stable, and honest about what was and was not computed.
    if !converged || !u.iter().all(|v| v.is_finite()) {
        return (vec![Vec3::ZERO; topo.ties.len()], iters);
    }

    let forces = topo
        .ties
        .iter()
        .map(|t| {
            let (a, b) = (t.a as usize, t.b as usize);
            if a >= n || b >= n || t.integrity <= 0.0 {
                return Vec3::ZERO;
            }
            let d = bodies[b].pos - bodies[a].pos;
            let len = d.norm();
            if len <= 0.0 {
                return Vec3::ZERO;
            }
            let axis = d.scale(1.0 / len);
            let radius = (t.area / std::f64::consts::PI).max(0.0).sqrt();
            let (ka, kt) = stiffness(t.area, radius, len);
            let rel = u[b] - u[a];
            let along = axis.scale(rel.dot(axis));
            let across = rel - along;
            along.scale(ka) + across.scale(kt)
        })
        .collect();
    (forces, iters)
}

fn dot(a: &[Vec3], b: &[Vec3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x.dot(*y)).sum()
}

/// Break every joint over its limit, then find what is no longer held.
pub fn apply_failures(
    bodies: &[Body],
    topo: &mut Topology,
    loads: &[JointLoad],
    field: &LoadField,
) -> FailureReport {
    let n = bodies.len().min(topo.support.len());
    let mut report = FailureReport::default();
    report.energy_delivered = field.energy_delivered;

    for i in 0..n {
        // Mechanisms can destroy a part outright, before any stress is involved.
        if field.destroyed.get(i).copied().unwrap_or(false) {
            if let Some(b) = topo.bonds.get_mut(i) {
                if b.radius > 0.0 {
                    b.integrity = 0.0;
                    if topo.support[i] != NO_SUPPORT {
                        topo.support[i] = NO_SUPPORT;
                        report.broken_sites.push(topo.site[i]);
                    }
                    if topo.material.combustible {
                        report.consumed_mass += bodies[i].mass;
                    }
                }
            }
            continue;
        }
        let load = match loads.get(i) {
            Some(l) => l,
            None => continue,
        };
        if load.utilisation > report.peak_utilisation {
            report.peak_utilisation = load.utilisation;
            report.peak_at = i as u32;
        }
        if load.utilisation >= 1.0 {
            if let Some(b) = topo.bonds.get_mut(i) {
                if b.parent != NO_SUPPORT {
                    b.integrity = 0.0;
                    topo.support[i] = NO_SUPPORT;
                    report.broken_sites.push(topo.site[i]);
                }
            }
        }
    }

    let mut grounded = vec![false; n];
    for i in 0..n {
        let mut cur = i;
        let mut guard = 0;
        loop {
            guard += 1;
            if guard > n {
                break;
            }
            let s = topo.support[cur];
            if s == NO_SUPPORT {
                grounded[i] = !report.broken_sites.contains(&topo.site[cur]);
                break;
            }
            let si = s as usize;
            if si >= n {
                break;
            }
            if grounded[si] {
                grounded[i] = true;
                break;
            }
            cur = si;
        }
    }

    for i in 0..n {
        let structural = topo.bonds.get(i).map(|b| b.radius > 0.0).unwrap_or(false)
            && topo.site[i] != NO_SUPPORT;
        if structural && !grounded[i] {
            report.detached.push(i as u32);
            report.detached_mass += bodies[i].mass;
        }
    }
    report
}

/// Deflection of a cantilever of length `l` carrying a tip load.
pub fn tip_deflection(force: f64, length: f64, radius: f64, material: &crate::topology::Material) -> f64 {
    let i = std::f64::consts::PI * radius.powi(4) / 4.0;
    if i <= 0.0 {
        return 0.0;
    }
    force * length.powi(3) / (3.0 * material.stiffness * i)
}

/// Named load cases, built out of the mechanisms above.
///
/// Everything here is a convenience. The solver has never heard of snow.
pub mod weather {
    use super::*;

    /// Falling snow of a given depth and density.
    ///
    /// Adhesion capacity scales steeply with density: dry powder barely sticks
    /// and blows off, wet snow near freezing bonds to the surface. There is no
    /// single capacity that is right for both — treating all snow alike makes
    /// powder lethal or wet snow harmless.
    pub fn snow(depth: f64, density: f64, footprint: f64) -> Mechanism {
        let capacity_depth = 0.06 * (density / 200.0).powf(1.5);
        Mechanism::SurfaceAccretion {
            areal_mass: depth * density,
            capacity: capacity_depth * density,
            footprint,
            material_density: density,
        }
    }

    /// Rime or freezing rain: thin, dense and it does not shed.
    pub fn ice(thickness: f64, footprint: f64) -> Mechanism {
        Mechanism::SurfaceAccretion {
            areal_mass: thickness * 917.0,
            capacity: 1.0e9,
            footprint,
            material_density: 917.0,
        }
    }

    /// Volcanic ash fall.
    pub fn ash(depth: f64, footprint: f64) -> Mechanism {
        Mechanism::SurfaceAccretion {
            areal_mass: depth * 700.0,
            capacity: 0.05 * 700.0,
            footprint,
            material_density: 700.0,
        }
    }

    /// Wind in ordinary air.
    pub fn wind(speed: f64, direction: Vec3) -> Mechanism {
        Mechanism::FlowDrag {
            velocity: direction.unit().scale(speed),
            fluid_density: 1.225,
            drag_coefficient: 1.2,
        }
    }

    /// A river, a storm surge, a tsunami — the same drag law, 800 times the
    /// fluid density.
    pub fn current(speed: f64, direction: Vec3) -> Mechanism {
        Mechanism::FlowDrag {
            velocity: direction.unit().scale(speed),
            fluid_density: 1000.0,
            drag_coefficient: 1.4,
        }
    }

    /// A lightning strike entering at a member.
    pub fn lightning(joules: f64, entry: u32) -> Mechanism {
        Mechanism::ConductedEnergy { joules, entry }
    }

    /// A fire front of a given flame temperature and residence time, reaching
    /// `flame_height` above the base.
    pub fn fire(temperature: f64, flame_height: f64, duration: f64) -> Mechanism {
        Mechanism::ThermalField {
            temperature,
            ceiling: flame_height,
            duration,
            coupling: 50.0,
        }
    }

    /// Gravity.
    pub fn gravity() -> Mechanism {
        Mechanism::BodyAcceleration(G_EARTH)
    }
}
