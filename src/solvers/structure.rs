//! Structural loading and failure.
//!
//! Given a structure's joints and whatever is pushing on it, work out the
//! internal forces, decide what breaks, and account for what falls off.
//!
//! # The method
//!
//! Loads accumulate from the leaves inward. For each part the solver carries
//! two running sums over that part's whole subtree: the total force, and the
//! total moment about the origin. The moment about any particular joint is then
//! `T - r_joint x F`, which is one cross product, and the peak fibre stress is
//! that moment divided by the section modulus.
//!
//! Because the support graph is a tree and the parts are emitted parents-first,
//! a single reverse pass over the array does the whole accumulation. There is
//! no linear system, no iteration and no convergence criterion — the answer is
//! exact for a determinate structure and conservative for a redundant one.
//!
//! # What this buys
//!
//! Damage stops being scripted. Nothing in the engine says "lightning destroys
//! a tree" or "heavy snow breaks branches". Snow adds mass to upward-facing
//! surfaces, wind adds drag to projected area, lightning deposits enthalpy
//! along a conduction path, fire raises temperature and consumes material —
//! and then the *same* stress calculation decides what survives. A branch
//! breaks because the moment at its base exceeded what its cross-section could
//! carry, which is also why real branches break.

use crate::math::Vec3;
use crate::morph::NO_SUPPORT;
use crate::state::Body;
use crate::topology::{Material, Topology};

/// Standard gravity, pointing down the z axis. Structures are generated with
/// z up, so this is the load case that matters.
pub const G_EARTH: Vec3 = Vec3 { x: 0.0, y: 0.0, z: -9.80665 };

/// Depth of snow a canopy retains once interception has saturated, metres.
/// Measured canopy interception tops out at a few kilograms per square metre;
/// at settled-snow density this reproduces that.
pub const SNOW_RETAINED_DEPTH: f64 = 0.06;

/// Something being done to a structure.
#[derive(Debug, Clone, Copy)]
pub enum Insult {
    /// Snow lying on upward-facing surfaces. Depth in metres of *fallen* snow;
    /// fresh is ~100 kg/m^3, settled ~200, wet ~400. `crown_area` is the
    /// silhouette the structure presents from above — pass the same crown
    /// projection the growth model uses for light capture, or zero to have it
    /// estimated from the geometry.
    Snow { depth: f64, density: f64, crown_area: f64 },
    /// Wind blowing along `direction` at `speed` m/s.
    Wind { speed: f64, direction: Vec3 },
    /// A lightning strike entering at a part and conducting to ground.
    Lightning { joules: f64, entry: u32 },
    /// A fire front: parts within reach are heated for `duration` seconds.
    Fire { temperature: f64, duration: f64, height: f64 },
    /// A localised impact.
    Impact { at: u32, impulse: Vec3 },
}

/// Per-joint state after a loading analysis.
#[derive(Debug, Clone, Copy)]
pub struct JointLoad {
    /// Peak fibre stress at the joint, pascals.
    pub stress: f64,
    /// Stress divided by the strength the joint still has. At or above 1 it
    /// fails.
    pub utilisation: f64,
    /// Force carried through the joint, newtons.
    pub force: Vec3,
    /// Bending moment at the joint, newton-metres.
    pub moment: Vec3,
    /// Mass supported by this joint, kilograms.
    pub carried: f64,
}

/// Outcome of loading a structure.
#[derive(Debug, Clone, Default)]
pub struct FailureReport {
    /// Sites, in the program's own naming, whose joints failed.
    pub broken_sites: Vec<u32>,
    /// Body indices that are no longer connected to the ground.
    pub detached: Vec<u32>,
    /// Mass that left the structure, kilograms.
    pub detached_mass: f64,
    /// Highest utilisation anywhere, whether or not anything broke.
    pub peak_utilisation: f64,
    /// Which body was working hardest.
    pub peak_at: u32,
    /// Mass consumed outright — burned, vaporised — as opposed to detached.
    pub consumed_mass: f64,
    /// Energy the insult delivered into the structure, joules.
    pub energy_delivered: f64,
}

/// Compute the internal loads throughout a structure.
///
/// `extra` is an optional per-body external force, used for snow, wind and
/// impacts. Returns one entry per body, indexed identically.
pub fn analyse(
    bodies: &[Body],
    topo: &Topology,
    gravity: Vec3,
    extra: Option<&[Vec3]>,
    temperature: &[f64],
) -> Vec<JointLoad> {
    let n = bodies.len().min(topo.support.len());
    let mut force = vec![Vec3::ZERO; n];
    let mut torque = vec![Vec3::ZERO; n];
    let mut carried = vec![0.0f64; n];

    // Own contribution.
    for i in 0..n {
        let mut f = gravity.scale(bodies[i].mass);
        if let Some(e) = extra {
            if i < e.len() {
                f += e[i];
            }
        }
        force[i] = f;
        torque[i] = bodies[i].pos.cross(f);
        carried[i] = bodies[i].mass;
    }

    // Accumulate leaves-inward. Parts are emitted parents-first, so iterating
    // backwards visits every child before its parent — the whole tree in one
    // pass, no sorting and no recursion.
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
        // Moment about this joint, from the subtree's total force and moment.
        let moment = torque[i] - at.cross(force[i]);
        // Only the component perpendicular to the member bends it; the parallel
        // component is torsion, which timber is far better at resisting.
        let axis = (topo.tip.get(i).copied().unwrap_or(Vec3::ZERO)
            - topo.base.get(i).copied().unwrap_or(Vec3::ZERO))
        .unit();
        let bending = if axis.norm2() > 0.0 {
            (moment - axis.scale(moment.dot(axis))).norm()
        } else {
            moment.norm()
        };

        let section = std::f64::consts::PI * radius.powi(3) / 4.0;
        let area = std::f64::consts::PI * radius * radius;
        let sigma_bend = if section > 0.0 { bending / section } else { 0.0 };
        let sigma_axial = if area > 0.0 {
            force[i].dot(axis).abs() / area
        } else {
            0.0
        };
        let stress = sigma_bend + sigma_axial;

        let t = temperature.get(i).copied().unwrap_or(290.0);
        let integrity = bond.map(|b| b.integrity).unwrap_or(1.0);
        let strength = topo.material.strength() * topo.material.strength_at(t) * integrity;
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

/// Break every joint that is over its limit, then find what is no longer held.
///
/// Failure cascades: a branch that sheds a limb is *relieved*, but a trunk that
/// loses its support takes everything above it. Both fall out of running the
/// connectivity check after the breaks rather than before.
pub fn apply_failures(
    bodies: &[Body],
    topo: &mut Topology,
    loads: &[JointLoad],
) -> FailureReport {
    let n = bodies.len().min(topo.support.len());
    let mut report = FailureReport::default();

    for (i, load) in loads.iter().enumerate().take(n) {
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

    // Anything whose chain of supports no longer reaches an anchored part is
    // detached. A part is anchored if it was generated with no support at all
    // and its own joint survived.
    let mut grounded = vec![false; n];
    for i in 0..n {
        // Walk up. Parents precede children, so this terminates.
        let mut cur = i;
        let mut guard = 0;
        loop {
            guard += 1;
            if guard > n {
                break;
            }
            let s = topo.support[cur];
            if s == NO_SUPPORT {
                // Root of its own chain: grounded only if it never broke away,
                // which for an originally-anchored part means its site is still
                // absent from the break list.
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
        // Loose litter was never part of the structure and does not "detach".
        let structural = topo.bonds.get(i).map(|b| b.radius > 0.0).unwrap_or(false)
            && topo.site[i] != NO_SUPPORT;
        if structural && !grounded[i] {
            report.detached.push(i as u32);
            report.detached_mass += bodies[i].mass;
        }
    }
    report
}

/// Turn an insult into per-body forces, temperatures and outright destruction.
///
/// Returns the external force field and the temperature field to feed to
/// [`analyse`], plus whatever the insult destroyed directly.
pub fn apply_insult(
    bodies: &[Body],
    topo: &mut Topology,
    insult: Insult,
    ambient: f64,
) -> (Vec<Vec3>, Vec<f64>, FailureReport) {
    let n = bodies.len().min(topo.support.len());
    let mut extra = vec![Vec3::ZERO; n];
    let mut temperature = vec![ambient; n];
    let mut report = FailureReport::default();

    match insult {
        Insult::Snow { depth, density, crown_area } => {
            // Snow settles on the upward-facing projected area of each member —
            // a near-vertical trunk collects almost nothing, a near-horizontal
            // limb collects along its whole length, which is exactly why limbs
            // and not trunks come down under snow.
            //
            // But the total is bounded by the crown's *silhouette*, not by the
            // sum of the members' areas. Branches shade one another and snow
            // falling between them reaches the ground. Summing per-member areas
            // over-counts by the crown's area index — a factor of three to ten
            // for a mature tree — and made 100 mm of snow destroy a tree that
            // in reality would not notice it.
            let mut per_member = vec![0.0f64; n];
            let mut total_projected = 0.0;
            let mut footprint_r2: f64 = 0.0;
            for i in 0..n {
                let axis = (topo.tip[i] - topo.base[i]).unit();
                let len = (topo.tip[i] - topo.base[i]).norm();
                let horizontal = (1.0 - axis.z.abs()).max(0.0);
                let a = 2.0 * topo.bonds[i].radius * len * horizontal;
                per_member[i] = a;
                total_projected += a;
                let r2 = bodies[i].pos.x * bodies[i].pos.x + bodies[i].pos.y * bodies[i].pos.y;
                footprint_r2 = footprint_r2.max(r2);
            }
            // The crown's silhouette, not the furthest twig's reach. Using the
            // maximum horizontal extent treats the crown as a full disc of that
            // radius and roughly doubles the load.
            let footprint = if crown_area > 0.0 {
                crown_area
            } else {
                std::f64::consts::PI * footprint_r2 * 0.45
            };
            // Interception saturates: a canopy holds only so much before the
            // rest slides off or falls through, which is why a metre of snow
            // is not twenty times as damaging as fifty millimetres.
            // Capacity depends strongly on how wet the snow is. Dry powder
            // barely adheres and blows off; wet snow near freezing bonds to
            // itself and to the bark, and that is the snow that brings limbs
            // down. Treating all snow alike either makes powder lethal or makes
            // wet snow harmless — there is no single value that is right.
            let capacity = SNOW_RETAINED_DEPTH * (density / 200.0).powf(1.5);
            let retained = capacity * (1.0 - (-depth / capacity.max(1e-6)).exp());
            let intercepted = footprint * retained * density;
            let share = if total_projected > 0.0 {
                intercepted / total_projected
            } else {
                0.0
            };
            for i in 0..n {
                extra[i] += G_EARTH.scale(per_member[i] * share);
            }
            report.energy_delivered = 0.0;
        }
        Insult::Wind { speed, direction } => {
            let d = direction.unit();
            const AIR_DENSITY: f64 = 1.225;
            const DRAG_COEFFICIENT: f64 = 1.2;
            for i in 0..n {
                let axis = (topo.tip[i] - topo.base[i]).unit();
                let len = (topo.tip[i] - topo.base[i]).norm();
                // Projected area is the part of the member across the flow.
                let across = (1.0 - axis.dot(d).abs()).max(0.0);
                let area = 2.0 * topo.bonds[i].radius * len * across;
                let f = 0.5 * AIR_DENSITY * speed * speed * DRAG_COEFFICIENT * area;
                extra[i] += d.scale(f);
            }
        }
        Insult::Lightning { joules, entry } => {
            // The channel conducts from the entry point to ground along the
            // support chain — the same path the loads travel, because both
            // follow the structure. Energy is deposited in proportion to each
            // member's resistance, which for a uniform resistivity goes as
            // length over area; a thin twig therefore takes far more energy per
            // kilogram than the trunk does, and is destroyed outright.
            // Only a structural member can carry the channel; a strike aimed at
            // loose litter has nothing to conduct through.
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
            let resistance: f64 = path
                .iter()
                .map(|&i| {
                    let len = (topo.tip[i] - topo.base[i]).norm().max(1e-9);
                    let a = topo.bonds[i].area().max(1e-12);
                    len / a
                })
                .sum();
            report.energy_delivered = joules;
            for &i in &path {
                let len = (topo.tip[i] - topo.base[i]).norm().max(1e-9);
                let a = topo.bonds[i].area().max(1e-12);
                let share = if resistance > 0.0 {
                    (len / a) / resistance
                } else {
                    0.0
                };
                let deposited = joules * share;
                // Enough enthalpy to boil the sap and blow the member apart?
                let needed = bodies[i].mass * topo.material.destruction_enthalpy();
                if deposited >= needed && needed > 0.0 {
                    topo.bonds[i].integrity = 0.0;
                    topo.support[i] = NO_SUPPORT;
                    report.broken_sites.push(topo.site[i]);
                    report.consumed_mass += 0.0; // the wood is shattered, not gone
                    temperature[i] = 1500.0;
                } else if needed > 0.0 {
                    // Not destroyed, but heated and weakened.
                    let fraction = deposited / needed;
                    temperature[i] = ambient + fraction * (1500.0 - ambient);
                    topo.bonds[i].integrity *= (1.0 - fraction).clamp(0.0, 1.0);
                }
            }
        }
        Insult::Fire { temperature: front, duration, height } => {
            // Everything below the flame height is heated towards the front
            // temperature; the approach is exponential with a timescale set by
            // the member's thermal mass, so thin twigs reach it and a thick
            // trunk barely warms. That is why a ground fire kills the understory
            // and scorches but does not fell the big trees.
            let ground = topo
                .base
                .iter()
                .zip(topo.bonds.iter())
                .filter(|(_, b)| b.radius > 0.0)
                .map(|(p, _)| p.z)
                .fold(f64::INFINITY, f64::min);
            for i in 0..n {
                // Only the structure burns. The node's unstructured remainder
                // is soil, air and water, and letting the fire consume it made
                // a three-tonne tree release fifty-eight tonnes of smoke.
                if topo.bonds[i].radius <= 0.0 {
                    continue;
                }
                let z = bodies[i].pos.z - ground;
                if z > height {
                    continue;
                }
                let r = topo.bonds[i].radius.max(1e-6);
                // Time constant ~ rho c r / h, with h a convective coefficient.
                let tau = 600.0 * r / 0.05;
                let reached = 1.0 - (-duration / tau.max(1e-6)).exp();
                temperature[i] = ambient + reached * (front - ambient);
                report.energy_delivered +=
                    bodies[i].mass * 1700.0 * (temperature[i] - ambient).max(0.0);
                // Above pyrolysis the material is consumed, not merely weakened.
                let (_, gone) = topo.material.thermal_limits();
                if temperature[i] >= gone {
                    topo.bonds[i].integrity = 0.0;
                    if topo.support[i] != NO_SUPPORT {
                        topo.support[i] = NO_SUPPORT;
                        report.broken_sites.push(topo.site[i]);
                    }
                    report.consumed_mass += bodies[i].mass;
                } else {
                    topo.bonds[i].integrity *=
                        topo.material.strength_at(temperature[i]).max(0.0);
                }
            }
        }
        Insult::Impact { at, impulse } => {
            if (at as usize) < n {
                extra[at as usize] += impulse;
                report.energy_delivered = impulse.norm2() / (2.0 * bodies[at as usize].mass.max(1e-12));
            }
        }
    }
    (extra, temperature, report)
}

/// Deflection of a cantilever of length `l` carrying a tip load `f`.
/// Not used for failure — it is what makes a loaded branch visibly droop.
pub fn tip_deflection(force: f64, length: f64, radius: f64, material: Material) -> f64 {
    let i = std::f64::consts::PI * radius.powi(4) / 4.0;
    if i <= 0.0 {
        return 0.0;
    }
    force * length.powi(3) / (3.0 * material.stiffness() * i)
}
