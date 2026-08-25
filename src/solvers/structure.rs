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
//! [`Topology::is_determinate`] decides, and the redundant case goes to
//! [`crate::solvers::frame`] — real 3D beam elements with rotational degrees of
//! freedom, solved matrix-free by preconditioned conjugate gradient, with Euler
//! buckling and plastic redistribution on top.
//!
//! The redundant path was a network of axial and transverse springs before. It
//! is not enough: a spring pair carries no moment between its ends, so a
//! member's rotation is invisible to its neighbours, and a fixed-fixed beam
//! comes out sixty-four times too flexible. Anything braced, portalised or
//! continuous is exactly the case where that error matters, which is the same
//! case that needs the redundant solver in the first place.

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
    /// Compressive load as a fraction of the Euler critical load. At or above 1
    /// the member buckles, however far its stress is from rupture.
    pub buckling: f64,
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

    if topo.is_determinate() {
        // Statically determinate: statics alone fixes the internal forces, and
        // one reverse pass over the array computes them exactly.
        return (accumulate(bodies, topo, &external, loads, n), false, 0);
    }

    // Redundant: how the load divides depends on relative stiffness, so it
    // needs a solve. Real beam elements with rotational degrees of freedom —
    // see `solvers::frame` for why a spring network is not good enough.
    match frame_analyse(bodies, topo, &external, loads, n) {
        Some((joints, iters)) => (joints, true, iters),
        // A solve that did not converge is not a result. Falling back to the
        // determinate answer ignores the alternative load paths, which is
        // conservative — it over-predicts the load on the primary one.
        None => (accumulate(bodies, topo, &external, loads, n), true, 0),
    }
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

        // Buckling is a stability failure, not a strength one: a slender member
        // in compression goes at a load a stress check never notices. The
        // determinate path has the axial force in hand, so it costs one
        // comparison to include.
        let member_length = (topo.tip.get(i).copied().unwrap_or(Vec3::ZERO)
            - topo.base.get(i).copied().unwrap_or(Vec3::ZERO))
        .norm();
        let buckling = if axial < 0.0 && member_length > 0.0 && radius > 0.0 {
            let inertia = std::f64::consts::PI * radius.powi(4) / 4.0;
            let critical = std::f64::consts::PI.powi(2) * topo.material.stiffness * inertia
                / (0.85 * member_length).powi(2);
            if critical > 0.0 {
                -axial / critical
            } else {
                f64::INFINITY
            }
        } else {
            0.0
        };

        let t = loads.temperature.get(i).copied().unwrap_or(loads.ambient);
        let integrity = bond.map(|b| b.integrity).unwrap_or(1.0);
        // Tension is the weak direction for brittle materials, and it is what
        // decides whether masonry topples or merely settles.
        let tensile = axial > 0.0;
        let ratio = if tensile { topo.material.tensile_ratio } else { 1.0 };
        let strength = topo.material.rupture * topo.material.strength_at(t) * integrity * ratio;
        let by_stress = if strength > 0.0 {
            stress / strength
        } else if stress > 0.0 {
            f64::INFINITY
        } else {
            0.0
        };
        // Two independent failure modes; whichever arrives first governs.
        let utilisation = by_stress.max(buckling);

        out.push(JointLoad {
            stress,
            buckling,
            utilisation,
            force: force[i],
            moment,
            carried: carried[i],
        });
    }
    out
}

/// Analyse a redundant structure with real beam elements.
///
/// Joints are the members' shared endpoints: a member's base *is* its parent's
/// tip, so the node count is the member count plus one anchor per ground
/// connection. Distributed member loads are lumped half to each end, which is
/// the standard consistent treatment for a uniformly loaded element.
/// The frame a topology describes: joints, members, and what holds it down.
///
/// Shared by the static analysis and by [`crate::solvers::dynamics`], which is
/// the point — a structure that is analysed and a structure that is animated
/// must be the same structure, down to which joints are welded together and
/// which members are pin-jointed.
pub struct BuiltFrame {
    pub frame: crate::solvers::frame::Frame,
    /// Node at each member's far end, or `u32::MAX` if the member has none.
    pub tip_node: Vec<u32>,
    /// Node at each member's supported end.
    pub base_node: Vec<u32>,
    /// Element index of each member, or `usize::MAX` if it has none.
    pub element_of: Vec<usize>,
}

/// Build the frame for a topology.
///
/// Joints are the members' shared endpoints: a member's base *is* its parent's
/// tip. Members that meet at a point without a support relation between them —
/// three bars to a common apex, a strut closing a triangle — weld into one
/// joint too, because giving each its own coincident node would leave the joint
/// free to come apart. Free and fixed nodes weld separately: an anchor that
/// happens to sit where a free joint is must not drag that joint into the
/// ground.
pub fn build_frame(topo: &Topology, n: usize) -> BuiltFrame {
    use crate::solvers::frame::Frame;

    let mut frame = Frame::new(topo.material);
    let mut tip_node = vec![u32::MAX; n];
    let mut base_node = vec![u32::MAX; n];
    let mut element_of = vec![usize::MAX; n];
    let mut weld = Weld::new(topo, n);

    for i in 0..n {
        if topo.bonds[i].radius <= 0.0 || topo.bonds[i].integrity <= 0.0 {
            continue;
        }
        tip_node[i] = weld.node(&mut frame, topo.tip[i], false);
    }
    for i in 0..n {
        if tip_node[i] == u32::MAX {
            continue;
        }
        let p = topo.support[i];
        base_node[i] = if p != NO_SUPPORT && (p as usize) < n && tip_node[p as usize] != u32::MAX {
            tip_node[p as usize]
        } else {
            weld.node(&mut frame, topo.base[i], true)
        };
    }
    for i in 0..n {
        if tip_node[i] == u32::MAX || base_node[i] == tip_node[i] {
            continue;
        }
        element_of[i] = frame.add_beam(base_node[i], tip_node[i], topo.bonds[i].radius);
        frame.elements[element_of[i]].integrity = topo.bonds[i].integrity;
    }
    for t in &topo.ties {
        let (a, b) = (t.a as usize, t.b as usize);
        if a >= n || b >= n || t.integrity <= 0.0 {
            continue;
        }
        if tip_node[a] == u32::MAX || tip_node[b] == u32::MAX || tip_node[a] == tip_node[b] {
            // Either end missing, or the tie's ends welded into one joint and
            // it has nothing left to hold.
            continue;
        }
        let radius = (t.area / std::f64::consts::PI).max(0.0).sqrt();
        let e = frame.add_tie(tip_node[a], tip_node[b], radius);
        frame.elements[e].integrity = t.integrity;
    }

    BuiltFrame { frame, tip_node, base_node, element_of }
}

fn frame_analyse(
    bodies: &[Body],
    topo: &Topology,
    external: &[Vec3],
    loads: &LoadField,
    n: usize,
) -> Option<(Vec<JointLoad>, u32)> {
    use crate::solvers::frame::Dof;

    let BuiltFrame { frame, tip_node, base_node, element_of } = build_frame(topo, n);
    if frame.elements.is_empty() {
        return None;
    }
    let mut load = vec![Dof::default(); frame.nodes.len()];
    for i in 0..n {
        if tip_node[i] == u32::MAX {
            continue;
        }
        let half = external.get(i).copied().unwrap_or(Vec3::ZERO).scale(0.5);
        load[base_node[i] as usize].t += half;
        load[tip_node[i] as usize].t += half;
    }

    let solution = frame.solve(&load);
    if !solution.converged {
        return None;
    }

    // Back to member index space, and re-derive the temperature-dependent
    // strength here so the two paths apply the same failure criteria.
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let e = element_of.get(i).copied().unwrap_or(usize::MAX);
        if e == usize::MAX || e >= solution.forces.len() {
            out.push(JointLoad {
                stress: 0.0,
                buckling: 0.0,
                utilisation: 0.0,
                force: Vec3::ZERO,
                moment: Vec3::ZERO,
                carried: bodies[i].mass,
            });
            continue;
        }
        let f = solution.forces[e];
        let (axis, length) = member_axis(topo, i);

        // Recover the axial force at the *base* section from the element's.
        //
        // A member's own load acts along its length, and lumping half of it to
        // each end is what makes the end moments right. It does not make the
        // axial force right: the base section carries all of that load and the
        // tip section carries none, so the element's single constant value is
        // the average, half of what the critical section actually sees. This
        // adds the missing half back. Without it a hanging twig reports two
        // thirds of its true axial stress, and the two solvers disagree by up
        // to 10% on exactly the members statics gets exactly right.
        let own = external.get(i).copied().unwrap_or(Vec3::ZERO).dot(axis) * 0.5;
        let axial = f.axial + own;

        let area = topo.bonds[i].area().max(1e-30);
        let section = topo.bonds[i].section_modulus().max(1e-30);
        let stress = f.moment / section + axial.abs() / area;

        let t = loads.temperature.get(i).copied().unwrap_or(loads.ambient);
        let integrity = topo.bonds[i].integrity;
        let tensile = axial > 0.0;
        let ratio = if tensile { topo.material.tensile_ratio } else { 1.0 };
        let strength = topo.material.rupture * topo.material.strength_at(t) * integrity * ratio;
        let by_stress = if strength > 0.0 {
            stress / strength
        } else if stress > 0.0 {
            f64::INFINITY
        } else {
            0.0
        };

        // Buckling likewise judges the section under the most compression.
        let buckling = if axial < 0.0 && length > 0.0 && topo.bonds[i].radius > 0.0 {
            let inertia = std::f64::consts::PI * topo.bonds[i].radius.powi(4) / 4.0;
            let critical = std::f64::consts::PI.powi(2) * topo.material.stiffness * inertia
                / (0.85 * length).powi(2);
            if critical > 0.0 { -axial / critical } else { f64::INFINITY }
        } else {
            0.0
        };

        out.push(JointLoad {
            stress,
            buckling,
            utilisation: by_stress.max(buckling),
            force: axis.scale(axial),
            moment: axis.scale(f.torsion),
            carried: bodies[i].mass,
        });
    }
    Some((out, solution.iterations))
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

/// Welds coincident points into shared frame nodes.
///
/// Independent members that meet in space are one joint, not two points that
/// happen to be at the same coordinates. Without this a truss whose bars all
/// reach a common apex would be three disconnected bars, and a solver asked to
/// hold them together with a zero-length tie has to invert a stiffness that is
/// infinite by construction.
///
/// Positions are quantised onto a grid a millionth of the structure's own size,
/// and a candidate is matched against that cell and its 26 neighbours, so points
/// that agree to within round-off weld whether or not their bits agree.
struct Weld {
    cells: std::collections::HashMap<(i64, i64, i64, bool), Vec<u32>>,
    eps: f64,
}

impl Weld {
    fn new(topo: &Topology, n: usize) -> Weld {
        // Scale from the structure's own extent, so the tolerance means the
        // same thing for a bacterium and for a bridge.
        let mut extent: f64 = 0.0;
        for i in 0..n {
            extent = extent.max(topo.tip[i].norm()).max(topo.base[i].norm());
        }
        let eps = if extent > 0.0 { extent * 1.0e-6 } else { 1.0e-9 };
        Weld { cells: std::collections::HashMap::new(), eps }
    }

    #[inline]
    fn cell(&self, p: Vec3, fixed: bool) -> (i64, i64, i64, bool) {
        (
            (p.x / self.eps).round() as i64,
            (p.y / self.eps).round() as i64,
            (p.z / self.eps).round() as i64,
            fixed,
        )
    }

    fn node(&mut self, frame: &mut crate::solvers::frame::Frame, p: Vec3, fixed: bool) -> u32 {
        let (cx, cy, cz, f) = self.cell(p, fixed);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(bucket) = self.cells.get(&(cx + dx, cy + dy, cz + dz, f)) {
                        for &id in bucket {
                            if (frame.nodes[id as usize] - p).norm() <= self.eps {
                                return id;
                            }
                        }
                    }
                }
            }
        }
        let id = frame.add_node(p, fixed);
        self.cells.entry((cx, cy, cz, f)).or_default().push(id);
        id
    }
}

// ---------------------------------------------------------------------------
// Dynamics
// ---------------------------------------------------------------------------

/// A structure with mass, ready to be integrated through time.
///
/// [`analyse`] answers "does this stand up under that load". This answers "what
/// does it *do* when that load arrives", which is a different question with a
/// different answer: a gust that a quasi-static check passes at 60% utilisation
/// can break the same member outright, because a load that arrives suddenly
/// deflects a structure about twice as far as the same load standing still.
pub struct DynamicStructure {
    pub dynamics: crate::solvers::dynamics::Dynamics,
    /// Node at each member's far end.
    pub tip_node: Vec<u32>,
    /// Node at each member's supported end.
    pub base_node: Vec<u32>,
    /// Element index of each member.
    pub element_of: Vec<usize>,
}

/// Give a topology mass and inertia.
///
/// The structural mass comes from the members' own geometry and the material's
/// density, which is what sets the rotational inertia at a joint. Anything a
/// part weighs *beyond* that — foliage, cladding, accreted snow, a floor's
/// contents — is added as point mass at the member's ends: it is carried, and
/// it changes what the structure does, but it does not stiffen anything.
pub fn dynamic_structure(bodies: &[Body], topo: &Topology) -> Option<DynamicStructure> {
    use crate::solvers::dynamics::Dynamics;

    let n = bodies.len().min(topo.support.len());
    let BuiltFrame { frame, tip_node, base_node, element_of } = build_frame(topo, n);
    if frame.elements.is_empty() {
        return None;
    }
    let density = frame.material.density;
    let mut dynamics = Dynamics::new(frame);
    for i in 0..n {
        let (tip, base) = (tip_node[i], base_node[i]);
        if tip == u32::MAX || base == u32::MAX {
            continue;
        }
        let (axis, len) = member_axis(topo, i);
        let _ = axis;
        let structural = density * topo.bonds[i].area() * len;
        let carried = (bodies[i].mass - structural).max(0.0);
        if carried > 0.0 {
            dynamics.add_point_mass(tip, carried * 0.5);
            dynamics.add_point_mass(base, carried * 0.5);
        }
    }
    Some(DynamicStructure { dynamics, tip_node, base_node, element_of })
}

impl DynamicStructure {
    /// Turn a load field into nodal forces.
    ///
    /// A member's load acts along its length, so half of it goes to each end —
    /// the same consistent lumping the static path uses, for the same reason:
    /// so that a structure held at a steady load settles to exactly the
    /// deflection [`analyse`] predicts for it.
    pub fn nodal_loads(&self, field: &LoadField) -> Vec<crate::solvers::frame::Dof> {
        use crate::solvers::frame::Dof;
        let mut load = vec![Dof::default(); self.dynamics.frame.nodes.len()];
        for i in 0..self.tip_node.len().min(field.len()) {
            let (tip, base) = (self.tip_node[i], self.base_node[i]);
            if tip == u32::MAX || base == u32::MAX {
                continue;
            }
            let half = field.force[i].scale(0.5);
            load[base as usize].t += half;
            load[tip as usize].t += half;
        }
        load
    }

    /// Advance by `h` seconds under a load field.
    pub fn advance(
        &mut self,
        field: &LoadField,
        h: f64,
    ) -> crate::solvers::dynamics::StepReport {
        let load = self.nodal_loads(field);
        self.dynamics.step(&load, h)
    }

    /// Move the bodies onto the deformed structure.
    ///
    /// Parts sit at their members' midpoints, so a part's position and velocity
    /// are the mean of its member's two ends. Writing back rather than
    /// integrating the bodies directly is deliberate: the frame's degrees of
    /// freedom are the joints, and a part is a view of the two joints it lies
    /// between, not an independent thing that could disagree with them.
    pub fn write_back(&self, bodies: &mut [Body], h: f64) {
        let deformed = self.dynamics.deformed();
        let vel = &self.dynamics.velocity;
        let _ = h;
        for i in 0..bodies.len().min(self.tip_node.len()) {
            let (tip, base) = (self.tip_node[i], self.base_node[i]);
            if tip == u32::MAX || base == u32::MAX {
                continue;
            }
            bodies[i].pos = (deformed[tip as usize] + deformed[base as usize]).scale(0.5);
            bodies[i].vel = (vel[tip as usize].t + vel[base as usize].t).scale(0.5);
        }
    }

    /// Members that have failed, in member-index space.
    pub fn failed_members(&self, broken: &[usize]) -> Vec<u32> {
        let mut out = Vec::new();
        for (member, &e) in self.element_of.iter().enumerate() {
            if e != usize::MAX && broken.contains(&e) {
                out.push(member as u32);
            }
        }
        out
    }

    /// Deformed geometry of every member, as `(base, tip)` pairs — what a
    /// renderer needs and nothing more.
    pub fn deformed_members(&self) -> Vec<(Vec3, Vec3)> {
        let d = self.dynamics.deformed();
        self.tip_node
            .iter()
            .zip(&self.base_node)
            .map(|(&t, &b)| {
                if t == u32::MAX || b == u32::MAX {
                    (Vec3::ZERO, Vec3::ZERO)
                } else {
                    (d[b as usize], d[t as usize])
                }
            })
            .collect()
    }
}
