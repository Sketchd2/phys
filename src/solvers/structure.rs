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
    /// Member indices that lost their support, as opposed to the program-stable
    /// site names in `broken_sites`. Indices are what [`detach`] needs; sites
    /// are what the event log needs, and they are not the same thing.
    pub broken_members: Vec<u32>,
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
    build_frame_with(topo, n, true)
}

/// As [`build_frame`], with a choice about what an unsupported member means.
///
/// For a structure that is standing, a member with no support is bolted to the
/// ground. For a piece that has come away, the very same member is *free* — it
/// is the break, and there is nothing on the other side of it any more. The
/// distinction is not visible anywhere in the topology, which describes what is
/// connected to what and has no opinion about the planet; it belongs to whoever
/// knows whether this object is standing or falling.
///
/// Getting it wrong is silent and total: a fragment built with anchored roots
/// hangs in the air exactly where it broke, and every test of falling debris
/// reports that nothing ever reached the ground.
pub fn build_frame_with(topo: &Topology, n: usize, anchored: bool) -> BuiltFrame {
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
            weld.node(&mut frame, topo.base[i], anchored)
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
                        report.broken_members.push(i as u32);
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
                    report.broken_members.push(i as u32);
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
    dynamic_structure_with(bodies, topo, true)
}

/// As [`dynamic_structure`], with a choice about whether unsupported members
/// are anchored to the ground or free. See [`build_frame_with`].
pub fn dynamic_structure_with(
    bodies: &[Body],
    topo: &Topology,
    anchored: bool,
) -> Option<DynamicStructure> {
    use crate::solvers::dynamics::Dynamics;

    let n = bodies.len().min(topo.support.len());
    let BuiltFrame { frame, tip_node, base_node, element_of } =
        build_frame_with(topo, n, anchored);
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

// ---------------------------------------------------------------------------
// Design
// ---------------------------------------------------------------------------

/// What an optimisation pass achieved.
#[derive(Debug, Clone, Copy, Default)]
pub struct DesignReport {
    pub passes: u32,
    /// Peak utilisation before and after.
    pub peak_before: f64,
    pub peak_after: f64,
    /// Standard deviation of utilisation across loaded members, before and
    /// after. This is the number the method is actually minimising: a design is
    /// efficient when every member is working equally hard.
    pub spread_before: f64,
    pub spread_after: f64,
    /// Structural volume before and after. These must be equal — the members
    /// are re-proportioned, not fed.
    pub volume_before: f64,
    pub volume_after: f64,
}

impl DesignReport {
    /// Fractional change in structural volume. Anything but zero means the
    /// optimiser bought its improvement with material it invented.
    pub fn volume_error(&self) -> f64 {
        if self.volume_before > 0.0 {
            (self.volume_after - self.volume_before).abs() / self.volume_before
        } else {
            0.0
        }
    }
}

/// Passes of fully-stressed sizing. Five is well past where the spread stops
/// falling for the structures the engine generates.
pub const DESIGN_PASSES: u32 = 5;

/// How far a member may be re-proportioned from what the generator produced.
///
/// The generator's proportions are a prior, not noise: they encode what the
/// program is *for*, and a member carrying nothing in every design case is
/// still holding the geometry together. Without a floor, fully-stressed design
/// happily reduces such a member to a thread, and the structure then fails at
/// the first load its design cases did not include — which is every real load.
pub const DESIGN_BOUNDS: (f64, f64) = (0.55, 2.6);

/// Re-proportion a structure's members for the loads they actually carry.
///
/// # Why a generated structure needs this
///
/// The generator decides where members go and how they branch. It has no way to
/// know what any of them will end up carrying, so the radii it produces are a
/// shape — scaled as a group to match the structural mass, but not related to
/// the loads member by member. The result is a structure where a few members
/// are at the point of failure while most of the material sits in members doing
/// nothing, and both problems have the same cause.
///
/// # Fully stressed design
///
/// Analyse; then give every member the section that would bring it to the same
/// utilisation as every other; then rescale the whole set back to the volume it
/// started with. Repeat. It converges in a handful of passes and it is the
/// oldest structural optimisation there is, because it is what the answer looks
/// like: at the optimum of a mass-constrained design, no member is idle.
///
/// Bending governs a slender member, and bending stress goes as `M/r^3`, so the
/// radius correction is the cube root of the utilisation ratio. Each pass is
/// damped and clamped, because a member that is briefly carrying nothing would
/// otherwise be sized out of existence and never recover.
///
/// # Why this is not a cheat
///
/// The structural mass does not change. `volume_error` is reported so that
/// "does not change" is something a test checks rather than something this
/// comment claims. What improves is *where the material is*, which is the one
/// thing a real tree spends its life adjusting — a trunk that lays down wood
/// where the wind bends it is running this loop, slowly, with the same
/// objective.
pub fn optimise(
    bodies: &mut [Body],
    topo: &mut Topology,
    cases: &[LoadField],
    passes: u32,
) -> DesignReport {
    let n = bodies
        .len()
        .min(topo.support.len())
        .min(cases.iter().map(|f| f.len()).min().unwrap_or(0));
    let mut report = DesignReport::default();
    if n == 0 || passes == 0 || cases.is_empty() {
        return report;
    }
    let original: Vec<f64> = topo.bonds.iter().take(n).map(|b| b.radius).collect();

    let volume = |t: &Topology| -> f64 {
        (0..n)
            .map(|i| {
                let len = (t.tip[i] - t.base[i]).norm();
                std::f64::consts::PI * t.bonds[i].radius * t.bonds[i].radius * len
            })
            .sum::<f64>()
    };
    report.volume_before = volume(topo);
    if report.volume_before <= 0.0 {
        return report;
    }

    // The envelope: for each member, the worst it does in any design case.
    //
    // Sizing against a single case produces a structure that is optimal for
    // that case and brittle in every other. A tree proportioned only for a wind
    // from the west is a tree that falls over in an easterly, and one
    // proportioned only for wind has nothing left for the winter it spends
    // under snow. Real design does this too and calls it load combinations.
    let envelope = |bodies: &[Body], topo: &Topology| -> Vec<f64> {
        let mut worst = vec![0.0f64; n];
        for field in cases {
            let loads = analyse(bodies, topo, field);
            for i in 0..n {
                let u = loads.get(i).map(|l| l.utilisation).unwrap_or(0.0);
                if u.is_finite() {
                    worst[i] = worst[i].max(u);
                }
            }
        }
        worst
    };

    let stats = |worst: &[f64]| -> (f64, f64) {
        let live: Vec<f64> = worst.iter().copied().filter(|u| *u > 0.0).collect();
        if live.is_empty() {
            return (0.0, 0.0);
        }
        let peak = live.iter().cloned().fold(0.0f64, f64::max);
        let mean = live.iter().sum::<f64>() / live.len() as f64;
        let var = live.iter().map(|u| (u - mean) * (u - mean)).sum::<f64>() / live.len() as f64;
        (peak, var.sqrt())
    };

    let mut worst = envelope(bodies, topo);
    let (peak, spread) = stats(&worst);
    report.peak_before = peak;
    report.spread_before = spread;
    report.peak_after = peak;
    report.spread_after = spread;

    for pass in 0..passes {
        // The target every member is aimed at: the mean of what they are all
        // doing now. Aiming at the peak would shrink everything; aiming at zero
        // is not a target at all.
        let live: Vec<f64> = worst.iter().copied().filter(|u| *u > 0.0).collect();
        if live.is_empty() {
            break;
        }
        let target = live.iter().sum::<f64>() / live.len() as f64;
        if target <= 0.0 {
            break;
        }

        for i in 0..n {
            let r = topo.bonds[i].radius;
            if r <= 0.0 {
                continue;
            }
            let u = worst[i];
            // A member carrying nothing measurable keeps most of its section
            // rather than vanishing.
            let ratio = if u > 1e-6 { u / target } else { 0.5 };
            let step = ratio.cbrt().clamp(0.7, 1.5);
            topo.bonds[i].radius =
                (r * step).clamp(original[i] * DESIGN_BOUNDS.0, original[i] * DESIGN_BOUNDS.1);
        }

        // Back to the volume it started with. Radii enter volume squared, so
        // the correction is a square root. The bounds are re-applied after,
        // which is why the volume is checked at the end rather than assumed.
        let after = volume(topo);
        if after <= 0.0 {
            break;
        }
        let renormalise = (report.volume_before / after).sqrt();
        for i in 0..n {
            topo.bonds[i].radius = (topo.bonds[i].radius * renormalise)
                .clamp(original[i] * DESIGN_BOUNDS.0, original[i] * DESIGN_BOUNDS.1);
            bodies[i].radius = topo.bonds[i].radius;
        }

        worst = envelope(bodies, topo);
        let (peak, spread) = stats(&worst);
        report.passes = pass + 1;
        report.peak_after = peak;
        report.spread_after = spread;
    }

    // A final pass at the volume alone, so the structure ends with exactly the
    // material it started with even where the bounds bit.
    let after = volume(topo);
    if after > 0.0 {
        let renormalise = (report.volume_before / after).sqrt();
        for i in 0..n {
            topo.bonds[i].radius *= renormalise;
            bodies[i].radius = topo.bonds[i].radius;
        }
    }
    report.volume_after = volume(topo);
    report
}

// ---------------------------------------------------------------------------
// Fragmentation
// ---------------------------------------------------------------------------

/// A structure separated into what still stands and what does not.
#[derive(Debug, Clone, Default)]
pub struct Detachment {
    /// Members that are still held, in their original indices.
    pub standing: Vec<u32>,
    /// One list of original member indices per detached piece. A piece is a
    /// connected component of the support forest once the broken members have
    /// been cut from their parents.
    pub pieces: Vec<Vec<u32>>,
}

/// Cut a structure at the members that failed, and find what comes away.
///
/// A break is not a member disappearing; it is a member ceasing to be
/// *supported*. Everything hanging off it comes with it, and what comes away is
/// no longer part of the structure it fell from — it is its own object, with its
/// own roots, its own centre of mass and its own reason to be analysed.
///
/// That last part is the whole point. The static analysis walks the support
/// forest from the leaves inward and needs roots to walk towards. A branch
/// still carrying its old support index would be analysed as though the trunk
/// were holding it up, which is exactly what stopped being true. So the piece is
/// re-rooted: the broken member becomes an anchor of the new object, and the
/// forest is recomputed from there.
pub fn detach(topo: &Topology, broken: &[u32]) -> Detachment {
    let n = topo.support.len();
    let mut out = Detachment::default();
    if n == 0 || broken.is_empty() {
        out.standing = (0..n as u32).collect();
        return out;
    }
    let cut: std::collections::HashSet<u32> = broken.iter().copied().collect();

    // Children, so the forest can be walked outwards from a break.
    let mut children: Vec<Vec<u32>> = vec![Vec::new(); n];
    for i in 0..n {
        let p = topo.support[i];
        if p != NO_SUPPORT && (p as usize) < n && p as usize != i {
            children[p as usize].push(i as u32);
        }
    }

    let mut fallen = vec![false; n];
    let mut piece_of = vec![usize::MAX; n];
    for &root in broken {
        let r = root as usize;
        if r >= n || fallen[r] {
            continue;
        }
        let piece = out.pieces.len();
        let mut stack = vec![root];
        let mut members = Vec::new();
        while let Some(m) = stack.pop() {
            let mi = m as usize;
            if fallen[mi] {
                continue;
            }
            fallen[mi] = true;
            piece_of[mi] = piece;
            members.push(m);
            for &c in &children[mi] {
                // A child that broke on its own account starts its own piece —
                // unless it is inside this one already, in which case it simply
                // becomes a free joint within it.
                if !fallen[c as usize] {
                    stack.push(c);
                }
            }
        }
        members.sort_unstable();
        out.pieces.push(members);
    }
    let _ = cut;
    out.standing = (0..n as u32).filter(|&i| !fallen[i as usize]).collect();
    out
}

/// One piece that has come away and is now falling.
///
/// It is an ordinary structure in every respect except that nothing holds it
/// down: the same topology, the same solver, the same failure criteria. What
/// makes it a fragment is only that its roots are free, so the dynamics has no
/// fixed nodes to react against and the whole thing accelerates.
pub struct Fragment {
    pub bodies: Vec<Body>,
    pub topo: Topology,
    pub dynamics: DynamicStructure,
    /// Seconds since it came away.
    pub age: f64,
    /// Whether any part of it is resting on the ground.
    pub grounded: bool,
    /// Original member indices in the structure it fell from, for reporting.
    pub came_from: Vec<u32>,
}

/// Build a standalone structure out of a set of members.
///
/// Indices are remapped into the new object's own space, and any support that
/// pointed outside the set becomes [`NO_SUPPORT`] — those are the new roots.
/// Ties with an end outside the set go, because there is nothing at the other
/// end of them any more.
pub fn extract(
    bodies: &[Body],
    topo: &Topology,
    members: &[u32],
) -> Option<(Vec<Body>, Topology)> {
    let n = topo.support.len();
    if members.is_empty() {
        return None;
    }
    let mut remap = vec![u32::MAX; n];
    for (new, &old) in members.iter().enumerate() {
        if (old as usize) < n {
            remap[old as usize] = new as u32;
        }
    }
    let mut out_bodies = Vec::with_capacity(members.len());
    let mut piece = Topology {
        material: topo.material,
        ..Topology::default()
    };
    for &old in members {
        let o = old as usize;
        if o >= n {
            return None;
        }
        out_bodies.push(bodies.get(o).copied().unwrap_or_default());
        let parent = topo.support[o];
        let mapped = if parent == NO_SUPPORT || (parent as usize) >= n {
            NO_SUPPORT
        } else {
            remap[parent as usize]
        };
        piece.support.push(if mapped == u32::MAX { NO_SUPPORT } else { mapped });
        piece.site.push(topo.site.get(o).copied().unwrap_or(old));
        piece.base.push(topo.base[o]);
        piece.tip.push(topo.tip[o]);
        let mut bond = topo.bonds[o];
        bond.child = remap[o];
        bond.parent = *piece.support.last().unwrap();
        piece.bonds.push(bond);
    }
    for t in &topo.ties {
        let (a, b) = (t.a as usize, t.b as usize);
        if a >= n || b >= n {
            continue;
        }
        if remap[a] == u32::MAX || remap[b] == u32::MAX {
            continue;
        }
        piece.ties.push(crate::topology::Tie {
            a: remap[a],
            b: remap[b],
            area: t.area,
            integrity: t.integrity,
        });
    }
    Some((out_bodies, piece))
}

impl Fragment {
    /// Make a falling piece out of a set of members, with the velocity they
    /// already had.
    pub fn new(bodies: &[Body], topo: &Topology, members: &[u32]) -> Option<Fragment> {
        let (mut piece_bodies, piece_topo) = extract(bodies, topo, members)?;
        // Free roots: nothing holds this up any more, which is why it is here.
        let dynamics = dynamic_structure_with(&piece_bodies, &piece_topo, false)?;
        for b in piece_bodies.iter_mut() {
            b.vel = Vec3::ZERO;
        }
        Some(Fragment {
            bodies: piece_bodies,
            topo: piece_topo,
            dynamics,
            age: 0.0,
            grounded: false,
            came_from: members.to_vec(),
        })
    }

    pub fn mass(&self) -> f64 {
        self.bodies.iter().map(|b| b.mass).sum()
    }

    /// Whether it has stopped moving enough to be litter rather than debris.
    ///
    /// Both conditions matter. A piece still in the air is not finished however
    /// slowly it is drifting, and one on the ground still sliding is not
    /// finished either.
    pub fn at_rest(&self) -> bool {
        let m = self.mass().max(1e-9);
        // Half a metre a second, whatever the piece weighs. An absolute energy
        // would call a twig settled and a trunk still falling at the same speed.
        self.grounded && self.dynamics.dynamics.kinetic_energy() < 0.5 * m * REST_SPEED * REST_SPEED
    }

    /// Lowest point of the piece.
    ///
    /// Skipping the members that have no node. `deformed_members` returns a
    /// pair of zeroes for those, and a minimum taken over them reports every
    /// piece as sitting at the origin however high it actually is.
    pub fn lowest(&self) -> f64 {
        self.dynamics
            .deformed_members()
            .iter()
            .filter(|(a, b)| a != b)
            .map(|(a, b)| a.z.min(b.z))
            .fold(f64::INFINITY, f64::min)
    }
}

/// A fragment touching something on its way down.
#[derive(Debug, Clone, Copy)]
pub struct Contact {
    /// Member of the falling piece.
    pub falling: u32,
    /// Member of the structure it hit, or `NO_SUPPORT` for the ground.
    pub struck: u32,
    /// Where, in world coordinates.
    pub at: Vec3,
    /// Unit normal, pointing from the struck member towards the falling one.
    pub normal: Vec3,
    /// Closing speed along the normal, m/s. Positive means they are still
    /// approaching.
    pub closing: f64,
}

/// Closest approach between two segments, as parameters along each.
fn segment_closest(a0: Vec3, a1: Vec3, b0: Vec3, b1: Vec3) -> (f64, f64) {
    let u = a1 - a0;
    let v = b1 - b0;
    let w = a0 - b0;
    let (a, b, c) = (u.dot(u), u.dot(v), v.dot(v));
    let (d, e) = (u.dot(w), v.dot(w));
    let denom = a * c - b * b;
    // Parallel segments: pick the projection of one endpoint and be done. The
    // pair is degenerate for contact purposes either way.
    if denom.abs() < 1e-18 {
        let t = if c > 0.0 { (e / c).clamp(0.0, 1.0) } else { 0.0 };
        return (0.0, t);
    }
    let s = ((b * e - c * d) / denom).clamp(0.0, 1.0);
    let t = ((a * e - b * d) / denom).clamp(0.0, 1.0);
    (s, t)
}

impl Fragment {
    /// Everything this piece is touching: members of the structure below it,
    /// and the ground.
    ///
    /// Members are capsules — a segment with a radius — so contact is the
    /// closest approach between two segments against the sum of their radii.
    /// The search is bounded by a grid over the struck structure, because a
    /// falling limb near a thousand-member crown must not cost a thousand
    /// tests per step.
    pub fn contacts(
        &self,
        struck_bodies: &[Body],
        struck: &Topology,
        ground: f64,
    ) -> Vec<Contact> {
        let mine = self.dynamics.deformed_members();
        let vel = &self.dynamics.dynamics.velocity;
        let mut out = Vec::new();
        let n = struck.support.len().min(struck_bodies.len());
        // A piece that has just torn free is still occupying the space it
        // occupied a moment ago, touching the parent it broke from and the
        // siblings beside it. Testing those as collisions catches the limb
        // before it has moved at all and pins it where it broke. It has to be
        // allowed to clear its own footprint first — which is what tearing
        // free is, and takes a limb about a tenth of a second.
        let tearing_free = self.age < TEAR_AWAY;

        for (i, &(a0, a1)) in mine.iter().enumerate() {
            if a0 == a1 {
                continue;
            }
            let ra = self.topo.bonds[i].radius;
            let node = self.dynamics.tip_node[i];
            let v = if node == u32::MAX {
                Vec3::ZERO
            } else {
                vel[node as usize].t
            };

            // The ground first: it is what most of a fallen limb ends up on.
            let low = a0.z.min(a1.z);
            if low - ra <= ground {
                if v.z < 0.0 {
                    let at = if a0.z < a1.z { a0 } else { a1 };
                    out.push(Contact {
                        falling: i as u32,
                        struck: NO_SUPPORT,
                        at,
                        normal: Vec3 { x: 0.0, y: 0.0, z: 1.0 },
                        closing: -v.z,
                    });
                }
                continue;
            }

            if tearing_free {
                continue;
            }
            for j in 0..n {
                let rb = struck.bonds[j].radius;
                if rb <= 0.0 || struck.bonds[j].integrity <= 0.0 {
                    continue;
                }
                let (b0, b1) = (struck.base[j], struck.tip[j]);
                // Cheap rejection on the segment midpoints before the real test.
                let reach = ra + rb + (a1 - a0).norm() * 0.5 + (b1 - b0).norm() * 0.5;
                if ((a0 + a1).scale(0.5) - (b0 + b1).scale(0.5)).norm2() > reach * reach {
                    continue;
                }
                let (s, t) = segment_closest(a0, a1, b0, b1);
                let pa = a0 + (a1 - a0).scale(s);
                let pb = b0 + (b1 - b0).scale(t);
                let d = pa - pb;
                let dist = d.norm();
                if dist >= ra + rb || dist <= 0.0 {
                    continue;
                }
                let normal = d.scale(1.0 / dist);
                let closing = -v.dot(normal);
                if closing <= 0.0 {
                    continue;
                }
                out.push(Contact {
                    falling: i as u32,
                    struck: j as u32,
                    at: pb,
                    normal,
                    closing,
                });
            }
        }
        out
    }

    /// Resolve a set of contacts: stop the piece, and hand the reaction to what
    /// it hit.
    ///
    /// Returns the impulses to deliver into the struck structure, as
    /// `(member, impulse)`. The piece keeps a fraction of its approach speed —
    /// wood on wood does not bounce — and the kinetic energy that disappears in
    /// the collision is what the struck member has to absorb, which is why the
    /// impulse is delivered through the ordinary mechanism vocabulary rather
    /// than applied to the geometry directly. What happens next is the same
    /// stress calculation that decides everything else.
    pub fn resolve(&mut self, contacts: &[Contact], restitution: f64) -> Vec<(u32, Vec3)> {
        let mut delivered: Vec<(u32, Vec3)> = Vec::new();
        if contacts.is_empty() {
            return delivered;
        }
        let hit_ground = contacts.iter().any(|c| c.struck == NO_SUPPORT);

        // The whole piece is behind the contact, not just the joint that
        // touched. A limb is stiff over its own length, so what a branch
        // underneath it has to stop is the limb's momentum and not the few
        // kilograms nearest the point of contact. Using the local lumped mass
        // instead under-reads the force by the number of joints in the piece —
        // fifty times, for a limb of any size — and a hundred and fifty
        // kilograms of falling wood then settles onto a twig without marking
        // it.
        let mass = self.mass().max(1e-9);
        let share = 1.0 / contacts.len() as f64;
        let bounce = 1.0 + restitution.clamp(0.0, 1.0);
        let mut change = Vec3::ZERO;
        for c in contacts {
            let j = mass * c.closing * bounce * share;
            let impulse = c.normal.scale(j);
            change += impulse.scale(1.0 / mass);
            if c.struck != NO_SUPPORT {
                delivered.push((c.struck, impulse.scale(-1.0)));
            } else {
                self.grounded = true;
            }
        }
        // Applied as a rigid-body velocity change, because that is what it is:
        // pushing one joint and letting the piece fold around it is a different
        // event, and not the one that happened.
        for v in self.dynamics.dynamics.velocity.iter_mut() {
            v.t += change;
        }

        // A branch hitting the forest floor thuds; it does not ring. The
        // impulse alone acts on one joint and leaves the rest of the piece
        // swinging about it, which is what an elastic body with free roots
        // does and is not what a limb in leaf litter does. Damping the whole
        // piece is the ground absorbing it.
        if hit_ground {
            for v in self.dynamics.dynamics.velocity.iter_mut() {
                v.t = v.t.scale(GROUND_FRICTION);
                v.r = v.r.scale(GROUND_FRICTION);
            }
            for a in self.dynamics.dynamics.acceleration.iter_mut() {
                *a = crate::solvers::frame::Dof::default();
            }
        }
        delivered
    }
}


/// How much motion survives hitting the ground. Litter does not bounce, and it
/// does not slide.
pub const GROUND_FRICTION: f64 = 0.25;

/// Speed below which a grounded piece is litter rather than debris, m/s.
pub const REST_SPEED: f64 = 0.5;


/// How long a piece is left alone after it comes away, seconds.
pub const TEAR_AWAY: f64 = 0.12;
