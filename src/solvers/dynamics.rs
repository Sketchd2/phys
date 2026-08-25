//! Structural dynamics: what the bonds make things *do*.
//!
//! # The gap this closes
//!
//! [`crate::topology`] knew what a bond could take before it failed, and
//! [`crate::solvers::frame`] knew how load divided between bonds. Neither knew
//! what a bond makes its neighbours *do*. The consequence was a structure that
//! could be analysed and could break, but could not move: a tree in a gale
//! either stood exactly still or snapped, with nothing in between. Swaying,
//! ringing down after a gust, whipping back when a limb goes — all of that is
//! the mass and the stiffness arguing with each other, and there was no mass.
//!
//! # Why implicit, and why on the same operator
//!
//! A slender member is stiff along its axis and soft across it; the ratio of
//! those two frequencies is `L/r`, so an explicit integrator has to resolve the
//! axial mode to simulate the bending one. For a 10 m trunk that is a timestep
//! near a microsecond — twenty thousand steps per frame to animate a sway with
//! a period of seconds. Explicit integration is not slow here, it is
//! *unusable*.
//!
//! Newmark-beta removes the stability limit entirely. Writing the step for the
//! displacement increment `d = x_{n+1} - x_n`:
//!
//! ```text
//!     [ M/(B h^2) + g C/(B h) + K ] d
//!         = f - K x_n + M[ v/(B h) + (1/2B - 1) a ] - C[ (1 - g/B) v + h(1 - g/2B) a ]
//! ```
//!
//! With Rayleigh damping `C = a M + b K` the left-hand side is a multiple of
//! the mass plus a multiple of the *same* stiffness operator the static
//! analysis uses. So [`Frame`] grew two coefficients and a lumped mass vector,
//! and everything else is shared: the same conjugate gradient, the same Jacobi
//! preconditioner, the same element forces, the same buckling and rupture
//! criteria. A structure cannot move according to one stiffness and break
//! according to another, because there is only one.
//!
//! # Why not backward Euler
//!
//! Because it does not oscillate for long enough to watch. Backward Euler is
//! the `g = 1, B = 1` corner of the same family and is unconditionally stable,
//! but its amplitude decays by `1/sqrt(1 + (w h)^2)` every step for reasons
//! that have nothing to do with the material: at sixty steps per cycle it
//! removes 93% of a structure's energy in four cycles. A tree under that
//! integrator deflects into the wind and stops dead. Trapezoidal Newmark
//! (`g = 1/2, B = 1/4`) is the member of the family that conserves energy
//! exactly for a linear system, so what damping there is, is the damping that
//! was asked for. [`Dynamics::numerical_damping`] can put dissipation back when
//! it is wanted — after a fracture, when the released energy arrives as
//! high-frequency ringing the timestep cannot resolve — and it is off by
//! default so that it is never doing so silently.
//!
//! # What this is not
//!
//! Small-displacement linear dynamics about the reference geometry. The
//! restoring force is exact for deflections small against member length, which
//! covers a building in a storm and a trunk in a gale; a sapling bent double is
//! outside it, and [`Dynamics::displacement_ratio`] reports how close the
//! current state is to leaving the regime rather than leaving the user to guess.

use crate::math::Vec3;
use crate::solvers::frame::{Dof, ElementForces, Frame};

/// A structure with mass, integrated through time.
#[derive(Debug, Clone)]
pub struct Dynamics {
    /// Reference geometry, elements, supports and material.
    pub frame: Frame,
    /// Displacement of each node from the reference geometry.
    pub displacement: Vec<Dof>,
    /// Velocity of each node, translational and angular.
    pub velocity: Vec<Dof>,
    /// Mass-proportional Rayleigh damping, 1/s. Drag on the structure as a
    /// whole; damps the low modes.
    pub mass_damping: f64,
    /// Stiffness-proportional Rayleigh damping, s. Internal friction in the
    /// material; damps the high modes, which is where it belongs.
    pub stiff_damping: f64,
    /// Numerical dissipation, 0..1. Zero is the trapezoidal rule and conserves
    /// energy exactly for a linear system; larger values damp the modes the
    /// timestep cannot resolve, at the cost of first-order accuracy.
    pub numerical_damping: f64,
    /// Acceleration carried between steps. Newmark needs it; it is state, not a
    /// derived quantity, and recomputing it from the forces each step would
    /// quietly turn the scheme into something else.
    pub acceleration: Vec<Dof>,
    /// Per-element stiffness, reduced where a member has yielded and zero where
    /// one has failed.
    pub stiffness: Vec<f64>,
}

/// What one dynamic step did.
#[derive(Debug, Clone, Default)]
pub struct StepReport {
    pub iterations: u32,
    pub converged: bool,
    /// Members that failed this step, by element index.
    pub broken: Vec<usize>,
    /// Kinetic energy after the step, J.
    pub kinetic: f64,
    /// Elastic strain energy stored after the step, J.
    pub strain: f64,
    /// Energy removed by damping and by the integrator's own dissipation, J.
    /// Positive means energy left the structure.
    pub dissipated: f64,
    /// Strain energy that was stored in members that failed this step, J. It
    /// is not dissipation: it went into breaking the section and into throwing
    /// the pieces, and booking it as damping would hide a fracture inside a
    /// friction term.
    pub released: f64,
    /// Largest nodal displacement as a fraction of the shortest member it is
    /// attached to. Above about 0.1 the linear restoring force is no longer
    /// trustworthy.
    pub displacement_ratio: f64,
}

impl Dynamics {
    /// Attach mass to a frame. Lumped mass and rotational inertia come from the
    /// members' own geometry and the material's density, so nothing here is a
    /// tuning parameter.
    pub fn new(frame: Frame) -> Dynamics {
        let n = frame.nodes.len();
        let stiffness = vec![frame.material.stiffness; frame.elements.len()];
        let mut d = Dynamics {
            frame,
            displacement: vec![Dof::default(); n],
            velocity: vec![Dof::default(); n],
            // Defaults chosen to give a wooden structure a few cycles of
            // visible sway before it settles, which is what a real one does.
            mass_damping: 0.4,
            stiff_damping: 4.0e-4,
            numerical_damping: 0.0,
            acceleration: vec![Dof::default(); n],
            stiffness,
        };
        d.relump();
        d
    }

    /// Recompute lumped mass and rotational inertia from the current members.
    ///
    /// Half of each member's mass goes to each of its ends — the standard
    /// lumping, and exact for the translational modes. The rotational entries
    /// use `m L^2 / 24` about the transverse axes, which is a half-member's
    /// moment about its own end, and `m r^2 / 4` about the member's own axis.
    /// A node with no rotational inertia would make the mass matrix singular
    /// in the rotational block, and the solve would be free to spin joints at
    /// unbounded rate for nothing.
    pub fn relump(&mut self) {
        let n = self.frame.nodes.len();
        self.frame.lumped = vec![Dof::default(); n];
        let rho = self.frame.material.density;
        for e in &self.frame.elements {
            let (a, b) = (e.a as usize, e.b as usize);
            let axis = self.frame.nodes[b] - self.frame.nodes[a];
            let l = axis.norm();
            if l <= 0.0 {
                continue;
            }
            let half = 0.5 * rho * e.area() * l;
            let (e1, e2, e3) = orthonormal(axis);
            let transverse = half * l * l / 24.0;
            let axial = half * e.radius * e.radius / 4.0;
            let sq = |v: Vec3| Vec3 { x: v.x * v.x, y: v.y * v.y, z: v.z * v.z };
            let rot = sq(e1).scale(axial) + (sq(e2) + sq(e3)).scale(transverse);
            for node in [a, b] {
                self.frame.lumped[node].t += Vec3 { x: half, y: half, z: half };
                self.frame.lumped[node].r += rot;
            }
        }
    }

    /// Add mass that is carried but does not itself resist — foliage, snow, a
    /// load on a floor. It changes what the structure does without changing
    /// what it is.
    pub fn add_point_mass(&mut self, node: u32, mass: f64) {
        if let Some(l) = self.frame.lumped.get_mut(node as usize) {
            l.t += Vec3 { x: mass, y: mass, z: mass };
        }
    }

    /// Advance by `h` seconds under nodal loads `f`.
    ///
    /// `f` is the *total* external load — gravity, wind, an impulse divided by
    /// `h`. Members whose internal force passes their limit are cut, and the
    /// strain energy they were holding is reported as released rather than
    /// quietly deleted.
    pub fn step(&mut self, f: &[Dof], h: f64) -> StepReport {
        let n = self.frame.nodes.len();
        let mut report = StepReport::default();
        if n == 0 || h <= 0.0 || self.frame.elements.is_empty() {
            return report;
        }
        let energy_before = self.kinetic_energy() + self.strain_energy();

        // Newmark coefficients. `xi = 0` is trapezoidal: gamma = 1/2, beta =
        // 1/4, and no numerical dissipation at all.
        let xi = self.numerical_damping.clamp(0.0, 1.0);
        let gamma = 0.5 + xi;
        let beta = 0.25 * (gamma + 0.5) * (gamma + 0.5);
        let (a_m, a_c) = (self.mass_damping, self.stiff_damping);

        // Operator: M/(B h^2) + gamma C/(B h) + K, with C = a_m M + a_c K.
        self.frame.mass_scale = 1.0 / (beta * h * h) + gamma * a_m / (beta * h);
        self.frame.stiff_scale = 1.0 + gamma * a_c / (beta * h);

        // Right-hand side. The stiffness part is -K x_n - C_k v_n', which is a
        // single application of the unshifted operator to a combined vector
        // rather than two element loops.
        let c1 = 1.0 - gamma / beta;
        let c2 = h * (1.0 - gamma / (2.0 * beta));
        let combined: Vec<Dof> = (0..n)
            .map(|i| {
                self.displacement[i]
                    .add(self.velocity[i].scale(a_c * c1))
                    .add(self.acceleration[i].scale(a_c * c2))
            })
            .collect();
        let (ms, ks) = (self.frame.mass_scale, self.frame.stiff_scale);
        self.frame.mass_scale = 0.0;
        self.frame.stiff_scale = 1.0;
        let mut k_term = vec![Dof::default(); n];
        self.frame.apply_operator(&combined, &self.stiffness, &mut k_term);
        self.frame.mass_scale = ms;
        self.frame.stiff_scale = ks;

        let m_v = 1.0 / (beta * h) - a_m * c1;
        let m_a = 1.0 / (2.0 * beta) - 1.0 - a_m * c2;
        let rhs: Vec<Dof> = (0..n)
            .map(|i| {
                let m = self.frame.lumped[i];
                let inertial = m
                    .mul(self.velocity[i].scale(m_v).add(self.acceleration[i].scale(m_a)));
                f.get(i)
                    .copied()
                    .unwrap_or_default()
                    .add(inertial)
                    .sub(k_term[i])
            })
            .collect();

        let (delta, iters, converged) =
            self.frame.solve_elastic_opt(&rhs, &self.stiffness, true);
        report.iterations = iters;
        report.converged = converged;
        self.frame.mass_scale = 0.0;
        self.frame.stiff_scale = 1.0;
        if !converged {
            // Never integrate an unconverged state: a residual that large is
            // not an approximate answer, it is a different structure.
            return report;
        }

        for i in 0..n {
            if self.frame.fixed[i] {
                self.velocity[i] = Dof::default();
                self.acceleration[i] = Dof::default();
                continue;
            }
            let d = delta[i];
            let next_a = d
                .scale(1.0 / (beta * h * h))
                .sub(self.velocity[i].scale(1.0 / (beta * h)))
                .sub(self.acceleration[i].scale(1.0 / (2.0 * beta) - 1.0));
            let next_v = self.velocity[i]
                .add(self.acceleration[i].scale(h * (1.0 - gamma)))
                .add(next_a.scale(h * gamma));
            self.velocity[i] = next_v;
            self.acceleration[i] = next_a;
            self.displacement[i] = self.displacement[i].add(d);
        }

        // Failure, on the current total displacement and by the same criteria
        // the static analysis uses.
        let forces = self.frame.element_forces(&self.displacement, &self.stiffness);
        let mut released = 0.0;
        for (i, fe) in forces.iter().enumerate() {
            if self.stiffness[i] <= 0.0 {
                continue;
            }
            if self.exceeds_limit(i, fe) {
                released += self.element_strain_energy(i, fe);
                self.stiffness[i] = 0.0;
                report.broken.push(i);
            }
        }

        report.kinetic = self.kinetic_energy();
        report.strain = self.strain_energy();
        report.released = released;
        let work: f64 = (0..n)
            .map(|i| f.get(i).copied().unwrap_or_default().dot(delta[i]))
            .sum();
        report.dissipated =
            energy_before + work - report.kinetic - report.strain - report.released;
        report.displacement_ratio = self.displacement_ratio();
        report
    }

    /// Strain energy held in one member, from its own internal forces.
    fn element_strain_energy(&self, i: usize, f: &ElementForces) -> f64 {
        let e = &self.frame.elements[i];
        let l = (self.frame.nodes[e.b as usize] - self.frame.nodes[e.a as usize]).norm();
        let k = self.stiffness[i];
        if l <= 0.0 || k <= 0.0 {
            return 0.0;
        }
        let ea = k * e.area();
        let ei = k * e.inertia();
        f.axial * f.axial * l / (2.0 * ea) + f.moment * f.moment * l / (2.0 * ei)
    }

    /// Whether a member has failed, by rupture or by buckling. The same test
    /// the static path applies, so a structure that survives a gust in the
    /// dynamic solver survives the same load held steady.
    fn exceeds_limit(&self, i: usize, f: &ElementForces) -> bool {
        let e = &self.frame.elements[i];
        let m = self.frame.material;
        let ratio = if f.axial > 0.0 { m.tensile_ratio } else { 1.0 };
        let strength = m.rupture * e.integrity * ratio;
        (strength > 0.0 && f.stress > strength) || f.buckling >= 1.0
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.frame
            .lumped
            .iter()
            .zip(&self.velocity)
            .map(|(m, v)| 0.5 * m.mul(*v).dot(*v))
            .sum()
    }

    /// Elastic energy stored in the members, `1/2 x^T K x`.
    pub fn strain_energy(&self) -> f64 {
        let n = self.frame.nodes.len();
        let mut kx = vec![Dof::default(); n];
        let mut f = self.frame.clone();
        f.mass_scale = 0.0;
        f.stiff_scale = 1.0;
        f.apply_operator(&self.displacement, &self.stiffness, &mut kx);
        0.5 * (0..n).map(|i| kx[i].dot(self.displacement[i])).sum::<f64>()
    }

    /// Current position of each node.
    pub fn deformed(&self) -> Vec<Vec3> {
        self.frame
            .nodes
            .iter()
            .zip(&self.displacement)
            .map(|(p, d)| *p + d.t)
            .collect()
    }

    /// Largest nodal displacement measured against the shortest member reaching
    /// that node. Beyond roughly 0.1 the small-displacement assumption is
    /// spent and the restoring force is being overestimated.
    pub fn displacement_ratio(&self) -> f64 {
        let mut scale = vec![f64::INFINITY; self.frame.nodes.len()];
        for e in &self.frame.elements {
            let l = (self.frame.nodes[e.b as usize] - self.frame.nodes[e.a as usize]).norm();
            if l <= 0.0 {
                continue;
            }
            for node in [e.a as usize, e.b as usize] {
                scale[node] = scale[node].min(l);
            }
        }
        (0..self.frame.nodes.len())
            .filter(|&i| scale[i].is_finite() && scale[i] > 0.0)
            .map(|i| self.displacement[i].t.norm() / scale[i])
            .fold(0.0f64, f64::max)
    }

    /// Free vibration period of the dominant mode, estimated by Rayleigh
    /// quotient on the current deformed shape.
    ///
    /// `T = 2 pi sqrt(x^T M x / x^T K x)`. Given a shape near the fundamental
    /// mode this is accurate to a fraction of a percent, and it costs one
    /// operator application rather than an eigensolve — which is what makes it
    /// affordable to ask every frame whether the timestep still resolves the
    /// motion.
    pub fn dominant_period(&self, shape: &[Dof]) -> f64 {
        let n = self.frame.nodes.len();
        if shape.len() < n {
            return 0.0;
        }
        let mut kx = vec![Dof::default(); n];
        let mut f = self.frame.clone();
        f.mass_scale = 0.0;
        f.stiff_scale = 1.0;
        f.apply_operator(shape, &self.stiffness, &mut kx);
        let num: f64 = (0..n)
            .filter(|&i| !self.frame.fixed[i])
            .map(|i| self.frame.lumped[i].mul(shape[i]).dot(shape[i]))
            .sum();
        let den: f64 = (0..n).map(|i| kx[i].dot(shape[i])).sum();
        if den <= 0.0 || num <= 0.0 {
            return 0.0;
        }
        2.0 * std::f64::consts::PI * (num / den).sqrt()
    }
}

/// An orthonormal basis with `e1` along `axis`. Duplicated from `frame` rather
/// than exported from it: it is three lines, and the alternative is a public
/// name that says nothing about frames.
fn orthonormal(axis: Vec3) -> (Vec3, Vec3, Vec3) {
    let e1 = axis.unit();
    let seed = if e1.z.abs() < 0.9 {
        Vec3 { x: 0.0, y: 0.0, z: 1.0 }
    } else {
        Vec3 { x: 1.0, y: 0.0, z: 0.0 }
    };
    let e2 = e1.cross(seed).unit();
    let e3 = e1.cross(e2);
    (e1, e2, e3)
}
