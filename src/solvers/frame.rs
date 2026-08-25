//! Three-dimensional frame analysis: real beam elements, buckling, and plastic
//! redistribution.
//!
//! # What was wrong with springs
//!
//! The first redundant solver modelled each connection as an anisotropic
//! spring: stiff along its axis, soft across it, with the transverse stiffness
//! set to `3EI/L^3`. That is the tip stiffness of a cantilever, and it is
//! exactly right for a cantilever — which is why the obvious validation, a
//! cantilever's `PL^3/3EI` tip deflection, passes and proves nothing.
//!
//! It is wrong for everything else, because it has no rotational degrees of
//! freedom. A real beam carries moments through its joints; a spring network
//! cannot, so it cannot tell a fixed end from a pinned one. The discriminating
//! case is a beam built in at both ends, whose midspan deflection is
//! `PL^3/192EI` — a factor of 64 stiffer than the cantilever, entirely because
//! the fixed ends resist rotation. A translation-only model gets that badly
//! wrong, and `tests/frame.rs` checks it.
//!
//! # The element
//!
//! Standard Euler-Bernoulli 3D frame element: 6 degrees of freedom per node,
//! 12 per element, with axial, torsional and biaxial bending terms. Every
//! member the engine generates has a circular section, so `Iy == Iz` and the
//! choice of principal axes is free — which removes the one genuinely fiddly
//! part of assembling a frame element.
//!
//! Ties are `truss` elements: axial only, no moment transfer. That is what
//! bracing physically is, and modelling it as a full beam would over-stiffen
//! the structure.
//!
//! # Why it is still matrix-free
//!
//! A structure may be rebuilt every frame, so assembling a sparse `6n x 6n`
//! matrix costs more than the solve. Conjugate gradient needs only the product
//! `K u`, which is a loop over elements. The addition here over the previous
//! solver is a Jacobi preconditioner, which matters far more with rotational
//! degrees of freedom than without: translations and rotations differ in units
//! and by orders of magnitude in scale, and unpreconditioned CG spends its
//! iterations discovering that.

use crate::math::{v3, Vec3};
use crate::topology::Material;

/// Poisson's ratio, for the shear modulus `G = E / 2(1+nu)`.
pub const POISSON: f64 = 0.3;

/// One member of a frame.
#[derive(Debug, Clone, Copy)]
pub struct Element {
    pub a: u32,
    pub b: u32,
    /// Circular section radius, m.
    pub radius: f64,
    /// Axial only, no moment transfer. Bracing and ties are pin-jointed.
    pub truss: bool,
    /// Remaining fraction of the section's strength, 0..1.
    pub integrity: f64,
}

impl Element {
    #[inline]
    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.radius * self.radius
    }
    /// Second moment of area of a circular section.
    #[inline]
    pub fn inertia(&self) -> f64 {
        std::f64::consts::PI * self.radius.powi(4) / 4.0
    }
    /// Polar second moment, which for a circle is twice the bending value.
    #[inline]
    pub fn polar(&self) -> f64 {
        2.0 * self.inertia()
    }
    /// Section modulus `I/c`.
    #[inline]
    pub fn section_modulus(&self) -> f64 {
        std::f64::consts::PI * self.radius.powi(3) / 4.0
    }
}

/// A structure ready to analyse: joints, members, and what is held down.
#[derive(Debug, Clone)]
pub struct Frame {
    pub nodes: Vec<Vec3>,
    pub elements: Vec<Element>,
    /// Nodes held against translation *and* rotation — a foundation, not a pin.
    pub fixed: Vec<bool>,
    pub material: Material,

    /// Lumped mass and rotational inertia per node. Empty for a static
    /// analysis, where inertia is by definition irrelevant.
    ///
    /// This is what turns the operator from `K` into `s_m M + s_k K`, and it is
    /// the whole of what [`crate::solvers::dynamics`] needs from this module:
    /// an implicit dynamic step is a static solve against a shifted operator,
    /// so the conjugate gradient, the preconditioner, the plastic
    /// redistribution and the failure criteria are all *the same code*. A
    /// structure cannot then move according to one stiffness and break
    /// according to another.
    pub lumped: Vec<Dof>,
    /// Coefficient on `lumped` in the operator. Zero — the default — leaves a
    /// purely static problem however much mass is attached.
    pub mass_scale: f64,
    /// Coefficient on the stiffness in the operator. One for a static solve;
    /// backward Euler with stiffness-proportional damping uses `1 + beta/h`.
    /// It scales the operator only, never the reported member forces.
    pub stiff_scale: f64,
}

/// Internal forces in one member, in its own local frame.
#[derive(Debug, Clone, Copy, Default)]
pub struct ElementForces {
    /// Positive in tension, N.
    pub axial: f64,
    /// Resultant transverse force, N.
    pub shear: f64,
    /// Largest end moment, N·m.
    pub moment: f64,
    pub torsion: f64,
    /// Peak fibre stress from bending plus axial, Pa.
    pub stress: f64,
    /// Compressive load as a fraction of the Euler critical load. At or above 1
    /// the member buckles, however far its stress is from rupture.
    pub buckling: f64,
}

#[derive(Debug, Clone, Default)]
pub struct FrameSolution {
    pub translation: Vec<Vec3>,
    pub rotation: Vec<Vec3>,
    pub forces: Vec<ElementForces>,
    pub iterations: u32,
    pub converged: bool,
    /// Elements that yielded and shed load to their neighbours.
    pub yielded: usize,
    /// Passes of plastic redistribution performed.
    pub redistributions: u32,
}

/// Six degrees of freedom at a joint.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Dof {
    pub t: Vec3,
    pub r: Vec3,
}

impl Dof {
    #[inline]
    pub fn add(self, o: Dof) -> Dof {
        Dof { t: self.t + o.t, r: self.r + o.r }
    }
    #[inline]
    pub fn sub(self, o: Dof) -> Dof {
        Dof { t: self.t - o.t, r: self.r - o.r }
    }
    #[inline]
    pub fn scale(self, s: f64) -> Dof {
        Dof { t: self.t.scale(s), r: self.r.scale(s) }
    }
    #[inline]
    pub fn dot(self, o: Dof) -> f64 {
        self.t.dot(o.t) + self.r.dot(o.r)
    }
    #[inline]
    pub fn is_finite(self) -> bool {
        self.t.is_finite() && self.r.is_finite()
    }
    /// Componentwise product, for applying a diagonal preconditioner.
    #[inline]
    pub fn mul(self, o: Dof) -> Dof {
        Dof {
            t: v3(self.t.x * o.t.x, self.t.y * o.t.y, self.t.z * o.t.z),
            r: v3(self.r.x * o.r.x, self.r.y * o.r.y, self.r.z * o.r.z),
        }
    }
}

/// An orthonormal basis with `e1` along the member.
fn basis(axis: Vec3) -> (Vec3, Vec3, Vec3) {
    let e1 = axis.unit();
    // Any perpendicular will do: the section is circular, so there are no
    // principal axes to align with.
    let seed = if e1.z.abs() < 0.9 { v3(0.0, 0.0, 1.0) } else { v3(1.0, 0.0, 0.0) };
    let e2 = e1.cross(seed).unit();
    let e3 = e1.cross(e2);
    (e1, e2, e3)
}

impl Frame {
    pub fn new(material: Material) -> Frame {
        Frame {
            nodes: Vec::new(),
            elements: Vec::new(),
            fixed: Vec::new(),
            material,
            lumped: Vec::new(),
            mass_scale: 0.0,
            stiff_scale: 1.0,
        }
    }

    pub fn add_node(&mut self, p: Vec3, fixed: bool) -> u32 {
        self.nodes.push(p);
        self.fixed.push(fixed);
        (self.nodes.len() - 1) as u32
    }

    pub fn add_beam(&mut self, a: u32, b: u32, radius: f64) -> usize {
        self.elements.push(Element { a, b, radius, truss: false, integrity: 1.0 });
        self.elements.len() - 1
    }

    pub fn add_tie(&mut self, a: u32, b: u32, radius: f64) -> usize {
        self.elements.push(Element { a, b, radius, truss: true, integrity: 1.0 });
        self.elements.len() - 1
    }

    #[inline]
    fn length(&self, e: &Element) -> f64 {
        (self.nodes[e.b as usize] - self.nodes[e.a as usize]).norm()
    }

    /// Element stiffness applied to a displacement state, in global
    /// coordinates. Returns the internal force and moment at each end.
    fn element_apply(&self, e: &Element, d1: Dof, d2: Dof, stiffness: f64) -> (Dof, Dof) {
        let d = self.nodes[e.b as usize] - self.nodes[e.a as usize];
        let l = d.norm();
        if l <= 0.0 || stiffness <= 0.0 {
            return (Dof::default(), Dof::default());
        }
        let (e1, e2, e3) = basis(d);
        let ea = stiffness * e.area();
        let ei = stiffness * e.inertia();
        let gj = stiffness / (2.0 * (1.0 + POISSON)) * e.polar();

        // Local components.
        let (u1, v1, w1) = (d1.t.dot(e1), d1.t.dot(e2), d1.t.dot(e3));
        let (u2, v2, w2) = (d2.t.dot(e1), d2.t.dot(e2), d2.t.dot(e3));

        // Axial: present in every element.
        let fx1 = ea / l * (u1 - u2);
        let fx2 = -fx1;

        if e.truss {
            // Pin-jointed: axial only, and no moment reaches the joints.
            return (
                Dof { t: e1.scale(fx1), r: Vec3::ZERO },
                Dof { t: e1.scale(fx2), r: Vec3::ZERO },
            );
        }

        let (tx1, ty1, tz1) = (d1.r.dot(e1), d1.r.dot(e2), d1.r.dot(e3));
        let (tx2, ty2, tz2) = (d2.r.dot(e1), d2.r.dot(e2), d2.r.dot(e3));

        // Torsion about the member axis.
        let mx1 = gj / l * (tx1 - tx2);
        let mx2 = -mx1;

        let l2 = l * l;
        let l3 = l2 * l;
        // Bending in the e1-e2 plane: transverse displacement v, rotation about
        // e3. This is the standard Euler-Bernoulli beam sub-matrix.
        let fy1 = ei * (12.0 * v1 / l3 + 6.0 * tz1 / l2 - 12.0 * v2 / l3 + 6.0 * tz2 / l2);
        let mz1 = ei * (6.0 * v1 / l2 + 4.0 * tz1 / l - 6.0 * v2 / l2 + 2.0 * tz2 / l);
        let fy2 = -fy1;
        let mz2 = ei * (6.0 * v1 / l2 + 2.0 * tz1 / l - 6.0 * v2 / l2 + 4.0 * tz2 / l);

        // Bending in the e1-e3 plane: rotation about e2, with the coupling
        // terms sign-flipped relative to the plane above.
        let fz1 = ei * (12.0 * w1 / l3 - 6.0 * ty1 / l2 - 12.0 * w2 / l3 - 6.0 * ty2 / l2);
        let my1 = ei * (-6.0 * w1 / l2 + 4.0 * ty1 / l + 6.0 * w2 / l2 + 2.0 * ty2 / l);
        let fz2 = -fz1;
        let my2 = ei * (-6.0 * w1 / l2 + 2.0 * ty1 / l + 6.0 * w2 / l2 + 4.0 * ty2 / l);

        (
            Dof {
                t: e1.scale(fx1) + e2.scale(fy1) + e3.scale(fz1),
                r: e1.scale(mx1) + e2.scale(my1) + e3.scale(mz1),
            },
            Dof {
                t: e1.scale(fx2) + e2.scale(fy2) + e3.scale(fz2),
                r: e1.scale(mx2) + e2.scale(my2) + e3.scale(mz2),
            },
        )
    }

    /// The operator `s_m M + s_k K` applied to a displacement state.
    ///
    /// Public because a dynamic solver needs the residual `f - K x` as well as
    /// the solve itself, and re-deriving the element loop to get it would be
    /// the exact duplication this module exists to avoid. Pass
    /// `mass_scale = 0`, `stiff_scale = 1` for the plain stiffness product.
    pub fn apply_operator(&self, u: &[Dof], stiff: &[f64], out: &mut [Dof]) {
        self.apply(u, stiff, out)
    }

    fn apply(&self, u: &[Dof], stiff: &[f64], out: &mut [Dof]) {
        for o in out.iter_mut() {
            *o = Dof::default();
        }
        for (i, e) in self.elements.iter().enumerate() {
            let (a, b) = (e.a as usize, e.b as usize);
            let (f1, f2) = self.element_apply(e, u[a], u[b], stiff[i] * self.stiff_scale);
            out[a] = out[a].add(f1);
            out[b] = out[b].add(f2);
        }
        if self.mass_scale != 0.0 {
            for (i, o) in out.iter_mut().enumerate() {
                let m = self.lumped.get(i).copied().unwrap_or_default();
                *o = o.add(m.mul(u[i]).scale(self.mass_scale));
            }
        }
        for (i, o) in out.iter_mut().enumerate() {
            if self.fixed[i] {
                *o = Dof::default();
            }
        }
    }

    /// Diagonal of the stiffness matrix, for Jacobi preconditioning.
    ///
    /// Without this, CG has to reconcile translational stiffnesses of order
    /// `EA/L` with rotational ones of order `EI/L` — quantities in different
    /// units whose ratio is `L^2/r^2`, four orders of magnitude for a slender
    /// member. It converged eventually and spent hundreds of iterations doing
    /// what one division per degree of freedom does here.
    fn diagonal(&self, stiff: &[f64]) -> Vec<Dof> {
        let mut d = vec![Dof::default(); self.nodes.len()];
        for (i, e) in self.elements.iter().enumerate() {
            let l = self.length(e);
            if l <= 0.0 || stiff[i] <= 0.0 {
                continue;
            }
            let (e1, e2, e3) = basis(self.nodes[e.b as usize] - self.nodes[e.a as usize]);
            let k = stiff[i] * self.stiff_scale;
            let ea = k * e.area() / l;
            let ei = k * e.inertia();
            let gj = k / (2.0 * (1.0 + POISSON)) * e.polar() / l;
            let sq = |v: Vec3| v3(v.x * v.x, v.y * v.y, v.z * v.z);
            let trans = sq(e1).scale(ea)
                + if e.truss { Vec3::ZERO } else { (sq(e2) + sq(e3)).scale(12.0 * ei / l.powi(3)) };
            let rot = if e.truss {
                Vec3::ZERO
            } else {
                sq(e1).scale(gj) + (sq(e2) + sq(e3)).scale(4.0 * ei / l)
            };
            for n in [e.a as usize, e.b as usize] {
                d[n].t += trans;
                d[n].r += rot;
            }
        }
        if self.mass_scale != 0.0 {
            for (i, di) in d.iter_mut().enumerate() {
                let m = self.lumped.get(i).copied().unwrap_or_default();
                *di = di.add(m.scale(self.mass_scale));
            }
        }
        d
    }

    /// Solve for nodal displacements under the given loads.
    ///
    /// `load` is force and moment per node. Returns displacements plus the
    /// internal forces in every member, including its buckling utilisation.
    pub fn solve(&self, load: &[Dof]) -> FrameSolution {
        self.solve_with(load, self.material.ductility > 0.0)
    }

    /// As [`solve`], with plastic redistribution optionally disabled.
    pub fn solve_with(&self, load: &[Dof], plastic: bool) -> FrameSolution {
        let n = self.nodes.len();
        let mut out = FrameSolution::default();
        if n == 0 || self.elements.is_empty() {
            return out;
        }
        // Per-element stiffness, reduced where a member has yielded.
        let mut stiff = vec![self.material.stiffness; self.elements.len()];
        let mut solution = self.solve_elastic(load, &stiff);
        out.iterations = solution.1;
        out.converged = solution.2;

        if plastic && out.converged {
            // Elastic-perfectly-plastic, by secant stiffness. A member past
            // yield cannot carry more, so its stiffness is reduced until the
            // force it attracts matches what it can hold, and the surplus goes
            // to its neighbours. Brittle materials skip this entirely — they
            // fracture rather than redistribute, which is the whole difference
            // between a steel frame sagging and a masonry wall dropping.
            for pass in 0..PLASTIC_PASSES {
                let forces = self.element_forces(&solution.0, &stiff);
                let mut yielded = 0;
                let mut changed = false;
                for (i, f) in forces.iter().enumerate() {
                    let cap = self.material.rupture
                        * self.material.ductility
                        * self.elements[i].integrity;
                    if cap > 0.0 && f.stress > cap {
                        let secant = (cap / f.stress).clamp(1e-4, 1.0);
                        let next = stiff[i] * secant;
                        if (next - stiff[i]).abs() > 1e-9 * stiff[i] {
                            stiff[i] = next;
                            changed = true;
                        }
                        yielded += 1;
                    }
                }
                out.yielded = yielded;
                out.redistributions = pass + 1;
                if !changed {
                    break;
                }
                solution = self.solve_elastic(load, &stiff);
                out.iterations += solution.1;
                if !solution.2 {
                    break;
                }
            }
        }

        out.forces = self.element_forces(&solution.0, &stiff);
        let (t, r): (Vec<Vec3>, Vec<Vec3>) =
            solution.0.iter().map(|d| (d.t, d.r)).unzip();
        out.translation = t;
        out.rotation = r;
        out
    }

    /// Jacobi-preconditioned conjugate gradient. Returns `(u, iterations,
    /// converged)`.
    fn solve_elastic(&self, load: &[Dof], stiff: &[f64]) -> (Vec<Dof>, u32, bool) {
        self.solve_elastic_opt(load, stiff, true)
    }

    /// As above, with preconditioning optionally disabled — which exists so
    /// that the claim "preconditioning matters here" can be measured rather
    /// than asserted against a magic iteration count.
    pub fn solve_elastic_opt(
        &self,
        load: &[Dof],
        stiff: &[f64],
        precondition: bool,
    ) -> (Vec<Dof>, u32, bool) {
        let n = self.nodes.len();
        let diag = self.diagonal(stiff);
        // Nodes touched by no element, and fixed nodes, are removed from the
        // system: they have an empty row but may carry a load, which is
        // inconsistent and which CG can never satisfy.
        let inv: Vec<Dof> = (0..n)
            .map(|i| {
                if self.fixed[i] {
                    return Dof::default();
                }
                let d = diag[i];
                // Unpreconditioned still masks out degrees of freedom with no
                // stiffness behind them — that is a well-posedness question,
                // not a convergence one.
                let f = |x: f64| {
                    if x > 0.0 {
                        if precondition { 1.0 / x } else { 1.0 }
                    } else {
                        0.0
                    }
                };
                Dof {
                    t: v3(f(d.t.x), f(d.t.y), f(d.t.z)),
                    r: v3(f(d.r.x), f(d.r.y), f(d.r.z)),
                }
            })
            .collect();

        let b: Vec<Dof> = (0..n)
            .map(|i| {
                if self.fixed[i] {
                    Dof::default()
                } else {
                    // Zero out any component with no stiffness behind it.
                    let mask = Dof {
                        t: v3(
                            (inv[i].t.x > 0.0) as u8 as f64,
                            (inv[i].t.y > 0.0) as u8 as f64,
                            (inv[i].t.z > 0.0) as u8 as f64,
                        ),
                        r: v3(
                            (inv[i].r.x > 0.0) as u8 as f64,
                            (inv[i].r.y > 0.0) as u8 as f64,
                            (inv[i].r.z > 0.0) as u8 as f64,
                        ),
                    };
                    load.get(i).copied().unwrap_or_default().mul(mask)
                }
            })
            .collect();

        // Preconditioner. The tree factorisation is exact for a structure
        // whose load path is a tree, so conjugate gradient finishes in one
        // iteration; a braced one keeps only the braces outside the factor and
        // converges in a handful. Jacobi remains the fallback, and remains what
        // `precondition = false` turns off, so the comparison the tests make
        // stays a comparison of preconditioners rather than of code paths.
        // Applying the factorisation costs roughly an order of magnitude more
        // than a diagonal division, so it has to remove about that much from
        // the iteration count to be worth it. What predicts that is the
        // structure's degree of static indeterminacy: on a tree the factor is
        // exact and removes everything, and on a frame carrying nearly a brace
        // per joint it removes a factor of two and loses on the exchange.
        // Above this ratio the diagonal is the better preconditioner, which is
        // why it stays.
        let factor = if precondition && self.redundancy() * FACTOR_REDUNDANCY_LIMIT <= n {
            Some(self.factor(stiff))
        } else {
            None
        };
        let mut u = vec![Dof::default(); n];
        let mut r = b.clone();
        let mut z = vec![Dof::default(); n];
        let precondition_into = |r: &[Dof], z: &mut [Dof]| match &factor {
            Some(f) if !f.is_empty() => {
                f.solve(r, z);
                // The factor knows nothing about degrees of freedom that carry
                // no stiffness at all; the diagonal mask does, and applying it
                // after keeps them out of the search directions.
                for i in 0..n {
                    if inv[i] == Dof::default() {
                        z[i] = Dof::default();
                    }
                }
            }
            _ => {
                for i in 0..n {
                    z[i] = r[i].mul(inv[i]);
                }
            }
        };
        precondition_into(&r, &mut z);
        let mut p = z.clone();
        let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a.dot(*b)).sum();
        let r0: f64 = r.iter().map(|a| a.dot(*a)).sum();
        if r0 <= 0.0 {
            return (u, 0, true);
        }
        let tol2 = r0 * 1e-18;
        let mut ap = vec![Dof::default(); n];
        let mut iters = 0u32;
        let mut converged = false;
        let max_iter = (n * 6).clamp(64, 6000);

        for _ in 0..max_iter {
            iters += 1;
            self.apply(&p, stiff, &mut ap);
            let denom: f64 = p.iter().zip(&ap).map(|(a, b)| a.dot(*b)).sum();
            if !(denom.abs() > 1e-300) || !denom.is_finite() {
                break;
            }
            let alpha = rz / denom;
            if !alpha.is_finite() {
                break;
            }
            for i in 0..n {
                u[i] = u[i].add(p[i].scale(alpha));
                r[i] = r[i].sub(ap[i].scale(alpha));
            }
            let rr: f64 = r.iter().map(|a| a.dot(*a)).sum();
            if !rr.is_finite() {
                break;
            }
            if rr <= tol2 {
                converged = true;
                break;
            }
            precondition_into(&r, &mut z);
            let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a.dot(*b)).sum();
            let beta = rz_new / rz;
            if !beta.is_finite() {
                break;
            }
            for i in 0..n {
                p[i] = z[i].add(p[i].scale(beta));
            }
            rz = rz_new;
        }
        if !u.iter().all(|d| d.is_finite()) {
            return (vec![Dof::default(); n], iters, false);
        }
        (u, iters, converged)
    }

    /// Internal forces in every member, from a displacement state.
    pub fn element_forces(&self, u: &[Dof], stiff: &[f64]) -> Vec<ElementForces> {
        self.elements
            .iter()
            .enumerate()
            .map(|(i, e)| {
                let (a, b) = (e.a as usize, e.b as usize);
                let l = self.length(e);
                if l <= 0.0 {
                    return ElementForces::default();
                }
                let (f1, _f2) = self.element_apply(e, u[a], u[b], stiff[i]);
                let (e1, _, _) = basis(self.nodes[b] - self.nodes[a]);
                // Sign convention: positive axial is tension.
                let axial = -f1.t.dot(e1);
                let shear = (f1.t - e1.scale(f1.t.dot(e1))).norm();
                let torsion = f1.r.dot(e1).abs();
                let bending = (f1.r - e1.scale(f1.r.dot(e1))).norm();
                // Both ends of a beam carry a moment; the larger governs.
                let moment = bending.max(if e.truss { 0.0 } else { shear * l * 0.5 });

                let area = e.area().max(1e-30);
                let section = e.section_modulus().max(1e-30);
                let stress = moment / section + axial.abs() / area;

                // Euler buckling. A slender member in compression fails at a
                // load far below its material strength, and nothing in a stress
                // check will ever notice — it is a stability failure, not a
                // strength one. Omitting it means a hundred-metre column
                // reports 57% utilised while it is folding up.
                let buckling = if axial < 0.0 {
                    let k = effective_length_factor(e.truss);
                    let critical =
                        std::f64::consts::PI.powi(2) * stiff[i] * e.inertia() / (k * l).powi(2);
                    if critical > 0.0 {
                        -axial / critical
                    } else {
                        f64::INFINITY
                    }
                } else {
                    0.0
                };

                ElementForces { axial, shear, moment, torsion, stress, buckling }
            })
            .collect()
    }

    /// Total elastic strain energy, useful for comparing load paths.
    pub fn strain_energy(&self, solution: &FrameSolution) -> f64 {
        solution
            .forces
            .iter()
            .enumerate()
            .map(|(i, f)| {
                let e = &self.elements[i];
                let l = self.length(e);
                if l <= 0.0 {
                    return 0.0;
                }
                let ea = self.material.stiffness * e.area();
                let ei = self.material.stiffness * e.inertia();
                f.axial * f.axial * l / (2.0 * ea) + f.moment * f.moment * l / (2.0 * ei)
            })
            .sum()
    }
}

/// Effective length factor for buckling.
///
/// A pin-jointed tie buckles as a pinned-pinned strut (`K = 1`); a member built
/// into a frame at both ends is stiffer (`K` nearer 0.7 in theory, but joint
/// flexibility eats most of that, so 0.85 is the honest number for a real
/// connection).
fn effective_length_factor(truss: bool) -> f64 {
    if truss {
        1.0
    } else {
        0.85
    }
}

const PLASTIC_PASSES: u32 = 6;

// ---------------------------------------------------------------------------
// The tree preconditioner
// ---------------------------------------------------------------------------

/// A six-by-six block: one joint's coupling to another.
#[derive(Debug, Clone, Copy)]
pub struct Mat6([[f64; 6]; 6]);

impl Mat6 {
    pub const ZERO: Mat6 = Mat6([[0.0; 6]; 6]);

    #[inline]
    fn sub_assign(&mut self, o: &Mat6) {
        for i in 0..6 {
            for j in 0..6 {
                self.0[i][j] -= o.0[i][j];
            }
        }
    }

    #[inline]
    fn add_assign(&mut self, o: &Mat6) {
        for i in 0..6 {
            for j in 0..6 {
                self.0[i][j] += o.0[i][j];
            }
        }
    }

    fn mul(&self, o: &Mat6) -> Mat6 {
        let mut out = Mat6::ZERO;
        for i in 0..6 {
            for k in 0..6 {
                let a = self.0[i][k];
                if a == 0.0 {
                    continue;
                }
                for j in 0..6 {
                    out.0[i][j] += a * o.0[k][j];
                }
            }
        }
        out
    }

    fn apply(&self, v: Dof) -> Dof {
        let x = [v.t.x, v.t.y, v.t.z, v.r.x, v.r.y, v.r.z];
        let mut y = [0.0f64; 6];
        for i in 0..6 {
            let mut s = 0.0;
            for j in 0..6 {
                s += self.0[i][j] * x[j];
            }
            y[i] = s;
        }
        Dof { t: v3(y[0], y[1], y[2]), r: v3(y[3], y[4], y[5]) }
    }

    /// Gauss-Jordan inverse with partial pivoting.
    ///
    /// A degree of freedom with no stiffness and no mass behind it — a node
    /// orphaned by a break — gives a singular block. Its row and column come
    /// back zero rather than infinite, which says "this direction is not
    /// determined" and is the same treatment the diagonal masking gives it.
    fn inverse(&self) -> Mat6 {
        let mut a = self.0;
        let mut inv = [[0.0f64; 6]; 6];
        for i in 0..6 {
            inv[i][i] = 1.0;
        }
        // Scale for deciding what counts as a zero pivot, from the block itself.
        let mut scale = 0.0f64;
        for row in &a {
            for v in row {
                scale = scale.max(v.abs());
            }
        }
        let tiny = scale * 1e-14;
        for col in 0..6 {
            let mut pivot = col;
            for r in col + 1..6 {
                if a[r][col].abs() > a[pivot][col].abs() {
                    pivot = r;
                }
            }
            if a[pivot][col].abs() <= tiny {
                // Undetermined direction: zero the row so it contributes
                // nothing rather than an unbounded correction.
                for j in 0..6 {
                    inv[col][j] = 0.0;
                }
                continue;
            }
            a.swap(col, pivot);
            inv.swap(col, pivot);
            let d = a[col][col];
            for j in 0..6 {
                a[col][j] /= d;
                inv[col][j] /= d;
            }
            for r in 0..6 {
                if r == col {
                    continue;
                }
                let f = a[r][col];
                if f == 0.0 {
                    continue;
                }
                for j in 0..6 {
                    a[r][j] -= f * a[col][j];
                    inv[r][j] -= f * inv[col][j];
                }
            }
        }
        Mat6(inv)
    }
}

/// An exact block factorisation of the structure's spanning forest.
///
/// # Why this exists
///
/// The load path of every structure the engine generates is a *tree*, and that
/// is what makes the static analysis an O(n) reverse pass with no solve at all.
/// The same fact does something for the dynamic solve that is easy to miss: a
/// tree has a **perfect elimination ordering**. Eliminate the leaves first and
/// each one's only surviving coupling is to its parent, so the factorisation
/// produces no fill-in whatever and costs one six-by-six inverse per joint.
///
/// For a determinate structure that makes this an *exact* solve, and
/// preconditioned conjugate gradient converges in a single iteration. For a
/// braced one the ties are the only edges outside the forest, so it remains a
/// very good approximation and the iteration count collapses.
///
/// This is not an optimisation bolted on afterwards. It is the same structural
/// fact the O(n) static pass rests on, applied to a different question.
///
/// # Why it was needed
///
/// Jacobi preconditioning is fine for a squat braced frame and hopeless for a
/// tree: a slender chain of `n` beam elements has a condition number growing
/// like `n^4`, so a 2000-member tree took 3645 iterations and 700 ms for one
/// substep. At twenty updates a second that is not a solver, it is a
/// screensaver.
#[derive(Debug, Clone, Default)]
pub struct TreeFactor {
    /// Joints in elimination order: every child before its parent.
    order: Vec<u32>,
    /// Parent of each joint in the spanning forest, or `u32::MAX` for a root.
    parent: Vec<u32>,
    /// The coupling block to that parent.
    to_parent: Vec<Mat6>,
    /// Inverse of the eliminated diagonal block.
    inverse: Vec<Mat6>,
    /// Joints outside the system: held down, or reached by nothing.
    excluded: Vec<bool>,
}

impl TreeFactor {
    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// Solve the factorised system for `r`, in place of a Jacobi division.
    ///
    /// `A = L D L^T` with `L_{p,c} = A_pc D_c^-1` for each child `c` of `p` —
    /// the only off-diagonal entries a tree has. Forward-substitute towards the
    /// roots, apply the inverse diagonal, back-substitute towards the leaves.
    pub fn solve(&self, r: &[Dof], out: &mut [Dof]) {
        for i in 0..self.parent.len() {
            out[i] = if self.excluded[i] { Dof::default() } else { r[i] };
        }
        // Forward: children first, each pushing its eliminated load onto its
        // parent.
        for &node in &self.order {
            let i = node as usize;
            let p = self.parent[i];
            if p == u32::MAX {
                continue;
            }
            let lifted = self.to_parent[i].transposed().apply(self.inverse[i].apply(out[i]));
            out[p as usize] = out[p as usize].sub(lifted);
        }
        // Back: parents first, each correcting its children.
        for &node in self.order.iter().rev() {
            let i = node as usize;
            let p = self.parent[i];
            let mut y = out[i];
            if p != u32::MAX {
                y = y.sub(self.to_parent[i].apply(out[p as usize]));
            }
            out[i] = self.inverse[i].apply(y);
        }
    }
}

impl Frame {
    /// Factorise the structure's spanning forest.
    ///
    /// Members outside the forest — the braces, which are exactly what makes a
    /// structure redundant — contribute their diagonal stiffness but not their
    /// coupling. Dropping them is what keeps the factorisation fill-free; it
    /// costs nothing in correctness because this is only ever a preconditioner.
    pub fn factor(&self, stiff: &[f64]) -> TreeFactor {
        let n = self.nodes.len();
        if n == 0 {
            return TreeFactor::default();
        }
        let mut diag = vec![Mat6::ZERO; n];
        // Adjacency, with the element index so blocks can be recovered.
        let mut adjacent: Vec<Vec<(u32, usize)>> = vec![Vec::new(); n];
        let mut blocks: Vec<(Mat6, Mat6, Mat6, Mat6)> = Vec::with_capacity(self.elements.len());
        for (i, e) in self.elements.iter().enumerate() {
            let k = stiff.get(i).copied().unwrap_or(0.0) * self.stiff_scale;
            let (aa, ab, ba, bb) = self.element_blocks(e, k);
            let (a, b) = (e.a as usize, e.b as usize);
            diag[a].add_assign(&aa);
            diag[b].add_assign(&bb);
            adjacent[a].push((e.b, i));
            adjacent[b].push((e.a, i));
            blocks.push((aa, ab, ba, bb));
        }
        if self.mass_scale != 0.0 {
            for (i, d) in diag.iter_mut().enumerate() {
                let m = self.lumped.get(i).copied().unwrap_or_default();
                let v = [m.t.x, m.t.y, m.t.z, m.r.x, m.r.y, m.r.z];
                for k in 0..6 {
                    d.0[k][k] += v[k] * self.mass_scale;
                }
            }
        }

        // The *maximum* spanning forest, by axial stiffness. Which members are
        // left out is the whole question: what is dropped is what the
        // preconditioner cannot see, and leaving out a column costs far more
        // than leaving out a slender stay.
        //
        // A breadth-first forest is not good enough, and fails in a way that
        // looks harmless. A stay skipping five joints along a chain is a
        // shortcut, so breadth-first discovery reaches the far joint through
        // the stay and drops a *chain* element instead — severing the load path
        // in the factorisation to keep a member a hundred times softer. On a
        // 200-element chain with ten stays that took the iteration count from
        // one to eleven hundred. Prim's algorithm takes the stiffest edge
        // available at every step and cannot make that mistake.
        let weight: Vec<f64> = self
            .elements
            .iter()
            .enumerate()
            .map(|(i, e)| {
                let l = self.length(e).max(1e-30);
                stiff.get(i).copied().unwrap_or(0.0) * e.area() / l
            })
            .collect();

        let mut parent = vec![u32::MAX; n];
        let mut to_parent = vec![Mat6::ZERO; n];
        let mut seen: Vec<bool> = (0..n).map(|i| self.fixed[i]).collect();
        let mut order: Vec<u32> = Vec::with_capacity(n);
        // Ordered by weight; `Reverse` is unnecessary because the heap is a
        // max-heap and the stiffest edge is the one to take.
        let mut frontier: std::collections::BinaryHeap<(OrderedWeight, u32, u32, u32)> =
            std::collections::BinaryHeap::new();
        for root in 0..n {
            if seen[root] {
                continue;
            }
            seen[root] = true;
            order.push(root as u32);
            frontier.clear();
            let mut push = |frontier: &mut std::collections::BinaryHeap<_>, node: u32| {
                for &(other, element) in &adjacent[node as usize] {
                    frontier.push((OrderedWeight(weight[element]), other, node, element as u32));
                }
            };
            push(&mut frontier, root as u32);
            while let Some((_, other, from, element)) = frontier.pop() {
                let j = other as usize;
                if seen[j] {
                    continue;
                }
                seen[j] = true;
                parent[j] = from;
                let (_aa, ab, ba, _bb) = blocks[element as usize];
                // The block on row `j`, column `from`.
                to_parent[j] = if self.elements[element as usize].a as usize == j { ab } else { ba };
                order.push(other);
                push(&mut frontier, other);
            }
        }
        // Children before parents. Prim discovers every joint after its parent,
        // so reversing the discovery order is exactly that.
        order.reverse();

        let mut inverse = vec![Mat6::ZERO; n];
        for &node in &order {
            let i = node as usize;
            inverse[i] = diag[i].inverse();
            let p = parent[i];
            if p != u32::MAX {
                let t = to_parent[i];
                let update = t.transposed().mul(&inverse[i].mul(&t));
                let mut d = diag[p as usize];
                d.sub_assign(&update);
                diag[p as usize] = d;
            }
        }

        let excluded: Vec<bool> = (0..n).map(|i| self.fixed[i]).collect();
        TreeFactor { order, parent, to_parent, inverse, excluded }
    }

    /// Members outside the spanning forest: the structure's degree of static
    /// indeterminacy, and what the factorisation cannot see.
    ///
    /// Members reaching a foundation do not count. A held-down joint has no
    /// unknowns, so it is eliminated from the system entirely and its members
    /// contribute stiffness to the free end and nothing else — there is no
    /// coupling for the factorisation to drop.
    pub fn redundancy(&self) -> usize {
        let n = self.nodes.len();
        let mut adjacent: Vec<Vec<u32>> = vec![Vec::new(); n];
        for e in &self.elements {
            adjacent[e.a as usize].push(e.b);
            adjacent[e.b as usize].push(e.a);
        }
        let mut seen: Vec<bool> = (0..n).map(|i| self.fixed[i]).collect();
        let mut in_forest = 0usize;
        let mut queue: std::collections::VecDeque<u32> = std::collections::VecDeque::new();
        for root in 0..n {
            if seen[root] {
                continue;
            }
            seen[root] = true;
            queue.push_back(root as u32);
            while let Some(node) = queue.pop_front() {
                for &other in &adjacent[node as usize] {
                    if seen[other as usize] {
                        continue;
                    }
                    seen[other as usize] = true;
                    in_forest += 1;
                    queue.push_back(other);
                }
            }
        }
        let free_members = self
            .elements
            .iter()
            .filter(|e| !self.fixed[e.a as usize] && !self.fixed[e.b as usize])
            .count();
        free_members.saturating_sub(in_forest)
    }

    /// The four six-by-six blocks of one element's stiffness, by applying it to
    /// each unit degree of freedom in turn. Twelve applications of a routine
    /// that already exists, rather than a second derivation of the same algebra
    /// that could drift away from it.
    fn element_blocks(&self, e: &Element, stiff: f64) -> (Mat6, Mat6, Mat6, Mat6) {
        let (mut aa, mut ab, mut ba, mut bb) =
            (Mat6::ZERO, Mat6::ZERO, Mat6::ZERO, Mat6::ZERO);
        let unit = |k: usize| {
            let mut d = Dof::default();
            match k {
                0 => d.t.x = 1.0,
                1 => d.t.y = 1.0,
                2 => d.t.z = 1.0,
                3 => d.r.x = 1.0,
                4 => d.r.y = 1.0,
                _ => d.r.z = 1.0,
            }
            d
        };
        let put = |m: &mut Mat6, col: usize, d: Dof| {
            let v = [d.t.x, d.t.y, d.t.z, d.r.x, d.r.y, d.r.z];
            for row in 0..6 {
                m.0[row][col] = v[row];
            }
        };
        for k in 0..6 {
            let (f1, f2) = self.element_apply(e, unit(k), Dof::default(), stiff);
            put(&mut aa, k, f1);
            put(&mut ba, k, f2);
            let (g1, g2) = self.element_apply(e, Dof::default(), unit(k), stiff);
            put(&mut ab, k, g1);
            put(&mut bb, k, g2);
        }
        (aa, ab, ba, bb)
    }
}

impl Mat6 {
    fn transposed(&self) -> Mat6 {
        let mut out = Mat6::ZERO;
        for i in 0..6 {
            for j in 0..6 {
                out.0[i][j] = self.0[j][i];
            }
        }
        out
    }
}


/// How redundant a structure may be before the tree factorisation stops paying
/// for itself, as a divisor of the joint count. Measured: a frame with a brace
/// per joint runs faster on the diagonal.
const FACTOR_REDUNDANCY_LIMIT: usize = 4;


/// A total order on stiffness, for the spanning-forest heap. Stiffnesses are
/// finite and non-negative by construction, so `total_cmp` is a total order
/// here without any of the caveats it usually carries.
#[derive(Debug, Clone, Copy, PartialEq)]
struct OrderedWeight(f64);

impl Eq for OrderedWeight {}
impl Ord for OrderedWeight {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}
impl PartialOrd for OrderedWeight {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
