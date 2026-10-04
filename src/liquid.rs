//! Liquids: rigid molecules in a periodic box, moved by a law between their
//! sites — `docs/PLAY.md` Phase 6, E8.
//!
//! # What this is for
//!
//! A boiling point is where a liquid's vapour presses at one atmosphere. The
//! owner chose to find it by direct coexistence: a slab of liquid in a long
//! periodic box with its own vapour either side, simulated in time, the
//! vapour's density read off at each temperature. This module is the
//! simulation; the law between molecules comes from the electronic structure
//! (route A's fitted site-site form, or E6's derived one) and is handed in.
//!
//! # Rigid molecules, by the owner's decision
//!
//! Each molecule keeps the shape E5 relaxed it to. Its state is a centre of
//! mass and its momentum, and an orientation (a quaternion taking its
//! principal axes to the box's) with its angular momentum about those axes.
//! Rotation is integrated by the splitting of Dullweber, Leimkuhler and
//! McLachlan (1997): a free rotor's step as exact rotations about one
//! principal axis at a time, `x(h/2) y(h/2) z(h) y(h/2) x(h/2)`, which is
//! symplectic and keeps the orientation a rotation to round-off — so energy
//! does not drift with the step count, which is what a coexistence run of a
//! nanosecond needs.
//!
//! # The box and the cut-off
//!
//! Molecules interact through the nearest periodic image of their centres,
//! every site of one with every site of the other at that one shift, so a
//! molecule's charges stay together. Between `r_on` and `r_cut` the pair's
//! whole energy is switched off smoothly on the centres' separation (a
//! quintic with continuous first and second derivatives), so forces are
//! continuous and energy is conserved. How far `r_cut` must reach is a
//! convergence question for each liquid — dispersion's tail is long — and is
//! measured, not assumed.
//!
//! # Temperature
//!
//! The engine's rule for thermostats: deterministic, its noise drawn from a
//! seeded stream, so a run replays exactly. A Langevin step on each
//! molecule's momentum and angular momentum, Ornstein-Uhlenbeck exact for the
//! step, between the two half-kicks.
//!
//! Atomic units throughout: bohr, hartree, electron masses, and time in
//! hbar / hartree (2.4189e-17 s).

use crate::math::{Quat, Vec3};
use crate::rng::{Purpose, Stream};

/// Boltzmann's constant in hartree per kelvin.
pub const K_B: f64 = 3.166811563e-6;

/// One unified atomic mass unit in electron masses.
pub const AMU: f64 = 1822.888486;

/// A kind of rigid molecule: its sites in its own principal frame, about its
/// centre of mass, each with a type the law reads; its mass and principal
/// moments of inertia.
#[derive(Debug, Clone)]
pub struct Kind {
    pub sites: Vec<Vec3>,
    pub types: Vec<usize>,
    pub mass: f64,
    pub inertia: [f64; 3],
}

impl Kind {
    /// From atoms (positions in bohr, masses in electron masses, a type for
    /// each): centred on the centre of mass and turned onto the principal
    /// axes of inertia.
    pub fn from_atoms(positions: &[[f64; 3]], masses: &[f64], types: &[usize]) -> Kind {
        let mass: f64 = masses.iter().sum();
        let mut com = Vec3::ZERO;
        for (p, m) in positions.iter().zip(masses) {
            com += Vec3 { x: p[0], y: p[1], z: p[2] }.scale(m / mass);
        }
        let rel: Vec<Vec3> = positions.iter().map(|p| Vec3 { x: p[0], y: p[1], z: p[2] } - com).collect();
        let mut t = crate::electrons::linalg::Matrix::zeros(3);
        for (r, m) in rel.iter().zip(masses) {
            let c = [r.x, r.y, r.z];
            let r2 = r.norm2();
            for i in 0..3 {
                for j in 0..3 {
                    let v = t.get(i, j) + m * (if i == j { r2 } else { 0.0 } - c[i] * c[j]);
                    t.set(i, j, v);
                }
            }
        }
        let (moments, axes) = crate::electrons::linalg::eigh(&t);
        // Columns of `axes` are the principal axes in the input frame; a
        // site's principal-frame coordinates are its projections on them.
        // A right-handed frame, so that the orientation is a rotation.
        let col = |k: usize| Vec3 { x: axes.get(0, k), y: axes.get(1, k), z: axes.get(2, k) };
        let (e0, e1) = (col(0), col(1));
        let e2 = e0.cross(e1);
        let sites = rel.iter().map(|r| Vec3 { x: r.dot(e0), y: r.dot(e1), z: r.dot(e2) }).collect();
        Kind { sites, types: types.to_vec(), mass, inertia: [moments[0], moments[1], moments[2]] }
    }
}

/// What the molecules feel: an energy between two sites of given types at a
/// distance, and its derivative with respect to that distance.
pub trait SiteLaw: Sync {
    fn site_pair(&self, type_a: usize, type_b: usize, r: f64) -> (f64, f64);
}

/// The quintic switch: 1 below `on`, 0 above `cut`, with continuous first and
/// second derivatives. Returns the switch and its derivative.
pub fn switch(r: f64, on: f64, cut: f64) -> (f64, f64) {
    if r <= on {
        return (1.0, 0.0);
    }
    if r >= cut {
        return (0.0, 0.0);
    }
    let x = (r - on) / (cut - on);
    let s = 1.0 - x * x * x * (10.0 - 15.0 * x + 6.0 * x * x);
    let ds = -30.0 * x * x * (1.0 - x) * (1.0 - x) / (cut - on);
    (s, ds)
}

/// Molecules of one kind in a periodic box.
#[derive(Debug, Clone)]
pub struct Liquid {
    pub kind: Kind,
    /// The box's edges, bohr; periodic in all three.
    pub cell: [f64; 3],
    pub com: Vec<Vec3>,
    /// Each molecule's orientation: principal frame to box.
    pub orientation: Vec<Quat>,
    pub momentum: Vec<Vec3>,
    /// Angular momentum about the principal axes.
    pub spin: Vec<Vec3>,
    pub r_on: f64,
    pub r_cut: f64,
}

/// Forces and torques on every molecule, and the potential energy.
pub struct Forces {
    pub energy: f64,
    pub force: Vec<Vec3>,
    /// Torque in the box frame.
    pub torque: Vec<Vec3>,
}

impl Liquid {
    /// The nearest image of `d`.
    fn image(&self, mut d: Vec3) -> Vec3 {
        let c = self.cell;
        d.x -= c[0] * (d.x / c[0]).round();
        d.y -= c[1] * (d.y / c[1]).round();
        d.z -= c[2] * (d.z / c[2]).round();
        d
    }

    /// Every site of molecule `i`, in the box frame, relative to its centre.
    fn arms(&self, i: usize) -> Vec<Vec3> {
        self.kind.sites.iter().map(|s| self.orientation[i].rotate(*s)).collect()
    }

    /// The potential energy, forces and torques, through `law`, every pair of
    /// molecules within the cut-off of each other's nearest image. Molecules
    /// are spread across threads by row; each row's sums are its own and are
    /// added in order, so the result does not depend on the thread count.
    pub fn forces(&self, law: &dyn SiteLaw) -> Forces {
        let n = self.com.len();
        let arms: Vec<Vec<Vec3>> = (0..n).map(|i| self.arms(i)).collect();
        let types = &self.kind.types;
        let job = |idx: &mut dyn Iterator<Item = usize>| {
            idx.map(|i| {
                let (mut e, mut f, mut t) = (0.0, Vec3::ZERO, Vec3::ZERO);
                for j in 0..n {
                    if j == i {
                        continue;
                    }
                    let d = self.image(self.com[i] - self.com[j]);
                    let dist = d.norm();
                    if dist >= self.r_cut {
                        continue;
                    }
                    let (s, ds) = switch(dist, self.r_on, self.r_cut);
                    let mut u = 0.0;
                    let (mut fi, mut ti) = (Vec3::ZERO, Vec3::ZERO);
                    for (a, (pa, ta)) in arms[i].iter().zip(types).enumerate() {
                        let _ = a;
                        for (pb, tb) in arms[j].iter().zip(types) {
                            let r = d + *pa - *pb;
                            let rr = r.norm();
                            let (v, dv) = law.site_pair(*ta, *tb, rr);
                            u += v;
                            // Force on site a of molecule i.
                            let fa = r.scale(-dv / rr);
                            fi += fa;
                            ti += pa.cross(fa);
                        }
                    }
                    // Half the pair's energy to each row.
                    e += 0.5 * s * u;
                    f += fi.scale(s) - d.scale(u * ds / dist);
                    t += ti.scale(s);
                }
                (i, e, f, t)
            }).collect::<Vec<_>>()
        };
        let mut energy_rows = vec![0.0; n];
        let mut force = vec![Vec3::ZERO; n];
        let mut torque = vec![Vec3::ZERO; n];
        for part in crate::electrons::scf::parallel_interleaved(n, &job) {
            for (i, e, f, t) in part {
                energy_rows[i] = e;
                force[i] = f;
                torque[i] = t;
            }
        }
        Forces { energy: energy_rows.iter().sum(), force, torque }
    }

    /// Kinetic energy: translational and rotational.
    pub fn kinetic(&self) -> (f64, f64) {
        let tr: f64 = self.momentum.iter().map(|p| p.norm2() / (2.0 * self.kind.mass)).sum();
        let i = self.kind.inertia;
        let rot: f64 = self.spin.iter().map(|l| l.x * l.x / (2.0 * i[0]) + l.y * l.y / (2.0 * i[1]) + l.z * l.z / (2.0 * i[2])).sum();
        (tr, rot)
    }

    /// Kick: momenta by forces and torques for `h`.
    fn kick(&mut self, forces: &Forces, h: f64) {
        for i in 0..self.com.len() {
            self.momentum[i] += forces.force[i].scale(h);
            // The torque in the molecule's principal frame.
            self.spin[i] += self.orientation[i].conjugate().rotate(forces.torque[i]).scale(h);
        }
    }

    /// Drift: centres by momenta, orientations by the free-rotor splitting.
    fn drift(&mut self, h: f64) {
        let inertia = self.kind.inertia;
        for i in 0..self.com.len() {
            self.com[i] += self.momentum[i].scale(h / self.kind.mass);
            let (mut q, mut l) = (self.orientation[i], self.spin[i]);
            for (axis, frac) in [(0usize, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)] {
                let comp = [l.x, l.y, l.z][axis];
                let angle = comp / inertia[axis] * h * frac;
                let e = [Vec3 { x: 1.0, y: 0.0, z: 0.0 }, Vec3 { x: 0.0, y: 1.0, z: 0.0 }, Vec3 { x: 0.0, y: 0.0, z: 1.0 }][axis];
                // The body turns by `angle` about its own axis; its angular
                // momentum, seen from the body, turns the other way.
                let r = Quat::from_axis_angle(e, angle);
                q = q.then(r);
                l = r.conjugate().rotate(l);
            }
            self.orientation[i] = q.unit();
            self.spin[i] = l;
        }
    }

    /// One step of `h`: velocity Verlet with the rotational splitting, and,
    /// with `bath = Some((temperature, friction, stream))`, a Langevin step on
    /// every momentum between the half-kicks. Returns the forces at the end,
    /// for the next step.
    pub fn step(&mut self, law: &dyn SiteLaw, forces: Forces, h: f64, bath: Option<(f64, f64, &mut Stream)>) -> Forces {
        self.kick(&forces, 0.5 * h);
        self.drift(0.5 * h);
        if let Some((temperature, friction, stream)) = bath {
            let c = (-friction * h).exp();
            let kt = K_B * temperature;
            let m = self.kind.mass;
            let inertia = self.kind.inertia;
            for i in 0..self.com.len() {
                let xi = stream.normal3();
                self.momentum[i] = self.momentum[i].scale(c) + xi.scale(((1.0 - c * c) * m * kt).sqrt());
                let zeta = stream.normal3();
                let l = self.spin[i];
                self.spin[i] = Vec3 {
                    x: c * l.x + zeta.x * ((1.0 - c * c) * inertia[0] * kt).sqrt(),
                    y: c * l.y + zeta.y * ((1.0 - c * c) * inertia[1] * kt).sqrt(),
                    z: c * l.z + zeta.z * ((1.0 - c * c) * inertia[2] * kt).sqrt(),
                };
            }
        }
        self.drift(0.5 * h);
        let next = self.forces(law);
        self.kick(&next, 0.5 * h);
        next
    }

    /// A stream for this run's thermal noise, seeded from a number so a run
    /// replays exactly.
    pub fn noise(seed: u64) -> Stream {
        Stream::at(seed, 0, 0, Purpose::ThermalNoise)
    }
}
