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

    /// If the law has induced dipoles (`induction.rs`): the polarisability of
    /// each site type (0 for a site that does not polarise) and the charge of
    /// each, which are the fields' sources. The pair terms above carry the
    /// permanent charges; this adds what the molecules do to each other's.
    fn induction(&self) -> Option<(&[f64], &[f64])> {
        None
    }
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
    /// The molecular virial `sum_(i<j) d_ij . F_ij` over the pairs' centres,
    /// `d_ij` the nearest-image separation and `F_ij` the force on `i` from
    /// `j`: with the translational kinetic energy it gives the pressure
    /// (see [`Liquid::pressure`]).
    pub virial: f64,
    /// The induced dipoles at the end of the evaluation, per molecule per
    /// polarisable site, to start the next from (empty without induction).
    pub dipoles: Vec<Vec<Vec3>>,
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
        self.forces_from(law, None)
    }

    /// [`Liquid::forces`], the induced dipoles (if the law has them) started
    /// from `warm`: the last step's, which the next step's differ from by a
    /// little, so the solve takes a few iterations and not twenty.
    pub fn forces_from(&self, law: &dyn SiteLaw, warm: Option<&[Vec<Vec3>]>) -> Forces {
        let n = self.com.len();
        let arms: Vec<Vec<Vec3>> = (0..n).map(|i| self.arms(i)).collect();
        let types = &self.kind.types;
        let job = |idx: &mut dyn Iterator<Item = usize>| {
            idx.map(|i| {
                let (mut e, mut f, mut t, mut w) = (0.0, Vec3::ZERO, Vec3::ZERO, 0.0);
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
                    let fij = fi.scale(s) - d.scale(u * ds / dist);
                    f += fij;
                    t += ti.scale(s);
                    w += d.dot(fij);
                }
                (i, e, f, t, w)
            }).collect::<Vec<_>>()
        };
        let mut energy_rows = vec![0.0; n];
        let mut virial_rows = vec![0.0; n];
        let mut force = vec![Vec3::ZERO; n];
        let mut torque = vec![Vec3::ZERO; n];
        for part in crate::electrons::scf::parallel_interleaved(n, &job) {
            for (i, e, f, t, w) in part {
                energy_rows[i] = e;
                force[i] = f;
                torque[i] = t;
                virial_rows[i] = w;
            }
        }
        // Every pair is in two rows, with the same product `d . F`: half the
        // rows' sum is the sum over pairs.
        let mut forces = Forces { energy: energy_rows.iter().sum(), force, torque, virial: 0.5 * virial_rows.iter().sum::<f64>(), dipoles: Vec::new() };
        if let Some((alpha, charge)) = law.induction() {
            self.add_induction(&mut forces, &arms, alpha, charge, warm);
        }
        forces
    }

    /// What the molecules do to each other's charges beyond the pair terms:
    /// the induced dipoles of the polarisable sites, in the field of every
    /// other molecule's permanent charges and of one another (`induction.rs`),
    /// each pair of molecules within the cut-off carrying the same switch as
    /// the pair terms. Adds the energy, forces, torques and virial.
    fn add_induction(&self, forces: &mut Forces, arms: &[Vec<Vec3>], alpha: &[f64], charge: &[f64], warm: Option<&[Vec<Vec3>]>) {
        use crate::induction::{solve, Charge, Cluster, Link, PolSite};
        let n = self.com.len();
        let types = &self.kind.types;
        let pol_sites: Vec<usize> = (0..types.len()).filter(|&k| alpha[types[k]] > 0.0).collect();
        let charged: Vec<usize> = (0..types.len()).filter(|&k| charge[types[k]] != 0.0).collect();
        let mut cl = Cluster::default();
        for i in 0..n {
            cl.pol.push(pol_sites.iter().map(|&k| PolSite { pos: self.com[i] + arms[i][k], alpha: alpha[types[k]] }).collect());
            cl.charges.push(charged.iter().map(|&k| Charge { pos: self.com[i] + arms[i][k], q: charge[types[k]] }).collect());
        }
        let mut links = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                let raw = self.com[i] - self.com[j];
                let d = self.image(raw);
                let dist = d.norm();
                if dist >= self.r_cut {
                    continue;
                }
                let (w, dw) = switch(dist, self.r_on, self.r_cut);
                links.push(Link { i, j, shift: raw - d, weight: w, dweight: dw, d });
            }
        }
        let out = solve(&cl, &links, warm, 1e-9);
        forces.energy += out.energy;
        forces.virial += out.virial;
        for i in 0..n {
            for (a, &k) in pol_sites.iter().enumerate() {
                let f = out.force_pol[i][a];
                forces.force[i] += f;
                forces.torque[i] += arms[i][k].cross(f);
            }
            for (c, &k) in charged.iter().enumerate() {
                let f = out.force_charge[i][c];
                forces.force[i] += f;
                forces.torque[i] += arms[i][k].cross(f);
            }
        }
        for (k, l) in links.iter().enumerate() {
            forces.force[l.i] += out.link_force[k];
            forces.force[l.j] -= out.link_force[k];
        }
        forces.dipoles = out.dipoles;
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
        let next = self.forces_from(law, if forces.dipoles.is_empty() { None } else { Some(&forces.dipoles) });
        self.kick(&next, 0.5 * h);
        next
    }

    /// A stream for this run's thermal noise, seeded from a number so a run
    /// replays exactly.
    pub fn noise(seed: u64) -> Stream {
        Stream::at(seed, 0, 0, Purpose::ThermalNoise)
    }
}

/// Tang and Toennies' damping of a dispersion term `C_n / r^n`:
/// `f_n(x) = 1 - e^-x sum_(k=0..n) x^k / k!`, which takes it smoothly to zero
/// where the clouds overlap. Returns `f_n` and its derivative in `x`.
pub fn tang_toennies(n: usize, x: f64) -> (f64, f64) {
    let (mut sum, mut term) = (1.0, 1.0);
    for k in 1..=n {
        term *= x / k as f64;
        sum += term;
    }
    let e = (-x).exp();
    // d/dx [1 - e^-x S_n] = e^-x (S_n - S_(n-1)) = e^-x x^n / n!.
    (1.0 - e * sum, e * term)
}

/// A site-site law with induced dipoles added (`induction.rs`): the law's pair
/// terms unchanged, and each site type's polarisability, so that the liquid
/// and the fit add what the molecules do to each other's charges. The atoms
/// that polarise are the types whose `alpha` is not zero.
pub struct Polarisable<'a> {
    pub law: &'a SiteSite,
    pub alpha: Vec<f64>,
}

impl SiteLaw for Polarisable<'_> {
    fn site_pair(&self, type_a: usize, type_b: usize, r: f64) -> (f64, f64) {
        self.law.site_pair(type_a, type_b, r)
    }

    fn induction(&self) -> Option<(&[f64], &[f64])> {
        Some((&self.alpha, &self.law.charge))
    }
}

/// Route A's site-site law (PLAY.md E8): between sites of types `a` and `b`,
/// `q_a q_b / r + A_ab exp(-B_ab r) - f_6(B_ab r) C6_ab / r^6 - f_8(B_ab r) C8_ab / r^8`
/// — electrostatics, exchange repulsion and damped dispersion, each a
/// physical mechanism; its numbers are fitted only to the engine's own pair
/// energies.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteSite {
    /// A charge per site type.
    pub charge: Vec<f64>,
    /// Per pair of types, `[A, B, C6, C8]`, indexed `a * types + b`, symmetric.
    pub pair: Vec<[f64; 4]>,
    /// The width of each type's charge, bohr: a Gaussian cloud and not a point,
    /// so two charges closer than their clouds feel `erf(r / s) / r` with
    /// `s^2 = sigma_a^2 + sigma_b^2` (charge penetration: a charge inside
    /// another's cloud sees less of it). Empty, or zero, is a point charge.
    /// Only the permanent electrostatics are smeared; the induced dipoles'
    /// field still comes from the point charges.
    pub sigma: Vec<f64>,
}

impl SiteLaw for SiteSite {
    fn site_pair(&self, ta: usize, tb: usize, r: f64) -> (f64, f64) {
        let types = self.charge.len();
        let [a, b, c6, c8] = self.pair[ta * types + tb];
        let qq = self.charge[ta] * self.charge[tb];
        let s2 = match (self.sigma.get(ta), self.sigma.get(tb)) {
            (Some(x), Some(y)) => x * x + y * y,
            _ => 0.0,
        };
        // qq times 1/r, smeared: the value and its derivative with respect to r.
        let (coul, dcoul) = if s2 > 0.0 && qq != 0.0 {
            let sg = s2.sqrt();
            let g = crate::induction::erf(r / sg);
            let c = 2.0 / (std::f64::consts::PI.sqrt() * sg) * (-r * r / s2).exp();
            (qq * g / r, qq * (c / r - g / (r * r)))
        } else {
            (qq / r, -qq / (r * r))
        };
        let rep = a * (-b * r).exp();
        let (f6, df6) = tang_toennies(6, b * r);
        let (f8, df8) = tang_toennies(8, b * r);
        let (r6, r8) = (r.powi(6), r.powi(8));
        let u = coul + rep - f6 * c6 / r6 - f8 * c8 / r8;
        let du = dcoul - b * rep - (b * df6 * c6 / r6 - 6.0 * f6 * c6 / (r6 * r)) - (b * df8 * c8 / r8 - 8.0 * f8 * c8 / (r8 * r));
        (u, du)
    }
}

impl Liquid {
    /// Molecules at the box-frame density profile's centre, as bins of
    /// molecules per unit volume along `z` from `-cell_z / 2` to `cell_z / 2`
    /// about the slab's centre of mass — found as a circular mean, since the
    /// box is periodic and the slab wanders.
    pub fn density_profile(&self, bins: usize) -> Vec<f64> {
        let lz = self.cell[2];
        let (mut c, mut s) = (0.0, 0.0);
        for r in &self.com {
            let th = 2.0 * std::f64::consts::PI * r.z / lz;
            c += th.cos();
            s += th.sin();
        }
        let centre = s.atan2(c) / (2.0 * std::f64::consts::PI) * lz;
        let mut hist = vec![0.0; bins];
        let vol = self.cell[0] * self.cell[1] * lz / bins as f64;
        for r in &self.com {
            let mut z = r.z - centre;
            z -= lz * (z / lz).round();
            let k = (((z / lz) + 0.5) * bins as f64).floor().clamp(0.0, (bins - 1) as f64) as usize;
            hist[k] += 1.0 / vol;
        }
        hist
    }
}

/// One computed pair energy for fitting a law: the two molecules' sites (box
/// frame, bohr) with their types, and the interaction energy (hartree).
#[derive(Debug, Clone)]
pub struct PairEnergy {
    pub a: Vec<(Vec3, usize)>,
    pub b: Vec<(Vec3, usize)>,
    pub energy: f64,
}

/// A law's energy for one pair configuration.
pub fn pair_energy(law: &dyn SiteLaw, pair: &PairEnergy) -> f64 {
    let mut e = 0.0;
    for (pa, ta) in &pair.a {
        for (pb, tb) in &pair.b {
            e += law.site_pair(*ta, *tb, (*pa - *pb).norm()).0;
        }
    }
    if let Some((alpha, charge)) = law.induction() {
        e += induction_pair_energy(alpha, charge, pair);
    }
    e
}

/// The induction energy of two molecules alone: their polarisable sites'
/// dipoles in each other's permanent charges, from `induction::solve` with one
/// link of weight one. Nothing for a molecule by itself (no field, no dipoles),
/// so this is the whole of induction's part of the pair's interaction energy.
pub fn induction_pair_energy(alpha: &[f64], charge: &[f64], pair: &PairEnergy) -> f64 {
    use crate::induction::{solve, Charge, Cluster, Link, PolSite};
    let mut cl = Cluster::default();
    for m in [&pair.a, &pair.b] {
        cl.pol.push(m.iter().filter(|(_, t)| alpha[*t] > 0.0).map(|(p, t)| PolSite { pos: *p, alpha: alpha[*t] }).collect());
        cl.charges.push(m.iter().filter(|(_, t)| charge[*t] != 0.0).map(|(p, t)| Charge { pos: *p, q: charge[*t] }).collect());
    }
    let links = [Link { i: 0, j: 1, shift: Vec3::ZERO, weight: 1.0, dweight: 0.0, d: Vec3::ZERO }];
    solve(&cl, &links, None, 1e-12).energy
}

/// What a fit produced: the law, and its weighted and unweighted root-mean-
/// square residuals over the data it was fitted to.
pub struct Fitted {
    pub law: SiteSite,
    pub weighted_rms: f64,
    pub rms: f64,
    pub iterations: usize,
}

/// Fit route A's site-site law to computed pair energies by Levenberg-
/// Marquardt, every number from the energies alone (PLAY.md E8): charges per
/// type with the molecule kept neutral (`multiplicity[t]` sites of type `t`
/// in one molecule), and per pair of types `A, B, C6, C8`, each fitted
/// through its logarithm so it stays positive. Each configuration weighs
/// `exp(-(E - E_min) / (k T))`, held at least at `floor`, so the fit spends
/// itself where the liquid goes while the repulsive wall still counts.
pub fn fit_site_site(data: &[PairEnergy], multiplicity: &[usize], start: &SiteSite, temperature: f64, floor: f64, max_iterations: usize) -> Fitted {
    fit_site_site_held(data, multiplicity, start, temperature, floor, max_iterations, &Held::default())
}

/// A molecule's sites with one more added: of type `site_type`, `distance`
/// bohr from the first site along the bisector of its bonds to the second and
/// third (the lone pairs' charge of a water, which `Kind::of_molecule` places
/// the same way). For laws fitted and read with a `bisector` line.
pub fn with_bisector_site(molecule: &[(Vec3, usize)], site_type: usize, distance: f64) -> Vec<(Vec3, usize)> {
    let o = molecule[0].0;
    let (u, v) = ((molecule[1].0 - o).unit(), (molecule[2].0 - o).unit());
    let mut out = molecule.to_vec();
    out.push((o + (u + v).unit().scale(distance), site_type));
    out
}

/// Numbers of a law a fit leaves at the values it was started with.
///
/// A site that carries charge and nothing else (an off-atom site standing in
/// for a lone pair's charge) has no repulsion and no dispersion of its own:
/// its pairs are held at zero. A dispersion coefficient the engine derived
/// elsewhere (E6's partition of the E7 kernel) is held rather than refitted,
/// when the question is what the energies say about everything else.
#[derive(Default)]
pub struct Held<'a> {
    /// Type pairs whose `A, B, C6, C8` are all held.
    pub pairs: &'a [(usize, usize)],
    /// Type pairs whose `C6` alone is held.
    pub dispersion: &'a [(usize, usize)],
    /// Fit a charge width for every type (see [`SiteSite::sigma`]).
    pub sigma: bool,
    /// Keep the starting charges (the last still follows from neutrality):
    /// charges derived from the molecule's own density (`phys-esp`), so that
    /// the pair energies fit only what the density does not say.
    pub charges: bool,
}

/// [`fit_site_site`] with some numbers held at their starting values.
pub fn fit_site_site_held(data: &[PairEnergy], multiplicity: &[usize], start: &SiteSite, temperature: f64, floor: f64, max_iterations: usize, held: &Held) -> Fitted {
    fit_site_site_polarised(data, multiplicity, start, temperature, floor, max_iterations, held, None)
}

/// [`fit_site_site_held`] with induced dipoles: `alpha` is each site type's
/// polarisability, derived elsewhere (`phys-polar`) and not fitted, and every
/// pair's energy in the fit is the pair law's plus the induction of the two
/// molecules (`induction_pair_energy`) in the law's own charges as they are
/// fitted. The law that comes out is the pair law; it is a law with induction
/// only together with `alpha`.
pub fn fit_site_site_polarised(data: &[PairEnergy], multiplicity: &[usize], start: &SiteSite, temperature: f64, floor: f64, max_iterations: usize, held: &Held, alpha: Option<&[f64]>) -> Fitted {
    let types = start.charge.len();
    let pairs: Vec<(usize, usize)> = (0..types).flat_map(|a| (a..types).map(move |b| (a, b))).collect();
    // Parameters: charges of types 0..types-1 (the last is set by
    // neutrality), then per unique pair ln A, ln B, ln C6, ln C8.
    let pack = |law: &SiteSite| -> Vec<f64> {
        let mut x: Vec<f64> = law.charge[..types - 1].to_vec();
        for &(a, b) in &pairs {
            for v in law.pair[a * types + b] {
                x.push(v.max(1e-300).ln());
            }
        }
        if held.sigma {
            for t in 0..types {
                x.push(law.sigma.get(t).copied().filter(|v| *v > 0.0).unwrap_or(0.5).ln());
            }
        }
        x
    };
    let unpack = |x: &[f64]| -> SiteSite {
        let mut charge: Vec<f64> = x[..types - 1].to_vec();
        let partial: f64 = charge.iter().zip(multiplicity).map(|(q, m)| q * *m as f64).sum();
        charge.push(-partial / multiplicity[types - 1] as f64);
        let mut pair = vec![[0.0; 4]; types * types];
        for (k, &(a, b)) in pairs.iter().enumerate() {
            let o = types - 1 + 4 * k;
            let v = [x[o].exp(), x[o + 1].exp(), x[o + 2].exp(), x[o + 3].exp()];
            pair[a * types + b] = v;
            pair[b * types + a] = v;
        }
        let sigma = if held.sigma { (0..types).map(|t| x[types - 1 + 4 * pairs.len() + t].exp()).collect() } else { start.sigma.clone() };
        SiteSite { charge, pair, sigma }
    };
    // What is fitted: every parameter but the held ones, which keep what
    // `start` gave them.
    let base = pack(start);
    let is_held = |list: &[(usize, usize)], a: usize, b: usize| list.iter().any(|&(p, q)| (p, q) == (a, b) || (q, p) == (a, b));
    let mut free: Vec<usize> = if held.charges { Vec::new() } else { (0..types - 1).collect() };
    for (k, &(a, b)) in pairs.iter().enumerate() {
        let o = types - 1 + 4 * k;
        if is_held(held.pairs, a, b) {
            continue;
        }
        for j in 0..4 {
            if !(j == 2 && is_held(held.dispersion, a, b)) {
                free.push(o + j);
            }
        }
    }
    if held.sigma {
        free.extend((0..types).map(|t| types - 1 + 4 * pairs.len() + t));
    }
    let expand = |xf: &[f64]| -> Vec<f64> {
        let mut full = base.clone();
        for (j, &i) in free.iter().enumerate() {
            full[i] = xf[j];
        }
        full
    };
    let e_min = data.iter().map(|d| d.energy).fold(f64::INFINITY, f64::min);
    let kt = K_B * temperature;
    let weights: Vec<f64> = data.iter().map(|d| (-(d.energy - e_min) / kt).exp().max(floor)).collect();
    let energy_of = |law: &SiteSite, d: &PairEnergy| -> f64 {
        match alpha {
            Some(a) => pair_energy(&Polarisable { law, alpha: a.to_vec() }, d),
            None => pair_energy(law, d),
        }
    };
    let residuals = |x: &[f64]| -> Vec<f64> {
        let law = unpack(&expand(x));
        data.iter().zip(&weights).map(|(d, w)| w.sqrt() * (energy_of(&law, d) - d.energy)).collect()
    };
    let cost = |r: &[f64]| r.iter().map(|v| v * v).sum::<f64>();
    let mut x: Vec<f64> = free.iter().map(|&i| base[i]).collect();
    let np = x.len();
    let mut r = residuals(&x);
    let mut c = cost(&r);
    let mut lambda = 1e-3;
    let mut iterations = 0;
    for it in 0..max_iterations {
        iterations = it + 1;
        // Jacobian by central differences: a handful of parameters, cheap.
        let mut jac = vec![vec![0.0; np]; r.len()];
        for p in 0..np {
            let h = 1e-6 * x[p].abs().max(1e-3);
            let (mut xp, mut xm) = (x.clone(), x.clone());
            xp[p] += h;
            xm[p] -= h;
            let (rp, rm) = (residuals(&xp), residuals(&xm));
            for i in 0..r.len() {
                jac[i][p] = (rp[i] - rm[i]) / (2.0 * h);
            }
        }
        let mut jtj = vec![vec![0.0; np]; np];
        let mut jtr = vec![0.0; np];
        for i in 0..r.len() {
            for p in 0..np {
                jtr[p] += jac[i][p] * r[i];
                for q in 0..np {
                    jtj[p][q] += jac[i][p] * jac[i][q];
                }
            }
        }
        let mut improved = false;
        for _ in 0..30 {
            let mut m = crate::electrons::linalg::Matrix::zeros(np);
            for p in 0..np {
                for q in 0..np {
                    m.set(p, q, jtj[p][q] + if p == q { lambda * jtj[p][p].max(1e-30) } else { 0.0 });
                }
            }
            let Some(step) = solve_dense(&m, &jtr) else {
                lambda *= 10.0;
                continue;
            };
            let trial: Vec<f64> = x.iter().zip(&step).map(|(a, s)| a - s).collect();
            let rt = residuals(&trial);
            let ct = cost(&rt);
            if ct.is_finite() && ct < c {
                x = trial;
                r = rt;
                let gain = (c - ct) / c.max(1e-300);
                c = ct;
                lambda = (lambda / 3.0).max(1e-12);
                improved = true;
                if gain < 1e-10 {
                    lambda = f64::INFINITY;
                }
                break;
            }
            lambda *= 4.0;
        }
        if !improved || !lambda.is_finite() {
            break;
        }
    }
    let law = unpack(&expand(&x));
    let wsum: f64 = weights.iter().sum();
    let weighted_rms = (c / wsum).sqrt();
    let rms = (data.iter().map(|d| (energy_of(&law, d) - d.energy).powi(2)).sum::<f64>() / data.len() as f64).sqrt();
    Fitted { law, weighted_rms, rms, iterations }
}

/// Solve a small dense symmetric positive system by Cholesky; `None` if it
/// is not positive definite.
fn solve_dense(m: &crate::electrons::linalg::Matrix, b: &[f64]) -> Option<Vec<f64>> {
    let n = m.n;
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = m.get(i, j);
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i * n + i] = s.sqrt();
            } else {
                l[i * n + j] = s / l[j * n + j];
            }
        }
    }
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i * n + k] * y[k];
        }
        y[i] = s / l[i * n + i];
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in i + 1..n {
            s -= l[k * n + i] * x[k];
        }
        x[i] = s / l[i * n + i];
    }
    Some(x)
}

impl Liquid {
    /// A slab of `n` molecules of `kind` at `density` (molecules per bohr^3),
    /// on a simple cubic lattice filling the middle of a box `side x side x
    /// length`, with empty space either side for the vapour; orientations and
    /// velocities drawn from `seed` at `temperature`, the total momentum zero.
    pub fn slab(kind: Kind, n: usize, density: f64, side: f64, length: f64, temperature: f64, seed: u64, r_on: f64, r_cut: f64) -> Liquid {
        // As near cubic as the box's side allows: the spacing the density
        // gives, rounded to fit the side, and the layers spaced to keep the
        // density.
        let across = ((side * density.cbrt()).round() as usize).max(1);
        let spacing = side / across as f64;
        let layers = n.div_ceil(across * across);
        let dz = 1.0 / (density * spacing * spacing);
        let z0 = -0.5 * layers as f64 * dz;
        let mut s = Stream::at(seed, 0, 0, Purpose::Positions);
        let mut v = Stream::at(seed, 0, 0, Purpose::Velocities);
        let (mut com, mut orientation, mut momentum, mut spin) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        'fill: for k in 0..layers {
            for i in 0..across {
                for j in 0..across {
                    if com.len() == n {
                        break 'fill;
                    }
                    com.push(Vec3 { x: (i as f64 + 0.5) * spacing, y: (j as f64 + 0.5) * spacing, z: z0 + (k as f64 + 0.5) * dz });
                    let (u1, u2, u3) = (s.uniform(), s.uniform(), s.uniform());
                    let tau = std::f64::consts::TAU;
                    let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
                    orientation.push(Quat { w: a * (tau * u2).sin(), v: Vec3 { x: a * (tau * u2).cos(), y: b * (tau * u3).sin(), z: b * (tau * u3).cos() } });
                    momentum.push(v.normal3().scale((kind.mass * K_B * temperature).sqrt()));
                    let i3 = kind.inertia;
                    spin.push(Vec3 { x: v.normal() * (i3[0] * K_B * temperature).sqrt(), y: v.normal() * (i3[1] * K_B * temperature).sqrt(), z: v.normal() * (i3[2] * K_B * temperature).sqrt() });
                }
            }
        }
        let mean = momentum.iter().fold(Vec3::ZERO, |acc, p| acc + *p).scale(1.0 / momentum.len() as f64);
        for p in momentum.iter_mut() {
            *p -= mean;
        }
        Liquid { kind, cell: [side, side, length], com, orientation, momentum, spin, r_on, r_cut }
    }
}

/// Liquid and vapour densities from a profile about the slab's centre
/// (`Liquid::density_profile`): the mean over the middle `core` fraction of
/// the bins, and over the `far` fraction farthest from the centre on both
/// sides.
pub fn coexisting_densities(profile: &[f64], core: f64, far: f64) -> (f64, f64) {
    let n = profile.len();
    let mid = n / 2;
    let half_core = ((core * n as f64) / 2.0).round().max(1.0) as usize;
    let liquid = profile[mid - half_core..mid + half_core].iter().sum::<f64>() / (2 * half_core) as f64;
    let edge = ((far * n as f64) / 2.0).round().max(1.0) as usize;
    let vapour = (profile[..edge].iter().sum::<f64>() + profile[n - edge..].iter().sum::<f64>()) / (2 * edge) as f64;
    (liquid, vapour)
}

impl SiteSite {
    /// As text: a line of charges, then a line `a b A B C6 C8` for each pair
    /// of types with `a <= b`, every number to the last digit.
    pub fn to_text(&self) -> String {
        let types = self.charge.len();
        let mut s = format!("charges {}\n", self.charge.iter().map(|q| format!("{q:e}")).collect::<Vec<_>>().join(" "));
        for a in 0..types {
            for b in a..types {
                let p = self.pair[a * types + b];
                s += &format!("pair {a} {b} {:e} {:e} {:e} {:e}\n", p[0], p[1], p[2], p[3]);
            }
        }
        s
    }

    /// The `alpha <type> <value>` lines of a law's text: each site type's
    /// polarisability in bohr cubed (0 for a type that has no line), for the
    /// law's induced dipoles. Derived by `phys-polar` from the molecule's own
    /// electrons and added to the law's text by hand or script.
    pub fn alpha_from_text(text: &str, types: usize) -> Vec<f64> {
        let mut alpha = vec![0.0; types];
        for line in text.lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            if w.first() == Some(&"alpha") && w.len() == 3 {
                if let (Ok(t), Ok(a)) = (w[1].parse::<usize>(), w[2].parse::<f64>()) {
                    if t < types {
                        alpha[t] = a;
                    }
                }
            }
        }
        alpha
    }

    /// Read back what [`SiteSite::to_text`] wrote.
    pub fn from_text(text: &str) -> Option<SiteSite> {
        let mut charge = Vec::new();
        let mut sigma = Vec::new();
        let mut pairs = Vec::new();
        for line in text.lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            match w.first().copied() {
                Some("charges") => charge = w[1..].iter().map(|x| x.parse().ok()).collect::<Option<Vec<f64>>>()?,
                // Read by `bisector_from_text`: where an off-atom site goes, not
                // part of the pair law.
                Some("bisector") => {}
                // Read by `alpha_from_text`: the polarisability of a site type.
                Some("alpha") => {}
                Some("sigma") => sigma = w[1..].iter().map(|x| x.parse().ok()).collect::<Option<Vec<f64>>>()?,
                Some("pair") if w.len() == 7 => {
                    let (a, b): (usize, usize) = (w[1].parse().ok()?, w[2].parse().ok()?);
                    let v: Vec<f64> = w[3..].iter().map(|x| x.parse().ok()).collect::<Option<Vec<f64>>>()?;
                    pairs.push((a, b, [v[0], v[1], v[2], v[3]]));
                }
                None => {}
                _ => return None,
            }
        }
        let types = charge.len();
        let mut pair = vec![[0.0; 4]; types * types];
        for (a, b, v) in pairs {
            pair[a * types + b] = v;
            pair[b * types + a] = v;
        }
        Some(SiteSite { charge, pair, sigma })
    }
}

impl Liquid {
    /// The pressure, hartree per cubic bohr: `(2 K_translation + W) / (3 V)`
    /// with `W` the molecular virial. For rigid molecules the centres' forces
    /// and kinetic energy are all there is to it, since the constraints that
    /// hold a molecule together do no work between molecules.
    pub fn pressure(&self, forces: &Forces) -> f64 {
        let volume = self.cell[0] * self.cell[1] * self.cell[2];
        let (translation, _) = self.kinetic();
        (2.0 * translation + forces.virial) / (3.0 * volume)
    }
}

/// One hartree per cubic bohr in bar.
pub const BAR_PER_HARTREE_PER_BOHR3: f64 = 2.942_101_57e8;

impl SiteSite {
    /// The `bisector <type> <distance>` line of a law's text, if it has one: a
    /// site of that type, carrying electrostatics only, at that distance (bohr)
    /// from the first atom along the bisector of its bonds to the second and
    /// third — the TIP4P-like charge site the water fit found it needed.
    pub fn bisector_from_text(text: &str) -> Option<(usize, f64)> {
        text.lines().find_map(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            (w.first() == Some(&"bisector") && w.len() == 3).then(|| Some((w[1].parse().ok()?, w[2].parse().ok()?))).flatten()
        })
    }
}

impl Kind {
    /// A water-like molecule's kind from its element numbers, positions (bohr)
    /// and site types, atoms first: masses from the elements, and, if asked,
    /// one more massless site of the given type on the bisector of atom 0's
    /// bonds to atoms 1 and 2, at the given distance from atom 0. A massless
    /// site changes neither the centre of mass nor the inertia, and rides in
    /// the principal frame with the rest.
    pub fn of_molecule(z: &[u32], positions: &[[f64; 3]], types: &[usize], bisector: Option<(usize, f64)>) -> Kind {
        let mut pos = positions.to_vec();
        let mut ty = types.to_vec();
        let mut masses: Vec<f64> = z.iter().map(|&zz| crate::chem::elements::Element(zz as u8).mass_kg().expect("a mass") / 1.66053906660e-27 * AMU).collect();
        if let Some((site_type, d)) = bisector {
            let at = |k: usize| Vec3 { x: positions[k][0], y: positions[k][1], z: positions[k][2] };
            let (o, u, v) = (at(0), (at(1) - at(0)).unit(), (at(2) - at(0)).unit());
            let m = o + (u + v).unit().scale(d);
            pos.push([m.x, m.y, m.z]);
            ty.push(site_type);
            masses.push(0.0);
        }
        Kind::from_atoms(&pos, &masses, &ty)
    }
}
