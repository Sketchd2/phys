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

    /// The width of each site type's charge cloud, for the field the charges
    /// make on the induced dipoles (empty: points).
    fn charge_widths(&self) -> &[f64] {
        &[]
    }

    /// The site-site law underneath, where the law is (or wraps) one: a
    /// molecule whose charges and widths follow its geometry (`geolaw.rs`)
    /// needs the table, not just the answers.
    fn site_site(&self) -> Option<&SiteSite> {
        None
    }

    /// The Ewald sum that carries the long-range part of the charges and the
    /// induced dipoles (`ewald.rs`), where the law is wrapped in [`EwaldLaw`]:
    /// the pair terms then have only their short-range Coulomb remainder.
    fn ewald(&self) -> Option<&crate::ewald::Ewald> {
        None
    }
}

/// Everything the site-site pair loop of a flexible liquid needs, for an engine
/// that does it elsewhere (the GPU): the molecules' sites relative to their
/// centres, each site's charge and Coulomb width, and the pairs of molecules
/// inside the cut-off with the centres' separation and the switch
/// (`induction::Link`: `weight` is the switch, `dweight` its slope).
pub struct PairProblem<'a> {
    pub table: &'a SiteSite,
    /// The site type of each of a molecule's sites.
    pub types: &'a [usize],
    /// Per molecule, each site's place relative to the molecule's centre.
    pub rel: &'a [Vec<Vec3>],
    pub charge: &'a [Vec<f64>],
    pub sigma: &'a [Vec<f64>],
    pub links: &'a [crate::induction::Link],
    /// The Ewald splitting width, when the law has only the short-range remainder of its Coulomb
    /// term: the pair loop takes `erf(r / s') / r` off each charged pair (`EwaldLaw::site_pair`).
    pub ewald_width: Option<f64>,
}

/// What the pair loop gives back.
pub struct PairResult {
    pub energy: f64,
    pub virial: f64,
    /// Per molecule, the force on each site.
    pub fsite: Vec<Vec<Vec3>>,
    /// Per molecule, the switch's force on its centre.
    pub fcentre: Vec<Vec3>,
    /// Per molecule, per site, the energy's derivative with respect to the
    /// site's charge and to its Coulomb width.
    pub gq: Vec<Vec<f64>>,
    pub gs: Vec<Vec<f64>>,
}

/// A way of doing the pair loop faster than `Flex::forces` does, to the same
/// answer to round-off (its sums are added in another, fixed, order).
pub trait PairEngine: Send + Sync {
    fn pairs(&self, problem: &PairProblem) -> Option<PairResult>;
}

static PAIR_ENGINE: std::sync::Mutex<Option<std::sync::Arc<dyn PairEngine>>> = std::sync::Mutex::new(None);

/// Install a pair-loop engine for this process (`None` removes it).
pub fn set_pair_engine(engine: Option<std::sync::Arc<dyn PairEngine>>) {
    *PAIR_ENGINE.lock().expect("the pair engine") = engine;
}

/// The installed pair-loop engine, if any.
pub fn pair_engine() -> Option<std::sync::Arc<dyn PairEngine>> {
    PAIR_ENGINE.lock().expect("the pair engine").clone()
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

/// Molecules in a periodic box: of one kind, or of several (a solution), the law then a table whose
/// site types are disjoint between the kinds (`docs/SOLUTIONS.md` 4.1: a mixture's law is the pure
/// components' laws plus the cross entries of one table). `kind` is the first species; `other_kinds`
/// are the rest, and `species[i]` says which kind molecule `i` is (empty: all the first).
#[derive(Debug, Clone)]
pub struct Liquid {
    pub kind: Kind,
    pub other_kinds: Vec<Kind>,
    pub species: Vec<usize>,
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

    /// Which kind molecule `i` is: 0 for the first.
    pub fn species_of(&self, i: usize) -> usize {
        self.species.get(i).copied().unwrap_or(0)
    }

    /// The kind of molecule `i`.
    pub fn kind_of(&self, i: usize) -> &Kind {
        match self.species_of(i) {
            0 => &self.kind,
            s => &self.other_kinds[s - 1],
        }
    }

    /// Every kind in the box, the first first.
    pub fn all_kinds(&self) -> Vec<&Kind> {
        std::iter::once(&self.kind).chain(self.other_kinds.iter()).collect()
    }

    /// Whether every molecule is of the first kind.
    pub fn homogeneous(&self) -> bool {
        self.species.iter().all(|&s| s == 0)
    }

    /// Every site of molecule `i`, in the box frame, relative to its centre.
    fn arms(&self, i: usize) -> Vec<Vec3> {
        self.kind_of(i).sites.iter().map(|s| self.orientation[i].rotate(*s)).collect()
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
        let t_forces = std::time::Instant::now();
        let n = self.com.len();
        let arms: Vec<Vec<Vec3>> = (0..n).map(|i| self.arms(i)).collect();
        let types = &self.kind.types;
        // The pair loop on an installed engine (the GPU), if there is one and it takes this law and
        // the molecules are alike (the engines' layouts are for identical molecules).
        let engine_result = crate::liquid::pair_engine().filter(|_| law.site_site().is_some() && self.homogeneous()).and_then(|e| {
            let tab = law.site_site()?;
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
                    links.push(crate::induction::Link { i, j, shift: raw - d, weight: w, dweight: dw, d });
                }
            }
            let qs: Vec<Vec<f64>> = (0..n).map(|_| types.iter().map(|&t| tab.charge[t]).collect()).collect();
            let ss: Vec<Vec<f64>> = (0..n).map(|_| types.iter().map(|&t| tab.sigma.get(t).copied().unwrap_or(0.0)).collect()).collect();
            e.pairs(&PairProblem { table: tab, types, rel: &arms, charge: &qs, sigma: &ss, links: &links, ewald_width: law.ewald().map(|w| w.s) })
        });
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
                    for (a, (pa, ta)) in arms[i].iter().zip(&self.kind_of(i).types).enumerate() {
                        let _ = a;
                        for (pb, tb) in arms[j].iter().zip(&self.kind_of(j).types) {
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
        let mut forces = match engine_result {
            Some(r) => {
                let force: Vec<Vec3> = (0..n).map(|i| r.fsite[i].iter().fold(r.fcentre[i], |acc, f| acc + *f)).collect();
                let torque: Vec<Vec3> = (0..n).map(|i| arms[i].iter().zip(&r.fsite[i]).fold(Vec3::ZERO, |acc, (p, f)| acc + p.cross(*f))).collect();
                Forces { energy: r.energy, force, torque, virial: r.virial, dipoles: Vec::new() }
            }
            None => {
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
                Forces { energy: energy_rows.iter().sum(), force, torque, virial: 0.5 * virial_rows.iter().sum::<f64>(), dipoles: Vec::new() }
            }
        };
        let t_pairs = t_forces.elapsed();
        if let Some(ew) = law.ewald() {
            self.add_ewald(&mut forces, &arms, law, ew, warm);
            if std::env::var_os("PHYS_PROFILE").is_some() {
                eprintln!("  forces: pair loop {:.1} ms, whole {:.1} ms", t_pairs.as_secs_f64() * 1e3, t_forces.elapsed().as_secs_f64() * 1e3);
            }
        } else if let Some((alpha, charge)) = law.induction() {
            self.add_induction(&mut forces, &arms, alpha, charge, law.charge_widths(), warm);
        }
        forces
    }

    /// What the molecules do to each other's charges beyond the pair terms:
    /// the induced dipoles of the polarisable sites, in the field of every
    /// other molecule's permanent charges and of one another (`induction.rs`),
    /// each pair of molecules within the cut-off carrying the same switch as
    /// the pair terms. Adds the energy, forces, torques and virial.
    fn add_induction(&self, forces: &mut Forces, arms: &[Vec<Vec3>], alpha: &[f64], charge: &[f64], widths: &[f64], warm: Option<&[Vec<Vec3>]>) {
        use crate::induction::{solve, Charge, Cluster, Link, PolSite};
        let n = self.com.len();
        let kinds = self.all_kinds();
        let pol_of: Vec<Vec<usize>> = kinds.iter().map(|kd| (0..kd.types.len()).filter(|&k| alpha[kd.types[k]] > 0.0).collect()).collect();
        let ch_of: Vec<Vec<usize>> = kinds.iter().map(|kd| (0..kd.types.len()).filter(|&k| charge[kd.types[k]] != 0.0).collect()).collect();
        let mut cl = Cluster::default();
        for i in 0..n {
            let s = self.species_of(i);
            let types = &kinds[s].types;
            cl.pol.push(pol_of[s].iter().map(|&k| PolSite { pos: self.com[i] + arms[i][k], alpha: alpha[types[k]] }).collect());
            cl.charges.push(ch_of[s].iter().map(|&k| Charge { pos: self.com[i] + arms[i][k], q: charge[types[k]], sigma: widths.get(types[k]).copied().unwrap_or(0.0) }).collect());
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
            let s = self.species_of(i);
            for (a, &k) in pol_of[s].iter().enumerate() {
                let f = out.force_pol[i][a];
                forces.force[i] += f;
                forces.torque[i] += arms[i][k].cross(f);
            }
            for (c, &k) in ch_of[s].iter().enumerate() {
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

    /// The long-range part of the permanent charges' interactions and of the
    /// induced dipoles' (`ewald.rs`, `PLAY.md` E8d), on top of the pair loop's
    /// short-range remainder: the reciprocal sums, the dipoles solved with
    /// the real-space remainder and the reciprocal coupling together
    /// (`induction::solve_ewald`), the reciprocal forces on every charge and
    /// dipole site, the constant self and same-molecule terms taken off the
    /// energy, and the reciprocal virial. Adds to `forces`. Molecules of several kinds each
    /// have their own charged sites, merged charges, polarisable sites and constant.
    fn add_ewald(&self, forces: &mut Forces, arms: &[Vec<Vec3>], law: &dyn SiteLaw, ew: &crate::ewald::Ewald, warm: Option<&[Vec<Vec3>]>) {
        use crate::induction::{solve_ewald, Charge, Cluster, Link, PolSite};
        let n = self.com.len();
        let kinds = self.all_kinds();
        let ss = law.site_site().expect("the Ewald sum needs a site-site law underneath");
        let induced = law.induction();
        /// What one kind of molecule contributes to the sums.
        struct Species {
            charged: Vec<usize>,
            /// (a member site for the arm, the summed charge)
            groups: Vec<(usize, f64)>,
            pol_sites: Vec<usize>,
            /// The long-range part of the pairs inside one molecule (its own constant).
            intra: f64,
        }
        // The first molecule of each kind, whose arms the kind's coincident charges and constant are read from.
        let first_of: Vec<Option<usize>> = (0..kinds.len()).map(|s| (0..n).find(|&i| self.species_of(i) == s)).collect();
        let ix: Vec<Species> = kinds
            .iter()
            .enumerate()
            .map(|(s, kd)| {
                let types = &kd.types;
                let q_of = |k: usize| ss.charge[types[k]];
                let arm = |k: usize| match first_of[s] {
                    Some(i0) => arms[i0][k],
                    None => kd.sites[k],
                };
                let charged: Vec<usize> = (0..types.len()).filter(|&k| q_of(k) != 0.0).collect();
                // In reciprocal space every charge is a point, so charges at one place (a core and its cloud
                // on the same atom) are one charge there: the structure factor, and with it the energy, the
                // forces on the molecule and the virial, are unchanged, and there are fewer sites to sum.
                // The self term and the same-molecule terms below keep the original charges.
                let mut groups: Vec<(usize, f64)> = Vec::new();
                for &k in &charged {
                    match groups.iter_mut().find(|(k0, _)| (arm(*k0) - arm(k)).norm() < 1e-9) {
                        Some(g) => g.1 += q_of(k),
                        None => groups.push((k, q_of(k))),
                    }
                }
                let pol_sites: Vec<usize> = match induced {
                    Some((alpha, _)) => (0..types.len()).filter(|&k| alpha[types[k]] > 0.0).collect(),
                    None => Vec::new(),
                };
                let mut intra = 0.0;
                for (a, &ka) in charged.iter().enumerate() {
                    for &kb in &charged[a + 1..] {
                        let r = (arm(ka) - arm(kb)).norm();
                        intra += q_of(ka) * q_of(kb) * crate::ewald::erf_over_r(r, ew.s).0;
                    }
                }
                Species { charged, groups, pol_sites, intra }
            })
            .collect();
        let (mut cpos, mut cq, mut ccom, mut carm) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut cq_original: Vec<f64> = Vec::new();
        let (mut dpos, mut dcom, mut darm) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..n {
            let s = self.species_of(i);
            let sp = &ix[s];
            for &(k, q) in &sp.groups {
                cpos.push(self.com[i] + arms[i][k]);
                cq.push(q);
                ccom.push(self.com[i]);
                carm.push(arms[i][k]);
            }
            for &k in &sp.charged {
                cq_original.push(ss.charge[kinds[s].types[k]]);
            }
            for &k in &sp.pol_sites {
                dpos.push(self.com[i] + arms[i][k]);
                dcom.push(self.com[i]);
                darm.push(arms[i][k]);
            }
        }
        let t_a = std::time::Instant::now();
        // The dipoles' cluster and the pairs of molecules within the cut-off.
        let mut cl = Cluster::default();
        let mut links = Vec::new();
        if let Some((alpha, charge)) = induced {
            let widths = law.charge_widths();
            for i in 0..n {
                let s = self.species_of(i);
                let types = &kinds[s].types;
                cl.pol.push(ix[s].pol_sites.iter().map(|&k| PolSite { pos: self.com[i] + arms[i][k], alpha: alpha[types[k]] }).collect());
                cl.charges.push(ix[s].charged.iter().map(|&k| Charge { pos: self.com[i] + arms[i][k], q: charge[types[k]], sigma: widths.get(types[k]).copied().unwrap_or(0.0) }).collect());
            }
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
        }
        // A whole-step engine on the card, if one is installed and takes this liquid; else the CPU's
        // (or the reciprocal engine's) work in its parts.
        let stepped = if induced.is_some() && ix.iter().any(|s| !s.pol_sites.is_empty()) && !links.is_empty() {
            crate::ewald::step_engine().and_then(|e| e.step(ew, &cl, &links, &cpos, &cq, &carm, &dpos, &darm, warm, ind_tolerance()))
        } else {
            None
        };
        let (sq, out, fin) = match stepped {
            Some(s) => (s.sq, if induced.is_some() { Some(s.induced) } else { None }, s.finished),
            None => {
                // The reciprocal work (`ewald.rs`): the CPU's, or the installed engine's on the card.
                let rec = std::cell::RefCell::new(crate::ewald::begin_reciprocal(ew, &cpos, &cq, &carm, &dpos, &darm));
                let sq = rec.borrow_mut().charge_sum();
                let mut out = None;
                let mut mu: Vec<Vec3> = Vec::new();
                if induced.is_some() {
                    let e0_lr = rec.borrow_mut().field_at_dipoles();
                    let t_response = std::cell::Cell::new(std::time::Duration::ZERO);
                    let response = |x: &[Vec3]| {
                        let t0 = std::time::Instant::now();
                        let r = rec.borrow_mut().dipole_response(x);
                        t_response.set(t_response.get() + t0.elapsed());
                        r
                    };
                    let t_b = std::time::Instant::now();
                    let o = crate::induction::solve_ewald(&cl, &links, ew.s, &e0_lr, &response, warm, ind_tolerance());
                    let _ = solve_ewald;
                    if std::env::var_os("PHYS_PROFILE").is_some() {
                        eprintln!("  add_ewald: before solve {:.1} ms; solve {:.1} ms of which reciprocal response {:.1} ms", (t_b - t_a).as_secs_f64() * 1e3, t_b.elapsed().as_secs_f64() * 1e3, t_response.get().as_secs_f64() * 1e3);
                    }
                    mu = o.dipoles.iter().flatten().cloned().collect();
                    out = Some(o);
                }
                let fin = rec.borrow_mut().finish(&mu);
                (sq, out, fin)
            }
        };
        // The permanent charges: the reciprocal energy, less each charge's self
        // term and less the long-range part of the pairs inside a molecule (a
        // rigid molecule's own constant), which the liquid does not have.
        let intra_total = if self.homogeneous() { ix[0].intra * n as f64 } else { (0..n).map(|i| ix[self.species_of(i)].intra).sum::<f64>() };
        forces.energy += ew.energy(&sq) - ew.charge_self(&cq_original) - intra_total;
        let mut dipoles = Vec::new();
        let mut mu: Vec<Vec3> = Vec::new();
        if let Some(out) = out {
            forces.energy += out.energy;
            forces.virial += out.virial;
            mu = out.dipoles.iter().flatten().cloned().collect();
            for i in 0..n {
                let s = self.species_of(i);
                for (a, &k) in ix[s].pol_sites.iter().enumerate() {
                    let f = out.force_pol[i][a];
                    forces.force[i] += f;
                    forces.torque[i] += arms[i][k].cross(f);
                }
                for (c, &k) in ix[s].charged.iter().enumerate() {
                    let f = out.force_charge[i][c];
                    forces.force[i] += f;
                    forces.torque[i] += arms[i][k].cross(f);
                }
            }
            for (k, l) in links.iter().enumerate() {
                forces.force[l.i] += out.link_force[k];
                forces.force[l.j] -= out.link_force[k];
            }
            dipoles = out.dipoles;
        }
        // The reciprocal forces on every charge and dipole site, and the reciprocal virial.
        let stot: Vec<crate::ewald::C> = if fin.d.is_empty() { sq.clone() } else { sq.iter().zip(&fin.d).map(|(a, b)| crate::ewald::C { re: a.re + b.re, im: a.im + b.im }).collect() };
        let fch_rec: Vec<Vec3> = fin.grad_charge.iter().map(|g| -*g).collect();
        let virial_rec = if std::env::var_os("PHYS_EWALD_NUMERIC_VIRIAL").is_some() { ew.virial(&ccom, &carm, &cq, &dcom, &darm, &mu) } else { ew.virial_from_sums(&stot, &fin.tc, &fin.td, &fin.d) };
        let mut idx = 0;
        for i in 0..n {
            for &k in &ix[self.species_of(i)].pol_sites {
                let f = -fin.grad_dipole_pos[idx];
                forces.force[i] += f;
                forces.torque[i] += arms[i][k].cross(f);
                idx += 1;
            }
        }
        let mut idx = 0;
        for i in 0..n {
            for &(k, _) in &ix[self.species_of(i)].groups {
                let f = fch_rec[idx];
                forces.force[i] += f;
                forces.torque[i] += arms[i][k].cross(f);
                idx += 1;
            }
        }
        forces.virial += virial_rec;
        if std::env::var_os("PHYS_PROFILE").is_some() {
            eprintln!("  add_ewald: whole {:.1} ms", t_a.elapsed().as_secs_f64() * 1e3);
        }
        if !dipoles.is_empty() {
            forces.dipoles = dipoles;
        }
    }

    /// Kinetic energy: translational and rotational.
    pub fn kinetic(&self) -> (f64, f64) {
        let tr: f64 = self.momentum.iter().enumerate().map(|(k, p)| p.norm2() / (2.0 * self.kind_of(k).mass)).sum();
        let rot: f64 = self
            .spin
            .iter()
            .enumerate()
            .map(|(k, l)| {
                let i = self.kind_of(k).inertia;
                l.x * l.x / (2.0 * i[0]) + l.y * l.y / (2.0 * i[1]) + l.z * l.z / (2.0 * i[2])
            })
            .sum();
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
        for i in 0..self.com.len() {
            let (mass, inertia) = (self.kind_of(i).mass, self.kind_of(i).inertia);
            self.com[i] += self.momentum[i].scale(h / mass);
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
            for i in 0..self.com.len() {
                let (m, inertia) = (self.kind_of(i).mass, self.kind_of(i).inertia);
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

    fn charge_widths(&self) -> &[f64] {
        &self.law.sigma
    }

    fn site_site(&self) -> Option<&SiteSite> {
        Some(self.law)
    }
}

/// The relative tolerance of the induced-dipole solve in an Ewald run: 1e-10 unless
/// `PHYS_IND_TOL` says otherwise. Forces at fixed dipoles are linear in the dipoles, so
/// 1e-7 is ample for dynamics and saves iterations; the finite-difference checks want the
/// default.
fn ind_tolerance() -> f64 {
    static TOL: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *TOL.get_or_init(|| std::env::var("PHYS_IND_TOL").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-10))
}

/// A law with its Coulomb term cut down to the short-range remainder
/// `erf(r / s) - erf(r / s')` over `r` (`s'` the Ewald width), the rest being
/// the reciprocal sum `Liquid::add_ewald` adds (`ewald.rs`). Everything else the
/// inner law has - repulsion, dispersion, induction - is passed through.
pub struct EwaldLaw<'a> {
    pub inner: &'a dyn SiteLaw,
    pub ewald: crate::ewald::Ewald,
}

impl SiteLaw for EwaldLaw<'_> {
    fn site_pair(&self, ta: usize, tb: usize, r: f64) -> (f64, f64) {
        let (u, du) = self.inner.site_pair(ta, tb, r);
        if let Some(ss) = self.inner.site_site() {
            let qq = ss.charge[ta] * ss.charge[tb];
            if qq != 0.0 {
                let (g, dg) = crate::ewald::erf_over_r(r, self.ewald.s);
                return (u - qq * g, du - qq * dg);
            }
        }
        (u, du)
    }

    fn induction(&self) -> Option<(&[f64], &[f64])> {
        self.inner.induction()
    }

    fn charge_widths(&self) -> &[f64] {
        self.inner.charge_widths()
    }

    fn site_site(&self) -> Option<&SiteSite> {
        self.inner.site_site()
    }

    fn ewald(&self) -> Option<&crate::ewald::Ewald> {
        Some(&self.ewald)
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
    /// Per pair of types (`a * types + b`, symmetric), the exponent the
    /// dispersion's Tang-Toennies damping uses in place of the repulsion's
    /// `B`: where it is zero or the vector is empty the damping follows `B`.
    /// Hartree-Fock's repulsion and the correlation energy's dispersion are
    /// different physics and need not fall off together.
    pub damp: Vec<f64>,
    /// Per pair of types (`a * types + b`, symmetric), the strength `K` of a
    /// repulsion between two charge clouds, `K exp(-r^2 / s^2)` with `s^2 =
    /// sigma_a^2 + sigma_b^2`: the overlap of two Gaussian clouds, which is
    /// where a Pauli repulsion proportional to the overlap of valence densities
    /// (the AVDO model) comes to for the clouds `phys-esp --cloud` finds. Empty,
    /// or zero, is none.
    pub overlap: Vec<f64>,
}

impl SiteSite {
    /// As [`SiteLaw::site_pair`] with the two sites' charges `qa`, `qb` and
    /// Coulomb widths `sa`, `sb` given and not read from the table, which is
    /// what a molecule whose cloud follows its geometry needs. The overlap
    /// repulsion keeps the table's (base) widths. Returns the energy, its
    /// derivative in `r`, and the energy's derivatives with respect to
    /// `[qa, qb, sa, sb]`.
    pub fn site_pair_with(&self, ta: usize, tb: usize, r: f64, qa: f64, qb: f64, sa: f64, sb: f64) -> (f64, f64, [f64; 4]) {
        let types = self.charge.len();
        let [a, b, c6, c8] = self.pair[ta * types + tb];
        let qq = qa * qb;
        let s2 = sa * sa + sb * sb;
        let mut grad = [0.0; 4];
        let (coul, dcoul) = if s2 > 0.0 && qq != 0.0 {
            let sg = s2.sqrt();
            let x = r / sg;
            let g = crate::induction::erf(x);
            let e = (-x * x).exp();
            let c = 2.0 / (std::f64::consts::PI.sqrt() * sg) * e;
            // d/d sigma of erf(r / sg) / r is -(2 / sqrt(pi)) e^(-x^2) sigma / sg^3.
            let dsig = -qq * 2.0 / std::f64::consts::PI.sqrt() * e / (sg * sg * sg);
            grad = [qb * g / r, qa * g / r, dsig * sa, dsig * sb];
            (qq * g / r, qq * (c / r - g / (r * r)))
        } else {
            grad = [qb / r, qa / r, 0.0, 0.0];
            (qq / r, -qq / (r * r))
        };
        let (ov, dov) = {
            let k = self.overlap.get(ta * types + tb).copied().unwrap_or(0.0);
            match (k > 0.0, self.sigma.get(ta), self.sigma.get(tb)) {
                (true, Some(x), Some(y)) if x * x + y * y > 0.0 => {
                    let s2 = x * x + y * y;
                    let v = k * (-r * r / s2).exp();
                    (v, -2.0 * r / s2 * v)
                }
                _ => (0.0, 0.0),
            }
        };
        let rep = a * (-b * r).exp();
        let bd = self.damp.get(ta * types + tb).copied().filter(|v| *v > 0.0).unwrap_or(b);
        let (f6, df6) = tang_toennies(6, bd * r);
        let (f8, df8) = tang_toennies(8, bd * r);
        let (r6, r8) = (r.powi(6), r.powi(8));
        let u = coul + rep + ov - f6 * c6 / r6 - f8 * c8 / r8;
        let du = dcoul + dov - b * rep - (bd * df6 * c6 / r6 - 6.0 * f6 * c6 / (r6 * r)) - (bd * df8 * c8 / r8 - 8.0 * f8 * c8 / (r8 * r));
        (u, du, grad)
    }
}

impl SiteLaw for SiteSite {
    fn site_site(&self) -> Option<&SiteSite> {
        Some(self)
    }

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
        let (ov, dov) = {
            let k = self.overlap.get(ta * types + tb).copied().unwrap_or(0.0);
            match (k > 0.0, self.sigma.get(ta), self.sigma.get(tb)) {
                (true, Some(x), Some(y)) if x * x + y * y > 0.0 => {
                    let s2 = x * x + y * y;
                    let v = k * (-r * r / s2).exp();
                    (v, -2.0 * r / s2 * v)
                }
                _ => (0.0, 0.0),
            }
        };
        let rep = a * (-b * r).exp();
        let bd = self.damp.get(ta * types + tb).copied().filter(|v| *v > 0.0).unwrap_or(b);
        let (f6, df6) = tang_toennies(6, bd * r);
        let (f8, df8) = tang_toennies(8, bd * r);
        let (r6, r8) = (r.powi(6), r.powi(8));
        let u = coul + rep + ov - f6 * c6 / r6 - f8 * c8 / r8;
        let du = dcoul + dov - b * rep - (bd * df6 * c6 / r6 - 6.0 * f6 * c6 / (r6 * r)) - (bd * df8 * c8 / r8 - 8.0 * f8 * c8 / (r8 * r));
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
        e += induction_pair_energy(alpha, charge, law.charge_widths(), pair);
    }
    e
}

/// The induction energy of two molecules alone: their polarisable sites'
/// dipoles in each other's permanent charges, from `induction::solve` with one
/// link of weight one. Nothing for a molecule by itself (no field, no dipoles),
/// so this is the whole of induction's part of the pair's interaction energy.
pub fn induction_pair_energy(alpha: &[f64], charge: &[f64], widths: &[f64], pair: &PairEnergy) -> f64 {
    use crate::induction::{solve, Charge, Cluster, Link, PolSite};
    let mut cl = Cluster::default();
    for m in [&pair.a, &pair.b] {
        cl.pol.push(m.iter().filter(|(_, t)| alpha[*t] > 0.0).map(|(p, t)| PolSite { pos: *p, alpha: alpha[*t] }).collect());
        cl.charges.push(m.iter().filter(|(_, t)| charge[*t] != 0.0).map(|(p, t)| Charge { pos: *p, q: charge[*t], sigma: widths.get(*t).copied().unwrap_or(0.0) }).collect());
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

/// A site of a rigid molecule that is not an atom, placed in the molecule's own
/// frame so that it follows the atoms: from the first atom `o`, with `u` and
/// `v` the unit vectors to the second and third and `n` the unit normal of
/// their plane, at `o + a u + b v + c n` (bohr). Where Boys-localised orbitals
/// put a bond or a lone pair (`phys-esp --sites`) it is one of these; the
/// law's text carries one `site <type> <a> <b> <c>` line for each.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtraSite {
    pub ty: usize,
    pub frame: [f64; 3],
}

impl ExtraSite {
    pub fn position(&self, o: Vec3, p1: Vec3, p2: Vec3) -> Vec3 {
        let (u, v) = ((p1 - o).unit(), (p2 - o).unit());
        let n = u.cross(v).unit();
        o + u.scale(self.frame[0]) + v.scale(self.frame[1]) + n.scale(self.frame[2])
    }
}

/// The `site` lines of a law's text.
pub fn extra_sites_from_text(text: &str) -> Vec<ExtraSite> {
    text.lines()
        .filter_map(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            (w.first() == Some(&"site") && w.len() == 5).then(|| Some(ExtraSite { ty: w[1].parse().ok()?, frame: [w[2].parse().ok()?, w[3].parse().ok()?, w[4].parse().ok()?] })).flatten()
        })
        .collect()
}

/// A molecule's sites with `sites` added, in the order given.
pub fn with_frame_sites(molecule: &[(Vec3, usize)], sites: &[ExtraSite]) -> Vec<(Vec3, usize)> {
    let mut out = molecule.to_vec();
    for s in sites {
        out.push((s.position(molecule[0].0, molecule[1].0, molecule[2].0), s.ty));
    }
    out
}

/// A pair's energy at two geometries a little apart along the line between
/// the molecules' centres, and the reference's difference: the force a liquid's
/// pressure is made of, for a fit to match (`phys-es-gpu --deriv`).
#[derive(Debug, Clone)]
pub struct Slope {
    pub plus: PairEnergy,
    pub minus: PairEnergy,
    /// Reference energy at `plus` less at `minus`, hartree.
    pub target: f64,
    pub weight: f64,
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
    /// Fit a dispersion damping exponent for every pair of types.
    pub damp: bool,
    /// Keep every pair's `A` and `B`: the repulsion fitted elsewhere.
    pub repulsion: bool,
    /// Keep `C6` and `C8` at zero: a fit to Hartree-Fock's energies, which
    /// have no dispersion in them.
    pub no_dispersion: bool,
    /// Type pairs whose `C6` and `C8` are held at zero (an off-atom site that
    /// repels but does not disperse).
    pub no_dispersion_pairs: &'a [(usize, usize)],
    /// Keep every `C8` at its starting value (zero when the start has none):
    /// a free C8 with a weak damping becomes a long-range term the fit's
    /// weights (a floor of 1e-3 beyond the well) do not see.
    pub no_c8: bool,
    /// Fit the overlap repulsion's strength for every pair of types that both
    /// have a charge cloud.
    pub overlap: bool,
    /// Keep every `A` where it started (zero for a law with no Born-Mayer
    /// repulsion of its own), with `B` left as the dispersion's damping.
    pub no_born_mayer: bool,
    /// If positive, the Born-Mayer exponent of every pair is held at this value
    /// (bohr^-1) and only the amplitudes are fitted: the exponent of exchange
    /// repulsion is set by the density's tail, `2 sqrt(2 IP)` (`PLAY.md` E8d).
    pub bm_exponent: f64,
    /// A soft weight for the repulsive wall (hartree; 0 is off): a pair more than 25
    /// kcal/mol above the lowest weighs `0.05 / (1 + (E / soft_wall)^2)` and not the
    /// Boltzmann weight's floor, so that a wall pair counts for its *relative* error
    /// (for E well above `soft_wall` the weighted residual is about `0.22 soft_wall`
    /// times the fractional miss) and a pair of 2000 kcal/mol does not outweigh the
    /// whole well. The floor alone weighs every wall pair the same, so the few
    /// highest decide.
    pub soft_wall: f64,
    /// The soft weight's size at zero energy (0 means 0.05): the smaller it is, the less the
    /// wall pulls on the well.
    pub soft_weight: f64,
    /// Differences of energies along the centres' line to be matched as well.
    pub slopes: &'a [Slope],
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
    let has_cloud = |a: usize, b: usize| -> bool { start.sigma.get(a).copied().unwrap_or(0.0) > 0.0 && start.sigma.get(b).copied().unwrap_or(0.0) > 0.0 };
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
        if held.damp {
            for &(a, b) in &pairs {
                x.push(law.damp.get(a * types + b).copied().filter(|v| *v > 0.0).unwrap_or(law.pair[a * types + b][1]).max(1e-300).ln());
            }
        }
        if held.overlap {
            for &(a, b) in &pairs {
                if has_cloud(a, b) {
                    x.push(law.overlap.get(a * types + b).copied().filter(|v| *v > 0.0).unwrap_or(0.05).ln());
                }
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
        let d0 = types - 1 + 4 * pairs.len() + if held.sigma { types } else { 0 };
        let damp = if held.damp {
            let mut d = vec![0.0; types * types];
            for (k, &(a, b)) in pairs.iter().enumerate() {
                d[a * types + b] = x[d0 + k].exp();
                d[b * types + a] = d[a * types + b];
            }
            d
        } else {
            start.damp.clone()
        };
        let o0 = d0 + if held.damp { pairs.len() } else { 0 };
        let overlap = if held.overlap {
            let mut o = vec![0.0; types * types];
            let mut k = 0;
            for &(a, b) in &pairs {
                if has_cloud(a, b) {
                    o[a * types + b] = x[o0 + k].exp();
                    o[b * types + a] = o[a * types + b];
                    k += 1;
                }
            }
            o
        } else {
            start.overlap.clone()
        };
        SiteSite { charge, pair, sigma, damp, overlap }
    };
    // What is fitted: every parameter but the held ones, which keep what
    // `start` gave them.
    let mut base = pack(start);
    if held.bm_exponent > 0.0 {
        for k in 0..pairs.len() {
            base[types - 1 + 4 * k + 1] = held.bm_exponent.ln();
        }
    }
    let is_held = |list: &[(usize, usize)], a: usize, b: usize| list.iter().any(|&(p, q)| (p, q) == (a, b) || (q, p) == (a, b));
    let mut free: Vec<usize> = if held.charges { Vec::new() } else { (0..types - 1).collect() };
    for (k, &(a, b)) in pairs.iter().enumerate() {
        let o = types - 1 + 4 * k;
        if is_held(held.pairs, a, b) {
            continue;
        }
        for j in 0..4 {
            if (j == 1 && held.bm_exponent > 0.0) || (held.repulsion && j < 2) || (held.no_born_mayer && j < 2) || (held.no_dispersion && j >= 2) || (is_held(held.no_dispersion_pairs, a, b) && j >= 2) || (held.no_c8 && j == 3) {
                continue;
            }
            if !(j == 2 && is_held(held.dispersion, a, b)) {
                free.push(o + j);
            }
        }
    }
    if held.sigma {
        free.extend((0..types).map(|t| types - 1 + 4 * pairs.len() + t));
    }
    if held.overlap {
        let o0 = types - 1 + 4 * pairs.len() + if held.sigma { types } else { 0 } + if held.damp { pairs.len() } else { 0 };
        let count = pairs.iter().filter(|&&(a, b)| has_cloud(a, b)).count();
        free.extend(o0..o0 + count);
    }
    if held.damp {
        let d0 = types - 1 + 4 * pairs.len() + if held.sigma { types } else { 0 };
        free.extend(pairs.iter().enumerate().filter(|(_, &(a, b))| !is_held(held.pairs, a, b)).map(|(k, _)| d0 + k));
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
    let weights: Vec<f64> = data
        .iter()
        .map(|d| {
            let above = d.energy - e_min;
            if held.soft_wall > 0.0 && above > 25.0 / 627.509474 {
                (if held.soft_weight > 0.0 { held.soft_weight } else { 0.05 }) / (1.0 + (d.energy / held.soft_wall).powi(2))
            } else {
                (-above / kt).exp().max(floor)
            }
        })
        .collect();
    let energy_of = |law: &SiteSite, d: &PairEnergy| -> f64 {
        match alpha {
            Some(a) => pair_energy(&Polarisable { law, alpha: a.to_vec() }, d),
            None => pair_energy(law, d),
        }
    };
    let residuals = |x: &[f64]| -> Vec<f64> {
        let law = unpack(&expand(x));
        let mut r: Vec<f64> = data.iter().zip(&weights).map(|(d, w)| w.sqrt() * (energy_of(&law, d) - d.energy)).collect();
        for sl in held.slopes {
            r.push(sl.weight.sqrt() * ((energy_of(&law, &sl.plus) - energy_of(&law, &sl.minus)) - sl.target));
        }
        r
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
        Liquid::mixture(vec![(kind, n)], density, side, length, temperature, seed, r_on, r_cut)
    }

    /// A slab of several kinds of molecule, `(kind, count)` each, mixed at random over the lattice
    /// sites (a seeded shuffle, its own stream, so a single kind draws exactly what [`Liquid::slab`]
    /// always did). The first kind is `kind`; the law's site types must be disjoint between kinds.
    pub fn mixture(kinds: Vec<(Kind, usize)>, density: f64, side: f64, length: f64, temperature: f64, seed: u64, r_on: f64, r_cut: f64) -> Liquid {
        let n: usize = kinds.iter().map(|k| k.1).sum();
        let mut species: Vec<usize> = kinds.iter().enumerate().flat_map(|(s, k)| std::iter::repeat_n(s, k.1)).collect();
        if kinds.len() > 1 {
            let mut shuffle = Stream::at(seed, 1, 0, Purpose::Positions);
            for i in (1..species.len()).rev() {
                let j = ((shuffle.uniform() * (i + 1) as f64) as usize).min(i);
                species.swap(i, j);
            }
        }
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
                    let kind = &kinds[species[com.len()]].0;
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
        let mut kinds = kinds.into_iter().map(|k| k.0);
        let kind = kinds.next().expect("a kind");
        let other_kinds: Vec<Kind> = kinds.collect();
        let species = if other_kinds.is_empty() { Vec::new() } else { species };
        Liquid { kind, other_kinds, species, cell: [side, side, length], com, orientation, momentum, spin, r_on, r_cut }
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
        for a in 0..types {
            for b in a..types {
                let d = self.damp.get(a * types + b).copied().unwrap_or(0.0);
                if d > 0.0 {
                    s += &format!("damp {a} {b} {d:e}
");
                }
            }
        }
        for a in 0..types {
            for b in a..types {
                let k = self.overlap.get(a * types + b).copied().unwrap_or(0.0);
                if k > 0.0 {
                    s += &format!("overlap {a} {b} {k:e}\n");
                }
            }
        }
        if self.sigma.iter().any(|v| *v > 0.0) {
            s += &format!("sigma {}
", self.sigma.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "));
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
        let mut damps: Vec<(usize, usize, f64)> = Vec::new();
        let mut overlaps: Vec<(usize, usize, f64)> = Vec::new();
        let mut pairs = Vec::new();
        for line in text.lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            match w.first().copied() {
                Some("charges") => charge = w[1..].iter().map(|x| x.parse().ok()).collect::<Option<Vec<f64>>>()?,
                // Read by `bisector_from_text`: where an off-atom site goes, not
                // part of the pair law.
                Some("bisector") => {}
                // Read by `extra_sites_from_text`.
                Some("site") => {}
                // Read by `alpha_from_text`: the polarisability of a site type.
                Some("alpha") => {}
                // Read by `geolaw::GeoCloud::from_text`: a cloud that follows the angle.
                Some("geo") | Some("geopoly") => {}
                Some("overlap") if w.len() == 4 => overlaps.push((w[1].parse().ok()?, w[2].parse().ok()?, w[3].parse().ok()?)),
                Some("damp") if w.len() == 4 => damps.push((w[1].parse().ok()?, w[2].parse().ok()?, w[3].parse().ok()?)),
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
        let mut damp = Vec::new();
        if !damps.is_empty() {
            damp = vec![0.0; types * types];
            for (a, b, v) in damps {
                damp[a * types + b] = v;
                damp[b * types + a] = v;
            }
        }
        let mut overlap = Vec::new();
        if !overlaps.is_empty() {
            overlap = vec![0.0; types * types];
            for (a, b, v) in overlaps {
                overlap[a * types + b] = v;
                overlap[b * types + a] = v;
            }
        }
        Some(SiteSite { charge, pair, sigma, damp, overlap })
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
        Kind::of_molecule_sites(z, positions, types, bisector, &[])
    }

    /// [`Kind::of_molecule`] with frame sites ([`ExtraSite`]) as well, each
    /// massless like the bisector site.
    pub fn of_molecule_sites(z: &[u32], positions: &[[f64; 3]], types: &[usize], bisector: Option<(usize, f64)>, frame: &[ExtraSite]) -> Kind {
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
        for s in frame {
            let at = |k: usize| Vec3 { x: positions[k][0], y: positions[k][1], z: positions[k][2] };
            let m = s.position(at(0), at(1), at(2));
            pos.push([m.x, m.y, m.z]);
            ty.push(s.ty);
            masses.push(0.0);
        }
        Kind::from_atoms(&pos, &masses, &ty)
    }
}
