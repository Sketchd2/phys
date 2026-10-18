//! Water that moves atom by atom — `docs/PLAY.md` E8c.
//!
//! [`liquid`](crate::liquid) moves rigid bodies: a centre, an orientation, and
//! the forces and torques on them. A molecule that bends and stretches needs the
//! force on each atom, and this module is that: the same intermolecular law
//! (`SiteLaw`, the induced dipoles of `induction.rs`, the same switch on the
//! molecules' centres) with its sites attached to atoms, or placed from them,
//! and the molecule's own energy from the fitted surface `phys-pes` and
//! `phys-vib` made, on which the atoms move.
//!
//! **Sites.** A site of the law is either on an atom (the cores, the clouds on
//! the oxygen and the hydrogens) or the one virtual site, on the bisector of the
//! oxygen's two bonds. The virtual site has no mass and no force of its own: the
//! force on it is handed back to the three atoms through the derivative of where
//! it is, [`bisector_transfer`].
//!
//! **Two ways to move.** Fully flexible, every atom free; or with the two O-H
//! lengths held (SHAKE and RATTLE), the bend free, which is the treatment E8c
//! calls the second. The integrator is BAOAB, a Langevin thermostat on the atoms:
//! it samples the same ensemble as any other and the friction only sets how fast.
//!
//! **What is not here.** The law's own numbers (charges, widths, sizes) do not
//! change with the molecule's shape: they are those of the shape the pair
//! energies were computed at, and only the places of the sites move with the
//! atoms. A molecule's dipole therefore changes as its bonds lengthen through the
//! cloud sites moving and the induced dipoles, not through the clouds' charges
//! changing; making them follow the geometry is a stage of its own.

use crate::induction::{solve, Charge, Cluster, Link, PolSite};
use crate::liquid::{switch, SiteLaw};
use crate::math::Vec3;
use crate::rng::{Purpose, Stream};

/// Where a site is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Anchor {
    Atom(usize),
    /// On the bisector of the first atom's bonds to the other two, this far
    /// (bohr) from it.
    Bisector(f64),
}

/// The law's sites of one molecule, in order: each a type and an anchor.
#[derive(Debug, Clone, PartialEq)]
pub struct FlexSites {
    pub types: Vec<usize>,
    pub anchors: Vec<Anchor>,
}

impl FlexSites {
    /// From a law's text: the atoms with their types, then the `bisector` line,
    /// then each `site` line, which must lie on an atom (frame `(0,0,0)` the
    /// first, `(a,0,0)` the second, `(0,b,0)` the third; the cloud model's are).
    pub fn from_text(text: &str, atom_types: [usize; 3]) -> FlexSites {
        let mut types: Vec<usize> = atom_types.to_vec();
        let mut anchors = vec![Anchor::Atom(0), Anchor::Atom(1), Anchor::Atom(2)];
        if let Some((t, d)) = crate::liquid::SiteSite::bisector_from_text(text) {
            types.push(t);
            anchors.push(Anchor::Bisector(d));
        }
        for s in crate::liquid::extra_sites_from_text(text) {
            let [a, b, c] = s.frame;
            let near = |x: f64| x.abs() < 1e-9;
            let atom = if near(a) && near(b) && near(c) {
                0
            } else if !near(a) && near(b) && near(c) {
                1
            } else if near(a) && !near(b) && near(c) {
                2
            } else {
                panic!("the site of type {} at frame ({a}, {b}, {c}) is not on an atom, and a flexible molecule has no frame to hold it in", s.ty)
            };
            types.push(s.ty);
            anchors.push(Anchor::Atom(atom));
        }
        FlexSites { types, anchors }
    }
}

/// The molecule's own energy as a polynomial in its two bond lengths and its
/// angle, the fit `phys-vib` writes to `surface-water.txt`.
#[derive(Debug, Clone)]
pub struct Surface {
    r0: f64,
    th0: f64,
    e_ref: f64,
    terms: Vec<(i32, i32, i32, f64)>,
    /// The energy at the surface's minimum, the zero of [`Surface::eval`].
    pub minimum: f64,
}

impl Surface {
    pub fn from_text(text: &str) -> Surface {
        let (mut r0, mut th0, mut e_ref) = (0.0, 0.0, 0.0);
        let mut terms = Vec::new();
        for l in text.lines() {
            let w: Vec<&str> = l.split_whitespace().collect();
            match w.first().copied() {
                Some("ref") => {
                    r0 = w[1].parse().expect("r0");
                    th0 = w[2].parse().expect("theta0");
                    e_ref = w[3].parse().expect("e_ref");
                }
                Some("term") => terms.push((w[1].parse().expect("a"), w[2].parse().expect("b"), w[3].parse().expect("c"), w[4].parse().expect("coefficient"))),
                _ => {}
            }
        }
        let mut s = Surface { r0, th0, e_ref, terms, minimum: 0.0 };
        // Newton on the surface for its minimum.
        let mut q = [r0, r0, th0];
        for _ in 0..40 {
            let h = 1e-4;
            let g = s.raw(q[0], q[1], q[2]).1;
            let mut hess = [[0.0; 3]; 3];
            for j in 0..3 {
                let (mut p, mut m) = (q, q);
                p[j] += h;
                m[j] -= h;
                let (gp, gm) = (s.raw(p[0], p[1], p[2]).1, s.raw(m[0], m[1], m[2]).1);
                for i in 0..3 {
                    hess[i][j] = (gp[i] - gm[i]) / (2.0 * h);
                }
            }
            let step = solve3(hess, g);
            for i in 0..3 {
                q[i] -= step[i];
            }
            if step.iter().map(|x| x.abs()).fold(0.0, f64::max) < 1e-13 {
                break;
            }
        }
        s.minimum = s.raw(q[0], q[1], q[2]).0;
        s
    }

    fn raw(&self, r1: f64, r2: f64, th: f64) -> (f64, [f64; 3]) {
        let (x, y, z) = ((r1 - self.r0) / 0.1, (r2 - self.r0) / 0.1, (th - self.th0) / 0.2);
        let mut e = self.e_ref;
        let mut g = [0.0; 3];
        let pw = |v: f64, n: i32| if n < 0 { 0.0 } else { v.powi(n) };
        for &(a, b, c, k) in &self.terms {
            e += k * pw(x, a) * pw(y, b) * pw(z, c);
            g[0] += k * (a as f64) * pw(x, a - 1) * pw(y, b) * pw(z, c) / 0.1;
            g[1] += k * (b as f64) * pw(x, a) * pw(y, b - 1) * pw(z, c) / 0.1;
            g[2] += k * (c as f64) * pw(x, a) * pw(y, b) * pw(z, c - 1) / 0.2;
        }
        (e, g)
    }

    /// The energy above the minimum and its gradient with respect to
    /// `(r1, r2, theta)`.
    pub fn eval(&self, r1: f64, r2: f64, th: f64) -> (f64, [f64; 3]) {
        let (e, g) = self.raw(r1, r2, th);
        (e - self.minimum, g)
    }
}

fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    let cof = |r: usize, c: usize| {
        let (r1, r2) = ((r + 1) % 3, (r + 2) % 3);
        let (c1, c2) = ((c + 1) % 3, (c + 2) % 3);
        a[r1][c1] * a[r2][c2] - a[r1][c2] * a[r2][c1]
    };
    let mut x = [0.0; 3];
    for i in 0..3 {
        x[i] = (0..3).map(|j| cof(j, i) * b[j]).sum::<f64>() / det;
    }
    x
}

/// The force on each atom of a molecule from the surface, and its energy.
pub fn intramolecular(surface: &Surface, pos: &[Vec3; 3]) -> (f64, [Vec3; 3]) {
    let a = pos[1] - pos[0];
    let b = pos[2] - pos[0];
    let (r1, r2) = (a.norm(), b.norm());
    let (u, v) = (a.scale(1.0 / r1), b.scale(1.0 / r2));
    let c = u.dot(v).clamp(-1.0, 1.0);
    let th = c.acos();
    let s = (1.0 - c * c).sqrt();
    let (e, [d1, d2, dt]) = surface.eval(r1, r2, th);
    // d theta / d H1 = -(v - c u) / (r1 sin theta), and the mirror for H2.
    let dth1 = (v - u.scale(c)).scale(-1.0 / (r1 * s));
    let dth2 = (u - v.scale(c)).scale(-1.0 / (r2 * s));
    let f1 = -(u.scale(d1) + dth1.scale(dt));
    let f2 = -(v.scale(d2) + dth2.scale(dt));
    (e, [-(f1 + f2), f1, f2])
}

/// The gradient of the H-O-H angle with respect to the three atoms.
pub fn angle_gradient(pos: &[Vec3; 3]) -> [Vec3; 3] {
    let a = pos[1] - pos[0];
    let b = pos[2] - pos[0];
    let (r1, r2) = (a.norm(), b.norm());
    let (u, v) = (a.scale(1.0 / r1), b.scale(1.0 / r2));
    let c = u.dot(v).clamp(-1.0, 1.0);
    let sn = (1.0 - c * c).sqrt();
    let dth1 = (v - u.scale(c)).scale(-1.0 / (r1 * sn));
    let dth2 = (u - v.scale(c)).scale(-1.0 / (r2 * sn));
    [-(dth1 + dth2), dth1, dth2]
}

/// The force `f` on the virtual bisector site `d` from the first atom handed to
/// the three atoms: `M = O + d w / |w|` with `w = u + v`, the unit vectors of
/// the two bonds.
pub fn bisector_transfer(pos: &[Vec3; 3], d: f64, f: Vec3) -> [Vec3; 3] {
    let a = pos[1] - pos[0];
    let b = pos[2] - pos[0];
    let (r1, r2) = (a.norm(), b.norm());
    let (u, v) = (a.scale(1.0 / r1), b.scale(1.0 / r2));
    let w = u + v;
    let wn = w.norm();
    let wh = w.scale(1.0 / wn);
    let g = (f - wh.scale(wh.dot(f))).scale(d / wn);
    let f1 = (g - u.scale(u.dot(g))).scale(1.0 / r1);
    let f2 = (g - v.scale(v.dot(g))).scale(1.0 / r2);
    [f - f1 - f2, f1, f2]
}

pub struct Flex {
    pub sites: FlexSites,
    /// Atomic masses (atomic units) of the molecule's three atoms.
    pub masses: [f64; 3],
    pub cell: [f64; 3],
    pub pos: Vec<[Vec3; 3]>,
    pub vel: Vec<[Vec3; 3]>,
    pub r_on: f64,
    pub r_cut: f64,
    pub surface: Surface,
    /// Hold both O-H lengths at this (bohr), or leave them free.
    pub bond: Option<f64>,
    /// A cloud that follows each molecule's angle (`geolaw.rs`): the law must be
    /// a `SiteSite` (or wrap one), whose table gives the pair terms and the base
    /// widths, while each molecule's charges, Coulomb widths and bisector
    /// distance come from here at its own angle.
    pub geo: Option<crate::geolaw::GeoCloud>,
}

pub struct FlexForces {
    /// The law's energy (pairs and induction), and the molecules' own.
    pub inter: f64,
    pub intra: f64,
    pub force: Vec<[Vec3; 3]>,
    /// The molecular virial over the molecules' centres of mass, as `liquid`'s.
    pub virial: f64,
    pub dipoles: Vec<Vec<Vec3>>,
}

impl Flex {
    fn total_mass(&self) -> f64 {
        self.masses.iter().sum()
    }

    pub fn centre(&self, i: usize) -> Vec3 {
        let m = self.total_mass();
        (0..3).fold(Vec3::ZERO, |acc, a| acc + self.pos[i][a].scale(self.masses[a] / m))
    }

    fn image(&self, mut d: Vec3) -> Vec3 {
        let c = self.cell;
        d.x -= c[0] * (d.x / c[0]).round();
        d.y -= c[1] * (d.y / c[1]).round();
        d.z -= c[2] * (d.z / c[2]).round();
        d
    }

    pub fn site_positions(&self, i: usize) -> Vec<Vec3> {
        let p = &self.pos[i];
        let geo_d = self.geo.as_ref().map(|g| g.params(crate::geolaw::GeoCloud::angle(p)).d);
        self.sites
            .anchors
            .iter()
            .map(|a| match *a {
                Anchor::Atom(k) => p[k],
                Anchor::Bisector(d) => {
                    let d = geo_d.unwrap_or(d);
                    let u = (p[1] - p[0]).unit();
                    let v = (p[2] - p[0]).unit();
                    p[0] + (u + v).unit().scale(d)
                }
            })
            .collect()
    }

    /// Forces on every site folded down to the atoms.
    fn to_atoms(&self, i: usize, site_forces: &[Vec3]) -> [Vec3; 3] {
        let mut out = [Vec3::ZERO; 3];
        let geo_d = self.geo.as_ref().map(|g| g.params(crate::geolaw::GeoCloud::angle(&self.pos[i])).d);
        for (k, a) in self.sites.anchors.iter().enumerate() {
            match *a {
                Anchor::Atom(j) => out[j] += site_forces[k],
                Anchor::Bisector(d) => {
                    let d = geo_d.unwrap_or(d);
                    let t = bisector_transfer(&self.pos[i], d, site_forces[k]);
                    for j in 0..3 {
                        out[j] += t[j];
                    }
                }
            }
        }
        out
    }

    pub fn forces(&self, law: &dyn SiteLaw, warm: Option<&[Vec<Vec3>]>) -> FlexForces {
        let n = self.pos.len();
        let ns = self.sites.types.len();
        let centres: Vec<Vec3> = (0..n).map(|i| self.centre(i)).collect();
        let sites: Vec<Vec<Vec3>> = (0..n).map(|i| self.site_positions(i)).collect();
        let types = &self.sites.types;
        let mtot = self.total_mass();
        // A cloud that follows each molecule's angle: its parameters, and the law's table.
        let geo_params: Option<Vec<crate::geolaw::CloudParams>> = self.geo.as_ref().map(|g| self.pos.iter().map(|p| g.params(crate::geolaw::GeoCloud::angle(p))).collect());
        let table: Option<&crate::liquid::SiteSite> = if self.geo.is_some() { Some(law.site_site().expect("a geometry-following cloud needs a site-site law")) } else { None };
        // The pairs of molecules inside the cut-off, for the induction and for a pair-loop engine.
        let pair_engine = crate::liquid::pair_engine().filter(|_| law.site_site().is_some());
        let mut links: Vec<Link> = Vec::new();
        if law.induction().is_some() || pair_engine.is_some() {
            for i in 0..n {
                for j in i + 1..n {
                    let raw = centres[i] - centres[j];
                    let d = self.image(raw);
                    let dist = d.norm();
                    if dist >= self.r_cut {
                        continue;
                    }
                    let (wgt, dw) = switch(dist, self.r_on, self.r_cut);
                    links.push(Link { i, j, shift: raw - d, weight: wgt, dweight: dw, d });
                }
            }
        }
        let engine_result = pair_engine.as_ref().and_then(|e| {
            let tab = law.site_site()?;
            let rel: Vec<Vec<Vec3>> = (0..n).map(|i| sites[i].iter().map(|x| *x - centres[i]).collect()).collect();
            let (qs, ss): (Vec<Vec<f64>>, Vec<Vec<f64>>) = match &geo_params {
                Some(gp) => ((0..n).map(|i| types.iter().map(|&t| gp[i].charge[t]).collect()).collect(), (0..n).map(|i| types.iter().map(|&t| gp[i].sigma[t]).collect()).collect()),
                None => ((0..n).map(|_| types.iter().map(|&t| tab.charge[t]).collect()).collect(), (0..n).map(|_| types.iter().map(|&t| tab.sigma.get(t).copied().unwrap_or(0.0)).collect()).collect()),
            };
            e.pairs(&crate::liquid::PairProblem { table: tab, types, rel: &rel, charge: &qs, sigma: &ss, links: &links, ewald_width: None })
        });
        let job = |idx: &mut dyn Iterator<Item = usize>| {
            idx.map(|i| {
                let mut fsite = vec![Vec3::ZERO; ns];
                let (mut gq, mut gs) = (vec![0.0; ns], vec![0.0; ns]);
                let (mut e, mut fc, mut w) = (0.0, Vec3::ZERO, 0.0);
                for j in 0..n {
                    if j == i {
                        continue;
                    }
                    let d = self.image(centres[i] - centres[j]);
                    let dist = d.norm();
                    if dist >= self.r_cut {
                        continue;
                    }
                    let (s, ds) = switch(dist, self.r_on, self.r_cut);
                    let mut u = 0.0;
                    let mut net = Vec3::ZERO;
                    for a in 0..ns {
                        let pa = sites[i][a] - centres[i];
                        for b in 0..ns {
                            let pb = sites[j][b] - centres[j];
                            let r = d + pa - pb;
                            let rr = r.norm();
                            let (v, dv) = match (&geo_params, table) {
                                (Some(gp), Some(t)) => {
                                    let (v, dv, g) = t.site_pair_with(types[a], types[b], rr, gp[i].charge[types[a]], gp[j].charge[types[b]], gp[i].sigma[types[a]], gp[j].sigma[types[b]]);
                                    gq[a] += s * g[0];
                                    gs[a] += s * g[2];
                                    (v, dv)
                                }
                                _ => law.site_pair(types[a], types[b], rr),
                            };
                            u += v;
                            let f = r.scale(-dv / rr);
                            fsite[a] += f.scale(s);
                            net += f.scale(s);
                        }
                    }
                    e += 0.5 * s * u;
                    let switch_force = d.scale(-u * ds / dist);
                    fc += switch_force;
                    w += d.dot(net + switch_force);
                }
                (i, e, fsite, fc, w, gq, gs)
            })
            .collect::<Vec<_>>()
        };
        let mut energy = 0.0;
        let mut virial = 0.0;
        let mut force = vec![[Vec3::ZERO; 3]; n];
        // Per molecule, per site: the energy's derivative with respect to the site's charge and
        // Coulomb width, and the force on the bisector site, for the angle's part of the force.
        let mut gq = vec![vec![0.0; ns]; n];
        let mut gs = vec![vec![0.0; ns]; n];
        let mut fbis = vec![Vec3::ZERO; n];
        let bis_site: Option<usize> = self.sites.anchors.iter().position(|a| matches!(a, Anchor::Bisector(_)));
        let parts = match engine_result {
            Some(r) => {
                energy = r.energy;
                virial = r.virial;
                vec![(0..n).map(|i| (i, 0.0, r.fsite[i].clone(), r.fcentre[i], 0.0, r.gq[i].clone(), r.gs[i].clone())).collect::<Vec<_>>()]
            }
            None => crate::electrons::scf::parallel_interleaved(n, &job),
        };
        for part in parts {
            for (i, e, fsite, fc, w, q, sg) in part {
                energy += e;
                virial += 0.5 * w;
                if let Some(kb) = bis_site {
                    fbis[i] += fsite[kb];
                }
                gq[i] = q;
                gs[i] = sg;
                let mut atoms = self.to_atoms(i, &fsite);
                for a in 0..3 {
                    atoms[a] += fc.scale(self.masses[a] / mtot);
                }
                for a in 0..3 {
                    force[i][a] += atoms[a];
                }
            }
        }
        let mut dipoles = Vec::new();
        if let Some((alpha, charge)) = law.induction() {
            let widths = law.charge_widths();
            let pol: Vec<usize> = (0..ns).filter(|&k| alpha[types[k]] > 0.0).collect();
            let charged: Vec<usize> = (0..ns).filter(|&k| charge[types[k]] != 0.0).collect();
            let mut cl = Cluster::default();
            for i in 0..n {
                cl.pol.push(pol.iter().map(|&k| PolSite { pos: sites[i][k], alpha: alpha[types[k]] }).collect());
                cl.charges.push(
                    charged
                        .iter()
                        .map(|&k| match &geo_params {
                            Some(gp) => Charge { pos: sites[i][k], q: gp[i].charge[types[k]], sigma: gp[i].sigma[types[k]] },
                            None => Charge { pos: sites[i][k], q: charge[types[k]], sigma: widths.get(types[k]).copied().unwrap_or(0.0) },
                        })
                        .collect(),
                );
            }
            let out = solve(&cl, &links, warm, 1e-9);
            energy += out.energy;
            virial += out.virial;
            if geo_params.is_some() {
                let (dq, ds) = crate::induction::charge_gradients(&cl, &links, &out.dipoles);
                for i in 0..n {
                    for (c, &k) in charged.iter().enumerate() {
                        gq[i][k] += dq[i][c];
                        gs[i][k] += ds[i][c];
                    }
                }
            }
            for i in 0..n {
                let mut fsite = vec![Vec3::ZERO; ns];
                for (a, &k) in pol.iter().enumerate() {
                    fsite[k] += out.force_pol[i][a];
                }
                for (c, &k) in charged.iter().enumerate() {
                    fsite[k] += out.force_charge[i][c];
                }
                if let Some(kb) = bis_site {
                    fbis[i] += fsite[kb];
                }
                let atoms = self.to_atoms(i, &fsite);
                for a in 0..3 {
                    force[i][a] += atoms[a];
                }
            }
            for (k, l) in links.iter().enumerate() {
                for a in 0..3 {
                    force[l.i][a] += out.link_force[k].scale(self.masses[a] / mtot);
                    force[l.j][a] -= out.link_force[k].scale(self.masses[a] / mtot);
                }
            }
            dipoles = out.dipoles;
        }
        // The angle enters the energy twice more: through each charge and width, and through the
        // bisector cloud's distance. With theta the molecule's angle, the generalised force
        // `F_bis . w_hat * d'(theta) - dE/dtheta` acts along `grad theta`.
        if let Some(g) = &self.geo {
            for i in 0..n {
                let th = crate::geolaw::GeoCloud::angle(&self.pos[i]);
                let dp = g.dparams(th);
                let mut de_dth = 0.0;
                for k in 0..ns {
                    de_dth += gq[i][k] * dp.charge[types[k]] + gs[i][k] * dp.sigma[types[k]];
                }
                let (u, v) = ((self.pos[i][1] - self.pos[i][0]).unit(), (self.pos[i][2] - self.pos[i][0]).unit());
                let kappa = fbis[i].dot((u + v).unit()) * dp.d - de_dth;
                let grad = angle_gradient(&self.pos[i]);
                for a in 0..3 {
                    force[i][a] += grad[a].scale(kappa);
                }
            }
        }
        let mut intra = 0.0;
        for i in 0..n {
            let (e, f) = intramolecular(&self.surface, &self.pos[i]);
            intra += e;
            for a in 0..3 {
                force[i][a] += f[a];
            }
        }
        FlexForces { inter: energy, intra, force, virial, dipoles }
    }

    /// Kinetic energy of the molecules' centres of mass, and in all.
    pub fn kinetic(&self) -> (f64, f64) {
        let m = self.total_mass();
        let mut com = 0.0;
        let mut all = 0.0;
        for i in 0..self.pos.len() {
            let p = (0..3).fold(Vec3::ZERO, |acc, a| acc + self.vel[i][a].scale(self.masses[a]));
            com += p.norm2() / (2.0 * m);
            all += (0..3).map(|a| 0.5 * self.masses[a] * self.vel[i][a].norm2()).sum::<f64>();
        }
        (com, all)
    }

    /// Pressure, hartree per cubic bohr, from the centres' kinetic energy and
    /// the molecular virial.
    pub fn pressure(&self, f: &FlexForces) -> f64 {
        let volume = self.cell[0] * self.cell[1] * self.cell[2];
        (2.0 * self.kinetic().0 + f.virial) / (3.0 * volume)
    }

    /// Degrees of freedom: nine an atom triple, less the constraints, less the
    /// centre of mass of the whole.
    pub fn degrees_of_freedom(&self) -> f64 {
        let per = if self.bond.is_some() { 7.0 } else { 9.0 };
        per * self.pos.len() as f64 - 3.0
    }

    /// Put each molecule's bonds back at their length after a move, SHAKE:
    /// along the bonds as they were before it.
    fn shake(&mut self, before: &[[Vec3; 3]], h: f64) {
        let Some(len) = self.bond else { return };
        let (m0, m1, m2) = (self.masses[0], self.masses[1], self.masses[2]);
        for i in 0..self.pos.len() {
            let old = before[i];
            let pairs = [(0usize, 1usize), (0, 2)];
            for _ in 0..200 {
                let mut worst: f64 = 0.0;
                for &(a, b) in &pairs {
                    let r = self.pos[i][b] - self.pos[i][a];
                    let diff = r.norm2() - len * len;
                    worst = worst.max(diff.abs() / (len * len));
                    let r_old = old[b] - old[a];
                    let (ma, mb) = (self.masses[a], self.masses[b]);
                    let g = diff / (2.0 * (1.0 / ma + 1.0 / mb) * r.dot(r_old));
                    self.pos[i][a] += r_old.scale(g / ma);
                    self.pos[i][b] -= r_old.scale(g / mb);
                    // The velocity follows the correction.
                    self.vel[i][a] += r_old.scale(g / (ma * h));
                    self.vel[i][b] -= r_old.scale(g / (mb * h));
                }
                if worst < 1e-12 {
                    break;
                }
            }
            let _ = (m0, m1, m2);
        }
    }

    /// Take the part of the velocities along the bonds out, RATTLE.
    fn rattle(&mut self) {
        if self.bond.is_none() {
            return;
        }
        for i in 0..self.pos.len() {
            for _ in 0..200 {
                let mut worst: f64 = 0.0;
                for &(a, b) in &[(0usize, 1usize), (0, 2)] {
                    let r = self.pos[i][b] - self.pos[i][a];
                    let vrel = self.vel[i][b] - self.vel[i][a];
                    let (ma, mb) = (self.masses[a], self.masses[b]);
                    let k = vrel.dot(r) / ((1.0 / ma + 1.0 / mb) * r.norm2());
                    self.vel[i][a] += r.scale(k / ma);
                    self.vel[i][b] -= r.scale(k / mb);
                    worst = worst.max(k.abs());
                }
                if worst < 1e-14 {
                    break;
                }
            }
        }
    }

    /// One BAOAB step of length `h`: kick, drift, thermostat, drift, kick. With
    /// `bath` (temperature in hartree, friction per atomic time unit, stream) the
    /// thermostat acts; without it the step is plain velocity Verlet.
    pub fn step(&mut self, law: &dyn SiteLaw, mut f: FlexForces, h: f64, mut bath: Option<(f64, f64, &mut Stream)>) -> FlexForces {
        let n = self.pos.len();
        let kick = |s: &mut Flex, f: &FlexForces, dt: f64| {
            for i in 0..n {
                for a in 0..3 {
                    s.vel[i][a] += f.force[i][a].scale(dt / s.masses[a]);
                }
            }
        };
        kick(self, &f, 0.5 * h);
        self.rattle();
        for half in 0..2 {
            let before = self.pos.clone();
            for i in 0..n {
                for a in 0..3 {
                    self.pos[i][a] += self.vel[i][a].scale(0.5 * h);
                }
            }
            self.shake(&before, 0.5 * h);
            if half == 0 {
                if let Some((kt, friction, stream)) = bath.as_mut() {
                    let c = (-*friction * h).exp();
                    let s2 = 1.0 - c * c;
                    for i in 0..n {
                        for a in 0..3 {
                            let sd = (s2 * *kt / self.masses[a]).sqrt();
                            let mut g = [0.0; 3];
                            for x in g.iter_mut() {
                                let (u1, u2) = (stream.uniform().max(1e-300), stream.uniform());
                                *x = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
                            }
                            self.vel[i][a] = self.vel[i][a].scale(c) + Vec3 { x: g[0], y: g[1], z: g[2] }.scale(sd);
                        }
                    }
                    self.rattle();
                }
            }
        }
        // Keep each molecule's centre in the box (the molecule whole).
        for i in 0..n {
            let c = self.centre(i);
            let w = c - self.image(c);
            if w.norm2() > 0.0 {
                for a in 0..3 {
                    self.pos[i][a] -= w;
                }
            }
        }
        let warm = std::mem::take(&mut f.dipoles);
        f = self.forces(law, if warm.is_empty() { None } else { Some(&warm) });
        kick(self, &f, 0.5 * h);
        self.rattle();
        f
    }
}
