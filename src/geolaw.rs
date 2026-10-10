//! A water law whose cloud follows the molecule's angle — `docs/PLAY.md` E8c,
//! treatment (2), and the first form of the geometry-following law.
//!
//! The core-and-cloud model (`phys-esp --cloud`): the oxygen and hydrogen nuclei
//! less their 1s electrons as point charges (+6, +1), and the eight valence
//! electrons as Gaussian clouds on the oxygen, on each hydrogen and on the
//! bisector of the oxygen's bonds. Its six numbers are the three clouds' widths
//! (as logarithms), the bisector cloud's distance from the oxygen, and the
//! populations of the bisector and the oxygen clouds (the hydrogens' follow from
//! the eight electrons). Here each is **a quadratic in the molecule's H-O-H
//! angle**, `x_p(theta) = c0 + c1 u + c2 u^2`, `u = (theta - theta0) / scale`:
//! the monomer-derived terms of E8c's hybrid law, fitted to the exact
//! electrostatic energies of pairs whose molecules each have their own angle.
//!
//! What does not follow the angle, by choice: the repulsion's overlap form uses
//! the *base* widths (`SiteSite::sigma`), so that the repulsion fitted to a
//! pair's remaining energy does not depend on the cloud being fitted; and the
//! polarisabilities are the atoms'.
//!
//! The energy of a liquid of such molecules depends on each molecule's angle
//! twice: through the atoms' positions, and through these parameters. The
//! derivative through the parameters is [`GeoCloud::dparams`] times the energy's
//! derivative with respect to each charge and width, which the pair loop and
//! `induction::charge_gradients` supply; `flex.rs` assembles the force.

use crate::induction::{solve, Charge, Cluster, Link, PolSite};
use crate::liquid::{extra_sites_from_text, with_bisector_site, with_frame_sites, ExtraSite, SiteSite};
use crate::math::Vec3;

/// Site types of the core-and-cloud law: 0 and 1 the cores, 2 the bisector
/// cloud, 3 the oxygen's, 4 the hydrogens'.
pub const TYPES: usize = 5;
/// Parameters of the cloud: ln width of the bisector, oxygen and hydrogen
/// clouds, the bisector cloud's distance, the bisector and oxygen populations.
pub const NPAR: usize = 6;

/// One molecule's cloud at one angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CloudParams {
    /// A charge per site type (the cores', which do not move, included).
    pub charge: [f64; TYPES],
    /// A width per site type (zero for the cores).
    pub sigma: [f64; TYPES],
    /// The bisector cloud's distance from the oxygen, bohr.
    pub d: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GeoCloud {
    /// The core charges (types 0 and 1), which are the elements' and do not move.
    pub cores: [f64; 2],
    pub theta0: f64,
    pub scale: f64,
    /// `coef[p]` the quadratic of parameter `p`.
    pub coef: [[f64; 3]; NPAR],
}

impl GeoCloud {
    /// A cloud that does not depend on the angle: the numbers `x`, as
    /// `phys-diag --fit-es` writes them.
    pub fn constant(cores: [f64; 2], x: [f64; NPAR], theta0: f64) -> GeoCloud {
        let mut coef = [[0.0; 3]; NPAR];
        for p in 0..NPAR {
            coef[p][0] = x[p];
        }
        GeoCloud { cores, theta0, scale: 0.2, coef }
    }

    pub fn x(&self, theta: f64) -> [f64; NPAR] {
        let u = (theta - self.theta0) / self.scale;
        let mut out = [0.0; NPAR];
        for p in 0..NPAR {
            out[p] = self.coef[p][0] + self.coef[p][1] * u + self.coef[p][2] * u * u;
        }
        out
    }

    /// `d x / d theta`.
    pub fn dx(&self, theta: f64) -> [f64; NPAR] {
        let u = (theta - self.theta0) / self.scale;
        let mut out = [0.0; NPAR];
        for p in 0..NPAR {
            out[p] = (self.coef[p][1] + 2.0 * self.coef[p][2] * u) / self.scale;
        }
        out
    }

    fn params_of(&self, x: &[f64; NPAR]) -> CloudParams {
        let (n_m, n_o) = (x[4], x[5]);
        let n_h = (8.0 - n_m - n_o) / 2.0;
        CloudParams { charge: [self.cores[0], self.cores[1], -n_m, -n_o, -n_h], sigma: [0.0, 0.0, x[0].exp(), x[1].exp(), x[2].exp()], d: x[3].abs().max(1e-3) }
    }

    pub fn params(&self, theta: f64) -> CloudParams {
        self.params_of(&self.x(theta))
    }

    /// The derivative of every parameter with respect to the angle.
    pub fn dparams(&self, theta: f64) -> CloudParams {
        let x = self.x(theta);
        let dx = self.dx(theta);
        let sigma = |p: usize| x[p].exp() * dx[p];
        CloudParams {
            charge: [0.0, 0.0, -dx[4], -dx[5], 0.5 * (dx[4] + dx[5])],
            sigma: [0.0, 0.0, sigma(0), sigma(1), sigma(2)],
            d: if x[3].abs() > 1e-3 { dx[3] * x[3].signum() } else { 0.0 },
        }
    }

    /// The text of a law that has one: `geo <theta0> <scale> <core O> <core H>`
    /// and a `geopoly <p> <c0> <c1> <c2>` line for each parameter.
    pub fn to_text(&self) -> String {
        let mut t = format!("geo {:e} {:e} {:e} {:e}\n", self.theta0, self.scale, self.cores[0], self.cores[1]);
        for p in 0..NPAR {
            t += &format!("geopoly {p} {:e} {:e} {:e}\n", self.coef[p][0], self.coef[p][1], self.coef[p][2]);
        }
        t
    }

    pub fn from_text(text: &str) -> Option<GeoCloud> {
        let mut head: Option<(f64, f64, [f64; 2])> = None;
        let mut coef = [[0.0; 3]; NPAR];
        let mut seen = 0;
        for l in text.lines() {
            let w: Vec<&str> = l.split_whitespace().collect();
            match w.first().copied() {
                Some("geo") if w.len() >= 5 => head = Some((w[1].parse().ok()?, w[2].parse().ok()?, [w[3].parse().ok()?, w[4].parse().ok()?])),
                Some("geopoly") if w.len() >= 5 => {
                    let p: usize = w[1].parse().ok()?;
                    for k in 0..3 {
                        coef[p][k] = w[2 + k].parse().ok()?;
                    }
                    seen += 1;
                }
                _ => {}
            }
        }
        let (theta0, scale, cores) = head?;
        (seen == NPAR).then_some(GeoCloud { cores, theta0, scale, coef })
    }

    /// The oxygen-centred angle of a molecule's three atoms.
    pub fn angle(pos: &[Vec3; 3]) -> f64 {
        let (a, b) = ((pos[1] - pos[0]).unit(), (pos[2] - pos[0]).unit());
        a.dot(b).clamp(-1.0, 1.0).acos()
    }
}

/// A molecule's sites for the pair energy: its three atoms (types 0, 1, 1),
/// the bisector cloud at `d`, and the frame sites (the oxygen cloud on the
/// oxygen, a hydrogen cloud on each hydrogen at the bond length).
pub fn placed(atoms: &[(Vec3, usize)], d: f64, frame: &[ExtraSite]) -> Vec<(Vec3, usize)> {
    with_frame_sites(&with_bisector_site(atoms, 2, d), frame)
}

/// The electrostatic energy of two molecules' clouds (hartree), site against
/// site as the law has it: `q_a q_b erf(r / s) / r`, `s^2 = sigma_a^2 + sigma_b^2`.
pub fn pair_electrostatics(sa: &[(Vec3, usize)], pa: &CloudParams, sb: &[(Vec3, usize)], pb: &CloudParams) -> f64 {
    let mut e = 0.0;
    for (xa, ta) in sa {
        for (xb, tb) in sb {
            let qq = pa.charge[*ta] * pb.charge[*tb];
            if qq == 0.0 {
                continue;
            }
            let r = (*xa - *xb).norm();
            let s2 = pa.sigma[*ta].powi(2) + pb.sigma[*tb].powi(2);
            e += if s2 > 0.0 { qq * crate::induction::erf(r / s2.sqrt()) / r } else { qq / r };
        }
    }
    e
}

/// The induction energy of two molecules (hartree): each atom a polarisable
/// site of `alpha[type]`, each site with a charge a source, by the same solve a
/// liquid uses.
pub fn pair_induction(alpha: &[f64], sa: &[(Vec3, usize)], pa: &CloudParams, sb: &[(Vec3, usize)], pb: &CloudParams) -> f64 {
    let build = |s: &[(Vec3, usize)], p: &CloudParams| -> (Vec<PolSite>, Vec<Charge>) {
        let pol = s.iter().filter(|(_, t)| alpha.get(*t).copied().unwrap_or(0.0) > 0.0 && *t < 2).map(|(x, t)| PolSite { pos: *x, alpha: alpha[*t] }).collect();
        let ch = s.iter().filter(|(_, t)| p.charge[*t] != 0.0).map(|(x, t)| Charge { pos: *x, q: p.charge[*t], sigma: p.sigma[*t] }).collect();
        (pol, ch)
    };
    let (pol_a, ch_a) = build(sa, pa);
    let (pol_b, ch_b) = build(sb, pb);
    let cl = Cluster { pol: vec![pol_a, pol_b], charges: vec![ch_a, ch_b] };
    let centre = |s: &[(Vec3, usize)]| s.iter().take(3).fold(Vec3::ZERO, |acc, (x, t)| acc + x.scale(if *t == 0 { 15.999 } else { 1.008 })).scale(1.0 / 18.015);
    let d = centre(sa) - centre(sb);
    let link = Link { i: 0, j: 1, shift: Vec3::ZERO, weight: 1.0, dweight: 0.0, d };
    solve(&cl, &[link], None, 1e-12).energy
}

/// A pair's electrostatics and induction with each molecule's own cloud: the
/// part of its energy that follows the angles. `a` and `b` are the molecules'
/// atoms (types 0, 1, 1).
pub fn pair_es_ind(geo: &GeoCloud, frame: &[ExtraSite], alpha: &[f64], a: &[(Vec3, usize)], b: &[(Vec3, usize)], follow: bool) -> f64 {
    let angle = |m: &[(Vec3, usize)]| if follow { GeoCloud::angle(&[m[0].0, m[1].0, m[2].0]) } else { geo.theta0 };
    let (pa, pb) = (geo.params(angle(a)), geo.params(angle(b)));
    let (sa, sb) = (placed(a, pa.d, frame), placed(b, pb.d, frame));
    pair_electrostatics(&sa, &pa, &sb, &pb) + pair_induction(alpha, &sa, &pa, &sb, &pb)
}

/// The correction a fit of the pair terms applies to a pair's energy: the part
/// of electrostatics and induction that the angles change from their value at
/// the reference angle, so that a fitter whose model has the reference cloud
/// finds the remaining terms. `E_target - correction` is what it fits.
pub fn correction(geo: &GeoCloud, frame: &[ExtraSite], alpha: &[f64], a: &[(Vec3, usize)], b: &[(Vec3, usize)]) -> f64 {
    pair_es_ind(geo, frame, alpha, a, b, true) - pair_es_ind(geo, frame, alpha, a, b, false)
}

/// The geo cloud's law text for a program that reads a law: the base law
/// (types 0..5, the cloud at the reference angle), the frame sites and the
/// bisector line, then the `geo` lines.
pub fn read_frame(text: &str) -> Vec<ExtraSite> {
    extra_sites_from_text(text)
}

/// Reference cloud of a law file as the numbers `phys-diag --fit-es` fits.
pub fn x_of_law(law: &SiteSite, bisector_distance: f64) -> [f64; NPAR] {
    [law.sigma[2].ln(), law.sigma[3].ln(), law.sigma[4].ln(), bisector_distance, -law.charge[2], -law.charge[3]]
}
