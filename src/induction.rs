//! Induced dipoles between rigid molecules — `docs/PLAY.md` Phase 6, E8a.
//!
//! # Why this is here
//!
//! A pair-additive law cannot make water: the first liquid built from one,
//! fitted to the engine's own Hartree-Fock + MP2 pair energies, has six
//! neighbours inside 3.5 A where water has four and no second shell at 4.5 A
//! (`PLAY.md` E8a). What a pair law lacks is that a molecule's charge cloud
//! answers the field of its neighbours, and the answer changes what the
//! neighbours feel: cooperativity. The simplest physical term that has it is
//! induction: each atom carries a polarisability, and the dipole it takes on
//! minimises its own cost plus its interaction with every other induced dipole
//! and with every other molecule's permanent charges.
//!
//! # The model
//!
//! ```text
//! E(mu) = sum_a mu_a^2 / 2 alpha_a + sum_(a<b) mu_a . T_ab . mu_b - sum_a mu_a . E0_a
//! ```
//!
//! minimised over the dipoles, so `(1/alpha + T) mu = E0` and
//! `E = - 1/2 sum mu . E0`. `E0_a` is the field at atom `a` of the permanent
//! charges of the *other* molecules (a molecule's own charges are already in
//! its own electronic structure), and `T` couples every pair of dipoles,
//! including those of one molecule, which is what makes a molecule's
//! polarisability less than its atoms' sum.
//!
//! **Each induced dipole is a Gaussian-smeared one.** Point dipoles at atomic
//! separations diverge (the polarisation catastrophe); the usual cure, Thole's
//! damping, brings a constant fitted to molecular polarisabilities. Here each
//! dipole is a Gaussian charge cloud, and the interaction of two is the
//! point-dipole one with `1/r` replaced by `erf(r / s) / r`, `s = sqrt(2
//! (sigma_a^2 + sigma_b^2))`; a permanent point charge acts on a smeared dipole
//! through the field of its Gaussian, `q r_vec / r^3 [erf(x) - (2/sqrt(pi)) x
//! e^(-x^2)]`, `x = r / (sqrt(2) sigma_a)`. The width is not a parameter:
//! `sigma = (sqrt(2/pi) alpha / 3)^(1/3)`, the relation under which a smeared
//! dipole has the polarisability it is given (Applequist, Mayer), so the only
//! inputs are the atoms' polarisabilities and the charges the law already has.
//!
//! **In a liquid** the pair interactions are switched on the molecules'
//! centres' separation exactly as the law's are (`liquid::switch`), the
//! dipoles solved by conjugate gradients from the previous step's, and the
//! force is the derivative at fixed dipoles, which is exact at the minimum.

use crate::math::Vec3;

const PI: f64 = std::f64::consts::PI;

/// The error function, by its series `2/sqrt(pi) e^(-x^2) sum 2^n x^(2n+1) /
/// (2n+1)!!` for `|x| < 4` and 1 beyond (`erfc(4)` is 1.5e-8, and the caller's
/// terms there multiply `e^(-x^2)`, 1e-7 of what they would be at 1).
pub fn erf(x: f64) -> f64 {
    let ax = x.abs();
    if ax == 0.0 {
        return 0.0;
    }
    if ax >= 5.5 {
        return x.signum();
    }
    let x2 = ax * ax;
    let mut term = ax;
    let mut sum = ax;
    let mut n = 0.0;
    loop {
        n += 1.0;
        term *= 2.0 * x2 / (2.0 * n + 1.0);
        sum += term;
        if term <= 1e-17 * sum {
            break;
        }
    }
    x.signum() * 2.0 / PI.sqrt() * (-x2).exp() * sum
}

/// The width of a smeared dipole of polarisability `alpha`.
pub fn sigma_of(alpha: f64) -> f64 {
    ((2.0 / PI).sqrt() * alpha / 3.0).cbrt()
}

/// One polarisable site.
#[derive(Debug, Clone, Copy)]
pub struct PolSite {
    pub pos: Vec3,
    pub alpha: f64,
}

/// One permanent point charge.
#[derive(Debug, Clone, Copy)]
pub struct Charge {
    pub pos: Vec3,
    pub q: f64,
}

/// Two molecules' interaction: indices, the shift added to the second's
/// positions (its nearest image), the weight the pair's terms carry and the
/// weight's derivative with respect to the centres' separation `d` (first
/// minus second, nearest image).
#[derive(Debug, Clone, Copy)]
pub struct Link {
    pub i: usize,
    pub j: usize,
    pub shift: Vec3,
    pub weight: f64,
    pub dweight: f64,
    pub d: Vec3,
}

/// The molecules: their polarisable sites and charges, in one frame.
#[derive(Debug, Clone, Default)]
pub struct Cluster {
    pub pol: Vec<Vec<PolSite>>,
    pub charges: Vec<Vec<Charge>>,
}

/// What a solve returns.
#[derive(Debug, Clone)]
pub struct Induced {
    pub energy: f64,
    /// The dipoles, per molecule per polarisable site.
    pub dipoles: Vec<Vec<Vec3>>,
    /// The force on every polarisable site and on every charge, from the
    /// pair terms at fixed dipoles.
    pub force_pol: Vec<Vec<Vec3>>,
    pub force_charge: Vec<Vec<Vec3>>,
    /// Per link, the force on the first molecule's centre from the weight's
    /// dependence on the separation (the second gets minus it).
    pub link_force: Vec<Vec3>,
    /// `sum_links d . F_link`, `F_link` the whole force on the first
    /// molecule from the second in this link.
    pub virial: f64,
    pub iterations: usize,
}

/// `erf(r/s)` and `c = 2/(sqrt(pi) s) e^(-r^2/s^2)` (its derivative).
fn screen(r: f64, s: f64) -> (f64, f64) {
    let x = r / s;
    (erf(x), 2.0 / (PI.sqrt() * s) * (-x * x).exp())
}

/// The screened dipole-dipole tensor `T = A I - Bp r r^T` between two dipoles
/// a distance `r` apart with screening length `s`, and `A'`, `Bp'`.
fn dd_terms(r: f64, s: f64) -> (f64, f64, f64, f64) {
    let (g, c) = screen(r, s);
    let r2 = r * r;
    let r3 = r2 * r;
    let a = g / r3 - c / r2;
    let bp = 3.0 * g / (r3 * r2) - c * (3.0 / (r2 * r2) + 2.0 / (s * s * r2));
    let da = 3.0 * c / r3 - 3.0 * g / (r2 * r2) + 2.0 * c / (s * s * r);
    let dbp = 15.0 * c / (r3 * r2) - 15.0 * g / (r3 * r3) + 10.0 * c / (s * s * r3) + 4.0 * c / (s * s * s * s * r);
    (a, bp, da, dbp)
}

/// The field of a point charge on a smeared dipole of width scale `s_a`:
/// `E = q h(r) r_vec`, `h = [erf(x) - (2/sqrt(pi)) x e^(-x^2)] / r^3`; and `h'`.
fn charge_terms(r: f64, s_a: f64) -> (f64, f64) {
    let x = r / s_a;
    let ex = (-x * x).exp();
    let gg = erf(x) - 2.0 / PI.sqrt() * x * ex;
    let dgg = 4.0 * r * r / (PI.sqrt() * s_a * s_a * s_a) * ex;
    let r3 = r * r * r;
    (gg / r3, dgg / r3 - 3.0 * gg / (r3 * r))
}

fn tensor_apply(a: f64, bp: f64, r: Vec3, v: Vec3) -> Vec3 {
    v.scale(a) - r.scale(bp * r.dot(v))
}

/// Solve a cluster: the dipoles, the energy, the forces. `warm` is the dipoles
/// to start from (the same shape as the result's), `None` for zero.
pub fn solve(cl: &Cluster, links: &[Link], warm: Option<&[Vec<Vec3>]>, tolerance: f64) -> Induced {
    let nm = cl.pol.len();
    // Flat index of molecule i's first polarisable site.
    let mut first = vec![0usize; nm + 1];
    for i in 0..nm {
        first[i + 1] = first[i] + cl.pol[i].len();
    }
    let n = first[nm];
    let sigma: Vec<Vec<f64>> = cl.pol.iter().map(|m| m.iter().map(|p| sigma_of(p.alpha)).collect()).collect();
    // The blocks of T and the right-hand side.
    let mut blocks: Vec<(usize, usize, f64, f64, Vec3)> = Vec::new();
    let mut e0 = vec![Vec3::ZERO; n];
    for i in 0..nm {
        for a in 0..cl.pol[i].len() {
            for b in a + 1..cl.pol[i].len() {
                let r = cl.pol[i][a].pos - cl.pol[i][b].pos;
                let s = (2.0 * (sigma[i][a].powi(2) + sigma[i][b].powi(2))).sqrt();
                let (ta, tbp, _, _) = dd_terms(r.norm(), s);
                blocks.push((first[i] + a, first[i] + b, ta, tbp, r));
            }
        }
    }
    for l in links {
        for (a, pa) in cl.pol[l.i].iter().enumerate() {
            for (b, pb) in cl.pol[l.j].iter().enumerate() {
                let r = pa.pos - (pb.pos + l.shift);
                let s = (2.0 * (sigma[l.i][a].powi(2) + sigma[l.j][b].powi(2))).sqrt();
                let (ta, tbp, _, _) = dd_terms(r.norm(), s);
                blocks.push((first[l.i] + a, first[l.j] + b, l.weight * ta, l.weight * tbp, r));
            }
            for c in &cl.charges[l.j] {
                let r = pa.pos - (c.pos + l.shift);
                let (h, _) = charge_terms(r.norm(), (2.0f64).sqrt() * sigma[l.i][a]);
                e0[first[l.i] + a] += r.scale(l.weight * c.q * h);
            }
        }
        for (b, pb) in cl.pol[l.j].iter().enumerate() {
            for c in &cl.charges[l.i] {
                let r = (pb.pos + l.shift) - c.pos;
                let (h, _) = charge_terms(r.norm(), (2.0f64).sqrt() * sigma[l.j][b]);
                e0[first[l.j] + b] += r.scale(l.weight * c.q * h);
            }
        }
    }
    let inv_alpha: Vec<f64> = cl.pol.iter().flat_map(|m| m.iter().map(|p| 1.0 / p.alpha)).collect();
    let matvec = |x: &[Vec3]| -> Vec<Vec3> {
        let mut y: Vec<Vec3> = x.iter().zip(&inv_alpha).map(|(v, ia)| v.scale(*ia)).collect();
        for &(a, b, ta, tbp, r) in &blocks {
            y[a] += tensor_apply(ta, tbp, r, x[b]);
            y[b] += tensor_apply(ta, tbp, r, x[a]);
        }
        y
    };
    let mut mu: Vec<Vec3> = match warm {
        Some(w) if w.iter().map(|m| m.len()).sum::<usize>() == n => w.iter().flatten().cloned().collect(),
        _ => vec![Vec3::ZERO; n],
    };
    // Conjugate gradients, Jacobi-preconditioned.
    let ax = matvec(&mu);
    let mut r: Vec<Vec3> = e0.iter().zip(&ax).map(|(b, a)| *b - *a).collect();
    let mut z: Vec<Vec3> = r.iter().zip(&inv_alpha).map(|(v, ia)| v.scale(1.0 / ia)).collect();
    let mut p = z.clone();
    let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a.dot(*b)).sum();
    let bnorm = e0.iter().map(|v| v.norm2()).sum::<f64>().sqrt().max(1e-300);
    let mut iterations = 0;
    for it in 0..200 {
        let rn = r.iter().map(|v| v.norm2()).sum::<f64>().sqrt();
        if rn <= tolerance * bnorm {
            break;
        }
        iterations = it + 1;
        let ap = matvec(&p);
        let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a.dot(*b)).sum();
        if pap <= 0.0 {
            break;
        }
        let alpha = rz / pap;
        for k in 0..n {
            mu[k] += p[k].scale(alpha);
            r[k] -= ap[k].scale(alpha);
        }
        z = r.iter().zip(&inv_alpha).map(|(v, ia)| v.scale(1.0 / ia)).collect();
        let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a.dot(*b)).sum();
        let beta = rz_new / rz;
        rz = rz_new;
        for k in 0..n {
            p[k] = z[k] + p[k].scale(beta);
        }
    }
    let energy = -0.5 * mu.iter().zip(&e0).map(|(m, e)| m.dot(*e)).sum::<f64>();
    // Forces at fixed dipoles.
    let mut force_pol: Vec<Vec<Vec3>> = cl.pol.iter().map(|m| vec![Vec3::ZERO; m.len()]).collect();
    let mut force_charge: Vec<Vec<Vec3>> = cl.charges.iter().map(|m| vec![Vec3::ZERO; m.len()]).collect();
    let mut link_force = vec![Vec3::ZERO; links.len()];
    let mut virial = 0.0;
    // Intramolecular dipole-dipole. Equal and opposite on a molecule's pair of
    // sites, so they add nothing to its net force, but they are not central, so
    // at fixed dipoles they do turn it: the torque of a rigid molecule is the
    // derivative of the minimised energy with respect to its rotation, and
    // those dipoles do not rotate with it until they have been re-minimised.
    for i in 0..nm {
        for a in 0..cl.pol[i].len() {
            for b in a + 1..cl.pol[i].len() {
                let (ma, mb) = (mu[first[i] + a], mu[first[i] + b]);
                let r = cl.pol[i][a].pos - cl.pol[i][b].pos;
                let rr = r.norm();
                let s = (2.0 * (sigma[i][a].powi(2) + sigma[i][b].powi(2))).sqrt();
                let (_, tbp, da, dbp) = dd_terms(rr, s);
                let grad = r.scale(da / rr * ma.dot(mb)) - r.scale(dbp / rr * ma.dot(r) * mb.dot(r)) - ma.scale(tbp * mb.dot(r)) - mb.scale(tbp * ma.dot(r));
                force_pol[i][a] -= grad;
                force_pol[i][b] += grad;
            }
        }
    }
    for (k, l) in links.iter().enumerate() {
        let mut u = 0.0;
        // Forces on i's sites from j, and on j's from i, before the weight.
        let mut f_i = Vec3::ZERO;
        for (a, pa) in cl.pol[l.i].iter().enumerate() {
            let ma = mu[first[l.i] + a];
            for (b, pb) in cl.pol[l.j].iter().enumerate() {
                let mb = mu[first[l.j] + b];
                let r = pa.pos - (pb.pos + l.shift);
                let rr = r.norm();
                let s = (2.0 * (sigma[l.i][a].powi(2) + sigma[l.j][b].powi(2))).sqrt();
                let (ta, tbp, da, dbp) = dd_terms(rr, s);
                u += ma.dot(mb) * ta - tbp * ma.dot(r) * mb.dot(r);
                // d/dr_vec of u.
                let grad = r.scale(da / rr * ma.dot(mb)) - r.scale(dbp / rr * ma.dot(r) * mb.dot(r)) - ma.scale(tbp * mb.dot(r)) - mb.scale(tbp * ma.dot(r));
                let f = grad.scale(-l.weight);
                force_pol[l.i][a] += f;
                force_pol[l.j][b] -= f;
                f_i += f;
            }
            for (ci, c) in cl.charges[l.j].iter().enumerate() {
                let r = pa.pos - (c.pos + l.shift);
                let rr = r.norm();
                let (h, dh) = charge_terms(rr, (2.0f64).sqrt() * sigma[l.i][a]);
                u -= c.q * h * ma.dot(r);
                let grad = r.scale(-c.q * dh / rr * ma.dot(r)) - ma.scale(c.q * h);
                let f = grad.scale(-l.weight);
                force_pol[l.i][a] += f;
                force_charge[l.j][ci] -= f;
                f_i += f;
            }
        }
        for (b, pb) in cl.pol[l.j].iter().enumerate() {
            let mb = mu[first[l.j] + b];
            for (ci, c) in cl.charges[l.i].iter().enumerate() {
                let r = (pb.pos + l.shift) - c.pos;
                let rr = r.norm();
                let (h, dh) = charge_terms(rr, (2.0f64).sqrt() * sigma[l.j][b]);
                u -= c.q * h * mb.dot(r);
                let grad = r.scale(-c.q * dh / rr * mb.dot(r)) - mb.scale(c.q * h);
                let f = grad.scale(-l.weight);
                force_pol[l.j][b] += f;
                force_charge[l.i][ci] -= f;
                // The force on molecule i from this pair is what its charge
                // feels: minus what the dipole feels.
                f_i -= f;
            }
        }
        // The weight's own dependence on the separation.
        let dist = l.d.norm();
        let fw = if dist > 0.0 { l.d.scale(-u * l.dweight / dist) } else { Vec3::ZERO };
        link_force[k] = fw;
        virial += l.d.dot(f_i + fw);
    }
    let dipoles: Vec<Vec<Vec3>> = (0..nm).map(|i| mu[first[i]..first[i + 1]].to_vec()).collect();
    Induced { energy, dipoles, force_pol, force_charge, link_force, virial, iterations }
}

/// The polarisability tensor of one molecule from its atoms' polarisabilities,
/// by this model: the response of the summed dipole to a uniform field, with
/// the atoms' dipole-dipole coupling and their widths, nothing else.
pub fn molecular_polarisability(sites: &[PolSite]) -> [[f64; 3]; 3] {
    let cl = Cluster { pol: vec![sites.to_vec()], charges: vec![Vec::new()] };
    let n = sites.len();
    let sigma: Vec<f64> = sites.iter().map(|p| sigma_of(p.alpha)).collect();
    // (1/alpha + T) mu = E, by solving the 3n x 3n system directly.
    let dim = 3 * n;
    let mut m = vec![0.0f64; dim * dim];
    for a in 0..n {
        for d in 0..3 {
            m[(3 * a + d) * dim + 3 * a + d] = 1.0 / sites[a].alpha;
        }
        for b in 0..n {
            if a == b {
                continue;
            }
            let r = sites[a].pos - sites[b].pos;
            let s = (2.0 * (sigma[a].powi(2) + sigma[b].powi(2))).sqrt();
            let (ta, tbp, _, _) = dd_terms(r.norm(), s);
            let rv = [r.x, r.y, r.z];
            for d in 0..3 {
                for e in 0..3 {
                    m[(3 * a + d) * dim + 3 * b + e] += (if d == e { ta } else { 0.0 }) - tbp * rv[d] * rv[e];
                }
            }
        }
    }
    let _ = cl;
    let mut alpha = [[0.0f64; 3]; 3];
    for e in 0..3 {
        let mut rhs = vec![0.0f64; dim];
        for a in 0..n {
            rhs[3 * a + e] = 1.0;
        }
        let x = solve_dense(&mut m.clone(), &mut rhs, dim);
        for d in 0..3 {
            alpha[d][e] = (0..n).map(|a| x[3 * a + d]).sum();
        }
    }
    alpha
}

/// Gaussian elimination with partial pivoting on a small dense system.
fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> Vec<f64> {
    for col in 0..n {
        let piv = (col..n).max_by(|&x, &y| a[x * n + col].abs().total_cmp(&a[y * n + col].abs())).expect("a pivot");
        if piv != col {
            for k in 0..n {
                a.swap(piv * n + k, col * n + k);
            }
            b.swap(piv, col);
        }
        for r in col + 1..n {
            let f = a[r * n + col] / a[col * n + col];
            for k in col..n {
                a[r * n + k] -= f * a[col * n + k];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in r + 1..n {
            s -= a[r * n + k] * x[k];
        }
        x[r] = s / a[r * n + r];
    }
    x
}
