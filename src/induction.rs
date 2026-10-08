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

/// Run `job` over a range of items split into a fixed number of contiguous
/// chunks, across threads when there is enough work to pay for them, and
/// return the chunks' results in order. The chunks do not depend on the number
/// of threads, so sums over the results are the same on any machine.
fn chunked<T: Send>(len: usize, job: &(dyn Fn(std::ops::Range<usize>) -> T + Sync)) -> Vec<T> {
    const CHUNKS: usize = 32;
    const MIN_WORK: usize = 400;
    if len < MIN_WORK {
        return vec![job(0..len)];
    }
    let chunks = CHUNKS.min(len);
    let parts = crate::electrons::scf::parallel_interleaved(chunks, &|idx: &mut dyn Iterator<Item = usize>| idx.map(|c| (c, job(c * len / chunks..(c + 1) * len / chunks))).collect::<Vec<_>>());
    let mut all: Vec<(usize, T)> = parts.into_iter().flatten().collect();
    all.sort_by_key(|(c, _)| *c);
    all.into_iter().map(|(_, t)| t).collect()
}

/// One screened dipole-dipole coupling between polarisable sites `a` and `b`
/// (flat indices), `r` the separation of `a` from `b`'s image, and the tensor's
/// pieces for it.
#[derive(Clone, Copy)]
struct Dd {
    a: usize,
    b: usize,
    ta: f64,
    tbp: f64,
    da: f64,
    dbp: f64,
    r: Vec3,
}

/// One permanent charge's field on one polarisable site: `r` runs from the
/// charge to the site (the image of whichever is in the second molecule), `h`
/// and `dh` are the field function and its slope, and `pol_on_first` says
/// whether the dipole is on the link's first molecule.
#[derive(Clone, Copy)]
struct Cd {
    a: usize,
    c: usize,
    q: f64,
    h: f64,
    dh: f64,
    r: Vec3,
    pol_on_first: bool,
}

/// Everything one link contributes, worked out once.
#[derive(Default)]
struct LinkItems {
    dd: Vec<Dd>,
    cd: Vec<Cd>,
}

/// Solve a cluster: the dipoles, the energy, the forces. `warm` is the dipoles
/// to start from (the same shape as the result's), `None` for zero.
///
/// The screened tensors are worked out once per link (across threads, in
/// chunks fixed by the number of links) and used by the conjugate-gradient
/// iterations and by the forces alike; the iterations and the forces run over
/// the same chunks, their partial sums added in chunk order.
pub fn solve(cl: &Cluster, links: &[Link], warm: Option<&[Vec<Vec3>]>, tolerance: f64) -> Induced {
    let nm = cl.pol.len();
    // Flat index of molecule i's first polarisable site, and first charge.
    let mut first = vec![0usize; nm + 1];
    let mut firstc = vec![0usize; nm + 1];
    for i in 0..nm {
        first[i + 1] = first[i] + cl.pol[i].len();
        firstc[i + 1] = firstc[i] + cl.charges[i].len();
    }
    let (n, nc) = (first[nm], firstc[nm]);
    let sigma: Vec<Vec<f64>> = cl.pol.iter().map(|m| m.iter().map(|p| sigma_of(p.alpha)).collect()).collect();
    let root2 = (2.0f64).sqrt();
    // Within a molecule.
    let mut intra: Vec<Dd> = Vec::new();
    for i in 0..nm {
        for a in 0..cl.pol[i].len() {
            for b in a + 1..cl.pol[i].len() {
                let r = cl.pol[i][a].pos - cl.pol[i][b].pos;
                let s = (2.0 * (sigma[i][a].powi(2) + sigma[i][b].powi(2))).sqrt();
                let (ta, tbp, da, dbp) = dd_terms(r.norm(), s);
                intra.push(Dd { a: first[i] + a, b: first[i] + b, ta, tbp, da, dbp, r });
            }
        }
    }
    // Between molecules: per link, the pairs of dipoles and the charges' fields.
    let items: Vec<LinkItems> = chunked(links.len(), &|range: std::ops::Range<usize>| -> Vec<LinkItems> {
        range
            .map(|k| {
                let l = &links[k];
                let mut it = LinkItems::default();
                for (a, pa) in cl.pol[l.i].iter().enumerate() {
                    for (b, pb) in cl.pol[l.j].iter().enumerate() {
                        let r = pa.pos - (pb.pos + l.shift);
                        let s = (2.0 * (sigma[l.i][a].powi(2) + sigma[l.j][b].powi(2))).sqrt();
                        let (ta, tbp, da, dbp) = dd_terms(r.norm(), s);
                        it.dd.push(Dd { a: first[l.i] + a, b: first[l.j] + b, ta, tbp, da, dbp, r });
                    }
                    for (ci, c) in cl.charges[l.j].iter().enumerate() {
                        let r = pa.pos - (c.pos + l.shift);
                        let (h, dh) = charge_terms(r.norm(), root2 * sigma[l.i][a]);
                        it.cd.push(Cd { a: first[l.i] + a, c: firstc[l.j] + ci, q: c.q, h, dh, r, pol_on_first: true });
                    }
                }
                for (b, pb) in cl.pol[l.j].iter().enumerate() {
                    for (ci, c) in cl.charges[l.i].iter().enumerate() {
                        let r = (pb.pos + l.shift) - c.pos;
                        let (h, dh) = charge_terms(r.norm(), root2 * sigma[l.j][b]);
                        it.cd.push(Cd { a: first[l.j] + b, c: firstc[l.i] + ci, q: c.q, h, dh, r, pol_on_first: false });
                    }
                }
                it
            })
            .collect()
    })
    .into_iter()
    .flatten()
    .collect();
    let mut e0 = vec![Vec3::ZERO; n];
    for (l, it) in links.iter().zip(&items) {
        for c in &it.cd {
            e0[c.a] += c.r.scale(l.weight * c.q * c.h);
        }
    }
    let inv_alpha: Vec<f64> = cl.pol.iter().flat_map(|m| m.iter().map(|p| 1.0 / p.alpha)).collect();
    let matvec = |x: &[Vec3]| -> Vec<Vec3> {
        let mut y: Vec<Vec3> = x.iter().zip(&inv_alpha).map(|(v, ia)| v.scale(*ia)).collect();
        for d in &intra {
            y[d.a] += tensor_apply(d.ta, d.tbp, d.r, x[d.b]);
            y[d.b] += tensor_apply(d.ta, d.tbp, d.r, x[d.a]);
        }
        let parts = chunked(links.len(), &|range: std::ops::Range<usize>| -> Vec<Vec3> {
            let mut part = vec![Vec3::ZERO; n];
            for k in range {
                let w = links[k].weight;
                for d in &items[k].dd {
                    part[d.a] += tensor_apply(w * d.ta, w * d.tbp, d.r, x[d.b]);
                    part[d.b] += tensor_apply(w * d.ta, w * d.tbp, d.r, x[d.a]);
                }
            }
            part
        });
        for part in parts {
            for (yy, p) in y.iter_mut().zip(part) {
                *yy += p;
            }
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
    let mut fpol = vec![Vec3::ZERO; n];
    let mut fch = vec![Vec3::ZERO; nc];
    // Intramolecular dipole-dipole. Equal and opposite on a molecule's pair of
    // sites, so they add nothing to its net force, but they are not central, so
    // at fixed dipoles they do turn it: the torque of a rigid molecule is the
    // derivative of the minimised energy with respect to its rotation, and
    // those dipoles do not rotate with it until they have been re-minimised.
    for d in &intra {
        let (ma, mb) = (mu[d.a], mu[d.b]);
        let rr = d.r.norm();
        let grad = d.r.scale(d.da / rr * ma.dot(mb)) - d.r.scale(d.dbp / rr * ma.dot(d.r) * mb.dot(d.r)) - ma.scale(d.tbp * mb.dot(d.r)) - mb.scale(d.tbp * ma.dot(d.r));
        fpol[d.a] -= grad;
        fpol[d.b] += grad;
    }
    struct Part {
        fpol: Vec<Vec3>,
        fch: Vec<Vec3>,
        link_force: Vec<Vec3>,
        virial: f64,
    }
    let parts = chunked(links.len(), &|range: std::ops::Range<usize>| -> Part {
        let mut part = Part { fpol: vec![Vec3::ZERO; n], fch: vec![Vec3::ZERO; nc], link_force: Vec::with_capacity(range.len()), virial: 0.0 };
        for k in range {
            let l = &links[k];
            let mut u = 0.0;
            let mut f_i = Vec3::ZERO;
            for d in &items[k].dd {
                let (ma, mb) = (mu[d.a], mu[d.b]);
                let rr = d.r.norm();
                u += ma.dot(mb) * d.ta - d.tbp * ma.dot(d.r) * mb.dot(d.r);
                let grad = d.r.scale(d.da / rr * ma.dot(mb)) - d.r.scale(d.dbp / rr * ma.dot(d.r) * mb.dot(d.r)) - ma.scale(d.tbp * mb.dot(d.r)) - mb.scale(d.tbp * ma.dot(d.r));
                let f = grad.scale(-l.weight);
                part.fpol[d.a] += f;
                part.fpol[d.b] -= f;
                f_i += f;
            }
            for c in &items[k].cd {
                let m = mu[c.a];
                let rr = c.r.norm();
                u -= c.q * c.h * m.dot(c.r);
                let grad = c.r.scale(-c.q * c.dh / rr * m.dot(c.r)) - m.scale(c.q * c.h);
                let f = grad.scale(-l.weight);
                part.fpol[c.a] += f;
                part.fch[c.c] -= f;
                // The force on the first molecule from this pair: what its
                // dipole feels, or minus what the dipole of the other feels
                // (its charge feels the opposite).
                if c.pol_on_first {
                    f_i += f;
                } else {
                    f_i -= f;
                }
            }
            // The weight's own dependence on the separation.
            let dist = l.d.norm();
            let fw = if dist > 0.0 { l.d.scale(-u * l.dweight / dist) } else { Vec3::ZERO };
            part.link_force.push(fw);
            part.virial += l.d.dot(f_i + fw);
        }
        part
    });
    let mut link_force = Vec::with_capacity(links.len());
    let mut virial = 0.0;
    for part in parts {
        for (a, b) in fpol.iter_mut().zip(&part.fpol) {
            *a += *b;
        }
        for (a, b) in fch.iter_mut().zip(&part.fch) {
            *a += *b;
        }
        link_force.extend(part.link_force);
        virial += part.virial;
    }
    let dipoles: Vec<Vec<Vec3>> = (0..nm).map(|i| mu[first[i]..first[i + 1]].to_vec()).collect();
    let force_pol: Vec<Vec<Vec3>> = (0..nm).map(|i| fpol[first[i]..first[i + 1]].to_vec()).collect();
    let force_charge: Vec<Vec<Vec3>> = (0..nm).map(|i| fch[firstc[i]..firstc[i + 1]].to_vec()).collect();
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
