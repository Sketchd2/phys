//! Ewald sums for the liquid's Gaussian charges and induced dipoles
//! (`docs/PLAY.md` E8d).
//!
//! # Why
//!
//! `phys-bulk` cuts every pair term off at 16 bohr with a switch on the
//! molecules' centres. Measured on the same configuration with the law fixed
//! (`runs/many/E8d-notes.md` 6.8) the pressure changes by thousands of bar with
//! the cut-off and is still doing so at 32 bohr: the switched sum of dipolar
//! interactions does not converge, and every pressure the project has
//! reported carries that error. This is the remedy that is standard and
//! does not depend on the law's numbers.
//!
//! # The split
//!
//! Every object the liquid carries - a permanent charge, an induced dipole - is
//! a Gaussian of some width `w`, and any two of them interact through
//! derivatives of `erf(r / s) / r` with `s^2 = w_a^2 + w_b^2` (charge-charge
//! `sigma_a^2 + sigma_b^2`, charge-dipole `2 sigma_d^2 + sigma_c^2`,
//! dipole-dipole `2 sigma_a^2 + 2 sigma_b^2`; `induction.rs`). One **common**
//! width `s' >= s` for every pair splits each of them as
//!
//! ```text
//! [ kernel(s) - kernel(s') ] + kernel(s')
//! ```
//!
//! The first bracket decays like `erfc(r / s')` and is summed in real space over
//! the nearest images, using the very functions `induction.rs` already has
//! (evaluated at `s` and at `s'` and differenced). The second is the same
//! operator on `erf(r / s') / r` for every pair whatever its own widths, so it
//! is the interaction of **point** multipoles smoothed by one Gaussian, and its
//! Fourier transform is `4 pi / k^2 exp(-k^2 s'^2 / 4)`:
//!
//! ```text
//! E_rec = (2 pi / V) sum_(k != 0) exp(-k^2 s'^2 / 4) / k^2 |S(k)|^2,
//! S(k)  = sum_i q_i e^(i k.r_i) + sum_a i (k.mu_a) e^(i k.r_a).
//! ```
//!
//! It includes each object's interaction with itself and with the other objects
//! of its own molecule, which the liquid does not have (the permanent charges
//! of one molecule do not act on its own dipoles, and the charges do not act on
//! one another), so those are subtracted; the dipoles' own mutual coupling
//! inside a molecule *is* in the model, and takes its short-range part in real
//! space like any other pair. Tin-foil boundaries (no `k = 0` term).
//!
//! The reciprocal virial is taken as the derivative of `E_rec` under rigid
//! scaling of the molecules' centres (their shapes and the dipoles fixed),
//! which is what the molecular virial of a liquid of rigid molecules is.

use crate::math::Vec3;

const PI: f64 = std::f64::consts::PI;

/// How many threads the reciprocal loops use. Each thread owns a disjoint slice
/// of the output and adds in a fixed order, so the numbers do not depend on it.
fn threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).min(12)
}

/// A complex number, as little as this needs.
#[derive(Clone, Copy, Debug, Default)]
pub struct C {
    pub re: f64,
    pub im: f64,
}

impl C {
    fn mul(self, o: C) -> C {
        C { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
    fn conj(self) -> C {
        C { re: self.re, im: -self.im }
    }
}

/// The reciprocal-space machinery for one box and one splitting width.
#[derive(Clone, Debug)]
pub struct Ewald {
    /// The common splitting width `s'` (bohr).
    pub s: f64,
    pub cell: [f64; 3],
    pub volume: f64,
    nmax: [i32; 3],
    /// The integer vectors of the half space (`k` and `-k` count once).
    idx: Vec<[i32; 3]>,
    kv: Vec<[f64; 3]>,
    fk: Vec<f64>,
    /// `k` at which `exp(-k^2 s'^2 / 4)` falls to the accuracy asked.
    pub kcut: f64,
}

impl Ewald {
    /// The vectors out to where the Gaussian factor is `accuracy` (1e-7 is
    /// ample: the sum's error is of that order of the energy scale).
    pub fn new(s: f64, cell: [f64; 3], accuracy: f64) -> Ewald {
        let kcut = 2.0 * (1.0 / accuracy).ln().sqrt() / s;
        let nmax = [0, 1, 2].map(|a| (kcut * cell[a] / (2.0 * PI)).ceil() as i32);
        let mut idx = Vec::new();
        for nx in 0..=nmax[0] {
            for ny in -nmax[1]..=nmax[1] {
                for nz in -nmax[2]..=nmax[2] {
                    if nx == 0 && (ny < 0 || (ny == 0 && nz <= 0)) {
                        continue;
                    }
                    let k = [2.0 * PI * nx as f64 / cell[0], 2.0 * PI * ny as f64 / cell[1], 2.0 * PI * nz as f64 / cell[2]];
                    if k[0] * k[0] + k[1] * k[1] + k[2] * k[2] <= kcut * kcut {
                        idx.push([nx, ny, nz]);
                    }
                }
            }
        }
        let (kv, fk) = Self::kdata(&idx, s, cell, 1.0);
        Ewald { s, cell, volume: cell[0] * cell[1] * cell[2], nmax, idx, kv, fk, kcut }
    }

    fn kdata(idx: &[[i32; 3]], s: f64, cell: [f64; 3], lambda: f64) -> (Vec<[f64; 3]>, Vec<f64>) {
        let mut kv = Vec::with_capacity(idx.len());
        let mut fk = Vec::with_capacity(idx.len());
        for n in idx {
            let k = [2.0 * PI * n[0] as f64 / (cell[0] * lambda), 2.0 * PI * n[1] as f64 / (cell[1] * lambda), 2.0 * PI * n[2] as f64 / (cell[2] * lambda)];
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            kv.push(k);
            fk.push((-k2 * s * s / 4.0).exp() / k2);
        }
        (kv, fk)
    }

    /// How many `k` vectors.
    pub fn nk(&self) -> usize {
        self.idx.len()
    }

    /// `e^(i k.r)` for every site and every `k`, site-major. `lambda` scales the
    /// box (and so every `k`) for the virial; positions are taken as given.
    pub fn phases(&self, pos: &[Vec3]) -> Vec<C> {
        let nk = self.idx.len();
        let mut out = vec![C::default(); pos.len() * nk];
        if pos.is_empty() {
            return out;
        }
        let n = [self.nmax[0] as usize, self.nmax[1] as usize, self.nmax[2] as usize];
        let per = pos.len().div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in out.chunks_mut(per * nk).enumerate() {
                sc.spawn(move || {
                    for (si, row) in chunk.chunks_mut(nk).enumerate() {
                        let p = pos[ci * per + si];
                        // e^(i n theta) for n = 0..=nmax on each axis; negative n by conjugate.
                        let table = |x: f64, a: usize| -> Vec<C> {
                            let th = 2.0 * PI * x / self.cell[a];
                            let step = C { re: th.cos(), im: th.sin() };
                            let mut t = Vec::with_capacity(n[a] + 1);
                            let mut cur = C { re: 1.0, im: 0.0 };
                            for m in 0..=n[a] {
                                // Recompute exactly now and then so the error does not creep.
                                if m % 8 == 0 {
                                    let a2 = th * m as f64;
                                    cur = C { re: a2.cos(), im: a2.sin() };
                                }
                                t.push(cur);
                                cur = cur.mul(step);
                            }
                            t
                        };
                        let (tx, ty, tz) = (table(p.x, 0), table(p.y, 1), table(p.z, 2));
                        let get = |t: &Vec<C>, m: i32| if m >= 0 { t[m as usize] } else { t[(-m) as usize].conj() };
                        for (j, m) in self.idx.iter().enumerate() {
                            row[j] = get(&tx, m[0]).mul(get(&ty, m[1])).mul(get(&tz, m[2]));
                        }
                    }
                });
            }
        });
        out
    }

    /// `S_q(k) = sum q e^(ik.r)`.
    pub fn charge_sum(&self, ph: &[C], q: &[f64]) -> Vec<C> {
        let nk = self.idx.len();
        let mut s = vec![C::default(); nk];
        let per = nk.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in s.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    let k0 = ci * per;
                    for (i, qi) in q.iter().enumerate() {
                        let row = &ph[i * nk + k0..i * nk + k0 + chunk.len()];
                        for (sk, p) in chunk.iter_mut().zip(row) {
                            sk.re += qi * p.re;
                            sk.im += qi * p.im;
                        }
                    }
                });
            }
        });
        s
    }

    /// `D(k) = sum i (k.mu) e^(ik.r)`.
    pub fn dipole_sum(&self, ph: &[C], mu: &[Vec3]) -> Vec<C> {
        let nk = self.idx.len();
        let mut s = vec![C::default(); nk];
        let per = nk.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in s.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    let k0 = ci * per;
                    for (a, m) in mu.iter().enumerate() {
                        let row = &ph[a * nk + k0..a * nk + k0 + chunk.len()];
                        for (j, (sk, p)) in chunk.iter_mut().zip(row).enumerate() {
                            let k = self.kv[k0 + j];
                            let km = k[0] * m.x + k[1] * m.y + k[2] * m.z;
                            // i (k.mu) (p.re + i p.im) = (k.mu)(-p.im + i p.re)
                            sk.re += -km * p.im;
                            sk.im += km * p.re;
                        }
                    }
                });
            }
        });
        s
    }

    /// `(2 pi / V) sum_(all k) f |S|^2 = (4 pi / V) sum_(half) f |S|^2`.
    pub fn energy(&self, s: &[C]) -> f64 {
        let mut e = 0.0;
        for (sk, f) in s.iter().zip(&self.fk) {
            e += f * (sk.re * sk.re + sk.im * sk.im);
        }
        4.0 * PI / self.volume * e
    }

    /// `dE/dmu_a = (8 pi / V) sum f Re[ S* i k e^(ik.r_a) ]` for every dipole.
    pub fn grad_dipole(&self, s: &[C], ph: &[C], n: usize) -> Vec<Vec3> {
        let nk = self.idx.len();
        let pref = 8.0 * PI / self.volume;
        let mut out = vec![Vec3::ZERO; n];
        if n == 0 {
            return out;
        }
        let per = n.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in out.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    for (ai, o) in chunk.iter_mut().enumerate() {
                        let a = ci * per + ai;
                        let row = &ph[a * nk..(a + 1) * nk];
                        let (mut gx, mut gy, mut gz) = (0.0, 0.0, 0.0);
                        for j in 0..nk {
                            // S* i e^(ikr): S* = (sr - i si), i e = i(pr + i pi) = (-pi + i pr)
                            let (sr, si) = (s[j].re, s[j].im);
                            let (pr, pi) = (row[j].re, row[j].im);
                            let re = sr * (-pi) + si * pr; // Re[(sr - i si)(-pi + i pr)]
                            let w = self.fk[j] * re;
                            gx += w * self.kv[j][0];
                            gy += w * self.kv[j][1];
                            gz += w * self.kv[j][2];
                        }
                        *o = Vec3 { x: gx * pref, y: gy * pref, z: gz * pref };
                    }
                });
            }
        });
        out
    }

    /// `dE/dr_i` of the charges: `(8 pi / V) sum f Re[ S* i k q_i e^(ik.r_i) ]`.
    pub fn grad_charge_pos(&self, s: &[C], ph: &[C], q: &[f64]) -> Vec<Vec3> {
        let g = self.grad_dipole(s, ph, q.len());
        g.into_iter().zip(q).map(|(v, qi)| v.scale(*qi)).collect()
    }

    /// `dE/dr_a` of the dipoles: `(8 pi / V) sum f Re[ -S* (k.mu_a) k e^(ik.r_a) ]`.
    pub fn grad_dipole_pos(&self, s: &[C], ph: &[C], mu: &[Vec3]) -> Vec<Vec3> {
        let nk = self.idx.len();
        let pref = 8.0 * PI / self.volume;
        let mut out = vec![Vec3::ZERO; mu.len()];
        if mu.is_empty() {
            return out;
        }
        let per = mu.len().div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in out.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    for (ai, o) in chunk.iter_mut().enumerate() {
                        let a = ci * per + ai;
                        let m = mu[a];
                        let row = &ph[a * nk..(a + 1) * nk];
                        let (mut gx, mut gy, mut gz) = (0.0, 0.0, 0.0);
                        for j in 0..nk {
                            let k = self.kv[j];
                            let km = k[0] * m.x + k[1] * m.y + k[2] * m.z;
                            // Re[ -(sr - i si)(pr + i pi) ] = -(sr pr + si pi)
                            let re = -(s[j].re * row[j].re + s[j].im * row[j].im);
                            let w = self.fk[j] * re * km;
                            gx += w * k[0];
                            gy += w * k[1];
                            gz += w * k[2];
                        }
                        *o = Vec3 { x: gx * pref, y: gy * pref, z: gz * pref };
                    }
                });
            }
        });
        out
    }

    /// The long-range field of the permanent charges at the dipoles' sites,
    /// `E0_LR = -dE_(q,mu)/dmu = -(8 pi / V) sum f Re[ S_q* i k e^(ik.r_a) ]`.
    pub fn field_at_dipoles(&self, sq: &[C], ph_dip: &[C], n: usize) -> Vec<Vec3> {
        self.grad_dipole(sq, ph_dip, n).into_iter().map(|v| -v).collect()
    }

    /// The long-range dipole-dipole coupling `T_LR mu`, less each dipole's
    /// coupling to itself (`+4 / (3 sqrt(pi) s'^3)` too large in `T`, since the
    /// reciprocal sum includes the pair `(a, a)`).
    pub fn dipole_response(&self, ph_dip: &[C], mu: &[Vec3]) -> Vec<Vec3> {
        let d = self.dipole_sum(ph_dip, mu);
        let g = self.grad_dipole(&d, ph_dip, mu.len());
        let own = 4.0 / (3.0 * PI.sqrt() * self.s.powi(3));
        g.into_iter().zip(mu).map(|(v, m)| v - m.scale(own)).collect()
    }

    /// `E_rec` with the box scaled by `lambda` and every molecule's centre with
    /// it, shapes and dipoles fixed: the sites are given as their molecule's
    /// centre and their arm from it.
    pub fn energy_scaled(&self, lambda: f64, ccom: &[Vec3], carm: &[Vec3], q: &[f64], dcom: &[Vec3], darm: &[Vec3], mu: &[Vec3]) -> f64 {
        // k.r with k -> k / lambda and r -> com lambda + arm is k.com + (k/lambda).arm:
        // the phase of (com + arm / lambda) in the unscaled box.
        let pc: Vec<Vec3> = ccom.iter().zip(carm).map(|(c, a)| *c + a.scale(1.0 / lambda)).collect();
        let pd: Vec<Vec3> = dcom.iter().zip(darm).map(|(c, a)| *c + a.scale(1.0 / lambda)).collect();
        let (kv, fk) = Self::kdata(&self.idx, self.s, self.cell, lambda);
        let phc = self.phases(&pc);
        let nk = self.idx.len();
        let mut s = self.charge_sum(&phc, q);
        if !mu.is_empty() {
            let phd = self.phases(&pd);
            for (a, m) in mu.iter().enumerate() {
                let row = &phd[a * nk..(a + 1) * nk];
                for j in 0..nk {
                    let k = kv[j];
                    let km = k[0] * m.x + k[1] * m.y + k[2] * m.z;
                    s[j].re += -km * row[j].im;
                    s[j].im += km * row[j].re;
                }
            }
        }
        let mut e = 0.0;
        for (sk, f) in s.iter().zip(&fk) {
            e += f * (sk.re * sk.re + sk.im * sk.im);
        }
        4.0 * PI / (self.volume * lambda.powi(3)) * e
    }

    /// `-dE_rec/d ln(lambda)`, the reciprocal part of the molecular virial, by a
    /// central difference.
    pub fn virial(&self, ccom: &[Vec3], carm: &[Vec3], q: &[f64], dcom: &[Vec3], darm: &[Vec3], mu: &[Vec3]) -> f64 {
        let eps = 1e-4;
        let up = self.energy_scaled(1.0 + eps, ccom, carm, q, dcom, darm, mu);
        let dn = self.energy_scaled(1.0 - eps, ccom, carm, q, dcom, darm, mu);
        // d/d ln(lambda) = lambda d/d lambda; the points are 1 +- eps.
        -(up - dn) / ((1.0 + eps).ln() - (1.0 - eps).ln())
    }

    /// `sum_i w_i (k.a_i) e^(ik.r_i)` for every `k` (a site's weight and its arm from its
    /// molecule's centre): what the derivative of a site's phase under rigid scaling needs.
    fn arm_sum(&self, ph: &[C], w: &[f64], arm: &[Vec3]) -> Vec<C> {
        let nk = self.idx.len();
        let mut s = vec![C::default(); nk];
        let per = nk.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in s.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    let k0 = ci * per;
                    for (i, (wi, a)) in w.iter().zip(arm).enumerate() {
                        let row = &ph[i * nk + k0..i * nk + k0 + chunk.len()];
                        for (j, (sk, p)) in chunk.iter_mut().zip(row).enumerate() {
                            let k = self.kv[k0 + j];
                            let c = wi * (k[0] * a.x + k[1] * a.y + k[2] * a.z);
                            sk.re += c * p.re;
                            sk.im += c * p.im;
                        }
                    }
                });
            }
        });
        s
    }

    /// `sum_a (k.mu_a)(k.b_a) e^(ik.r_a)` for every `k`.
    fn dipole_arm_sum(&self, ph: &[C], mu: &[Vec3], arm: &[Vec3]) -> Vec<C> {
        let nk = self.idx.len();
        let mut s = vec![C::default(); nk];
        let per = nk.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in s.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    let k0 = ci * per;
                    for (a, (m, b)) in mu.iter().zip(arm).enumerate() {
                        let row = &ph[a * nk + k0..a * nk + k0 + chunk.len()];
                        for (j, (sk, p)) in chunk.iter_mut().zip(row).enumerate() {
                            let k = self.kv[k0 + j];
                            let c = (k[0] * m.x + k[1] * m.y + k[2] * m.z) * (k[0] * b.x + k[1] * b.y + k[2] * b.z);
                            sk.re += c * p.re;
                            sk.im += c * p.im;
                        }
                    }
                });
            }
        });
        s
    }

    /// The reciprocal part of the molecular virial, `-dE_rec/d ln(lambda)` under rigid scaling,
    /// analytically, from the structure factor already in hand. With `k -> k / lambda`,
    /// `V -> V lambda^3` and a site's phase `k.R + (k / lambda).a`:
    ///
    /// ```text
    /// dE/dlnL = -3E + (4 pi / V) sum f (k^2 s^2 / 2 + 2) |S|^2 + (8 pi / V) sum f Re[ S* dS/dlnL ],
    /// dS/dlnL = -i sum q (k.a) e^(ik.r) - D + sum (k.mu)(k.b) e^(ik.r),
    /// ```
    ///
    /// `D` the dipoles' own sum. Agrees with [`Ewald::virial`] (the central difference) to
    /// rounding; it costs one pass over the sites where that costs two full sets of phases.
    pub fn virial_analytic(&self, s: &[C], ph_charge: &[C], q: &[f64], carm: &[Vec3], ph_dip: &[C], mu: &[Vec3], darm: &[Vec3], d: &[C]) -> f64 {
        let tc = self.arm_sum(ph_charge, q, carm);
        let td = if mu.is_empty() { vec![C::default(); self.idx.len()] } else { self.dipole_arm_sum(ph_dip, mu, darm) };
        let mut e = 0.0;
        let mut damp = 0.0;
        let mut cross = 0.0;
        for j in 0..self.idx.len() {
            let k = self.kv[j];
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            let sj = s[j];
            let s2 = sj.re * sj.re + sj.im * sj.im;
            e += self.fk[j] * s2;
            damp += self.fk[j] * (0.5 * k2 * self.s * self.s + 2.0) * s2;
            // dS/dlnL = -i tc - d + td ; -i (a + ib) = b - i a
            let ds = C { re: tc[j].im - if d.is_empty() { 0.0 } else { d[j].re } + td[j].re, im: -tc[j].re - if d.is_empty() { 0.0 } else { d[j].im } + td[j].im };
            // Re[S* dS] = sr dsr + si dsi
            cross += self.fk[j] * (sj.re * ds.re + sj.im * ds.im);
        }
        let pref = 4.0 * PI / self.volume;
        let de = -3.0 * pref * e + pref * damp + 2.0 * pref * cross;
        -de
    }

    /// The charges' self energy in the reciprocal sum, to subtract: `sum q^2 / (sqrt(pi) s')`.
    pub fn charge_self(&self, q: &[f64]) -> f64 {
        q.iter().map(|x| x * x).sum::<f64>() / (PI.sqrt() * self.s)
    }
}

/// `erf(r / s) / r` and its derivative in `r`.
pub fn erf_over_r(r: f64, s: f64) -> (f64, f64) {
    // Two sites at the same point (a core and a cloud on one atom): the limit.
    if r < 1e-9 {
        return (2.0 / (PI.sqrt() * s), 0.0);
    }
    let x = r / s;
    let g = crate::induction::erf(x);
    let c = 2.0 / (PI.sqrt() * s) * (-x * x).exp();
    (g / r, c / r - g / (r * r))
}
