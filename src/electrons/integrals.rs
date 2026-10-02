//! Integrals over Gaussian shells: McMurchie and Davidson.
//!
//! The product of two Gaussians is a Gaussian about a point between them, and
//! McMurchie-Davidson expands that product in *Hermite* Gaussians, whose
//! derivatives are easy. Overlap and kinetic energy are then products of
//! one-dimensional factors; anything Coulombic — attraction to a nucleus, the
//! repulsion of two electron distributions — is a sum of Hermite Coulomb
//! integrals `R_tuv`, built by recursion from the Boys function.
//!
//! Nothing here is specific to any element or angular momentum: the same code
//! handles an `s` shell on hydrogen and a `g` shell on uranium.

use super::basis::{components, component_scale, Basis, Shell};
use super::boys::boys;
use super::linalg::Matrix;

/// Hermite expansion coefficients in one dimension:
/// `x_A^i x_B^j exp(..) = sum_t E[i][j][t] Lambda_t`.
///
/// Stored flat as `[(i * (jmax+1) + j) * (imax+jmax+1) + t]`.
struct Hermite {
    jmax: usize,
    tmax: usize,
    e: Vec<f64>,
}

impl Hermite {
    fn new(imax: usize, jmax: usize, a: f64, b: f64, ab: f64) -> Hermite {
        let p = a + b;
        let mu = a * b / p;
        let tmax = imax + jmax;
        let w = tmax + 1;
        let mut e = vec![0.0; (imax + 1) * (jmax + 1) * w];
        let idx = |i: usize, j: usize, t: usize| (i * (jmax + 1) + j) * w + t;
        // x_PA and x_PB, with P the product centre and ab = A - B.
        let xpa = -b / p * ab;
        let xpb = a / p * ab;
        let half = 0.5 / p;
        e[idx(0, 0, 0)] = (-mu * ab * ab).exp();
        for i in 0..=imax {
            for j in 0..=jmax {
                if i == 0 && j == 0 {
                    continue;
                }
                for t in 0..=(i + j) {
                    let v = if i > 0 {
                        let (ii, jj) = (i - 1, j);
                        let low = if t > 0 { e[idx(ii, jj, t - 1)] } else { 0.0 };
                        let mid = if t <= ii + jj { e[idx(ii, jj, t)] } else { 0.0 };
                        let high = if t + 1 <= ii + jj { e[idx(ii, jj, t + 1)] } else { 0.0 };
                        half * low + xpa * mid + (t + 1) as f64 * high
                    } else {
                        let (ii, jj) = (i, j - 1);
                        let low = if t > 0 { e[idx(ii, jj, t - 1)] } else { 0.0 };
                        let mid = if t <= ii + jj { e[idx(ii, jj, t)] } else { 0.0 };
                        let high = if t + 1 <= ii + jj { e[idx(ii, jj, t + 1)] } else { 0.0 };
                        half * low + xpb * mid + (t + 1) as f64 * high
                    };
                    e[idx(i, j, t)] = v;
                }
            }
        }
        Hermite { jmax, tmax, e }
    }

    #[inline]
    fn get(&self, i: usize, j: usize, t: usize) -> f64 {
        if t > i + j {
            return 0.0;
        }
        self.e[(i * (self.jmax + 1) + j) * (self.tmax + 1) + t]
    }
}

/// Hermite Coulomb integrals `R^0_tuv(p, PC)` for `t + u + v <= lmax`, flat as
/// `[(t * (lmax+1) + u) * (lmax+1) + v]`.
pub(crate) fn hermite_coulomb(lmax: usize, p: f64, pc: [f64; 3]) -> Vec<f64> {
    let mut scratch = Vec::new();
    let mut out = Vec::new();
    hermite_coulomb_into(lmax, p, pc, &mut scratch, &mut out);
    out
}

/// As [`hermite_coulomb`], into `out`, with `scratch` reused between calls:
/// this runs once per primitive quartet, and allocating two vectors each time
/// was a measurable part of building the integrals.
pub(crate) fn hermite_coulomb_into(lmax: usize, p: f64, pc: [f64; 3], scratch: &mut Vec<f64>, out: &mut Vec<f64>) {
    let n1 = lmax + 1;
    let r2 = pc[0] * pc[0] + pc[1] * pc[1] + pc[2] * pc[2];
    let mut f = [0.0f64; 32];
    boys(lmax, p * r2, &mut f[..n1]);
    let at = |n: usize, t: usize, u: usize, v: usize| ((n * n1 + t) * n1 + u) * n1 + v;
    let need = n1 * n1 * n1 * n1;
    scratch.clear();
    scratch.resize(need, 0.0);
    let r = scratch;
    let mut m2p = 1.0;
    for n in 0..=lmax {
        r[at(n, 0, 0, 0)] = m2p * f[n];
        m2p *= -2.0 * p;
    }
    for n in (0..lmax).rev() {
        let top = lmax - n;
        for t in 0..=top {
            for u in 0..=(top - t) {
                for v in 0..=(top - t - u) {
                    if t + u + v == 0 {
                        continue;
                    }
                    let val = if t > 0 {
                        let a = if t > 1 { (t - 1) as f64 * r[at(n + 1, t - 2, u, v)] } else { 0.0 };
                        a + pc[0] * r[at(n + 1, t - 1, u, v)]
                    } else if u > 0 {
                        let a = if u > 1 { (u - 1) as f64 * r[at(n + 1, t, u - 2, v)] } else { 0.0 };
                        a + pc[1] * r[at(n + 1, t, u - 1, v)]
                    } else {
                        let a = if v > 1 { (v - 1) as f64 * r[at(n + 1, t, u, v - 2)] } else { 0.0 };
                        a + pc[2] * r[at(n + 1, t, u, v - 1)]
                    };
                    r[at(n, t, u, v)] = val;
                }
            }
        }
    }
    out.clear();
    out.extend_from_slice(&r[..n1 * n1 * n1]);
}

/// A primitive pair: the Hermite tables and the product's exponent and centre.
pub(crate) struct Pair {
    p: f64,
    /// The second shell's exponent, which the kinetic integral needs.
    beta: f64,
    centre: [f64; 3],
    coef: f64,
    hx: Hermite,
    hy: Hermite,
    hz: Hermite,
}

pub(crate) fn pairs(a: &Shell, b: &Shell, extra_j: usize) -> Vec<Pair> {
    let ab = [a.centre[0] - b.centre[0], a.centre[1] - b.centre[1], a.centre[2] - b.centre[2]];
    let mut out = Vec::with_capacity(a.exponents.len() * b.exponents.len());
    for (ea, ca) in a.exponents.iter().zip(&a.coefficients) {
        for (eb, cb) in b.exponents.iter().zip(&b.coefficients) {
            let p = ea + eb;
            let centre = [
                (ea * a.centre[0] + eb * b.centre[0]) / p,
                (ea * a.centre[1] + eb * b.centre[1]) / p,
                (ea * a.centre[2] + eb * b.centre[2]) / p,
            ];
            out.push(Pair {
                p,
                beta: *eb,
                centre,
                coef: ca * cb,
                hx: Hermite::new(a.l, b.l + extra_j, *ea, *eb, ab[0]),
                hy: Hermite::new(a.l, b.l + extra_j, *ea, *eb, ab[1]),
                hz: Hermite::new(a.l, b.l + extra_j, *ea, *eb, ab[2]),
            });
        }
    }
    out
}

/// Overlap, kinetic and nuclear-attraction matrices. `nuclei` are
/// `(charge, position)` in atomic units.
pub fn one_electron(basis: &Basis, nuclei: &[(f64, [f64; 3])]) -> (Matrix, Matrix, Matrix) {
    let n = basis.size;
    let mut s = Matrix::zeros(n);
    let mut t = Matrix::zeros(n);
    let mut v = Matrix::zeros(n);
    let pi = std::f64::consts::PI;
    for (ia, a) in basis.shells.iter().enumerate() {
        for (ib, b) in basis.shells.iter().enumerate().take(ia + 1) {
            let ca = components(a.l);
            let cb = components(b.l);
            let prs = pairs(a, b, 2);
            for (ka, ka3) in ca.iter().enumerate() {
                for (kb, kb3) in cb.iter().enumerate() {
                    let scale = component_scale(a.l, *ka3) * component_scale(b.l, *kb3);
                    let (mut ss, mut tt, mut vv) = (0.0, 0.0, 0.0);
                    for pr in &prs {
                        let (i, j, k) = (ka3[0], ka3[1], ka3[2]);
                        let (l, m, nn) = (kb3[0], kb3[1], kb3[2]);
                        let ex = |d: &Hermite, i: usize, j: usize| d.get(i, j, 0);
                        let sx = ex(&pr.hx, i, l);
                        let sy = ex(&pr.hy, j, m);
                        let sz = ex(&pr.hz, k, nn);
                        let norm = (pi / pr.p).powf(1.5);
                        ss += pr.coef * sx * sy * sz * norm;
                        // Kinetic, one dimension at a time:
                        // T_d = b(2l+1) S(i,l) - 2b^2 S(i,l+2) - l(l-1)/2 S(i,l-2).
                        let kin = |h: &Hermite, i: usize, l: usize, beta: f64| {
                            let mut r = beta * (2 * l + 1) as f64 * h.get(i, l, 0) - 2.0 * beta * beta * h.get(i, l + 2, 0);
                            if l >= 2 {
                                r -= 0.5 * (l * (l - 1)) as f64 * h.get(i, l - 2, 0);
                            }
                            r
                        };
                        let beta = pr.beta;
                        let tx = kin(&pr.hx, i, l, beta);
                        let ty = kin(&pr.hy, j, m, beta);
                        let tz = kin(&pr.hz, k, nn, beta);
                        tt += pr.coef * (tx * sy * sz + sx * ty * sz + sx * sy * tz) * norm;
                        // Nuclear attraction.
                        let lsum = a.l + b.l;
                        for (zc, c) in nuclei {
                            let pc = [pr.centre[0] - c[0], pr.centre[1] - c[1], pr.centre[2] - c[2]];
                            let r = hermite_coulomb(lsum, pr.p, pc);
                            let w = lsum + 1;
                            let mut sum = 0.0;
                            for tt_ in 0..=(i + l) {
                                let etx = pr.hx.get(i, l, tt_);
                                if etx == 0.0 {
                                    continue;
                                }
                                for u in 0..=(j + m) {
                                    let ety = pr.hy.get(j, m, u);
                                    if ety == 0.0 {
                                        continue;
                                    }
                                    for vv_ in 0..=(k + nn) {
                                        sum += etx * ety * pr.hz.get(k, nn, vv_) * r[(tt_ * w + u) * w + vv_];
                                    }
                                }
                            }
                            vv -= zc * 2.0 * pi / pr.p * sum * pr.coef;
                        }
                    }
                    let (row, col) = (basis.offsets[ia] + ka, basis.offsets[ib] + kb);
                    for (m_, val) in [(&mut s, ss), (&mut t, tt), (&mut v, vv)] {
                        m_.set(row, col, val * scale);
                        m_.set(col, row, val * scale);
                    }
                }
            }
        }
    }
    (s, t, v)
}

/// The electron-repulsion block `(ab|cd)` for four shells, in chemists'
/// notation, flat over `[ka][kb][kc][kd]` components.
pub fn eri_block(a: &Shell, b: &Shell, c: &Shell, d: &Shell) -> Vec<f64> {
    eri_from_pairs(&pairs(a, b, 0), &pairs(c, d, 0), [a.l, b.l, c.l, d.l])
}

/// As [`eri_block`], from primitive pairs already built — so that one bra pair
/// against many kets builds its Hermite tables once, not once per ket.
pub(crate) fn eri_from_pairs(bra: &[Pair], ket: &[Pair], ls: [usize; 4]) -> Vec<f64> {
    let pi = std::f64::consts::PI;
    let (ca, cb, cc, cd) = (components(ls[0]), components(ls[1]), components(ls[2]), components(ls[3]));
    let ltot = ls[0] + ls[1] + ls[2] + ls[3];
    let w = ltot + 1;
    let mut out = vec![0.0; ca.len() * cb.len() * cc.len() * cd.len()];
    let mut scratch = Vec::new();
    let mut r = Vec::new();
    for p in bra {
        for q in ket {
            let alpha = p.p * q.p / (p.p + q.p);
            let pq = [p.centre[0] - q.centre[0], p.centre[1] - q.centre[1], p.centre[2] - q.centre[2]];
            hermite_coulomb_into(ltot, alpha, pq, &mut scratch, &mut r);
            let pre = 2.0 * pi.powf(2.5) / (p.p * q.p * (p.p + q.p).sqrt()) * p.coef * q.coef;
            let mut o = 0;
            for x in &ca {
                for y in &cb {
                    for z in &cc {
                        for wv in &cd {
                            let mut sum = 0.0;
                            for t in 0..=(x[0] + y[0]) {
                                let e1 = p.hx.get(x[0], y[0], t);
                                if e1 == 0.0 {
                                    continue;
                                }
                                for u in 0..=(x[1] + y[1]) {
                                    let e2 = e1 * p.hy.get(x[1], y[1], u);
                                    if e2 == 0.0 {
                                        continue;
                                    }
                                    for v in 0..=(x[2] + y[2]) {
                                        let e3 = e2 * p.hz.get(x[2], y[2], v);
                                        if e3 == 0.0 {
                                            continue;
                                        }
                                        for tau in 0..=(z[0] + wv[0]) {
                                            let f1 = q.hx.get(z[0], wv[0], tau);
                                            if f1 == 0.0 {
                                                continue;
                                            }
                                            for nu in 0..=(z[1] + wv[1]) {
                                                let f2 = f1 * q.hy.get(z[1], wv[1], nu);
                                                if f2 == 0.0 {
                                                    continue;
                                                }
                                                for phi in 0..=(z[2] + wv[2]) {
                                                    let f3 = f2 * q.hz.get(z[2], wv[2], phi);
                                                    if f3 == 0.0 {
                                                        continue;
                                                    }
                                                    let sign = if (tau + nu + phi) % 2 == 0 { 1.0 } else { -1.0 };
                                                    sum += e3 * sign * f3 * r[((t + tau) * w + (u + nu)) * w + (v + phi)];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            out[o] += pre * sum;
                            o += 1;
                        }
                    }
                }
            }
        }
    }
    let mut o = 0;
    for x in &ca {
        for y in &cb {
            for z in &cc {
                for wv in &cd {
                    out[o] *= component_scale(ls[0], *x) * component_scale(ls[1], *y) * component_scale(ls[2], *z) * component_scale(ls[3], *wv);
                    o += 1;
                }
            }
        }
    }
    out
}
