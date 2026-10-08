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
    /// The first shell's exponent, which a derivative with respect to its
    /// centre needs.
    alpha: f64,
    /// The second shell's exponent, which the kinetic integral needs.
    beta: f64,
    centre: [f64; 3],
    coef: f64,
    hx: Hermite,
    hy: Hermite,
    hz: Hermite,
}

pub(crate) fn pairs(a: &Shell, b: &Shell, extra_j: usize) -> Vec<Pair> {
    pairs_ext(a, b, 0, extra_j)
}

/// Primitive pairs with Hermite tables reaching `extra_i` above the first
/// shell's angular momentum and `extra_j` above the second's.
pub(crate) fn pairs_ext(a: &Shell, b: &Shell, extra_i: usize, extra_j: usize) -> Vec<Pair> {
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
                alpha: *ea,
                beta: *eb,
                centre,
                coef: ca * cb,
                hx: Hermite::new(a.l + extra_i, b.l + extra_j, *ea, *eb, ab[0]),
                hy: Hermite::new(a.l + extra_i, b.l + extra_j, *ea, *eb, ab[1]),
                hz: Hermite::new(a.l + extra_i, b.l + extra_j, *ea, *eb, ab[2]),
            });
        }
    }
    out
}

/// Overlap, kinetic and nuclear-attraction matrices. `nuclei` are
/// `(charge, position)` in atomic units.
pub fn one_electron(basis: &Basis, nuclei: &[(f64, [f64; 3])]) -> (Matrix, Matrix, Matrix) {
    one_electron_where(basis, nuclei, &|_, _| true)
}

/// As [`one_electron`], for only the shell pairs `(a, b)`, `b <= a`, that
/// `keep` accepts; the rest of each matrix is left zero. A basis grown for a
/// molecule estimates each candidate against the current functions and itself
/// and never against another candidate, and the blocks between candidates
/// were most of the work.
///
/// Shell pairs are spread across threads, and the nuclear attraction's Hermite
/// Coulomb table is built once per primitive pair and nucleus rather than once
/// per component pair as well: it does not depend on the components, and for
/// a d against an f it was being built sixty times over. Each element is the
/// same sum in the same order as before, so the matrices are bit-identical.
pub fn one_electron_where(basis: &Basis, nuclei: &[(f64, [f64; 3])], keep: &(dyn Fn(usize, usize) -> bool + Sync)) -> (Matrix, Matrix, Matrix) {
    let n = basis.size;
    let pi = std::f64::consts::PI;
    let shell_pairs: Vec<(usize, usize)> = (0..basis.shells.len()).flat_map(|a| (0..=a).map(move |b| (a, b))).filter(|&(a, b)| keep(a, b)).collect();
    let job = |idx: &mut dyn Iterator<Item = usize>| {
        let mut out: Vec<(usize, usize, Vec<[f64; 3]>)> = Vec::new();
        for (ia, ib) in idx.map(|i| shell_pairs[i]) {
            let (a, b) = (&basis.shells[ia], &basis.shells[ib]);
            let ca = components(a.l);
            let cb = components(b.l);
            let prs = pairs(a, b, 2);
            let lsum = a.l + b.l;
            let w = lsum + 1;
            // Per primitive pair, per nucleus: the Hermite Coulomb table.
            let rs: Vec<Vec<Vec<f64>>> = prs.iter().map(|pr| {
                nuclei.iter().map(|(_, c)| {
                    let pc = [pr.centre[0] - c[0], pr.centre[1] - c[1], pr.centre[2] - c[2]];
                    hermite_coulomb(lsum, pr.p, pc)
                }).collect()
            }).collect();
            let mut block = vec![[0.0; 3]; ca.len() * cb.len()];
            for (ka, ka3) in ca.iter().enumerate() {
                for (kb, kb3) in cb.iter().enumerate() {
                    let scale = component_scale(a.l, *ka3) * component_scale(b.l, *kb3);
                    let (mut ss, mut tt, mut vv) = (0.0, 0.0, 0.0);
                    for (pi_, pr) in prs.iter().enumerate() {
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
                        for (ci, (zc, _)) in nuclei.iter().enumerate() {
                            let r = &rs[pi_][ci];
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
                    block[ka * cb.len() + kb] = [ss * scale, tt * scale, vv * scale];
                }
            }
            out.push((ia, ib, block));
        }
        out
    };
    let mut s = Matrix::zeros(n);
    let mut t = Matrix::zeros(n);
    let mut v = Matrix::zeros(n);
    for part in super::scf::parallel_interleaved(shell_pairs.len(), &job) {
        for (ia, ib, block) in part {
            let nb = basis.shells[ib].size();
            for ka in 0..basis.shells[ia].size() {
                for kb in 0..nb {
                    let (row, col) = (basis.offsets[ia] + ka, basis.offsets[ib] + kb);
                    let [x, y, z] = block[ka * nb + kb];
                    for (m_, val) in [(&mut s, x), (&mut t, y), (&mut v, z)] {
                        m_.set(row, col, val);
                        m_.set(col, row, val);
                    }
                }
            }
        }
    }
    (s, t, v)
}

/// The dipole matrices about `origin`: `<a| (x - C_x) |b>`, then y, then z, over
/// the basis. With overlap, kinetic and nuclear attraction these are the one-
/// electron integrals of a uniform electric field, which is what lets the
/// engine ask a molecule how it answers one (its polarisability, by finite
/// field). In McMurchie-Davidson form the one-dimensional moment is
/// `(E_1 + (P - C) E_0) sqrt(pi / p)`: the product's Hermite expansion has
/// its first term at the product's centre `P` and its second integrates to one.
pub fn dipole(basis: &Basis, origin: [f64; 3]) -> [Matrix; 3] {
    let n = basis.size;
    let pi = std::f64::consts::PI;
    let mut out = [Matrix::zeros(n), Matrix::zeros(n), Matrix::zeros(n)];
    for ia in 0..basis.shells.len() {
        for ib in 0..=ia {
            let (a, b) = (&basis.shells[ia], &basis.shells[ib]);
            let (ca, cb) = (components(a.l), components(b.l));
            let prs = pairs(a, b, 0);
            for (ka, ka3) in ca.iter().enumerate() {
                for (kb, kb3) in cb.iter().enumerate() {
                    let scale = component_scale(a.l, *ka3) * component_scale(b.l, *kb3);
                    let mut m = [0.0f64; 3];
                    for pr in &prs {
                        let (i, j, k) = (ka3[0], ka3[1], ka3[2]);
                        let (l, mm, nn) = (kb3[0], kb3[1], kb3[2]);
                        let s0 = [pr.hx.get(i, l, 0), pr.hy.get(j, mm, 0), pr.hz.get(k, nn, 0)];
                        let s1 = [pr.hx.get(i, l, 1), pr.hy.get(j, mm, 1), pr.hz.get(k, nn, 1)];
                        let norm = (pi / pr.p).powf(1.5);
                        for d in 0..3 {
                            let moment = s1[d] + (pr.centre[d] - origin[d]) * s0[d];
                            let others: f64 = (0..3).filter(|&e| e != d).map(|e| s0[e]).product();
                            m[d] += pr.coef * moment * others * norm;
                        }
                    }
                    let (row, col) = (basis.offsets[ia] + ka, basis.offsets[ib] + kb);
                    for d in 0..3 {
                        out[d].set(row, col, m[d] * scale);
                        out[d].set(col, row, m[d] * scale);
                    }
                }
            }
        }
    }
    out
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


/// The Hermite sum of one electron-repulsion integral over *unnormalised*
/// Cartesian components `x y | z w`, for one primitive quartet, given the
/// Hermite Coulomb table `r` of width `w`.
#[allow(clippy::too_many_arguments)]
#[inline]
fn eri_raw(p: &Pair, q: &Pair, x: [usize; 3], y: [usize; 3], z: [usize; 3], wv: [usize; 3], r: &[f64], w: usize) -> f64 {
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
    sum
}

/// `d/dA (ab|cd)` along x, y and z, `A` the centre of `a`: three blocks laid
/// out as [`eri_block`]'s.
///
/// Moving a Gaussian's centre is shifting its power:
/// `d/dA_x [x^i exp(-a x^2)] = 2a x^(i+1) exp(..) - i x^(i-1) exp(..)`
/// (with `x` measured from `A`), so the derivative is two integrals one
/// angular momentum either side, with the same normalisation as the original.
pub fn eri_derivative(a: &Shell, b: &Shell, c: &Shell, d: &Shell) -> [Vec<f64>; 3] {
    let pi = std::f64::consts::PI;
    let (ca, cb, cc, cd) = (components(a.l), components(b.l), components(c.l), components(d.l));
    let bra = pairs_ext(a, b, 1, 0);
    let ket = pairs(c, d, 0);
    let ltot = a.l + 1 + b.l + c.l + d.l;
    let w = ltot + 1;
    let size = ca.len() * cb.len() * cc.len() * cd.len();
    let mut out = [vec![0.0; size], vec![0.0; size], vec![0.0; size]];
    let mut scratch = Vec::new();
    let mut r = Vec::new();
    for p in &bra {
        for q in &ket {
            let alpha = p.p * q.p / (p.p + q.p);
            let pq = [p.centre[0] - q.centre[0], p.centre[1] - q.centre[1], p.centre[2] - q.centre[2]];
            hermite_coulomb_into(ltot, alpha, pq, &mut scratch, &mut r);
            let pre = 2.0 * pi.powf(2.5) / (p.p * q.p * (p.p + q.p).sqrt()) * p.coef * q.coef;
            let mut o = 0;
            for x in &ca {
                for y in &cb {
                    for z in &cc {
                        for wv in &cd {
                            for dir in 0..3 {
                                let mut up = *x;
                                up[dir] += 1;
                                let mut val = 2.0 * p.alpha * eri_raw(p, q, up, *y, *z, *wv, &r, w);
                                if x[dir] > 0 {
                                    let mut down = *x;
                                    down[dir] -= 1;
                                    val -= x[dir] as f64 * eri_raw(p, q, down, *y, *z, *wv, &r, w);
                                }
                                out[dir][o] += pre * val;
                            }
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
                    let sc = component_scale(a.l, *x) * component_scale(b.l, *y) * component_scale(c.l, *z) * component_scale(d.l, *wv);
                    for dir in 0..3 {
                        out[dir][o] *= sc;
                    }
                    o += 1;
                }
            }
        }
    }
    out
}

/// Hermite potential of one auxiliary component at every bra Hermite index:
/// `g[t,u,v] = sum_(tau nu phi) (-1)^(tau+nu+phi) f_(tau nu phi) R_(t+tau, u+nu, v+phi)`
/// for `t + u + v <= lab`, laid out `(t * (lab+1) + u) * (lab+1) + v`. `f` is
/// the component's own Hermite expansion, `z` its powers, and `scale` is
/// folded in as a weight so that several components can be summed into one.
#[inline]
fn add_hermite_potential(q: &Pair, z: [usize; 3], scale: f64, r: &[f64], w: usize, lab: usize, g: &mut [f64]) {
    let h = lab + 1;
    for tau in 0..=z[0] {
        let f1 = q.hx.get(z[0], 0, tau);
        if f1 == 0.0 {
            continue;
        }
        for nu in 0..=z[1] {
            let f2 = f1 * q.hy.get(z[1], 0, nu);
            if f2 == 0.0 {
                continue;
            }
            for phi in 0..=z[2] {
                let f3 = f2 * q.hz.get(z[2], 0, phi);
                if f3 == 0.0 {
                    continue;
                }
                let sign = if (tau + nu + phi) % 2 == 0 { scale } else { -scale };
                let c = sign * f3;
                for t in 0..=lab {
                    for u in 0..=(lab - t) {
                        let row = ((t + tau) * w + (u + nu)) * w + phi;
                        let base = (t * h + u) * h;
                        for v in 0..=(lab - t - u) {
                            g[base + v] += c * r[row + v];
                        }
                    }
                }
            }
        }
    }
}

/// `sum_tuv E^x_t E^y_u E^z_v g[t,u,v]` for bra powers `x` on the first
/// centre and `y` on the second: the bra's half of a Hermite contraction.
#[inline]
fn bra_contract(p: &Pair, x: [usize; 3], y: [usize; 3], g: &[f64], lab: usize) -> f64 {
    let h = lab + 1;
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
            let base = (t * h + u) * h;
            let mut s3 = 0.0;
            for v in 0..=(x[2] + y[2]) {
                s3 += p.hz.get(x[2], y[2], v) * g[base + v];
            }
            sum += e2 * s3;
        }
    }
    sum
}

/// Three-centre integrals `(ab|P)` of a bra pair against an auxiliary shell
/// (paired with the unit function), laid out as [`eri_block`]'s
/// `[ka][kb][kp]` — the same numbers as `eri_from_pairs` with
/// `ls = [la, lb, lp, 0]`.
///
/// Organised so the auxiliary side is summed once per primitive quartet: each
/// auxiliary component becomes a Hermite *potential* over the bra's Hermite
/// indices, which every bra component then reads. Done the other way round —
/// the full six-deep sum per pair of components — the same work was repeated
/// for every bra component, and was most of the time it took to set up a
/// fitted Coulomb matrix and nearly all of its gradient.
pub(crate) fn eri_three(bra: &[Pair], ket: &[Pair], la: usize, lb: usize, lp: usize) -> Vec<f64> {
    let mut scratch = ThreeScratch::new();
    let mut out = Vec::new();
    eri_three_into(bra, ket, la, lb, lp, &mut scratch, &mut out);
    out
}

/// Working space for [`eri_three_into`], kept by one thread across calls.
///
/// Every call used to allocate seven small vectors, and a solve makes around
/// ten million calls. Spread across threads that is slower than one thread on
/// Windows, where the threads queue on the allocator: measured, an
/// allocation-heavy job ran 0.5x on twelve threads against 10.0x for the same
/// job's arithmetic, and the three-centre integrals of water's round-9 basis
/// ran on about 1.4 cores.
pub(crate) struct ThreeScratch {
    scratch: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
    /// `components(l)` for every `l` up to the highest asked for so far.
    comps: Vec<Vec<[usize; 3]>>,
    /// Per-component scale factors, auxiliary and bra.
    ps: Vec<f64>,
    ab: Vec<f64>,
}

impl ThreeScratch {
    pub(crate) fn new() -> ThreeScratch {
        ThreeScratch { scratch: Vec::new(), r: Vec::new(), g: Vec::new(), comps: Vec::new(), ps: Vec::new(), ab: Vec::new() }
    }

    fn reach(&mut self, l: usize) {
        while self.comps.len() <= l {
            let next = self.comps.len();
            self.comps.push(components(next));
        }
    }
}

/// [`eri_three`] into `out`, using `s` for its working space. The same sums in
/// the same order, so the values are bit-identical.
pub(crate) fn eri_three_into(bra: &[Pair], ket: &[Pair], la: usize, lb: usize, lp: usize, s: &mut ThreeScratch, out: &mut Vec<f64>) {
    let pi = std::f64::consts::PI;
    s.reach(la.max(lb).max(lp));
    let ThreeScratch { scratch, r, g, comps, .. } = s;
    let (ca, cb, cp) = (&comps[la], &comps[lb], &comps[lp]);
    let lab = la + lb;
    let ltot = lab + lp;
    let w = ltot + 1;
    let h = lab + 1;
    let (na, nb, np) = (ca.len(), cb.len(), cp.len());
    out.clear();
    out.resize(na * nb * np, 0.0);
    g.clear();
    g.resize(h * h * h, 0.0);
    for p in bra {
        for q in ket {
            let alpha = p.p * q.p / (p.p + q.p);
            let pq = [p.centre[0] - q.centre[0], p.centre[1] - q.centre[1], p.centre[2] - q.centre[2]];
            hermite_coulomb_into(ltot, alpha, pq, scratch, r);
            let pre = 2.0 * pi.powf(2.5) / (p.p * q.p * (p.p + q.p).sqrt()) * p.coef * q.coef;
            for (kk, z) in cp.iter().enumerate() {
                g.iter_mut().for_each(|x| *x = 0.0);
                add_hermite_potential(q, *z, 1.0, r, w, lab, g);
                for (ia, x) in ca.iter().enumerate() {
                    for (ib, y) in cb.iter().enumerate() {
                        out[(ia * nb + ib) * np + kk] += pre * bra_contract(p, *x, *y, g, lab);
                    }
                }
            }
        }
    }
    for (ia, x) in ca.iter().enumerate() {
        for (ib, y) in cb.iter().enumerate() {
            for (kk, z) in cp.iter().enumerate() {
                out[(ia * nb + ib) * np + kk] *= component_scale(la, *x) * component_scale(lb, *y) * component_scale(lp, *z);
            }
        }
    }
}

/// The derivatives, with respect to the two bra centres, of
/// `sum_(mu nu k) D_(mu nu) c_k (mu nu|k)` for one bra shell pair and one
/// auxiliary shell: `d` is the density over the pair's components
/// (`[ka][kb]`) and `c` the fitted coefficients over the auxiliary shell's.
///
/// `bra` must be built with Hermite tables one higher on both centres
/// (`pairs_ext(a, b, 1, 1)`). The fitted coefficients are summed into one
/// Hermite potential before anything touches the bra, so the auxiliary side
/// costs one pass per primitive quartet whatever its angular momentum.
pub(crate) fn eri_three_gradient(bra: &[Pair], ket: &[Pair], la: usize, lb: usize, lp: usize, d: &[f64], c: &[f64]) -> ([f64; 3], [f64; 3]) {
    eri_three_gradient_with(bra, ket, la, lb, lp, d, c, &mut ThreeScratch::new())
}

/// [`eri_three_gradient`] using `s` for its working space (see
/// [`ThreeScratch`]). The same sums in the same order.
#[allow(clippy::too_many_arguments)]
pub(crate) fn eri_three_gradient_with(bra: &[Pair], ket: &[Pair], la: usize, lb: usize, lp: usize, d: &[f64], c: &[f64], s: &mut ThreeScratch) -> ([f64; 3], [f64; 3]) {
    let pi = std::f64::consts::PI;
    s.reach(la.max(lb).max(lp));
    let ThreeScratch { scratch, r, g, comps, ps: pscale, ab: abscale } = s;
    let (ca, cb, cp) = (&comps[la], &comps[lb], &comps[lp]);
    let lab = la + lb + 1;
    let ltot = lab + lp;
    let w = ltot + 1;
    let h = lab + 1;
    let nb = cb.len();
    let mut ga = [0.0; 3];
    let mut gb = [0.0; 3];
    g.clear();
    g.resize(h * h * h, 0.0);
    pscale.clear();
    pscale.extend(cp.iter().map(|z| component_scale(lp, *z)));
    abscale.clear();
    abscale.extend(ca.iter().flat_map(|x| cb.iter().map(move |y| component_scale(la, *x) * component_scale(lb, *y))));
    for p in bra {
        for q in ket {
            let alpha = p.p * q.p / (p.p + q.p);
            let pq = [p.centre[0] - q.centre[0], p.centre[1] - q.centre[1], p.centre[2] - q.centre[2]];
            hermite_coulomb_into(ltot, alpha, pq, scratch, r);
            let pre = 2.0 * pi.powf(2.5) / (p.p * q.p * (p.p + q.p).sqrt()) * p.coef * q.coef;
            g.iter_mut().for_each(|x| *x = 0.0);
            for (kk, z) in cp.iter().enumerate() {
                if c[kk] != 0.0 {
                    add_hermite_potential(q, *z, c[kk] * pscale[kk], r, w, lab, g);
                }
            }
            for (ia, x) in ca.iter().enumerate() {
                for (ib, y) in cb.iter().enumerate() {
                    let dm = d[ia * nb + ib];
                    if dm == 0.0 {
                        continue;
                    }
                    let wt = pre * dm * abscale[ia * nb + ib];
                    for dir in 0..3 {
                        let mut up = *x;
                        up[dir] += 1;
                        let mut va = 2.0 * p.alpha * bra_contract(p, up, *y, g, lab);
                        if x[dir] > 0 {
                            let mut down = *x;
                            down[dir] -= 1;
                            va -= x[dir] as f64 * bra_contract(p, down, *y, g, lab);
                        }
                        let mut up = *y;
                        up[dir] += 1;
                        let mut vb = 2.0 * p.beta * bra_contract(p, *x, up, g, lab);
                        if y[dir] > 0 {
                            let mut down = *y;
                            down[dir] -= 1;
                            vb -= y[dir] as f64 * bra_contract(p, *x, down, g, lab);
                        }
                        ga[dir] += wt * va;
                        gb[dir] += wt * vb;
                    }
                }
            }
        }
    }
    (ga, gb)
}

/// `sum_k c_k (ab|k)` over one auxiliary shell's components: the potential of
/// a fitted density between two functions, laid out `[ka][kb]`. The fitted
/// coefficients are summed into one Hermite potential first, as in
/// [`eri_three_gradient`], so the auxiliary side costs one pass per primitive
/// quartet.
pub(crate) fn eri_three_contracted(bra: &[Pair], ket: &[Pair], la: usize, lb: usize, lp: usize, c: &[f64]) -> Vec<f64> {
    let mut scratch = ThreeScratch::new();
    let mut out = Vec::new();
    eri_three_contracted_into(bra, ket, la, lb, lp, c, &mut scratch, &mut out);
    out
}

/// [`eri_three_contracted`] into `out`, using `s` for its working space (see
/// [`ThreeScratch`]). The same sums in the same order.
#[allow(clippy::too_many_arguments)]
pub(crate) fn eri_three_contracted_into(bra: &[Pair], ket: &[Pair], la: usize, lb: usize, lp: usize, c: &[f64], s: &mut ThreeScratch, out: &mut Vec<f64>) {
    let pi = std::f64::consts::PI;
    s.reach(la.max(lb).max(lp));
    let ThreeScratch { scratch, r, g, comps, ps: pscale, .. } = s;
    let (ca, cb, cp) = (&comps[la], &comps[lb], &comps[lp]);
    let lab = la + lb;
    let ltot = lab + lp;
    let w = ltot + 1;
    let h = lab + 1;
    let nb = cb.len();
    out.clear();
    out.resize(ca.len() * nb, 0.0);
    g.clear();
    g.resize(h * h * h, 0.0);
    pscale.clear();
    pscale.extend(cp.iter().map(|z| component_scale(lp, *z)));
    for p in bra {
        for q in ket {
            let alpha = p.p * q.p / (p.p + q.p);
            let pq = [p.centre[0] - q.centre[0], p.centre[1] - q.centre[1], p.centre[2] - q.centre[2]];
            hermite_coulomb_into(ltot, alpha, pq, scratch, r);
            let pre = 2.0 * pi.powf(2.5) / (p.p * q.p * (p.p + q.p).sqrt()) * p.coef * q.coef;
            g.iter_mut().for_each(|x| *x = 0.0);
            for (kk, z) in cp.iter().enumerate() {
                if c[kk] != 0.0 {
                    add_hermite_potential(q, *z, c[kk] * pscale[kk], r, w, lab, g);
                }
            }
            for (ia, x) in ca.iter().enumerate() {
                for (ib, y) in cb.iter().enumerate() {
                    out[ia * nb + ib] += pre * bra_contract(p, *x, *y, g, lab);
                }
            }
        }
    }
    for (ia, x) in ca.iter().enumerate() {
        for (ib, y) in cb.iter().enumerate() {
            out[ia * nb + ib] *= component_scale(la, *x) * component_scale(lb, *y);
        }
    }
}

/// [`eri_three_contracted`] for three shells.
pub fn eri_three_contracted_block(a: &Shell, b: &Shell, aux: &Shell, c: &[f64]) -> Vec<f64> {
    eri_three_contracted(&pairs(a, b, 0), &pairs(aux, &Shell::unit(), 0), a.l, b.l, aux.l, c)
}

/// [`eri_three`] for three shells: `(ab|P)`.
pub fn eri_three_block(a: &Shell, b: &Shell, aux: &Shell) -> Vec<f64> {
    eri_three(&pairs(a, b, 0), &pairs(aux, &Shell::unit(), 0), a.l, b.l, aux.l)
}

/// [`eri_three_gradient`] for three shells.
pub fn eri_three_gradient_block(a: &Shell, b: &Shell, aux: &Shell, d: &[f64], c: &[f64]) -> ([f64; 3], [f64; 3]) {
    eri_three_gradient(&pairs_ext(a, b, 1, 1), &pairs(aux, &Shell::unit(), 0), a.l, b.l, aux.l, d, c)
}

/// Gradient contributions of the one-electron energy, given the total density
/// `d` and energy-weighted density `w`: for each nucleus,
/// `sum D dH/dR - sum W dS/dR`, where `H` is kinetic plus nuclear attraction
/// and the derivative includes both the moving basis functions and the
/// nucleus's own pull (Hellmann-Feynman). `owner[s]` is the nucleus shell `s`
/// sits on.
pub fn one_electron_gradient(basis: &Basis, nuclei: &[(f64, [f64; 3])], owner: &[usize], d: &Matrix, w: &Matrix) -> Vec<[f64; 3]> {
    let pi = std::f64::consts::PI;
    let n = basis.size;
    let mut g = vec![[0.0; 3]; nuclei.len()];
    for (ia, a) in basis.shells.iter().enumerate() {
        for (ib, b) in basis.shells.iter().enumerate() {
            let ca = components(a.l);
            let cb = components(b.l);
            let prs = pairs_ext(a, b, 1, 2);
            let lsum = a.l + 1 + b.l;
            for (ka, x) in ca.iter().enumerate() {
                for (kb, y) in cb.iter().enumerate() {
                    let (mu, nu) = (basis.offsets[ia] + ka, basis.offsets[ib] + kb);
                    let (dmn, wmn) = (d.a[mu * n + nu], w.a[mu * n + nu]);
                    let sc = component_scale(a.l, *x) * component_scale(b.l, *y);
                    for pr in &prs {
                        let norm = (pi / pr.p).powf(1.5);
                        let s1 = |h: &Hermite, i: usize, j: usize| h.get(i, j, 0);
                        let kin = |h: &Hermite, i: usize, l: usize| {
                            let beta = pr.beta;
                            let mut r = beta * (2 * l + 1) as f64 * h.get(i, l, 0) - 2.0 * beta * beta * h.get(i, l + 2, 0);
                            if l >= 2 {
                                r -= 0.5 * (l * (l - 1)) as f64 * h.get(i, l - 2, 0);
                            }
                            r
                        };
                        // Raw overlap and kinetic for a bra of powers `c`.
                        let st = |c: [usize; 3]| {
                            let (sx, sy, sz) = (s1(&pr.hx, c[0], y[0]), s1(&pr.hy, c[1], y[1]), s1(&pr.hz, c[2], y[2]));
                            let (tx, ty, tz) = (kin(&pr.hx, c[0], y[0]), kin(&pr.hy, c[1], y[1]), kin(&pr.hz, c[2], y[2]));
                            (sx * sy * sz * norm, (tx * sy * sz + sx * ty * sz + sx * sy * tz) * norm)
                        };
                        // Raw nuclear attraction to every nucleus.
                        let vn = |c: [usize; 3], r: &[f64], wd: usize| {
                            let mut sum = 0.0;
                            for t in 0..=(c[0] + y[0]) {
                                let ex = pr.hx.get(c[0], y[0], t);
                                if ex == 0.0 {
                                    continue;
                                }
                                for u in 0..=(c[1] + y[1]) {
                                    let ey = pr.hy.get(c[1], y[1], u);
                                    if ey == 0.0 {
                                        continue;
                                    }
                                    for v in 0..=(c[2] + y[2]) {
                                        sum += ex * ey * pr.hz.get(c[2], y[2], v) * r[(t * wd + u) * wd + v];
                                    }
                                }
                            }
                            sum
                        };
                        let tables: Vec<Vec<f64>> = nuclei.iter().map(|(_, cpos)| {
                            hermite_coulomb(lsum + 1, pr.p, [pr.centre[0] - cpos[0], pr.centre[1] - cpos[1], pr.centre[2] - cpos[2]])
                        }).collect();
                        let wd = lsum + 2;
                        for dir in 0..3 {
                            // Basis-function derivative on the bra (a's nucleus).
                            let mut up = *x;
                            up[dir] += 1;
                            let (mut ds, mut dt) = st(up);
                            ds *= 2.0 * pr.alpha;
                            dt *= 2.0 * pr.alpha;
                            let mut dv = 0.0;
                            for (k, (z, _)) in nuclei.iter().enumerate() {
                                dv -= z * 2.0 * pi / pr.p * 2.0 * pr.alpha * vn(up, &tables[k], wd);
                            }
                            if x[dir] > 0 {
                                let mut down = *x;
                                down[dir] -= 1;
                                let (s0, t0) = st(down);
                                ds -= x[dir] as f64 * s0;
                                dt -= x[dir] as f64 * t0;
                                for (k, (z, _)) in nuclei.iter().enumerate() {
                                    dv += z * 2.0 * pi / pr.p * x[dir] as f64 * vn(down, &tables[k], wd);
                                }
                            }
                            // Both the bra and the ket move with their atoms; the
                            // full double loop over shells counts each once as the
                            // bra, so the factor is 2 for a symmetric density.
                            g[owner[ia]][dir] += 2.0 * pr.coef * sc * (dmn * (dt + dv) - wmn * ds);
                        }
                        // Hellmann-Feynman: the nucleus's own pull on the pair,
                        // dV/dC = +Z (2 pi / p) sum E E E R_(t+1).
                        for (k, (z, _)) in nuclei.iter().enumerate() {
                            for dir in 0..3 {
                                let mut sum = 0.0;
                                for t in 0..=(x[0] + y[0]) {
                                    let ex = pr.hx.get(x[0], y[0], t);
                                    if ex == 0.0 {
                                        continue;
                                    }
                                    for u in 0..=(x[1] + y[1]) {
                                        let ey = pr.hy.get(x[1], y[1], u);
                                        if ey == 0.0 {
                                            continue;
                                        }
                                        for v in 0..=(x[2] + y[2]) {
                                            let ez = pr.hz.get(x[2], y[2], v);
                                            let (tt, uu, vv) = match dir {
                                                0 => (t + 1, u, v),
                                                1 => (t, u + 1, v),
                                                _ => (t, u, v + 1),
                                            };
                                            sum += ex * ey * ez * tables[k][(tt * wd + uu) * wd + vv];
                                        }
                                    }
                                }
                                g[k][dir] += dmn * pr.coef * sc * z * 2.0 * pi / pr.p * sum;
                            }
                        }
                    }
                }
            }
        }
    }
    g
}
