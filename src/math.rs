//! Minimal deterministic vector/tensor math.
//!
//! Everything here is `f64` with a *fixed* evaluation order. We never use
//! `mul_add` (its fusion is target-dependent) and never reduce in parallel
//! without a deterministic tree, because bit-identical replay is a hard
//! requirement of the engine (see `docs/DESIGN.md`, "Determinism contract").

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

pub const ZERO: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };

#[inline]
pub const fn v3(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3 { x, y, z }
}

impl Vec3 {
    pub const ZERO: Vec3 = ZERO;

    #[inline]
    pub const fn splat(a: f64) -> Self {
        Vec3 { x: a, y: a, z: a }
    }

    #[inline]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    #[inline]
    pub fn cross(self, o: Self) -> Self {
        v3(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    #[inline]
    pub fn norm2(self) -> f64 {
        self.dot(self)
    }

    #[inline]
    pub fn norm(self) -> f64 {
        self.norm2().sqrt()
    }

    #[inline]
    pub fn scale(self, s: f64) -> Self {
        v3(self.x * s, self.y * s, self.z * s)
    }

    /// Normalised, or `ZERO` for a zero-length vector (never NaN).
    #[inline]
    pub fn unit(self) -> Self {
        let n = self.norm();
        if n > 0.0 {
            self.scale(1.0 / n)
        } else {
            ZERO
        }
    }

    #[inline]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    #[inline]
    pub fn max_abs(self) -> f64 {
        self.x.abs().max(self.y.abs()).max(self.z.abs())
    }

    /// Outer product `self ⊗ other` as a symmetric-capable 3x3.
    pub fn outer(self, o: Self) -> Mat3 {
        Mat3::new([
            [self.x * o.x, self.x * o.y, self.x * o.z],
            [self.y * o.x, self.y * o.y, self.y * o.z],
            [self.z * o.x, self.z * o.y, self.z * o.z],
        ])
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    #[inline]
    fn add(self, o: Vec3) -> Vec3 {
        v3(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl Sub for Vec3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, o: Vec3) -> Vec3 {
        v3(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, s: f64) -> Vec3 {
        self.scale(s)
    }
}
impl Div<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn div(self, s: f64) -> Vec3 {
        self.scale(1.0 / s)
    }
}
impl Neg for Vec3 {
    type Output = Vec3;
    #[inline]
    fn neg(self) -> Vec3 {
        v3(-self.x, -self.y, -self.z)
    }
}
impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, o: Vec3) {
        *self = *self + o;
    }
}
impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, o: Vec3) {
        *self = *self - o;
    }
}

/// Row-major 3x3.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Mat3(pub [[f64; 3]; 3]);

impl Mat3 {
    pub const fn new(m: [[f64; 3]; 3]) -> Self {
        Mat3(m)
    }

    pub fn identity() -> Self {
        Mat3::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    }

    pub fn zero() -> Self {
        Mat3::new([[0.0; 3]; 3])
    }

    pub fn add(self, o: Self) -> Self {
        let mut r = Mat3::zero();
        for i in 0..3 {
            for j in 0..3 {
                r.0[i][j] = self.0[i][j] + o.0[i][j];
            }
        }
        r
    }

    pub fn scaled(self, s: f64) -> Self {
        let mut r = Mat3::zero();
        for i in 0..3 {
            for j in 0..3 {
                r.0[i][j] = self.0[i][j] * s;
            }
        }
        r
    }

    pub fn mul_vec(self, v: Vec3) -> Vec3 {
        v3(
            self.0[0][0] * v.x + self.0[0][1] * v.y + self.0[0][2] * v.z,
            self.0[1][0] * v.x + self.0[1][1] * v.y + self.0[1][2] * v.z,
            self.0[2][0] * v.x + self.0[2][1] * v.y + self.0[2][2] * v.z,
        )
    }

    pub fn trace(self) -> f64 {
        self.0[0][0] + self.0[1][1] + self.0[2][2]
    }

    /// Solve `self * x = b` by Gauss-Jordan with partial pivoting.
    ///
    /// Returns `None` when the matrix is singular to working precision, which
    /// happens for genuinely degenerate mass distributions (all mass on a line,
    /// or a single point). Callers must handle it — for the inertia tensor a
    /// singular solve means "this configuration cannot carry angular momentum
    /// about that axis", which is a physical statement, not an error.
    pub fn solve(self, b: Vec3) -> Option<Vec3> {
        let mut a = self.0;
        let mut r = [b.x, b.y, b.z];
        let scale = a
            .iter()
            .flat_map(|row| row.iter())
            .fold(0.0f64, |m, v| m.max(v.abs()));
        if scale == 0.0 || !scale.is_finite() {
            return None;
        }
        for col in 0..3 {
            // partial pivot
            let mut piv = col;
            for row in (col + 1)..3 {
                if a[row][col].abs() > a[piv][col].abs() {
                    piv = row;
                }
            }
            if a[piv][col].abs() < 1e-13 * scale {
                return None;
            }
            a.swap(col, piv);
            r.swap(col, piv);
            let d = a[col][col];
            for j in 0..3 {
                a[col][j] /= d;
            }
            r[col] /= d;
            for row in 0..3 {
                if row != col {
                    let f = a[row][col];
                    if f != 0.0 {
                        for j in 0..3 {
                            a[row][j] -= f * a[col][j];
                        }
                        r[row] -= f * r[col];
                    }
                }
            }
        }
        let out = v3(r[0], r[1], r[2]);
        if out.is_finite() {
            Some(out)
        } else {
            None
        }
    }
}

/// Deterministic pairwise (tree) summation.
///
/// Plain sequential summation of 10^7 terms loses ~7 digits; pairwise loses
/// ~log2(n). More importantly this is *associativity-fixed*, so a CPU run and a
/// GPU run that use the same tree shape produce bit-identical results.
pub fn det_sum(xs: &[f64]) -> f64 {
    match xs.len() {
        0 => 0.0,
        1 => xs[0],
        2 => xs[0] + xs[1],
        n => {
            let h = n / 2;
            det_sum(&xs[..h]) + det_sum(&xs[h..])
        }
    }
}

pub fn det_sum_v3(xs: &[Vec3]) -> Vec3 {
    match xs.len() {
        0 => ZERO,
        1 => xs[0],
        2 => xs[0] + xs[1],
        n => {
            let h = n / 2;
            det_sum_v3(&xs[..h]) + det_sum_v3(&xs[h..])
        }
    }
}

/// Pairwise sum of `f(i)` over `0..n` without materialising the slice.
pub fn det_sum_by<F: Fn(usize) -> f64>(n: usize, f: &F) -> f64 {
    fn go<F: Fn(usize) -> f64>(a: usize, b: usize, f: &F) -> f64 {
        match b - a {
            0 => 0.0,
            1 => f(a),
            2 => f(a) + f(a + 1),
            n => {
                let m = a + n / 2;
                go(a, m, f) + go(m, b, f)
            }
        }
    }
    go(0, n, f)
}

pub fn det_sum_v3_by<F: Fn(usize) -> Vec3>(n: usize, f: &F) -> Vec3 {
    fn go<F: Fn(usize) -> Vec3>(a: usize, b: usize, f: &F) -> Vec3 {
        match b - a {
            0 => ZERO,
            1 => f(a),
            2 => f(a) + f(a + 1),
            n => {
                let m = a + n / 2;
                go(a, m, f) + go(m, b, f)
            }
        }
    }
    go(0, n, f)
}

/// A rotation, as a unit quaternion.
///
/// # Why a quaternion and not an angle-axis vector
///
/// Because rotation between updates has to be *exact*, not approximate. A node
/// that is only re-solved every few minutes still has to be drawn, and asked
/// where its surface is, at every instant in between — and rigid rotation is a
/// closed-form solution, so there is no reason for the answer to drift.
/// Composing angle-axis vectors is only correct to first order in the angle;
/// composing quaternions is correct for any angle, which is what lets a planet
/// be integrated on its own slow cadence and still turn smoothly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quat {
    pub w: f64,
    pub v: Vec3,
}

impl Default for Quat {
    fn default() -> Quat {
        Quat::IDENTITY
    }
}

impl Quat {
    pub const IDENTITY: Quat = Quat { w: 1.0, v: Vec3::ZERO };

    /// A rotation of `angle` radians about `axis`.
    pub fn from_axis_angle(axis: Vec3, angle: f64) -> Quat {
        let n = axis.norm();
        if n <= 0.0 || !angle.is_finite() {
            return Quat::IDENTITY;
        }
        let half = angle * 0.5;
        Quat { w: half.cos(), v: axis.scale(half.sin() / n) }
    }

    /// The rotation an angular velocity produces in `dt` seconds. Exact for
    /// constant `omega`, which is what a rigid body has between torques.
    pub fn from_rate(omega: Vec3, dt: f64) -> Quat {
        Quat::from_axis_angle(omega, omega.norm() * dt)
    }

    pub fn conjugate(self) -> Quat {
        Quat { w: self.w, v: -self.v }
    }

    pub fn norm(self) -> f64 {
        (self.w * self.w + self.v.dot(self.v)).sqrt()
    }

    /// Renormalise. Repeated composition drifts off the unit sphere by a part
    /// in 10^16 per step, which over a long-lived node is worth spending a
    /// square root on.
    pub fn unit(self) -> Quat {
        let n = self.norm();
        if n <= 0.0 {
            return Quat::IDENTITY;
        }
        Quat { w: self.w / n, v: self.v.scale(1.0 / n) }
    }

    /// Compose: `self` applied after `other`.
    pub fn then(self, other: Quat) -> Quat {
        Quat {
            w: self.w * other.w - self.v.dot(other.v),
            v: other.v.scale(self.w) + self.v.scale(other.w) + self.v.cross(other.v),
        }
    }

    /// Rotate a vector.
    pub fn rotate(self, p: Vec3) -> Vec3 {
        // p + 2 v x (v x p + w p), the standard form: no matrix, no trig.
        let t = self.v.cross(p).scale(2.0);
        p + t.scale(self.w) + self.v.cross(t)
    }

    /// Angle of this rotation, radians.
    pub fn angle(self) -> f64 {
        2.0 * self.v.norm().atan2(self.w.abs())
    }

    pub fn is_finite(self) -> bool {
        self.w.is_finite() && self.v.is_finite()
    }
}
