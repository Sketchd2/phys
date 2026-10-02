//! Exchange-correlation functionals: LDA and PBE, spin-resolved.
//!
//! # What is derived and what is stored
//!
//! * **Exchange** of the uniform electron gas is exact and has no parameter:
//!   `-(3/4)(6/pi)^(1/3) (rho_a^(4/3) + rho_b^(4/3))`.
//! * **Correlation** of the uniform gas is not known in closed form. It was
//!   computed by quantum Monte Carlo (Ceperley and Alder, 1980) and Perdew and
//!   Wang (1992) gave it the analytic form used here, with its high- and
//!   low-density limits built in. Its constants are a *fit to that simulation*:
//!   a derived law of a model system, stored as a shortcut, which is the input
//!   every non-empirical functional is built on. The digits are libxc's
//!   (`PW_MOD`), so that this agrees with the reference implementation.
//! * **PBE** (Perdew, Burke and Ernzerhof, 1996) adds the density's gradient.
//!   Its constants are not fitted to anything chemical: `kappa = 0.804` is the
//!   Lieb-Oxford bound, `mu` makes exchange cancel correlation's gradient term
//!   to second order, `beta` is the high-density gradient expansion of
//!   correlation and `gamma = (1 - ln 2) / pi^2` its logarithmic limit.
//!
//! # Derivatives by arithmetic, not by hand
//!
//! A Kohn-Sham potential needs the energy density's derivatives with respect
//! to both spin densities and the three gradient invariants. PBE's are long,
//! and transcribing them is where implementations go wrong. Here the energy is
//! evaluated once in forward-mode dual numbers ([`D`]) carrying all five
//! partial derivatives, so the potential is exactly the derivative of the
//! energy that is evaluated — by construction rather than by care.

use std::ops::{Add, Div, Mul, Neg, Sub};

/// A value and its derivatives with respect to
/// `(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct D {
    pub v: f64,
    pub d: [f64; 5],
}

impl D {
    pub fn c(v: f64) -> D {
        D { v, d: [0.0; 5] }
    }
    fn var(v: f64, k: usize) -> D {
        let mut d = [0.0; 5];
        d[k] = 1.0;
        D { v, d }
    }
    fn chain(self, f: f64, df: f64) -> D {
        let mut d = self.d;
        for x in &mut d {
            *x *= df;
        }
        D { v: f, d }
    }
    pub fn ln(self) -> D {
        self.chain(self.v.ln(), 1.0 / self.v)
    }
    pub fn exp(self) -> D {
        let e = self.v.exp();
        self.chain(e, e)
    }
    pub fn sqrt(self) -> D {
        let s = self.v.sqrt();
        self.chain(s, 0.5 / s)
    }
    pub fn powf(self, p: f64) -> D {
        let f = self.v.powf(p);
        self.chain(f, if self.v == 0.0 { 0.0 } else { p * f / self.v })
    }
    pub fn powi(self, p: i32) -> D {
        self.chain(self.v.powi(p), p as f64 * self.v.powi(p - 1))
    }
}

impl Add for D {
    type Output = D;
    fn add(self, o: D) -> D {
        let mut d = self.d;
        for k in 0..5 {
            d[k] += o.d[k];
        }
        D { v: self.v + o.v, d }
    }
}
impl Sub for D {
    type Output = D;
    fn sub(self, o: D) -> D {
        let mut d = self.d;
        for k in 0..5 {
            d[k] -= o.d[k];
        }
        D { v: self.v - o.v, d }
    }
}
impl Mul for D {
    type Output = D;
    fn mul(self, o: D) -> D {
        let mut d = [0.0; 5];
        for k in 0..5 {
            d[k] = self.d[k] * o.v + self.v * o.d[k];
        }
        D { v: self.v * o.v, d }
    }
}
impl Div for D {
    type Output = D;
    fn div(self, o: D) -> D {
        let inv = 1.0 / o.v;
        let mut d = [0.0; 5];
        for k in 0..5 {
            d[k] = (self.d[k] - self.v * inv * o.d[k]) * inv;
        }
        D { v: self.v * inv, d }
    }
}
impl Neg for D {
    type Output = D;
    fn neg(self) -> D {
        let mut d = self.d;
        for x in &mut d {
            *x = -*x;
        }
        D { v: -self.v, d }
    }
}
impl Add<f64> for D {
    type Output = D;
    fn add(self, o: f64) -> D {
        D { v: self.v + o, d: self.d }
    }
}
impl Sub<f64> for D {
    type Output = D;
    fn sub(self, o: f64) -> D {
        D { v: self.v - o, d: self.d }
    }
}
impl Mul<f64> for D {
    type Output = D;
    fn mul(self, o: f64) -> D {
        let mut d = self.d;
        for x in &mut d {
            *x *= o;
        }
        D { v: self.v * o, d }
    }
}
impl Div<f64> for D {
    type Output = D;
    fn div(self, o: f64) -> D {
        self * (1.0 / o)
    }
}
impl Add<D> for f64 {
    type Output = D;
    fn add(self, o: D) -> D {
        o + self
    }
}
impl Sub<D> for f64 {
    type Output = D;
    fn sub(self, o: D) -> D {
        -o + self
    }
}
impl Mul<D> for f64 {
    type Output = D;
    fn mul(self, o: D) -> D {
        o * self
    }
}
impl Div<D> for f64 {
    type Output = D;
    fn div(self, o: D) -> D {
        D::c(self) / o
    }
}

/// Which functional.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Functional {
    /// Local spin density: Slater exchange, PW92 correlation.
    Lda,
    /// Perdew-Burke-Ernzerhof.
    Pbe,
}

impl Functional {
    pub fn needs_gradient(self) -> bool {
        matches!(self, Functional::Pbe)
    }
}

/// Below this a spin density is treated as absent.
pub const DENSITY_FLOOR: f64 = 1e-14;

const PI: f64 = std::f64::consts::PI;

/// PW92's `G(rs)` for one parameter set.
fn pw_g(rs: D, a: f64, a1: f64, b1: f64, b2: f64, b3: f64, b4: f64) -> D {
    let srs = rs.sqrt();
    let den = 2.0 * a * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs);
    -2.0 * a * (1.0 + a1 * rs) * (1.0 + 1.0 / den).ln()
}

/// Uniform-gas correlation energy per electron, `(rs, zeta)`.
fn pw92(rs: D, zeta: D) -> D {
    let ec0 = pw_g(rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294);
    let ec1 = pw_g(rs, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517);
    let mac = pw_g(rs, 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671);
    let fz20 = 1.709920934161365617563962776245;
    let fz = ((1.0 + zeta).powf(4.0 / 3.0) + (1.0 - zeta).powf(4.0 / 3.0) - 2.0) / (2f64.powf(4.0 / 3.0) - 2.0);
    let z4 = zeta.powi(4);
    ec0 - mac * fz * (1.0 - z4) / fz20 + (ec1 - ec0) * fz * z4
}

/// Exchange energy density of a spin-unpolarised density `n` with
/// `|grad n|^2 = g2`, per volume: the spin-scaling relation builds the
/// polarised case from two of these.
fn exchange_unpolarised(n: D, g2: Option<D>) -> D {
    let ex_unif = -0.75 * (3.0 / PI).powf(1.0 / 3.0) * n.powf(4.0 / 3.0);
    match g2 {
        None => ex_unif,
        Some(g2) => {
            let kappa = 0.804;
            let mu = 0.2195149727645171;
            let kf = (3.0 * PI * PI * n).powf(1.0 / 3.0);
            let s2 = g2 / (4.0 * kf * kf * n * n);
            let fx = 1.0 + kappa - kappa / (1.0 + mu * s2 / kappa);
            ex_unif * fx
        }
    }
}

/// Energy density per volume and its five partial derivatives.
pub fn evaluate(f: Functional, ra: f64, rb: f64, saa: f64, sab: f64, sbb: f64) -> (f64, [f64; 5]) {
    let ra_ = ra.max(0.0);
    let rb_ = rb.max(0.0);
    if ra_ + rb_ < DENSITY_FLOOR {
        return (0.0, [0.0; 5]);
    }
    let a = D::var(ra_, 0);
    let b = D::var(rb_, 1);
    let gaa = D::var(saa.max(0.0), 2);
    let gab = D::var(sab, 3);
    let gbb = D::var(sbb.max(0.0), 4);
    let grad = f.needs_gradient();
    // Exchange, spin by spin: E_x[a, b] = (E_x[2a] + E_x[2b]) / 2.
    let mut e = D::c(0.0);
    if ra_ > DENSITY_FLOOR {
        e = e + 0.5 * exchange_unpolarised(2.0 * a, grad.then_some(4.0 * gaa));
    }
    if rb_ > DENSITY_FLOOR {
        e = e + 0.5 * exchange_unpolarised(2.0 * b, grad.then_some(4.0 * gbb));
    }
    // Correlation.
    let n = a + b;
    let zeta = ((a - b) / n).powf(1.0).clamp_unit();
    let rs = (3.0 / (4.0 * PI) / n).powf(1.0 / 3.0);
    let ec = pw92(rs, zeta);
    let mut ecorr = n * ec;
    if grad {
        let beta = 0.06672455060314922;
        let gamma = (1.0 - 2f64.ln()) / (PI * PI);
        let phi = ((1.0 + zeta).powf(2.0 / 3.0) + (1.0 - zeta).powf(2.0 / 3.0)) * 0.5;
        let kf = (3.0 * PI * PI * n).powf(1.0 / 3.0);
        let ks = (4.0 * kf / PI).sqrt();
        let g2 = gaa + 2.0 * gab + gbb;
        let t2 = g2 / (4.0 * phi * phi * ks * ks * n * n);
        let phi3 = phi.powi(3);
        let aa = (beta / gamma) / ((-ec / (gamma * phi3)).exp() - 1.0);
        let at2 = aa * t2;
        let h = gamma * phi3 * (1.0 + (beta / gamma) * t2 * (1.0 + at2) / (1.0 + at2 + at2 * at2)).ln();
        ecorr = ecorr + n * h;
    }
    e = e + ecorr;
    (e.v, e.d)
}

trait ClampUnit {
    fn clamp_unit(self) -> Self;
}
impl ClampUnit for D {
    /// Spin polarisation held inside `[-1, 1]`, where `(1 - zeta)^(4/3)` exists.
    fn clamp_unit(self) -> D {
        if self.v > 1.0 {
            D { v: 1.0, d: [0.0; 5] }
        } else if self.v < -1.0 {
            D { v: -1.0, d: [0.0; 5] }
        } else {
            self
        }
    }
}
