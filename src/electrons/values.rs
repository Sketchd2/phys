//! Basis functions evaluated at points, with their gradients — what a grid
//! integral and a gradient-corrected functional need.

use super::basis::{component_scale, components, Basis};

/// Values of the listed shells' functions at `p`, packed in list order into
/// `out` (and `grad`, three per function) — what a screened batch needs,
/// without evaluating shells that are zero there.
pub fn at_shells(basis: &Basis, shells: &[usize], p: [f64; 3], out: &mut [f64], mut grad: Option<&mut [f64]>) {
    let mut o = 0;
    for &is in shells {
        let sh = &basis.shells[is];
        let d = [p[0] - sh.centre[0], p[1] - sh.centre[1], p[2] - sh.centre[2]];
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        let (mut rad, mut drad) = (0.0, 0.0);
        for (a, c) in sh.exponents.iter().zip(&sh.coefficients) {
            let e = c * (-a * r2).exp();
            rad += e;
            drad -= 2.0 * a * e;
        }
        for c in components(sh.l).iter() {
            let scale = component_scale(sh.l, *c);
            let pw = |x: f64, n: usize| if n == 0 { 1.0 } else { x.powi(n as i32) };
            let (px, py, pz) = (pw(d[0], c[0]), pw(d[1], c[1]), pw(d[2], c[2]));
            let ang = px * py * pz;
            out[o] = scale * ang * rad;
            if let Some(g) = grad.as_deref_mut() {
                let dp = |x: f64, n: usize| if n == 0 { 0.0 } else { n as f64 * x.powi(n as i32 - 1) };
                g[3 * o] = scale * (dp(d[0], c[0]) * py * pz * rad + ang * drad * d[0]);
                g[3 * o + 1] = scale * (px * dp(d[1], c[1]) * pz * rad + ang * drad * d[1]);
                g[3 * o + 2] = scale * (px * py * dp(d[2], c[2]) * rad + ang * drad * d[2]);
            }
            o += 1;
        }
    }
}

/// Values of every basis function at `p`, into `out[0..basis.size]`, and the
/// gradient into `grad` (`3 * basis.size`, `[x0, y0, z0, x1, ...]`) if given.
pub fn at(basis: &Basis, p: [f64; 3], out: &mut [f64], mut grad: Option<&mut [f64]>) {
    for (is, sh) in basis.shells.iter().enumerate() {
        let d = [p[0] - sh.centre[0], p[1] - sh.centre[1], p[2] - sh.centre[2]];
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        // radial part and its derivative factor: sum c exp(-a r2), sum -2 a c exp(-a r2)
        let (mut rad, mut drad) = (0.0, 0.0);
        for (a, c) in sh.exponents.iter().zip(&sh.coefficients) {
            let e = c * (-a * r2).exp();
            rad += e;
            drad -= 2.0 * a * e;
        }
        let off = basis.offsets[is];
        for (k, c) in components(sh.l).iter().enumerate() {
            let scale = component_scale(sh.l, *c);
            let pw = |x: f64, n: usize| if n == 0 { 1.0 } else { x.powi(n as i32) };
            let (px, py, pz) = (pw(d[0], c[0]), pw(d[1], c[1]), pw(d[2], c[2]));
            let ang = px * py * pz;
            out[off + k] = scale * ang * rad;
            if let Some(g) = grad.as_deref_mut() {
                let dp = |x: f64, n: usize| if n == 0 { 0.0 } else { n as f64 * x.powi(n as i32 - 1) };
                let gx = dp(d[0], c[0]) * py * pz * rad + ang * drad * d[0];
                let gy = px * dp(d[1], c[1]) * pz * rad + ang * drad * d[1];
                let gz = px * py * dp(d[2], c[2]) * rad + ang * drad * d[2];
                g[3 * (off + k)] = scale * gx;
                g[3 * (off + k) + 1] = scale * gy;
                g[3 * (off + k) + 2] = scale * gz;
            }
        }
    }
}

/// Values, gradients and second derivatives of the listed shells' functions at
/// `p`, packed in list order: `hess` holds six per function
/// (`xx, yy, zz, xy, xz, yz`). What a gradient-corrected functional's force
/// needs, since its potential already involves first derivatives.
pub fn at_shells_hessian(basis: &Basis, shells: &[usize], p: [f64; 3], val: &mut [f64], grad: &mut [f64], hess: &mut [f64]) {
    let mut o = 0;
    for &is in shells {
        let sh = &basis.shells[is];
        let d = [p[0] - sh.centre[0], p[1] - sh.centre[1], p[2] - sh.centre[2]];
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        for c in components(sh.l).iter() {
            let scale = component_scale(sh.l, *c);
            let (mut v, mut g, mut h) = (0.0, [0.0; 3], [0.0; 6]);
            for (a, coef) in sh.exponents.iter().zip(&sh.coefficients) {
                let e = coef * (-a * r2).exp();
                // Per axis: f, f', f'' of x^n exp(-a x^2) without the common factor.
                let axis = |x: f64, n: usize| {
                    let pw = |k: i32| if k < 0 { 0.0 } else { x.powi(k) };
                    let n_ = n as i32;
                    let f0 = pw(n_);
                    let f1 = n as f64 * pw(n_ - 1) - 2.0 * a * pw(n_ + 1);
                    let f2 = (n * n.saturating_sub(1)) as f64 * pw(n_ - 2) - 2.0 * a * (2 * n + 1) as f64 * pw(n_) + 4.0 * a * a * pw(n_ + 2);
                    (f0, f1, f2)
                };
                let (x0, x1, x2) = axis(d[0], c[0]);
                let (y0, y1, y2) = axis(d[1], c[1]);
                let (z0, z1, z2) = axis(d[2], c[2]);
                v += e * x0 * y0 * z0;
                g[0] += e * x1 * y0 * z0;
                g[1] += e * x0 * y1 * z0;
                g[2] += e * x0 * y0 * z1;
                h[0] += e * x2 * y0 * z0;
                h[1] += e * x0 * y2 * z0;
                h[2] += e * x0 * y0 * z2;
                h[3] += e * x1 * y1 * z0;
                h[4] += e * x1 * y0 * z1;
                h[5] += e * x0 * y1 * z1;
            }
            val[o] = scale * v;
            for k in 0..3 {
                grad[3 * o + k] = scale * g[k];
            }
            for k in 0..6 {
                hess[6 * o + k] = scale * h[k];
            }
            o += 1;
        }
    }
}
