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
