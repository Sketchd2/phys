//! Non-local correlation: the van der Waals density functional's kernel —
//! `docs/PLAY.md` Phase 6, E7.
//!
//! # What it is
//!
//! A semilocal functional sees the density at a point and its gradient, and
//! so cannot see two molecules attract each other across a gap where neither
//! has any density: dispersion is a correlation between fluctuations in two
//! places at once. Dion, Rydberg, Schroder, Langreth and Lundqvist (2004)
//! wrote that correlation as a double integral over the density,
//!
//! ```text
//! E_c^nl = 1/2 ∫∫ n(r) phi(q0(r) |r - r'|, q0(r') |r - r'|) n(r') dr dr'
//! ```
//!
//! with a universal kernel `phi(d1, d2)` and a local wavevector `q0` set by
//! the density, its gradient and the uniform gas's correlation energy there.
//!
//! # What is derived and what is stored
//!
//! The kernel comes from a plasmon-pole model of the electron gas's response,
//! and nothing in it was fitted to a molecule: `gamma = 4 pi / 9` follows from
//! the model's small-wavevector limit, and the one constant in `q0`,
//! `Z_ab = -0.8491`, is the gradient expansion of the slowly varying gas's
//! exchange. It is a function of two numbers, computed here by numerical
//! integration of Dion's formula, once, and stored as a table — the shortcut —
//! rather than read from a published one. The exchange it is paired with is
//! the empirical part, and is chosen elsewhere (PLAY.md E7).
//!
//! The formula, eq. 14 of Dion et al., with `T`'s factor of one half folded
//! into the prefactor:
//!
//! ```text
//! phi(d1, d2) = 1/pi^2 ∫0^∞ ∫0^∞ a^2 b^2 W(a, b) T(nu(a), nu(b), nu'(a), nu'(b)) da db
//! W(a, b)     = 2 [(3 - a^2) b cos b sin a + (3 - b^2) a cos a sin b
//!                  + (a^2 + b^2 - 3) sin a sin b - 3 a b cos a cos b] / (a^3 b^3)
//! T(w,x,y,z)  = [1/(w + x) + 1/(y + z)] [1/((w + y)(x + z)) + 1/((w + z)(y + x))]
//! nu(y)       = y^2 / (2 h(y / d1)),  nu'(y) = y^2 / (2 h(y / d2)),
//! h(y)        = 1 - exp(-gamma y^2)
//! ```
//!
//! checked against Quantum ESPRESSO's kernel generator
//! (`generate_vdW_kernel_table.f90`), which writes the same integrand.

use std::f64::consts::PI;

/// `gamma` in the kernel's `h(y) = 1 - exp(-gamma y^2)`: `4 pi / 9`.
pub const GAMMA: f64 = 4.0 * PI / 9.0;

/// The gradient coefficient in `q0` for vdW-DF1, from the slowly varying
/// electron gas's exchange: `Z_ab = -0.8491`.
pub const Z_AB_DF1: f64 = -0.8491;

/// The kernel's large-separation limit: `phi -> -C / (d1^2 d2^2 (d1^2 + d2^2))`
/// with `C = 12 (4 pi / 9)^3`, which is what makes the energy between two
/// distant pieces of density fall as `R^-6`.
pub const ASYMPTOTE_C: f64 = 12.0 * GAMMA * GAMMA * GAMMA;

/// A quadrature for the kernel's double integral over `a` and `b` on
/// `[0, a_max]`: Gauss-Legendre in `atan(a)`, which crowds points near the
/// origin where `nu` changes fastest and thins them where the integrand only
/// oscillates and decays.
pub struct KernelQuadrature {
    a: Vec<f64>,
    /// `W(a_i, a_j) a_i^2 a_j^2` times both weights, the part of the integrand
    /// that does not depend on `d1` or `d2`.
    w: Vec<f64>,
}

impl KernelQuadrature {
    pub fn new(points: usize, a_max: f64) -> KernelQuadrature {
        let (x, wx) = super::grid::gauss_legendre(points);
        let (lo, hi) = (0f64.atan(), a_max.atan());
        let (mid, half) = (0.5 * (lo + hi), 0.5 * (hi - lo));
        let a: Vec<f64> = x.iter().map(|t| (mid + half * t).tan()).collect();
        // d a = (1 + a^2) d atan(a).
        let wa: Vec<f64> = wx.iter().zip(&a).map(|(w, a)| w * half * (1.0 + a * a)).collect();
        let n = a.len();
        let mut w = vec![0.0; n * n];
        for i in 0..n {
            let (ai, si, ci) = (a[i], a[i].sin(), a[i].cos());
            for j in 0..n {
                let (aj, sj, cj) = (a[j], a[j].sin(), a[j].cos());
                // a^2 b^2 W(a, b) = 2 [...] / (a b).
                let bracket = (3.0 - ai * ai) * aj * cj * si + (3.0 - aj * aj) * ai * ci * sj + (ai * ai + aj * aj - 3.0) * si * sj - 3.0 * ai * aj * ci * cj;
                w[i * n + j] = wa[i] * wa[j] * 2.0 * bracket / (ai * aj);
            }
        }
        KernelQuadrature { a, w }
    }

    /// `phi(d1, d2)` by this quadrature.
    pub fn phi(&self, d1: f64, d2: f64) -> f64 {
        if d1 == 0.0 && d2 == 0.0 {
            return 0.0;
        }
        let nu = |y: f64, d: f64| -> f64 {
            if d == 0.0 {
                // h -> 1 as d -> 0 at fixed y.
                0.5 * y * y
            } else {
                let h = 1.0 - (-GAMMA * y * y / (d * d)).exp();
                if h > 0.0 { 0.5 * y * y / h } else { 0.5 * d * d / GAMMA }
            }
        };
        let n = self.a.len();
        let v1: Vec<f64> = self.a.iter().map(|&y| nu(y, d1)).collect();
        let v2: Vec<f64> = self.a.iter().map(|&y| nu(y, d2)).collect();
        let mut sum = 0.0;
        for i in 0..n {
            let (w_, y) = (v1[i], v2[i]);
            for j in 0..n {
                let (x, z) = (v1[j], v2[j]);
                let t = (1.0 / (w_ + x) + 1.0 / (y + z)) * (1.0 / ((w_ + y) * (x + z)) + 1.0 / ((w_ + z) * (y + x)));
                sum += t * self.w[i * n + j];
            }
        }
        sum / (PI * PI)
    }
}

/// The kernel's large-separation form.
pub fn phi_asymptotic(d1: f64, d2: f64) -> f64 {
    -ASYMPTOTE_C / (d1 * d1 * d2 * d2 * (d1 * d1 + d2 * d2))
}

/// Where both arguments are at least this, the kernel is its asymptote:
/// measured with the converged quadrature, `phi / asymptote` is 1.0000 at
/// `(12, 12)` and 0.9999 at `(8, 12)` and `(10, 15)`, and the quadrature
/// itself turns noisy on the vanishing values beyond. With one argument
/// small it is nothing like the asymptote — 0.005 at `(23.8, 0.24)` — so this
/// is a test on the smaller argument.
pub const ASYMPTOTIC_FROM: f64 = 12.0;

/// The quadrature the stored kernel is computed with: 1024 points to
/// `a_max = 256`. Measured against 4096 to 512, the kernel agrees to about
/// 1e-9 for `d` up to 6; with too few points for its range a quadrature
/// under-resolves the oscillating integrand (1024 to 1024 is off by 1e-6),
/// so the points scale with the range.
pub fn converged_quadrature() -> KernelQuadrature {
    KernelQuadrature::new(1024, 256.0)
}

/// The kernel tabulated on a square grid in `u = d / (1 + d)`, for each
/// argument from 0 to `d_max`, with bicubic interpolation between nodes.
pub struct KernelTable {
    /// Nodes per side.
    pub m: usize,
    pub d_max: f64,
    /// `phi` at `(u_i, u_j)`, row-major, symmetric.
    pub values: Vec<f64>,
}

impl KernelTable {
    fn u_max(d_max: f64) -> f64 {
        d_max / (1.0 + d_max)
    }

    /// The argument at node `i`.
    pub fn node(m: usize, d_max: f64, i: usize) -> f64 {
        let u = Self::u_max(d_max) * i as f64 / (m - 1) as f64;
        u / (1.0 - u)
    }

    /// Compute the table: `m (m + 1) / 2` kernel values, by symmetry, spread
    /// across threads.
    pub fn build(m: usize, d_max: f64, q: &KernelQuadrature) -> KernelTable {
        let pairs: Vec<(usize, usize)> = (0..m).flat_map(|i| (0..=i).map(move |j| (i, j))).collect();
        let job = |idx: &mut dyn Iterator<Item = usize>| {
            idx.map(|k| {
                let (i, j) = pairs[k];
                (i, j, q.phi(Self::node(m, d_max, i), Self::node(m, d_max, j)))
            }).collect::<Vec<_>>()
        };
        let mut values = vec![0.0; m * m];
        for part in super::scf::parallel_interleaved(pairs.len(), &job) {
            for (i, j, v) in part {
                values[i * m + j] = v;
                values[j * m + i] = v;
            }
        }
        KernelTable { m, d_max, values }
    }

    /// `phi(d1, d2)`: the asymptote where both are at least
    /// [`ASYMPTOTIC_FROM`], otherwise the table, with an argument beyond
    /// `d_max` held at `d_max` (where the kernel is of order 1e-7; what that
    /// does to an energy is measured by moving `d_max`).
    pub fn phi(&self, d1: f64, d2: f64) -> f64 {
        if d1.min(d2) >= ASYMPTOTIC_FROM {
            return phi_asymptotic(d1, d2);
        }
        let (d1, d2) = (d1.min(self.d_max), d2.min(self.d_max));
        let scale = (self.m - 1) as f64 / Self::u_max(self.d_max);
        let (x, y) = (d1 / (1.0 + d1) * scale, d2 / (1.0 + d2) * scale);
        self.bicubic(x, y)
    }

    /// Catmull-Rom bicubic interpolation at fractional node position `(x, y)`,
    /// with the grid continued past its edges by linear extrapolation.
    fn bicubic(&self, x: f64, y: f64) -> f64 {
        let m = self.m;
        let last = (m - 2) as f64;
        let (x, y) = (x.clamp(0.0, (m - 1) as f64), y.clamp(0.0, (m - 1) as f64));
        let (i, j) = (x.floor().min(last) as isize, y.floor().min(last) as isize);
        let (tx, ty) = (x - i as f64, y - j as f64);
        let at = |a: isize, b: isize| -> f64 {
            // Linear extrapolation one node beyond either edge.
            let ext = |k: isize| -> (usize, usize, f64) {
                if k < 0 {
                    (0, 1, -1.0)
                } else if k as usize >= m {
                    (m - 1, m - 2, -1.0)
                } else {
                    (k as usize, k as usize, 0.0)
                }
            };
            let (a0, a1, fa) = ext(a);
            let (b0, b1, fb) = ext(b);
            let v = |p: usize, q: usize| self.values[p * m + q];
            let base = v(a0, b0);
            let da = if fa != 0.0 { base - v(a1, b0) } else { 0.0 };
            let db = if fb != 0.0 { base - v(a0, b1) } else { 0.0 };
            base + da + db
        };
        let w = |t: f64| -> [f64; 4] {
            let (t2, t3) = (t * t, t * t * t);
            [-0.5 * t3 + t2 - 0.5 * t, 1.5 * t3 - 2.5 * t2 + 1.0, -1.5 * t3 + 2.0 * t2 + 0.5 * t, 0.5 * t3 - 0.5 * t2]
        };
        let (wx, wy) = (w(tx), w(ty));
        let mut s = 0.0;
        for (p, wxp) in wx.iter().enumerate() {
            for (q, wyq) in wy.iter().enumerate() {
                s += wxp * wyq * at(i - 1 + p as isize, j - 1 + q as isize);
            }
        }
        s
    }
}
