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

    /// `phi(d1, d2)` and its partial derivatives with respect to `d1` and
    /// `d2`, of exactly the function [`KernelTable::phi`] evaluates: the
    /// interpolant's own derivatives inside the table, the asymptote's where
    /// it is used, and zero along an argument held at `d_max`.
    pub fn phi_and_slopes(&self, d1: f64, d2: f64) -> (f64, f64, f64) {
        if d1.min(d2) >= ASYMPTOTIC_FROM {
            let f = phi_asymptotic(d1, d2);
            let s = d1 * d1 + d2 * d2;
            return (f, f * (-2.0 / d1 - 2.0 * d1 / s), f * (-2.0 / d2 - 2.0 * d2 / s));
        }
        let (c1, c2) = (d1.min(self.d_max), d2.min(self.d_max));
        let scale = (self.m - 1) as f64 / Self::u_max(self.d_max);
        let (x, y) = (c1 / (1.0 + c1) * scale, c2 / (1.0 + c2) * scale);
        let (f, fx, fy) = self.bicubic_with_slopes(x, y);
        // dx / dd = scale / (1 + d)^2 inside the table, nothing past its edge.
        let gx = if d1 < self.d_max { scale / ((1.0 + c1) * (1.0 + c1)) } else { 0.0 };
        let gy = if d2 < self.d_max { scale / ((1.0 + c2) * (1.0 + c2)) } else { 0.0 };
        (f, fx * gx, fy * gy)
    }

    fn bicubic(&self, x: f64, y: f64) -> f64 {
        self.bicubic_with_slopes(x, y).0
    }

    /// Catmull-Rom bicubic interpolation at fractional node position `(x, y)`,
    /// with the grid continued past its edges by linear extrapolation.
    fn bicubic_with_slopes(&self, x: f64, y: f64) -> (f64, f64, f64) {
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
        let dw = |t: f64| -> [f64; 4] {
            let t2 = t * t;
            [-1.5 * t2 + 2.0 * t - 0.5, 4.5 * t2 - 5.0 * t, -4.5 * t2 + 4.0 * t + 0.5, 1.5 * t2 - t]
        };
        let (wx, wy, dwx, dwy) = (w(tx), w(ty), dw(tx), dw(ty));
        let (mut s, mut sx, mut sy) = (0.0, 0.0, 0.0);
        for p in 0..4 {
            for q in 0..4 {
                let v = at(i - 1 + p as isize, j - 1 + q as isize);
                s += wx[p] * wy[q] * v;
                sx += dwx[p] * wy[q] * v;
                sy += wx[p] * dwy[q] * v;
            }
        }
        (s, sx, sy)
    }
}

/// The electron density and the square of its gradient at every point of a
/// grid, for the total of both spins: what the non-local functional reads.
/// One entry per grid point, in the grid's order.
pub fn density_on(basis: &super::basis::Basis, batches: &super::scf::Batches, d_total: &super::linalg::Matrix, points: usize) -> Vec<(f64, f64)> {
    density_and_gradient_on(basis, batches, d_total, points).iter().map(|(n, g)| (*n, g[0] * g[0] + g[1] * g[1] + g[2] * g[2])).collect()
}

/// The total density and its gradient vector at every point of a grid.
pub fn density_and_gradient_on(basis: &super::basis::Basis, batches: &super::scf::Batches, d_total: &super::linalg::Matrix, points: usize) -> Vec<(f64, [f64; 3])> {
    let n = basis.size;
    let job = |idx: &mut dyn Iterator<Item = usize>| {
        let mut out: Vec<(usize, f64, [f64; 3])> = Vec::new();
        for bi in idx {
            let (pts, _, funcs, shells) = &batches.batches[bi];
            let k = funcs.len();
            let mut sub = vec![0.0; k * k];
            for (i, &fi) in funcs.iter().enumerate() {
                for (j, &fj) in funcs.iter().enumerate() {
                    sub[i * k + j] = d_total.a[fi * n + fj];
                }
            }
            let mut phi = vec![0.0; k];
            let mut g3 = vec![0.0; 3 * k];
            let mut x = vec![0.0; k];
            for (p, pt) in pts.iter().enumerate() {
                super::values::at_shells(basis, shells, *pt, &mut phi, Some(&mut g3));
                for i in 0..k {
                    x[i] = (0..k).map(|j| sub[i * k + j] * phi[j]).sum();
                }
                let rho: f64 = (0..k).map(|i| phi[i] * x[i]).sum();
                let mut g = [0.0; 3];
                for (dir, gd) in g.iter_mut().enumerate() {
                    *gd = 2.0 * (0..k).map(|i| g3[3 * i + dir] * x[i]).sum::<f64>();
                }
                out.push((batches.indices[bi][p], rho, g));
            }
        }
        out
    };
    let mut values = vec![(0.0, [0.0; 3]); points];
    for part in super::scf::parallel_interleaved(batches.batches.len(), &job) {
        for (i, rho, g) in part {
            values[i] = (rho, g);
        }
    }
    values
}

/// `q0`, the local wavevector the kernel's arguments are scaled by (Dion et
/// al. eq. 11-12): `q0 = -(4 pi / 3) eps_xc^0`, with `eps_xc^0` the uniform
/// gas's exchange-correlation energy per electron corrected by the gradient
/// expansion of exchange, `- eps_x^LDA (Z_ab / 9) s^2`, `s = |grad n| / (2 k_F n)`.
pub fn q0(n: f64, grad2: f64, z_ab: f64) -> f64 {
    let kf = (3.0 * PI * PI * n).powf(1.0 / 3.0);
    let s2 = grad2 / (4.0 * kf * kf * n * n);
    let ex = -3.0 * kf / (4.0 * PI);
    let ec = super::functional::lda_correlation_per_electron(n);
    -(4.0 * PI / 3.0) * (ex + ec - ex * (z_ab / 9.0) * s2)
}

/// [`q0`] and its derivatives with respect to `n` and to `|grad n|^2`.
/// Written out, `q0 = k_F - (4 pi / 3) eps_c - Z_ab |grad n|^2 / (36 k_F n^2)`.
pub fn q0_and_slopes(n: f64, grad2: f64, z_ab: f64) -> (f64, f64, f64) {
    let kf = (3.0 * PI * PI * n).powf(1.0 / 3.0);
    let (ec, dec) = super::functional::lda_correlation_per_electron_and_slope(n);
    let q = kf - (4.0 * PI / 3.0) * ec - z_ab * grad2 / (36.0 * kf * n * n);
    let dkf = kf / (3.0 * n);
    let dq_dn = dkf - (4.0 * PI / 3.0) * dec + 7.0 * z_ab * grad2 / (108.0 * kf * n * n * n);
    let dq_dg2 = -z_ab / (36.0 * kf * n * n);
    (q, dq_dn, dq_dg2)
}

/// One point of density as the non-local energy sees it: where it is, its
/// weight times its density, and its `q0`.
#[derive(Debug, Clone, Copy)]
pub struct Site {
    pub r: [f64; 3],
    pub wn: f64,
    pub q: f64,
}

/// The sites of a density on a grid, dropping points whose weighted density
/// is below `floor` — what that drops is measured by moving `floor`.
pub fn sites(grid: &super::grid::Grid, density: &[(f64, f64)], z_ab: f64, floor: f64) -> Vec<Site> {
    grid.points.iter().zip(&grid.weights).zip(density).filter_map(|((r, w), (n, g2))| {
        let wn = w * n;
        (*n > 0.0 && wn.abs() >= floor).then(|| Site { r: *r, wn, q: q0(*n, *g2, z_ab) })
    }).collect()
}

/// `E_c^nl = 1/2 sum_i sum_j w_i n_i w_j n_j phi(q_i r_ij, q_j r_ij)`, every
/// pair once and the diagonal (`r = 0`, `phi(0, 0)`) included as the product
/// quadrature has it. Rows are spread across threads, and each row's sum is
/// its own, added in row order, so the result does not depend on the thread
/// count.
pub fn nonlocal_energy(sites: &[Site], table: &KernelTable) -> f64 {
    let job = |idx: &mut dyn Iterator<Item = usize>| {
        idx.map(|i| {
            let a = &sites[i];
            let mut row = 0.5 * a.wn * a.wn * table.phi(0.0, 0.0);
            for b in &sites[i + 1..] {
                let d = [a.r[0] - b.r[0], a.r[1] - b.r[1], a.r[2] - b.r[2]];
                let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                row += a.wn * b.wn * table.phi(a.q * r, b.q * r);
            }
            (i, row)
        }).collect::<Vec<_>>()
    };
    let mut rows = vec![0.0; sites.len()];
    for part in super::scf::parallel_interleaved(sites.len(), &job) {
        for (i, v) in part {
            rows[i] = v;
        }
    }
    rows.iter().sum()
}

/// The `C6` the non-local energy gives two systems far apart:
/// `E -> -C6 / R^6` with `C6 = C sum_(i in A) sum_(j in B) w_i n_i w_j n_j /
/// (q_i^2 q_j^2 (q_i^2 + q_j^2))`, from the kernel's asymptote.
pub fn c6(a: &[Site], b: &[Site]) -> f64 {
    let job = |idx: &mut dyn Iterator<Item = usize>| {
        idx.map(|i| {
            let s = &a[i];
            let q2 = s.q * s.q;
            (i, b.iter().map(|t| { let p2 = t.q * t.q; s.wn * t.wn / (q2 * p2 * (q2 + p2)) }).sum::<f64>())
        }).collect::<Vec<_>>()
    };
    let mut rows = vec![0.0; a.len()];
    for part in super::scf::parallel_interleaved(a.len(), &job) {
        for (i, v) in part {
            rows[i] = v;
        }
    }
    ASYMPTOTE_C * rows.iter().sum::<f64>()
}

/// The non-local energy and its derivatives at every grid point: what a
/// self-consistent field and its forces need.
pub struct Nonlocal {
    pub energy: f64,
    /// `dE / dn` at each grid point, per unit weight.
    pub v_n: Vec<f64>,
    /// `dE / d|grad n|^2` at each grid point, per unit weight.
    pub v_g2: Vec<f64>,
}

/// [`nonlocal_energy`] with its derivatives. For each kept point `i`, over
/// every kept point `j` (itself included, at `phi(0, 0)`):
/// `A_i = sum_j w_j n_j phi(q_i r, q_j r)` and
/// `B_i = sum_j w_j n_j r d1phi(q_i r, q_j r)`; then `E = 1/2 sum_i w_i n_i A_i`,
/// `dE/dn_i = w_i (A_i + n_i B_i dq_i/dn)` and
/// `dE/d|grad n|^2_i = w_i n_i B_i dq_i/d|grad n|^2`. Each row is its own sum
/// over every point, twice the kernel evaluations of visiting each pair once,
/// so that rows can run on any thread and the result does not depend on how
/// many there are.
pub fn nonlocal(grid: &super::grid::Grid, density: &[(f64, [f64; 3])], z_ab: f64, floor: f64, table: &KernelTable) -> Nonlocal {
    struct Point { at: usize, r: [f64; 3], wn: f64, n: f64, q: f64, dq_dn: f64, dq_dg2: f64 }
    let points: Vec<Point> = grid.points.iter().zip(&grid.weights).zip(density).enumerate().filter_map(|(at, ((r, w), (n, g)))| {
        let wn = w * n;
        if *n <= 0.0 || wn.abs() < floor {
            return None;
        }
        let (q, dq_dn, dq_dg2) = q0_and_slopes(*n, g[0] * g[0] + g[1] * g[1] + g[2] * g[2], z_ab);
        Some(Point { at, r: *r, wn, n: *n, q, dq_dn, dq_dg2 })
    }).collect();
    let job = |idx: &mut dyn Iterator<Item = usize>| {
        idx.map(|i| {
            let a = &points[i];
            let (mut sa, mut sb) = (0.0, 0.0);
            for b in &points {
                let d = [a.r[0] - b.r[0], a.r[1] - b.r[1], a.r[2] - b.r[2]];
                let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                let (f, f1, _) = table.phi_and_slopes(a.q * r, b.q * r);
                sa += b.wn * f;
                sb += b.wn * r * f1;
            }
            (i, sa, sb)
        }).collect::<Vec<_>>()
    };
    let mut ab = vec![(0.0, 0.0); points.len()];
    for part in super::scf::parallel_interleaved(points.len(), &job) {
        for (i, sa, sb) in part {
            ab[i] = (sa, sb);
        }
    }
    let mut energy = 0.0;
    let mut v_n = vec![0.0; grid.points.len()];
    let mut v_g2 = vec![0.0; grid.points.len()];
    for (p, (sa, sb)) in points.iter().zip(&ab) {
        energy += 0.5 * p.wn * sa;
        v_n[p.at] = sa + p.n * sb * p.dq_dn;
        v_g2[p.at] = p.n * sb * p.dq_dg2;
    }
    Nonlocal { energy, v_n, v_g2 }
}

/// The non-local correlation's Kohn-Sham matrix, the same for both spins:
/// `V_mn = sum_p w_p [v_n phi_m phi_n + 2 v_g2 grad n . grad(phi_m phi_n)]`.
/// Batches run on any thread and are added in batch order.
pub fn nonlocal_matrix(basis: &super::basis::Basis, batches: &super::scf::Batches, density: &[(f64, [f64; 3])], nl: &Nonlocal) -> super::linalg::Matrix {
    let n = basis.size;
    let job = |idx: &mut dyn Iterator<Item = usize>| {
        idx.map(|bi| {
            let (pts, ws, funcs, shells) = &batches.batches[bi];
            let k = funcs.len();
            let mut block = vec![0.0; k * k];
            let mut phi = vec![0.0; k];
            let mut g3 = vec![0.0; 3 * k];
            let mut c = vec![0.0; k];
            for (p, pt) in pts.iter().enumerate() {
                let at = batches.indices[bi][p];
                let (vn, vg) = (nl.v_n[at], nl.v_g2[at]);
                if vn == 0.0 && vg == 0.0 {
                    continue;
                }
                let gn = density[at].1;
                super::values::at_shells(basis, shells, *pt, &mut phi, Some(&mut g3));
                // V = sum_p w_p (phi_m c_n + c_m phi_n), c = v_n phi / 2 + 2 v_g2 grad n . grad phi.
                for m in 0..k {
                    c[m] = 0.5 * vn * phi[m] + 2.0 * vg * (gn[0] * g3[3 * m] + gn[1] * g3[3 * m + 1] + gn[2] * g3[3 * m + 2]);
                }
                let w = ws[p];
                for m in 0..k {
                    for l in 0..k {
                        block[m * k + l] += w * (phi[m] * c[l] + c[m] * phi[l]);
                    }
                }
            }
            (bi, block)
        }).collect::<Vec<_>>()
    };
    let mut parts: Vec<(usize, Vec<f64>)> = super::scf::parallel_interleaved(batches.batches.len(), &job).into_iter().flatten().collect();
    parts.sort_by_key(|(bi, _)| *bi);
    let mut v = super::linalg::Matrix::zeros(n);
    for (bi, block) in parts {
        let funcs = &batches.batches[bi].2;
        let k = funcs.len();
        for (m, &fm) in funcs.iter().enumerate() {
            for (l, &fl) in funcs.iter().enumerate() {
                v.a[fm * n + fl] += block[m * k + l];
            }
        }
    }
    v
}
