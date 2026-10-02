//! A quadrature grid over a molecule, for the integrals that have no closed
//! form — above all the exchange-correlation energy.
//!
//! **Generated, not tabulated.** The usual angular rules (Lebedev's) are tables
//! of points and weights; the ones here are computed: Gauss-Legendre nodes in
//! `cos(theta)`, found by Newton's method on the Legendre polynomials, times an
//! even spread in `phi`. That product is exact for every spherical harmonic up
//! to degree `2 n_theta - 1` and costs about half as many points again as
//! Lebedev for the same degree, which is the price of not storing a table.
//! Radially, Gauss-Chebyshev of the second kind with Becke's map
//! `r = R (1 + x) / (1 - x)`. Each atom's grid is then weighted by Becke's fuzzy
//! partition so that the molecule is the sum of atomic integrals and nothing is
//! counted twice.

/// Points, bohr, and weights.
#[derive(Debug, Clone, Default)]
pub struct Grid {
    pub points: Vec<[f64; 3]>,
    pub weights: Vec<f64>,
    /// The atom each point belongs to and moves with.
    pub owner: Vec<usize>,
    /// Each point's radial-times-angular weight, before Becke's partition:
    /// `weights = raw * becke_share(owner, point)`.
    pub raw: Vec<f64>,
}

/// Atom `owner`'s share of space at `p`: Becke's fuzzy cell, `P_A = Z_A / sum_B
/// Z_B`, `Z_B = prod_(C != B) s(mu_BC)`.
pub fn becke_share(atoms: &[[f64; 3]], owner: usize, p: [f64; 3]) -> f64 {
    let n = atoms.len();
    if n == 1 {
        return 1.0;
    }
    let dist = |a: [f64; 3], b: [f64; 3]| ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
    let mut cell = vec![1.0; n];
    for i in 0..n {
        let ri = dist(p, atoms[i]);
        for j in 0..n {
            if i == j {
                continue;
            }
            let rj = dist(p, atoms[j]);
            cell[i] *= becke_s((ri - rj) / dist(atoms[i], atoms[j]));
        }
    }
    let total: f64 = cell.iter().sum();
    if total > 0.0 { cell[owner] / total } else { 0.0 }
}

/// How `becke_share` at a point of atom `owner` — a point that moves with that
/// atom, `p = R_owner + u` — changes as each atom moves: by central differences
/// of the share, which is a smooth function of positions alone.
pub fn becke_share_gradient(atoms: &[[f64; 3]], owner: usize, p: [f64; 3]) -> Vec<[f64; 3]> {
    let h = 1e-5;
    let n = atoms.len();
    let mut out = vec![[0.0; 3]; n];
    if n == 1 {
        return out;
    }
    let u = [p[0] - atoms[owner][0], p[1] - atoms[owner][1], p[2] - atoms[owner][2]];
    for d in 0..n {
        for k in 0..3 {
            let mut plus = atoms.to_vec();
            let mut minus = atoms.to_vec();
            plus[d][k] += h;
            minus[d][k] -= h;
            let at = |a: &[[f64; 3]]| [a[owner][0] + u[0], a[owner][1] + u[1], a[owner][2] + u[2]];
            out[d][k] = (becke_share(&plus, owner, at(&plus)) - becke_share(&minus, owner, at(&minus))) / (2.0 * h);
        }
    }
    out
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    let pi = std::f64::consts::PI;
    for i in 0..n.div_ceil(2) {
        // Tricomi's estimate of the root, then Newton.
        let mut z = (pi * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, z);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let pn = if n == 0 { 1.0 } else if n == 1 { z } else { p1 };
            let pnm1 = if n == 1 { 1.0 } else { p0 };
            dp = n as f64 * (z * pn - pnm1) / (z * z - 1.0);
            let dz = pn / dp;
            z -= dz;
            if dz.abs() < 1e-16 {
                break;
            }
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        let wi = 2.0 / ((1.0 - z * z) * dp * dp);
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
    (x, w)
}

/// Points on the unit sphere and weights summing to `4 pi`, exact for
/// spherical harmonics up to degree `2 n_theta - 1`.
pub fn sphere(n_theta: usize) -> (Vec<[f64; 3]>, Vec<f64>) {
    let (ct, wt) = gauss_legendre(n_theta);
    let n_phi = 2 * n_theta;
    let dphi = 2.0 * std::f64::consts::PI / n_phi as f64;
    let mut pts = Vec::with_capacity(n_theta * n_phi);
    let mut ws = Vec::with_capacity(n_theta * n_phi);
    for (c, w) in ct.iter().zip(&wt) {
        let s = (1.0 - c * c).max(0.0).sqrt();
        for k in 0..n_phi {
            let phi = (k as f64 + 0.5) * dphi;
            pts.push([s * phi.cos(), s * phi.sin(), *c]);
            ws.push(w * dphi);
        }
    }
    (pts, ws)
}

/// Radial points and weights for `integral_0^inf f(r) r^2 dr`, scale `scale`.
pub fn radial(n: usize, scale: f64) -> (Vec<f64>, Vec<f64>) {
    let pi = std::f64::consts::PI;
    let mut r = Vec::with_capacity(n);
    let mut w = Vec::with_capacity(n);
    for i in 1..=n {
        let a = i as f64 * pi / (n + 1) as f64;
        let x = a.cos();
        let s = a.sin();
        // Chebyshev second kind integrates g(x) sqrt(1 - x^2); divide it out.
        let wx = pi / (n + 1) as f64 * s * s / s;
        let ri = scale * (1.0 + x) / (1.0 - x);
        let drdx = 2.0 * scale / ((1.0 - x) * (1.0 - x));
        r.push(ri);
        w.push(wx * drdx * ri * ri);
    }
    (r, w)
}

/// Becke's cell function: 1 deep inside atom `i`'s cell, 0 deep inside
/// another's, smooth between.
fn becke_s(mu: f64) -> f64 {
    let mut f = mu;
    for _ in 0..3 {
        f = 1.5 * f - 0.5 * f * f * f;
    }
    0.5 * (1.0 - f)
}

/// A molecular grid. `atoms` are positions with a size each (bohr), which is
/// the scale of that atom's radial map.
///
/// **Becke's size adjustment is deliberately not used.** It moves the boundary
/// between unlike atoms towards the smaller one, and in doing so sharpens the
/// partition exactly where an angular rule integrates it worst: on water, the
/// adjustment made the grid's overlap error 27x larger at the same points
/// (2.5e-5 against 9.3e-7). The size is used where it helps, in the radial map.
pub fn molecular(atoms: &[([f64; 3], f64)], n_radial: usize, n_theta: usize) -> Grid {
    molecular_pruned(atoms, n_radial, n_theta, false)
}

/// The angular order used at radius `r` on an atom of size `size`, pruned near
/// the nucleus, where the density is nearly spherical.
///
/// **Not far out**, although the density is smooth there too: measured on
/// water with the derived basis, halving the order beyond five atomic sizes
/// moved the energy 1.3e-5 Ha, because hydrogen's diffuse d functions (down to
/// an exponent of 0.011) still have angular structure at three bohr. Pruning
/// the core alone moves it under 1e-6.
pub fn pruned_theta(n_theta: usize, r: f64, size: f64) -> usize {
    let x = r / size;
    let f = if x < 0.25 {
        1.0 / 3.0
    } else if x < 0.5 {
        0.5
    } else {
        1.0
    };
    ((n_theta as f64 * f).ceil() as usize).max(4).min(n_theta)
}

/// As [`molecular`], optionally pruned by [`pruned_theta`].
pub fn molecular_pruned(atoms: &[([f64; 3], f64)], n_radial: usize, n_theta: usize, prune: bool) -> Grid {
    let mut spheres: std::collections::HashMap<usize, (Vec<[f64; 3]>, Vec<f64>)> = std::collections::HashMap::new();
    let mut grid = Grid::default();
    let positions: Vec<[f64; 3]> = atoms.iter().map(|a| a.0).collect();
    for (ia, (centre, size)) in atoms.iter().enumerate() {
        let (rr, rw) = radial(n_radial, *size);
        for (r, wr) in rr.iter().zip(&rw) {
            let nt = if prune { pruned_theta(n_theta, *r, *size) } else { n_theta };
            let (sph, sw) = spheres.entry(nt).or_insert_with(|| sphere(nt));
            for (d, wa) in sph.iter().zip(sw.iter()) {
                let p = [centre[0] + r * d[0], centre[1] + r * d[1], centre[2] + r * d[2]];
                let w = becke_share(&positions, ia, p);
                let weight = wr * wa * w;
                if weight > 1e-20 {
                    grid.points.push(p);
                    grid.weights.push(weight);
                    grid.owner.push(ia);
                    grid.raw.push(wr * wa);
                }
            }
        }
    }
    grid
}
