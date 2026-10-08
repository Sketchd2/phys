//! Boys localisation of occupied orbitals — `docs/PLAY.md` E8b.
//!
//! The canonical orbitals of a molecule spread over all of it. Rotating the
//! occupied ones among themselves leaves the density, the energy and
//! everything measurable unchanged, and the rotation that makes each orbital
//! as compact as it can be (Foster and Boys: maximise the sum of the squared
//! distances of the orbitals' centroids from the origin) turns them into what
//! a chemist draws: a bond, a lone pair, each with a place. Those places are
//! derived from the molecule's own electrons, and a force field's extra sites
//! — the lone pairs a water's hydrogen bonds point at — can be put there
//! instead of being guessed at (`phys-esp --sites`).

use super::linalg::Matrix;

/// The rotated orbitals and where each one's centre of charge is.
pub struct Localised {
    /// All the orbitals in `c`'s layout (`c[ao * stride + orbital]`), those in
    /// the localised range rotated among themselves, the rest untouched.
    pub c: Vec<f64>,
    /// The centroid `<i|r|i>` of each localised orbital, bohr.
    pub centroids: Vec<[f64; 3]>,
    /// Sweeps it took.
    pub sweeps: usize,
}

/// Localise orbitals `range` of `c` (`n` basis functions, `stride` orbitals
/// to a row) with the position matrices `dip` (`<m|r_k|n>`, about any one
/// origin: the localisation does not depend on it).
pub fn boys(c: &[f64], stride: usize, n: usize, range: std::ops::Range<usize>, dip: &[Matrix; 3]) -> Localised {
    let m = range.len();
    let mut c = c.to_vec();
    // r[k][i][j] over the localised orbitals.
    let mut r = vec![vec![vec![0.0f64; m]; m]; 3];
    for k in 0..3 {
        for i in 0..m {
            for j in 0..m {
                let mut s = 0.0;
                for a in 0..n {
                    let ca = c[a * stride + range.start + i];
                    if ca == 0.0 {
                        continue;
                    }
                    for b in 0..n {
                        s += ca * c[b * stride + range.start + j] * dip[k].get(a, b);
                    }
                }
                r[k][i][j] = s;
            }
        }
    }
    let mut sweeps = 0;
    for sweep in 0..200 {
        sweeps = sweep + 1;
        let mut biggest: f64 = 0.0;
        for i in 0..m {
            for j in i + 1..m {
                let (mut a_sum, mut b_sum) = (0.0, 0.0);
                for k in 0..3 {
                    let d = r[k][i][i] - r[k][j][j];
                    a_sum += r[k][i][j].powi(2) - 0.25 * d * d;
                    b_sum += r[k][i][j] * d;
                }
                let theta = 0.25 * b_sum.atan2(-a_sum);
                if theta.abs() < 1e-14 {
                    continue;
                }
                biggest = biggest.max(theta.abs());
                let (s, co) = theta.sin_cos();
                // Orbitals i and j: i' = c i + s j, j' = -s i + c j.
                for a in 0..n {
                    let (x, y) = (c[a * stride + range.start + i], c[a * stride + range.start + j]);
                    c[a * stride + range.start + i] = co * x + s * y;
                    c[a * stride + range.start + j] = -s * x + co * y;
                }
                for k in 0..3 {
                    for l in 0..m {
                        let (x, y) = (r[k][i][l], r[k][j][l]);
                        r[k][i][l] = co * x + s * y;
                        r[k][j][l] = -s * x + co * y;
                    }
                    for l in 0..m {
                        let (x, y) = (r[k][l][i], r[k][l][j]);
                        r[k][l][i] = co * x + s * y;
                        r[k][l][j] = -s * x + co * y;
                    }
                }
            }
        }
        if biggest < 1e-12 {
            break;
        }
    }
    let centroids = (0..m).map(|i| [r[0][i][i], r[1][i][i], r[2][i][i]]).collect();
    Localised { c, centroids, sweeps }
}
