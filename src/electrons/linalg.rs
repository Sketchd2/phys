//! Dense symmetric linear algebra, without dependencies.
//!
//! The core crate has none and builds for wasm32, so the eigensolver an
//! electronic-structure calculation needs is written here: Householder
//! reduction to tridiagonal form followed by the implicit QL algorithm — the
//! EISPACK routines `tred2` and `tql2`, in the form JAMA made familiar. It is
//! `O(n^3)` and stable, which is all a self-consistent field asks of it.

/// A square matrix, row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    pub n: usize,
    pub a: Vec<f64>,
}

impl Matrix {
    pub fn zeros(n: usize) -> Matrix {
        Matrix { n, a: vec![0.0; n * n] }
    }

    pub fn identity(n: usize) -> Matrix {
        let mut m = Matrix::zeros(n);
        for i in 0..n {
            m.a[i * n + i] = 1.0;
        }
        m
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.a[i * self.n + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.n + j] = v;
    }

    #[inline]
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.n + j] += v;
    }

    pub fn transpose(&self) -> Matrix {
        let n = self.n;
        let mut t = Matrix::zeros(n);
        for i in 0..n {
            for j in 0..n {
                t.a[j * n + i] = self.a[i * n + j];
            }
        }
        t
    }

    /// `self * other`.
    pub fn mul(&self, other: &Matrix) -> Matrix {
        let n = self.n;
        let mut c = Matrix::zeros(n);
        for i in 0..n {
            for k in 0..n {
                let aik = self.a[i * n + k];
                if aik == 0.0 {
                    continue;
                }
                let row = &other.a[k * n..k * n + n];
                let out = &mut c.a[i * n..i * n + n];
                for j in 0..n {
                    out[j] += aik * row[j];
                }
            }
        }
        c
    }

    /// `trace(self * other)`, without forming the product.
    pub fn dot(&self, other: &Matrix) -> f64 {
        // trace(A B) = sum_ij A_ij B_ji
        let n = self.n;
        let mut s = 0.0;
        for i in 0..n {
            for j in 0..n {
                s += self.a[i * n + j] * other.a[j * n + i];
            }
        }
        s
    }

    /// Largest absolute difference from the transpose.
    pub fn asymmetry(&self) -> f64 {
        let n = self.n;
        let mut worst: f64 = 0.0;
        for i in 0..n {
            for j in 0..i {
                worst = worst.max((self.a[i * n + j] - self.a[j * n + i]).abs());
            }
        }
        worst
    }
}

/// Eigenvalues ascending, and eigenvectors as the *columns* of the returned
/// matrix, of a symmetric matrix.
pub fn eigh(m: &Matrix) -> (Vec<f64>, Matrix) {
    let n = m.n;
    if n == 0 {
        return (Vec::new(), Matrix::zeros(0));
    }
    // V holds the matrix and becomes the eigenvectors; d and e the diagonal
    // and off-diagonal of the tridiagonal form.
    let mut v: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| m.get(i, j)).collect()).collect();
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    tred2(&mut v, &mut d, &mut e);
    tql2(&mut v, &mut d, &mut e);
    // Sort ascending, carrying the vectors.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| d[a].total_cmp(&d[b]));
    let values: Vec<f64> = order.iter().map(|&k| d[k]).collect();
    let mut vectors = Matrix::zeros(n);
    for (col, &k) in order.iter().enumerate() {
        for row in 0..n {
            vectors.set(row, col, v[row][k]);
        }
    }
    (values, vectors)
}

fn tred2(v: &mut [Vec<f64>], d: &mut [f64], e: &mut [f64]) {
    let n = d.len();
    for j in 0..n {
        d[j] = v[n - 1][j];
    }
    for i in (1..n).rev() {
        let mut scale = 0.0;
        let mut h = 0.0;
        for k in 0..i {
            scale += d[k].abs();
        }
        if scale == 0.0 {
            e[i] = d[i - 1];
            for j in 0..i {
                d[j] = v[i - 1][j];
                v[i][j] = 0.0;
                v[j][i] = 0.0;
            }
        } else {
            for k in 0..i {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            let mut f = d[i - 1];
            let mut g = h.sqrt();
            if f > 0.0 {
                g = -g;
            }
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for j in 0..i {
                e[j] = 0.0;
            }
            for j in 0..i {
                f = d[j];
                v[j][i] = f;
                g = e[j] + v[j][j] * f;
                for k in (j + 1)..i {
                    g += v[k][j] * d[k];
                    e[k] += v[k][j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for j in 0..i {
                e[j] /= h;
                f += e[j] * d[j];
            }
            let hh = f / (h + h);
            for j in 0..i {
                e[j] -= hh * d[j];
            }
            for j in 0..i {
                f = d[j];
                g = e[j];
                for k in j..i {
                    v[k][j] -= f * e[k] + g * d[k];
                }
                d[j] = v[i - 1][j];
                v[i][j] = 0.0;
            }
        }
        d[i] = h;
    }
    for i in 0..n - 1 {
        v[n - 1][i] = v[i][i];
        v[i][i] = 1.0;
        let h = d[i + 1];
        if h != 0.0 {
            for k in 0..=i {
                d[k] = v[k][i + 1] / h;
            }
            for j in 0..=i {
                let mut g = 0.0;
                for k in 0..=i {
                    g += v[k][i + 1] * v[k][j];
                }
                for k in 0..=i {
                    v[k][j] -= g * d[k];
                }
            }
        }
        for k in 0..=i {
            v[k][i + 1] = 0.0;
        }
    }
    for j in 0..n {
        d[j] = v[n - 1][j];
        v[n - 1][j] = 0.0;
    }
    v[n - 1][n - 1] = 1.0;
    e[0] = 0.0;
}

fn tql2(v: &mut [Vec<f64>], d: &mut [f64], e: &mut [f64]) {
    let n = d.len();
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    let mut f = 0.0;
    let mut tst1: f64 = 0.0;
    let eps = f64::EPSILON;
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n {
            if e[m].abs() <= eps * tst1 {
                break;
            }
            m += 1;
        }
        if m > l {
            loop {
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for i in (l + 2)..n {
                    d[i] -= h;
                }
                f += h;
                p = d[m];
                let mut c = 1.0;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0;
                let mut s2 = 0.0;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for k in 0..n {
                        h = v[k][i + 1];
                        v[k][i + 1] = s * v[k][i] + c * h;
                        v[k][i] = c * v[k][i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= eps * tst1 {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }
}

/// A transformation `X` with `X^T S X = 1`, from the overlap `S` of a basis
/// that may be nearly linearly dependent.
///
/// Canonical orthogonalisation: eigenvectors of `S` scaled by `s^-1/2`, with
/// those whose eigenvalue is below `threshold` *dropped* rather than inverted.
/// An even-tempered basis grown towards completeness becomes nearly dependent
/// by construction, and inverting a 1e-10 eigenvalue amplifies round-off by
/// 1e5. Returns `X` as `n` rows by `m` columns (row-major), and `m`.
pub fn orthogonaliser(s: &Matrix, threshold: f64) -> (Vec<f64>, usize) {
    let n = s.n;
    let (vals, vecs) = eigh(s);
    let keep: Vec<usize> = (0..n).filter(|&k| vals[k] > threshold).collect();
    let m = keep.len();
    let mut x = vec![0.0; n * m];
    for (col, &k) in keep.iter().enumerate() {
        let f = 1.0 / vals[k].sqrt();
        for row in 0..n {
            x[row * m + col] = vecs.get(row, k) * f;
        }
    }
    (x, m)
}

/// `X^T A X` for `X` of `n` rows by `m` columns.
pub fn congruence(x: &[f64], m: usize, a: &Matrix) -> Matrix {
    let n = a.n;
    // T = A X  (n x m)
    let mut t = vec![0.0; n * m];
    for i in 0..n {
        for k in 0..n {
            let aik = a.a[i * n + k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..m {
                t[i * m + j] += aik * x[k * m + j];
            }
        }
    }
    let mut out = Matrix::zeros(m);
    for i in 0..m {
        for k in 0..n {
            let xki = x[k * m + i];
            if xki == 0.0 {
                continue;
            }
            for j in 0..m {
                out.a[i * m + j] += xki * t[k * m + j];
            }
        }
    }
    out
}

/// Solve `F C = S C e` for symmetric `F` and positive-definite `S`, given the
/// orthogonaliser of `S`: eigenvalues ascending and coefficients as the columns
/// of an `n x m` row-major array.
pub fn generalised(f: &Matrix, x: &[f64], m: usize) -> (Vec<f64>, Vec<f64>) {
    let n = f.n;
    let fp = congruence(x, m, f);
    let (e, cp) = eigh(&fp);
    let mut c = vec![0.0; n * m];
    for i in 0..n {
        for k in 0..m {
            let xik = x[i * m + k];
            if xik == 0.0 {
                continue;
            }
            for j in 0..m {
                c[i * m + j] += xik * cp.get(k, j);
            }
        }
    }
    (e, c)
}
