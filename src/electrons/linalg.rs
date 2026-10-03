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

/// `out[i][j] = sum_p a[i][p] b[j][p]`: `a` is `m x kd` and `b` is `n x kd`,
/// both rows contiguous, `out` is `m x n` and is overwritten.
///
/// Every dense product in a grid integral takes this shape once its operands
/// are laid out by function, and written as a plain triple loop it was most of
/// the exchange-correlation build. Here `b` is packed four rows at a time so
/// that four columns of the result load together, and the result is built
/// four by four so each value loaded is used four times: 10.6 GFLOP/s on one
/// 2.8 GHz core with two-wide vectors, which is that width's peak. Where the
/// processor has four-wide vectors (AVX) the same code is compiled for them
/// and chosen at run time. Each element is summed over `p` in order with a
/// separate multiply and add (no fused multiply-add), so the result is the
/// same to the last bit on either path and equal to the plain loop.
pub fn product_nt(a: &[f64], b: &[f64], m: usize, n: usize, kd: usize, out: &mut [f64]) {
    assert!(a.len() >= m * kd && b.len() >= n * kd && out.len() >= m * n);
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx") {
            // Safety: the processor has just been found to support AVX.
            unsafe { product_nt_avx(a, b, m, n, kd, out) };
            return;
        }
    }
    product_nt_body(a, b, m, n, kd, out);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx")]
unsafe fn product_nt_avx(a: &[f64], b: &[f64], m: usize, n: usize, kd: usize, out: &mut [f64]) {
    product_nt_body(a, b, m, n, kd, out);
}

#[inline(always)]
fn product_nt_body(a: &[f64], b: &[f64], m: usize, n: usize, kd: usize, out: &mut [f64]) {
    let (m4, n4) = (m - m % 4, n - n % 4);
    // b's rows four at a time, interleaved: panel[jb][p][c] = b[4 jb + c][p].
    let mut panel = vec![0.0; n4 * kd];
    for jb in 0..n4 / 4 {
        let dst = &mut panel[jb * 4 * kd..(jb + 1) * 4 * kd];
        for c in 0..4 {
            let row = &b[(4 * jb + c) * kd..(4 * jb + c + 1) * kd];
            for (p, &y) in row.iter().enumerate() {
                dst[4 * p + c] = y;
            }
        }
    }
    for i in (0..m4).step_by(4) {
        let a0 = &a[i * kd..(i + 1) * kd];
        let a1 = &a[(i + 1) * kd..(i + 2) * kd];
        let a2 = &a[(i + 2) * kd..(i + 3) * kd];
        let a3 = &a[(i + 3) * kd..(i + 4) * kd];
        for jb in 0..n4 / 4 {
            let pan = &panel[jb * 4 * kd..(jb + 1) * 4 * kd];
            let mut s = [[0.0f64; 4]; 4];
            for p in 0..kd {
                let y: &[f64; 4] = pan[4 * p..4 * p + 4].try_into().expect("four");
                let x = [a0[p], a1[p], a2[p], a3[p]];
                for r in 0..4 {
                    for c in 0..4 {
                        s[r][c] += x[r] * y[c];
                    }
                }
            }
            let j = 4 * jb;
            for (r, row) in s.iter().enumerate() {
                out[(i + r) * n + j..(i + r) * n + j + 4].copy_from_slice(row);
            }
        }
        for j in n4..n {
            let bj = &b[j * kd..(j + 1) * kd];
            for (r, ar) in [a0, a1, a2, a3].iter().enumerate() {
                let mut t = 0.0;
                for p in 0..kd {
                    t += ar[p] * bj[p];
                }
                out[(i + r) * n + j] = t;
            }
        }
    }
    for i in m4..m {
        let ai = &a[i * kd..(i + 1) * kd];
        for j in 0..n {
            let bj = &b[j * kd..(j + 1) * kd];
            let mut t = 0.0;
            for p in 0..kd {
                t += ai[p] * bj[p];
            }
            out[i * n + j] = t;
        }
    }
}

/// `m x n` row-major into `n x m` row-major.
pub fn transpose(a: &[f64], m: usize, n: usize, out: &mut [f64]) {
    for i in 0..m {
        for j in 0..n {
            out[j * m + i] = a[i * n + j];
        }
    }
}

/// A step of a pivoted Cholesky runs on one thread unless it has at least this
/// much work in it (rows left times columns so far): starting the threads for
/// a step costs about 0.7 ms, and a step smaller than this is done sooner on
/// one. Measured threading every step against none, on near-full-rank metrics:
/// 256 rows 0.175 s against 0.003, 2048 rows 1.53 against 1.37, 4096 rows 4.9
/// against 10.9 — the early steps of any factorisation are small, and the late
/// steps of a large one are not.
const PARALLEL_CHOLESKY_WORK: usize = 2_000_000;

/// Pivoted Cholesky of a symmetric positive semi-definite matrix: at each step
/// the row with the largest remaining diagonal is taken, and the factorisation
/// stops when that falls to `cutoff` times the largest diagonal there was.
///
/// Returns the rows taken, in the order taken, and the factor `L` over them
/// (`k x k`, row-major, lower triangular in that order), with
/// `v[kept][kept] = L L^T`. The rows not taken are the ones the others already
/// represent to within the cutoff — the near-dependent directions of a dense
/// set — so the factorisation selects a well-conditioned subset rather than
/// inverting through the dependence.
pub fn pivoted_cholesky(v: &Matrix, cutoff: f64) -> (Vec<usize>, Vec<f64>) {
    let n = v.n;
    let mut d: Vec<f64> = (0..n).map(|i| v.a[i * n + i]).collect();
    let top = d.iter().cloned().fold(0.0f64, f64::max);
    let mut taken = vec![false; n];
    let mut kept = Vec::new();
    // Columns of the factor over every row, row-major by row: l[i][k].
    let mut l: Vec<Vec<f64>> = vec![Vec::new(); n];
    while kept.len() < n {
        let (mut p, mut best) = (usize::MAX, -1.0);
        for i in 0..n {
            if !taken[i] && d[i] > best {
                best = d[i];
                p = i;
            }
        }
        if best <= cutoff * top {
            break;
        }
        let lkk = best.sqrt();
        taken[p] = true;
        let k = kept.len();
        kept.push(p);
        let lp = l[p].clone();
        // Each row's update reads only its own row, the pivot's and the
        // metric, so a large step's rows are split across threads; every
        // element is the same expression either way, so the factor is
        // bit-identical. Serial, a 5566-function metric (water's grown basis
        // at round 9) took 28 s of a 305 s solve.
        let update = |base: usize, rows: &mut [Vec<f64>], diag: &mut [f64]| {
            for (o, (li, di)) in rows.iter_mut().zip(diag.iter_mut()).enumerate() {
                let i = base + o;
                if taken[i] && i != p {
                    continue;
                }
                let x = if i == p {
                    lkk
                } else {
                    let s: f64 = li.iter().zip(&lp).map(|(a, b)| a * b).sum();
                    (v.a[i * n + p] - s) / lkk
                };
                li.push(x);
                if i != p {
                    *di -= x * x;
                }
            }
        };
        #[cfg(not(target_arch = "wasm32"))]
        let workers = if (n - k) * k >= PARALLEL_CHOLESKY_WORK { std::thread::available_parallelism().map(|w| w.get()).unwrap_or(1) } else { 1 };
        #[cfg(target_arch = "wasm32")]
        let workers = 1;
        if workers <= 1 {
            update(0, &mut l, &mut d);
        } else {
            #[cfg(not(target_arch = "wasm32"))]
            {
                let chunk = n.div_ceil(workers);
                let update = &update;
                std::thread::scope(|scope| {
                    for (c, (rows, diag)) in l.chunks_mut(chunk).zip(d.chunks_mut(chunk)).enumerate() {
                        scope.spawn(move || update(c * chunk, rows, diag));
                    }
                });
            }
        }
        debug_assert_eq!(l[p].len(), k + 1);
    }
    let k = kept.len();
    let mut factor = vec![0.0; k * k];
    for (r, &i) in kept.iter().enumerate() {
        factor[r * k..r * k + l[i].len().min(k)].copy_from_slice(&l[i][..l[i].len().min(k)]);
    }
    (kept, factor)
}

/// Solve `L L^T x = b` for a `k x k` lower-triangular `L`, row-major.
pub fn cholesky_solve(l: &[f64], k: usize, b: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; k];
    for i in 0..k {
        let s: f64 = (0..i).map(|j| l[i * k + j] * y[j]).sum();
        y[i] = (b[i] - s) / l[i * k + i];
    }
    let mut x = vec![0.0; k];
    for i in (0..k).rev() {
        let s: f64 = (i + 1..k).map(|j| l[j * k + i] * x[j]).sum();
        x[i] = (y[i] - s) / l[i * k + i];
    }
    x
}
