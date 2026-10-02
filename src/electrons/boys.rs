//! The Boys function, `F_n(T) = integral_0^1 t^(2n) exp(-T t^2) dt`.
//!
//! Every Coulomb integral over Gaussians reduces to it, so it is evaluated to
//! round-off rather than to a tolerance. The highest order is computed
//! directly — by its convergent series for moderate `T`, by the asymptotic form
//! where that is exact to double precision — and the lower orders come from the
//! downward recursion `F_(n-1) = (2T F_n + exp(-T)) / (2n - 1)`, which is stable
//! in that direction.

/// `F_0(T) ... F_nmax(T)` into `out[0..=nmax]`.
pub fn boys(nmax: usize, t: f64, out: &mut [f64]) {
    debug_assert!(out.len() > nmax);
    let et = (-t).exp();
    let top = if t < 1e-14 {
        1.0 / (2 * nmax + 1) as f64
    } else if t > 50.0 {
        // F_n(T) ~ (2n-1)!! / 2^(n+1) * sqrt(pi / T^(2n+1)); the neglected
        // term is of order exp(-T), below 2e-22 here.
        let mut df = 1.0; // (2n-1)!!
        for k in 1..=nmax {
            df *= (2 * k - 1) as f64;
        }
        df / 2f64.powi(nmax as i32 + 1) * (std::f64::consts::PI / t.powi(2 * nmax as i32 + 1)).sqrt()
    } else {
        // exp(-T) * sum_k (2T)^k / ((2n+1)(2n+3)...(2n+2k+1)); every term is
        // positive, so the sum is stopped when one no longer changes it.
        let mut term = 1.0 / (2 * nmax + 1) as f64;
        let mut sum = term;
        let mut k = 1;
        loop {
            term *= 2.0 * t / (2 * nmax + 2 * k + 1) as f64;
            let next = sum + term;
            if next == sum {
                break;
            }
            sum = next;
            k += 1;
            if k > 2000 {
                break;
            }
        }
        et * sum
    };
    out[nmax] = top;
    for n in (1..=nmax).rev() {
        out[n - 1] = (2.0 * t * out[n] + et) / (2 * n - 1) as f64;
    }
}
