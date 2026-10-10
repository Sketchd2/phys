// Shared by the GPU kernels that need the reference's special functions in
// double precision: WGSL has no f64 exp, so `gexp` is built here; `gerf` is the
// reference's own series. Included by gen_induction.py in each template at
// `@COMMON@`; edit it here.

const SQRT_PI = 1.7724538509055160273lf;
const LOG2E = 1.4426950408889634074lf;
const LN2_HI = 0.693147180369123816490lf;
const LN2_LO = 1.90821492927058770002e-10lf;
const INV_ODD = array<f64, 130>(@INV_ODD@);

// exp(x) for x <= 0 to full double precision: x = n ln2 + r, |r| <= ln2/2, a
// degree-14 polynomial for e^r, and 2^n by exact repeated squaring.
fn gexp(x: f64) -> f64 {
    if (x < -700.0lf) {
        return 0.0lf;
    }
    let n = floor(x * LOG2E + 0.5lf);
    let r = (x - n * LN2_HI) - n * LN2_LO;
    @HORNER@
    var m = i32(n);
    var base = 2.0lf;
    if (m < 0) {
        base = 0.5lf;
        m = -m;
    }
    var scale = 1.0lf;
    for (var b = 0; b < 11; b = b + 1) {
        if ((m & 1) != 0) {
            scale = scale * base;
        }
        base = base * base;
        m = m >> 1u;
    }
    return poly * scale;
}

// exp(-x^2), taken as 0 beyond x^2 = 40 (below 5e-18 of one, under the
// reference's own rounding at the places it is used).
fn gauss(x: f64) -> f64 {
    let x2 = x * x;
    if (x2 > 40.0lf) {
        return 0.0lf;
    }
    return gexp(-x2);
}

// The reference's erf: the series for |x| < 5.5, the sign beyond.
fn gerf(x: f64) -> f64 {
    let ax = abs(x);
    if (ax == 0.0lf) {
        return 0.0lf;
    }
    if (ax >= 5.5lf) {
        if (x < 0.0lf) {
            return -1.0lf;
        }
        return 1.0lf;
    }
    let x2 = ax * ax;
    let t2 = 2.0lf * x2;
    var term = ax;
    var sum = ax;
    for (var n = 1u; n < 129u; n = n + 1u) {
        term = term * t2 * INV_ODD[n];
        sum = sum + term;
        if (term <= 1.0e-17lf * sum) {
            break;
        }
    }
    var v = 2.0lf / SQRT_PI * gexp(-x2) * sum;
    if (x < 0.0lf) {
        v = -v;
    }
    return v;
}

