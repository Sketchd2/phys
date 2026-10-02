//! Gaussian shells: the functions electrons are expanded in.
//!
//! A shell is every Cartesian Gaussian `x^i y^j z^k exp(-a r^2)` with
//! `i + j + k = l` about one centre, sharing exponents and contraction
//! coefficients. Cartesian rather than spherical: a Cartesian `d` shell carries
//! an `s`-like combination along with its five `d` functions, which costs one
//! extra function and nothing in accuracy (a larger space can only lower a
//! variational energy). Atomic units throughout — bohr and hartree.

/// One shell.
#[derive(Debug, Clone, PartialEq)]
pub struct Shell {
    /// Bohr.
    pub centre: [f64; 3],
    /// Angular momentum.
    pub l: usize,
    pub exponents: Vec<f64>,
    /// Coefficients of each primitive *including* its normalisation for the
    /// `(l, 0, 0)` component; the other components differ only by
    /// [`component_scale`].
    pub coefficients: Vec<f64>,
}

/// `(2n - 1)!!`, with `(-1)!! = 1`.
pub fn double_factorial_odd(n: usize) -> f64 {
    let mut r = 1.0;
    let mut k = 2 * n as i64 - 1;
    while k > 1 {
        r *= k as f64;
        k -= 2;
    }
    r
}

/// Normalisation of `x^l exp(-a r^2)`.
pub fn primitive_norm(a: f64, l: usize) -> f64 {
    (2.0 * a / std::f64::consts::PI).powf(0.75) * (4.0 * a).powf(l as f64 / 2.0)
        / double_factorial_odd(l).sqrt()
}

/// The factor taking the `(l,0,0)` normalisation to the `(i,j,k)` one.
pub fn component_scale(l: usize, c: [usize; 3]) -> f64 {
    (double_factorial_odd(l)
        / (double_factorial_odd(c[0]) * double_factorial_odd(c[1]) * double_factorial_odd(c[2])))
    .sqrt()
}

/// The Cartesian components of a shell, in a fixed order.
pub fn components(l: usize) -> Vec<[usize; 3]> {
    let mut out = Vec::with_capacity((l + 1) * (l + 2) / 2);
    for i in (0..=l).rev() {
        for j in (0..=(l - i)).rev() {
            out.push([i, j, l - i - j]);
        }
    }
    out
}

impl Shell {
    /// A shell of one normalised primitive.
    pub fn primitive(centre: [f64; 3], l: usize, exponent: f64) -> Shell {
        Shell { centre, l, exponents: vec![exponent], coefficients: vec![primitive_norm(exponent, l)] }
    }

    /// A contracted shell, normalised as a whole.
    pub fn contracted(centre: [f64; 3], l: usize, exponents: Vec<f64>, coefficients: Vec<f64>) -> Shell {
        let raw: Vec<f64> = exponents.iter().zip(&coefficients).map(|(a, c)| c * primitive_norm(*a, l)).collect();
        // <phi|phi> for the (l,0,0) component of the contraction.
        let mut s = 0.0;
        for (a, ca) in exponents.iter().zip(&raw) {
            for (b, cb) in exponents.iter().zip(&raw) {
                let p = a + b;
                s += ca * cb * double_factorial_odd(l) / (2.0 * p).powi(l as i32) * (std::f64::consts::PI / p).powf(1.5);
            }
        }
        let f = 1.0 / s.sqrt();
        Shell { centre, l, exponents, coefficients: raw.iter().map(|c| c * f).collect() }
    }

    /// The constant function 1, as a shell: what turns a four-centre Coulomb
    /// integral into a three- or two-centre one.
    pub fn unit() -> Shell {
        Shell { centre: [0.0; 3], l: 0, exponents: vec![0.0], coefficients: vec![1.0] }
    }

    pub fn size(&self) -> usize {
        (self.l + 1) * (self.l + 2) / 2
    }
}

/// A basis: shells, and where each starts in the list of functions.
#[derive(Debug, Clone, Default)]
pub struct Basis {
    pub shells: Vec<Shell>,
    pub offsets: Vec<usize>,
    pub size: usize,
}

impl Basis {
    pub fn new(shells: Vec<Shell>) -> Basis {
        let mut offsets = Vec::with_capacity(shells.len());
        let mut size = 0;
        for s in &shells {
            offsets.push(size);
            size += s.size();
        }
        Basis { shells, offsets, size }
    }

    pub fn max_l(&self) -> usize {
        self.shells.iter().map(|s| s.l).max().unwrap_or(0)
    }
}
