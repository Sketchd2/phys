//! Deterministic, *stateless-addressable* randomness.
//!
//! The engine's core trick — regenerating unobserved detail instead of storing
//! it — only works if regeneration is exactly repeatable. So randomness here is
//! never a global stream you draw from in whatever order the scheduler happens
//! to run. It is a pure function
//!
//! ```text
//!     value = f(world_seed, path_key, epoch, purpose, index)
//! ```
//!
//! Any thread, any frame, any order: the same tuple yields the same bits. That
//! is what makes "coarsen a star, fly away, come back, refine it again" give
//! you back the same star.

use crate::math::{v3, Vec3};

/// SplitMix64. Fast, no state carried between calls, excellent avalanche.
#[inline]
pub const fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Mix two 64-bit words into one. Used to fold path components together.
#[inline]
pub const fn mix2(a: u64, b: u64) -> u64 {
    splitmix64(a ^ splitmix64(b.wrapping_add(0x2545_F491_4F6C_DD1D)))
}

/// What a random draw is *for*. Two different purposes at the same node draw
/// from decorrelated streams, so adding a new physical process to the engine
/// never disturbs the numbers an existing one gets — old scenarios keep
/// replaying identically.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u64)]
pub enum Purpose {
    Positions = 0x01,
    Velocities = 0x02,
    Masses = 0x03,
    Composition = 0x04,
    StellarIMF = 0x05,
    Decay = 0x06,
    Scattering = 0x07,
    QuantumMeasure = 0x08,
    ThermalNoise = 0x09,
    PhotonEmission = 0x0A,
    Structure = 0x0B,
    Spin = 0x0C,
}

/// A reproducible draw sequence. Constructed from an address, not from entropy.
#[derive(Debug, Clone, Copy)]
pub struct Stream {
    base: u64,
    counter: u64,
}

impl Stream {
    /// The canonical constructor: everything that seeds randomness in this
    /// engine goes through here.
    pub fn at(world_seed: u64, path_key: u128, epoch: u32, purpose: Purpose) -> Stream {
        let lo = path_key as u64;
        let hi = (path_key >> 64) as u64;
        let mut b = mix2(world_seed, lo);
        b = mix2(b, hi);
        b = mix2(b, epoch as u64);
        b = mix2(b, purpose as u64);
        Stream { base: b, counter: 0 }
    }

    /// Derive a decorrelated substream.
    ///
    /// Needed wherever the *same* address is drawn from repeatedly over time —
    /// a thermostat, a scattering kernel, anything called once per step. Making
    /// randomness a pure function of the address is the point of this module,
    /// but if the address does not include the step, "deterministic" quietly
    /// becomes "identical every step": the thermostat then applies the same
    /// kick to the same atom forever, which is a constant force, and the system
    /// heats without bound. The bug looks like a thermostat calibration error
    /// and is not one.
    pub fn split(&self, index: u64) -> Stream {
        Stream {
            base: mix2(self.base, index),
            counter: 0,
        }
    }

    /// Direct indexed access — order-independent, so a GPU kernel can compute
    /// draw `i` for particle `i` in any lane order and match the CPU exactly.
    #[inline]
    pub fn nth_u64(&self, i: u64) -> u64 {
        splitmix64(self.base ^ splitmix64(i.wrapping_add(0x1234_5678_9ABC_DEF0)))
    }

    #[inline]
    pub fn nth_f64(&self, i: u64) -> f64 {
        u64_to_unit(self.nth_u64(i))
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let v = self.nth_u64(self.counter);
        self.counter += 1;
        v
    }

    /// Uniform on [0,1).
    #[inline]
    pub fn uniform(&mut self) -> f64 {
        u64_to_unit(self.next_u64())
    }

    #[inline]
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }

    /// Standard normal, via the Box-Muller transform.
    ///
    /// Box-Muller rather than ziggurat: it is branch-free enough to be
    /// bit-identical between the CPU reference and a GPU kernel, which matters
    /// more here than the ~2x throughput a ziggurat would buy.
    pub fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(f64::MIN_POSITIVE);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    pub fn normal3(&mut self) -> Vec3 {
        v3(self.normal(), self.normal(), self.normal())
    }

    /// Isotropic unit vector (Marsaglia).
    pub fn direction(&mut self) -> Vec3 {
        let z = self.range(-1.0, 1.0);
        let phi = self.range(0.0, std::f64::consts::TAU);
        let r = (1.0 - z * z).max(0.0).sqrt();
        v3(r * phi.cos(), r * phi.sin(), z)
    }

    /// Uniform inside the unit ball.
    pub fn in_ball(&mut self) -> Vec3 {
        let d = self.direction();
        let r = self.uniform().cbrt();
        d.scale(r)
    }

    pub fn exponential(&mut self, rate: f64) -> f64 {
        if rate <= 0.0 {
            return f64::INFINITY;
        }
        -(1.0 - self.uniform()).max(f64::MIN_POSITIVE).ln() / rate
    }

    /// Poisson draw (Knuth for small means, normal approximation for large).
    pub fn poisson(&mut self, lambda: f64) -> u64 {
        if lambda <= 0.0 {
            return 0;
        }
        if lambda < 30.0 {
            let l = (-lambda).exp();
            let mut k = 0u64;
            let mut p = 1.0;
            loop {
                p *= self.uniform();
                if p <= l {
                    return k;
                }
                k += 1;
                if k > 1_000_000 {
                    return k;
                }
            }
        } else {
            let g = lambda + lambda.sqrt() * self.normal();
            g.max(0.0).round() as u64
        }
    }

    /// Maxwell-Boltzmann speed for mass `m` at temperature `t`: each Cartesian
    /// component is N(0, sqrt(kT/m)).
    pub fn maxwell(&mut self, mass: f64, temperature: f64) -> Vec3 {
        if mass <= 0.0 || temperature <= 0.0 {
            return Vec3::ZERO;
        }
        let sigma = (crate::units::K_B * temperature / mass).sqrt();
        self.normal3().scale(sigma)
    }

    /// Sample an index from unnormalised weights. Deterministic given the
    /// stream; used for Born-rule outcomes and branching ratios alike.
    pub fn weighted(&mut self, weights: &[f64]) -> usize {
        let total = crate::math::det_sum(weights);
        if !(total > 0.0) {
            return 0;
        }
        let mut u = self.uniform() * total;
        for (i, w) in weights.iter().enumerate() {
            u -= *w;
            if u <= 0.0 {
                return i;
            }
        }
        weights.len() - 1
    }

    /// Power-law sample on [lo, hi] with `dN/dx ∝ x^alpha` (alpha != -1).
    /// Used for the stellar IMF and for turbulent energy cascades.
    pub fn power_law(&mut self, lo: f64, hi: f64, alpha: f64) -> f64 {
        let u = self.uniform();
        if (alpha + 1.0).abs() < 1e-9 {
            lo * (hi / lo).powf(u)
        } else {
            let a1 = alpha + 1.0;
            let lo_a = lo.powf(a1);
            let hi_a = hi.powf(a1);
            (lo_a + u * (hi_a - lo_a)).powf(1.0 / a1)
        }
    }
}

/// 53-bit mantissa fill: uniform on [0,1) with no rounding bias.
#[inline]
pub fn u64_to_unit(x: u64) -> f64 {
    ((x >> 11) as f64) * (1.0 / (1u64 << 53) as f64)
}

/// A one-shot hash → unit float, for when you want a value without a stream.
#[inline]
pub fn hash_unit(a: u64, b: u64) -> f64 {
    u64_to_unit(mix2(a, b))
}
