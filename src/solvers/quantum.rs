//! The quantum tier — and the reason the whole engine is allowed to work.
//!
//! # Measurement is not a metaphor here
//!
//! The engine's central move is to leave detail unmaterialised until someone
//! looks, then generate it from a seeded distribution. At macroscopic scales
//! that is an engineering trick with a consistency proof attached. At subatomic
//! scales it is *literally the physics*: an unobserved system genuinely does
//! not have a definite position, and a measurement genuinely does sample an
//! outcome from a distribution and leave the system in that outcome.
//!
//! So below the decoherence scale the engine stops pretending to store
//! particles with trajectories and stores what actually exists — occupation
//! amplitudes over states. `measure` samples the Born rule and writes the
//! result to the ledger, and from then on that value is a fact about the world.
//! This is the same "generate on demand, commit on observation" machinery used
//! at every other tier, except that here it is not an approximation to
//! anything. The approximation runs the other way: classical trajectories are
//! the approximation, and `regime` says when they are safe.
//!
//! # The information bound
//!
//! The uncertainty principle puts a hard ceiling on how much detail any region
//! can be asked to have: a phase-space volume `V * p^3` contains at most
//! `V p^3 / h^3` distinguishable states. "Simulate down to the subatomic level"
//! is therefore a *finite* demand per unit volume, and `max_distinguishable_states`
//! computes it. An observer cannot request infinite detail no matter how far
//! they zoom, because there is no infinite detail there to request.

use crate::math::Vec3;
use crate::rng::{Purpose, Stream};
use crate::units::*;

/// Which description is valid here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Regime {
    /// Trajectories are meaningful. Everything above the molecular tier.
    Classical,
    /// Trajectories are meaningful but tunnelling and zero-point motion are
    /// not negligible. Molecular tier, cold condensed matter.
    Semiclassical,
    /// Only distributions are meaningful. Bound electrons, nucleons, anything
    /// below the thermal de Broglie wavelength.
    Quantum,
}

/// Thermal de Broglie wavelength — the length at which a particle's wave nature
/// stops being ignorable.
#[inline]
pub fn thermal_wavelength(mass: f64, temperature: f64) -> f64 {
    if mass <= 0.0 || temperature <= 0.0 {
        return f64::INFINITY;
    }
    H_PLANCK / (2.0 * std::f64::consts::PI * mass * K_B * temperature).sqrt()
}

#[inline]
pub fn de_broglie(mass: f64, speed: f64) -> f64 {
    if mass <= 0.0 || speed <= 0.0 {
        return f64::INFINITY;
    }
    H_PLANCK / (mass * speed)
}

/// Decide the regime from the particle's mass, temperature and the separation
/// being resolved.
pub fn regime(mass: f64, temperature: f64, separation: f64) -> Regime {
    let lambda = thermal_wavelength(mass, temperature);
    if separation > 100.0 * lambda {
        Regime::Classical
    } else if separation > lambda {
        Regime::Semiclassical
    } else {
        Regime::Quantum
    }
}

/// Maximum number of distinguishable states in a phase-space volume.
///
/// The engine consults this before honouring a refinement request: asking for
/// more detail than this is asking for information that does not exist, and the
/// correct response is to stop refining, not to invent it.
pub fn max_distinguishable_states(volume: f64, momentum_spread: f64) -> f64 {
    if volume <= 0.0 || momentum_spread <= 0.0 {
        return 0.0;
    }
    volume * momentum_spread.powi(3) / H_PLANCK.powi(3)
}

/// The smallest length scale at which it is meaningful to ask for a particle's
/// position, given how well its momentum is known.
#[inline]
pub fn position_resolution_limit(momentum_spread: f64) -> f64 {
    if momentum_spread <= 0.0 {
        f64::INFINITY
    } else {
        H_BAR / (2.0 * momentum_spread)
    }
}

/// Enforce `Δx Δp >= ħ/2` on a sampled pair. Returns the corrected momentum
/// spread. The engine calls this on every materialisation at the atomic tier
/// and below, so it can never hand an observer a state that violates it.
pub fn enforce_uncertainty(position_spread: f64, momentum_spread: f64) -> f64 {
    if position_spread <= 0.0 {
        return momentum_spread;
    }
    let floor = H_BAR / (2.0 * position_spread);
    momentum_spread.max(floor)
}

/// A discrete quantum state with an occupation amplitude.
#[derive(Debug, Clone, Copy)]
pub struct Level {
    /// Energy above the ground state, J.
    pub energy: f64,
    /// Degeneracy.
    pub g: f64,
    /// Occupation probability (not amplitude — the engine tracks populations,
    /// which is what survives decoherence and what an experiment measures).
    pub population: f64,
}

/// A statistical description of a subatomic system: populations over levels,
/// never trajectories.
#[derive(Debug, Clone)]
pub struct Ensemble {
    pub levels: Vec<Level>,
    pub temperature: f64,
}

impl Ensemble {
    /// Hydrogenic level structure: `E_n = -13.6 Z^2 / n^2` eV. Enough to give
    /// the right spectral lines, which is what an observer actually sees.
    pub fn hydrogenic(z: f64, n_max: usize, temperature: f64) -> Ensemble {
        let ground = -13.605_693_122_994 * z * z * EV;
        let mut levels = Vec::with_capacity(n_max);
        for n in 1..=n_max {
            let e = -13.605_693_122_994 * z * z / (n * n) as f64 * EV;
            levels.push(Level {
                energy: e - ground,
                g: 2.0 * (n * n) as f64,
                population: 0.0,
            });
        }
        let mut e = Ensemble {
            levels,
            temperature,
        };
        e.thermalise();
        e
    }

    /// Populate by the Boltzmann distribution at the ensemble's temperature.
    pub fn thermalise(&mut self) {
        let kt = K_B * self.temperature.max(1e-6);
        let mut total = 0.0;
        for l in self.levels.iter_mut() {
            l.population = l.g * (-l.energy / kt).exp();
            total += l.population;
        }
        if total > 0.0 {
            for l in self.levels.iter_mut() {
                l.population /= total;
            }
        }
    }

    pub fn mean_energy(&self) -> f64 {
        self.levels.iter().map(|l| l.population * l.energy).sum()
    }

    /// Sample a definite level — a measurement. The Born rule *is* the
    /// weighted draw; there is nothing extra to model.
    ///
    /// Deterministic in the stream, so a replay of the same session measures
    /// the same outcome. That is not a violation of quantum indeterminacy: the
    /// outcome is unpredictable to anyone inside the simulation, which is the
    /// only place the question can be asked.
    pub fn measure_level(&self, stream: &mut Stream) -> usize {
        let weights: Vec<f64> = self.levels.iter().map(|l| l.population).collect();
        stream.weighted(&weights)
    }

    /// Collapse onto a measured level. Subsequent measurements agree with this
    /// one until something re-thermalises the system — which is exactly the
    /// behaviour the ledger needs from a "committed" fact.
    pub fn collapse_to(&mut self, index: usize) {
        for (i, l) in self.levels.iter_mut().enumerate() {
            l.population = if i == index { 1.0 } else { 0.0 };
        }
    }

    /// Photon energy released by a transition between two levels.
    pub fn transition_energy(&self, from: usize, to: usize) -> f64 {
        self.levels[from].energy - self.levels[to].energy
    }
}

/// WKB tunnelling probability through a square barrier.
///
/// Needed at the nuclear tier — fusion in a stellar core is tunnelling, and
/// without it the Sun does not shine at the temperature it actually has.
pub fn tunnelling_probability(mass: f64, barrier: f64, energy: f64, width: f64) -> f64 {
    if energy >= barrier || mass <= 0.0 || width <= 0.0 {
        return 1.0;
    }
    let kappa = (2.0 * mass * (barrier - energy)).sqrt() / H_BAR;
    (-2.0 * kappa * width).exp().clamp(0.0, 1.0)
}

/// Zero-point energy of a harmonic mode. Matters for the heat capacity of cold
/// solids, and for why a crystal does not simply stop at absolute zero.
#[inline]
pub fn zero_point_energy(omega: f64) -> f64 {
    0.5 * H_BAR * omega
}

/// Sample a photon energy from a blackbody at temperature `t`.
///
/// Rejection sampling against the Planck distribution in the dimensionless
/// variable `x = hν/kT`, whose mode is at 2.821. This is what turns a node's
/// temperature into something an observer's spectrograph can see.
pub fn sample_blackbody_photon(temperature: f64, stream: &mut Stream) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    let peak = 1.421; // maximum of x^2/(e^x - 1)
    for _ in 0..64 {
        let x = stream.range(1e-6, 20.0);
        let p = x * x / ((x.exp() - 1.0).max(1e-300));
        if stream.uniform() * peak < p {
            return x * K_B * temperature;
        }
    }
    2.821 * K_B * temperature
}

/// Compton scattering: the wavelength shift for a photon deflected by `theta`.
pub fn compton_shift(theta: f64) -> f64 {
    LAMBDA_COMPTON_E * (1.0 - theta.cos())
}

/// A quantum measurement performed by an instrument on a node.
///
/// The pair `(value, disturbance)` is the honest output of any real
/// measurement: you learn something, and you change the thing you learned it
/// from. The engine applies `disturbance` as an actual energy injection, which
/// is why repeatedly measuring a small system in this simulation heats it — as
/// it would in a laboratory.
#[derive(Debug, Clone, Copy)]
pub struct Measurement {
    pub value: f64,
    pub uncertainty: f64,
    pub disturbance: f64,
}

/// Measure a position to a requested precision. The finer you look, the more
/// momentum you deposit — `Δp >= ħ/2Δx`, deposited as real kinetic energy.
pub fn measure_position(
    mass: f64,
    true_pos: Vec3,
    axis: Vec3,
    precision: f64,
    stream: &mut Stream,
) -> Measurement {
    let sigma = precision.max(1e-30);
    let value = true_pos.dot(axis.unit()) + sigma * stream.normal();
    let dp = H_BAR / (2.0 * sigma);
    let disturbance = if mass > 0.0 { dp * dp / (2.0 * mass) } else { 0.0 };
    Measurement {
        value,
        uncertainty: sigma,
        disturbance,
    }
}

/// Spontaneous emission lifetime of an excited hydrogenic level (rough scaling).
pub fn radiative_lifetime(n: usize, z: f64) -> f64 {
    if n <= 1 {
        return f64::INFINITY;
    }
    1.6e-9 * (n as f64).powi(5) / z.powi(4)
}

/// Sample whether a transition happened in time `dt`.
pub fn emitted(lifetime: f64, dt: f64, stream: &mut Stream) -> bool {
    if !lifetime.is_finite() || lifetime <= 0.0 {
        return false;
    }
    stream.uniform() < 1.0 - (-dt / lifetime).exp()
}

/// A convenience for the observation system: how much detail can be extracted
/// from a region before the uncertainty principle says "there is no more".
pub fn detail_ceiling(volume: f64, temperature: f64, mass: f64) -> f64 {
    let p = (3.0 * mass * K_B * temperature).max(0.0).sqrt();
    max_distinguishable_states(volume, p)
}

/// Stream for a measurement on a given node — routed through the standard
/// address so that measurements replay.
pub fn measurement_stream(world_seed: u64, path_key: u128, epoch: u32, index: u64) -> Stream {
    let s = Stream::at(world_seed, path_key, epoch, Purpose::QuantumMeasure);
    let _ = s.nth_u64(index);
    s
}
