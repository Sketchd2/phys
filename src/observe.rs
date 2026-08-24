//! Observation, measurement, and the ledger.
//!
//! # Observation is the engine's input, not its output
//!
//! In a conventional simulation you compute the world and then look at it. Here
//! the looking comes first: an observation is a *demand* for a particular
//! resolution over a particular solid angle, and it is that demand — nothing
//! else — that decides which parts of the world get materialised. Point a
//! telescope at a molecular cloud and the cloud acquires structure; look away
//! and the structure is released.
//!
//! # Two rules keep that honest
//!
//! **Retardation.** You never see a node's current state, only its state on
//! your past light cone. A star 8 kpc away is shown as it was 26,000 years ago.
//! This is not a stylistic choice — it is what makes the "only simulate what is
//! observed" strategy self-consistent, because the state you need is one the
//! engine has already computed and can no longer be asked to revise.
//!
//! **The ledger.** A measurement is a commitment. Once an outcome has been
//! reported to an observer it is a fact about the world, recorded, and every
//! later query returns it. Without this the deception is trivially detectable:
//! measure a decay time twice and get two different answers.

use crate::causal::RetardedView;
use crate::coords::{angular_size, doppler};
use crate::ids::{NodeIdx, PathKey};
use crate::math::Vec3;
use crate::rng::{Purpose, Stream};
use crate::units::*;
use std::collections::HashMap;

/// An observer: a position, a velocity, and — crucially — a *resolution*.
///
/// The resolution is what the engine budgets against. An observer with a wide
/// field and coarse resolution is nearly free; one with a 1 fm aperture is
/// enormously expensive but only over a nanometre-wide field. The engine will
/// serve either, and the cost model in `budget.rs` is what stops a user from
/// asking for both at once.
#[derive(Debug, Clone, Copy)]
pub struct Observer {
    /// Node whose frame the observer sits in.
    pub anchor: NodeIdx,
    pub offset: Vec3,
    pub velocity: Vec3,
    /// Direction of view.
    pub look: Vec3,
    /// Half-angle of the field of view, radians.
    pub field: f64,
    /// Smallest angle the observer can resolve, radians. Together with the
    /// distance to a target this fixes the tier that must be materialised.
    pub angular_resolution: f64,
    /// How long the observer integrates. Longer means fainter things become
    /// visible and fast things blur — both modelled.
    pub integration_time: f64,
    /// How far into the future the engine promises consistency for this
    /// observer. Sets the causal horizon.
    pub horizon: f64,
    /// Relative importance when the frame budget cannot serve everyone.
    pub priority: f64,
}

impl Default for Observer {
    fn default() -> Self {
        Observer {
            anchor: NodeIdx(0),
            offset: Vec3::ZERO,
            velocity: Vec3::ZERO,
            look: crate::math::v3(0.0, 0.0, 1.0),
            field: 0.5,
            angular_resolution: 1e-4,
            integration_time: 1.0,
            horizon: 1.0,
            priority: 1.0,
        }
    }
}

impl Observer {
    /// The linear resolution this observer demands at distance `d`.
    #[inline]
    pub fn linear_resolution(&self, d: f64) -> f64 {
        (d * self.angular_resolution).max(1e-18)
    }

    /// The tier that must be materialised to satisfy this observer at distance
    /// `d`. This single line is the entire level-of-detail policy.
    pub fn required_tier(&self, d: f64) -> Tier {
        Tier::for_resolution(self.linear_resolution(d))
    }

    /// Is `dir` inside the field of view?
    pub fn sees(&self, dir: Vec3) -> bool {
        let l = self.look.unit();
        let d = dir.unit();
        d.dot(l) >= self.field.cos()
    }
}

/// What an instrument reports. Every variant carries its own uncertainty,
/// because an instrument that reports a number without one is lying.
#[derive(Debug, Clone)]
pub enum Reading {
    /// Bolometric flux, W/m^2, plus the light travel time it arrived after.
    Flux { value: f64, uncertainty: f64, delay: f64 },
    /// Spectral bins: (wavelength m, specific flux).
    Spectrum { bins: Vec<(f64, f64)>, redshift: f64 },
    Temperature { kelvin: f64, uncertainty: f64 },
    /// Counts in an integration window, with Poisson noise already applied.
    ParticleCount { counts: u64, expected: f64 },
    /// Position along an axis, with the disturbance the measurement caused.
    Position { value: f64, uncertainty: f64, disturbance: f64 },
    Composition { fractions: [f64; NSPECIES] },
    /// Bulk properties an observer could infer without resolving structure.
    Bulk { mass: f64, radius: f64, velocity: Vec3 },
    /// The instrument could not resolve the target at all.
    Unresolved { angular_size: f64, needed: f64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Instrument {
    /// Bolometric imaging.
    Imager,
    Spectrometer,
    /// Direct contact thermometry — only valid if the observer is inside.
    Thermometer,
    /// Counts particles crossing a detector area.
    ParticleDetector,
    /// Position measurement, subject to the uncertainty principle.
    Interferometer,
    MassSpectrometer,
}

/// A fact the world has committed to.
///
/// Once an outcome is in here it is no longer derivable from the procedural
/// seed — it *is* the world. This is the only structure whose size grows with
/// what users do rather than with the size of the universe, and keeping it
/// small is why the engine can be interactive at all.
#[derive(Debug, Clone, Copy)]
pub struct Fact {
    pub value: f64,
    pub time: f64,
    pub quantity: Quantity,
    /// Which observation produced it, for provenance and replay.
    pub sequence: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Quantity {
    DecayTime,
    QuantumLevel,
    Position,
    Momentum,
    Spin,
    PhotonCount,
    Temperature,
}

/// The commitment store.
#[derive(Default)]
pub struct Ledger {
    facts: HashMap<(PathKey, Quantity), Fact>,
    pub sequence: u64,
    pub queries: u64,
    pub commits: u64,
}

impl Ledger {
    pub fn new() -> Ledger {
        Ledger::default()
    }

    /// Look up a committed value, or sample and commit one.
    ///
    /// This is the operation the whole design rests on. The first call draws
    /// from the physically correct distribution; every later call returns the
    /// same number. An observer inside the simulation cannot distinguish this
    /// from a world in which the value was there all along — which is the
    /// standard the engine holds itself to.
    pub fn get_or_sample<F: FnOnce(&mut Stream) -> f64>(
        &mut self,
        key: PathKey,
        quantity: Quantity,
        time: f64,
        world_seed: u64,
        epoch: u32,
        sample: F,
    ) -> Fact {
        self.queries += 1;
        if let Some(f) = self.facts.get(&(key, quantity)) {
            return *f;
        }
        let mut stream = Stream::at(world_seed, key.0, epoch, Purpose::QuantumMeasure);
        let value = sample(&mut stream);
        self.sequence += 1;
        self.commits += 1;
        let fact = Fact {
            value,
            time,
            quantity,
            sequence: self.sequence,
        };
        self.facts.insert((key, quantity), fact);
        fact
    }

    pub fn peek(&self, key: PathKey, quantity: Quantity) -> Option<Fact> {
        self.facts.get(&(key, quantity)).copied()
    }

    /// Record a value produced by an interaction rather than a measurement.
    pub fn commit(&mut self, key: PathKey, quantity: Quantity, value: f64, time: f64) {
        self.sequence += 1;
        self.commits += 1;
        self.facts.insert(
            (key, quantity),
            Fact {
                value,
                time,
                quantity,
                sequence: self.sequence,
            },
        );
    }

    pub fn len(&self) -> usize {
        self.facts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.facts.is_empty()
    }

    /// Bytes held. Reported in the debug overlay next to the procedural world
    /// size, because the ratio between them is the engine's headline number.
    pub fn bytes(&self) -> usize {
        self.facts.len() * (std::mem::size_of::<Fact>() + std::mem::size_of::<(PathKey, Quantity)>())
    }

    /// Release facts that have become indistinguishable from what the
    /// procedural generator would produce anyway — a measured value that agrees
    /// with the distribution's median to within its own uncertainty carries no
    /// information and need not be stored.
    ///
    /// This is the pressure valve that keeps the ledger from growing without
    /// bound in a long session. It is safe precisely because "indistinguishable"
    /// is checked against the observer's own resolution.
    pub fn compact<F: Fn(&Fact) -> bool>(&mut self, redundant: F) -> usize {
        let before = self.facts.len();
        self.facts.retain(|_, f| !redundant(f));
        before - self.facts.len()
    }
}

/// A resolved observation of one node.
#[derive(Debug, Clone)]
pub struct Sighting {
    pub node: NodeIdx,
    pub key: PathKey,
    /// Where and when it was, on the observer's past light cone.
    pub view: RetardedView,
    pub angular_size: f64,
    /// Whether the observer's resolution is fine enough to see structure.
    pub resolved: bool,
    /// Tier the engine must materialise to satisfy this sighting.
    pub required_tier: Tier,
    pub flux: f64,
    pub doppler: f64,
    pub reading: Reading,
}

/// Bolometric flux at the observer from a source of luminosity `l` at distance
/// `d`, including the relativistic beaming factor.
pub fn flux(l: f64, d: f64, doppler_factor: f64) -> f64 {
    if d <= 0.0 {
        return f64::INFINITY;
    }
    // The fourth power is not a typo: one factor from photon energy, one from
    // arrival rate, two from solid-angle beaming.
    l * doppler_factor.powi(4) / (4.0 * std::f64::consts::PI * d * d)
}

/// Apply an instrument to a retarded view of a node.
pub fn read(
    instrument: Instrument,
    obs: &Observer,
    view: &RetardedView,
    agg: &crate::state::Aggregate,
    stream: &mut Stream,
) -> Reading {
    let d = view.distance.max(1e-30);
    let sep_dir = (view.snapshot.offset - obs.offset).unit();
    let dop = doppler(view.snapshot.velocity - obs.velocity, sep_dir);
    let theta = angular_size(agg.radius, d);

    match instrument {
        Instrument::Imager => {
            let f = flux(agg.luminosity, d, dop);
            if theta < obs.angular_resolution {
                // Unresolved: the observer gets a point source and a bound on
                // its size, which is exactly what a real telescope gets.
                return Reading::Unresolved {
                    angular_size: theta,
                    needed: obs.angular_resolution,
                };
            }
            // Photon shot noise: the fundamental limit on any flux measurement.
            let photon_e = 2.7 * K_B * agg.temperature.max(2.725);
            let n = if photon_e > 0.0 {
                f * obs.integration_time / photon_e
            } else {
                0.0
            };
            let rel = if n > 0.0 { 1.0 / n.sqrt() } else { 1.0 };
            Reading::Flux {
                value: f,
                uncertainty: f * rel,
                delay: view.delay,
            }
        }
        Instrument::Spectrometer => {
            let mut bins = Vec::with_capacity(32);
            let t = agg.temperature.max(2.725);
            for i in 0..32 {
                let lambda = 1e-8 * (10.0f64).powf(i as f64 / 8.0);
                // Planck function, Doppler-shifted into the observer's frame.
                let l_shift = lambda / dop;
                let x = H_PLANCK * C / (l_shift * K_B * t);
                let b = if x < 500.0 && x > 1e-6 {
                    2.0 * H_PLANCK * C * C / l_shift.powi(5) / ((x.exp() - 1.0).max(1e-300))
                } else {
                    0.0
                };
                bins.push((lambda, b * dop.powi(3)));
            }
            Reading::Spectrum {
                bins,
                redshift: 1.0 / dop - 1.0,
            }
        }
        Instrument::Thermometer => Reading::Temperature {
            kelvin: agg.temperature,
            // A thermometer in contact with N particles cannot do better than
            // the thermodynamic fluctuation limit.
            uncertainty: agg.temperature / agg.particle_count().max(1.0).sqrt(),
        },
        Instrument::ParticleDetector => {
            let f = flux(agg.luminosity, d, dop);
            let photon_e = (2.7 * K_B * agg.temperature.max(2.725)).max(1e-30);
            let expected = f * obs.integration_time / photon_e;
            Reading::ParticleCount {
                counts: stream.poisson(expected.min(1e12)),
                expected,
            }
        }
        Instrument::Interferometer => {
            let precision = obs.linear_resolution(d);
            let m = crate::solvers::quantum::measure_position(
                agg.mass,
                view.snapshot.offset,
                sep_dir,
                precision,
                stream,
            );
            Reading::Position {
                value: m.value,
                uncertainty: m.uncertainty,
                disturbance: m.disturbance,
            }
        }
        Instrument::MassSpectrometer => Reading::Composition {
            fractions: agg.composition.0,
        },
    }
}

/// Ways a user can act on the world.
///
/// The set is deliberately small and physical. Every one of them is expressed
/// as a change to a conserved quantity delivered at a place and a time, which
/// means none of them can violate the invariants the rest of the engine
/// depends on — a user cannot "cheat" the conservation checks through the
/// interaction API, only through the explicit authoring path, which says so.
#[derive(Debug, Clone, Copy)]
pub enum Interaction {
    /// Apply an impulse to a node. Momentum is taken from the observer's own
    /// budget so the global total is unchanged.
    Impulse { target: NodeIdx, dp: Vec3 },
    /// Deposit energy — a laser, a beam, a detonation.
    Deposit { target: NodeIdx, joules: f64, radius: f64 },
    /// Remove energy: a heat sink, a cooling laser.
    Extract { target: NodeIdx, joules: f64 },
    /// Add matter with a given composition and velocity.
    Inject {
        target: NodeIdx,
        mass: f64,
        composition: crate::state::Composition,
        velocity: Vec3,
    },
    /// Measure — which is an interaction, because measurement disturbs.
    Measure {
        target: NodeIdx,
        instrument: Instrument,
        quantity: Quantity,
    },
    /// Pin a region so its detail persists even when nobody is looking.
    Pin { target: NodeIdx },
    /// Explicit authoring: set a bulk property directly. Flagged in the audit
    /// log because it is the one path that can break conservation, and the
    /// engine reports the discontinuity rather than hiding it.
    Author {
        target: NodeIdx,
        property: Property,
        value: f64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Property {
    Mass,
    Temperature,
    Radius,
    Charge,
    Luminosity,
}

/// Record of an authoring action that broke conservation, so the audit trail
/// can distinguish "the engine drifted" from "the user changed things".
#[derive(Debug, Clone, Copy)]
pub struct AuthorEvent {
    pub key: PathKey,
    pub property: Property,
    pub delta_energy: f64,
    pub time: f64,
}
