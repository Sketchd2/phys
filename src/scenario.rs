//! Starting points: what the world is, before anybody looks at it.
//!
//! A scenario is one [`Aggregate`] and one refinement policy at one tier. That
//! is genuinely all it takes, because everything else the engine does is
//! *derived*: the detail is generated on demand from the aggregate, the solver
//! is chosen by the tier, the timestep by the physics, and the whole ladder
//! below is produced by refining. A galaxy and an iron nucleus differ by
//! forty-five orders of magnitude in size and by nothing at all in structure.
//!
//! # Why these are worth having as data
//!
//! The engine could always be pointed at any scale — [`crate::engine::galaxy`]
//! is one of these, written by hand. What was missing was a way to *say* which
//! one without writing code, and a way for a viewer to offer the choice. Every
//! scenario below is a few lines of bulk state, and the interesting claim is
//! that nothing else is needed: no per-scale renderer, no per-scale solver
//! selection, no per-scale camera. Descending from a galaxy to a nucleus and
//! loading a nucleus directly land in the same place.
//!
//! # The energy budgets are not decoration
//!
//! An aggregate that says it is a star has to carry a star's internal energy,
//! binding energy and angular momentum, because those are what the
//! materialisation is constrained to reproduce. Give a nucleus a thermal energy
//! computed as though it were a gas and the projection will still close its
//! books — on a nucleus whose nucleons are moving at a hundredth of the speed
//! they should be. So each budget below is derived from what the object *is*:
//! the virial theorem for anything self-gravitating, Fermi motion for a
//! nucleus, the equipartition of a gas for a gas.

use crate::engine::{default_spec, galaxy};
use crate::math::v3;
use crate::prolong::{MassSpectrum, Profile, ProlongSpec};
use crate::state::{Aggregate, BodyKind, Composition};
use crate::tree::Tree;
use crate::units::*;

/// One thing the engine can be pointed at.
pub struct Scenario {
    pub name: &'static str,
    /// One line on what it is and what is interesting about it.
    pub blurb: &'static str,
    pub tier: Tier,
    /// Characteristic size, metres. For labelling, not for physics.
    pub scale: f64,
    build: fn(u64) -> Tree,
}

impl Scenario {
    pub fn build(&self, world_seed: u64) -> Tree {
        (self.build)(world_seed)
    }
}

/// Everything on the shelf, coarsest first.
pub const ALL: &[Scenario] = &[
    Scenario {
        name: "Spiral galaxy",
        blurb: "10^9 stars in a rotating disc, held together by a halo that is not made of the thing being refined.",
        tier: Tier::Galactic,
        scale: 15.0 * KPC,
        build: build_galaxy,
    },
    Scenario {
        name: "Molecular cloud",
        blurb: "10^5 solar masses of cold gas at 20 K, turbulent, on the edge of collapsing into stars.",
        tier: Tier::Stellar,
        scale: 20.0 * PARSEC,
        build: build_cloud,
    },
    Scenario {
        name: "The Sun",
        blurb: "One solar mass at hydrostatic equilibrium. Refine it and the pp chain lights up in the core.",
        tier: Tier::Planetary,
        scale: R_SUN,
        build: build_star,
    },
    Scenario {
        name: "Rocky planet",
        blurb: "An Earth: iron core, silicate mantle, and a virial budget that knows the difference.",
        tier: Tier::Planetary,
        scale: 6.371e6,
        build: build_planet,
    },
    Scenario {
        name: "Granite block",
        blurb: "A cubic metre of rock at room temperature — bulk matter, where thermodynamics is the whole physics.",
        tier: Tier::Continuum,
        scale: 0.5,
        build: build_rock,
    },
    Scenario {
        name: "Water vapour",
        blurb: "A cube of steam at 400 K. Molecular dynamics with covalent bonds that form and break.",
        tier: Tier::Molecular,
        scale: 3.0e-9,
        build: build_vapour,
    },
    Scenario {
        name: "Carbon atom",
        blurb: "One atom, and the twelve nucleons inside it. Below here the engine stops producing trajectories, because there is nothing else there to describe.",
        tier: Tier::Atomic,
        scale: 7.0e-11,
        build: build_atom,
    },
    Scenario {
        name: "Iron nucleus",
        blurb: "Fifty-six nucleons at nuclear density, moving at a fifth of light speed because the exclusion principle says so.",
        tier: Tier::Nuclear,
        scale: 4.6e-15,
        build: build_nucleus,
    },
];

pub fn by_index(i: usize) -> &'static Scenario {
    &ALL[i.min(ALL.len() - 1)]
}

fn build_galaxy(seed: u64) -> Tree {
    galaxy(seed, 1e9)
}

/// A giant molecular cloud: cold, turbulent, marginally bound.
///
/// The internal energy is *not* thermal. At 20 K the sound speed is 0.3 km/s
/// and the observed line widths are ten times that: a cloud's kinetic budget is
/// supersonic turbulence, and a cloud given only its thermal energy would
/// collapse in a free-fall time instead of living for tens of them.
fn build_cloud(seed: u64) -> Tree {
    let mass = 1.0e5 * M_SUN;
    let radius = 20.0 * PARSEC;
    let mut agg = Aggregate::neutral(mass, radius, 20.0, Composition::primordial());
    let sigma = (G * mass / radius).sqrt();
    agg.internal_energy = 0.5 * mass * sigma * sigma;
    agg.binding_energy = -0.6 * G * mass * mass / radius;
    agg.spin = v3(0.0, 0.0, 0.25 * mass * sigma * radius);
    let spec = ProlongSpec {
        count: 8_000,
        profile: Profile::Plummer,
        spectrum: MassSpectrum::Equal,
        kind: BodyKind::GasParcel,
        composition_scatter: 0.05,
        turbulent_fraction: 0.7,
    };
    Tree::new(seed, agg, Tier::Stellar, spec)
}

/// A star, sized by the virial theorem rather than by a number typed in.
///
/// For a self-gravitating ball of ideal gas the internal energy is exactly half
/// the magnitude of the binding energy. That single relation fixes the whole
/// budget from the mass and the radius, and it is why a star that is made
/// smaller gets *hotter* — which is the reason stars work at all.
fn build_star(seed: u64) -> Tree {
    let mass = M_SUN;
    let radius = R_SUN;
    let binding = -0.6 * G * mass * mass / radius;
    let mut agg = Aggregate::neutral(mass, radius, 5.8e3, Composition::solar());
    agg.internal_energy = -0.5 * binding;
    agg.binding_energy = binding;
    agg.luminosity = 3.828e26;
    agg.spin = v3(0.0, 0.0, 1.9e41);
    Tree::new(seed, agg, Tier::Planetary, default_spec(Tier::Planetary))
}

/// A rocky planet. Same virial relation, a thousandth the mass, and a
/// composition that is mostly iron and silicon rather than mostly hydrogen.
fn build_planet(seed: u64) -> Tree {
    let mass = 5.972e24;
    let radius = 6.371e6;
    let binding = -0.6 * G * mass * mass / radius;
    let mut comp = [0.0; NSPECIES];
    comp[Species::Iron as usize] = 0.32;
    comp[Species::Silicon as usize] = 0.15;
    comp[Species::Oxygen as usize] = 0.30;
    comp[Species::Other as usize] = 0.23;
    let mut agg = Aggregate::neutral(mass, radius, 2000.0, Composition(comp).normalised());
    // A planet is not an ideal gas: most of its binding is held by material
    // strength and electron degeneracy, not by heat. Booking the full virial
    // internal energy would have the Earth at 10^5 K throughout.
    agg.internal_energy = -0.1 * binding;
    agg.binding_energy = binding;
    agg.spin = v3(0.0, 0.0, 7.05e33);
    Tree::new(seed, agg, Tier::Planetary, default_spec(Tier::Planetary))
}

/// A cubic metre of granite. Cold, dense, and held together by chemistry rather
/// than by gravity — so the binding energy is the lattice, not `GM^2/R`.
fn build_rock(seed: u64) -> Tree {
    let density = 2650.0;
    let half: f64 = 0.5;
    let volume = (2.0 * half).powi(3);
    let mass = density * volume;
    let mut comp = [0.0; NSPECIES];
    comp[Species::Oxygen as usize] = 0.47;
    comp[Species::Silicon as usize] = 0.28;
    comp[Species::Iron as usize] = 0.05;
    comp[Species::Other as usize] = 0.20;
    let mut agg = Aggregate::neutral(mass, half * 3f64.sqrt(), 290.0, Composition(comp).normalised());
    // Cohesive energy of a silicate, a few electron volts per atom.
    agg.binding_energy = -mass * agg.composition.nucleons_per_kg() / 20.0 * 5.0 * EV;
    Tree::new(seed, agg, Tier::Continuum, default_spec(Tier::Continuum))
}

/// Water vapour: a box of molecules hot enough to be a gas and cool enough that
/// the bonds hold. This is the tier where chemistry happens.
fn build_vapour(seed: u64) -> Tree {
    let count = 512.0;
    let mass_per = 18.015 * AMU;
    let mass = count * mass_per;
    let radius = 3.0e-9;
    let mut comp = [0.0; NSPECIES];
    comp[Species::Hydrogen as usize] = 2.0 * 1.008 / 18.015;
    comp[Species::Oxygen as usize] = 15.999 / 18.015;
    let mut agg = Aggregate::neutral(mass, radius, 400.0, Composition(comp).normalised());
    // Bound, and by a lot: two O-H bonds per molecule at 4.8 eV each.
    agg.binding_energy = -count * 2.0 * 4.81 * EV;
    Tree::new(seed, agg, Tier::Molecular, default_spec(Tier::Molecular))
}

/// A carbon atom, which refines into its own nucleons rather than into more
/// atoms.
///
/// The refinement table's atomic entry is "a *molecule* resolves into atoms",
/// which is the right policy for a node the size of a molecule and the wrong
/// one for a node that is a single atom: there is exactly one atom in a carbon
/// atom, and refining it that way shows you one body. What is inside an atom is
/// its nucleus, so this scenario carries the nuclear policy at the atomic tier
/// — the one place where what a node *is* and how big it is disagree.
fn build_atom(seed: u64) -> Tree {
    let mass = 12.011 * AMU;
    let radius = 7.0e-11;
    let mut agg = Aggregate::neutral(mass, radius, 300.0, Composition::pure(Species::Carbon));
    // Total electronic binding of neutral carbon, about 1030 eV.
    agg.binding_energy = -1030.0 * EV;
    agg.internal_energy = 1030.0 * EV * 0.5;
    Tree::new(seed, agg, Tier::Atomic, default_spec(Tier::Nuclear))
}

/// An iron nucleus, at the top of the binding energy curve.
///
/// Its internal energy is Fermi motion, not heat. Nucleons in a nucleus are
/// confined to a few femtometres, and the exclusion principle then puts them at
/// around 30 MeV of kinetic energy each — a fifth of light speed — whatever the
/// temperature is. A nucleus given a thermal budget at 300 K would have its
/// nucleons sitting still, which is not a nucleus.
fn build_nucleus(seed: u64) -> Tree {
    let a: f64 = 56.0;
    let mass = 55.845 * AMU;
    let radius = 1.2e-15 * a.cbrt();
    let mut agg = Aggregate::neutral(mass, radius, 1.0e9, Composition::pure(Species::Iron));
    agg.binding_energy = -8.79 * MEV * a;
    agg.internal_energy = 0.6 * 33.0 * MEV * a;
    agg = agg.with_charge(26.0 * E_CHARGE);
    Tree::new(seed, agg, Tier::Nuclear, default_spec(Tier::Nuclear))
}
