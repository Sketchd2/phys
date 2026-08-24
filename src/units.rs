//! SI constants and the scale ladder.
//!
//! Everything internal is SI. The engine spans 36 orders of magnitude in
//! length (kiloparsecs to femtometres) and 51 in time (Gyr to zeptoseconds),
//! so *no single coordinate system is ever used globally* — see `coords.rs`.

#![allow(clippy::excessive_precision)]

// ---- fundamental constants (CODATA 2018) --------------------------------
pub const C: f64 = 299_792_458.0; // m/s, exact
pub const C2: f64 = C * C;
pub const G: f64 = 6.674_30e-11; // m^3 kg^-1 s^-2
pub const H_BAR: f64 = 1.054_571_817e-34; // J s
pub const H_PLANCK: f64 = 6.626_070_15e-34; // J s, exact
pub const K_B: f64 = 1.380_649e-23; // J/K, exact
pub const E_CHARGE: f64 = 1.602_176_634e-19; // C, exact
pub const EPS0: f64 = 8.854_187_812_8e-12; // F/m
pub const MU0: f64 = 1.256_637_062_12e-6; // N/A^2
pub const K_COULOMB: f64 = 8.987_551_792_3e9; // 1/(4 pi eps0)
pub const SIGMA_SB: f64 = 5.670_374_419e-8; // W m^-2 K^-4
pub const N_AVOGADRO: f64 = 6.022_140_76e23; // 1/mol
pub const A_RAD: f64 = 7.565_733e-16; // radiation constant, J m^-3 K^-4

// ---- masses --------------------------------------------------------------
pub const M_PROTON: f64 = 1.672_621_923_69e-27;
pub const M_NEUTRON: f64 = 1.674_927_498_04e-27;
pub const M_ELECTRON: f64 = 9.109_383_701_5e-31;
pub const AMU: f64 = 1.660_539_066_60e-27;
pub const M_SUN: f64 = 1.988_47e30;
pub const M_EARTH: f64 = 5.972_2e24;

// ---- lengths / times -----------------------------------------------------
pub const AU: f64 = 1.495_978_707e11;
pub const PARSEC: f64 = 3.085_677_581_49e16;
pub const KPC: f64 = 1e3 * PARSEC;
pub const LIGHT_YEAR: f64 = 9.460_730_472_580_8e15;
pub const R_SUN: f64 = 6.957e8;
pub const R_EARTH: f64 = 6.371e6;
pub const BOHR: f64 = 5.291_772_109_03e-11;
pub const FEMTOMETRE: f64 = 1e-15;
pub const YEAR: f64 = 3.155_695_2e7; // Julian year, s
pub const MYR: f64 = 1e6 * YEAR;
pub const GYR: f64 = 1e9 * YEAR;

// ---- energies ------------------------------------------------------------
pub const EV: f64 = E_CHARGE; // J
pub const KEV: f64 = 1e3 * EV;
pub const MEV: f64 = 1e6 * EV;
pub const GEV: f64 = 1e9 * EV;

/// Compton wavelength of the electron — below this, a "particle position" is
/// not a meaningful classical concept and the engine switches representation.
pub const LAMBDA_COMPTON_E: f64 = 2.426_310_238_67e-12;

/// The scale ladder.
///
/// Each tier owns a *representation*, a *solver*, and a *timestep law*. A tier
/// is not merely a zoom level: crossing a tier boundary changes what the state
/// vector means (a `Stellar` node's "temperature" is a sub-grid ISM
/// temperature; a `Nuclear` node has no temperature at all, it has excitation
/// levels). Prolongation and restriction operators (`prolong.rs`) translate
/// between adjacent tiers and are required to conserve the invariant set
/// exactly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Tier {
    /// 10^19 – 10^21 m. Collisionless N-body + dark matter halo. Myr steps.
    Galactic = 0,
    /// 10^13 – 10^19 m. GMCs, star clusters, ISM hydrodynamics. kyr steps.
    Stellar = 1,
    /// 10^6 – 10^13 m. Stars, planets, orbits, stellar structure. s–hr steps.
    Planetary = 2,
    /// 10^-3 – 10^6 m. Continuum: solids, fluids, thermodynamics. ms steps.
    Continuum = 3,
    /// 10^-9 – 10^-3 m. Molecular dynamics, force fields. fs steps.
    Molecular = 4,
    /// 10^-11 – 10^-9 m. Atomic/electronic structure. as steps.
    Atomic = 5,
    /// < 10^-14 m. Nuclear and subatomic. Statistical, not trajectorial. zs.
    Nuclear = 6,
}

pub const TIER_COUNT: usize = 7;

impl Tier {
    pub const ALL: [Tier; TIER_COUNT] = [
        Tier::Galactic,
        Tier::Stellar,
        Tier::Planetary,
        Tier::Continuum,
        Tier::Molecular,
        Tier::Atomic,
        Tier::Nuclear,
    ];

    pub fn index(self) -> usize {
        self as usize
    }

    pub fn from_index(i: usize) -> Tier {
        Tier::ALL[i.min(TIER_COUNT - 1)]
    }

    /// One tier finer, saturating at `Nuclear`.
    pub fn finer(self) -> Tier {
        Tier::from_index(self.index() + 1)
    }

    /// One tier coarser, saturating at `Galactic`.
    pub fn coarser(self) -> Tier {
        Tier::from_index(self.index().saturating_sub(1))
    }

    pub fn name(self) -> &'static str {
        match self {
            Tier::Galactic => "galactic",
            Tier::Stellar => "stellar",
            Tier::Planetary => "planetary",
            Tier::Continuum => "continuum",
            Tier::Molecular => "molecular",
            Tier::Atomic => "atomic",
            Tier::Nuclear => "nuclear",
        }
    }

    /// Characteristic length of the tier, in metres. Used to pick a tier for a
    /// given observation aperture and to set the tree's opening criterion.
    pub fn length(self) -> f64 {
        match self {
            Tier::Galactic => 1e20,
            Tier::Stellar => 1e16,
            Tier::Planetary => 1e9,
            Tier::Continuum => 1e0,
            Tier::Molecular => 1e-9,
            Tier::Atomic => 1e-10,
            Tier::Nuclear => 1e-15,
        }
    }

    /// Characteristic integration step of the tier, in seconds. These are the
    /// *natural* steps; the scheduler shrinks them further under the causal
    /// constraint (`causal.rs`).
    pub fn dt(self) -> f64 {
        match self {
            Tier::Galactic => 1e5 * YEAR,
            Tier::Stellar => 1e2 * YEAR,
            Tier::Planetary => 1e2,
            Tier::Continuum => 1e-3,
            Tier::Molecular => 1e-15,
            Tier::Atomic => 1e-18,
            Tier::Nuclear => 1e-21,
        }
    }

    /// Light-crossing time of the tier's characteristic length. This is the
    /// engine's lookahead: no node can influence another before this elapses.
    pub fn light_crossing(self) -> f64 {
        self.length() / C
    }

    /// Lower bound of the tier's length range, in metres.
    ///
    /// A tier is a *physics regime*, not a tree level. Many refinements happen
    /// within one tier — a molecular cloud splits into clumps into cores into
    /// protostars, all of it `Stellar` — and the tier only changes when the
    /// characteristic size crosses one of these boundaries and a different
    /// description becomes appropriate. Conflating tiers with tree depth is a
    /// tempting simplification and a badly wrong one: it gives you seven levels
    /// of refinement between a galaxy and a nucleus, when the mass ratio alone
    /// demands more than twenty.
    pub fn floor(self) -> f64 {
        match self {
            Tier::Galactic => 1e18,  // ~30 pc: below this, self-gravitating clouds
            Tier::Stellar => 1e12,   // ~7 AU: below this, individual bodies
            Tier::Planetary => 1e4,  // 10 km: below this, bulk material
            Tier::Continuum => 1e-8, // 10 nm: below this, discrete molecules
            Tier::Molecular => 3e-10,
            Tier::Atomic => 1e-14,
            Tier::Nuclear => 0.0,
        }
    }

    /// The tier appropriate to an object of this size. Derived from the
    /// physical scale, never from how deep in the tree the object happens to
    /// sit.
    pub fn containing(metres: f64) -> Tier {
        for t in Tier::ALL {
            if metres >= t.floor() {
                return t;
            }
        }
        Tier::Nuclear
    }

    /// The tier that must be materialised to resolve detail at `metres`.
    pub fn for_resolution(metres: f64) -> Tier {
        Tier::containing(metres)
    }
}

/// Chemical species tracked in every aggregate state.
///
/// Eight buckets is a deliberate compromise: it is enough to drive nuclear
/// burning, cooling curves, opacity and chemistry, and it fits a cache line
/// alongside the rest of the state. Finer speciation is materialised on demand
/// at `Molecular` and below.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Species {
    Hydrogen = 0,
    Helium = 1,
    Carbon = 2,
    Nitrogen = 3,
    Oxygen = 4,
    Silicon = 5,
    Iron = 6,
    /// Everything else: heavier nuclei plus, at galactic tier, cold dark
    /// matter (which is inert under every force but gravity).
    Other = 7,
}

pub const NSPECIES: usize = 8;

impl Species {
    pub const ALL: [Species; NSPECIES] = [
        Species::Hydrogen,
        Species::Helium,
        Species::Carbon,
        Species::Nitrogen,
        Species::Oxygen,
        Species::Silicon,
        Species::Iron,
        Species::Other,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Species::Hydrogen => "H",
            Species::Helium => "He",
            Species::Carbon => "C",
            Species::Nitrogen => "N",
            Species::Oxygen => "O",
            Species::Silicon => "Si",
            Species::Iron => "Fe",
            Species::Other => "Z",
        }
    }

    /// Atomic number.
    pub fn z(self) -> f64 {
        match self {
            Species::Hydrogen => 1.0,
            Species::Helium => 2.0,
            Species::Carbon => 6.0,
            Species::Nitrogen => 7.0,
            Species::Oxygen => 8.0,
            Species::Silicon => 14.0,
            Species::Iron => 26.0,
            Species::Other => 30.0,
        }
    }

    /// Mass number (mean, for the lumped bucket).
    pub fn a(self) -> f64 {
        match self {
            Species::Hydrogen => 1.0,
            Species::Helium => 4.0,
            Species::Carbon => 12.0,
            Species::Nitrogen => 14.0,
            Species::Oxygen => 16.0,
            Species::Silicon => 28.0,
            Species::Iron => 56.0,
            Species::Other => 65.0,
        }
    }

    pub fn mass_kg(self) -> f64 {
        self.a() * AMU
    }

    /// Binding energy per nucleon, MeV — drives fusion energetics and tells the
    /// engine that iron is the floor.
    pub fn binding_per_nucleon_mev(self) -> f64 {
        match self {
            Species::Hydrogen => 0.0,
            Species::Helium => 7.074,
            Species::Carbon => 7.680,
            Species::Nitrogen => 7.476,
            Species::Oxygen => 7.976,
            Species::Silicon => 8.447,
            Species::Iron => 8.790,
            Species::Other => 8.60,
        }
    }
}

/// Format a length with an appropriate astronomical or atomic unit.
pub fn fmt_length(m: f64) -> String {
    let a = m.abs();
    if a >= 0.1 * KPC {
        format!("{:.3} kpc", m / KPC)
    } else if a >= 0.1 * PARSEC {
        format!("{:.3} pc", m / PARSEC)
    } else if a >= 0.1 * AU {
        format!("{:.3} AU", m / AU)
    } else if a >= 1e3 {
        format!("{:.3} km", m / 1e3)
    } else if a >= 1e-3 {
        format!("{:.3} m", m)
    } else if a >= 1e-9 {
        format!("{:.3} um", m * 1e6)
    } else if a >= 1e-12 {
        format!("{:.3} nm", m * 1e9)
    } else {
        format!("{:.3} fm", m * 1e15)
    }
}

pub fn fmt_time(s: f64) -> String {
    let a = s.abs();
    if a >= 0.1 * GYR {
        format!("{:.3} Gyr", s / GYR)
    } else if a >= 0.1 * MYR {
        format!("{:.3} Myr", s / MYR)
    } else if a >= 0.1 * YEAR {
        format!("{:.3} yr", s / YEAR)
    } else if a >= 1.0 {
        format!("{:.3} s", s)
    } else if a >= 1e-9 {
        format!("{:.3} ns", s * 1e9)
    } else if a >= 1e-15 {
        format!("{:.3} fs", s * 1e15)
    } else {
        format!("{:.3e} s", s)
    }
}

pub fn fmt_mass(kg: f64) -> String {
    let a = kg.abs();
    if a >= 1e-3 * M_SUN {
        format!("{:.4} Msun", kg / M_SUN)
    } else if a >= 1e-6 * M_EARTH {
        format!("{:.4} Mearth", kg / M_EARTH)
    } else if a >= 1e-3 {
        format!("{:.4} kg", kg)
    } else if a >= 1e-24 {
        format!("{:.4} amu", kg / AMU)
    } else {
        format!("{:.4e} kg", kg)
    }
}
