//! What condensed matter does when it is squeezed.
//!
//! # Why this exists
//!
//! [`crate::state::Matter::pressure`] is an ideal gas plus radiation, and until
//! this module it was the only equation of state in the engine. Handed a bucket
//! of water it answered 4.0x10^8 Pa; handed a solid, it burst the solid from
//! the inside. `docs/PLAY.md` §4.2 measured the first and §3.3 the second, and
//! Phase 5's second piece is the fix: *a liquid equation of state, so
//! `pressure()` stops returning 4x10^8 Pa for a bucket of water*.
//!
//! # The law, and where its one number comes from
//!
//! One form for both condensed phases, Murnaghan's:
//!
//! ```text
//!   p(rho) = (K / K') ((rho / rho_0)^K' - 1)
//! ```
//!
//! which for a liquid with `K' = 7` is the Tait–Cole equation every
//! weakly-compressible SPH code uses, and for a solid with `K' = 4` is the
//! standard first-order account of compression. The exponents are the owner's
//! decision for Phase 5, taken with this module's alternatives on the table.
//!
//! **The bulk modulus `K` is derived, and from different energies for the two
//! phases**, because what holds a liquid together is not what holds a molecule
//! together:
//!
//! - **A liquid's is its vaporisation energy density.** Taking a liquid apart
//!   means separating its molecules, and the cost of that per unit volume is the
//!   latent heat of vaporisation times the density — which `chem::react` already
//!   derives from Trouton's rule. A liquid's bulk modulus is of the order of its
//!   cohesive energy density, and for water the two agree closely: 2.2x10^9 Pa
//!   measured against 2.1x10^9 of energy density.
//!
//!   The substance's `cohesive_energy` is **not** the right number, and was
//!   measured before this was written: for water it is 9.07 eV a molecule —
//!   the O–H bonds — which through `Material::of`'s stiffness rule gives
//!   2.4x10^11 Pa and a sound speed of 12 km/s, a hundred times too stiff. What
//!   melts and boils in water is the hydrogen bonding between molecules, which
//!   is exactly what Trouton prices.
//!
//! - **A solid's is its stiffness.** `Material::of` already derives Young's
//!   modulus from the cohesive energy density, and a bulk modulus follows at
//!   the Poisson's ratio of 0.3 the engine uses everywhere else
//!   (`Material::yield_stress`): `K = E / 3(1 - 2 nu) = E / 1.2`. Taken at full
//!   packing — the equation of state is of the *substance*, and how much of a
//!   node's volume it fills is the node's business, not the law's.
//!
//! Measured on the engine's own water: 4.07x10^9 Pa and 1569 m/s, against real
//! water's 2.2x10^9 and 1500. The modulus is high by the same 65% the engine's
//! estimate of water's *density* is high by (1653 kg/m^3), and the sound speed
//! — which divides one by the other — lands within five per cent.
//!
//! # What it is not
//!
//! **Not a function of temperature.** Tait and Murnaghan are barotropic: the
//! pressure depends on density alone. Thermal pressure in a condensed phase is
//! real (a Grüneisen term) and small against the elastic one at any
//! temperature where the phase exists, which is the standard reason to leave it
//! out. The energy that compression puts in is carried — see
//! [`Condensed::elastic_energy`] — so leaving it out of the *pressure* does not
//! leak it out of the books.
//!
//! **Not a tension carrier at a free surface.** Below its rest density a liquid
//! would, by the formula, pull: `-K/K'` at zero density, 5.8x10^8 Pa for water.
//! Real water cavitates long before that, and a node or a particle that is
//! *less* dense than its rest density is one that is partly empty rather than
//! stretched. [`Condensed::pressure`] floors at zero, which is what a free
//! surface is.

use crate::chem::{Mixture, Phase, Properties, Registry};
use crate::state::Matter;

/// A liquid's stiffening exponent: Tait–Cole.
pub const LIQUID_STIFFENING: f64 = 7.0;
/// A solid's stiffening exponent: Murnaghan's first-order value.
pub const SOLID_STIFFENING: f64 = 4.0;
/// Poisson's ratio for a dense solid, the same 0.3 `Material::yield_stress`
/// uses to turn Young's modulus into a shear modulus.
pub const POISSON: f64 = 0.3;

/// A condensed phase's equation of state: a rest density, a bulk modulus, and
/// how quickly the modulus rises under compression.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Condensed {
    /// Density at zero pressure, kg/m^3.
    pub rest_density: f64,
    /// Bulk modulus at zero pressure, Pa.
    pub bulk_modulus: f64,
    /// `dK/dp`, dimensionless.
    pub stiffening: f64,
}

impl Condensed {
    /// A substance's liquid, from its vaporisation energy density.
    pub fn liquid(props: &Properties) -> Option<Condensed> {
        let rho = props.density;
        let l = crate::chem::react::heat_of_vaporisation(props);
        let k = l * rho;
        (rho > 0.0 && k > 0.0 && k.is_finite()).then_some(Condensed {
            rest_density: rho,
            bulk_modulus: k,
            stiffening: LIQUID_STIFFENING,
        })
    }

    /// A substance's solid, from the stiffness `Material::of` derives.
    pub fn solid(props: &Properties) -> Option<Condensed> {
        let rho = props.density;
        let e = crate::material::dense_stiffness(props);
        let k = e / (3.0 * (1.0 - 2.0 * POISSON));
        (rho > 0.0 && k > 0.0 && k.is_finite()).then_some(Condensed {
            rest_density: rho,
            bulk_modulus: k,
            stiffening: SOLID_STIFFENING,
        })
    }

    /// The phase a substance is in, as a condensed equation of state, or `None`
    /// for a gas.
    pub fn of(props: &Properties, phase: Phase) -> Option<Condensed> {
        match phase {
            Phase::Liquid => Condensed::liquid(props),
            Phase::Solid => Condensed::solid(props),
            _ => None,
        }
    }

    /// Pressure at a density, Pa, floored at zero — see the module doc.
    pub fn pressure(&self, rho: f64) -> f64 {
        if !(rho > self.rest_density) {
            return 0.0;
        }
        let n = self.stiffening;
        self.bulk_modulus / n * ((rho / self.rest_density).powf(n) - 1.0)
    }

    /// `sqrt(dp/drho)` at a density, m/s — taken at the rest density for
    /// anything less dense, because the signal in a partly empty region is
    /// still carried by the material that is there.
    pub fn sound_speed(&self, rho: f64) -> f64 {
        let r = rho.max(self.rest_density);
        let n = self.stiffening;
        (self.bulk_modulus / self.rest_density * (r / self.rest_density).powf(n - 1.0)).sqrt()
    }

    /// Energy stored per kilogram by compressing from rest to `rho`, J/kg:
    /// `integral p / rho^2 drho`. Zero at or below rest.
    pub fn elastic_energy(&self, rho: f64) -> f64 {
        if !(rho > self.rest_density) {
            return 0.0;
        }
        let n = self.stiffening;
        let r0 = self.rest_density;
        let k = self.bulk_modulus;
        // integral_{r0}^{rho} K/n ((r/r0)^n - 1) / r^2 dr
        //   = K/n [ r^(n-1) / ((n-1) r0^n) + 1/r ]_{r0}^{rho}
        let f = |r: f64| r.powf(n - 1.0) / ((n - 1.0) * r0.powf(n)) + 1.0 / r;
        k / n * (f(rho) - f(r0))
    }

    /// Density at which this phase carries `p`, kg/m^3 — the inverse of
    /// [`Condensed::pressure`] for `p >= 0`.
    pub fn density_at(&self, p: f64) -> f64 {
        let n = self.stiffening;
        self.rest_density * (1.0 + n * p.max(0.0) / self.bulk_modulus).powf(1.0 / n)
    }
}

/// What a node's matter, or one of its bodies, answers to when squeezed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Eos {
    /// The ideal gas plus radiation `Matter::pressure` has always been. Also
    /// the answer for matter nobody has described, which is almost everything
    /// in a galaxy and genuinely is a gas.
    Gas,
    /// A condensed phase, or a blend of them.
    Condensed(Condensed),
}

impl Eos {
    /// The equation of state a mixture's matter answers to.
    ///
    /// **Condensed wherever `Matter::gas_law_applies` says the gas law does
    /// not**, and made of the condensed pools. Several condensed pools
    /// blend as a Reuss average — compliances add by volume, which is what
    /// pressure equal in every phase means — and the rest density is the mass
    /// over the volume the pools fill.
    pub fn of_mixture(mix: &Mixture, reg: &Registry) -> Eos {
        if mix.is_empty() {
            return Eos::Gas;
        }
        let mut mass = 0.0;
        let mut volume = 0.0;
        let mut compliance = 0.0;
        let mut stiffening = 0.0;
        let mut total = 0.0;
        for pool in mix.entries() {
            total += pool.fraction;
            let Some(s) = reg.get(pool.substance) else { continue };
            let Some(c) = Condensed::of(&s.props, pool.phase) else { continue };
            let v = pool.fraction / c.rest_density;
            mass += pool.fraction;
            volume += v;
            compliance += v / c.bulk_modulus;
            stiffening += pool.fraction * c.stiffening;
        }
        // The line `Matter::gas_law_applies` draws, so the two cannot disagree
        // about which matter the gas law describes.
        if !(total > 0.0) || mix.in_phase(Phase::Gas) >= 0.5 || !(volume > 0.0) || !(compliance > 0.0)
        {
            return Eos::Gas;
        }
        Eos::Condensed(Condensed {
            rest_density: mass / volume,
            bulk_modulus: volume / compliance,
            stiffening: stiffening / mass,
        })
    }

    /// The equation of state of what a structure is *not* made of: the node's
    /// liquid and gas pools, with its solids — which the structure is — left
    /// out. A patch of ground with water on it has loose contents that are
    /// water, and pricing them as the blend priced them at 2436 kg/m^3 of
    /// rest density where water rests at 1653.
    pub fn of_loose(mix: &Mixture, reg: &Registry) -> Eos {
        let mut loose = Mixture::new();
        for p in mix.entries() {
            if p.phase != Phase::Solid {
                loose.add(p.substance, p.phase, p.fraction);
            }
        }
        Eos::of_mixture(&loose, reg)
    }

    pub fn of_matter(m: &Matter, reg: &Registry) -> Eos {
        Eos::of_mixture(&m.mixture, reg)
    }

    pub fn condensed(&self) -> Option<Condensed> {
        match self {
            Eos::Condensed(c) => Some(*c),
            Eos::Gas => None,
        }
    }

    /// Signal speed through a node's matter, m/s.
    ///
    /// For a gas, exactly [`Matter::sound_speed`], cap and all. For condensed
    /// matter, the material's own — at the node's bulk density if it is
    /// compressed, and at rest if the node is partly empty.
    pub fn sound_speed(&self, m: &Matter) -> f64 {
        match self {
            Eos::Gas => m.sound_speed(),
            Eos::Condensed(c) => c.sound_speed(m.density()),
        }
    }
}
