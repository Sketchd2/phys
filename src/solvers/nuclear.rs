//! Nuclear processes: what makes stars shine and what makes matter change.
//!
//! Everything here is *rate-based*. There is no attempt to integrate the strong
//! interaction — that would be lattice QCD, which costs supercomputer-months
//! per femtometre. Instead the engine uses measured cross sections and decay
//! constants, and samples events from them. This is the same choice every
//! stellar evolution code makes, and it is not a compromise: the tabulated
//! rates *are* the experimental facts, and a first-principles calculation would
//! be reproducing them with less accuracy.
//!
//! The consistency requirement is that composition changes here must be
//! reflected in the conserved tuple. Fusing hydrogen to helium reduces the
//! rest mass and releases the difference as energy, and both sides of that
//! trade are booked.

use crate::rng::Stream;
use crate::state::Composition;
use crate::units::*;

/// Energy released by converting `mass` of `from` into `to`, from the binding
/// energy curve. Positive means exothermic.
///
/// This is the *only* place composition changes are allowed to create energy,
/// and it is exactly the difference in nuclear binding — which is why iron is
/// the end of the line, and why the engine gets that for free rather than
/// having it hard-coded.
pub fn fusion_energy(from: Species, to: Species, mass: f64) -> f64 {
    let nucleons = mass / from.mass_kg() * from.a();
    let delta_per_nucleon = (to.binding_per_nucleon_mev() - from.binding_per_nucleon_mev()) * MEV;
    nucleons * delta_per_nucleon
}

/// Gamow peak energy for a reaction between nuclei of charges `z1`, `z2` at
/// temperature `t`. The narrow window where tunnelling probability and the
/// Maxwell tail overlap — the reason stellar burning is so temperature
/// sensitive.
pub fn gamow_peak(z1: f64, z2: f64, reduced_mass: f64, t: f64) -> f64 {
    let kt = K_B * t;
    let b = std::f64::consts::PI * z1 * z2 * E_CHARGE * E_CHARGE / (EPS0 * H_PLANCK)
        * (reduced_mass / 2.0).sqrt();
    (b * b * kt * kt / 4.0).powf(1.0 / 3.0)
}

/// Proton-proton chain energy generation rate, W/kg.
///
/// The standard `epsilon ∝ rho X^2 T^4` scaling with the screening and
/// temperature-exponent fit used in stellar structure codes. At solar central
/// conditions (1.5e7 K, 1.5e5 kg/m^3, X=0.35) it returns a few times 10^-3
/// W/kg, which integrates to roughly a solar luminosity — the check in
/// `tests/nuclear.rs`.
pub fn pp_chain_rate(rho: f64, temperature: f64, x_hydrogen: f64) -> f64 {
    if temperature < 1e6 || x_hydrogen <= 0.0 {
        return 0.0;
    }
    let t6 = temperature / 1e6;
    let t9 = temperature / 1e9;
    // Clayton's form, converted to SI. In cgs it reads
    //     eps = 2.38e6 rho X^2 T6^(-2/3) exp(-33.80 T6^(-1/3))  erg/g/s
    // and the conversion carries three factors that are easy to drop:
    // rho g/cm^3 -> kg/m^3 (10^-3), erg/g -> J/kg (10^-4), giving 2.38e6 x
    // 10^-7 = 0.238. Getting this wrong by the 10^5 that the raw cgs
    // coefficient implies makes the Sun a hundred thousand times too bright,
    // which is at least an obvious failure.
    const PP_SI: f64 = 0.238;
    let g = 33.80 / t6.cbrt();
    PP_SI * rho * x_hydrogen * x_hydrogen / (t6 * t6).cbrt() * (-g).exp() * (1.0 + 0.012 * t9)
}

/// CNO cycle rate, W/kg. Takes over above ~1.8e7 K, and its far steeper
/// temperature dependence is why massive stars are convective in the core.
pub fn cno_rate(rho: f64, temperature: f64, x_hydrogen: f64, x_cno: f64) -> f64 {
    if temperature < 5e6 || x_hydrogen <= 0.0 || x_cno <= 0.0 {
        return 0.0;
    }
    let t6 = temperature / 1e6;
    // cgs 8.67e27 erg/g/s -> SI, same conversion as the pp chain.
    const CNO_SI: f64 = 8.67e20;
    let g = 152.28 / t6.cbrt();
    CNO_SI * rho * x_hydrogen * x_cno / (t6 * t6).cbrt() * (-g).exp()
}

/// Triple-alpha rate, W/kg. Helium burning, and the origin of carbon.
pub fn triple_alpha_rate(rho: f64, temperature: f64, y_helium: f64) -> f64 {
    if temperature < 1e8 || y_helium <= 0.0 {
        return 0.0;
    }
    let t8 = temperature / 1e8;
    // cgs 5.09e11 erg/g/s with rho^2 -> SI picks up 10^-6 from the density
    // squared and 10^-4 from the specific energy rate.
    const TRIPLE_ALPHA_SI: f64 = 5.09e1;
    TRIPLE_ALPHA_SI * rho * rho * y_helium.powi(3) / t8.powi(3) * (-44.027 / t8).exp()
}

/// Total nuclear energy generation for a parcel, W/kg, plus the composition
/// change it implies over `dt`.
pub struct BurnResult {
    pub power_per_kg: f64,
    pub new_composition: Composition,
    pub energy_released: f64,
}

pub fn burn(comp: Composition, rho: f64, temperature: f64, mass: f64, dt: f64) -> BurnResult {
    let x_h = comp.get(Species::Hydrogen);
    let y_he = comp.get(Species::Helium);
    let x_cno = comp.get(Species::Carbon) + comp.get(Species::Nitrogen) + comp.get(Species::Oxygen);

    let pp = pp_chain_rate(rho, temperature, x_h);
    let cno = cno_rate(rho, temperature, x_h, x_cno);
    let tri = triple_alpha_rate(rho, temperature, y_he);
    let total = pp + cno + tri;

    let mut c = comp.0;
    let energy = total * mass * dt;

    // Convert the released energy back into a mass transfer along the binding
    // curve, so composition and energy stay in step. H -> He releases
    // 0.7% of the rest mass; He -> C releases 0.07%.
    let h_to_he = (pp + cno) * mass * dt / (0.00712 * C2);
    let he_to_c = tri * mass * dt / (0.00076 * C2);
    let h_avail = c[Species::Hydrogen as usize] * mass;
    let he_avail = c[Species::Helium as usize] * mass;
    let dh = h_to_he.min(h_avail * 0.1);
    let dhe = he_to_c.min(he_avail * 0.1);
    if mass > 0.0 {
        c[Species::Hydrogen as usize] -= dh / mass;
        c[Species::Helium as usize] += dh / mass - dhe / mass;
        c[Species::Carbon as usize] += dhe / mass;
    }

    BurnResult {
        power_per_kg: total,
        new_composition: Composition(c).normalised(),
        energy_released: energy,
    }
}

/// A radioactive species the engine can follow explicitly when a user zooms in
/// on one. Half-lives in seconds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Isotope {
    Neutron,
    Tritium,
    Carbon14,
    Aluminium26,
    Iron60,
    Nickel56,
    Cobalt56,
    Uranium238,
    Potassium40,
}

impl Isotope {
    pub fn half_life(self) -> f64 {
        match self {
            Isotope::Neutron => 878.4,
            Isotope::Tritium => 3.888e8,
            Isotope::Carbon14 => 1.808e11,
            Isotope::Aluminium26 => 2.26e13,
            Isotope::Iron60 => 8.2e13,
            Isotope::Nickel56 => 5.27e5,
            Isotope::Cobalt56 => 6.67e6,
            Isotope::Uranium238 => 1.41e17,
            Isotope::Potassium40 => 3.938e16,
        }
    }

    /// Energy released per decay, J.
    pub fn q_value(self) -> f64 {
        match self {
            Isotope::Neutron => 0.782 * MEV,
            Isotope::Tritium => 0.0186 * MEV,
            Isotope::Carbon14 => 0.156 * MEV,
            Isotope::Aluminium26 => 4.004 * MEV,
            Isotope::Iron60 => 3.05 * MEV,
            Isotope::Nickel56 => 2.136 * MEV,
            Isotope::Cobalt56 => 4.566 * MEV,
            Isotope::Uranium238 => 4.27 * MEV,
            Isotope::Potassium40 => 1.311 * MEV,
        }
    }

    #[inline]
    pub fn decay_constant(self) -> f64 {
        std::f64::consts::LN_2 / self.half_life()
    }

    /// Expected decays from `n` nuclei in time `dt`, sampled as a Poisson
    /// process. For large `n` this is the deterministic exponential law; for
    /// small `n` — the case when a user is watching a handful of atoms — it is
    /// genuinely stochastic, which is the physically correct behaviour and the
    /// thing a deterministic rate equation gets visibly wrong.
    pub fn sample_decays(self, n: f64, dt: f64, stream: &mut Stream) -> f64 {
        if n <= 0.0 || dt <= 0.0 {
            return 0.0;
        }
        let p = 1.0 - (-self.decay_constant() * dt).exp();
        let expected = n * p;
        if n > 1e6 {
            expected
        } else {
            stream.poisson(expected) as f64
        }
    }

    /// Time until the next decay of a single nucleus. This is what the engine
    /// samples when the user is watching one atom — and it is the cleanest
    /// example of the whole architecture: the answer does not exist until
    /// asked for, is drawn from the correct distribution, and once drawn is a
    /// permanent fact recorded in the ledger.
    pub fn sample_lifetime(self, stream: &mut Stream) -> f64 {
        stream.exponential(self.decay_constant())
    }
}

/// Thomson scattering cross section — the opacity floor for ionised gas.
pub const SIGMA_THOMSON: f64 = 6.652_458_7321e-29;

/// Electron-scattering opacity, m^2/kg.
pub fn electron_opacity(comp: Composition) -> f64 {
    0.02 * (1.0 + comp.get(Species::Hydrogen)) * 10.0
}

/// Kramers' bound-free/free-free opacity, m^2/kg.
pub fn kramers_opacity(rho: f64, temperature: f64, metallicity: f64) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    4.34e21 * (1.0 + metallicity * 10.0) * rho * temperature.powf(-3.5) * 0.1
}

/// Nuclear radius from the mass number: `r = 1.2 A^(1/3)` fm.
#[inline]
pub fn nuclear_radius(a: f64) -> f64 {
    1.2e-15 * a.max(1.0).cbrt()
}

/// Coulomb barrier between two nuclei, J.
pub fn coulomb_barrier(z1: f64, z2: f64, a1: f64, a2: f64) -> f64 {
    let r = nuclear_radius(a1) + nuclear_radius(a2);
    K_COULOMB * z1 * z2 * E_CHARGE * E_CHARGE / r
}

/// Semi-empirical mass formula binding energy, MeV — used when the engine needs
/// a nucleus that is not in the lumped species table.
pub fn semf_binding(z: f64, a: f64) -> f64 {
    if a <= 0.0 {
        return 0.0;
    }
    let n = a - z;
    let a_v = 15.75;
    let a_s = 17.8;
    let a_c = 0.711;
    let a_a = 23.7;
    let delta = if (z as i64) % 2 == 0 && (n as i64) % 2 == 0 {
        11.18 / a.sqrt()
    } else if (z as i64) % 2 == 1 && (n as i64) % 2 == 1 {
        -11.18 / a.sqrt()
    } else {
        0.0
    };
    a_v * a - a_s * a.powf(2.0 / 3.0) - a_c * z * (z - 1.0) / a.cbrt()
        - a_a * (n - z) * (n - z) / a
        + delta
}
