//! How fast heat crosses matter: thermal conductivity, derived.
//!
//! # What this closes
//!
//! `docs/PLAY.md` D3 says conduction, diffusion and radiative exchange are
//! three calls to one function with different coefficients, and until this
//! module one of the three existed. There was no thermal conductivity anywhere
//! in the engine — `Material` carries an *electrical* resistivity and a
//! specific heat, `chem::Properties` has no transport quantity at all — so two
//! things touching only radiated at each other, which is wrong by orders of
//! magnitude for a hand on cold metal or a pan on a hob.
//!
//! # The law, per phase — the owner's decision for Phase 5
//!
//! A table of conductivities per material is what axiom one forbids. What is
//! derivable is a law per *phase*, which D17 made possible by putting a real
//! phase on every node's pools:
//!
//! - **Condensed: the Einstein–Cahill–Pohl minimum conductivity.** Heat in an
//!   insulator is carried by vibrations, and in the limit where every
//!   vibration is scattered within a wavelength the conductivity is
//!
//!   ```text
//!     k_min = (pi/6)^(1/3) k_B n^(2/3) (v_l + 2 v_t)
//!   ```
//!
//!   with `n` the number density of *atoms* and the `v` the longitudinal and
//!   two transverse sound speeds. The engine has one sound speed per condensed
//!   phase — `eos::Condensed::sound_speed`, from its bulk modulus — and uses it
//!   for all three modes: a liquid carries no transverse wave at all and a
//!   solid's is about half its longitudinal one, so this is high by up to
//!   half for a solid and is the same approximation Bridgman's liquid formula
//!   makes. It is a *lower bound* for a good crystal and right to within a few
//!   for a glass or a liquid, which is most of what a play space is made of.
//!
//! - **Gas: kinetic theory**, Chapman and Enskog's hard-sphere result:
//!
//!   ```text
//!     eta = (5/16) sqrt(pi m k_B T) / (pi d^2),    k = (5/2) eta c_v
//!   ```
//!
//!   with the molecular diameter `d` measured from the volume one molecule
//!   occupies in the condensed phase, `(m / rho)^(1/3)`. The elementary
//!   mean-free-path form `k = rho c_v v_mean lambda / 3` is the same physics
//!   with a coefficient that is known to be low by a factor of about four,
//!   and this is the version with that factor derived rather than guessed.
//!   Independent of density, which is the classic result: a gas conducts as
//!   well at a tenth of an atmosphere as at one, until the mean free path
//!   reaches the size of the container.
//!
//! **No electronic term.** A metal's heat is carried mostly by its electrons,
//! and Wiedemann–Franz would derive that from the electrical resistivity — but
//! the engine does not derive resistivity (`material.rs` says so in its own
//! table), and a law applied to a stated number is the stated number wearing a
//! law's clothes. So a metal is priced as an insulator of its own stiffness,
//! and is low by an order.

use crate::chem::{Mixture, Phase, Properties, Registry};
use crate::eos::Condensed;
use crate::units::{K_B, N_AVOGADRO};

/// Thermal conductivity of one substance in one phase, W/m/K.
pub fn conductivity(props: &Properties, phase: Phase, temperature: f64) -> f64 {
    match phase {
        Phase::Solid | Phase::Liquid => match Condensed::of(props, phase) {
            Some(c) => minimum_conductivity(props, &c),
            None => 0.0,
        },
        _ => gas_conductivity(props, temperature),
    }
}

/// Einstein–Cahill–Pohl, with one sound speed for all three modes.
pub fn minimum_conductivity(props: &Properties, c: &Condensed) -> f64 {
    let atoms = props.atoms_per_unit.max(1) as f64;
    if !(props.unit_mass > 0.0) {
        return 0.0;
    }
    let n = c.rest_density / props.unit_mass * atoms;
    let v = c.sound_speed(c.rest_density);
    (std::f64::consts::PI / 6.0).cbrt() * K_B * n.powf(2.0 / 3.0) * 3.0 * v
}

/// Chapman–Enskog, hard spheres of the size the molecule occupies condensed.
pub fn gas_conductivity(props: &Properties, temperature: f64) -> f64 {
    let m = props.unit_mass;
    let rho = props.density;
    if !(m > 0.0) || !(rho > 0.0) || !(temperature > 0.0) {
        return 0.0;
    }
    let d = (m / rho).cbrt();
    let viscosity =
        5.0 / 16.0 * (std::f64::consts::PI * m * K_B * temperature).sqrt() / (std::f64::consts::PI * d * d);
    // Translational heat capacity per kilogram — the same `3/2 k` per particle
    // the engine's thermal energy is written with.
    let c_v = 1.5 * K_B / m;
    2.5 * viscosity * c_v
}

/// Conductivity of a node's matter, W/m/K, or `None` for matter nobody has
/// described — which is not the same answer as "measured, and a poor
/// conductor".
///
/// A mixture's pools are combined by volume in parallel (Wiener's upper
/// bound). Which bound a real mixture sits near depends on how its phases are
/// arranged, and a node does not know that; the upper one is the right one for
/// a connected solid with something in its pores, which is the common case.
pub fn conductivity_of(mix: &Mixture, reg: &Registry, temperature: f64) -> Option<f64> {
    if mix.is_empty() {
        return None;
    }
    let mut volume = 0.0;
    let mut sum = 0.0;
    for pool in mix.entries() {
        let Some(s) = reg.get(pool.substance) else { continue };
        let rho = match Condensed::of(&s.props, pool.phase) {
            Some(c) => c.rest_density,
            // A gas pool's volume is whatever the node gives it; weight it by
            // the ideal-gas volume of its mass at this temperature and one
            // bar, which is what "mostly air" means in a room.
            None => {
                let molar = s.props.molar_mass.max(1e-30);
                1.0e5 * molar / (K_B * N_AVOGADRO * temperature.max(1.0))
            }
        };
        if !(rho > 0.0) {
            continue;
        }
        let v = pool.fraction / rho;
        volume += v;
        sum += v * conductivity(&s.props, pool.phase, temperature);
    }
    (volume > 0.0).then(|| sum / volume)
}

/// The conductance between two touching things, W/K.
///
/// Two lumps of radius `r_a` and `r_b` whose centres are `distance` apart
/// share a face of area `pi r_min^2` while they touch, and heat crosses it
/// through each one's half of the path in series:
/// `G = A / (r_a / k_a + r_b / k_b)`. Zero once they are apart — across a gap
/// the radiative coefficient is the one that applies, and it is already there.
pub fn conductive_conductance(k_a: f64, r_a: f64, k_b: f64, r_b: f64, distance: f64) -> f64 {
    if !(k_a > 0.0) || !(k_b > 0.0) || !(r_a > 0.0) || !(r_b > 0.0) {
        return 0.0;
    }
    if distance > r_a + r_b {
        return 0.0;
    }
    let area = std::f64::consts::PI * r_a.min(r_b).powi(2);
    area / (r_a / k_a + r_b / k_b)
}

/// The conductance between two neighbouring parcels of one medium, W/K.
///
/// Not two lumps meeting at their surfaces. A sampled parcel's `radius` is not
/// the cell it stands for — measured on a block of silicate drawn in 64
/// parcels, neighbours sit 1.35 to 1.6 diameters apart, so "touching" never
/// happens inside a continuum and a test of it found no pair at all. What two
/// neighbouring cells of one medium share is a face of about the square of
/// their spacing across a path of the spacing, which is the finite-volume
/// conductance `k d^2 / d = k d`, with the two half-paths in series.
pub fn continuum_conductance(k_a: f64, k_b: f64, distance: f64) -> f64 {
    if !(k_a > 0.0) || !(k_b > 0.0) || !(distance > 0.0) {
        return 0.0;
    }
    2.0 * k_a * k_b / (k_a + k_b) * distance
}
