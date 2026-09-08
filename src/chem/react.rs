//! Chemistry that happens while the world is running.
//!
//! # The shape this borrows
//!
//! `solvers::nuclear::burn` takes a composition, a density, a temperature and a
//! span, and gives back a changed composition and the energy it released. It
//! costs O(substances) rather than O(particles) and runs on a node's matter,
//! which is what lets a star burn while nobody is looking at its interior.
//!
//! This is the same shape one tier down. A node's [`Mixture`] says what it is
//! made of; its temperature says what state each of those should be in and how
//! much of one will dissolve in another; and one pass moves mass between pools
//! and books the heat. Salt dissolving in water is not a special case written
//! into the engine — it is what this pass does to a mixture that happens to
//! contain a soluble lattice and a liquid.
//!
//! # What it may and may not change
//!
//! It moves mass between *phases* of substances and nothing else. That single
//! summarising is what makes it safe to run on every node every frame: the
//! elemental account cannot move, because dissolving a salt does not transmute
//! anything, so `Matter::composition`, the baryon number and the mass are
//! all untouched by construction rather than by care.
//!
//! What does move is energy, and it is booked. Melting absorbs, freezing
//! releases, boiling absorbs, dissolving does a little of either.
//!
//! # Rates, not equilibria
//!
//! A pass does not jump to the answer. Each process relaxes toward its
//! equilibrium with a time constant, because a spoonful of salt does not
//! dissolve instantly and a lake does not freeze the moment it drops below
//! zero. The constant is the node's own mixing time, which the engine already
//! computes for a completely different reason — how long before its detail
//! stops meaning anything — and which is the right answer here too: it is how
//! long the node takes to stir itself.
//!
//! # Where the numbers come from
//!
//! Nothing here is tabulated per substance. Melting and boiling points, the
//! solubility of one substance in another, and the latent heats all come out of
//! `analyse`, which derives them from the arrangement. Latent heats use the two
//! classic entropy rules — Trouton for vaporisation, Richard for fusion — which
//! turn a temperature into an enthalpy with one constant each.

use super::analyse::solubility_in;
use super::registry::{Mixture, Phase, Registry, SubstanceId};
use crate::units::N_AVOGADRO;

/// Entropy of fusion, J/(mol K). Richard's rule: a solid loses about this much
/// order when it melts, whatever it is made of. The fusion counterpart of
/// Trouton's constant, and it turns a melting point into a latent heat.
pub const RICHARD: f64 = 8.3;
/// Richard's rule assumes a solid that is not associated. Ice is held together
/// by a hydrogen-bond network and is far more ordered than that, so it loses
/// much more on melting: water needs 334 kJ/kg where the plain rule predicts
/// 126. The same exception as Trouton's, for the same reason.
pub const RICHARD_ASSOCIATED: f64 = 22.0;
/// Entropy of vaporisation, J/(mol K).
pub const TROUTON: f64 = 88.0;
/// And its own associated-liquid exception.
pub const TROUTON_ASSOCIATED: f64 = 110.0;

/// What one pass of chemistry did to a node.
///
/// Every field is a mass *fraction* of the node except `heat`, which is joules
/// per kilogram of node. The caller multiplies by the node's mass, because this
/// pass never sees it — which is what keeps it usable on a gram and on a
/// planet.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ReactionReport {
    pub melted: f64,
    pub frozen: f64,
    pub boiled: f64,
    pub condensed: f64,
    pub dissolved: f64,
    pub precipitated: f64,
    /// Net energy absorbed by the matter, J/kg. Positive cools the node.
    pub heat: f64,
    /// Pools that had no substance behind them in the registry, so nothing
    /// could be decided about them. Should be zero.
    pub unresolved: u32,
}

impl ReactionReport {
    /// Whether anything happened at all.
    pub fn quiet(&self) -> bool {
        self.melted == 0.0
            && self.frozen == 0.0
            && self.boiled == 0.0
            && self.condensed == 0.0
            && self.dissolved == 0.0
            && self.precipitated == 0.0
    }

    pub fn moved(&self) -> f64 {
        self.melted + self.frozen + self.boiled + self.condensed + self.dissolved + self.precipitated
    }
}

/// Latent heat of fusion, J/kg, from Richard's rule.
pub fn heat_of_fusion(props: &super::Properties) -> f64 {
    if props.molar_mass <= 0.0 {
        return 0.0;
    }
    let entropy = if props.hydrogen_bonds > 0 { RICHARD_ASSOCIATED } else { RICHARD };
    entropy * props.melting_point / props.molar_mass
}

/// Latent heat of vaporisation, J/kg, from Trouton's rule.
pub fn heat_of_vaporisation(props: &super::Properties) -> f64 {
    if props.molar_mass <= 0.0 {
        return 0.0;
    }
    let entropy = if props.hydrogen_bonds > 0 { TROUTON_ASSOCIATED } else { TROUTON };
    entropy * props.boiling_point / props.molar_mass
}

/// The phase a pure substance would be in at this temperature.
fn phase_at(props: &super::Properties, temperature: f64) -> Phase {
    if temperature >= props.boiling_point {
        Phase::Gas
    } else if temperature >= props.melting_point {
        Phase::Liquid
    } else {
        Phase::Solid
    }
}

/// How much of the node's mass is liquid and can act as a solvent, and which
/// substance dominates it.
fn dominant_solvent(mix: &Mixture) -> (Option<SubstanceId>, f64) {
    let mut best: Option<SubstanceId> = None;
    let mut best_f = 0.0;
    let mut total = 0.0;
    for p in mix.entries() {
        if p.phase == Phase::Liquid {
            total += p.fraction;
            if p.fraction > best_f {
                best_f = p.fraction;
                best = Some(p.substance);
            }
        }
    }
    (best, total)
}

/// Run one pass of chemistry over a node's mixture.
///
/// `tau` is how long the node takes to stir itself — its mixing time. A pass
/// covering much longer than that reaches equilibrium; a much shorter one
/// barely moves. Pass zero or infinity for "settle immediately", which is what
/// authoring wants and what a frame never does.
///
/// Returns what moved and the net heat, in joules per kilogram of node.
pub fn react(
    mix: &mut Mixture,
    reg: &Registry,
    temperature: f64,
    dt: f64,
    tau: f64,
) -> ReactionReport {
    let mut report = ReactionReport::default();
    if !(dt > 0.0) || !temperature.is_finite() || mix.is_empty() {
        return report;
    }

    // How far toward equilibrium one pass gets. An exponential relaxation, so
    // the answer does not depend on how the span happened to be cut up.
    let approach = if !(tau > 0.0) || !tau.is_finite() {
        1.0
    } else {
        1.0 - (-dt / tau).exp()
    };

    // --- phase changes -----------------------------------------------------
    //
    // Melting, freezing, boiling and condensing, each against the substance's
    // own derived transition temperature. A dissolved pool is left alone: it is
    // already surrounded by solvent and has no lattice left to melt.
    let pools: Vec<(SubstanceId, Phase, f64)> =
        mix.entries().iter().map(|p| (p.substance, p.phase, p.fraction)).collect();
    for (id, phase, fraction) in pools {
        if phase == Phase::Dissolved {
            continue;
        }
        let Some(s) = reg.get(id) else {
            report.unresolved += 1;
            continue;
        };
        let want = phase_at(&s.props, temperature);
        if want == phase {
            continue;
        }
        let amount = fraction * approach;
        let moved = mix.convert(id, phase, want, amount);
        if moved <= 0.0 {
            continue;
        }
        let fusion = heat_of_fusion(&s.props);
        let vapour = heat_of_vaporisation(&s.props);
        // Energy is booked by what the transition *is*, not by which direction
        // the pass happened to walk: solid to gas absorbs both latent heats.
        let (heat, melted, frozen, boiled, condensed) = match (phase, want) {
            (Phase::Solid, Phase::Liquid) => (fusion, moved, 0.0, 0.0, 0.0),
            (Phase::Liquid, Phase::Solid) => (-fusion, 0.0, moved, 0.0, 0.0),
            (Phase::Liquid, Phase::Gas) => (vapour, 0.0, 0.0, moved, 0.0),
            (Phase::Gas, Phase::Liquid) => (-vapour, 0.0, 0.0, 0.0, moved),
            (Phase::Solid, Phase::Gas) => (fusion + vapour, moved, 0.0, moved, 0.0),
            (Phase::Gas, Phase::Solid) => (-(fusion + vapour), 0.0, moved, 0.0, moved),
            _ => (0.0, 0.0, 0.0, 0.0, 0.0),
        };
        report.heat += heat * moved;
        report.melted += melted;
        report.frozen += frozen;
        report.boiled += boiled;
        report.condensed += condensed;
    }

    // --- dissolution and precipitation -------------------------------------
    let (solvent_id, solvent_mass) = dominant_solvent(mix);
    let Some(solvent_id) = solvent_id else {
        // Nothing liquid, so nothing is in solution. Anything still dissolved
        // has lost its solvent and comes out.
        let stranded: Vec<(SubstanceId, f64)> = mix
            .entries()
            .iter()
            .filter(|p| p.phase == Phase::Dissolved)
            .map(|p| (p.substance, p.fraction))
            .collect();
        for (id, fraction) in stranded {
            let out = mix.convert(id, Phase::Dissolved, Phase::Solid, fraction * approach);
            report.precipitated += out;
        }
        return report;
    };
    let Some(solvent) = reg.get(solvent_id) else {
        report.unresolved += 1;
        return report;
    };

    // Each solute once, whatever phases it is spread across.
    let solutes: std::collections::BTreeSet<SubstanceId> = mix
        .entries()
        .iter()
        .filter(|p| p.substance != solvent_id && p.phase != Phase::Gas)
        .map(|p| p.substance)
        .collect();

    for id in &solutes {
        let Some(s) = reg.get(*id) else {
            report.unresolved += 1;
            continue;
        };
        // Saturation: kilograms of solute per kilogram of solvent, times how
        // much solvent there is.
        let ceiling = solubility_in(&s.props, &solvent.props) * solvent_mass;
        let held = mix.pool(*id, Phase::Dissolved);
        // Enthalpy of solution. A lattice must be pulled apart, which costs;
        // the pieces are then solvated, which pays. The balance is small
        // compared with either half, and its sign is what makes some salts
        // cool a glass of water and others warm it.
        let solution_heat = 0.08 * s.props.lattice_binding_ev * crate::units::E_CHARGE
            * N_AVOGADRO
            / s.props.molar_mass.max(1e-30);

        if held < ceiling {
            // Room for more: whatever solid or liquid of it is present goes in.
            let room = ceiling - held;
            for source in [Phase::Solid, Phase::Liquid] {
                let available = mix.pool(*id, source);
                if available <= 0.0 {
                    continue;
                }
                let want = (room - report.dissolved).min(available) * approach;
                if want <= 0.0 {
                    continue;
                }
                let moved = mix.convert(*id, source, Phase::Dissolved, want);
                report.dissolved += moved;
                report.heat += solution_heat * moved;
            }
        } else if held > ceiling {
            // Oversaturated — it was cooled, or the solvent boiled off.
            let excess = (held - ceiling) * approach;
            let back = phase_at(&s.props, temperature);
            let back = if back == Phase::Gas { Phase::Solid } else { back };
            let moved = mix.convert(*id, Phase::Dissolved, back, excess);
            report.precipitated += moved;
            report.heat -= solution_heat * moved;
        }
    }

    mix.compact();
    report
}

/// Settle a mixture to the state it would reach given unlimited time.
///
/// What authoring wants: a beaker of salt water described as "salt and water"
/// should *be* salt water without anyone having to step the world. Equivalent
/// to a `react` over an infinite span, which is what `tau = 0` means.
pub fn equilibrate(mix: &mut Mixture, reg: &Registry, temperature: f64) -> ReactionReport {
    let mut total = ReactionReport::default();
    // A handful of passes, because dissolving changes what is liquid, which
    // changes what can dissolve. It converges in two or three; ten is a bound,
    // not an expectation.
    for _ in 0..10 {
        let r = react(mix, reg, temperature, 1.0, 0.0);
        total.melted += r.melted;
        total.frozen += r.frozen;
        total.boiled += r.boiled;
        total.condensed += r.condensed;
        total.dissolved += r.dissolved;
        total.precipitated += r.precipitated;
        total.heat += r.heat;
        total.unresolved = total.unresolved.max(r.unresolved);
        if r.quiet() {
            break;
        }
    }
    total
}
