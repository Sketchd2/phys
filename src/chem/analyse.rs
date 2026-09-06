//! Deriving a substance's properties from its arrangement.
//!
//! # What "derived" means here, and what it does not
//!
//! The rule this engine works to is that only physical laws sit on the axiom
//! side. Caffeine is not an axiom: it is eight carbons, ten hydrogens, four
//! nitrogens and two oxygens joined in a particular way, and everything about
//! it should follow from that plus the constants of the elements involved.
//!
//! So the axioms here are per *element* — atomic weight, electronegativity,
//! covalent radius, one homonuclear bond energy — which are measured constants
//! of nature in the same sense as the gravitational constant. Everything per
//! *substance* is computed:
//!
//! * **Molar mass** is an exact sum. No model, no error.
//! * **Bond energy** for any pair comes from Pauling's relation, the extra
//!   stability an unequal pair gets from its ionic resonance:
//!   `D(A–B) = ½[D(A–A) + D(B–B)] + 96.5 (χA − χB)²` kJ/mol. One number per
//!   element gives every pair, which is the difference between a table that
//!   describes what somebody wrote down and one that describes every pair there
//!   is.
//! * **Bond length** is the sum of covalent radii with the Schomaker–Stevenson
//!   correction `− 9|χA − χB|` pm, which is why H–F is shorter than the radii
//!   alone would say.
//! * **Ionic character** is Pauling's `1 − exp(−¼(χA − χB)²)`. This is the one
//!   that matters most at play scale: it is what tells the engine that sodium
//!   chloride is a salt and methane is not, without either being named.
//! * **Cohesive energy** is the sum over bonds, which for a lattice is per
//!   formula unit and for a molecule is the energy to take it apart.
//! * **Polarity** follows from the bond dipoles and the geometry, and with
//!   ionic character it decides what dissolves in what.
//!
//! # Where it stops
//!
//! Melting and boiling points are *correlations* against cohesive energy, not
//! derivations, and they are marked as such. A protein's function, an enzyme's
//! specificity, a drug's effect: none of that is here, and the honest position
//! is that those are tabulated responses until they can be solved for.
//! [`Confidence`] is how a caller tells which it got, and
//! [`super::registry::Provenance`] records it against the substance for good,
//! so a later build that *can* derive one of these can find every substance
//! whose value came from a table and recompute it.

use super::arrange::{Arrangement, Lattice, Order};
use super::Element;
use crate::units::N_AVOGADRO;

/// kJ/mol per unit of squared electronegativity difference, in Pauling's
/// relation. His original constant, in his original units.
const PAULING_IONIC_KJ: f64 = 96.5;
/// Schomaker–Stevenson shortening, picometres per unit of electronegativity
/// difference.
const SS_SHORTENING_PM: f64 = 9.0;

/// How much weight to put on a derived number.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Confidence {
    /// Exact given the inputs: a sum of atomic masses, a count of atoms.
    Exact,
    /// A physical relation with measured element constants behind it.
    /// Pauling's bond energy, Schomaker–Stevenson lengths, ionic character.
    Derived,
    /// A correlation that reproduces the trend and not the number. Melting
    /// points, densities.
    Correlated,
    /// The arrangement contains an element this build has no data for, so the
    /// result leans on defaults. Anything at this level should be treated as
    /// an order of magnitude.
    Guessed,
}

/// Why an arrangement cannot be analysed.
///
/// A refusal, not a clamp. An arrangement that breaks valence is not a
/// substance with unusual properties, it is a description of something that
/// does not hold together, and giving it properties anyway would let it into
/// the world.
#[derive(Debug, Clone, PartialEq)]
pub enum Illegal {
    /// No atoms at all.
    Empty,
    /// An atom holds more bonds than it has valence slots. The check that
    /// stops a hydrogen acquiring five neighbours.
    OverBonded { atom: usize, element: Element, used: usize, allowed: usize },
    /// A bond names an atom that is not there.
    DanglingBond { bond: usize },
    /// A molecule in two or more pieces, which is two substances rather than
    /// one. Lattices are exempt: their connectivity is in the tiling.
    Disconnected { pieces: usize },
    /// An element outside the table, in a position where guessing would be
    /// dishonest.
    UnknownElement(Element),
}

impl std::fmt::Display for Illegal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Illegal::Empty => write!(f, "no atoms"),
            Illegal::OverBonded { atom, element, used, allowed } => write!(
                f,
                "atom {atom} ({element}) holds {used} bonds but has {allowed} valence slots"
            ),
            Illegal::DanglingBond { bond } => write!(f, "bond {bond} names an atom that is not there"),
            Illegal::Disconnected { pieces } => {
                write!(f, "a molecule in {pieces} pieces is {pieces} substances")
            }
            Illegal::UnknownElement(e) => write!(f, "no data for element Z={}", e.z()),
        }
    }
}

/// What one bond turned out to be.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondFacts {
    /// Dissociation energy, joules.
    pub energy: f64,
    /// Equilibrium length, metres.
    pub length: f64,
    /// Fraction of the bond that is charge transfer rather than sharing, 0..1.
    pub ionicity: f64,
}

/// Everything the engine needs to know about a substance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Properties {
    /// Mass of one formula unit, kg. Exact.
    pub unit_mass: f64,
    /// Molar mass, kg/mol. Exact.
    pub molar_mass: f64,
    /// Energy to take one formula unit apart, joules. Positive.
    pub cohesive_energy: f64,
    /// Mean fraction of charge transfer across the bonds, 0..1. Above about
    /// 0.5 the substance behaves as a salt.
    pub ionicity: f64,
    /// Where this substance sits on the "like dissolves like" axis: how much
    /// charge it separates, in arbitrary but consistent units.
    ///
    /// Computed here rather than by [`solubility_in`], because it needs to know
    /// whether the substance is a lattice. An ionic crystal's formula unit has
    /// an enormous nominal dipole — a sodium ion beside a chloride is 8.5
    /// debye — and treating that as if it were a molecular dipole made salt
    /// come out *less* soluble in water than methane.
    pub polarity: f64,
    /// Net dipole, in coulomb-metres, from the bond dipoles and the geometry.
    /// Zero for anything symmetric, which is why carbon dioxide is non-polar
    /// despite two strongly polar bonds.
    pub dipole: f64,
    /// Hydrogen-bond donors per formula unit: hydrogens bonded to nitrogen,
    /// oxygen or fluorine.
    ///
    /// Carried on the substance because two entropy rules need it. Trouton's
    /// and Richard's constants both assume a liquid and a solid that are not
    /// associated, and hydrogen-bonded substances are more ordered than either
    /// assumes — which is why water boils at 373 K rather than 462 and needs
    /// 334 kJ/kg to melt rather than 126.
    pub hydrogen_bonds: u8,
    /// Energy per atom holding this substance's crystal together, eV. Zero for
    /// a molecular substance, which has no lattice to break before it
    /// dissolves.
    pub lattice_binding_ev: f64,
    /// Estimated density, kg/m^3.
    pub density: f64,
    /// Estimated melting point, K.
    pub melting_point: f64,
    /// Estimated boiling point, K.
    pub boiling_point: f64,
    /// Solubility in water, kg per kg of water at room temperature. Derived
    /// from ionicity and polarity by "like dissolves like", which is a rule
    /// about the *relation* between two substances rather than a property of
    /// one — see [`solubility_in`].
    pub water_solubility: f64,
    /// The weakest link in the numbers above.
    pub confidence: Confidence,
}

/// Check an arrangement holds together, and say why if it does not.
pub fn legality(a: &Arrangement) -> Result<(), Illegal> {
    if a.atoms.is_empty() {
        return Err(Illegal::Empty);
    }
    for (i, b) in a.bonds.iter().enumerate() {
        if b.a as usize >= a.atoms.len() || b.b as usize >= a.atoms.len() {
            return Err(Illegal::DanglingBond { bond: i });
        }
    }
    let used = a.used_slots();
    for (i, e) in a.atoms.iter().enumerate() {
        let allowed = match e.valence() {
            Some(v) => v,
            None => return Err(Illegal::UnknownElement(*e)),
        };
        // An ion has given an electron away or taken one on, so its capacity
        // moves with its charge — which is exactly how a sodium ion in salt
        // holds a bond that neutral sodium's single valence slot would already
        // have spent.
        let allowed = if a.charge != 0 { allowed + a.charge.unsigned_abs() as usize } else { allowed };
        if used[i] > allowed {
            return Err(Illegal::OverBonded { atom: i, element: *e, used: used[i], allowed });
        }
    }
    if a.lattice == Lattice::Molecular && a.atoms.len() > 1 {
        let pieces = components(a);
        if pieces > 1 {
            return Err(Illegal::Disconnected { pieces });
        }
    }
    Ok(())
}

fn components(a: &Arrangement) -> usize {
    let adj = a.neighbours();
    let mut seen = vec![false; a.atoms.len()];
    let mut pieces = 0;
    for start in 0..a.atoms.len() {
        if seen[start] {
            continue;
        }
        pieces += 1;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(i) = stack.pop() {
            for (j, _) in &adj[i] {
                if !seen[*j] {
                    seen[*j] = true;
                    stack.push(*j);
                }
            }
        }
    }
    pieces
}

/// Fraction of a bond that is charge transfer rather than sharing.
///
/// Pauling's `1 − exp(−¼Δχ²)`. Zero for a homonuclear bond by construction,
/// and it is what separates a salt from a molecule without either being named.
pub fn ionicity(a: Element, b: Element) -> f64 {
    match (a.electronegativity(), b.electronegativity()) {
        (Some(x), Some(y)) if x > 0.0 && y > 0.0 => {
            let d = x - y;
            1.0 - (-0.25 * d * d).exp()
        }
        _ => 0.0,
    }
}

/// Dissociation energy of one bond, joules.
///
/// Pauling's relation: the geometric intuition is that an unequal pair is more
/// stable than the average of the two equal pairs, by an amount that grows with
/// the square of how unequal they are.
pub fn bond_energy(a: Element, b: Element, order: Order) -> Option<f64> {
    let (da, db) = (a.homonuclear_bond()?, b.homonuclear_bond()?);
    let (xa, xb) = (a.electronegativity()?, b.electronegativity()?);
    let extra = if xa > 0.0 && xb > 0.0 {
        PAULING_IONIC_KJ * 1e3 / N_AVOGADRO * (xa - xb).powi(2)
    } else {
        0.0
    };
    Some((0.5 * (da + db) + extra) * order.strength())
}

/// Equilibrium bond length, metres. Covalent radii, shortened by
/// Schomaker–Stevenson where the pair is unequal.
pub fn bond_length(a: Element, b: Element, order: Order) -> Option<f64> {
    let (ra, rb) = (a.covalent_radius()?, b.covalent_radius()?);
    let (xa, xb) = (a.electronegativity()?, b.electronegativity()?);
    let shorten = if xa > 0.0 && xb > 0.0 {
        SS_SHORTENING_PM * 1e-12 * (xa - xb).abs()
    } else {
        0.0
    };
    // A double bond is about 13% shorter than a single, a triple about 22%.
    let multiplicity = match order {
        Order::Single | Order::Ionic => 1.0,
        Order::Double => 0.87,
        Order::Triple => 0.78,
        Order::Hydrogen => 2.6,
    };
    Some(((ra + rb) * multiplicity - shorten).max(3e-11))
}

/// Everything about one bond.
pub fn bond_facts(a: Element, b: Element, order: Order) -> Option<BondFacts> {
    Some(BondFacts {
        energy: bond_energy(a, b, order)?,
        length: bond_length(a, b, order)?,
        ionicity: ionicity(a, b),
    })
}

/// Derive a substance's properties from its arrangement.
pub fn analyse(arr: &Arrangement) -> Result<Properties, Illegal> {
    legality(arr)?;
    let formula = arr.formula();
    let unit_mass = formula.mass().ok_or(Illegal::UnknownElement(arr.atoms[0]))?;

    let mut confidence = Confidence::Exact;
    let mut cohesive = 0.0;
    let mut ionic_sum = 0.0;
    let mut dipole_sum = 0.0;
    let mut counted = 0.0;

    for b in &arr.bonds {
        let (ea, eb) = (arr.atoms[b.a as usize], arr.atoms[b.b as usize]);
        match bond_facts(ea, eb, b.order) {
            Some(f) => {
                cohesive += f.energy;
                ionic_sum += f.ionicity;
                // A bond dipole is the transferred charge times the separation.
                dipole_sum += f.ionicity * crate::units::E_CHARGE * f.length;
                counted += 1.0;
                confidence = confidence.max(Confidence::Derived);
            }
            None => {
                confidence = Confidence::Guessed;
            }
        }
    }

    let ionicity = if counted > 0.0 { ionic_sum / counted } else { 0.0 };

    // The dipole needs the *shape*, not the graph. Two bond dipoles of equal
    // size either cancel or add depending on the angle between them, and that
    // is the entire difference between carbon dioxide, which is linear and
    // barely dissolves, and water, which is bent and dissolves everything. An
    // earlier version tried to decide it from which atoms were
    // interchangeable and gave water a dipole of exactly zero.
    let conformer = super::geometry::embed(arr);
    let dipole = super::geometry::dipole(arr, &conformer);
    let _ = dipole_sum;

    // Density from how much space the atoms take up. For a lattice the cell
    // says so directly; for a molecule, from the covalent volume with a packing
    // fraction that reproduces ordinary liquids and solids.
    let (density, density_conf) = match arr.lattice {
        Lattice::Cubic { a } if a > 0.0 => (unit_mass / (a * a * a), Confidence::Derived),
        Lattice::Hexagonal { a, c } if a > 0.0 && c > 0.0 => {
            let volume = 3.0f64.sqrt() * 0.5 * a * a * c;
            (unit_mass / volume, Confidence::Derived)
        }
        _ => {
            // Van der Waals radii, not covalent. How much room an atom takes up
            // when it is *not* bonded is what sets a condensed phase's density,
            // and the two radii differ by a factor of two to three — cubed,
            // that is the eightfold error that had water at 10,700 kg/m^3.
            //
            // Bonded atoms overlap, so the sum of free-atom spheres
            // overestimates; the factor below is what reproduces ordinary
            // liquids across the range this table covers.
            let mut volume = 0.0;
            for e in &arr.atoms {
                let r = e.vdw_radius().unwrap_or(1.7e-10);
                volume += 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
            }
            (unit_mass / (volume * 0.62).max(1e-45), Confidence::Correlated)
        }
    };
    confidence = confidence.max(density_conf);

    // Melting and boiling are set by how much energy holds one unit *in place*,
    // and for a molecular substance that is not its bonds — it is the much
    // weaker forces between whole molecules. Getting that distinction right is
    // why salt melts at 1074 K and sugar at 460 K.
    let (melting_point, boiling_point) = phase_points(arr, &conformer, unit_mass, dipole, cohesive);
    confidence = confidence.max(Confidence::Correlated);

    let polarity = polarity_of(arr, ionicity, dipole);
    // Zero for a molecule: there is no lattice to take apart before it can
    // dissolve, only itself to surround.
    let lattice_binding_ev = match arr.lattice {
        Lattice::Molecular => 0.0,
        _ if arr.atoms.is_empty() => 0.0,
        _ => {
            // Lattice energy grows as the product of the ionic charges, which
            // is why sodium chloride dissolves and uranium dioxide does not
            // despite both being ionic oxide-or-halide lattices. The charges
            // are not in the arrangement — nothing declares an oxidation state
            // — but each element's ordinary valence stands in for one, and it
            // is enough to separate a 1:1 salt from a 4:2 oxide.
            let mut charge_product = 1.0;
            let mut counted = 0.0;
            for b in &arr.bonds {
                let (za, zb) = (
                    arr.atoms[b.a as usize].valence().unwrap_or(1).max(1) as f64,
                    arr.atoms[b.b as usize].valence().unwrap_or(1).max(1) as f64,
                );
                charge_product += za * zb;
                counted += 1.0;
            }
            let mean = if counted > 0.0 { (charge_product - 1.0) / counted } else { 1.0 };
            cohesive / arr.atoms.len() as f64 / crate::units::E_CHARGE * mean.powf(0.4)
        }
    };
    let water_solubility =
        dissolves(polarity, unit_mass * N_AVOGADRO, lattice_binding_ev, WATER_POLARITY);

    // A lattice's "dipole" is an artefact of picking one formula unit out of an
    // infinite alternating solid; what makes it polar is that it is ionic.
    Ok(Properties {
        unit_mass,
        polarity,
        hydrogen_bonds: hydrogen_bond_donors(arr).min(255) as u8,
        lattice_binding_ev,
        molar_mass: unit_mass * N_AVOGADRO,
        cohesive_energy: cohesive,
        ionicity,
        dipole,
        density,
        melting_point,
        boiling_point,
        water_solubility,
        confidence,
    })
}

/// Trouton's constant: entropy of vaporisation, J/(mol K). A liquid whose
/// molecules are not associated loses about this much order when it boils,
/// whatever it is made of, which is what makes a boiling point derivable from
/// a vaporisation enthalpy at all.
const TROUTON: f64 = 88.0;
/// Hydrogen-bonded liquids are more ordered to begin with, so they lose more.
/// This is the standard exception to Trouton's rule, and without it water
/// boils at 462 K.
const TROUTON_ASSOCIATED: f64 = 110.0;
/// Vaporisation enthalpy per gram per mole, kJ/mol, from dispersion alone.
/// Calibrated on the non-associating substances in `tests/chem.rs`; dispersion
/// scales with polarisability, which scales with how many electrons there are,
/// which for ordinary matter tracks molar mass closely enough to use it.
const DISPERSION_KJ_PER_G: f64 = 0.51;
/// One hydrogen bond's worth, kJ/mol.
const HBOND_KJ: f64 = 15.5;

/// How many hydrogen-bond donors an arrangement has: hydrogens bonded to
/// nitrogen, oxygen or fluorine.
///
/// The single largest thing separating substances of similar size, and the
/// reason water is a liquid at room temperature while methane is a gas
/// eighteen times lighter than it has any right to be.
pub fn hydrogen_bond_donors(arr: &Arrangement) -> usize {
    let adj = arr.neighbours();
    let mut n = 0;
    for (i, e) in arr.atoms.iter().enumerate() {
        if *e != Element::HYDROGEN {
            continue;
        }
        for (j, _) in &adj[i] {
            let z = arr.atoms[*j].z();
            if z == 7 || z == 8 || z == 9 {
                n += 1;
            }
        }
    }
    n
}

/// Melting and boiling points, K.
///
/// For a molecular substance: estimate the vaporisation enthalpy from
/// dispersion (scaling with mass), the dipole, and hydrogen bonding, then
/// divide by Trouton's constant. For a lattice there is nothing to vaporise
/// as a unit — the bonds themselves must go — so it scales with the cohesive
/// energy instead.
///
/// The melting point is taken as a fixed fraction of the boiling point. That
/// fraction really does vary from 0.45 to 0.81 across ordinary substances, so
/// it is the weakest number this module produces and is marked accordingly.
fn phase_points(
    arr: &Arrangement,
    _conformer: &[crate::math::Vec3],
    unit_mass: f64,
    dipole: f64,
    cohesive: f64,
) -> (f64, f64) {
    match arr.lattice {
        Lattice::Molecular => {
            let grams_per_mol = unit_mass * N_AVOGADRO * 1e3;
            let donors = hydrogen_bond_donors(arr);
            // Keesom: dipole-dipole attraction, in the same units.
            let debye = dipole / 3.33564e-30;
            let vap_kj =
                DISPERSION_KJ_PER_G * grams_per_mol + HBOND_KJ * donors as f64 + 1.6 * debye * debye;
            let trouton = if donors > 0 { TROUTON_ASSOCIATED } else { TROUTON };
            let boiling = (vap_kj * 1e3 / trouton).clamp(2.0, 4000.0);
            (boiling * 0.6, boiling)
        }
        _ => {
            let per_atom = if arr.atoms.is_empty() {
                0.0
            } else {
                cohesive / arr.atoms.len() as f64
            };
            let melting = (0.028 * per_atom / crate::units::K_B).clamp(1.0, 6000.0);
            (melting, (melting * 1.7).clamp(2.0, 8000.0))
        }
    }
}

/// Where a substance sits on the "like dissolves like" axis.
///
/// Charge separation of every kind that a solvent can grip: ionic character,
/// the net dipole, and hydrogen bonding, which is the one that cannot be left
/// out. Ethanol and propane have comparable dipoles and ethanol is miscible
/// with water while propane is not, and the whole of that difference is one
/// hydroxyl group.
///
/// A lattice takes the ionic term alone. Its formula unit has an enormous
/// nominal dipole — a sodium ion beside a chloride is 8.5 debye — but that is
/// an artefact of cutting one unit out of an infinite alternating solid, and
/// treating it as a molecular dipole made salt come out less soluble in water
/// than methane.
pub fn polarity_of(arr: &Arrangement, ionicity: f64, dipole: f64) -> f64 {
    match arr.lattice {
        Lattice::Molecular => {
            let donors = hydrogen_bond_donors(arr) as f64;
            let acceptors = arr
                .atoms
                .iter()
                .filter(|e| e.z() == 7 || e.z() == 8 || e.z() == 9)
                .count() as f64;
            ionicity * 6.0 + (dipole / 3.33564e-30).min(4.0) + 1.5 * donors + 0.6 * acceptors
        }
        _ => ionicity * 9.0,
    }
}

/// Water's own position on that axis.
///
/// A constant rather than a call, because every solubility is measured against
/// water and analysing water to find out would be circular. `water_is_where_we_
/// think_it_is` in `tests/chem.rs` checks the two agree, so the constant cannot
/// drift away from the model that produced it.
pub const WATER_POLARITY: f64 = 7.12;

/// How much of one substance dissolves in another, kg per kg, at room
/// temperature.
///
/// # An order of magnitude, and it says so
///
/// Aqueous solubility spans about ten decades between methane and sodium
/// chloride, so this is a power of ten and is worth about a decade and a half
/// either way. That is enough for the questions the engine asks — does this
/// dissolve, roughly how much before it saturates — and not enough to design
/// anything with. `Confidence::Correlated` says so, and `tests/chem.rs`
/// measures the residuals against real values rather than claiming a fit.
///
/// Three terms. How unalike the two are on the polarity axis, which is
/// "like dissolves like" itself. How large the solute is, since a big molecule
/// costs more cavity than its polar groups can pay for — the reason methanol is
/// miscible and octanol is not. And, for a crystal, how strongly it is bound,
/// because a lattice has to be taken apart before any of it can be surrounded.
///
/// # Where the number stops being a number
///
/// It is calibrated on small molecules and singly-charged salts in water, which
/// is the case that matters at play scale, and it is good to a decade and a
/// half there. Outside it, read the answer as a *class* rather than a
/// quantity:
///
/// * Multiply-charged lattices — metal oxides, phosphates — have lattice
///   energies growing as the product of the ionic charges, and the valence
///   proxy for that is coarse. Uranium dioxide comes out at 10^-9 against a
///   real 10^-9, and that agreement is luckier than the model deserves.
/// * Network solids like quartz are limited by how fast they dissolve rather
///   than by whether they can, and nothing here models a rate.
///
/// In both cases the answer to take from it is "negligible", which is the
/// answer the engine needs, and not the digits.
pub fn dissolves(
    solute_polarity: f64,
    solute_molar_mass: f64,
    solute_lattice_ev_per_atom: f64,
    solvent_polarity: f64,
) -> f64 {
    let log10 = SOLUBILITY_INTERCEPT
        - 1.43 * (solute_polarity - solvent_polarity).abs()
        - 0.05 * (solute_molar_mass / 0.018)
        - 1.2 * solute_lattice_ev_per_atom;
    // Past a few kilograms per kilogram the distinction stops meaning anything:
    // the substance is miscible.
    10f64.powf(log10.min(1.0))
}

/// Fitted on methane, carbon dioxide, ethanol and sodium chloride — four
/// substances spanning five decades of solubility and every mechanism the axis
/// above describes.
const SOLUBILITY_INTERCEPT: f64 = 5.31;

/// How much of one substance dissolves in another, kg per kg.
///
/// The general form of [`Properties::water_solubility`], and the reason that
/// field is only a convenience: solubility is a property of a *pair*, not of a
/// substance, and an engine that stored it per substance could never answer
/// "does this dissolve in ethanol".
pub fn solubility_in(solute: &Properties, solvent: &Properties) -> f64 {
    dissolves(
        solute.polarity,
        solute.molar_mass,
        solute.lattice_binding_ev,
        solvent.polarity,
    )
}
