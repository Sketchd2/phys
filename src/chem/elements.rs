//! The elements, as data rather than as a lumped bucket.
//!
//! # Why this is not `units::Species`
//!
//! [`Species`](crate::units::Species) is eight buckets — H, He, C, N, O, Si,
//! Fe and everything-else — and for what it does that is not a compromise. It
//! is the account nuclear burning runs on, and burning cares about hydrogen,
//! helium, the CNO catalysts, silicon and the iron floor. Astrophysics uses
//! very nearly that list for the same reason.
//!
//! It cannot describe salt. Sodium and chlorine both fall into `Other`, whose
//! mass number is 65 because it stands in for zinc-ish heavy nuclei, so a
//! kilogram of sodium counted through that bucket contains 2.8 times too few
//! atoms. That is harmless at stellar tier, where nobody counts sodium atoms,
//! and useless at the bench.
//!
//! So the two coexist. `Species` stays the *elemental and nuclear* account,
//! carried by every aggregate and every body, conserved through prolongation
//! and restriction, unchanged by this module. [`Element`] is the *chemical*
//! account, used where arrangements are analysed, and it names real elements by
//! atomic number. A substance's formula is in elements; its contribution to a
//! node's bulk composition is that formula lumped back into species, which is
//! how the two stay reconcilable — see `chem::registry`, which tests it.
//!
//! # The table
//!
//! Standard atomic weights (CIAAW), Pauling electronegativities, and Cordero
//! covalent radii for Z = 1..36, which spans everything ordinary matter at
//! human, animal and vehicle scale is made of. Past krypton the table stops and
//! [`Element::known`] says so rather than inventing values: an estimated
//! electronegativity would propagate silently into every bond it touched.

/// An element, by atomic number.
///
/// A newtype rather than an enum on purpose. The set of elements is not a thing
/// this engine gets to choose, and an enum would make "the elements we happen
/// to have written down" into a type — which is exactly the mistake `Species`
/// makes deliberately and this module exists to avoid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Element(pub u8);

struct Entry {
    symbol: &'static str,
    /// Standard atomic weight, in unified atomic mass units.
    weight: f64,
    /// Pauling electronegativity. Zero where it is undefined (noble gases),
    /// which callers must treat as "no opinion" rather than "very electropositive".
    electronegativity: f64,
    /// Cordero single-bond covalent radius, picometres. How close two atoms
    /// sit when they are *bonded*.
    radius_pm: f64,
    /// Bondi van der Waals radius, picometres. How much room an atom takes up
    /// when it is *not* bonded, which is what sets the density of a condensed
    /// phase.
    ///
    /// These are two to three times the covalent radius, and volume goes as the
    /// cube, so using the wrong one is an eightfold error in density — which is
    /// exactly what it was until water came out at 10,700 kg/m^3.
    vdw_pm: f64,
    /// Valence electrons, for counting lone pairs. Group number for the main
    /// group; the transition metals are given the count that reproduces their
    /// usual coordination rather than their full d shell.
    valence_electrons: u8,
    /// Bonds it will ordinarily hold.
    valence: u8,
    /// Homonuclear single-bond dissociation energy, kJ/mol.
    ///
    /// One measured number per element, which with electronegativity is enough
    /// to give *every* pair a bond energy through Pauling's relation — see
    /// `analyse::bond_energy`. That is the difference between a table of bonds,
    /// which can only describe pairs somebody wrote down, and a table of
    /// elements, which describes every pair there is.
    ///
    /// Zero where a single covalent bond is not a meaningful object: the noble
    /// gases, and the elements whose condensed phase is metallic rather than
    /// covalent. Transition-metal entries are the least reliable numbers in
    /// this file and are marked as such by `Analysis::confidence`.
    bond_kj: f64,
}

#[allow(clippy::too_many_arguments)]
const fn e(
    symbol: &'static str,
    weight: f64,
    electronegativity: f64,
    radius_pm: f64,
    vdw_pm: f64,
    valence: u8,
    valence_electrons: u8,
    bond_kj: f64,
) -> Entry {
    Entry { symbol, weight, electronegativity, radius_pm, vdw_pm, valence, valence_electrons, bond_kj }
}

/// Indexed by atomic number; slot zero is a placeholder so `Z` indexes directly.
const TABLE: [Entry; 37] = [
    e("?", 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0),
    e("H", 1.008, 2.2, 31.0, 120.0, 1, 1, 436.0),
    e("He", 4.0026, 0.0, 28.0, 140.0, 0, 2, 0.0),
    e("Li", 6.94, 0.98, 128.0, 182.0, 1, 1, 105.0),
    e("Be", 9.0122, 1.57, 96.0, 153.0, 2, 2, 208.0),
    e("B", 10.81, 2.04, 84.0, 192.0, 3, 3, 293.0),
    e("C", 12.011, 2.55, 76.0, 170.0, 4, 4, 346.0),
    e("N", 14.007, 3.04, 71.0, 155.0, 3, 5, 167.0),
    e("O", 15.999, 3.44, 66.0, 152.0, 2, 6, 142.0),
    e("F", 18.998, 3.98, 57.0, 147.0, 1, 7, 155.0),
    e("Ne", 20.18, 0.0, 58.0, 154.0, 0, 8, 0.0),
    e("Na", 22.99, 0.93, 166.0, 227.0, 1, 1, 72.0),
    e("Mg", 24.305, 1.31, 141.0, 173.0, 2, 2, 129.0),
    e("Al", 26.982, 1.61, 121.0, 184.0, 3, 3, 167.0),
    e("Si", 28.085, 1.9, 111.0, 210.0, 4, 4, 222.0),
    e("P", 30.974, 2.19, 107.0, 180.0, 3, 5, 201.0),
    e("S", 32.06, 2.58, 105.0, 180.0, 2, 6, 226.0),
    e("Cl", 35.45, 3.16, 102.0, 175.0, 1, 7, 242.0),
    e("Ar", 39.95, 0.0, 106.0, 188.0, 0, 8, 0.0),
    e("K", 39.098, 0.82, 203.0, 275.0, 1, 1, 49.0),
    e("Ca", 40.078, 1.0, 176.0, 231.0, 2, 2, 105.0),
    e("Sc", 44.956, 1.36, 170.0, 211.0, 3, 3, 138.0),
    e("Ti", 47.867, 1.54, 160.0, 187.0, 4, 4, 118.0),
    e("V", 50.942, 1.63, 153.0, 179.0, 5, 5, 269.0),
    e("Cr", 51.996, 1.66, 139.0, 189.0, 3, 3, 155.0),
    e("Mn", 54.938, 1.55, 139.0, 197.0, 2, 2, 116.0),
    e("Fe", 55.845, 1.83, 132.0, 194.0, 3, 3, 118.0),
    e("Co", 58.933, 1.88, 126.0, 192.0, 2, 2, 167.0),
    e("Ni", 58.693, 1.91, 124.0, 163.0, 2, 2, 204.0),
    e("Cu", 63.546, 1.9, 132.0, 140.0, 2, 2, 177.0),
    e("Zn", 65.38, 1.65, 122.0, 139.0, 2, 2, 29.0),
    e("Ga", 69.723, 1.81, 122.0, 187.0, 3, 3, 113.0),
    e("Ge", 72.63, 2.01, 120.0, 211.0, 4, 4, 188.0),
    e("As", 74.922, 2.18, 119.0, 185.0, 3, 5, 146.0),
    e("Se", 78.971, 2.55, 120.0, 190.0, 2, 6, 172.0),
    e("Br", 79.904, 2.96, 120.0, 185.0, 1, 7, 193.0),
    e("Kr", 83.798, 3.0, 116.0, 202.0, 0, 8, 0.0),
];

/// Highest atomic number the table covers.
pub const HEAVIEST: u8 = 36;

impl Element {
    pub const HYDROGEN: Element = Element(1);
    pub const CARBON: Element = Element(6);
    pub const NITROGEN: Element = Element(7);
    pub const OXYGEN: Element = Element(8);
    pub const SODIUM: Element = Element(11);
    pub const SILICON: Element = Element(14);
    pub const SULFUR: Element = Element(16);
    pub const CHLORINE: Element = Element(17);
    pub const CALCIUM: Element = Element(20);
    pub const IRON: Element = Element(26);

    /// Atomic number.
    pub fn z(self) -> u8 {
        self.0
    }

    /// Whether this build has data for it. Everything else returns `None`
    /// rather than an estimate, because an invented electronegativity would
    /// propagate silently into every bond it touched.
    pub fn known(self) -> bool {
        self.0 >= 1 && self.0 <= HEAVIEST
    }

    fn entry(self) -> Option<&'static Entry> {
        if self.known() {
            Some(&TABLE[self.0 as usize])
        } else {
            None
        }
    }

    pub fn symbol(self) -> &'static str {
        self.entry().map(|e| e.symbol).unwrap_or("?")
    }

    /// Standard atomic weight, unified atomic mass units.
    pub fn weight(self) -> Option<f64> {
        self.entry().map(|e| e.weight)
    }

    /// Mass of one atom, kg.
    pub fn mass_kg(self) -> Option<f64> {
        self.weight().map(|w| w * crate::units::AMU)
    }

    /// Pauling electronegativity, or `None` where it is undefined. Note that a
    /// noble gas reports `Some(0.0)`: the value is defined as "does not apply",
    /// not missing, and lumping the two together would make argon look like the
    /// most electropositive element there is.
    pub fn electronegativity(self) -> Option<f64> {
        self.entry().map(|e| e.electronegativity)
    }

    /// Single-bond covalent radius, metres. Use for bond lengths.
    pub fn covalent_radius(self) -> Option<f64> {
        self.entry().map(|e| e.radius_pm * 1e-12)
    }

    /// Van der Waals radius, metres. Use for how much room an atom occupies,
    /// which is what sets a condensed phase's density.
    pub fn vdw_radius(self) -> Option<f64> {
        self.entry().map(|e| e.vdw_pm * 1e-12)
    }

    /// Valence electrons, for counting lone pairs — which is what makes water
    /// bent rather than linear, and therefore polar rather than not.
    pub fn valence_electrons(self) -> Option<usize> {
        self.entry().map(|e| e.valence_electrons as usize)
    }

    /// How many bonds it will ordinarily hold.
    pub fn valence(self) -> Option<usize> {
        self.entry().map(|e| e.valence as usize)
    }

    /// Homonuclear single-bond dissociation energy, joules per bond.
    ///
    /// `None` for an element outside the table; `Some(0.0)` where a single
    /// covalent bond is not a meaningful object, which the caller must treat
    /// as "this element does not bond this way" rather than "a free bond".
    pub fn homonuclear_bond(self) -> Option<f64> {
        self.entry()
            .map(|e| e.bond_kj * 1e3 / crate::units::N_AVOGADRO)
    }

    /// Parse a symbol. Case-sensitive, as chemistry is: `Co` is cobalt and
    /// `CO` is carbon monoxide.
    pub fn from_symbol(s: &str) -> Option<Element> {
        (1..=HEAVIEST).map(Element).find(|el| el.symbol() == s)
    }

    /// Which of the engine's eight lumped buckets this element counts toward.
    ///
    /// The bridge between the two accounts. Everything the coarse tiers cannot
    /// name individually lands in `Other`, which is what that bucket is for.
    pub fn species(self) -> crate::units::Species {
        use crate::units::Species;
        match self.0 {
            1 => Species::Hydrogen,
            2 => Species::Helium,
            6 => Species::Carbon,
            7 => Species::Nitrogen,
            8 => Species::Oxygen,
            14 => Species::Silicon,
            26 => Species::Iron,
            _ => Species::Other,
        }
    }
}

impl std::fmt::Display for Element {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.symbol())
    }
}
