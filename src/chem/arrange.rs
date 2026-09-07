//! What a substance is: which atoms, joined how.
//!
//! # Arrangement, not name
//!
//! Nothing here says "water". An [`Arrangement`] is two hydrogens and an oxygen
//! with two bonds between them, and every property the engine uses is derived
//! from that by `analyse` — mass exactly, bond energies and lengths from the
//! elements' own measured constants, polarity from the geometry. The name is a
//! label a person may attach afterwards and the physics never reads.
//!
//! That is the difference between a table of substances, which can only
//! describe what somebody wrote down, and this, which describes anything that
//! can be built. A player who assembles something nobody anticipated gets its
//! real properties, because the properties were never stored against a name in
//! the first place.
//!
//! # Two shapes of matter
//!
//! A [`Molecule`] is a finite graph: this many atoms, these bonds, and then it
//! ends. Water, caffeine, a protein.
//!
//! A [`Lattice`] repeats: a unit cell and the vectors that tile it. Salt, ice,
//! quartz, steel. The distinction is not cosmetic — a molecular substance melts
//! when its intermolecular forces give way and a lattice melts when its bonds
//! do, which is why sugar melts at 460 K and salt at 1074 K.
//!
//! # Canonical form
//!
//! Two arrangements that describe the same substance must be recognised as the
//! same however they were written down: the same atoms in a different order,
//! the same bonds listed backwards. [`Arrangement::canonical`] rewrites one
//! into a form that depends only on the graph, and [`Arrangement::fingerprint`]
//! hashes it, which is what lets the registry deduplicate.
//!
//! The canonicalisation is a Morgan-style refinement: each atom is labelled by
//! its element and degree, then repeatedly relabelled by the multiset of its
//! neighbours' labels until the labelling stops changing. This distinguishes
//! everything the engine can currently build. It is not a complete graph
//! canonicalisation — two different molecules can in principle refine to the
//! same labelling, and for such a pair the registry would treat them as one
//! substance. `analyse` therefore also compares formulas and bond multisets
//! before declaring a match, which closes every case short of a genuine
//! regular-graph collision.

use std::collections::BTreeMap;

use super::Element;

/// How strongly two atoms are joined.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Order {
    Single,
    Double,
    Triple,
    /// Charge transfer rather than sharing: the two atoms are ions held by
    /// electrostatics. Derived from electronegativity difference rather than
    /// declared — see `analyse::ionicity` — but recorded here when a caller
    /// builds a lattice it already knows to be ionic.
    Ionic,
    /// A hydrogen bond: too weak to be a bond and too strong to ignore, and the
    /// entire reason water behaves as it does.
    Hydrogen,
}

impl Order {
    /// Multiplier on the single-bond energy.
    pub fn strength(self) -> f64 {
        match self {
            Order::Single => 1.0,
            Order::Double => 1.9,
            Order::Triple => 2.7,
            Order::Ionic => 1.0,
            Order::Hydrogen => 0.05,
        }
    }

    /// How many of an atom's valence slots it occupies.
    pub fn slots(self) -> usize {
        match self {
            Order::Single | Order::Ionic => 1,
            Order::Double => 2,
            Order::Triple => 3,
            // A hydrogen bond is an association, not a shared pair, so it does
            // not compete for valence. Counting it would make ordinary water
            // illegal.
            Order::Hydrogen => 0,
        }
    }
}

/// One bond, by atom index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Bond {
    pub a: u16,
    pub b: u16,
    pub order: Order,
}

impl Bond {
    pub fn new(a: usize, b: usize, order: Order) -> Bond {
        Bond { a: a.min(b) as u16, b: a.max(b) as u16, order }
    }
}

/// How a unit cell repeats.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Lattice {
    /// It does not: a finite molecule.
    Molecular,
    /// Cubic, with the given cell edge in metres. Rock salt, diamond, iron.
    Cubic { a: f64 },
    /// Hexagonal close packed or similar, cell edges in metres.
    Hexagonal { a: f64, c: f64 },
}

/// A substance's structure.
#[derive(Debug, Clone, PartialEq)]
pub struct Arrangement {
    /// Atoms, by element. Index into this is what `Bond` refers to.
    pub atoms: Vec<Element>,
    pub bonds: Vec<Bond>,
    pub lattice: Lattice,
    /// Net charge in elementary units. Non-zero makes this an ion, which is
    /// what a dissolved salt actually consists of.
    pub charge: i8,
}

/// Elemental counts. What is conserved when a substance reacts.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Formula(pub BTreeMap<Element, u32>);

impl Formula {
    /// Hill notation: carbon first, then hydrogen, then everything else
    /// alphabetically. `H2O`, `C8H10N4O2`.
    pub fn hill(&self) -> String {
        let mut out = String::new();
        let mut rest: Vec<(&Element, &u32)> = self.0.iter().collect();
        rest.sort_by_key(|(e, _)| e.symbol());
        let push = |out: &mut String, e: Element, n: u32| {
            out.push_str(e.symbol());
            if n > 1 {
                out.push_str(&n.to_string());
            }
        };
        if let Some(n) = self.0.get(&Element::CARBON) {
            push(&mut out, Element::CARBON, *n);
            if let Some(h) = self.0.get(&Element::HYDROGEN) {
                push(&mut out, Element::HYDROGEN, *h);
            }
        }
        for (e, n) in rest {
            let carbon_block = self.0.contains_key(&Element::CARBON)
                && (*e == Element::CARBON || *e == Element::HYDROGEN);
            if carbon_block {
                continue;
            }
            push(&mut out, *e, *n);
        }
        out
    }

    /// Total mass of one formula unit, kg. `None` if any element is off the
    /// table, because a molar mass that silently omitted an atom would be
    /// worse than no answer.
    pub fn mass(&self) -> Option<f64> {
        let mut total = 0.0;
        for (e, n) in &self.0 {
            total += e.mass_kg()? * *n as f64;
        }
        Some(total)
    }

    pub fn atoms(&self) -> u32 {
        self.0.values().sum()
    }

    /// This formula's mass split across the engine's eight lumped buckets, as
    /// fractions summing to one.
    ///
    /// The bridge back to `Matter::composition`. A substance's contribution
    /// to a node's bulk elemental account is exactly this, scaled by how much
    /// of the node is that substance — which is the invariant the registry
    /// tests.
    pub fn as_composition(&self) -> Option<crate::state::Composition> {
        let mut c = [0.0f64; crate::units::COARSE_ELEMENTS];
        let mut total = 0.0;
        for (e, n) in &self.0 {
            let m = e.mass_kg()? * *n as f64;
            c[e.coarse_element() as usize] += m;
            total += m;
        }
        if total > 0.0 {
            for slot in c.iter_mut() {
                *slot /= total;
            }
        }
        Some(crate::state::Composition(c))
    }
}

impl Arrangement {
    /// A finite molecule.
    pub fn molecule(atoms: Vec<Element>, bonds: Vec<Bond>) -> Arrangement {
        Arrangement { atoms, bonds, lattice: Lattice::Molecular, charge: 0 }
    }

    /// A single atom or a monatomic ion.
    pub fn atom(e: Element, charge: i8) -> Arrangement {
        Arrangement { atoms: vec![e], bonds: Vec::new(), lattice: Lattice::Molecular, charge }
    }

    /// A repeating solid.
    pub fn crystal(atoms: Vec<Element>, bonds: Vec<Bond>, lattice: Lattice) -> Arrangement {
        Arrangement { atoms, bonds, lattice, charge: 0 }
    }

    pub fn formula(&self) -> Formula {
        let mut f = BTreeMap::new();
        for a in &self.atoms {
            *f.entry(*a).or_insert(0) += 1;
        }
        Formula(f)
    }

    /// Bonds at each atom, as (neighbour, order).
    pub fn neighbours(&self) -> Vec<Vec<(usize, Order)>> {
        let mut out = vec![Vec::new(); self.atoms.len()];
        for b in &self.bonds {
            let (a, c) = (b.a as usize, b.b as usize);
            if a < out.len() && c < out.len() && a != c {
                out[a].push((c, b.order));
                out[c].push((a, b.order));
            }
        }
        out
    }

    /// Valence slots used at each atom.
    pub fn used_slots(&self) -> Vec<usize> {
        let mut out = vec![0usize; self.atoms.len()];
        for b in &self.bonds {
            for i in [b.a as usize, b.b as usize] {
                if i < out.len() {
                    out[i] += b.order.slots();
                }
            }
        }
        out
    }

    /// A relabelling that depends only on the graph, not on the order the
    /// atoms happen to have been listed in.
    ///
    /// Morgan-style refinement: start from (element, degree), then repeatedly
    /// replace each atom's label with a hash of its own label and the sorted
    /// multiset of its neighbours' labels. Two atoms that are genuinely
    /// interchangeable keep the same label; everything else separates.
    pub fn invariants(&self) -> Vec<u64> {
        let n = self.atoms.len();
        let adj = self.neighbours();
        let mut label: Vec<u64> = (0..n)
            .map(|i| {
                let mut h = fnv(self.atoms[i].z() as u64);
                h = fnv(h ^ adj[i].len() as u64);
                h
            })
            .collect();
        // n rounds is enough for the label to have travelled the diameter of
        // any connected graph this size.
        for _ in 0..n.min(16) {
            let mut next = Vec::with_capacity(n);
            for i in 0..n {
                let mut around: Vec<u64> =
                    adj[i].iter().map(|(j, o)| fnv(label[*j] ^ (*o as u64 + 1))).collect();
                around.sort_unstable();
                let mut h = label[i];
                for a in around {
                    h = fnv(h ^ a);
                }
                next.push(h);
            }
            if next == label {
                break;
            }
            label = next;
        }
        label
    }

    /// The arrangement rewritten with its atoms in canonical order.
    ///
    /// Ties in the refined labelling are broken by the original index, so this
    /// is deterministic for every input and idempotent: canonicalising twice
    /// gives the same thing.
    pub fn canonical(&self) -> Arrangement {
        let inv = self.invariants();
        let mut order: Vec<usize> = (0..self.atoms.len()).collect();
        order.sort_by_key(|&i| (inv[i], self.atoms[i].z(), i));
        let mut position = vec![0usize; order.len()];
        for (new, &old) in order.iter().enumerate() {
            position[old] = new;
        }
        let atoms = order.iter().map(|&i| self.atoms[i]).collect();
        let mut bonds: Vec<Bond> = self
            .bonds
            .iter()
            .filter(|b| (b.a as usize) < position.len() && (b.b as usize) < position.len())
            .map(|b| Bond::new(position[b.a as usize], position[b.b as usize], b.order))
            .collect();
        bonds.sort();
        bonds.dedup();
        Arrangement { atoms, bonds, lattice: self.lattice, charge: self.charge }
    }

    /// A hash of the canonical form. Equal fingerprints mean the registry will
    /// treat two arrangements as the same substance.
    pub fn fingerprint(&self) -> u64 {
        let c = self.canonical();
        let mut h = fnv(c.atoms.len() as u64);
        for a in &c.atoms {
            h = fnv(h ^ a.z() as u64);
        }
        for b in &c.bonds {
            h = fnv(h ^ ((b.a as u64) << 24 | (b.b as u64) << 8 | b.order as u64));
        }
        h = fnv(h ^ (c.charge as i64 as u64));
        match c.lattice {
            Lattice::Molecular => fnv(h ^ 1),
            Lattice::Cubic { a } => fnv(h ^ 2 ^ a.to_bits()),
            Lattice::Hexagonal { a, c } => fnv(fnv(h ^ 3 ^ a.to_bits()) ^ c.to_bits()),
        }
    }

    /// Whether two arrangements describe the same substance.
    ///
    /// Compares the canonical forms outright rather than trusting the hash,
    /// and the formula and bond multiset alongside, so the one theoretical
    /// weakness of Morgan refinement — two distinct regular graphs refining to
    /// the same labelling — cannot merge two real substances.
    pub fn same_as(&self, other: &Arrangement) -> bool {
        if self.formula() != other.formula() || self.charge != other.charge {
            return false;
        }
        let (a, b) = (self.canonical(), other.canonical());
        a.atoms == b.atoms && a.bonds == b.bonds && a.lattice == b.lattice
    }
}

#[inline]
fn fnv(v: u64) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for byte in v.to_le_bytes() {
        h ^= byte as u64;
        h = h.wrapping_mul(0x1000_0000_01b3);
    }
    h
}
