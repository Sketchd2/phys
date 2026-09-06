//! The catalogue of substances, and what a node is made of.
//!
//! # Analysed once, then compiled
//!
//! `analyse` derives a substance's properties from its arrangement every time
//! it is called, and for a molecule of any size that is not free. The registry
//! is where the answer is kept: hand it an arrangement, get back a
//! [`SubstanceId`], and every later question about that substance is a lookup.
//!
//! Which is exactly the rule this project works to — derived and stored, never
//! pre-tabulated. Nothing in this file knows what water is. It knows how to
//! analyse an arrangement it has never seen, and it remembers what it found.
//!
//! # Deduplication, with fuzziness
//!
//! Two players who independently build the same molecule must end up with the
//! same substance, or the world accumulates a thousand indistinguishable kinds
//! of water. [`Arrangement::fingerprint`] handles the identical case.
//!
//! The fuzzy case is mixtures rather than molecules: a coffee at 7.1 grams per
//! litre and one at 7.2 are the same drink, and storing both is waste. So a
//! [`Mixture`] rounds its fractions to [`Mixture::TOLERANCE`] before it is
//! compared or hashed, and two recipes within that of each other are one
//! recipe. A caller who means the difference asks for it with
//! [`Registry::intern_exact`], which is the escape hatch for the case where
//! the ratio really is the point.
//!
//! # Two accounts, reconciled
//!
//! A node carries both `Aggregate::composition` — eight lumped species, which
//! is what nuclear burning and mass conservation run on — and a `Mixture`,
//! which says how much of that mass is in *which substances*. The second is a
//! speciation of the first and may be partial: at ten million kelvin there are
//! no molecules, so the mixture is empty and the composition is the whole
//! story.
//!
//! [`Mixture::composition`] projects a mixture back down to the eight buckets,
//! and it must agree with the aggregate's own composition to within round-off.
//! That is the invariant tying the two together, and `tests/chem.rs` asserts it
//! rather than trusting it.

use std::collections::HashMap;

use super::analyse::{analyse, Illegal, Properties};
use super::arrange::{Arrangement, Formula};
use crate::state::Composition;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SubstanceId(pub u32);

impl SubstanceId {
    /// Matter that has not been speciated: plasma, a monatomic gas, a stellar
    /// interior. Not "unknown" — there is genuinely no molecule there.
    pub const UNSPECIATED: SubstanceId = SubstanceId(u32::MAX);
}

/// Where a substance's properties came from.
///
/// Kept per substance and for good, so that a later build which can derive
/// something this one had to be told can find every substance whose value came
/// from a table and recompute it. That is the upgrade path for the biological
/// properties this engine cannot yet solve for.
#[derive(Debug, Clone, PartialEq)]
pub enum Provenance {
    /// Computed from the arrangement by `analyse`.
    Derived,
    /// Supplied from outside, overriding what was derived. Carries what it
    /// replaced, so the disagreement is visible rather than lost.
    Measured { replaced: Box<Properties> },
}

/// One substance, as the world knows it.
#[derive(Debug, Clone)]
pub struct Substance {
    pub id: SubstanceId,
    pub arrangement: Arrangement,
    pub formula: Formula,
    pub props: Properties,
    pub provenance: Provenance,
    /// A name somebody attached. The physics never reads this: two substances
    /// with the same arrangement are the same substance whatever they are
    /// called, and one with no name behaves exactly like one with a name.
    pub label: Option<String>,
}

/// Everything the world has ever analysed.
#[derive(Debug, Default)]
pub struct Registry {
    substances: Vec<Substance>,
    by_fingerprint: HashMap<u64, Vec<u32>>,
    /// Arrangements actually put through `analyse`.
    pub analyses: u64,
    /// Times an arrangement was recognised instead.
    pub hits: u64,
}

impl Registry {
    pub fn new() -> Registry {
        Registry::default()
    }

    pub fn len(&self) -> usize {
        self.substances.len()
    }

    pub fn is_empty(&self) -> bool {
        self.substances.is_empty()
    }

    pub fn get(&self, id: SubstanceId) -> Option<&Substance> {
        self.substances.get(id.0 as usize)
    }

    pub fn all(&self) -> &[Substance] {
        &self.substances
    }

    /// Analyse an arrangement, or recognise one already known.
    ///
    /// The compile step. An arrangement nobody has built before is analysed and
    /// kept; one that matches something already here returns that instead, so
    /// two players who independently assemble the same molecule get the same
    /// substance.
    pub fn intern(&mut self, arrangement: Arrangement) -> Result<SubstanceId, Illegal> {
        let print = arrangement.fingerprint();
        if let Some(candidates) = self.by_fingerprint.get(&print) {
            for i in candidates {
                if self.substances[*i as usize].arrangement.same_as(&arrangement) {
                    self.hits += 1;
                    return Ok(self.substances[*i as usize].id);
                }
            }
        }
        self.compile(arrangement, print)
    }

    /// Analyse an arrangement as a new substance even if an identical one is
    /// already known.
    ///
    /// For the case where the caller means the difference — a recipe kept apart
    /// on purpose, a sample being tracked separately from its own kind.
    pub fn intern_exact(&mut self, arrangement: Arrangement) -> Result<SubstanceId, Illegal> {
        let print = arrangement.fingerprint();
        self.compile(arrangement, print)
    }

    fn compile(&mut self, arrangement: Arrangement, print: u64) -> Result<SubstanceId, Illegal> {
        let props = analyse(&arrangement)?;
        let id = SubstanceId(self.substances.len() as u32);
        self.substances.push(Substance {
            id,
            formula: arrangement.formula(),
            arrangement,
            props,
            provenance: Provenance::Derived,
            label: None,
        });
        self.by_fingerprint.entry(print).or_default().push(id.0);
        self.analyses += 1;
        Ok(id)
    }

    /// Attach a name. Cosmetic by construction — nothing in the physics reads
    /// it, and two substances with the same arrangement stay the same substance
    /// whatever they are called.
    pub fn name(&mut self, id: SubstanceId, label: &str) {
        if let Some(s) = self.substances.get_mut(id.0 as usize) {
            s.label = Some(label.to_string());
        }
    }

    /// Override derived properties with measured ones.
    ///
    /// The honest route for anything this engine cannot yet solve for. What was
    /// derived is kept alongside, so the disagreement between the model and the
    /// measurement stays visible instead of being quietly replaced — and a
    /// later build that improves the model can find every substance that needed
    /// this and check whether it still does.
    pub fn measured(&mut self, id: SubstanceId, props: Properties) -> bool {
        match self.substances.get_mut(id.0 as usize) {
            Some(s) => {
                let replaced = Box::new(s.props);
                s.props = props;
                s.provenance = Provenance::Measured { replaced };
                true
            }
            None => false,
        }
    }

    /// Find a substance by the name somebody gave it. A convenience for
    /// authoring and for tests; the engine never looks anything up this way.
    pub fn by_name(&self, label: &str) -> Option<SubstanceId> {
        self.substances
            .iter()
            .find(|s| s.label.as_deref() == Some(label))
            .map(|s| s.id)
    }
}

/// How many substances one node can be made of before the smallest are lumped.
///
/// Eight, for the same reason the elemental account has eight buckets: a node
/// is a *bulk* description, and a bulk description that tracked forty trace
/// species would cost more than the detail it stands in for. What falls off
/// the end is not lost — its mass stays in the aggregate's composition, it
/// simply stops being attributed to a named substance.
pub const MIXTURE_SLOTS: usize = 8;

/// What a node is made of, by substance.
///
/// Fixed size and `Copy`, so it can sit in an `Aggregate` without an
/// allocation. Fractions are of the node's total mass and need not sum to one:
/// what is left over is matter with no molecular identity, which at ten million
/// kelvin is all of it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mixture {
    slots: [(SubstanceId, f64); MIXTURE_SLOTS],
    used: u8,
}

impl Default for Mixture {
    fn default() -> Mixture {
        Mixture { slots: [(SubstanceId::UNSPECIATED, 0.0); MIXTURE_SLOTS], used: 0 }
    }
}

impl Mixture {
    /// How close two fractions have to be to count as the same recipe.
    ///
    /// A coffee at 7.1 grams per litre and one at 7.2 are the same drink. One
    /// part in a thousand is fine enough that no two substances a person would
    /// call different collide, and coarse enough that a recipe re-entered by
    /// hand matches the one it was copied from.
    pub const TOLERANCE: f64 = 1e-3;

    pub fn new() -> Mixture {
        Mixture::default()
    }

    pub fn len(&self) -> usize {
        self.used as usize
    }

    pub fn is_empty(&self) -> bool {
        self.used == 0
    }

    pub fn entries(&self) -> &[(SubstanceId, f64)] {
        &self.slots[..self.used as usize]
    }

    /// Mass fraction of one substance.
    pub fn fraction_of(&self, id: SubstanceId) -> f64 {
        self.entries().iter().find(|(s, _)| *s == id).map(|(_, f)| *f).unwrap_or(0.0)
    }

    /// Fraction of the node's mass that has a molecular identity at all.
    pub fn speciated(&self) -> f64 {
        self.entries().iter().map(|(_, f)| *f).sum()
    }

    /// Add mass of a substance, merging with what is already there.
    ///
    /// Returns false if the mixture was full and this was smaller than
    /// everything in it, which is the case where the caller is told its trace
    /// species did not make the cut rather than silently losing it.
    pub fn add(&mut self, id: SubstanceId, fraction: f64) -> bool {
        if !fraction.is_finite() || fraction <= 0.0 || id == SubstanceId::UNSPECIATED {
            return false;
        }
        for slot in self.slots[..self.used as usize].iter_mut() {
            if slot.0 == id {
                slot.1 += fraction;
                return true;
            }
        }
        if (self.used as usize) < MIXTURE_SLOTS {
            self.slots[self.used as usize] = (id, fraction);
            self.used += 1;
            return true;
        }
        // Full. Displace the smallest, if this beats it.
        let (i, smallest) = self.slots[..MIXTURE_SLOTS]
            .iter()
            .enumerate()
            .map(|(i, (_, f))| (i, *f))
            .fold((0, f64::INFINITY), |a, b| if b.1 < a.1 { b } else { a });
        if fraction > smallest {
            self.slots[i] = (id, fraction);
            true
        } else {
            false
        }
    }

    /// Scale every fraction so they sum to `total`.
    pub fn normalise_to(&mut self, total: f64) {
        let sum = self.speciated();
        if sum > 0.0 && total >= 0.0 {
            let k = total / sum;
            for slot in self.slots[..self.used as usize].iter_mut() {
                slot.1 *= k;
            }
        }
    }

    /// Sorted, with fractions rounded to [`Mixture::TOLERANCE`].
    ///
    /// The form two mixtures are compared in, so that "the same recipe written
    /// twice" and "the same recipe measured twice" both come out equal.
    pub fn canonical(&self) -> Mixture {
        let mut out = *self;
        for slot in out.slots[..out.used as usize].iter_mut() {
            slot.1 = (slot.1 / Mixture::TOLERANCE).round() * Mixture::TOLERANCE;
        }
        out.slots[..out.used as usize].sort_by(|a, b| a.0.cmp(&b.0));
        out
    }

    /// Whether two mixtures are the same recipe, within tolerance.
    pub fn same_as(&self, other: &Mixture) -> bool {
        let (a, b) = (self.canonical(), other.canonical());
        a.used == b.used && a.entries() == b.entries()
    }

    /// A hash of the canonical form, for deduplicating recipes.
    pub fn fingerprint(&self) -> u64 {
        let c = self.canonical();
        let mut h = 0xcbf2_9ce4_8422_2325u64;
        for (id, f) in c.entries() {
            for byte in id.0.to_le_bytes() {
                h = (h ^ byte as u64).wrapping_mul(0x1000_0000_01b3);
            }
            for byte in ((f / Mixture::TOLERANCE).round() as i64).to_le_bytes() {
                h = (h ^ byte as u64).wrapping_mul(0x1000_0000_01b3);
            }
        }
        h
    }

    /// This mixture projected onto the engine's eight lumped species.
    ///
    /// The reconciliation between the two accounts. Returns the composition of
    /// the speciated part only, together with what fraction of the mass that
    /// was — a caller comparing against an aggregate has to know how much of it
    /// this claims to explain.
    pub fn composition(&self, reg: &Registry) -> (Composition, f64) {
        let mut acc = [0.0f64; crate::units::NSPECIES];
        let mut explained = 0.0;
        for (id, frac) in self.entries() {
            let Some(s) = reg.get(*id) else { continue };
            let Some(c) = s.formula.as_composition() else { continue };
            for (slot, part) in acc.iter_mut().zip(c.0.iter()) {
                *slot += part * frac;
            }
            explained += frac;
        }
        if explained > 0.0 {
            for slot in acc.iter_mut() {
                *slot /= explained;
            }
        }
        (Composition(acc), explained)
    }

    /// Mean molar mass of the speciated part, kg/mol.
    ///
    /// What the lumped species account cannot give: a kilogram of salt counted
    /// through `Species::Other` contains 2.8 times too few particles, and this
    /// is the number that fixes it.
    pub fn molar_mass(&self, reg: &Registry) -> Option<f64> {
        let mut moles = 0.0;
        let mut mass = 0.0;
        for (id, frac) in self.entries() {
            let s = reg.get(*id)?;
            if s.props.molar_mass > 0.0 {
                moles += frac / s.props.molar_mass;
                mass += frac;
            }
        }
        if moles > 0.0 {
            Some(mass / moles)
        } else {
            None
        }
    }
}
