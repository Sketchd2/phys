//! Saving and loading a world.
//!
//! # What is durable, and what is not
//!
//! The engine's whole argument is that most of the world is *derivable* — a
//! galaxy's stars are a maximum-entropy sample of its aggregate, and coarsening
//! them returns the aggregate exactly. So a save file does not contain the
//! world; it contains the part of the world that could not be worked out again:
//!
//! | Durable | Why |
//! |---|---|
//! | Node aggregates, frames, tiers, addresses | This *is* the world |
//! | Pinned detail (`Tree::persisted`) | Somebody touched it, so it is no longer a sample of anything |
//! | Morphology and topology | A tree that lost a branch has lost it |
//! | Ledger facts | A measurement made a value a fact; re-sampling it would be a lie |
//! | The audit log | The record of where conservation was deliberately broken |
//! | World clock, pace, environments | The scenario's own state |
//!
//! | Transient | Why |
//! |---|---|
//! | Materialised bodies of unpinned nodes | Regenerated from address and epoch, bit-for-bit |
//! | `last_report` | Diagnostic output of the last materialisation, reproduced by the next one |
//! | Histories and clocks | Rebuilt as the simulation runs |
//! | In-flight solver state (shaking, falling) | A save between frames has none |
//! | Observers, frame budget, causal gate | Session and machine configuration, not world state |
//!
//! Leaving the regenerable parts out is not an optimisation bolted on the side.
//! It is the same claim `tests/consistency.rs` already makes, stored: if a
//! reloaded world regenerates detail that differs from what it had, then either
//! prolongation is not deterministic or the epoch bookkeeping is wrong, and
//! `tests/persistence.rs` fails rather than the difference going unnoticed.
//!
//! # Determinism of the file itself
//!
//! Two saves of the same world produce identical bytes. Hash maps are written
//! in sorted key order for exactly this reason: a file that differed run to run
//! would make "did anything actually change" unanswerable, and content
//! addressing impossible later.

use crate::causal::{Influence, InfluenceKind, Mailbox};
use crate::coords::Motion;
use crate::ids::{NodeIdx, PathKey};
use crate::morph::{Environment, Event, EventKind, Morphology, Program};
use crate::observe::{AuthorEvent, Fact, Ledger, Property, Quantity};
use crate::prolong::{MassSpectrum, Profile, ProlongReport, ProlongSpec};
use crate::state::{Aggregate, Body, BodyKind, Composition};
use crate::topology::{Joint, Material, Tie, Topology};
use crate::tree::{Node, Residency, Tree, TreeStats};
use crate::units::{CoarseElement, Tier, COARSE_ELEMENTS};
use crate::wire::{Reader, Result, WireError, Writer};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// enums
// ---------------------------------------------------------------------------

fn put_tier(w: &mut Writer, t: Tier) {
    w.u8(t.index() as u8);
}
fn get_tier(r: &mut Reader) -> Result<Tier> {
    let t = r.tag("tier", Tier::ALL.len() as u8)?;
    Ok(Tier::ALL[t as usize])
}

const BODY_KINDS: [BodyKind; 12] = [
    BodyKind::Super,
    BodyKind::Star,
    BodyKind::CompactObject,
    BodyKind::Planet,
    BodyKind::GasParcel,
    BodyKind::Grain,
    BodyKind::Molecule,
    BodyKind::Atom,
    BodyKind::Nucleus,
    BodyKind::Nucleon,
    BodyKind::Electron,
    BodyKind::Photon,
];

fn put_body_kind(w: &mut Writer, k: BodyKind) {
    w.u8(k as u8);
}
fn get_body_kind(r: &mut Reader) -> Result<BodyKind> {
    let t = r.tag("body kind", BODY_KINDS.len() as u8)?;
    Ok(BODY_KINDS[t as usize])
}

const RESIDENCIES: [Residency; 4] = [
    Residency::Speculative,
    Residency::Observed,
    Residency::Causal,
    Residency::Pinned,
];

fn put_residency(w: &mut Writer, v: Residency) {
    let t = RESIDENCIES.iter().position(|&x| x == v).unwrap_or(0);
    w.u8(t as u8);
}
fn get_residency(r: &mut Reader) -> Result<Residency> {
    let t = r.tag("residency", RESIDENCIES.len() as u8)?;
    Ok(RESIDENCIES[t as usize])
}

const PROGRAMS: [Program; 4] = [Program::Tree, Program::Coral, Program::Tower, Program::Wall];

fn put_program(w: &mut Writer, p: Program) {
    let t = PROGRAMS.iter().position(|&x| x == p).unwrap_or(0);
    w.u8(t as u8);
}
fn get_program(r: &mut Reader) -> Result<Program> {
    let t = r.tag("program", PROGRAMS.len() as u8)?;
    Ok(PROGRAMS[t as usize])
}

const EVENT_KINDS: [EventKind; 4] = [
    EventKind::Severed,
    EventKind::Damaged,
    EventKind::Suppressed,
    EventKind::Completed,
];

fn put_event_kind(w: &mut Writer, k: EventKind) {
    let t = EVENT_KINDS.iter().position(|&x| x == k).unwrap_or(0);
    w.u8(t as u8);
}
fn get_event_kind(r: &mut Reader) -> Result<EventKind> {
    let t = r.tag("event kind", EVENT_KINDS.len() as u8)?;
    Ok(EVENT_KINDS[t as usize])
}

const QUANTITIES: [Quantity; 7] = [
    Quantity::DecayTime,
    Quantity::QuantumLevel,
    Quantity::Position,
    Quantity::Momentum,
    Quantity::Spin,
    Quantity::PhotonCount,
    Quantity::Temperature,
];

fn put_quantity(w: &mut Writer, q: Quantity) {
    let t = QUANTITIES.iter().position(|&x| x == q).unwrap_or(0);
    w.u8(t as u8);
}
/// Decode a quantity tag that arrived as a database column rather than from
/// the wire.
#[cfg(feature = "postgres")]
pub(crate) fn quantity_from(tag: u8) -> Result<Quantity> {
    QUANTITIES
        .get(tag as usize)
        .copied()
        .ok_or(crate::wire::WireError::BadTag { what: "quantity", tag: tag as u64 })
}

/// As `quantity_from`, for an authored property.
#[cfg(feature = "postgres")]
pub(crate) fn property_from(tag: u8) -> Result<Property> {
    PROPERTIES
        .get(tag as usize)
        .copied()
        .ok_or(crate::wire::WireError::BadTag { what: "property", tag: tag as u64 })
}

fn get_quantity(r: &mut Reader) -> Result<Quantity> {
    let t = r.tag("quantity", QUANTITIES.len() as u8)?;
    Ok(QUANTITIES[t as usize])
}

const PROPERTIES: [Property; 6] = [
    Property::Mass,
    Property::Temperature,
    Property::Radius,
    Property::Charge,
    Property::Luminosity,
    Property::TimeRate,
];

fn put_property(w: &mut Writer, p: Property) {
    let t = PROPERTIES.iter().position(|&x| x == p).unwrap_or(0);
    w.u8(t as u8);
}
fn get_property(r: &mut Reader) -> Result<Property> {
    let t = r.tag("property", PROPERTIES.len() as u8)?;
    Ok(PROPERTIES[t as usize])
}

// ---------------------------------------------------------------------------
// leaves
// ---------------------------------------------------------------------------

fn put_composition(w: &mut Writer, c: &Composition) {
    for s in CoarseElement::ALL {
        w.f64(c.get(s));
    }
}
fn get_composition(r: &mut Reader) -> Result<Composition> {
    let mut a = [0.0f64; COARSE_ELEMENTS];
    for slot in a.iter_mut() {
        *slot = r.f64()?;
    }
    Ok(Composition(a))
}

fn put_body(w: &mut Writer, b: &Body) {
    w.vec3(b.pos);
    w.vec3(b.vel);
    w.f64(b.mass);
    w.f64(b.radius);
    w.f64(b.charge);
    w.f64(b.internal_energy);
    w.f64(b.temperature);
    put_composition(w, &b.composition);
    w.vec3(b.spin);
    w.u32(b.slot);
    put_body_kind(w, b.kind);
}
fn get_body(r: &mut Reader) -> Result<Body> {
    Ok(Body {
        pos: r.vec3()?,
        vel: r.vec3()?,
        mass: r.f64()?,
        radius: r.f64()?,
        charge: r.f64()?,
        internal_energy: r.f64()?,
        temperature: r.f64()?,
        composition: get_composition(r)?,
        spin: r.vec3()?,
        slot: r.u32()?,
        kind: get_body_kind(r)?,
    })
}

/// Smallest number of bytes one body can occupy on the wire. Used to bound a
/// length prefix before allocating; it must never overstate the true size.
const BODY_MIN_BYTES: usize = 8 * (3 + 3 + 5 + COARSE_ELEMENTS + 3) + 4 + 1;

pub(crate) fn put_bodies_pub(w: &mut Writer, bodies: &[Body]) {
    w.seq(bodies.len());
    for b in bodies {
        put_body(w, b);
    }
}
pub(crate) fn get_bodies_pub(r: &mut Reader) -> Result<Vec<Body>> {
    let n = r.seq("bodies", BODY_MIN_BYTES)?;
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(get_body(r)?);
    }
    Ok(v)
}

pub(crate) fn put_aggregate(w: &mut Writer, a: &Aggregate) {
    w.f64(a.mass);
    w.vec3(a.com);
    w.vec3(a.momentum);
    w.vec3(a.spin);
    w.f64(a.internal_energy);
    w.f64(a.binding_energy);
    w.f64(a.external_potential);
    w.f64(a.radius);
    w.f64(a.temperature);
    put_composition(w, &a.composition);
    w.f64(a.charge);
    w.f64(a.baryon_number);
    w.f64(a.lepton_number);
    w.f64(a.entropy);
    w.f64(a.entropy_exported);
    w.f64(a.chemical_energy);
    w.f64(a.magnetic_energy);
    w.f64(a.luminosity);
}
pub(crate) fn get_aggregate(r: &mut Reader) -> Result<Aggregate> {
    Ok(Aggregate {
        mass: r.f64()?,
        com: r.vec3()?,
        momentum: r.vec3()?,
        spin: r.vec3()?,
        internal_energy: r.f64()?,
        binding_energy: r.f64()?,
        external_potential: r.f64()?,
        radius: r.f64()?,
        temperature: r.f64()?,
        composition: get_composition(r)?,
        charge: r.f64()?,
        baryon_number: r.f64()?,
        lepton_number: r.f64()?,
        entropy: r.f64()?,
        entropy_exported: r.f64()?,
        chemical_energy: r.f64()?,
        magnetic_energy: r.f64()?,
        luminosity: r.f64()?,
    })
}

fn put_motion(w: &mut Writer, f: &Motion) {
    w.vec3(f.offset);
    w.vec3(f.velocity);
    w.quat(f.orientation);
    w.vec3(f.spin_rate);
    w.f64(f.proper_time);
}
fn get_motion(r: &mut Reader) -> Result<Motion> {
    Ok(Motion {
        offset: r.vec3()?,
        velocity: r.vec3()?,
        orientation: r.quat()?,
        spin_rate: r.vec3()?,
        proper_time: r.f64()?,
    })
}

fn put_profile(w: &mut Writer, p: &Profile) {
    match p {
        Profile::Uniform => w.u8(0),
        Profile::Plummer => w.u8(1),
        Profile::Disk { scale_height_ratio } => {
            w.u8(2);
            w.f64(*scale_height_ratio);
        }
        Profile::Shell => w.u8(3),
        Profile::WoodsSaxon => w.u8(4),
        Profile::Lattice => w.u8(5),
    }
}
fn get_profile(r: &mut Reader) -> Result<Profile> {
    Ok(match r.tag("profile", 6)? {
        0 => Profile::Uniform,
        1 => Profile::Plummer,
        2 => Profile::Disk { scale_height_ratio: r.f64()? },
        3 => Profile::Shell,
        4 => Profile::WoodsSaxon,
        _ => Profile::Lattice,
    })
}

fn put_spectrum(w: &mut Writer, s: &MassSpectrum) {
    match s {
        MassSpectrum::Equal => w.u8(0),
        MassSpectrum::Kroupa { min_msun, max_msun } => {
            w.u8(1);
            w.f64(*min_msun);
            w.f64(*max_msun);
        }
        MassSpectrum::PowerLaw { alpha, ratio } => {
            w.u8(2);
            w.f64(*alpha);
            w.f64(*ratio);
        }
        MassSpectrum::CoarseElement => w.u8(3),
    }
}
fn get_spectrum(r: &mut Reader) -> Result<MassSpectrum> {
    Ok(match r.tag("mass spectrum", 4)? {
        0 => MassSpectrum::Equal,
        1 => MassSpectrum::Kroupa { min_msun: r.f64()?, max_msun: r.f64()? },
        2 => MassSpectrum::PowerLaw { alpha: r.f64()?, ratio: r.f64()? },
        _ => MassSpectrum::CoarseElement,
    })
}

pub(crate) fn put_spec(w: &mut Writer, s: &ProlongSpec) {
    w.u32(s.count as u32);
    put_profile(w, &s.profile);
    put_spectrum(w, &s.spectrum);
    put_body_kind(w, s.kind);
    w.f64(s.composition_scatter);
    w.f64(s.turbulent_fraction);
}
pub(crate) fn get_spec(r: &mut Reader) -> Result<ProlongSpec> {
    Ok(ProlongSpec {
        count: r.u32()? as usize,
        profile: get_profile(r)?,
        spectrum: get_spectrum(r)?,
        kind: get_body_kind(r)?,
        composition_scatter: r.f64()?,
        turbulent_fraction: r.f64()?,
    })
}


// ---------------------------------------------------------------------------
// chemistry
// ---------------------------------------------------------------------------
//
// The substance catalogue is world state, not scenery. A node's mixture names
// substances by id, so a world reloaded without its registry would be pointing
// at nothing — and the ids are positions in the catalogue, so the order it is
// written in *is* the identity. It is written as a list rather than a map for
// exactly that reason.

const SUBSTANCE_MIN_BYTES: usize = 4 + 4 + 1 + 8 * 12 + 1 + 1;
const MIXTURE_MIN_BYTES: usize = 16 + 4;

fn put_element(w: &mut Writer, e: crate::chem::Element) {
    w.u8(e.z());
}
fn get_element(r: &mut Reader) -> Result<crate::chem::Element> {
    Ok(crate::chem::Element(r.u8()?))
}

const ORDERS: [crate::chem::Order; 5] = [
    crate::chem::Order::Single,
    crate::chem::Order::Double,
    crate::chem::Order::Triple,
    crate::chem::Order::Ionic,
    crate::chem::Order::Hydrogen,
];

const PHASES: [crate::chem::Phase; 4] = [
    crate::chem::Phase::Solid,
    crate::chem::Phase::Liquid,
    crate::chem::Phase::Gas,
    crate::chem::Phase::Dissolved,
];

fn put_arrangement(w: &mut Writer, a: &crate::chem::Arrangement) {
    w.seq(a.atoms.len());
    for e in &a.atoms {
        put_element(w, *e);
    }
    w.seq(a.bonds.len());
    for b in &a.bonds {
        w.u16(b.a);
        w.u16(b.b);
        w.u8(ORDERS.iter().position(|o| *o == b.order).unwrap_or(0) as u8);
    }
    match a.lattice {
        crate::chem::Lattice::Molecular => w.u8(0),
        crate::chem::Lattice::Cubic { a } => {
            w.u8(1);
            w.f64(a);
        }
        crate::chem::Lattice::Hexagonal { a, c } => {
            w.u8(2);
            w.f64(a);
            w.f64(c);
        }
    }
    w.u8(a.charge as u8);
}

fn get_arrangement(r: &mut Reader) -> Result<crate::chem::Arrangement> {
    let n = r.seq("atoms", 1)?;
    let mut atoms = Vec::with_capacity(n);
    for _ in 0..n {
        atoms.push(get_element(r)?);
    }
    let n = r.seq("bonds", 5)?;
    let mut bonds = Vec::with_capacity(n);
    for _ in 0..n {
        let a = r.u16()?;
        let b = r.u16()?;
        let order = ORDERS[r.tag("bond order", 5)? as usize];
        bonds.push(crate::chem::Bond { a, b, order });
    }
    let lattice = match r.tag("lattice", 3)? {
        0 => crate::chem::Lattice::Molecular,
        1 => crate::chem::Lattice::Cubic { a: r.f64()? },
        _ => crate::chem::Lattice::Hexagonal { a: r.f64()?, c: r.f64()? },
    };
    let charge = r.u8()? as i8;
    Ok(crate::chem::Arrangement { atoms, bonds, lattice, charge })
}

const CONFIDENCES: [crate::chem::Confidence; 4] = [
    crate::chem::Confidence::Exact,
    crate::chem::Confidence::Derived,
    crate::chem::Confidence::Correlated,
    crate::chem::Confidence::Guessed,
];

fn put_properties(w: &mut Writer, p: &crate::chem::Properties) {
    for v in [
        p.unit_mass,
        p.molar_mass,
        p.cohesive_energy,
        p.ionicity,
        p.polarity,
        p.lattice_binding_ev,
        p.dipole,
        p.density,
        p.melting_point,
        p.boiling_point,
        p.water_solubility,
    ] {
        w.f64(v);
    }
    w.u8(p.hydrogen_bonds);
    w.u8(CONFIDENCES.iter().position(|c| *c == p.confidence).unwrap_or(3) as u8);
}

fn get_properties(r: &mut Reader) -> Result<crate::chem::Properties> {
    Ok(crate::chem::Properties {
        unit_mass: r.f64()?,
        molar_mass: r.f64()?,
        cohesive_energy: r.f64()?,
        ionicity: r.f64()?,
        polarity: r.f64()?,
        lattice_binding_ev: r.f64()?,
        dipole: r.f64()?,
        density: r.f64()?,
        melting_point: r.f64()?,
        boiling_point: r.f64()?,
        water_solubility: r.f64()?,
        hydrogen_bonds: r.u8()?,
        confidence: CONFIDENCES[r.tag("confidence", 4)? as usize],
    })
}

pub fn put_registry(w: &mut Writer, reg: &crate::chem::Registry) {
    w.seq(reg.len());
    for s in reg.all() {
        put_arrangement(w, &s.arrangement);
        put_properties(w, &s.props);
        match &s.provenance {
            crate::chem::Provenance::Derived => w.u8(0),
            crate::chem::Provenance::Measured { replaced } => {
                w.u8(1);
                put_properties(w, replaced);
            }
        }
        match &s.label {
            Some(l) => {
                w.bool(true);
                w.str(l);
            }
            None => w.bool(false),
        }
    }
}

pub fn get_registry(r: &mut Reader) -> Result<crate::chem::Registry> {
    let n = r.seq("substances", SUBSTANCE_MIN_BYTES)?;
    let mut reg = crate::chem::Registry::new();
    for i in 0..n {
        let arrangement = get_arrangement(r)?;
        let props = get_properties(r)?;
        let provenance = match r.tag("provenance", 2)? {
            0 => crate::chem::Provenance::Derived,
            _ => crate::chem::Provenance::Measured { replaced: Box::new(get_properties(r)?) },
        };
        let label = if r.bool()? { Some(r.str()?) } else { None };
        // `intern_exact`, not `intern`: the ids are positions in the catalogue
        // and every mixture in the file refers to them, so deduplicating on the
        // way back in would renumber everything. A saved registry is restored
        // exactly as it was saved, including any deliberate duplicates.
        let id = reg
            .intern_exact(arrangement)
            .map_err(|e| WireError::Io(format!("substance {i} does not analyse: {e}")))?;
        reg.restore(id, props, provenance, label);
    }
    Ok(reg)
}

pub fn put_mixture(w: &mut Writer, m: &crate::chem::Mixture) {
    w.seq(m.len());
    for p in m.entries() {
        w.u32(p.substance.0);
        w.u8(PHASES.iter().position(|x| *x == p.phase).unwrap_or(0) as u8);
        w.f64(p.fraction);
    }
}

pub fn get_mixture(r: &mut Reader) -> Result<crate::chem::Mixture> {
    let n = r.seq("pools", 13)?;
    let mut m = crate::chem::Mixture::new();
    for _ in 0..n {
        let id = crate::chem::SubstanceId(r.u32()?);
        let phase = PHASES[r.tag("phase", 4)? as usize];
        m.add(id, phase, r.f64()?);
    }
    Ok(m)
}

// ---------------------------------------------------------------------------
// morphology
// ---------------------------------------------------------------------------

pub(crate) fn put_morphology(w: &mut Writer, m: &Morphology) {
    put_program(w, m.program);
    for g in m.genome {
        w.u32(g.to_bits());
    }
    w.f64(m.age);
    w.f64(m.built);
    w.f64(m.progress);
    w.seq(m.events.len());
    for e in &m.events {
        w.f64(e.at);
        put_event_kind(w, e.kind);
        w.u32(e.site);
        w.f64(e.magnitude);
    }
    w.f64(m.checkpoint_age);
    w.f64(m.design_mass);
}

/// 8 + 1 + 4 + 8 for one event.
const EVENT_MIN_BYTES: usize = 21;

pub(crate) fn get_morphology(r: &mut Reader) -> Result<Morphology> {
    let program = get_program(r)?;
    let mut genome = [0.0f32; 8];
    for g in genome.iter_mut() {
        *g = f32::from_bits(r.u32()?);
    }
    let age = r.f64()?;
    let built = r.f64()?;
    let progress = r.f64()?;
    let n = r.seq("events", EVENT_MIN_BYTES)?;
    let mut events = Vec::with_capacity(n);
    for _ in 0..n {
        events.push(Event {
            at: r.f64()?,
            kind: get_event_kind(r)?,
            site: r.u32()?,
            magnitude: r.f64()?,
        });
    }
    Ok(Morphology {
        program,
        genome,
        age,
        built,
        progress,
        events,
        checkpoint_age: r.f64()?,
        design_mass: r.f64()?,
    })
}

// ---------------------------------------------------------------------------
// topology
// ---------------------------------------------------------------------------

fn put_material(w: &mut Writer, m: &Material) {
    w.str(m.name);
    w.f64(m.density);
    w.f64(m.rupture);
    w.f64(m.tensile_ratio);
    w.f64(m.stiffness);
    w.f64(m.thermal_onset);
    w.f64(m.thermal_gone);
    w.f64(m.destruction_enthalpy);
    w.f64(m.specific_heat);
    w.f64(m.resistivity);
    w.bool(m.combustible);
    w.f64(m.ductility);
}
fn get_material(r: &mut Reader) -> Result<Material> {
    // The name is the one field that cannot come back from bytes: it is a
    // `&'static str`. Recover it from the preset table, or fall back to a
    // static label. Every number below round-trips exactly.
    let name = Material::static_name(&r.str()?);
    Ok(Material {
        name,
        density: r.f64()?,
        rupture: r.f64()?,
        tensile_ratio: r.f64()?,
        stiffness: r.f64()?,
        thermal_onset: r.f64()?,
        thermal_gone: r.f64()?,
        destruction_enthalpy: r.f64()?,
        specific_heat: r.f64()?,
        resistivity: r.f64()?,
        combustible: r.bool()?,
        ductility: r.f64()?,
    })
}

fn put_topology(w: &mut Writer, t: &Topology) {
    w.seq(t.joints.len());
    for b in &t.joints {
        w.u32(b.child);
        w.u32(b.parent);
        w.vec3(b.at);
        w.f64(b.radius);
        w.f64(b.integrity);
    }
    w.seq(t.support.len());
    for &s in &t.support {
        w.u32(s);
    }
    w.seq(t.site.len());
    for &s in &t.site {
        w.u32(s);
    }
    w.seq(t.base.len());
    for &v in &t.base {
        w.vec3(v);
    }
    w.seq(t.tip.len());
    for &v in &t.tip {
        w.vec3(v);
    }
    put_material(w, &t.material);
    w.seq(t.ties.len());
    for tie in &t.ties {
        w.u32(tie.a);
        w.u32(tie.b);
        w.f64(tie.area);
        w.f64(tie.integrity);
    }
}

const BOND_MIN_BYTES: usize = 4 + 4 + 24 + 8 + 8;
const TIE_MIN_BYTES: usize = 4 + 4 + 8 + 8;

fn get_topology(r: &mut Reader) -> Result<Topology> {
    let n = r.seq("bonds", BOND_MIN_BYTES)?;
    let mut joints = Vec::with_capacity(n);
    for _ in 0..n {
        joints.push(Joint {
            child: r.u32()?,
            parent: r.u32()?,
            at: r.vec3()?,
            radius: r.f64()?,
            integrity: r.f64()?,
        });
    }
    let n = r.seq("support", 4)?;
    let mut support = Vec::with_capacity(n);
    for _ in 0..n {
        support.push(r.u32()?);
    }
    let n = r.seq("site", 4)?;
    let mut site = Vec::with_capacity(n);
    for _ in 0..n {
        site.push(r.u32()?);
    }
    let n = r.seq("base", 24)?;
    let mut base = Vec::with_capacity(n);
    for _ in 0..n {
        base.push(r.vec3()?);
    }
    let n = r.seq("tip", 24)?;
    let mut tip = Vec::with_capacity(n);
    for _ in 0..n {
        tip.push(r.vec3()?);
    }
    let material = get_material(r)?;
    let n = r.seq("ties", TIE_MIN_BYTES)?;
    let mut ties = Vec::with_capacity(n);
    for _ in 0..n {
        ties.push(Tie { a: r.u32()?, b: r.u32()?, area: r.f64()?, integrity: r.f64()? });
    }
    Ok(Topology { joints, support, site, base, tip, material, ties })
}

// ---------------------------------------------------------------------------
// nodes and the tree
// ---------------------------------------------------------------------------

fn put_option<T, F: FnOnce(&mut Writer, &T)>(w: &mut Writer, v: &Option<T>, f: F) {
    match v {
        None => w.bool(false),
        Some(x) => {
            w.bool(true);
            f(w, x);
        }
    }
}

pub(crate) fn put_node_payload(w: &mut Writer, n: &Node) {
    w.u128(n.key.0);
    w.u32(n.parent.0);
    w.u32(n.slot);
    w.u32(n.depth);
    put_tier(w, n.tier);
    put_aggregate(w, &n.agg);
    put_motion(w, &n.motion);
    // Only pinned detail is written. Everything else is regenerated from the
    // node's address and epoch, and `tests/persistence.rs` checks that the
    // regenerated bodies match what was discarded.
    if n.pinned {
        put_bodies_pub(w, &n.bodies);
    } else {
        w.seq(0);
    }
    w.f64(n.potential);
    w.seq(n.children.len());
    for c in &n.children {
        w.u32(c.0);
    }
    put_spec(w, &n.spec);
    w.u32(n.epoch);
    w.f64(n.time);
    w.f64(n.last_disturbed);
    w.f64(n.last_solved);
    w.f64(n.last_grown);
    put_residency(w, n.residency);
    w.bool(n.pinned);
    w.f64(n.bubble);
    w.bool(n.alive);
    put_option(w, &n.morphology, put_morphology);
    put_option(w, &n.topology, put_topology);
    w.u64(n.steps_taken);
}

pub(crate) fn get_node_payload(r: &mut Reader) -> Result<Node> {
    let key = PathKey(r.u128()?);
    let parent = NodeIdx(r.u32()?);
    let slot = r.u32()?;
    let depth = r.u32()?;
    let tier = get_tier(r)?;
    let agg = get_aggregate(r)?;
    let motion = get_motion(r)?;
    let bodies = get_bodies_pub(r)?;
    let potential = r.f64()?;
    let n = r.seq("children", 4)?;
    let mut children = Vec::with_capacity(n);
    for _ in 0..n {
        children.push(NodeIdx(r.u32()?));
    }
    let spec = get_spec(r)?;
    Ok(Node {
        key,
        parent,
        slot,
        depth,
        tier,
        agg,
        motion,
        bodies,
        potential,
        children,
        spec,
        epoch: r.u32()?,
        time: r.f64()?,
        last_disturbed: r.f64()?,
        last_solved: r.f64()?,
        last_grown: r.f64()?,
        residency: get_residency(r)?,
        pinned: r.bool()?,
        bubble: r.f64()?,
        alive: r.bool()?,
        morphology: if r.bool()? { Some(get_morphology(r)?) } else { None },
        topology: if r.bool()? { Some(get_topology(r)?) } else { None },
        steps_taken: r.u64()?,
        // Diagnostic output of the last materialisation. Regenerated by the
        // next one; storing it would be storing a derived value.
        last_report: ProlongReport::default(),
    })
}

pub(crate) fn put_tree_stats(w: &mut Writer, s: &TreeStats) {
    w.u64(s.materialisations);
    w.u64(s.coarsenings);
    w.u64(s.idempotent_coarsenings);
    w.u64(s.structures);
    w.u64(s.growth_steps);
    w.u64(s.damage_events);
    w.f64(s.external_energy_absorbed);
    w.u64(s.bodies_created);
    w.u64(s.bodies_discarded);
    w.u64(s.promotions);
    w.u64(s.persisted_bodies);
    w.f64(s.worst_conservation_error);
}
pub(crate) fn get_tree_stats(r: &mut Reader) -> Result<TreeStats> {
    Ok(TreeStats {
        materialisations: r.u64()?,
        coarsenings: r.u64()?,
        idempotent_coarsenings: r.u64()?,
        structures: r.u64()?,
        growth_steps: r.u64()?,
        damage_events: r.u64()?,
        external_energy_absorbed: r.f64()?,
        bodies_created: r.u64()?,
        bodies_discarded: r.u64()?,
        promotions: r.u64()?,
        persisted_bodies: r.u64()?,
        worst_conservation_error: r.f64()?,
    })
}

/// Smallest a node can be: the fixed fields, with every optional part absent.
const NODE_MIN_BYTES: usize = 16 + 4 + 4 + 4 + 1 + 8 * 24 + 8 * 11 + 4 + 4 + 1 + 1 + 1 + 1 + 1 + 8;

fn put_tree(w: &mut Writer, t: &Tree) {
    w.u64(t.world_seed);
    w.u32(t.root.0);
    w.seq(t.nodes.len());
    for n in &t.nodes {
        put_node_payload(w, n);
    }
    // Sorted, so two saves of the same world are byte-identical.
    let mut keys: Vec<&PathKey> = t.persisted.keys().collect();
    keys.sort_by_key(|k| k.0);
    w.seq(keys.len());
    for k in keys {
        w.u128(k.0);
        put_bodies_pub(w, &t.persisted[k]);
    }
    put_tree_stats(w, &t.stats);
}

fn get_tree(r: &mut Reader) -> Result<Tree> {
    let world_seed = r.u64()?;
    let root = NodeIdx(r.u32()?);
    let n = r.seq("nodes", NODE_MIN_BYTES)?;
    let mut nodes = Vec::with_capacity(n);
    for _ in 0..n {
        nodes.push(get_node_payload(r)?);
    }
    let n = r.seq("persisted", 16 + 4)?;
    let mut persisted = HashMap::with_capacity(n);
    for _ in 0..n {
        let k = PathKey(r.u128()?);
        persisted.insert(k, get_bodies_pub(r)?);
    }
    let stats = get_tree_stats(r)?;
    Ok(Tree::restore(nodes, root, world_seed, persisted, stats))
}

// ---------------------------------------------------------------------------
// the world
// ---------------------------------------------------------------------------

/// Everything durable about a world, independent of how it is stored.
///
/// A `World` also carries observers, a frame budget and a causal gate. Those
/// are session and machine configuration rather than facts about the world, so
/// they are not saved: a world reloaded on another machine gets that machine's
/// budget, and whoever opens it supplies their own observers.
pub struct Snapshot {
    /// Reconstructed mailbox, so `view()` has something to borrow.
    pub mailbox_view: Mailbox,
    pub tree: Tree,
    pub ledger: Ledger,
    pub time: f64,
    pub pace: f64,
    pub time_rate: f64,
    pub time_throttle: f64,
    pub paced_to: NodeIdx,
    pub pace_mode: crate::engine::PaceMode,
    pub labour_rate: f64,
    pub rejected_transactions: u64,
    pub environments: HashMap<PathKey, Environment>,
    pub substances: crate::chem::Registry,
    pub mixtures: HashMap<PathKey, crate::chem::Mixture>,
    pub audit: Vec<AuthorEvent>,
    /// Influences posted and not yet arrived. Durable: an impulse in the
    /// light-delay between the act and its landing is an action somebody took,
    /// and a save that dropped it would quietly undo them.
    pub in_flight: Vec<Influence>,
    pub delivered: u64,
    pub in_flight_peak: usize,
}

/// A borrowed view of the same thing, for writing.
///
/// Saving must not clone the tree. A world large enough to be worth saving is
/// large enough that doubling it in memory to serialise it is the wrong shape,
/// and the borrow costs nothing.
pub struct WorldView<'a> {
    pub tree: &'a Tree,
    pub ledger: &'a Ledger,
    pub time: f64,
    pub pace: f64,
    pub time_rate: f64,
    pub time_throttle: f64,
    pub paced_to: NodeIdx,
    pub pace_mode: crate::engine::PaceMode,
    pub labour_rate: f64,
    pub rejected_transactions: u64,
    pub environments: &'a HashMap<PathKey, Environment>,
    pub substances: &'a crate::chem::Registry,
    pub mixtures: &'a HashMap<PathKey, crate::chem::Mixture>,
    pub audit: &'a [AuthorEvent],
    pub mailbox: &'a Mailbox,
}

impl Snapshot {
    pub fn view(&self) -> WorldView<'_> {
        WorldView {
            tree: &self.tree,
            ledger: &self.ledger,
            time: self.time,
            pace: self.pace,
            time_rate: self.time_rate,
            time_throttle: self.time_throttle,
            paced_to: self.paced_to,
            pace_mode: self.pace_mode,
            labour_rate: self.labour_rate,
            rejected_transactions: self.rejected_transactions,
            environments: &self.environments,
            substances: &self.substances,
            mixtures: &self.mixtures,
            audit: &self.audit,
            mailbox: &self.mailbox_view,
        }
    }
}

pub(crate) fn put_environment_pub(w: &mut Writer, e: &Environment) {
    w.f64(e.light_flux);
    w.f64(e.temperature);
    w.f64(e.water);
    w.f64(e.crowding);
    w.f64(e.reservoir_mass);
    w.f64(e.labour);
}
pub(crate) fn get_environment_pub(r: &mut Reader) -> Result<Environment> {
    Ok(Environment {
        light_flux: r.f64()?,
        temperature: r.f64()?,
        water: r.f64()?,
        crowding: r.f64()?,
        reservoir_mass: r.f64()?,
        labour: r.f64()?,
    })
}

/// Serialise a snapshot. Deterministic: the same world always gives the same
/// bytes.
pub fn encode(s: WorldView<'_>) -> Vec<u8> {
    let mut w = Writer::new();
    w.header();
    put_tree(&mut w, s.tree);

    let mut facts: Vec<(PathKey, Quantity, Fact)> = s.ledger.entries().collect();
    facts.sort_by_key(|(k, q, _)| (k.0, *q as u8));
    w.seq(facts.len());
    for (k, q, f) in &facts {
        w.u128(k.0);
        put_quantity(&mut w, *q);
        w.f64(f.value);
        w.f64(f.time);
        w.u64(f.sequence);
    }
    w.u64(s.ledger.sequence);
    w.u64(s.ledger.queries);
    w.u64(s.ledger.commits);

    w.f64(s.time);
    w.f64(s.pace);
    w.f64(s.time_rate);
    w.f64(s.time_throttle);
    w.u32(s.paced_to.0);
    w.bool(s.pace_mode == crate::engine::PaceMode::Fixed);
    w.f64(s.labour_rate);
    w.u64(s.rejected_transactions);

    let mut envs: Vec<(&PathKey, &Environment)> = s.environments.iter().collect();
    envs.sort_by_key(|(k, _)| k.0);
    w.seq(envs.len());
    for (k, e) in envs {
        w.u128(k.0);
        put_environment_pub(&mut w, e);
    }

    put_registry(&mut w, s.substances);
    let mut mixes: Vec<(&PathKey, &crate::chem::Mixture)> = s.mixtures.iter().collect();
    mixes.sort_by_key(|(k, _)| k.0);
    w.seq(mixes.len());
    for (k, m) in mixes {
        w.u128(k.0);
        put_mixture(&mut w, m);
    }

    w.seq(s.audit.len());
    for a in s.audit {
        w.u128(a.key.0);
        put_property(&mut w, a.property);
        w.f64(a.delta_energy);
        w.f64(a.time);
    }

    // Sorted by arrival, so two saves of one world are the same bytes even
    // though a binary heap does not iterate in order.
    let mut flying: Vec<&Influence> = s.mailbox.in_flight().collect();
    flying.sort_by(|a, b| {
        a.arrives
            .partial_cmp(&b.arrives)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.target.0.cmp(&b.target.0))
    });
    w.seq(flying.len());
    for i in flying {
        w.f64(i.arrives);
        w.u32(i.target.0);
        put_influence_kind(&mut w, i.kind);
        w.f64(i.energy);
        w.vec3(i.momentum);
        w.f64(i.source_distance);
    }
    w.u64(s.mailbox.delivered);
    w.u32(s.mailbox.in_flight_peak as u32);
    w.finish()
}

const FACT_MIN_BYTES: usize = 16 + 1 + 8 + 8 + 8;
const ENV_MIN_BYTES: usize = 16 + 8 * 6;
const AUDIT_MIN_BYTES: usize = 16 + 1 + 8 + 8;
const INFLUENCE_MIN_BYTES: usize = 8 + 4 + 1 + 8 + 24 + 8;

const INFLUENCE_KINDS: [InfluenceKind; 5] = [
    InfluenceKind::Radiation,
    InfluenceKind::Blast,
    InfluenceKind::Impact,
    InfluenceKind::Probe,
    InfluenceKind::UserImpulse,
];

pub(crate) fn influence_kind_tag(k: InfluenceKind) -> u8 {
    INFLUENCE_KINDS.iter().position(|&x| x == k).unwrap_or(0) as u8
}

/// Decode an influence-kind tag that arrived as a database column.
#[cfg(feature = "postgres")]
pub(crate) fn influence_kind_from(tag: u8) -> Result<InfluenceKind> {
    INFLUENCE_KINDS
        .get(tag as usize)
        .copied()
        .ok_or(crate::wire::WireError::BadTag { what: "influence kind", tag: tag as u64 })
}

fn put_influence_kind(w: &mut Writer, k: InfluenceKind) {
    w.u8(influence_kind_tag(k));
}
fn get_influence_kind(r: &mut Reader) -> Result<InfluenceKind> {
    let t = r.tag("influence kind", INFLUENCE_KINDS.len() as u8)?;
    Ok(INFLUENCE_KINDS[t as usize])
}

/// Parse a snapshot. Never panics, whatever the bytes contain.
pub fn decode(bytes: &[u8]) -> Result<Snapshot> {
    let mut r = Reader::new(bytes);
    r.header()?;
    let tree = get_tree(&mut r)?;

    let n = r.seq("facts", FACT_MIN_BYTES)?;
    let mut facts = Vec::with_capacity(n);
    for _ in 0..n {
        let key = PathKey(r.u128()?);
        let quantity = get_quantity(&mut r)?;
        facts.push((
            key,
            quantity,
            Fact { value: r.f64()?, time: r.f64()?, quantity, sequence: r.u64()? },
        ));
    }
    let sequence = r.u64()?;
    let queries = r.u64()?;
    let commits = r.u64()?;
    let ledger = Ledger::restore(facts, sequence, queries, commits);

    let time = r.f64()?;
    let pace = r.f64()?;
    let time_rate = r.f64()?;
    let time_throttle = r.f64()?;
    let paced_to = NodeIdx(r.u32()?);
    let pace_mode = if r.bool()? {
        crate::engine::PaceMode::Fixed
    } else {
        crate::engine::PaceMode::Follow
    };
    let labour_rate = r.f64()?;
    let rejected_transactions = r.u64()?;

    let n = r.seq("environments", ENV_MIN_BYTES)?;
    let mut environments = HashMap::with_capacity(n);
    for _ in 0..n {
        let k = PathKey(r.u128()?);
        environments.insert(k, get_environment_pub(&mut r)?);
    }

    let substances = get_registry(&mut r)?;
    let n = r.seq("mixtures", MIXTURE_MIN_BYTES)?;
    let mut mixtures = HashMap::with_capacity(n);
    for _ in 0..n {
        let k = PathKey(r.u128()?);
        mixtures.insert(k, get_mixture(&mut r)?);
    }

    let n = r.seq("audit", AUDIT_MIN_BYTES)?;
    let mut audit = Vec::with_capacity(n);
    for _ in 0..n {
        audit.push(AuthorEvent {
            key: PathKey(r.u128()?),
            property: get_property(&mut r)?,
            delta_energy: r.f64()?,
            time: r.f64()?,
        });
    }

    let n = r.seq("in flight", INFLUENCE_MIN_BYTES)?;
    let mut in_flight = Vec::with_capacity(n);
    for _ in 0..n {
        in_flight.push(Influence {
            arrives: r.f64()?,
            target: NodeIdx(r.u32()?),
            kind: get_influence_kind(&mut r)?,
            energy: r.f64()?,
            momentum: r.vec3()?,
            source_distance: r.f64()?,
        });
    }
    let delivered = r.u64()?;
    let in_flight_peak = r.u32()? as usize;

    r.finish()?;
    Ok(Snapshot {
        mailbox_view: Mailbox::restore(in_flight.clone(), delivered, in_flight_peak),
        in_flight,
        delivered,
        in_flight_peak,
        tree,
        ledger,
        time,
        pace,
        time_rate,
        time_throttle,
        paced_to,
        pace_mode,
        labour_rate,
        rejected_transactions,
        environments,
        substances,
        mixtures,
        audit,
    })
}

// ---------------------------------------------------------------------------
// stores
// ---------------------------------------------------------------------------

/// Where a world lives.
///
/// A trait rather than a function because the point of Phase 1 is that the
/// *storage* decision stays reversible. `FileStore` is what exists today; a
/// SpacetimeDB-backed implementation is what §03 of the review proposes, and if
/// that bet goes badly it costs one implementation rather than the project.
pub trait WorldStore {
    /// Write the whole world.
    fn save(&mut self, view: WorldView<'_>) -> Result<()>;

    /// Write only what has changed since the world instant `since`.
    ///
    /// This is the operation that decides whether a store can carry a large
    /// world, and the reason it can is the scheduler: a node that was coasted
    /// has not been re-derived, so there is nothing new to say about it. Writes
    /// are proportional to *events* — solves and disturbances — rather than to
    /// how much world there is.
    ///
    /// The default implementation writes everything, so a blob store is still a
    /// perfectly valid store; it just does more work than it needs to. That
    /// default is what keeps the trait honest as a swap point rather than a
    /// Postgres-shaped hole.
    fn save_since(&mut self, view: WorldView<'_>, since: f64) -> Result<Flushed> {
        let total = view.tree.nodes.len();
        self.save(view)?;
        let _ = since;
        Ok(Flushed { nodes_written: total, nodes_skipped: 0, incremental: false })
    }

    fn load(&mut self) -> Result<Snapshot>;

    /// Bytes the store is holding, for reporting.
    fn size(&self) -> usize;
}

/// What a flush actually did.
#[derive(Debug, Clone, Copy, Default)]
pub struct Flushed {
    pub nodes_written: usize,
    pub nodes_skipped: usize,
    /// False when the store fell back to writing everything.
    pub incremental: bool,
}

/// Which nodes have something new to say since the instant `since`.
///
/// A node is dirty when its dynamics were re-derived (`last_solved`) or
/// something happened to it (`last_disturbed`). A node that was merely carried
/// forward is not: its position at any instant is a closed-form function of the
/// state already stored, so the reader can reconstruct it by coasting — which is
/// exactly what [`Snapshot::catch_up`] does on load.
pub fn dirty_nodes(tree: &Tree, since: f64) -> Vec<usize> {
    (0..tree.nodes.len())
        .filter(|&i| {
            let n = &tree.nodes[i];
            n.last_solved > since || n.last_disturbed > since
        })
        .collect()
}

impl Snapshot {
    /// Bring every node up to the world instant.
    ///
    /// A store that wrote only the dirty nodes left the rest at whatever instant
    /// they were last written at. Carrying them forward is the same closed-form
    /// step the scheduler performs every frame, and it is exact — which is what
    /// makes skipping the write safe rather than lossy.
    pub fn catch_up(&mut self) {
        let instant = self.time;
        for n in self.tree.nodes.iter_mut() {
            if !n.alive {
                continue;
            }
            let dt = instant - n.time;
            if dt > 0.0 {
                n.motion.advance(dt);
                n.time = instant;
            }
        }
    }
}

/// In memory. For tests, and for a world that has not been given a home yet.
#[derive(Default)]
pub struct MemoryStore {
    bytes: Vec<u8>,
}

impl MemoryStore {
    pub fn new() -> MemoryStore {
        MemoryStore::default()
    }
    pub fn raw(&self) -> &[u8] {
        &self.bytes
    }
}

impl WorldStore for MemoryStore {
    fn save(&mut self, view: WorldView<'_>) -> Result<()> {
        self.bytes = encode(view);
        Ok(())
    }
    fn load(&mut self) -> Result<Snapshot> {
        decode(&self.bytes)
    }
    fn size(&self) -> usize {
        self.bytes.len()
    }
}

/// One file on disk.
///
/// The write goes to a temporary file and is then renamed over the target, so
/// an interrupted save leaves the previous world intact rather than a truncated
/// one. A world file is the only copy of everything a player has done to it.
pub struct FileStore {
    path: std::path::PathBuf,
}

impl FileStore {
    pub fn new<P: Into<std::path::PathBuf>>(path: P) -> FileStore {
        FileStore { path: path.into() }
    }
    pub fn path(&self) -> &std::path::Path {
        &self.path
    }
}

impl WorldStore for FileStore {
    fn save(&mut self, view: WorldView<'_>) -> Result<()> {
        use std::io::Write;
        let bytes = encode(view);
        let tmp = self.path.with_extension("tmp");
        {
            let mut f = std::fs::File::create(&tmp).map_err(|e| WireError::Io(e.to_string()))?;
            f.write_all(&bytes).map_err(|e| WireError::Io(e.to_string()))?;
            f.sync_all().map_err(|e| WireError::Io(e.to_string()))?;
        }
        std::fs::rename(&tmp, &self.path).map_err(|e| WireError::Io(e.to_string()))?;
        Ok(())
    }
    fn load(&mut self) -> Result<Snapshot> {
        let bytes = std::fs::read(&self.path).map_err(|e| WireError::Io(e.to_string()))?;
        decode(&bytes)
    }
    fn size(&self) -> usize {
        std::fs::metadata(&self.path).map(|m| m.len() as usize).unwrap_or(0)
    }
}
