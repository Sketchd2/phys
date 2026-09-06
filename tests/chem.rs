//! Chemistry derived rather than tabulated.
//!
//! The rule this engine works to is that only physical laws are axioms.
//! Caffeine is not an axiom: it is eight carbons, ten hydrogens, four nitrogens
//! and two oxygens joined a particular way, and everything about it should
//! follow from that. So the axioms here are per *element* — atomic weight,
//! electronegativity, two radii, one bond energy — and every substance property
//! is computed.
//!
//! These tests are therefore mostly *measurements*. A test that asserted
//! "water's boiling point is derived correctly" would be worth nothing without
//! saying against what; each one below prints the derived value against the
//! real one and holds a tolerance that says how good the model actually is.
//! Where the model is only good to a factor of two, the tolerance says two.

use phys::chem::analyse::{
    bond_energy, bond_length, ionicity, solubility_in, Illegal, WATER_POLARITY,
};
use phys::chem::arrange::{Bond, Lattice, Order};
use phys::chem::*;

fn el(z: u8) -> Element {
    Element(z)
}

fn water() -> Arrangement {
    Arrangement::molecule(
        vec![el(8), el(1), el(1)],
        vec![Bond::new(0, 1, Order::Single), Bond::new(0, 2, Order::Single)],
    )
}
fn methane() -> Arrangement {
    Arrangement::molecule(
        vec![el(6), el(1), el(1), el(1), el(1)],
        (1..5).map(|i| Bond::new(0, i, Order::Single)).collect(),
    )
}
fn carbon_dioxide() -> Arrangement {
    Arrangement::molecule(
        vec![el(6), el(8), el(8)],
        vec![Bond::new(0, 1, Order::Double), Bond::new(0, 2, Order::Double)],
    )
}
fn rock_salt() -> Arrangement {
    // One formula unit of the rock-salt lattice: the cell edge is 564 pm and
    // holds four of them.
    Arrangement::crystal(
        vec![el(11), el(17)],
        vec![Bond::new(0, 1, Order::Ionic)],
        Lattice::Cubic { a: 5.64e-10 / 4f64.cbrt() },
    )
}
fn ethanol() -> Arrangement {
    Arrangement::molecule(
        vec![el(6), el(6), el(8), el(1), el(1), el(1), el(1), el(1), el(1)],
        vec![
            Bond::new(0, 1, Order::Single),
            Bond::new(1, 2, Order::Single),
            Bond::new(0, 3, Order::Single),
            Bond::new(0, 4, Order::Single),
            Bond::new(0, 5, Order::Single),
            Bond::new(1, 6, Order::Single),
            Bond::new(1, 7, Order::Single),
            Bond::new(2, 8, Order::Single),
        ],
    )
}
/// Glycine, which nothing in this engine has ever heard of.
///
/// `H2N–CH2–COOH`. Built by hand rather than transcribed, so what is being
/// tested is that the analyser handles a molecule nobody anticipated — three
/// elements, a double bond, a hydroxyl — and not that somebody typed in a
/// structure file correctly.
fn glycine() -> Arrangement {
    //          0:N  1:C  2:C  3:O(=)  4:O(H)  5..9:H
    let atoms = vec![el(7), el(6), el(6), el(8), el(8), el(1), el(1), el(1), el(1), el(1)];
    let bonds = vec![
        Bond::new(0, 1, Order::Single),  // N-C
        Bond::new(1, 2, Order::Single),  // C-C
        Bond::new(2, 3, Order::Double),  // C=O
        Bond::new(2, 4, Order::Single),  // C-O
        Bond::new(0, 5, Order::Single),  // N-H
        Bond::new(0, 6, Order::Single),  // N-H
        Bond::new(1, 7, Order::Single),  // C-H
        Bond::new(1, 8, Order::Single),  // C-H
        Bond::new(4, 9, Order::Single),  // O-H
    ];
    Arrangement::molecule(atoms, bonds)
}

// ---------------------------------------------------------------------------
// the exact part
// ---------------------------------------------------------------------------

/// Molar mass is a sum of measured atomic weights and has no model in it, so
/// it should be right to the last digit the table carries.
#[test]
fn molar_mass_is_exact() {
    for (name, arr, real) in [
        ("water", water(), 18.015),
        ("methane", methane(), 16.043),
        ("carbon dioxide", carbon_dioxide(), 44.009),
        ("sodium chloride", rock_salt(), 58.44),
        ("ethanol", ethanol(), 46.069),
    ] {
        let p = analyse(&arr).unwrap_or_else(|e| panic!("{name}: {e}"));
        let got = p.molar_mass * 1000.0;
        println!("  {name:<16} {:<8} {got:>8.3} g/mol against {real:>8.3}", arr.formula().hill());
        assert!(
            (got - real).abs() < 0.01,
            "{name}: {got} g/mol against a real {real}"
        );
    }
}

/// Hill notation, so a formula reads the way a chemist writes it.
#[test]
fn formulas_read_correctly() {
    assert_eq!(water().formula().hill(), "H2O");
    assert_eq!(methane().formula().hill(), "CH4");
    assert_eq!(carbon_dioxide().formula().hill(), "CO2");
    assert_eq!(ethanol().formula().hill(), "C2H6O");
    // Hill puts carbon and hydrogen first *only* when there is carbon; salt is
    // alphabetical, so "ClNa" is right and "NaCl" is the common name.
    assert_eq!(rock_salt().formula().hill(), "ClNa");
    assert_eq!(glycine().formula().hill(), "C2H5NO2");
}

// ---------------------------------------------------------------------------
// the derived part
// ---------------------------------------------------------------------------

/// Pauling's ionic character is the discriminator that matters most at play
/// scale: it is what tells the engine a salt from a molecule with neither being
/// named.
#[test]
fn ionic_character_separates_salts_from_molecules() {
    let cases = [
        ("Na-Cl", el(11), el(17), 0.70),
        ("H-O", el(1), el(8), 0.33),
        ("C-H", el(6), el(1), 0.03),
        ("C-C", el(6), el(6), 0.0),
    ];
    for (name, a, b, expect) in cases {
        let got = ionicity(a, b);
        println!("  {name:<6} ionic character {got:.3} (Pauling {expect:.2})");
        assert!((got - expect).abs() < 0.05, "{name}: {got} against {expect}");
    }
    assert!(
        ionicity(el(11), el(17)) > 0.5 && ionicity(el(6), el(1)) < 0.1,
        "the salt/molecule split has to be unambiguous or nothing downstream works"
    );
}

/// Bond energies for pairs nobody tabulated, from one number per element.
#[test]
fn bond_energies_come_out_of_paulings_relation() {
    let kj = |a, b| bond_energy(a, b, Order::Single).unwrap() * phys::units::N_AVOGADRO / 1e3;
    for (name, a, b, real) in [
        ("H-H", el(1), el(1), 436.0),
        ("C-C", el(6), el(6), 346.0),
        ("O-H", el(8), el(1), 463.0),
        ("C-H", el(6), el(1), 413.0),
        ("H-Cl", el(1), el(17), 431.0),
        ("C-O", el(6), el(8), 358.0),
    ] {
        let got = kj(a, b);
        println!("  {name:<6} {got:>7.1} kJ/mol against a measured {real:>7.1}");
        assert!(
            (got - real).abs() / real < 0.25,
            "{name}: {got:.1} kJ/mol against {real}"
        );
    }
}

/// Schomaker–Stevenson lengths: the sum of two covalent radii, shortened where
/// the pair is unequal.
///
/// Hydrogen is the known exception and is kept in the list rather than quietly
/// dropped. Cordero's radii are fitted to crystal structures, where hydrogen
/// appears in X–H bonds; its effective radius in H2 is about 37 pm rather than
/// the 31 pm that fit everything else, so the additive model gives 62 pm for a
/// bond that is really 74. Every other bond here is within a few per cent, and
/// a model that is good to 5% except for one diatomic is worth having as long
/// as it says which one.
#[test]
fn bond_lengths_come_out_of_the_radii() {
    let cases = [
        ("O-H", el(8), el(1), Order::Single, 96.0),
        ("C-C", el(6), el(6), Order::Single, 154.0),
        ("C=O", el(6), el(8), Order::Double, 123.0),
        ("C-Cl", el(6), el(17), Order::Single, 177.0),
        ("N-H", el(7), el(1), Order::Single, 101.0),
        ("C-N", el(6), el(7), Order::Single, 147.0),
    ];
    let mut worst: f64 = 0.0;
    for (name, a, b, order, real_pm) in cases {
        let got = bond_length(a, b, order).unwrap() * 1e12;
        let err = (got - real_pm).abs() / real_pm;
        worst = worst.max(err);
        println!("  {name:<6} {got:>6.1} pm against a measured {real_pm:>6.1}  ({:+.1}%)",
            (got / real_pm - 1.0) * 100.0);
    }
    assert!(worst < 0.12, "worst bond-length error {:.1}%", worst * 100.0);

    let h2 = bond_length(el(1), el(1), Order::Single).unwrap() * 1e12;
    println!("  H-H    {h2:>6.1} pm against a measured   74.0  — the known exception");
    assert!((h2 - 74.0).abs() / 74.0 < 0.20, "even the exception should be within 20%");
}

/// The shape, from VSEPR and lone-pair counting. This is what separates carbon
/// dioxide from water, and nothing about the *graph* can do it.
#[test]
fn geometry_reproduces_the_angles_that_matter() {
    for (name, arr, real) in [
        ("water", water(), 104.5),
        ("carbon dioxide", carbon_dioxide(), 180.0),
        ("methane", methane(), 109.47),
    ] {
        let pos = geometry::embed(&arr);
        let adj = arr.neighbours();
        let centre = (0..arr.atoms.len()).max_by_key(|i| adj[*i].len()).unwrap();
        let dirs: Vec<_> =
            adj[centre].iter().map(|(j, _)| (pos[*j] - pos[centre]).unit()).collect();
        let mut worst: f64 = 0.0;
        for i in 0..dirs.len() {
            for j in i + 1..dirs.len() {
                let a = dirs[i].dot(dirs[j]).clamp(-1.0, 1.0).acos().to_degrees();
                worst = worst.max((a - real).abs());
            }
        }
        println!("  {name:<16} every bond angle within {worst:.2} deg of {real}");
        assert!(worst < 0.5, "{name}: off by {worst:.2} degrees");
    }
}

/// Bond dipoles add or cancel depending on the shape, and that single fact is
/// most of what "does it dissolve" means.
#[test]
fn shape_decides_polarity() {
    let w = analyse(&water()).unwrap();
    let c = analyse(&carbon_dioxide()).unwrap();
    let m = analyse(&methane()).unwrap();
    let debye = |p: &Properties| p.dipole / 3.33564e-30;
    println!(
        "  water {:.2} D (real 1.85), carbon dioxide {:.2e} D (real 0), methane {:.2e} D (real 0)",
        debye(&w),
        debye(&c),
        debye(&m)
    );
    assert!(
        (debye(&w) - 1.85).abs() < 0.5,
        "water should be polar: {:.2} D",
        debye(&w)
    );
    // Both have strongly polar bonds and neither is polar, because they are
    // symmetric. This is the case an earlier version got wrong by asking
    // whether atoms were interchangeable instead of building the shape.
    assert!(debye(&c) < 1e-6, "carbon dioxide is linear: {:.3e} D", debye(&c));
    assert!(debye(&m) < 1e-6, "methane is tetrahedral: {:.3e} D", debye(&m));
}

// ---------------------------------------------------------------------------
// the correlated part, measured honestly
// ---------------------------------------------------------------------------

/// Boiling points from a vaporisation enthalpy over Trouton's constant.
///
/// The tolerance here is 25 per cent and that is the honest number: this is a
/// correlation across four decades of intermolecular strength, not a
/// derivation. It reproduces the *ordering* and the rough magnitude, which is
/// what the engine needs to know whether something is a gas at room
/// temperature.
#[test]
fn boiling_points_are_right_to_about_a_quarter() {
    let cases = [
        ("methane", methane(), 111.7),
        ("water", water(), 373.1),
        ("ethanol", ethanol(), 351.5),
        ("sodium chloride", rock_salt(), 1738.0),
    ];
    let mut worst: f64 = 0.0;
    for (name, arr, real) in cases {
        let p = analyse(&arr).unwrap();
        let err = (p.boiling_point - real).abs() / real;
        worst = worst.max(err);
        println!(
            "  {name:<16} boils at {:>7.1} K against a real {real:>7.1}  ({:+.0}%)",
            p.boiling_point,
            (p.boiling_point / real - 1.0) * 100.0
        );
    }
    assert!(worst < 0.25, "worst boiling-point error {:.0}%", worst * 100.0);
}

/// Melting points are the weakest number this module produces, and the
/// tolerance says so.
#[test]
fn melting_points_are_only_the_right_order() {
    for (name, arr, real) in [
        ("water", water(), 273.1),
        ("methane", methane(), 90.7),
        ("sodium chloride", rock_salt(), 1074.0),
    ] {
        let p = analyse(&arr).unwrap();
        println!(
            "  {name:<16} melts at {:>7.1} K against a real {real:>7.1}  ({:+.0}%)",
            p.melting_point,
            (p.melting_point / real - 1.0) * 100.0
        );
        assert!(
            p.melting_point > real * 0.4 && p.melting_point < real * 2.2,
            "{name}: {:.1} K against {real}",
            p.melting_point
        );
        assert_eq!(p.confidence, Confidence::Correlated);
    }
}

/// The whole point at play scale: salt dissolves and methane does not, and
/// neither was named anywhere.
#[test]
fn salt_dissolves_and_methane_does_not() {
    let salt = analyse(&rock_salt()).unwrap();
    let gas = analyse(&methane()).unwrap();
    let fizz = analyse(&carbon_dioxide()).unwrap();
    println!(
        "  salt {:.3e} kg/kg (real 3.6e-1), carbon dioxide {:.3e} (real 1.5e-3), methane {:.3e} (real 2.2e-5)",
        salt.water_solubility, fizz.water_solubility, gas.water_solubility
    );
    // Ordering, which must be exactly right.
    assert!(salt.water_solubility > fizz.water_solubility);
    assert!(fizz.water_solubility > gas.water_solubility);
    // Magnitude, to within a decade and a half — this spans five decades and is
    // reported as a power of ten for that reason.
    for (name, got, real) in [
        ("salt", salt.water_solubility, 0.36),
        ("carbon dioxide", fizz.water_solubility, 1.5e-3),
        ("methane", gas.water_solubility, 2.2e-5),
    ] {
        let decades = (got / real).log10().abs();
        assert!(decades < 1.5, "{name} off by {decades:.2} decades");
    }
}

/// Every solubility is measured against water, so water's own position on the
/// polarity axis is a constant rather than a call — analysing water to find out
/// would be circular. This is what stops the constant drifting away from the
/// model that produced it.
#[test]
fn water_is_where_we_think_it_is() {
    let w = analyse(&water()).unwrap();
    println!("  water's polarity {:.3}, the constant says {WATER_POLARITY:.3}", w.polarity);
    assert!(
        (w.polarity - WATER_POLARITY).abs() < 0.05,
        "WATER_POLARITY is {WATER_POLARITY} but analysing water gives {:.3}",
        w.polarity
    );
}

/// Solubility is a property of a *pair*, and the general form has to agree
/// with the water-specific one — they used to be two models giving answers
/// sixty times apart for the same question.
#[test]
fn like_dissolves_like() {
    let w = analyse(&water()).unwrap();
    let salt = analyse(&rock_salt()).unwrap();
    let gas = analyse(&methane()).unwrap();
    let oil = analyse(&ethanol()).unwrap();

    let in_water = |p: &Properties| solubility_in(p, &w);
    println!(
        "  in water: salt {:.3e}, ethanol {:.3e}, methane {:.3e}",
        in_water(&salt),
        in_water(&oil),
        in_water(&gas)
    );
    println!("  methane in ethanol: {:.3e}", solubility_in(&gas, &oil));
    assert!(in_water(&salt) > in_water(&gas), "salt is more soluble in water than methane");
    // The convenience field and the general function are one model now, so
    // they must give the same answer. Not bit-identical: `water_solubility`
    // uses the `WATER_POLARITY` constant, since analysing water to find out
    // how soluble things are in water would be circular, while this route uses
    // the analysed value. The gap between them is the gap in the previous
    // assertion, and a hundredth of a decade is what that is worth.
    for (name, p) in [("salt", &salt), ("methane", &gas), ("ethanol", &oil)] {
        let a = p.water_solubility;
        let b = in_water(p);
        let decades = (a / b).log10().abs();
        println!("  {name:<8} water_solubility {a:.3e} against solubility_in {b:.3e} ({decades:.4} decades apart)");
        assert!(
            decades < 0.05,
            "{name}: water_solubility {a:.3e} against solubility_in {b:.3e}"
        );
    }
    assert!(
        solubility_in(&gas, &oil) > solubility_in(&gas, &w),
        "a non-polar gas prefers a less polar solvent"
    );
}

// ---------------------------------------------------------------------------
// legality: an arrangement that cannot exist gets no properties
// ---------------------------------------------------------------------------

#[test]
fn an_impossible_arrangement_is_refused_not_approximated() {
    // Five bonds on a hydrogen.
    let bad = Arrangement::molecule(
        vec![el(1), el(1), el(1), el(1), el(1), el(1)],
        (1..6).map(|i| Bond::new(0, i, Order::Single)).collect(),
    );
    match analyse(&bad) {
        Err(Illegal::OverBonded { element, used, allowed, .. }) => {
            println!("  refused: {element} holding {used} bonds with {allowed} slots");
            assert_eq!(used, 5);
            assert_eq!(allowed, 1);
        }
        other => panic!("a hydrogen with five bonds must be refused, got {other:?}"),
    }

    // Two separate molecules described as one substance.
    let split = Arrangement::molecule(
        vec![el(1), el(1), el(8), el(8)],
        vec![Bond::new(0, 1, Order::Single), Bond::new(2, 3, Order::Double)],
    );
    assert!(matches!(analyse(&split), Err(Illegal::Disconnected { pieces: 2 })));

    // Nothing at all.
    assert!(matches!(analyse(&Arrangement::molecule(vec![], vec![])), Err(Illegal::Empty)));

    // An element the table has no data for, rather than a guess.
    // Past the end of the table — 94 is plutonium and perfectly ordinary.
    let exotic = Arrangement::atom(Element(120), 0);
    assert!(matches!(analyse(&exotic), Err(Illegal::UnknownElement(_))));
}

// ---------------------------------------------------------------------------
// analysed once, then compiled
// ---------------------------------------------------------------------------

#[test]
fn the_same_substance_built_twice_is_one_substance() {
    let mut reg = Registry::new();
    let a = reg.intern(water()).unwrap();

    // The same molecule, written down in a different order and with its bonds
    // the other way round. Two players will not agree on atom ordering.
    let same = Arrangement::molecule(
        vec![el(1), el(8), el(1)],
        vec![Bond::new(1, 2, Order::Single), Bond::new(0, 1, Order::Single)],
    );
    let b = reg.intern(same).unwrap();

    println!("  {} substances after interning water twice, {} analyses, {} hits",
        reg.len(), reg.analyses, reg.hits);
    assert_eq!(a, b, "the same molecule written differently must be one substance");
    assert_eq!(reg.len(), 1);
    assert_eq!(reg.analyses, 1, "the second one must not have been analysed again");
    assert_eq!(reg.hits, 1);

    // And a different substance is a different substance.
    let c = reg.intern(methane()).unwrap();
    assert_ne!(a, c);
    assert_eq!(reg.len(), 2);
}

#[test]
fn a_caller_can_insist_on_a_separate_entry() {
    let mut reg = Registry::new();
    let a = reg.intern(water()).unwrap();
    let b = reg.intern_exact(water()).unwrap();
    assert_ne!(a, b, "intern_exact is the escape hatch for a recipe kept apart on purpose");
    assert_eq!(reg.len(), 2);
    // And it is still recognised as the same arrangement afterwards.
    assert!(reg.get(a).unwrap().arrangement.same_as(&reg.get(b).unwrap().arrangement));
}

#[test]
fn a_name_is_cosmetic() {
    let mut reg = Registry::new();
    let id = reg.intern(water()).unwrap();
    let before = reg.get(id).unwrap().props;
    reg.name(id, "water");
    assert_eq!(reg.by_name("water"), Some(id));
    assert_eq!(reg.get(id).unwrap().props, before, "naming must not change any physics");
    // And a substance nobody named behaves identically.
    let mut other = Registry::new();
    let anon = other.intern(water()).unwrap();
    assert_eq!(other.get(anon).unwrap().props, before);
}

#[test]
fn a_measurement_overrides_but_does_not_erase() {
    let mut reg = Registry::new();
    let id = reg.intern(water()).unwrap();
    let derived = reg.get(id).unwrap().props;
    assert_eq!(reg.get(id).unwrap().provenance, Provenance::Derived);

    let mut measured = derived;
    measured.boiling_point = 373.15;
    measured.melting_point = 273.15;
    assert!(reg.measured(id, measured));

    let s = reg.get(id).unwrap();
    assert_eq!(s.props.boiling_point, 373.15);
    match &s.provenance {
        Provenance::Measured { replaced } => {
            println!(
                "  derived {:.1} K, measured {:.1} K — the disagreement stays visible",
                replaced.boiling_point, s.props.boiling_point
            );
            assert_eq!(replaced.boiling_point, derived.boiling_point);
        }
        other => panic!("expected a measurement record, got {other:?}"),
    }
}

/// Something the engine has never heard of, analysed anyway.
#[test]
fn an_unknown_molecule_still_gets_properties() {
    let mut reg = Registry::new();
    let id = reg.intern(glycine()).expect("a novel molecule must still analyse");
    let s = reg.get(id).unwrap();
    println!(
        "  {} : M {:.2} g/mol, {} atoms, boils at {:.0} K, solubility {:.2e} kg/kg [{:?}]",
        s.formula.hill(),
        s.props.molar_mass * 1000.0,
        s.arrangement.atoms.len(),
        s.props.boiling_point,
        s.props.water_solubility,
        s.props.confidence
    );
    assert!(s.props.molar_mass > 0.0);
    assert!(s.label.is_none(), "nothing named it and nothing needed to");
}

// ---------------------------------------------------------------------------
// mixtures, and the two accounts agreeing
// ---------------------------------------------------------------------------

#[test]
fn a_mixture_reconciles_with_the_elemental_account() {
    let mut reg = Registry::new();
    let salt = reg.intern(rock_salt()).unwrap();
    let h2o = reg.intern(water()).unwrap();

    // Seawater, roughly: 3.5% salt by mass.
    let mut brine = Mixture::new();
    assert!(brine.add(salt, 0.035));
    assert!(brine.add(h2o, 0.965));

    let (comp, explained) = brine.composition(&reg);
    println!(
        "  brine explains {:.1}% of the mass; as species H {:.4} O {:.4} Other {:.4}",
        explained * 100.0,
        comp.get(phys::units::Species::Hydrogen),
        comp.get(phys::units::Species::Oxygen),
        comp.get(phys::units::Species::Other),
    );
    assert!((explained - 1.0).abs() < 1e-12);
    // Sodium and chlorine both lump into `Other`, and their combined mass
    // fraction has to come back out.
    assert!((comp.get(phys::units::Species::Other) - 0.035).abs() < 1e-9);
    let sum: f64 = phys::units::Species::ALL.iter().map(|s| comp.get(*s)).sum();
    assert!((sum - 1.0).abs() < 1e-12, "a composition must sum to one, got {sum}");

    // The molar mass the lumped account cannot give.
    let m = brine.molar_mass(&reg).unwrap();
    println!("  mean molar mass {:.2} g/mol", m * 1000.0);
    assert!(m > 0.0);
}

#[test]
fn similar_recipes_are_one_recipe() {
    let mut reg = Registry::new();
    let salt = reg.intern(rock_salt()).unwrap();
    let h2o = reg.intern(water()).unwrap();

    let brew = |s: f64| {
        let mut m = Mixture::new();
        m.add(salt, s);
        m.add(h2o, 1.0 - s);
        m
    };
    let a = brew(0.0350);
    let b = brew(0.0351);
    let c = brew(0.0500);

    println!(
        "  3.50% and 3.51% same: {}; 3.50% and 5.00% same: {}",
        a.same_as(&b),
        a.same_as(&c)
    );
    assert!(a.same_as(&b), "a difference under the tolerance is not a different recipe");
    assert_eq!(a.fingerprint(), b.fingerprint());
    assert!(!a.same_as(&c), "a real difference must survive");
    assert_ne!(a.fingerprint(), c.fingerprint());

    // Order of assembly must not matter either.
    let mut reversed = Mixture::new();
    reversed.add(h2o, 0.965);
    reversed.add(salt, 0.035);
    assert!(a.same_as(&reversed));
}

#[test]
fn a_mixture_is_bounded_and_says_what_it_dropped() {
    let mut reg = Registry::new();
    let mut m = Mixture::new();
    let mut ids = Vec::new();
    for i in 0..12 {
        // Twelve distinguishable substances: chains of carbon of different
        // lengths, none of which anything has ever heard of.
        let n = i + 2;
        let atoms = vec![el(6); n];
        let bonds = (0..n - 1).map(|k| Bond::new(k, k + 1, Order::Single)).collect();
        ids.push(reg.intern(Arrangement::molecule(atoms, bonds)).unwrap());
    }
    for (k, id) in ids.iter().enumerate() {
        m.add(*id, 0.01 * (k + 1) as f64);
    }
    println!(
        "  {} substances offered, {} kept, {:.3} of the mass speciated",
        ids.len(),
        m.len(),
        m.speciated()
    );
    assert_eq!(m.len(), phys::chem::registry::MIXTURE_SLOTS);
    // A trace species that does not make the cut is told so.
    let tiny = ids[0];
    assert!(!m.add(tiny, 1e-9), "a trace below everything held must be refused, not silently lost");
}

// ---------------------------------------------------------------------------
// the table itself
// ---------------------------------------------------------------------------

/// Hydrogen to plutonium, with nothing missing and nothing invented.
#[test]
fn the_table_runs_to_plutonium() {
    assert_eq!(Element::from_symbol("Pu"), Some(Element(94)));
    assert_eq!(Element(94).symbol(), "Pu");
    assert!(!Element(95).known(), "the table must stop where the data does");

    // Every element in range answers every question.
    for z in 1..=elements::HEAVIEST {
        let e = Element(z);
        assert!(e.known(), "Z={z} is in range and must be known");
        for (what, present) in [
            ("weight", e.weight().is_some()),
            ("electronegativity", e.electronegativity().is_some()),
            ("covalent radius", e.covalent_radius().is_some()),
            ("van der Waals radius", e.vdw_radius().is_some()),
            ("valence", e.valence().is_some()),
            ("valence electrons", e.valence_electrons().is_some()),
            ("bond energy", e.homonuclear_bond().is_some()),
        ] {
            assert!(present, "{} (Z={z}) has no {what}", e.symbol());
        }
        assert!(e.weight().unwrap() > 0.0, "{} has no mass", e.symbol());
        assert!(e.vdw_radius().unwrap() > e.covalent_radius().unwrap(),
            "{}: a van der Waals radius must exceed the covalent one", e.symbol());
    }

    // Symbols are unique, or `from_symbol` would be ambiguous.
    let mut seen = std::collections::BTreeSet::new();
    for z in 1..=elements::HEAVIEST {
        assert!(seen.insert(Element(z).symbol()), "duplicate symbol at Z={z}");
    }
    println!("  {} elements, H to Pu, complete", elements::HEAVIEST);
}

/// Atomic weights against the values a chemist would look up.
#[test]
fn atomic_weights_are_right_across_the_table() {
    for (sym, real) in [
        ("H", 1.008), ("C", 12.011), ("O", 15.999), ("Na", 22.990), ("Cl", 35.45),
        ("Fe", 55.845), ("Ag", 107.87), ("I", 126.90), ("W", 183.84), ("Au", 196.97),
        ("Pb", 207.2), ("U", 238.03), ("Pu", 244.0),
    ] {
        let e = Element::from_symbol(sym).unwrap_or_else(|| panic!("{sym} missing"));
        let got = e.weight().unwrap();
        println!("  {sym:<3} Z={:<3} {got:>8.3} u", e.z());
        assert!((got - real).abs() < 0.01, "{sym}: {got} against {real}");
    }
}

/// The heavy elements are usable, not merely present: a uranium compound has to
/// analyse like any other.
#[test]
fn a_heavy_element_analyses_like_any_other() {
    let u = Element::from_symbol("U").unwrap();
    let o = Element::from_symbol("O").unwrap();
    // Uranium dioxide, the ceramic reactor fuel is made of. Fluorite lattice,
    // cell edge 547 pm holding four formula units.
    let uo2 = Arrangement::crystal(
        vec![u, o, o],
        vec![Bond::new(0, 1, Order::Ionic), Bond::new(0, 2, Order::Ionic)],
        Lattice::Cubic { a: 5.47e-10 / 4f64.cbrt() },
    );
    let p = analyse(&uo2).expect("uranium dioxide must analyse");
    println!(
        "  {} : M {:.2} g/mol (real 270.03), density {:.0} kg/m3 (real 10970), \
         ionicity {:.2}, melts at {:.0} K (real 3138)",
        uo2.formula().hill(),
        p.molar_mass * 1000.0,
        p.density,
        p.ionicity,
        p.melting_point
    );
    assert!((p.molar_mass * 1000.0 - 270.03).abs() < 0.05, "molar mass is exact or nothing is");
    // Density comes from the lattice, so it should be close.
    assert!(
        (p.density / 10970.0 - 1.0).abs() < 0.15,
        "density {:.0} against a real 10970",
        p.density
    );
    assert!(p.ionicity > 0.4, "a metal oxide should be substantially ionic");
    // And it is not soluble in water, which is what makes it a fuel pellet
    // rather than a hazard the moment it rains. A multiply-charged lattice is
    // outside the range the solubility model is calibrated on, so what is
    // asserted here is the class — negligible — and not the digits.
    println!("  water solubility {:.3e} kg/kg (real ~1e-9: negligible)", p.water_solubility);
    assert!(p.water_solubility < 1e-6, "a fuel pellet must not dissolve in rain");
}

/// Plutonium metal: the case the table was extended for.
#[test]
fn plutonium_is_a_substance_like_any_other() {
    let mut reg = Registry::new();
    let pu = Element::from_symbol("Pu").unwrap();
    let metal = Arrangement::crystal(
        vec![pu],
        Vec::new(),
        // Delta-phase plutonium: face-centred cubic, 463 pm, four atoms a cell.
        Lattice::Cubic { a: 4.63e-10 / 4f64.cbrt() },
    );
    let id = reg.intern(metal).expect("plutonium must analyse");
    let s = reg.get(id).unwrap();
    println!(
        "  {} : M {:.1} g/mol (real 244.0), density {:.0} kg/m3 (real 15920)",
        s.formula.hill(),
        s.props.molar_mass * 1000.0,
        s.props.density
    );
    assert!((s.props.molar_mass * 1000.0 - 244.0).abs() < 0.1);
    assert!((s.props.density / 15920.0 - 1.0).abs() < 0.15);
}
