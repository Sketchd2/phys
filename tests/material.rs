//! What a thing is made of, measured rather than looked up.
//!
//! `docs/PLAY.md` D13 and D14. A material used to be thirteen numbers somebody
//! typed, reached by asking what *generated* a node — `topology.material` for a
//! structure, `morphology.material()` for one not yet materialised, and `None`
//! for a rock. §2A measured what that cost: a boulder had no surface and could
//! not collide, which is a provenance test standing in for a state measurement.
//!
//! Now every number comes out of `chem::Properties`, which comes out of an
//! `Arrangement` — atoms and bonds — and the one thing that is not chemistry is
//! the **flaw scale**, which comes out of how the thing was made.
//!
//! These tests are therefore mostly *measurements*, in the style of
//! `tests/chem.rs`: each prints the derived value against the one the retired
//! table held and holds a tolerance that says how good the model actually is.

use phys::chem::{Mixture, Phase, Registry};
use phys::material::{substances, Formation, Material};

/// The eight stresses the retired `rupture` column held.
///
/// **This is the only place they survive**, and the derivation never sees them:
/// nothing in `src/` reads this list, and `Material::of` is given an
/// arrangement and a formation history and nothing else. `docs/PLAY.md` D14 is
/// explicit that reproducing them *exactly* would mean the derivation had been
/// fitted to them, so what is asserted below is the **ordering** and the
/// **order of magnitude**.
const RETIRED: &[(&str, f64)] = &[
    ("steel", 400.0e6),
    ("reinforced frame", 180.0e6),
    ("bedrock", 1.3e8),
    ("dry timber", 70.0e6),
    ("green wood", 45.0e6),
    ("aragonite", 12.0e6),
    ("masonry", 2.0e6),
    ("ice", 1.7e6),
];

fn retired(name: &str) -> f64 {
    RETIRED.iter().find(|(n, _)| *n == name).map(|(_, v)| *v).expect("a named preset")
}

/// A piece of a stated size, at the temperature it is used at.
///
/// Two things a bare `strength()` cannot express and the retired table never
/// distinguished. A table of stresses is implicitly *for a piece*, at
/// laboratory size and laboratory temperature, and the derivation now says
/// which — a metre across, at 290 K, except for a material that is not a solid
/// at 290 K, which is compared at 60% of its own melting point instead. That is
/// a rule rather than a column: it takes the same branch for all eight, and
/// only ice goes the other way, because only ice melts below room temperature.
fn service(m: &Material) -> f64 {
    290.0f64.min(0.6 * m.thermal_gone)
}

/// Every preset lands within a factor of three of the stress it replaces.
///
/// The weaker half of D14's honesty test, and the one that catches a derivation
/// that has gone somewhere silly. **Bedrock is no longer excepted.** It used to
/// be, by name and with the reason measured — its grain comes out at 21 cm and
/// a granite at 130 MPa implies a 34 um crack — and closing that is Phase 4's
/// first inherited measurement. See `Material::worst_flaw`.
#[test]
fn the_derived_strengths_are_the_right_size() {
    let mut worst = (0.0f64, String::new());
    for m in Material::presets() {
        let want = retired(m.name);
        let t = service(m);
        let got = m.strength_of(1.0, t);
        let ratio = got / want;
        println!(
            "  {:<18} r* {:>9.3e}  grain {:>9.3e}  worst flaw {:>9.3e} m   \
             derived {:>9.3e} Pa   retired {:>9.3e} Pa   ratio {ratio:>6.3}",
            m.name,
            m.nucleus,
            m.flaw_size,
            m.worst_flaw_at(1.0, t).min(1.0).min(m.flaw_size),
            got,
            want
        );
        let off = if ratio > 1.0 { ratio } else { 1.0 / ratio };
        if off > worst.0 {
            worst = (off, m.name.to_string());
        }
    }
    assert!(
        worst.0 < 3.0,
        "{} is off by {:.1}x, and the whole table is within 2.7x",
        worst.1,
        worst.0
    );
}

/// The same rock is stronger as a pebble than as a cliff face.
///
/// `docs/PLAY.md` Phase 4's third done-when. Two effects that point opposite
/// ways and are the same sentence about where the worst crack is: a small piece
/// cannot contain a large flaw, and a large one has more chances of containing
/// one. Weibull's size effect falls out with modulus `2b = 12`, which is what
/// fixes the exponent — see `Material::worst_flaw`.
#[test]
fn a_boulder_and_a_cliff_face_of_the_same_rock_differ_by_size() {
    let rock = Material::bedrock();
    let mut last = f64::INFINITY;
    let mut sizes = Vec::new();
    for size in [0.01, 0.1, 1.0, 10.0, 100.0] {
        let s = rock.strength_of(size, 290.0);
        println!(
            "  bedrock {size:>7} m across: worst flaw {:>9.3e} m, strength {s:>9.4e} Pa",
            rock.worst_flaw(size * size * size).min(size).min(rock.flaw_size)
        );
        assert!(s < last, "a bigger piece of the same rock is weaker");
        last = s;
        sizes.push(s);
    }
    // Weibull with modulus 2b: a decade of size is a factor of 10^(3/2b) in
    // strength, which for b = 6 is 1.778. Asserted as the law rather than as
    // the five numbers, so a change to the exponent fails here and is read.
    let decade = 10.0f64.powf(3.0 / (2.0 * Material::flaw_exponent()));
    for pair in sizes.windows(2) {
        assert!(
            (pair[0] / pair[1] - decade).abs() < 1e-9 * decade,
            "a decade of size should be {decade:.4}x of strength, and is {:.4}",
            pair[0] / pair[1]
        );
    }
    // And a boulder is not a cliff.
    println!(
        "  a 1 m boulder is {:.2}x the strength of a 100 m cliff face",
        rock.strength_of(1.0, 290.0) / rock.strength_of(100.0, 290.0)
    );
    assert!(rock.strength_of(1.0, 290.0) / rock.strength_of(100.0, 290.0) > 3.0);
}

/// The ordering survives, and it was never shown the values it reproduces.
///
/// The stronger half of D14's honesty test. Sorted by derived strength at each
/// material's service condition, the eight come out in the retired table's own
/// order except for two adjacent pairs, and both are recorded rather than
/// tolerated:
///
/// - **steel and a reinforced frame swap.** The cause is measured and is not in
///   this law: a metal's cohesive energy derives about 2.3x low, because iron
///   really has eight nearest neighbours and the valence model allows three, so
///   steel's melting point comes out 596 K against 1811 and it is *already
///   creeping* at 290 K. `docs/BACKLOG.md` carries the root cause.
/// - **ice and masonry swap**, as they did before, and they are adjacent: the
///   table separates them by 18%, which no derivation from first principles
///   should be expected to resolve.
///
/// **Bedrock is third, which is where the table puts it.** It used to derive as
/// the weakest of the eight, and the assertion below used to be written in that
/// direction with a note saying a future fix should fail here and be read. This
/// is that reading.
#[test]
fn the_ordering_survives_without_the_table() {
    let mut derived: Vec<(&str, f64)> = Material::presets()
        .iter()
        .map(|m| (m.name, m.strength_of(1.0, service(m))))
        .collect();
    derived.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("finite strengths"));
    for (i, (name, s)) in derived.iter().enumerate() {
        println!("  {}. {name:<18} {s:.3e} Pa", i + 1);
    }
    let order: Vec<&str> = derived.iter().map(|(n, _)| *n).collect();

    // The pairs the derivation is held to, each one a claim about physics
    // rather than about a table.
    let rank = |name: &str| order.iter().position(|n| *n == name).expect("present");
    assert!(
        rank("reinforced frame") < rank("bedrock"),
        "a reinforced frame carries more than bedrock"
    );
    assert!(rank("bedrock") < rank("dry timber"), "rock carries more than timber");
    assert!(
        rank("dry timber") < rank("green wood"),
        "slow-grown seasoned timber is stronger than fast-grown green wood, and \
         the only thing that differs between them is how thick a layer the tree \
         laid down"
    );
    assert!(rank("green wood") < rank("aragonite"), "timber beats a coral skeleton");
    assert!(rank("aragonite") < rank("masonry"), "a coral skeleton beats masonry");
    assert!(
        rank("steel") < rank("dry timber"),
        "steel carries more than timber even while it is creeping"
    );

    // And the two the table asserts that this does not reproduce, asserted in
    // the direction they actually come out, so that a future fix *fails* here
    // and gets read rather than quietly changing the answer.
    assert!(
        rank("reinforced frame") < rank("steel"),
        "steel derives below a reinforced frame because its melting point \
         derives 3x low and it is creeping at room temperature; if that has \
         changed, the valence model has been fixed and this test is what should \
         say so"
    );
    assert!(
        rank("ice") < rank("masonry"),
        "ice and masonry swap, and they are 18% apart in the table"
    );
}

/// The same substance, two histories, two materials.
///
/// D14's central claim and the reason a *stress* could never have expressed it:
/// masonry and bedrock are both silicate and green wood and dry timber are both
/// cellulose, so the retired table recorded four materials where there are two
/// substances. Nothing here says anything about silicate or cellulose that is
/// not in their atoms; everything that separates each pair is in the second
/// argument.
#[test]
fn one_substance_and_two_histories_make_two_materials() {
    // Silicate: fired and laid in courses, against crystallised in the crust.
    let brick = Material::masonry();
    let rock = Material::bedrock();
    println!(
        "  silicate — masonry a {:.3e} m -> {:.3e} Pa, bedrock a {:.3e} m -> {:.3e} Pa",
        brick.flaw_size,
        brick.strength(),
        rock.flaw_size,
        rock.strength()
    );
    // Surface energy is a property of the bonds alone, so it is identical:
    // both are silicate, and nothing about how they were made changes what it
    // costs to part one.
    assert_eq!(
        brick.surface_energy, rock.surface_energy,
        "the same substance costs the same to part"
    );
    // Stiffness is not, and the difference is entirely how much of the volume
    // each one fills: mortar leaves voids and a pluton does not. Gibson and
    // Ashby's square.
    let packing = 1900.0 / 2644.0;
    assert!(
        (brick.stiffness / rock.stiffness - packing * packing).abs() < 0.02,
        "masonry is {:.3} of bedrock's stiffness where its packing squared is {:.3}",
        brick.stiffness / rock.stiffness,
        packing * packing
    );
    assert!(
        brick.flaw_size != rock.flaw_size,
        "two histories should not give one flaw scale"
    );

    // Cellulose: fast summer growth against slow growth, seasoned.
    let green = Material::green_wood();
    let dry = Material::dry_timber();
    println!(
        "  cellulose — green a {:.3e} m -> {:.3e} Pa, dry a {:.3e} m -> {:.3e} Pa",
        green.flaw_size,
        green.strength(),
        dry.flaw_size,
        dry.strength()
    );
    assert!(dry.strength() > green.strength(), "slow growth makes stronger wood");
    assert!(dry.density < green.density, "and drier wood is lighter");
}

/// Strength is size-dependent, which a tabulated stress cannot be.
///
/// D14: "a thin fibre is genuinely stronger than a thick bar of the same
/// material, because a small piece cannot contain a large flaw. That is real,
/// measurable, currently inexpressible, and free — every member already carries
/// a radius."
#[test]
fn a_small_piece_cannot_hold_a_large_flaw() {
    let m = Material::masonry();
    let a = m.flaw_size;
    let big = m.strength_at_size(a * 100.0);
    let small = m.strength_at_size(a / 100.0);
    println!(
        "  masonry flaw {a:.3e} m: a piece 100x larger fails at {big:.3e} Pa, one \
         100x smaller at {small:.3e} Pa"
    );
    assert_eq!(big, m.strength(), "a piece bigger than the flaw sees the whole flaw");
    // Griffith goes as 1/sqrt(a), so a hundredfold smaller flaw is ten times
    // the stress.
    assert!(
        (small / big - 10.0).abs() < 0.2,
        "a hundredth of the flaw should be ten times the strength, and gave {:.3}",
        small / big
    );
}

/// A rock that no `Program` made has a material, which is the point of D13.
///
/// The measurement §2A made was that `surface_of` returned `None` for anything
/// without a `Topology` or a `Morphology`, so a boulder could not collide. A
/// node that carries a mixture now answers from it.
#[test]
fn matter_that_nothing_generated_still_has_a_material() {
    use phys::state::{Composition, Matter};

    let mut reg = Registry::new();
    let silica = reg.intern(substances::silica_arrangement()).expect("silica analyses");
    let mut mix = Mixture::new();
    mix.add(silica, Phase::Solid, 1.0);

    let matter = Matter::neutral(2650.0, 0.5, 290.0, Composition::primordial());
    let formation = Formation::of_matter(&matter, &reg.get(silica).unwrap().props);
    let m = Material::measured(&mix, &reg, formation).expect("a solid has a material");
    println!(
        "  a boulder of silicate: rho {:.1} kg/m^3, E {:.3e} Pa, a {:.3e} m, sigma {:.3e} Pa",
        m.density,
        m.stiffness,
        m.flaw_size,
        m.strength()
    );
    assert!(m.density > 0.0 && m.stiffness > 0.0 && m.strength() > 0.0);

    // And matter with no solid in it has none, which is D13's own line: a
    // liquid's surface belongs to its container and a gas has none at all.
    let mut liquid = Mixture::new();
    liquid.add(silica, Phase::Liquid, 1.0);
    assert!(
        Material::measured(&liquid, &reg, formation).is_none(),
        "a melt is not a solid and has no surface of its own"
    );
    assert!(
        Material::measured(&Mixture::new(), &reg, formation).is_none(),
        "matter nobody has described has no material either"
    );
}

/// Densities come out right, which is the part with no free parameter at all.
///
/// Nothing in the flaw-scale argument touches density: it is the substance's
/// own, from its unit cell or its van der Waals volume, times how much of the
/// space the process filled. Four of the eight state a packing and come out
/// exact by construction; the other four state nothing and are the measurement.
#[test]
fn the_densities_are_derived_and_land() {
    let expected: &[(&str, f64, f64)] = &[
        // name, what it should be, how far off it may be
        ("steel", 7850.0, 0.05),
        ("bedrock", 2700.0, 0.10),
        ("ice", 917.0, 1.0),
        ("reinforced frame", 2400.0, 0.40),
    ];
    for (name, want, tolerance) in expected {
        let m = Material::presets().iter().find(|m| m.name == *name).expect("a preset");
        let off = (m.density - want).abs() / want;
        println!("  {name:<18} {:>8.1} kg/m^3 against {want:>8.1}, {:.1}% off", m.density, off * 100.0);
        assert!(
            off <= *tolerance,
            "{name} derives at {:.1} kg/m^3 against {want}",
            m.density
        );
    }
}
