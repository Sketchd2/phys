//! The engine's central claim, stated as executable assertions.
//!
//! If these pass, no experiment performed at a coarse scale can detect that the
//! fine detail was generated rather than simulated.

use phys::math::v3;
use phys::prolong::*;
use phys::state::*;
use phys::units::*;

fn sample_aggregates() -> Vec<(&'static str, Aggregate)> {
    let mut out = Vec::new();

    let mut cloud = Aggregate::neutral(3.0e4 * M_SUN, 12.0 * PARSEC, 20.0, Composition::solar());
    cloud.momentum = v3(1e33, -2e32, 5e32);
    cloud.spin = v3(0.0, 0.0, 4e50);
    out.push(("molecular cloud", cloud));

    let mut star = Aggregate::neutral(M_SUN, R_SUN, 5.8e6, Composition::solar());
    star.spin = v3(1e41, 0.0, 1.1e42);
    star.binding_energy = -0.6 * G * M_SUN * M_SUN / R_SUN;
    out.push(("star", star));

    let mut planet = Aggregate::neutral(M_EARTH, R_EARTH, 3000.0, Composition::pure(Species::Silicon));
    planet.spin = v3(0.0, 0.0, 7.05e33);
    out.push(("planet", planet));

    let grain = Aggregate::neutral(1e-9, 1e-5, 150.0, Composition::pure(Species::Carbon));
    out.push(("dust grain", grain));

    let ion = Aggregate::neutral(56.0 * AMU, 4.6e-15, 1e7, Composition::pure(Species::Iron))
        .with_charge(26.0 * E_CHARGE);
    assert_eq!(ion.validate(), 0.0, "test fixture must be self-consistent");
    out.push(("iron nucleus (fully stripped)", ion));

    let hot = Aggregate::neutral(1e-20, 1e-9, 1e6, Composition::primordial());
    out.push(("hot plasma parcel", hot));

    out
}

fn specs() -> Vec<(&'static str, ProlongSpec)> {
    vec![
        ("uniform/equal", ProlongSpec::new(1000, Profile::Uniform, MassSpectrum::Equal, BodyKind::GasParcel)),
        ("plummer/kroupa", ProlongSpec::new(1000, Profile::Plummer, MassSpectrum::Kroupa { min_msun: 0.08, max_msun: 60.0 }, BodyKind::Star)),
        ("disk/powerlaw", ProlongSpec { count: 2000, profile: Profile::Disk { scale_height_ratio: 0.1 }, spectrum: MassSpectrum::PowerLaw { alpha: -1.8, ratio: 50.0 }, kind: BodyKind::Super, composition_scatter: 0.2, turbulent_fraction: 0.5 }),
        ("shell/equal", ProlongSpec::new(500, Profile::Shell, MassSpectrum::Equal, BodyKind::GasParcel)),
        ("woods-saxon", ProlongSpec::new(56, Profile::WoodsSaxon, MassSpectrum::Equal, BodyKind::Nucleon)),
        ("lattice/species", ProlongSpec::new(64, Profile::Lattice, MassSpectrum::Species, BodyKind::Atom)),
        ("tiny", ProlongSpec::new(2, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain)),
    ]
}

/// The invariant: restriction after prolongation returns the original state.
#[test]
fn round_trip_conserves_everything() {
    let mut worst = 0.0f64;
    let mut worst_case = String::new();
    for (aname, agg) in sample_aggregates() {
        for (sname, spec) in specs() {
            for seed in [1u64, 0xDEAD_BEEF, 0x5EED_5EED] {
                let (bodies, report) = prolong(&agg, spec, seed, 0xABCD_1234, 0);
                assert!(!bodies.is_empty(), "{aname}/{sname} produced nothing");
                let mut back = restrict(&bodies, report.potential);
                back.external_potential = agg.external_potential;
                let scales = Scales::of(&bodies);
                let err = back.conserved().error_against(&agg.conserved(), &scales);
                if err > worst {
                    worst = err;
                    worst_case = format!("{aname} / {sname} / seed {seed:x}");
                }
                assert!(
                    err < 1e-9,
                    "{aname}/{sname} seed {seed:x}: conservation error {err:.3e}"
                );
                    assert!(back.mass > 0.0 && back.is_finite_state());
                assert!(
                    back.validate() < 1e-9,
                    "{aname}/{sname}: restricted state inconsistent ({:.3e})",
                    back.validate()
                );
            }
        }
    }
    println!("worst conservation error {worst:.3e} at {worst_case}");
}

/// Mass, charge, baryon and lepton number are exactly additive, with no
/// tolerance for scale-relative fudging: these are counted quantities.
#[test]
fn counted_quantities_are_exact() {
    for (name, agg) in sample_aggregates() {
        for (sname, spec) in specs() {
            let (bodies, report) = prolong(&agg, spec, 7, 0x1111, 0);
            let back = restrict(&bodies, report.potential);
            let dm = (back.mass - agg.mass).abs() / agg.mass;
            assert!(dm < 1e-14, "{name}/{sname}: mass drift {dm:.3e}");
            let dq = (back.charge - agg.charge).abs() / agg.charge.abs().max(E_CHARGE);
            assert!(dq < 1e-10, "{name}/{sname}: charge drift {dq:.3e}");
            let db = (back.baryon_number - agg.baryon_number).abs()
                / agg.baryon_number.abs().max(1.0);
            assert!(db < 1e-10, "{name}/{sname}: baryon drift {db:.3e}");
        }
    }
}

/// Regeneration is a pure function of the address. Same seed, same key, same
/// epoch — same bodies, bit for bit, however many times you ask.
#[test]
fn regeneration_is_bit_identical() {
    let agg = sample_aggregates()[0].1;
    let spec = specs()[2].1;
    let a = prolong(&agg, spec, 99, 0xFEED_FACE, 3).0;
    let b = prolong(&agg, spec, 99, 0xFEED_FACE, 3).0;
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(x.pos, y.pos);
        assert_eq!(x.vel, y.vel);
        assert_eq!(x.mass, y.mass);
        assert_eq!(x.charge, y.charge);
        assert_eq!(x.spin, y.spin);
    }
}

/// Different addresses must give genuinely different worlds — otherwise the
/// galaxy is a tiling of one cloud.
#[test]
fn different_addresses_decorrelate() {
    let agg = sample_aggregates()[0].1;
    let spec = specs()[0].1;
    let a = prolong(&agg, spec, 1, 0x1, 0).0;
    let b = prolong(&agg, spec, 1, 0x2, 0).0;
    let c = prolong(&agg, spec, 2, 0x1, 0).0;
    let d = prolong(&agg, spec, 1, 0x1, 1).0;
    let differs = |x: &Vec<Body>, y: &Vec<Body>| x.iter().zip(y).filter(|(p, q)| p.pos != q.pos).count();
    assert!(differs(&a, &b) > a.len() * 9 / 10, "path key must decorrelate");
    assert!(differs(&a, &c) > a.len() * 9 / 10, "world seed must decorrelate");
    assert!(differs(&a, &d) > a.len() * 9 / 10, "epoch must decorrelate");
}

/// A tree walked down and back up returns to exactly where it started.
#[test]
fn tree_round_trip_is_idempotent() {
    use phys::engine::{default_spec, galaxy};
    use phys::tree::Tree;
    let mut tree: Tree = galaxy(0xC0FFEE, 1e9);
    let root = tree.root;
    let before_agg = tree.nodes[root.get()].agg;
    let first: Vec<Body> = tree.refine(root).to_vec();
    let err = tree.coarsen(root);
    let after_agg = tree.nodes[root.get()].agg;
    assert!(err < 1e-9, "coarsening error {err:.3e}");
    assert_eq!(before_agg.mass, after_agg.mass, "mass must be untouched");
    assert_eq!(
        before_agg.total_energy(),
        after_agg.total_energy(),
        "energy must be untouched"
    );
    let second: Vec<Body> = tree.refine(root).to_vec();
    assert_eq!(first.len(), second.len());
    for (a, b) in first.iter().zip(&second) {
        assert_eq!(a.pos, b.pos, "regenerated body moved");
        assert_eq!(a.vel, b.vel);
    }
    let _ = default_spec(Tier::Galactic);
}

/// Descending 20+ levels must not accumulate error.
#[test]
fn deep_descent_stays_exact() {
    use phys::engine::{default_spec, galaxy, World};
    let mut w = World::new(galaxy(0x5EED, 1e9), 20.0);
    let root = w.tree.root;
    let path = w.drill(root, Tier::Nuclear, &default_spec);
    assert!(path.len() >= 15, "expected a deep descent, got {}", path.len());
    for &idx in &path {
        let n = &w.tree.nodes[idx.get()];
        assert!(
            n.last_report.conservation_error < 1e-9,
            "{} at depth {}: error {:.3e}",
            n.tier.name(),
            n.depth,
            n.last_report.conservation_error
        );
        assert!(n.agg.mass > 0.0 && n.agg.mass.is_finite());
        assert!(n.agg.radius > 0.0 && n.agg.radius.is_finite());
    }
    // The descent must actually span the scales it claims to.
    let top = w.tree.nodes[path[0].get()].agg.radius;
    let bottom = w.tree.nodes[path[path.len() - 1].get()].agg.radius;
    assert!(
        top / bottom > 1e30,
        "descent spanned only {:.1e} in scale",
        top / bottom
    );
    assert!(w.tree.stats.worst_conservation_error < 1e-9);
}

trait FiniteState {
    fn is_finite_state(&self) -> bool;
}
impl FiniteState for Aggregate {
    fn is_finite_state(&self) -> bool {
        self.is_finite() && self.total_energy().is_finite()
    }
}
