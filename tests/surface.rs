//! What a solid presents to the world. `docs/PLAY.md` D18.
//!
//! §2A's finding was that **the engine has no representation for the shape of a
//! solid, at any scale**. What existed was a point and a radius for a body, one
//! scalar for a whole node, and a member list that was a *structural analysis*
//! artefact — it exists so the stress solver knows where bending peaks, and
//! collision was improvised out of it.
//!
//! D18 replaces the improvisation with a stated boundary: a union of solid
//! convex primitives, each carrying its own material, emitted by whatever
//! generated the thing. These tests hold it to the three claims that make it
//! work — the primitives are solid, the generator states them, and they
//! reconcile with what the node is made of.

use phys::chem::{Mixture, Phase};
use phys::engine::World;
use phys::math::{v3, Quat, Vec3};
use phys::morph::NO_SUPPORT;
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::topology::{Joint, Material, Topology};
use phys::tree::Tree;
use phys::units::Tier;

/// A node of wood, with nothing else said about it.
fn wooden(seed: u64, mass: f64, radius: f64, count: usize) -> (World, phys::ids::NodeIdx) {
    let spec = SampleSpec::new(count, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let matter = Matter::neutral(mass, radius, 291.0, Composition::primordial());
    let mut w = World::new(Tree::new(seed, matter, Tier::Continuum, spec), 20.0);
    let root = w.tree.root;
    let cellulose = w
        .substances
        .intern(phys::material::substances::cellulose_arrangement())
        .expect("cellulose analyses");
    let mut mix = Mixture::new();
    mix.add(cellulose, Phase::Solid, 1.0);
    w.set_mixture(root, mix);
    (w, root)
}

/// A rock that no `Program` made presents a surface.
///
/// §2A measured the old behaviour: `surface_of` read `topology.material` or
/// `morphology.material()` and returned `None` for everything else, so a rock,
/// a boulder and a ball of wood had no surface at all. That is a *provenance
/// test standing in for a state measurement*, and §3.3 had already made the
/// same call correctly one layer down.
#[test]
fn a_thing_nobody_built_presents_a_surface() {
    let (mut w, root) = wooden(0x5011D, 900.0, 2.0, 64);
    let s = w.surface_of_node(root).clone();
    println!(
        "  an unmaterialised ball of wood presents {} piece(s), bound {:.3} m",
        s.len(),
        s.pieces().first().map(|p| p.hull.bound()).unwrap_or(0.0)
    );
    assert_eq!(s.len(), 1, "a rock is one filled solid, not a cloud of them");
    let p = &s.pieces()[0];
    assert!(p.material.density > 0.0 && p.material.strength() > 0.0);
    assert!(
        (p.hull.bound() - 2.0).abs() < 1e-9,
        "unmaterialised, the piece is the sphere the matter claims: {}",
        p.hull.bound()
    );
}

/// Matter that is not solid presents nothing, and that is measured.
///
/// D13 is explicit that it is a decision about solids and only about solids: a
/// liquid's surface is a property of its container and the field it is in, so
/// it cannot be derived once and kept, and a cloud of gas has none for the same
/// reason it has no shape. The worked case is a cup of water — the ceramic has
/// a boundary, the water in it has phase fractions and no geometry at all.
#[test]
fn a_liquid_presents_no_surface_and_says_so() {
    let (mut w, root) = wooden(0x1191D, 900.0, 2.0, 64);
    assert!(!w.surface_of_node(root).is_empty(), "precondition: solid wood has one");

    // The same substance, melted. Nothing else about the node changes.
    let cellulose = w.mixture_of(root).entries()[0].substance;
    let mut melted = Mixture::new();
    melted.add(cellulose, Phase::Liquid, 1.0);
    w.set_mixture(root, melted);
    w.tree.nodes[root.get()].epoch += 1;

    let s = w.surface_of_node(root);
    println!("  melted, it presents {} piece(s)", s.len());
    assert!(
        s.is_empty(),
        "a liquid's surface belongs to its container, and this one claimed {} pieces",
        s.len()
    );
}

/// A wooden box of six walls, stated by whatever built it — here the test,
/// standing in for a recipe.
fn boxed(seed: u64) -> (World, phys::ids::NodeIdx) {
    const L: f64 = 3.0;
    const R: f64 = 0.3;
    let (mut w, root) = wooden(seed, 4800.0, L * 3f64.sqrt(), 6);
    w.tree.refine(root);
    let members: Vec<(Vec3, Vec3)> = (0..3)
        .flat_map(|axis| {
            [1.0f64, -1.0].into_iter().map(move |sign| match axis {
                0 => (v3(sign * L, -L, 0.0), v3(sign * L, L, 0.0)),
                1 => (v3(-L, sign * L, 0.0), v3(L, sign * L, 0.0)),
                _ => (v3(-L, 0.0, sign * L), v3(L, 0.0, sign * L)),
            })
        })
        .collect();
    let nd = &mut w.tree.nodes[root.get()];
    let n = members.len();
    let (mut joints, mut base, mut tip) = (Vec::new(), Vec::new(), Vec::new());
    for (i, (b, t)) in members.iter().enumerate() {
        if let Some(body) = nd.bodies.get_mut(i) {
            body.pos = (*b + *t).scale(0.5);
            body.radius = R;
            body.mass = 4800.0 / n as f64;
        }
        joints.push(Joint {
            child: i as u32,
            parent: if i == 0 { NO_SUPPORT } else { 0 },
            at: *b,
            radius: R,
            integrity: 1.0,
        });
        base.push(*b);
        tip.push(*t);
    }
    nd.topology = Some(Topology {
        joints,
        support: (0..n).map(|i| if i == 0 { NO_SUPPORT } else { 0 }).collect(),
        site: (0..n as u32).collect(),
        base,
        tip,
        material: Material::green_wood(),
        ties: Vec::new(),
    });
    nd.epoch += 1;
    (w, root)
}

/// A box is six solid slabs and not one hull over the whole thing.
///
/// **The rule the whole decision rests on**, and the general statement of a
/// failure this engine already measured: one convex hull over a box *encloses
/// its own cavity*, so anything inside reads as deeply interpenetrating on
/// every frame. A shape language that can express "hollow" invites exactly that
/// mistake; one that cannot, cannot.
#[test]
fn a_hollow_box_is_its_walls_and_not_one_hull() {
    const L: f64 = 3.0;
    let (mut w, root) = boxed(0xB0FFED);
    let s = w.surface_of_node(root).clone();
    println!("  the box presents {} pieces", s.len());
    assert_eq!(s.len(), 6, "six walls, stated by what built it");

    // The cavity is where no primitive is, which is the whole of how a void is
    // represented. A point at the centre is outside every piece.
    let centre = Vec3::ZERO;
    for (i, p) in s.pieces().iter().enumerate() {
        let inside = (p.hull.centre() - centre).norm() < p.hull.bound() * 0.5;
        assert!(!inside, "piece {i} swallows the cavity it is supposed to enclose");
    }

    // And one hull over all six *would*, which is what makes this worth
    // asserting rather than assuming.
    let whole = phys::shape::Hull::of_spheres(
        s.pieces().iter().flat_map(|p| p.hull.spheres().iter().copied()),
    );
    let reach = whole.support(v3(0.0, 0.0, 1.0));
    println!(
        "  one hull over all six reaches {:.3} m up the axis and would enclose the cavity",
        reach.z
    );
    assert!(
        reach.z > L,
        "the single-hull comparison is not a comparison unless it really does cover the box"
    );
}

/// A surface is derived once and kept, and regenerated when the arrangement
/// changes.
///
/// D13 read literally — "deriving *once* and storing the result" — applied to
/// shape, which is the one place the engine never applied it. The old
/// `collision_shape` rebuilt a proxy from the member list on every frame and
/// persisted nothing, which is the half of the third axiom that saves nothing.
#[test]
fn a_surface_is_baked_once_and_kept_until_the_epoch_moves() {
    let (mut w, root) = wooden(0xBA4ED, 900.0, 2.0, 64);
    let before = w.stats.surfaces_baked;
    for _ in 0..50 {
        w.surface_of_node(root);
    }
    let after_first = w.stats.surfaces_baked;
    println!("  fifty asks cost {} bake(s)", after_first - before);
    assert_eq!(after_first - before, 1, "an undisturbed node bakes once");

    // Something happened to it.
    w.tree.nodes[root.get()].epoch += 1;
    w.surface_of_node(root);
    assert_eq!(
        w.stats.surfaces_baked - after_first,
        1,
        "a node whose arrangement changed must bake again"
    );
}

/// The surface turns with the node, and is stored in the node's own frame.
#[test]
fn a_surface_is_stored_in_the_nodes_own_frame() {
    // **A box and not a ball**, deliberately: a sphere turned about its own
    // centre is unchanged, so a rotation test on one asserts nothing. The box
    // has a wall whose ends a quarter turn visibly moves.
    let (mut w, root) = boxed(0x7024ED);
    let s = w.surface_of_node(root).clone();
    assert!(!s.is_empty());

    let quarter = Quat::from_axis_angle(v3(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
    let here = s.placed(Quat::IDENTITY, v3(10.0, 0.0, 0.0));
    let turned = s.placed(quarter, v3(10.0, 0.0, 0.0));
    assert_eq!(here.len(), s.len());
    // The *pieces* move, and the box as a whole is symmetric under a quarter
    // turn about z, so the extents agree while the individual walls do not.
    let reach = |hs: &[phys::shape::Hull], d: Vec3| {
        hs.iter().map(|h| h.support(d).dot(d)).fold(f64::MIN, f64::max)
    };
    let moved = here
        .iter()
        .zip(&turned)
        .map(|(a, b)| (a.centre() - b.centre()).norm())
        .fold(0.0f64, f64::max);
    println!(
        "  turned a quarter about z: extents {:.3} -> {:.3} along x, and the \
         furthest wall centre moved {moved:.3} m",
        reach(&here, v3(1.0, 0.0, 0.0)),
        reach(&turned, v3(1.0, 0.0, 0.0))
    );
    assert!(
        moved > 1.0,
        "a quarter turn should carry the walls with it, and moved them {moved} m"
    );
    assert!(
        (reach(&here, v3(0.0, 0.0, 1.0)) - reach(&turned, v3(0.0, 0.0, 1.0))).abs() < 1e-9,
        "a turn about z should not move the z extent"
    );
    // And the surface itself is untouched: it is stored in the node's own frame
    // and `placed` returns a view, because one with a position baked into it
    // would be wrong on the next frame.
    assert_eq!(
        w.surface_of_node(root).pieces()[0].hull,
        s.pieces()[0].hull,
        "placing a surface must not move the stored one"
    );
}

/// The surface materials reconcile with what the node is made of.
///
/// **D18's invariant**, and it is a conservation statement about material in
/// the same family as `summarise(sample(m)) == m`. Nothing in D17 or D18 alone
/// stops the two answers drifting: a node whose mixture is all water by mass
/// while its primitives present steel would be telling a contact something its
/// own bulk contradicts.
#[test]
fn a_surface_cannot_claim_what_the_node_is_not_made_of() {
    use phys::shape::{Mismatch, Piece, Surface};

    let (mut w, root) = wooden(0x5A1E, 900.0, 2.0, 16);
    let mixture = w.mixture_of(root);
    let wood = mixture.entries()[0].substance;

    // What the engine bakes reconciles.
    let baked = w.surface_of_node(root).clone();
    assert!(
        baked.reconcile(&mixture, 900.0, 1e-6).is_ok(),
        "the engine's own bake should agree with the node's pools"
    );
    assert_eq!(w.stats.surface_mismatches, 0);

    // A surface claiming a substance the node has no solid pool of is refused.
    let steel = w
        .substances
        .intern(phys::material::substances::iron_arrangement())
        .expect("iron analyses");
    let lying = Surface::new(vec![Piece {
        hull: phys::shape::Hull::sphere(Vec3::ZERO, 1.0),
        material: Material::steel(),
        substance: steel,
        mass: 900.0,
    }]);
    assert_eq!(
        lying.reconcile(&mixture, 900.0, 1e-6),
        Err(Mismatch::NotHeld { substance: steel }),
        "a wooden node presenting steel is two answers to one question"
    );

    // And one claiming more of a substance than the node holds is refused too.
    let greedy = Surface::new(vec![Piece {
        hull: phys::shape::Hull::sphere(Vec3::ZERO, 1.0),
        material: Material::green_wood(),
        substance: wood,
        mass: 1800.0,
    }]);
    assert!(
        matches!(
            greedy.reconcile(&mixture, 900.0, 1e-6),
            Err(Mismatch::OverClaimed { .. })
        ),
        "a 900 kg node cannot present 1800 kg of wood"
    );
}
