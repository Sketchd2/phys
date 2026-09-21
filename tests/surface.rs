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
            bond: phys::chem::SubstanceId::UNSPECIATED,
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
        bonds: Vec::new(),
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

/// Something set down on a coursed wall rests on it, rather than in it.
///
/// `docs/PLAY.md` Phase 4's third done-when. A grown or coursed structure
/// emitted one capsule per member, so a wall was a row of beads: the top of a
/// block is its own round shoulder, the seam between two blocks is a valley
/// between two shoulders, and a thing set down on the wall finds a different
/// height depending where along it you put it. Measured before this, on the
/// same wall: **0.169 m** of scallop on a 1.2 m panel.
///
/// What is measured here is the height a small sphere comes to rest at, swept
/// along the top of the wall. A flat top gives one height everywhere; a row of
/// beads gives a sawtooth, and the tooth is what something walking along the
/// wall falls into.
#[test]
fn something_set_down_on_a_coursed_wall_rests_on_top_of_it() {
    use phys::morph::Program;
    let (mut w, root) = wooden(0xC0_0453D, 4000.0, 1.2, 400);
    w.emplace(root, Program::Wall, 4000.0, None);
    w.tree.refine(root);
    let surface = w.surface_of_node(root).clone();
    assert!(surface.len() > 8, "a coursed wall is many blocks: {}", surface.len());

    let hulls: Vec<&phys::shape::Hull> = surface.hulls().collect();
    // The wall's own bounding box, measured off the pieces.
    let (mut lo_x, mut hi_x) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut lo_y, mut hi_y, mut hi_z) = (f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for h in &hulls {
        for sp in h.spheres() {
            lo_x = lo_x.min(sp.centre.x - sp.radius);
            hi_x = hi_x.max(sp.centre.x + sp.radius);
            lo_y = lo_y.min(sp.centre.y - sp.radius);
            hi_y = hi_y.max(sp.centre.y + sp.radius);
            hi_z = hi_z.max(sp.centre.z + sp.radius);
        }
    }
    let width = hi_x - lo_x;
    let thickness = hi_y - lo_y;
    let probe = width * 0.01;

    // Where a small sphere comes to rest, directly above each sample point.
    // The search runs between a height that is inside the wall — its own centre
    // plane — and one clear above it.
    let rest_at = |x: f64, y: f64| -> Option<f64> {
        let touches = |z: f64| {
            hulls.iter().any(|h| {
                phys::shape::closest(&phys::shape::Hull::sphere(v3(x, y, z), probe), h)
                    .map(|c| c.gap <= 0.0)
                    .unwrap_or(false)
            })
        };
        let (mut lo, mut hi) = (0.0, hi_z + 4.0 * probe);
        if !touches(lo) {
            return None;
        }
        for _ in 0..48 {
            let mid = 0.5 * (lo + hi);
            if touches(mid) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(lo)
    };

    // A grid over the middle of the top face, not a line along it. A line
    // along the courses is exactly the sweep a row of capsules passes: the
    // blocks of one course are collinear, so their tubes union into a smooth
    // cylinder and the top *along* the wall is flat. The fall is off the
    // shoulder — across the thickness — and between one course and the next.
    let mid_x = 0.5 * (lo_x + hi_x);
    let mid_y = 0.5 * (lo_y + hi_y);
    let mut heights = Vec::new();
    for i in 0..8 {
        for j in 0..8 {
            let x = mid_x + width * 0.3 * (2.0 * (i as f64 + 0.5) / 8.0 - 1.0);
            let y = mid_y + thickness * 0.3 * (2.0 * (j as f64 + 0.5) / 8.0 - 1.0);
            if let Some(z) = rest_at(x, y) {
                heights.push(z);
            }
        }
    }
    assert!(heights.len() > 32, "the sweep found the wall: {} of 64", heights.len());
    let hi = heights.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let lo = heights.iter().cloned().fold(f64::INFINITY, f64::min);
    println!(
        "  a wall {:.3} m across and {:.3} m thick: rest height varies {:.4} m \
         over {} samples",
        width,
        thickness,
        hi - lo,
        heights.len()
    );
    assert!(
        hi - lo < 0.01 * thickness,
        "something set down on the top of the wall drops {:.4} m, against a \
         wall {:.3} m thick",
        hi - lo,
        thickness
    );
}

/// A patch of ground is a surface, not a bed of posts.
///
/// The same measurement as the wall, on the thing `docs/PLAY.md` Phase 4 says
/// it first matters for: "Terrain is where it is first load-bearing rather than
/// cosmetic, because something has to *stand* on it."
///
/// Two defects met here, and the second is upstream of the first. The columns
/// were round on a square grid, so even touching at their midlines they left a
/// hole at every corner of it. And `render_terrain` drew a *cube*: the columns
/// stood on `z = -1`, one full half-width down, so a patch that says it is
/// 19.4 m across and 2.4 m deep drew one 25.1 m across and 22.3 m deep — nine
/// times its own volume at its own density — and the sampler's cross-section
/// correction thinned every column to **0.2848** of its cell to make the
/// numbers agree. The ground covered less than a third of the ground.
#[test]
fn a_patch_of_ground_has_no_holes_in_it() {
    use phys::morph::Program;
    let (mut w, root) = wooden(0x6204_11D, 2.4e6, 12.0, 400);
    w.emplace(root, Program::Terrain, 2.4e6, None);
    w.tree.refine(root);

    // Every column fills its cell exactly: the correction that used to thin
    // them has nothing left to do, because the patch now draws the slab it says
    // it is and the cell is stated rather than free.
    let bodies = &w.tree.nodes[root.get()].bodies;
    let boxed: Vec<&phys::state::Body> = bodies.iter().filter(|b| b.is_boxed()).collect();
    assert!(boxed.len() > 16, "a patch is a grid of columns: {}", boxed.len());
    let mut xs: Vec<f64> = boxed.iter().map(|b| b.pos.x).collect();
    xs.sort_by(|a, b| a.total_cmp(b));
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    let pitch = xs[1] - xs[0];
    let covering = 2.0 * boxed[0].half.x / pitch;
    println!("  cells {:.4} m apart, columns {:.4} m wide, covering {covering:.4}", pitch, 2.0 * boxed[0].half.x);
    assert!(
        (covering - 1.0).abs() < 1e-9,
        "the columns tile their cells, and cover {covering}"
    );

    let surface = w.surface_of_node(root).clone();
    let hulls: Vec<&phys::shape::Hull> = surface.hulls().collect();
    let (mut hi_z, mut lo_z, mut wide) = (f64::NEG_INFINITY, f64::INFINITY, 0.0f64);
    for h in &hulls {
        for sp in h.spheres() {
            hi_z = hi_z.max(sp.centre.z + sp.radius);
            lo_z = lo_z.min(sp.centre.z - sp.radius);
            wide = wide.max(sp.centre.x.abs() + sp.radius);
        }
    }
    let depth = hi_z - lo_z;
    let probe = pitch * 0.05;
    let ground_at = |x: f64, y: f64| -> Option<f64> {
        let touches = |z: f64| {
            hulls.iter().any(|h| {
                phys::shape::closest(&phys::shape::Hull::sphere(v3(x, y, z), probe), h)
                    .map(|c| c.gap <= 0.0)
                    .unwrap_or(false)
            })
        };
        let (mut lo, mut hi) = (lo_z + 0.05 * depth, hi_z + 4.0 * probe);
        if !touches(lo) {
            return None;
        }
        for _ in 0..48 {
            let mid = 0.5 * (lo + hi);
            if touches(mid) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(lo)
    };

    // Walk a line across the patch, through the corners of the grid rather than
    // down the middle of a row — the holes a round column leaves are exactly
    // where four cells meet.
    let n = 40;
    let span = wide * 0.5;
    let mut found = 0;
    let mut worst_step = 0.0f64;
    let mut last: Option<f64> = None;
    for k in 0..n {
        let t = -span + 2.0 * span * (k as f64 + 0.5) / n as f64;
        match ground_at(t, t) {
            Some(z) => {
                found += 1;
                if let Some(p) = last {
                    worst_step = worst_step.max((z - p).abs());
                }
                last = Some(z);
            }
            None => last = None,
        }
    }
    println!(
        "  a diagonal of {n} steps found ground {found} times, worst step \
         {worst_step:.4} m over a patch {depth:.3} m deep"
    );
    assert_eq!(found, n, "every step of the walk is on ground");
    assert!(
        worst_step < 0.25 * depth,
        "a step across the patch drops {worst_step:.4} m, against a patch \
         {depth:.3} m deep"
    );
}
