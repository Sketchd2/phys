//! Breaking apart: re-rooting, falling, and what the falling hits.
//!
//! A break is not a member disappearing. It is a member ceasing to be
//! *supported*, and everything hanging off it comes with it — as its own
//! object, with its own roots, its own centre of mass, and a reason to be
//! analysed all over again.

use phys::engine::{galaxy, World};
use phys::math::v3;
use phys::morph::{Environment, Program, NO_SUPPORT};
use phys::engine::default_spec;
use phys::sampler::sample_structured;
use phys::solvers::structure::*;
use phys::state::{Composition, Matter};
use phys::topology::{Material, Member, Topology};
use phys::units::{Tier, YEAR};

/// A node on the surface of an Earth, ready to have something planted in it.
///
/// `docs/PLAY.md` D6 made gravity derived, and that changed what these tests
/// have to set up. A node promoted straight out of a galaxy root feels the
/// galaxy's field — about 10^-13 m/s^2 — so nothing falls, no snow load breaks
/// anything, and a gale takes limbs off that then hang in the air. That is
/// correct: there is nothing underneath it. It used to feel Earth's surface
/// gravity wherever it was, because the load came from
/// `solvers::structure::G_EARTH` rather than from anything the scene contained.
///
/// So the scene contains a planet. `Tree::gravity_at` derives 9.82 m/s^2 from
/// its mass and radius, with no mention of Earth anywhere.
///
/// Placed on the **+z axis** deliberately. The derived field points at the
/// planet's centre, and a structure's own geometry is generated with `+z` up,
/// so putting the node anywhere else would load a tree sideways — `Motion`
/// carries an orientation and `gravity_at` does not yet compose it. See the
/// note on `Tree::gravity_at`.
fn on_an_earth(seed: u64, mass: f64, radius: f64, count: usize) -> (World, phys::ids::NodeIdx) {
    const EARTH_MASS: f64 = 5.972e24;
    const EARTH_RADIUS: f64 = 6.371e6;
    let mut w = World::new(galaxy(seed, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let planet = w.tree.promote(root, 7, phys::engine::default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.spec.count = 8;
    }
    w.tree.refine(planet);
    let node = w.tree.promote(planet, 0, phys::engine::default_spec(Tier::Continuum));
    {
        let n = &mut w.tree.nodes[node.get()];
        n.matter = Matter::neutral(mass, radius, 291.0, Program::Tree.substrate());
        n.spec.count = count;
        n.motion.offset = v3(0.0, 0.0, EARTH_RADIUS);
    }
    (w, node)
}


/// Earth's surface gravity, stated rather than assumed.
///
/// `docs/PLAY.md` D6 retired `solvers::structure::G_EARTH`: a structure is now
/// proportioned and loaded in the field it is actually in, derived by
/// `Tree::gravity_at` from whatever it sits inside. These tests are about
/// structure generation rather than about gravity, so they name a field and
/// hold it fixed — which is what the constant was doing for them, said out loud.
const SURFACE_G: phys::math::Vec3 = phys::math::Vec3 { x: 0.0, y: 0.0, z: -9.80665 };


fn tree(mass: f64, budget: usize) -> (Vec<phys::state::Body>, Topology) {
    let mut m = phys::morph::Morphology::new(Program::Tree, 0xACE, 0x1234, 0);
    m.built = mass;
    m.age = 45.0 * YEAR;
    let mut matter = Matter::neutral(mass, m.extent(), 291.0, Program::Tree.substrate());
    matter.chemical_energy = m.stored_energy();
    let (b, t, _) = sample_structured(&matter, &m, budget, 7, 0x1234, 0, SURFACE_G);
    (b, t)
}

/// Cutting a structure has to produce a piece whose support graph makes sense
/// on its own.
///
/// This is the part that is easy to get wrong and impossible to see: a branch
/// that keeps its old support index is still, as far as any analysis is
/// concerned, being held up by the trunk it fell off.
#[test]
fn a_severed_piece_is_re_rooted() {
    // A chain of five members, anchored at the first. Cut the middle one.
    let members: Vec<Member> = (0..5)
        .map(|i| {
            Member::new(
                v3(0.0, 0.0, i as f64),
                v3(0.0, 0.0, i as f64 + 1.0),
                0.05,
                if i == 0 { NO_SUPPORT } else { i as u32 - 1 },
            )
        })
        .collect();
    let mut topo = Topology::from_parts(&members, &[], Material::dry_timber());
    let bodies: Vec<phys::state::Body> = (0..5)
        .map(|i| phys::state::Body {
            pos: v3(0.0, 0.0, i as f64 + 0.5),
            mass: 10.0,
            radius: 0.05,
            ..Default::default()
        })
        .collect();

    // Member 2 loses its support, exactly as `apply_failures` would leave it.
    topo.support[2] = NO_SUPPORT;
    let cut = detach(&topo, &[2]);
    println!("  standing {:?}, pieces {:?}", cut.standing, cut.pieces);
    assert_eq!(cut.standing, vec![0, 1], "the base should still be standing");
    assert_eq!(cut.pieces.len(), 1, "one break, one piece");
    assert_eq!(cut.pieces[0], vec![2, 3, 4], "everything above the break comes with it");

    let (piece_bodies, piece) = extract(&bodies, &topo, &cut.pieces[0]).expect("a piece");
    assert_eq!(piece_bodies.len(), 3);
    // Re-rooted: the break is now an anchor, and the two above it hang off it
    // in the piece's own index space.
    assert_eq!(piece.support, vec![NO_SUPPORT, 0, 1], "support was not remapped");
    assert!(piece.is_determinate(), "a cut branch is still a tree");
    // And it can be analysed on its own terms, which is the point.
    let mut field = LoadField::new(piece_bodies.len(), 290.0);
    field.apply(&weather::gravity(SURFACE_G), &piece_bodies, &piece);
    let loads = analyse(&piece_bodies, &piece, &field);
    println!(
        "  the piece analysed alone: root carries {:.1} kg, tip carries {:.1} kg",
        loads[0].carried, loads[2].carried
    );
    assert!(
        loads[0].carried > loads[2].carried,
        "the piece's own root should carry the most"
    );
    assert!(
        (loads[0].carried - 30.0).abs() < 1e-6,
        "the piece's root should carry all three members, not {:.3} kg",
        loads[0].carried
    );
}

/// A break in a real tree must produce a real piece.
#[test]
fn breaking_a_tree_produces_falling_pieces() {
    let (bodies, mut topo) = tree(900.0, 1500);
    let mut field = LoadField::new(bodies.len(), 290.0);
    field.apply(&weather::wind(85.0, v3(1.0, 0.0, 0.0)), &bodies, &topo);
    field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
    let loads = analyse(&bodies, &topo, &field);
    let failures = apply_failures(&bodies, &mut topo, &loads, &field);
    assert!(
        !failures.broken_members.is_empty(),
        "a 48 m/s wind should break something"
    );

    let cut = detach(&topo, &failures.broken_members);
    let mut heaviest = 0.0f64;
    let mut total = 0.0f64;
    let mut pieces = 0;
    for members in &cut.pieces {
        let Some(frag) = Fragment::new(&bodies, &topo, members) else {
            continue;
        };
        pieces += 1;
        heaviest = heaviest.max(frag.mass());
        total += frag.mass();
    }
    println!(
        "  {} joints broke into {pieces} pieces totalling {total:.1} kg; heaviest {heaviest:.1} kg; \
         {} members still standing",
        failures.broken_members.len(),
        cut.standing.len()
    );
    assert!(pieces > 0, "nothing came away");
    assert!(total > 0.0, "the pieces weigh nothing");
    // Every member is either standing or in exactly one piece.
    let mut seen = vec![0u32; topo.support.len()];
    for &m in &cut.standing {
        seen[m as usize] += 1;
    }
    for members in &cut.pieces {
        for &m in members {
            seen[m as usize] += 1;
        }
    }
    assert!(
        seen.iter().all(|&c| c == 1),
        "a member was counted twice or not at all"
    );
}

/// The whole chain: a limb comes off, falls, hits what is under it, and what it
/// hits has to answer for it.
#[test]
fn a_falling_limb_damages_what_it_lands_on() {
    let (mut world, node) = on_an_earth(0x5EED, 4000.0, 6.0, 900);
    world.plant(node, Program::Tree, Some(Environment::default()));
    for _ in 0..70 {
        world.grow_node(node, YEAR);
    }

    // A gale hard enough to take limbs off, and no harder.
    //
    // **It was 38 m/s and is now 42.** `docs/PLAY.md` D14 made the material
    // measured rather than tabulated, and green wood comes out 3.2x stiffer
    // than the table held — a factor that is the accuracy the derivation claims
    // and not a change in behaviour. This tree fails by *buckling*, which
    // scales with stiffness, so the wind that prunes it rises with its square
    // root.
    //
    // The window is narrow and worth stating: measured, 40 m/s breaks 166
    // joints and the debris hits 19 members, 42 sheds cleanly, and 44 breaks
    // 281 and shreds the crown so that what comes off falls through what is
    // left of it. That last is true and is not what this test is about.
    let out = world.damage(node, &[weather::wind(42.0, v3(1.0, 0.0, 0.0))]);
    println!(
        "  the gale broke {} joints into {} falling pieces",
        out.broken_joints, out.detached_pieces
    );
    assert!(out.broken_joints > 0, "the gale broke nothing");
    assert!(out.detached_pieces > 0, "nothing came away as its own object");
    assert!(!world.falling().is_empty(), "nothing is falling");

    let start: f64 = world.falling().iter().map(|(_, f)| f.mass()).sum();
    let start_height = world
        .falling()
        .iter()
        .map(|(_, f)| f.lowest())
        .fold(f64::NEG_INFINITY, f64::max);

    let mut contacts = 0;
    let mut struck = 0;
    let mut secondary = 0;
    let mut secondary_mass = 0.0;
    let mut settled = 0;
    let mut peak = 0.0f64;
    for _ in 0..400 {
        let r = world.drop_fragments(0.02);
        contacts += r.contacts;
        struck += r.struck_members;
        secondary += r.secondary_breaks;
        secondary_mass += r.secondary_mass;
        settled += r.settled;
        peak = peak.max(r.peak_utilisation);
        if world.falling().is_empty() {
            break;
        }
    }
    let end_height = world
        .falling()
        .iter()
        .map(|(_, f)| f.lowest())
        .fold(f64::INFINITY, f64::min);

    println!(
        "  {start:.1} kg of debris fell from up to {start_height:.2} m; {contacts} contacts, \
         {struck} members struck, {secondary} joints broken by the impacts \
         ({secondary_mass:.1} kg), peak utilisation under impact {peak:.2}, {settled} pieces \
         came to rest; lowest still airborne {}",
        if end_height.is_finite() { format!("{end_height:.2} m") } else { "none".into() }
    );
    assert!(start > 0.0, "nothing came away");
    assert!(start_height > 1.0, "the debris started at ground level");
    assert!(struck > 0, "the debris fell through everything");
    // The point of the exercise: a limb on its way down is a load on whatever
    // is beneath it, and the ordinary stress calculation decides the rest.
    assert!(
        peak > 1.0,
        "no impact loaded a member past failure; the worst reached {peak:.2}"
    );
    assert!(
        secondary > 0,
        "the impacts broke nothing at all, which is not a fall, it is a fade"
    );
    assert_eq!(settled, 24, "the debris never came to rest");
}

/// Debris has to stop. A simulation that never lets go of a fallen branch
/// spends its whole budget on litter.
#[test]
fn debris_comes_to_rest() {
    let (bodies, mut topo) = tree(900.0, 800);
    let mut field = LoadField::new(bodies.len(), 290.0);
    field.apply(&weather::wind(95.0, v3(1.0, 0.0, 0.0)), &bodies, &topo);
    field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
    let loads = analyse(&bodies, &topo, &field);
    let failures = apply_failures(&bodies, &mut topo, &loads, &field);
    let cut = detach(&topo, &failures.broken_members);
    let Some(members) = cut.pieces.iter().max_by_key(|p| p.len()) else {
        panic!("nothing broke");
    };
    let mut frag = Fragment::new(&bodies, &topo, members).expect("a piece");

    use phys::solvers::frame::Dof;
    let n = frag.dynamics.dynamics.frame.joints.len();
    // Where the ground is in this structure's own frame: a generated tree is
    // recentred on its centre of mass, so its foundations are not at zero.
    let ground = (0..topo.support.len())
        .filter(|&i| topo.support[i] == NO_SUPPORT && topo.joints[i].radius > 0.0)
        .map(|i| topo.base[i].z.min(topo.tip[i].z))
        .fold(f64::INFINITY, f64::min);
    let start = frag.lowest() - ground;
    let mut steps = 0;
    for _ in 0..600 {
        let mut load = vec![Dof::default(); n];
        for i in 0..n {
            let m = frag.dynamics.dynamics.frame.lumped[i].t.z;
            load[i].t = SURFACE_G.scale(m);
        }
        frag.dynamics.dynamics.step(&load, 0.01);
        steps += 1;
        let hits = frag.contacts(&[], &Topology::default(), ground);
        frag.resolve(&hits, 0.15);
        if frag.at_rest() {
            break;
        }
    }
    println!(
        "  a {:.1} kg piece fell from {start:.2} m and landed after {steps} steps ({:.2} s)",
        frag.mass(),
        steps as f64 * 0.01
    );
    assert!(frag.at_rest(), "the piece never came to rest");
    // Free fall from that height, as a sanity check that it fell rather than
    // being teleported: t = sqrt(2h/g).
    let expected = (2.0 * start.max(0.0) / 9.80665).sqrt();
    let took = steps as f64 * 0.01;
    assert!(
        took > expected * 0.5 && took < expected * 3.0 + 1.0,
        "fell for {took:.2} s; free fall from {start:.2} m is {expected:.2} s"
    );
}

/// A branch lands on the next tree.
///
/// Phase 1's last done-when, and the one that needed the adjacency relation
/// rather than another mechanism. Until this, a falling piece could strike only
/// the structure it came off: `drop_fragments` keyed its candidates by the
/// piece's own node and tested nothing else, so a limb passed through a
/// neighbouring tree as though it were not there.
///
/// Measured before the fix, on exactly this arrangement — two trees 8 m apart,
/// each grown to a 10.4 m radius, so the crowns overlap by two metres: nine
/// limbs came off one in a gale, made 225 contacts and struck 160 members,
/// every one of them in the tree they fell from.
#[test]
fn a_branch_lands_on_the_next_tree() {
    let (mut w, a) = on_an_earth(0x5EED, 4000.0, 6.0, 900);
    w.plant(a, Program::Tree, Some(Environment::default()));
    for _ in 0..70 {
        w.grow_node(a, YEAR);
    }

    // A second tree beside the first, close enough that their crowns overlap.
    let planet = w.tree.nodes[a.get()].parent;
    let b = w.tree.promote(planet, 1, default_spec(Tier::Continuum));
    assert!(!b.is_none(), "the planet should have another body to promote");
    let beside = w.tree.nodes[a.get()].motion.offset + v3(8.0, 0.0, 0.0);
    {
        let n = &mut w.tree.nodes[b.get()];
        n.matter = Matter::neutral(4000.0, 6.0, 291.0, Program::Tree.substrate());
        n.spec.count = 900;
        n.motion.offset = beside;
    }
    w.plant(b, Program::Tree, Some(Environment::default()));
    for _ in 0..70 {
        w.grow_node(b, YEAR);
    }
    w.tree.refine(b);

    let (ra, rb) = (
        w.tree.nodes[a.get()].matter.radius,
        w.tree.nodes[b.get()].matter.radius,
    );
    assert!(
        ra + rb > 8.0,
        "the crowns have to overlap for this to be a test: {ra} and {rb} at 8 m"
    );

    // And a third, far enough away that nothing should reach it. This is the
    // control, and it is what makes the *offset* part of the test: a neighbour
    // whose geometry is read without being moved into the falling piece's frame
    // sits on top of the piece instead of beside it, so it gets hit harder
    // rather than not at all. Asserting only that something was struck passed
    // with the offset deleted.
    let far = w.tree.promote(planet, 2, default_spec(Tier::Continuum));
    assert!(!far.is_none());
    let away = w.tree.nodes[a.get()].motion.offset + v3(60.0, 0.0, 0.0);
    {
        let n = &mut w.tree.nodes[far.get()];
        n.matter = Matter::neutral(4000.0, 6.0, 291.0, Program::Tree.substrate());
        n.spec.count = 900;
        n.motion.offset = away;
    }
    w.plant(far, Program::Tree, Some(Environment::default()));
    for _ in 0..70 {
        w.grow_node(far, YEAR);
    }
    w.tree.refine(far);
    assert!(
        w.tree.nodes[far.get()].matter.radius * 2.0 < 60.0,
        "the far tree has to be out of reach for this to be a control"
    );

    // 38 m/s before the material was measured; see
    // `a_falling_limb_damages_what_it_lands_on` for the arithmetic.
    let out = w.damage(a, &[weather::wind(70.0, v3(1.0, 0.0, 0.15))]);
    assert!(out.detached_pieces > 0, "the gale took nothing off the first tree");
    let before = w.tree.nodes[b.get()].morphology.as_ref().unwrap().built;
    let before_far = w.tree.nodes[far.get()].morphology.as_ref().unwrap().built;

    for _ in 0..200 {
        w.drop_fragments(0.05);
        if w.falling().is_empty() {
            break;
        }
    }

    let after = w.tree.nodes[b.get()].morphology.as_ref().unwrap().built;
    assert!(
        after < before,
        "{} limbs fell through a tree whose crown they were inside: the second \
         tree still has all {before} kg of it",
        out.detached_pieces
    );
    let after_far = w.tree.nodes[far.get()].morphology.as_ref().unwrap().built;
    assert_eq!(
        after_far, before_far,
        "a tree sixty metres away lost {} kg to limbs falling off another one",
        before_far - after_far
    );
}
