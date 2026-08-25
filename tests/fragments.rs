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
use phys::prolong::prolong_structured;
use phys::solvers::structure::*;
use phys::state::Aggregate;
use phys::topology::{Material, Member, Topology};
use phys::units::{Tier, YEAR};

fn tree(mass: f64, budget: usize) -> (Vec<phys::state::Body>, Topology) {
    let mut m = phys::morph::Morphology::new(Program::Tree, 0xACE, 0x1234, 0);
    m.built = mass;
    m.age = 45.0 * YEAR;
    let mut agg = Aggregate::neutral(mass, m.extent(), 291.0, Program::Tree.substrate());
    agg.chemical_energy = m.stored_energy();
    let (b, t, _) = prolong_structured(&agg, &m, budget, 7, 0x1234, 0);
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
    let mut topo = Topology::from_parts(&members, &[], Material::DRY_TIMBER);
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
    field.apply(&weather::gravity(), &piece_bodies, &piece);
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
    field.apply(&weather::wind(48.0, v3(1.0, 0.0, 0.0)), &bodies, &topo);
    field.apply(&weather::gravity(), &bodies, &topo);
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
    let mut world = World::new(galaxy(0x5EED, 1e9), 20.0);
    let root = world.tree.root;
    world.tree.refine(root);
    let node = world.tree.promote(root, 3, default_spec(Tier::Stellar));
    {
        let n = &mut world.tree.nodes[node.get()];
        n.agg = Aggregate::neutral(4000.0, 6.0, 291.0, Program::Tree.substrate());
        n.spec.count = 900;
    }
    world.plant(node, Program::Tree, Environment::default());
    for _ in 0..70 {
        world.grow_node(node, YEAR);
    }

    // A gale hard enough to take limbs off, and no harder. Past about 44 m/s
    // this tree is shredded rather than pruned, and debris from a shredded
    // crown mostly falls through what is left of it — which is true, and not
    // what this test is about.
    let out = world.damage(node, &[weather::wind(38.0, v3(1.0, 0.0, 0.0))]);
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
    field.apply(&weather::wind(55.0, v3(1.0, 0.0, 0.0)), &bodies, &topo);
    field.apply(&weather::gravity(), &bodies, &topo);
    let loads = analyse(&bodies, &topo, &field);
    let failures = apply_failures(&bodies, &mut topo, &loads, &field);
    let cut = detach(&topo, &failures.broken_members);
    let Some(members) = cut.pieces.iter().max_by_key(|p| p.len()) else {
        panic!("nothing broke");
    };
    let mut frag = Fragment::new(&bodies, &topo, members).expect("a piece");

    use phys::solvers::frame::Dof;
    let n = frag.dynamics.dynamics.frame.nodes.len();
    // Where the ground is in this structure's own frame: a generated tree is
    // recentred on its centre of mass, so its foundations are not at zero.
    let ground = (0..topo.support.len())
        .filter(|&i| topo.support[i] == NO_SUPPORT && topo.bonds[i].radius > 0.0)
        .map(|i| topo.base[i].z.min(topo.tip[i].z))
        .fold(f64::INFINITY, f64::min);
    let start = frag.lowest() - ground;
    let mut steps = 0;
    for _ in 0..600 {
        let mut load = vec![Dof::default(); n];
        for i in 0..n {
            let m = frag.dynamics.dynamics.frame.lumped[i].t.z;
            load[i].t = G_EARTH.scale(m);
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
