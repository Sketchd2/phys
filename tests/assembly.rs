//! A composite is one node with a recipe. `docs/PLAY.md` D15.

use phys::assembly::{Assembly, Part, UNJOINED};
use phys::chem::{Arrangement, Bond, Element, Order, Phase, SubstanceId};
use phys::engine::{galaxy, World};
use phys::math::{v3, Vec3};
use phys::morph::Program;
use phys::state::Matter;
use phys::units::*;

/// A wooden box: six panels, five of them joined to the floor.
///
/// 1.2 m on a side, 25 mm planks. The floor is the anchor and everything else
/// is joined to it, which is a star and not a chain — a crate is not a tower,
/// and nothing about a wall depends on the wall next to it.
fn box_parts(oak: SubstanceId, glue: SubstanceId, density: f64) -> Assembly {
    let half = 0.6;
    let t = 0.0125;
    let faces: [(Vec3, Vec3); 6] = [
        (v3(0.0, 0.0, -half), v3(half, half, t)),   // floor
        (v3(0.0, 0.0, half), v3(half, half, t)),    // lid
        (v3(-half, 0.0, 0.0), v3(t, half, half)),   // four sides
        (v3(half, 0.0, 0.0), v3(t, half, half)),
        (v3(0.0, -half, 0.0), v3(half, t, half)),
        (v3(0.0, half, 0.0), v3(half, t, half)),
    ];
    let parts = faces
        .iter()
        .enumerate()
        .map(|(i, (at, h))| {
            let mass = 8.0 * h.x * h.y * h.z * density;
            let p = Part::new(*at, *h, mass, oak, i as u32);
            if i == 0 {
                p
            } else {
                // Glued along the edge it meets the floor on: the seam is the
                // panel's thickness by its width, not its whole face.
                p.joined(0, glue, 2.0 * half * 2.0 * t)
            }
        })
        .collect();
    Assembly::new(parts)
}

/// Cellulose and a glue, built by hand: nothing in the engine knows what wood
/// is, which is the point of the chemistry layer.
fn woods(w: &mut World) -> (SubstanceId, SubstanceId) {
    let oak = w
        .substances
        .intern(Arrangement::molecule(
            vec![Element(6), Element(1), Element(1), Element(8)],
            vec![Bond::new(0, 1, Order::Single), Bond::new(0, 2, Order::Single), Bond::new(0, 3, Order::Single)],
        ))
        .expect("cellulose analyses");
    w.substances.name(oak, "cellulose");
    let glue = w
        .substances
        .intern(Arrangement::molecule(
            vec![Element(6), Element(1), Element(8)],
            vec![Bond::new(0, 1, Order::Single), Bond::new(0, 2, Order::Double)],
        ))
        .expect("glue analyses");
    w.substances.name(glue, "glue");
    (oak, glue)
}

fn a_box() -> (World, phys::ids::NodeIdx, SubstanceId, SubstanceId) {
    let mut w = World::new(galaxy(0xB0C5, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let node = w.tree.promote(root, 3, phys::engine::default_spec(Tier::Continuum));
    let (oak, glue) = woods(&mut w);
    let parts = box_parts(oak, glue, 700.0);
    let mass = parts.mass();
    {
        let n = &mut w.tree.nodes[node.get()];
        n.matter = Matter::neutral(mass, 1.0, 291.0, Program::Tree.substrate());
        let mut mix = phys::chem::Mixture::new();
        mix.add(oak, Phase::Solid, 1.0);
        n.matter.mixture = mix;
        n.spec.count = 64;
    }
    w.assemble(node, Program::Tree, parts, None);
    (w, node, oak, glue)
}

#[test]
fn a_box_is_one_node_with_six_walls() {
    let (mut w, node, _, _) = a_box();
    let n = &w.tree.nodes[node.get()];
    let m = n.morphology.as_ref().expect("the box has a recipe");
    let a = m.assembly.as_ref().expect("and the recipe is an assembly");
    println!("  recipe: {} parts, {} bytes, extent {:.3} m, {:.1} kg",
        a.len(), m.state_bytes(), m.extent(), a.mass());
    assert_eq!(a.len(), 6);
    assert_eq!(a.parts.iter().filter(|p| p.joined_to != UNJOINED).count(), 5);

    // Sampled, the parts land where the recipe put them.
    let bodies = w.tree.refine(node).to_vec();
    let structural = bodies.iter().filter(|b| b.radius > 0.0).count();
    println!("  sampled: {} bodies, {} structural", bodies.len(), structural);
    let n = &w.tree.nodes[node.get()];
    let topo = n.topology.as_ref().expect("an assembly has a topology");
    for i in 0..6 {
        let want = n.morphology.as_ref().unwrap().assembly.as_ref().unwrap().parts[i].at;
        let got = bodies[i].pos;
        println!(
            "    part {i}: recipe {:?} sampled ({:.4}, {:.4}, {:.4}) joint r {:.4}",
            (want.x, want.y, want.z), got.x, got.y, got.z, topo.joints[i].radius
        );
    }
}

/// The parts land exactly where the recipe put them.
///
/// Sharper than it looks. `matter.radius` means *equivalent uniform sphere*
/// everywhere in this engine — `summarise` reports it from an rms and `sample`
/// scales what it draws until it comes back — and a recipe that stated its
/// bounding radius instead put this box's walls 1.87x too far out, silently,
/// because the sampler helpfully scaled the whole thing up until the two
/// numbers agreed. A recipe that says where a wall is and a sampler that puts
/// it somewhere else is the failure this whole phase is about, so it is
/// asserted rather than eyeballed.
#[test]
fn the_parts_land_where_the_recipe_put_them() {
    let (mut w, node, _, _) = a_box();
    let recipe: Vec<Vec3> = w.tree.nodes[node.get()]
        .morphology
        .as_ref()
        .unwrap()
        .assembly
        .as_ref()
        .unwrap()
        .parts
        .iter()
        .map(|p| p.at)
        .collect();
    let bodies = w.tree.refine(node).to_vec();
    let mut worst = 0.0f64;
    for (want, got) in recipe.iter().zip(&bodies) {
        worst = worst.max((*want - got.pos).norm());
    }
    println!("  worst part displacement: {worst:.3e} m over {} parts", recipe.len());
    assert!(worst < 1e-9, "the sampler moved the parts by {worst:.3e} m");

    // And the seams are the seams somebody glued, not the seams a
    // fully-stressed optimiser would have wanted.
    let topo = w.tree.nodes[node.get()].topology.as_ref().unwrap();
    let seam = (0.03f64 / std::f64::consts::PI).sqrt();
    for i in 1..6 {
        assert!(
            (topo.joints[i].radius - seam).abs() / seam < 1e-9,
            "joint {i} is {:.5} m, not the {:.5} m seam the recipe states",
            topo.joints[i].radius,
            seam
        );
    }
}

/// Six panels present six surfaces, each with its own material.
///
/// D18 says material attaches per primitive, and until an assembly existed
/// nothing in the engine had more than one: a tree is wood everywhere and a
/// sampled rock is one substance throughout. A box is the first thing that can
/// have an oak lid on a steel frame, and the contact has to read the material
/// of the piece it actually struck.
#[test]
fn a_box_presents_one_solid_per_wall() {
    let (mut w, node, _, _) = a_box();
    w.tree.refine(node);
    let surface = w.surface_of_node(node).clone();
    println!(
        "  surface: {} pieces, {} bytes, {} spheres in the first",
        surface.len(),
        surface.bytes(),
        surface.pieces()[0].hull.len()
    );
    assert_eq!(surface.len(), 6, "six panels should present six solids");
    // Filled, never hollow: each piece is a slab of eight corner spheres, so
    // nothing inside the box reads as interpenetrating anything.
    for p in surface.pieces() {
        assert_eq!(p.hull.len(), 8, "a panel is a filled slab, not a shell");
    }

    // A point in the middle of the box is inside *no* piece, which is what "a
    // void is not represented at all" means when you go and look.
    let inside_a_piece = surface
        .pieces()
        .iter()
        .any(|p| phys::shape::closest(&p.hull, &phys::shape::Hull::sphere(Vec3::ZERO, 1e-6)).is_none_or(|c| c.gap <= 0.0));
    assert!(!inside_a_piece, "the cavity of the box is inside one of its own walls");
}

/// Joining is an act. A plank welded on is a part of the recipe afterwards, and
/// the node it arrived as is gone.
#[test]
fn joining_absorbs_a_node_into_the_recipe() {
    let (mut w, node, oak, glue) = a_box();
    w.tree.refine(node);
    // Promote the lid out, then put it back.
    let spec = w.tree.nodes[node.get()].spec;
    let lid = w.tree.promote(node, 1, spec);
    assert!(!lid.is_none(), "the lid should promote");
    let before = w.tree.nodes[node.get()].morphology.as_ref().unwrap().assembly.as_ref().unwrap().len();
    let live_before = w.tree.live_count();

    let site = w.join(node, lid, glue, 0.03).expect("the lid rejoins");
    let m = w.tree.nodes[node.get()].morphology.as_ref().unwrap();
    let a = m.assembly.as_ref().unwrap();
    println!(
        "  {} parts -> {} after joining at site {site}; {} live nodes -> {}",
        before,
        a.len(),
        live_before,
        w.tree.live_count()
    );
    assert_eq!(a.len(), before + 1, "the join did not add a part");
    assert!(!w.tree.nodes[lid.get()].alive, "the joined node is still its own object");
    assert!(w.tree.live_count() < live_before, "joining cost no nodes, which is the point of it");
    assert_eq!(a.parts.last().unwrap().substance, oak);
    assert_eq!(a.parts.last().unwrap().join, glue);
    assert!(w.tree.nodes[node.get()].contains_edit, "the join is not recorded as an edit");
}

/// Breaking is the same transform run backwards: the wall becomes a node of its
/// own, carrying its own recipe, and is still a wall.
#[test]
fn a_wall_that_comes_off_is_still_a_wall() {
    let (mut w, node, _, _) = a_box();
    w.tree.refine(node);
    let live_before = w.tree.live_count();

    let wall = w.detach(node, 3);
    assert!(!wall.is_none(), "nothing came off");
    let child = &w.tree.nodes[wall.get()];
    let a = child.morphology.as_ref().expect("the wall has a recipe of its own");
    let parts = a.assembly.as_ref().expect("and it is an assembly");
    println!(
        "  detached: {} live nodes -> {}, the piece is {} part(s), {:.1} kg, {:.3} m",
        live_before,
        w.tree.live_count(),
        parts.len(),
        child.matter.mass,
        child.matter.radius
    );
    assert_eq!(parts.len(), 1, "a wall that comes off should be one part, not a lump");
    assert!(child.matter.mass > 0.0);

    // The box still describes six parts, five of them joined: the recipe now
    // says five walls plus a break.
    let box_parts = w.tree.nodes[node.get()].morphology.as_ref().unwrap().assembly.as_ref().unwrap();
    assert_eq!(box_parts.len(), 6);
    assert_eq!(box_parts.parts.iter().filter(|p| p.joined_to != UNJOINED).count(), 4);
    let events = &w.tree.nodes[node.get()].morphology.as_ref().unwrap().events;
    assert_eq!(events.len(), 1, "the break is not in the recipe's event log");
    assert_eq!(events[0].site, 3);
}

/// Struck hard enough, a wall comes away — and only then.
///
/// The done-when's middle clause, and the reason the seam is stated rather than
/// optimised: a box whose joins are sized to carry whatever they are asked to
/// carry never comes apart, which is a box that is not made of anything.
#[test]
fn a_hard_enough_strike_takes_a_wall_off() {
    for (speed, expect) in [(20.0f64, false), (700.0, true)] {
        let (mut w, node, _, _) = a_box();
        let out = w.damage(node, &[phys::solvers::structure::weather::wind(speed, v3(1.0, 0.0, 0.0))]);
        let a = w.tree.nodes[node.get()].morphology.as_ref().unwrap().assembly.as_ref().unwrap();
        let still_on = a.parts.iter().filter(|p| p.joined_to != UNJOINED).count();
        println!(
            "  {speed:>5.0} m/s: peak utilisation {:.2}, {} joints broke, {} pieces away, {still_on}/5 still joined",
            out.peak_utilisation, out.broken_joints, out.detached_pieces
        );
        if expect {
            assert!(out.detached_pieces > 0, "{speed} m/s took nothing off");
            assert!(still_on < 5, "the recipe still says every wall is attached");
            assert!(out.detached_mass > 0.0);
        } else {
            assert_eq!(out.detached_pieces, 0, "{speed} m/s should not dismantle a crate");
            assert_eq!(still_on, 5);
        }
    }
}

/// With nothing watching, the box is its recipe — and comes back the same box.
#[test]
fn a_box_nobody_is_watching_is_its_recipe() {
    let (mut w, node, _, _) = a_box();
    let before: Vec<Vec3> = w.tree.refine(node).iter().map(|b| b.pos).collect();
    let detail = w.tree.nodes[node.get()].detail_bytes();
    let recipe = w.tree.nodes[node.get()].morphology.as_ref().unwrap().state_bytes();

    assert!(w.collapsible(node), "a described box must be able to release its detail");
    let filed = w.tree.persisted.len();
    w.tree.coarsen(node);
    assert!(w.tree.nodes[node.get()].bodies.is_empty(), "it kept its bodies");
    assert_eq!(w.tree.persisted.len(), filed, "it filed a body list it can regenerate");

    let after: Vec<Vec3> = w.tree.refine(node).iter().map(|b| b.pos).collect();
    let worst = before
        .iter()
        .zip(&after)
        .map(|(a, b)| (*a - *b).norm())
        .fold(0.0f64, f64::max);
    println!("  detail {detail} B -> recipe {recipe} B; regenerated to {worst:.3e} m");
    assert_eq!(before.len(), after.len());
    assert_eq!(worst, 0.0, "the box came back a different box");
    assert!(recipe < detail, "the recipe is not smaller than the detail it replaces");
}
