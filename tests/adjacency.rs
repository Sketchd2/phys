//! Which things are next to which.
//!
//! `docs/BACKLOG.md`'s coupling audit found that contact, heat conduction, mass
//! diffusion, sibling-to-sibling radiation, friction and debris landing on
//! anything but its own parent are one missing primitive rather than six
//! missing features. `docs/PLAY.md` D3 is the design; this is what holds it to
//! it.
//!
//! The property that matters most here is that there is no per-tier variant.
//! Everything is in the node's own frame and its own units, so the same code
//! has to index a galaxy's arms and a nucleus's nucleons — and a test that only
//! ever ran at one scale would not notice the day someone added a constant with
//! metres baked into it.

use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::neighbourhood::{Neighbourhood, Occupant};
use phys::state::{Body, Composition};
use phys::units::*;

/// A ring of `n` bodies of radius `r` on a circle of radius `ring`, in a node
/// frame. Deliberately parameterised by scale so one test can run at galactic
/// and nuclear sizes without a second copy.
fn ring(n: usize, ring: f64, r: f64) -> Vec<Body> {
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            Body {
                pos: v3(ring * a.cos(), ring * a.sin(), 0.0),
                vel: v3(0.0, 0.0, 0.0),
                mass: 1.0,
                radius: r,
                charge: 0.0,
                internal_energy: 0.0,
                temperature: 100.0,
                composition: Composition::primordial(),
                spin: v3(0.0, 0.0, 0.0),
                ..Default::default()
            }
        })
        .collect()
}

fn no_children(_: NodeIdx) -> Option<(phys::math::Vec3, f64)> {
    None
}

// ---------------------------------------------------------------------------
// scale
// ---------------------------------------------------------------------------

/// The same arrangement, at eight scales spanning thirty-five orders of
/// magnitude, must give the same answer *and* partition the space the same way.
///
/// This is the claim D3 rests on — "there is no per-tier variant and there must
/// not be one" — and the first version of this test could not fail. It checked
/// only the answer, and the answer is filtered by true distance after the grid
/// lookup, so a length baked into the spacing degrades efficiency without
/// changing results. Verified: clamping the spacing to one metre left every
/// assertion passing while collapsing thirty orders of magnitude into a single
/// cell. So this asserts on the partition too.
#[test]
fn adjacency_is_the_same_at_every_scale() {
    let mut answers = Vec::new();
    for exponent in [-15.0, -10.0, -6.0, -2.0, 0.0, 6.0, 12.0, 20.0] {
        let scale = 10f64.powf(exponent);
        let bodies = ring(12, scale, 0.05 * scale);
        // Adjacent bodies on a twelve-ring of radius R are 2R sin(pi/12) apart,
        // which is 0.518 R; their surfaces are 0.468 R apart. The next ones out
        // are 1.0 R apart. So a query at 0.6 R finds each body's two immediate
        // neighbours and nothing further.
        let n = Neighbourhood::build(&bodies, &[], 1.0 * scale, 0, no_children);
        let mut found = n.near(bodies[0].pos, 0.6 * scale).expect("query is within reach");
        // The grid answers in cell-visit order, which is deterministic but not
        // slot order; `the_answer_is_deterministic` covers the order itself.
        found.sort();

        // Brute force, at this scale, from the same geometry. If the grid ever
        // returns something different the index is wrong rather than merely
        // slow — and this is what a spacing too *small* for the scale breaks.
        let mut expected: Vec<_> = bodies
            .iter()
            .enumerate()
            .filter(|(_, b)| (b.pos - bodies[0].pos).norm() - b.radius <= 0.6 * scale)
            .map(|(i, _)| Occupant::Body(i as u32))
            .collect();
        expected.sort();
        assert_eq!(
            found, expected,
            "at 10^{exponent} the grid disagreed with brute force"
        );

        answers.push((exponent, found.len(), n.cells()));
    }

    let (_, first_found, first_cells) = answers[0];
    assert_eq!(first_found, 3, "each body should find itself and its two neighbours");
    assert!(first_cells > 1, "the grid should partition, not hold one bucket");
    for (exponent, found, cells) in &answers {
        assert_eq!(
            *found, first_found,
            "scale 10^{exponent} found {found} where 10^{} found {first_found}",
            answers[0].0
        );
        assert_eq!(
            *cells, first_cells,
            "scale 10^{exponent} partitioned into {cells} cells where 10^{} used \
             {first_cells} — a length that is only right at one tier has crept into \
             the spacing",
            answers[0].0
        );
    }
}

// ---------------------------------------------------------------------------
// what a node holds
// ---------------------------------------------------------------------------

/// A promoted child stands in for the body in its slot, and is not counted
/// twice.
///
/// `children` runs parallel to `bodies`, so indexing both would make a thing
/// its own neighbour and double every quantity that crossed a boundary.
#[test]
fn a_promoted_child_replaces_its_body_rather_than_joining_it() {
    let mut w = World::new(galaxy(0xADDA, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 64;
    let root = w.tree.root;
    w.tree.refine(root);
    let bodies = w.tree.nodes[root.get()].bodies.len();

    let before = w.tree.neighbourhood(root);
    assert_eq!(before.len(), bodies, "every body is an occupant");
    assert!(
        before.occupants().iter().all(|o| matches!(o, Occupant::Body(_))),
        "nothing is promoted yet"
    );

    let tier = w.tree.nodes[root.get()].tier;
    let child = w.tree.promote(root, 3, default_spec(tier.finer()));
    assert!(!child.is_none());

    let after = w.tree.neighbourhood(root);
    assert_eq!(
        after.len(),
        bodies,
        "promotion changed the resolution of one slot, not the number of things"
    );
    let promoted: Vec<_> = after
        .occupants()
        .iter()
        .filter(|o| matches!(o, Occupant::Child(_)))
        .collect();
    assert_eq!(promoted.len(), 1, "exactly the promoted slot is a child");
    assert!(
        after.occupants().contains(&Occupant::Child(child)),
        "and it is the child that was promoted"
    );
    assert!(
        !after.occupants().contains(&Occupant::Body(3)),
        "slot 3's body must not appear alongside its own child"
    );
}

/// An unmaterialised node has nothing next to anything, and says so rather than
/// inventing an empty answer that looks like a full one.
#[test]
fn a_node_with_no_contents_has_an_empty_neighbourhood() {
    let w = World::new(galaxy(0xADDB, 1e9), 20.0);
    let n = w.tree.neighbourhood(w.tree.root);
    assert!(n.is_empty());
    assert_eq!(n.near(v3(0.0, 0.0, 0.0), 0.0).map(|v| v.len()), Some(0));
}

// ---------------------------------------------------------------------------
// the query
// ---------------------------------------------------------------------------

/// A query wider than the index can answer is refused, not silently truncated.
///
/// The 27-cell walk finds everything within one cell and nothing beyond it, so
/// a wider question would come back short. A short answer to a question about
/// what is nearby is how a thing falls through a floor.
#[test]
fn a_query_wider_than_the_index_is_refused() {
    let bodies = ring(20, 1.0, 0.01);
    let n = Neighbourhood::build(&bodies, &[], 0.2, 0, no_children);
    assert!(n.near(bodies[0].pos, n.reach() * 0.5).is_some());
    assert!(
        n.near(bodies[0].pos, n.reach() * 1.5).is_none(),
        "a query beyond the grid's reach must be refused rather than answered short"
    );
}

/// Size counts on the occupant's side: a large thing is found by a query that
/// would have missed its centre.
///
/// The spacing has to be at least twice the widest occupant, or the 27-cell
/// walk cannot reach from a small thing to a large one whose surface is close
/// but whose centre is cells away. The geometry here is chosen so the widest
/// term is what carries it: the node's resolution alone would put these five
/// cells apart and find nothing.
#[test]
fn a_large_occupant_is_found_by_its_surface_not_its_centre() {
    let mut bodies = ring(2, 1.0, 0.01);
    bodies[0].pos = v3(0.0, 0.0, 0.0);
    bodies[1].pos = v3(2.5, 0.0, 0.0);
    bodies[1].radius = 2.4;

    let resolution = 0.5;
    let n = Neighbourhood::build(&bodies, &[], resolution, 0, no_children);
    assert!(
        n.reach() > resolution * 2.0,
        "the widest occupant, not the resolution, should have set the spacing: \
         reach {} against resolution {resolution}",
        n.reach()
    );

    let from = bodies[0].pos;
    let centre_gap = (bodies[1].pos - from).norm();
    let surface_gap = centre_gap - bodies[1].radius;
    assert!(centre_gap > 2.0 && surface_gap < 0.2, "precondition: far centre, near surface");

    let found = n.near(from, 0.2).expect("within reach");
    assert!(
        found.contains(&Occupant::Body(1)),
        "the wide body's surface is {surface_gap:.3} away even though its centre is \
         {centre_gap:.3}"
    );
}

/// Overlap is distinct from nearness, because contact is impulsive and
/// conduction is not.
#[test]
fn touching_finds_overlap_and_not_mere_proximity() {
    let mut bodies = ring(3, 1.0, 0.1);
    // Slot 1 moved to overlap slot 0; slot 2 left where it is.
    bodies[1].pos = bodies[0].pos + v3(0.15, 0.0, 0.0);
    let n = Neighbourhood::build(&bodies, &[], 2.0, 0, no_children);

    let touching = n.touching(0);
    assert!(touching.contains(&Occupant::Body(1)), "0 and 1 interpenetrate");
    assert!(!touching.contains(&Occupant::Body(2)), "2 is near but not touching");
    assert!(!touching.contains(&Occupant::Body(0)), "nothing touches itself");
}

/// The index knows when the node it describes has moved on.
#[test]
fn a_neighbourhood_knows_when_it_is_stale() {
    let bodies = ring(4, 1.0, 0.01);
    let n = Neighbourhood::build(&bodies, &[], 1.0, 7, no_children);
    assert!(n.is_current(7));
    assert!(!n.is_current(8), "a node whose epoch moved has different contents");
}

/// Determinism, which the whole engine rests on: the same contents give the
/// same answer in the same order, however the grid's buckets happen to hash.
#[test]
fn the_answer_is_deterministic() {
    let bodies = ring(64, 1.0, 0.02);
    let first = {
        let n = Neighbourhood::build(&bodies, &[], 0.5, 0, no_children);
        n.near(bodies[0].pos, 0.4).unwrap()
    };
    for _ in 0..8 {
        let n = Neighbourhood::build(&bodies, &[], 0.5, 0, no_children);
        assert_eq!(n.near(bodies[0].pos, 0.4).unwrap(), first);
    }
}

/// A node whose contents are all at one point still answers, rather than
/// dividing by a zero spread.
#[test]
fn coincident_contents_do_not_break_the_index() {
    let mut bodies = ring(8, 1.0, 0.01);
    for b in bodies.iter_mut() {
        b.pos = v3(0.0, 0.0, 0.0);
    }
    let n = Neighbourhood::build(&bodies, &[], 1.0, 0, no_children);
    let found = n.near(v3(0.0, 0.0, 0.0), 0.5).expect("within reach");
    assert_eq!(found.len(), 8, "all eight are at the query point");
}

/// A node whose radius is zero, or whose bodies carry a non-finite position,
/// must not take the index with them. Both have happened.
#[test]
fn degenerate_geometry_is_survived() {
    let mut bodies = ring(4, 1.0, 0.01);
    bodies[2].pos = v3(f64::NAN, 0.0, 0.0);
    bodies[3].pos = v3(f64::INFINITY, 0.0, 0.0);
    let n = Neighbourhood::build(&bodies, &[], 0.0, 0, no_children);
    let found = n.near(bodies[0].pos, n.reach()).expect("within reach");
    assert!(
        !found.contains(&Occupant::Body(2)) && !found.contains(&Occupant::Body(3)),
        "a body with no position is nobody's neighbour"
    );
}
