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

// ---------------------------------------------------------------------------
// Exchange: the transport half of D3.
// ---------------------------------------------------------------------------

use phys::neighbourhood::{exchange, radiative_area, radiative_conductance, Reservoir};

/// Heat moves the way it is asked to, and the two sides never cross.
///
/// The crossing is the whole reason `exchange` solves the pair rather than
/// multiplying `G * delta * dt`: the engine's steps are seconds long and two
/// small things in contact equilibrate in microseconds, so the explicit form
/// overshoots and then oscillates with growing amplitude.
#[test]
fn an_exchange_never_carries_the_two_sides_past_each_other() {
    let (ca, cb) = (10.0, 4.0);
    // A conductance and a span whose product dwarfs both capacities: the
    // explicit form would move about 30,000 times the whole difference.
    for dt in [1e-6, 1e-3, 1.0, 1e3, 1e9] {
        let a = Reservoir::new(1000.0, ca);
        let b = Reservoir::new(100.0, cb);
        let q = exchange(a, b, 1e5, dt);
        let (ta, tb) = (a.potential - q / ca, b.potential + q / cb);
        assert!(
            ta >= tb - 1e-9,
            "dt={dt}: the sides crossed — a ended at {ta} and b at {tb}, \
             so {q} J moved where at most the difference should have"
        );
        assert!(ta <= 1000.0 + 1e-9 && tb >= 100.0 - 1e-9, "dt={dt}: {ta} {tb}");
    }
}

/// Given long enough, they meet at the capacity-weighted mean, which is where
/// energy conservation says they have to.
#[test]
fn an_exchange_settles_at_the_weighted_mean() {
    let (ca, cb) = (10.0, 4.0);
    let a = Reservoir::new(1000.0, ca);
    let b = Reservoir::new(100.0, cb);
    let q = exchange(a, b, 1e5, 1e9);
    let (ta, tb) = (a.potential - q / ca, b.potential + q / cb);
    let mean = (1000.0 * ca + 100.0 * cb) / (ca + cb);
    assert!((ta - mean).abs() < 1e-6, "a settled at {ta}, not {mean}");
    assert!((tb - mean).abs() < 1e-6, "b settled at {tb}, not {mean}");
}

/// For a short step it *is* `G * delta * dt`. The exact solution has to agree
/// with the law it solves in the limit where the law is unambiguous.
#[test]
fn a_short_exchange_is_the_explicit_rate() {
    let a = Reservoir::new(400.0, 1e6);
    let b = Reservoir::new(300.0, 1e6);
    let (g, dt) = (2.0, 1e-3);
    let q = exchange(a, b, g, dt);
    let explicit = g * (a.potential - b.potential) * dt;
    let rel = (q - explicit).abs() / explicit.abs();
    assert!(rel < 1e-6, "short-step exchange moved {q} J, explicit says {explicit}");
}

/// A bath is a reservoir the exchange cannot move, and the pair law has to
/// degrade to the one-body law rather than to zero or to a division by
/// infinity.
#[test]
fn a_bath_does_not_move() {
    let a = Reservoir::new(400.0, 100.0);
    let sky = Reservoir::bath(2.725);
    let q = exchange(a, sky, 1.0, 10.0);
    let one_body = 100.0 * (400.0 - 2.725) * (1.0 - (-1.0 * 10.0 / 100.0f64).exp());
    assert!(q > 0.0, "nothing left the warm side towards a cold bath");
    assert!(
        (q - one_body).abs() / one_body < 1e-9,
        "a bath exchange moved {q} J, the one-body law says {one_body}"
    );
}

/// Whichever way round the pair is named, the same heat crosses the same
/// boundary. A boundary that answered differently depending on which side
/// asked would create or destroy energy the moment both sides asked.
#[test]
fn an_exchange_is_the_same_boundary_from_either_side() {
    let a = Reservoir::new(700.0, 3.0);
    let b = Reservoir::new(120.0, 11.0);
    let there = exchange(a, b, 0.5, 20.0);
    let back = exchange(b, a, 0.5, 20.0);
    assert!(
        (there + back).abs() < 1e-12 * there.abs().max(1.0),
        "a->b moved {there} J but b->a moved {back} J"
    );
    assert!(radiative_area(3.0, 7.0, 40.0) == radiative_area(7.0, 3.0, 40.0));
}

/// The radiative conductance is Stefan-Boltzmann exactly, not linearised about
/// a mean. `T^4 - T^4` factors, so there is no reason to approximate it.
#[test]
fn the_radiative_conductance_reproduces_stefan_boltzmann() {
    for (ta, tb) in [(300.0, 290.0), (6000.0, 3.0), (3.0, 6000.0), (1e4, 1.0)] {
        let area = 2.5;
        let g: f64 = radiative_conductance(ta, tb, area);
        let linear = g * (ta - tb);
        let quartic: f64 =
            phys::units::SIGMA_SB * area * (ta.powi(4) - tb.powi(4));
        let rel = (linear - quartic).abs() / quartic.abs().max(1e-300);
        assert!(
            rel < 1e-12,
            "at {ta} K against {tb} K the conductance gives {linear} W where \
             Stefan-Boltzmann gives {quartic} W"
        );
    }
}

/// The exchange area is the *reciprocal* one, and reciprocity is what makes
/// the transfer conserve: `A_a F_ab` has to equal `A_b F_ba` or the two sides
/// disagree about how big their shared boundary is.
#[test]
fn the_radiative_area_falls_off_as_the_inverse_square() {
    let (ra, rb) = (0.3, 0.7);
    let near = radiative_area(ra, rb, 10.0);
    let far = radiative_area(ra, rb, 20.0);
    assert!(
        (near / far - 4.0).abs() < 1e-9,
        "doubling the distance changed the area by {}x, not 4x",
        near / far
    );
    // And it is bounded: two spheres in contact cannot see more of each other
    // than a hemisphere of the smaller.
    let touching = radiative_area(ra, rb, ra + rb);
    assert!(
        touching <= 2.0 * std::f64::consts::PI * ra * ra + 1e-12,
        "touching spheres exchange over {touching} m^2, more than the smaller has"
    );
    assert!(touching > 0.0);
}

/// A rocky planet refined once, with two of its bodies promoted, the second
/// placed 2.2 radii from the first, and each given a temperature.
///
/// Planetary rather than Continuum deliberately. A Continuum node whose matter
/// carries a chemical cohesive energy has its sampled geometry inflated by
/// about 4.3e5 before `sampler::sample` gives up relaxing it — see the note in
/// `docs/BACKLOG.md` — so nothing at that tier is currently next to anything.
/// This is the finest tier where the geometry is trustworthy today.
fn adjacent_pair(
    ta: f64,
    tb: f64,
) -> (phys::engine::World, phys::ids::NodeIdx, phys::ids::NodeIdx) {
    use phys::engine::{default_spec, World};
    let sc = phys::scenario::ALL
        .iter()
        .find(|s| s.name == "Rocky planet")
        .expect("the scenario shelf has a rocky planet");
    let mut w = World::new(sc.build(0xC0FFEE), 1.0);
    let root = w.tree.root;
    w.tree.nodes[0].spec.count = 64;
    // Paced to its subject. A world runs at one second per second (`PLAY.md`
    // D1) and two planet-sized chunks take centuries to equilibrate radiatively
    // — the conductance is right, the span is what has to be large.
    w.pace_to(root);
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    let a = w.tree.promote(root, 0, default_spec(tier.finer()));
    let b = w.tree.promote(root, 1, default_spec(tier.finer()));
    assert!(!a.is_none() && !b.is_none(), "both promotions should succeed");
    let ra = w.tree.nodes[a.get()].matter.radius;
    let at = w.tree.nodes[a.get()].motion.offset;
    w.tree.nodes[b.get()].motion.offset = at + phys::math::v3(2.2 * ra, 0.0, 0.0);
    w.tree.nodes[a.get()].matter.temperature = ta;
    w.tree.nodes[b.get()].matter.temperature = tb;
    w.tree.pin(a);
    w.tree.pin(b);
    (w, a, b)
}

fn settled(ta: f64, tb: f64) -> (f64, f64) {
    let (mut w, a, b) = adjacent_pair(ta, tb);
    for _ in 0..60 {
        w.step_frame(50_000.0);
    }
    (
        w.tree.nodes[a.get()].matter.temperature,
        w.tree.nodes[b.get()].matter.temperature,
    )
}

/// The Phase 1 headline: a hot node beside a cold one equilibrates without
/// either being told the other exists.
///
/// Measured against a control, and the *choice* of control is the whole design
/// of this test. Both nodes also run `evolve_matter`, which radiates to the sky
/// and absorbs the parent's light, so "the cold one got warmer" would have
/// passed with no coupling at all. The first control tried here moved the
/// partner out of range instead — which also moved it to a different distance
/// from the parent, changed its illumination, and reproduced the exact
/// signature the exchange was supposed to produce. That version passed with
/// `exchange` stubbed to return zero, which is how it was caught.
///
/// So the control holds the geometry fixed to the last bit and varies only the
/// partner's temperature. Everything positional — illumination, the solver, the
/// sampler's streams — is then identical between the two runs by construction,
/// and the only thing that can separate them is heat crossing the boundary.
#[test]
fn a_hot_node_beside_a_cold_one_equilibrates() {
    // The cold node, beside a hot partner and beside a cold one.
    let (_, warmed) = settled(6000.0, 50.0);
    let (_, alone) = settled(50.0, 50.0);
    // And the assertion is on the *ratio*, not on the sign of the difference,
    // because the control is not perfectly clean and saying so is cheaper than
    // pretending. Raising the partner to 6000 K also changes what the partner
    // hands its parent, so `environment_at` gives this node marginally more
    // light in one run than in the other, and that leak has the same sign as
    // the effect being measured. With `exchange` stubbed to return zero it is
    // worth 2.4e-8 K against the exchange's 2.5e-3 K — five orders of magnitude
    // down, but strictly positive, so `warmed > alone` passed with no coupling
    // at all. This test was written that way first, and that is how it was
    // caught.
    //
    // Against the background warming both runs get from the parent, the
    // measured separation is 4.95x with the exchange and 1.07x without it.
    let (gained, background) = (warmed - 50.0, alone - 50.0);
    let ratio = gained / background.max(1e-300);
    assert!(
        ratio > 2.0,
        "a 50 K node warmed by {gained} K beside a 6000 K neighbour and by \
         {background} K beside a 50 K one, a factor of {ratio} — where the same \
         pair with no coupling at all comes out at 1.07x. Nothing crossed."
    );

    // The hot side is deliberately *not* tested by the mirror of this control.
    // Raising the partner from 50 K to 6000 K does not hold the rest of the
    // system fixed: a 6000 K node radiates through `evolve_matter` into
    // everything around it, and the parent's state comes back changed. The run
    // was measured and the hot node cooled *more* beside a hot partner than
    // beside a cold one — a real effect of that feedback, and nothing to do
    // with the boundary being tested. What leaves the hot side is measured
    // directly instead, in `heat_crosses_a_boundary_and_the_books_close`, where
    // there is no second mechanism in the way.
}

/// Heat leaves one side, arrives at the other, and the sum does not move.
///
/// Between two *bodies* in one node rather than two promoted children, which is
/// what makes it a clean measurement: `evolve_matter` acts on a node's matter
/// and not on the bodies inside it, so there is no second path by which energy
/// can enter or leave the list being summed. Whatever the two temperatures do
/// to each other, they did to each other.
///
/// The conservation bound is the real assertion here. An exchange that moved
/// heat in the right direction while quietly minting or losing some of it would
/// satisfy every other test in this file.
#[test]
fn heat_crosses_a_boundary_and_the_books_close() {
    use phys::engine::World;
    let sc = phys::scenario::ALL
        .iter()
        .find(|s| s.name == "Rocky planet")
        .expect("the scenario shelf has a rocky planet");
    let mut w = World::new(sc.build(0xC0FFEE), 1.0);
    let root = w.tree.root;
    w.tree.nodes[0].spec.count = 64;
    // Paced to its subject, as above: a second of radiative exchange between
    // planet-sized bodies moves nothing measurable.
    w.pace_to(root);
    w.tree.refine(root);
    // Without a pin the bodies are discarded the moment nobody is looking, and
    // there is nothing left to measure.
    w.tree.pin(root);

    let nb = w.tree.neighbourhood(root);
    let pairs = nb
        .pairs(nb.resolution())
        .expect("the node's own resolution is within its index's reach");
    assert!(!pairs.is_empty(), "a refined planet's bodies are next to nothing");
    let (i, j) = pairs[0];

    {
        let b = &mut w.tree.nodes[root.get()].bodies;
        b[i].temperature = 9000.0;
        b[j].temperature = 30.0;
    }
    let before: f64 = w.tree.nodes[root.get()].bodies.iter().map(|b| b.internal_energy).sum();

    for _ in 0..5 {
        w.step_frame(50_000.0);
    }

    let b = &w.tree.nodes[root.get()].bodies;
    assert_eq!(b.len(), 64, "the pinned node lost its detail");
    let (hot, cold) = (b[i].temperature, b[j].temperature);
    assert!(hot < 9000.0, "the hot body never cooled: still {hot} K");
    assert!(cold > 30.0, "the cold body never warmed: still {cold} K");
    assert!(hot > cold, "the pair crossed: {hot} K against {cold} K");

    let after: f64 = b.iter().map(|b| b.internal_energy).sum();
    let drift = (after - before).abs() / before.abs();
    assert!(
        drift < 1e-12,
        "the exchange moved heat but the books did not close: total internal          energy went from {before} J to {after} J, a relative drift of {drift}"
    );
}

/// An exchange does not pin the node it lands on.
///
/// `Tree::pin` is one-way and pins the whole ancestry, so pinning on ordinary
/// transport would pin every node with a warm neighbour, permanently. That is
/// axiom four — detail exists where something is happening — exactly inverted.
#[test]
fn an_exchange_does_not_pin_what_it_touches() {
    use phys::causal::InfluenceKind;
    use phys::engine::{default_spec, World};
    use phys::math::Vec3;
    let sc = phys::scenario::ALL
        .iter()
        .find(|s| s.name == "Rocky planet")
        .expect("the scenario shelf has a rocky planet");
    let mut w = World::new(sc.build(0xBEEF), 1.0);
    let root = w.tree.root;
    w.tree.nodes[0].spec.count = 64;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    let child = w.tree.promote(root, 0, default_spec(tier.finer()));
    assert!(!w.tree.nodes[child.get()].pinned, "a fresh promotion should not be pinned");

    let before = w.tree.nodes[child.get()].matter.internal_energy;
    w.mailbox
        .post(child, w.time, 0.0, InfluenceKind::Exchange, 1.0e20, Vec3::ZERO);
    w.step_frame(50_000.0);
    assert!(
        w.tree.nodes[child.get()].matter.internal_energy > before,
        "the exchange was never delivered"
    );
    assert!(
        !w.tree.nodes[child.get()].pinned,
        "an exchange pinned the node it landed on"
    );

    // The contrast, and the reason the kind is separate: a user's impulse is
    // information from outside that no re-sampling reproduces, and it does pin.
    w.mailbox
        .post(child, w.time, 0.0, InfluenceKind::UserImpulse, 1.0, Vec3::ZERO);
    w.step_frame(50_000.0);
    assert!(
        w.tree.nodes[child.get()].pinned,
        "a user impulse should still pin, and now does not"
    );
}

// ---------------------------------------------------------------------------
// Contact: the impulsive half of D3.
// ---------------------------------------------------------------------------

use phys::neighbourhood::{contact, restitution, yield_velocity, Side, Surface};
use phys::topology::Material;

fn side(x: f64, vx: f64, m: f64, r: f64, mat: &Material) -> Side {
    Side {
        pos: v3(x, 0.0, 0.0),
        velocity: v3(vx, 0.0, 0.0),
        mass: m,
        radius: r,
        heat_capacity: 1000.0,
        surface: Surface::of(mat),
    }
}

/// Restitution is a function of the impact, not a property of the material.
///
/// This is the whole reason D3 asks for it to be derived. A tabulated
/// coefficient is right at one speed; the same pair of surfaces returns a third
/// of a slow approach and a seventh of a fast one, and no single number is both.
#[test]
fn restitution_falls_with_the_speed_of_the_impact() {
    let wood = Surface::of(&Material::GREEN_WOOD);
    let slow = restitution(&wood, &wood, 1.0);
    let fast = restitution(&wood, &wood, 20.0);
    assert!(
        slow > fast,
        "wood returned {slow} of a 1 m/s approach and {fast} of a 20 m/s one"
    );
    assert!((0.0..=1.0).contains(&slow) && (0.0..=1.0).contains(&fast));
    // Below the yield velocity nothing is lost, because nothing has yielded.
    let v_y = yield_velocity(&wood, &wood);
    assert!(v_y > 0.0, "wood has no yield velocity at all");
    assert_eq!(restitution(&wood, &wood, v_y * 0.5), 1.0);
}

/// Masonry does not bounce and a steel frame does, at the same speed, without
/// either having been told so.
#[test]
fn a_stiff_strong_surface_returns_more_than_a_weak_one() {
    let frame = Surface::of(&Material::REINFORCED_FRAME);
    let masonry = Surface::of(&Material::MASONRY);
    let (a, b) = (
        restitution(&frame, &frame, 5.0),
        restitution(&masonry, &masonry, 5.0),
    );
    assert!(
        a > 4.0 * b,
        "a reinforced frame returned {a} and masonry {b} at the same 5 m/s"
    );
}

/// Nothing without a surface collides with anything.
///
/// A gas parcel and a star cluster are things a node holds, and giving them a
/// surface so that the collision code has something to read would be the engine
/// being *told* they are solid.
#[test]
fn a_thing_with_no_surface_does_not_collide() {
    let steel = Material::STEEL;
    let mut ghost = side(1.9, -3.0, 125.0, 1.0, &steel);
    ghost.surface = Surface { density: 0.0, stiffness: 0.0, strength: 0.0 };
    assert!(contact(&side(0.0, 3.0, 125.0, 1.0, &steel), &ghost).is_none());
    // And a pair already separating is left alone, or two overlapping things
    // buzz against each other forever.
    assert!(contact(
        &side(0.0, -3.0, 125.0, 1.0, &steel),
        &side(1.9, 3.0, 125.0, 1.0, &steel)
    )
    .is_none());
}

/// The contact conserves momentum and energy, and the energy it does not
/// return is accounted for as heat rather than dropped.
#[test]
fn a_contact_closes_its_books() {
    let m = Material::STEEL;
    let (a, b) = (side(0.0, 3.0, 125.0, 1.0, &m), side(1.9, -3.0, 200.0, 1.0, &m));
    let c = contact(&a, &b).expect("an approaching overlap is a contact");

    // One impulse, applied with both signs.
    let total = c.normal + c.friction;
    let (va, vb) = (
        a.velocity - total.scale(1.0 / a.mass),
        b.velocity + total.scale(1.0 / b.mass),
    );
    let before = a.velocity.scale(a.mass) + b.velocity.scale(b.mass);
    let after = va.scale(a.mass) + vb.scale(b.mass);
    assert!(
        (after - before).norm() < 1e-9,
        "momentum went from {before:?} to {after:?}"
    );

    let ke_before = 0.5 * a.mass * a.velocity.norm2() + 0.5 * b.mass * b.velocity.norm2();
    let ke_after = 0.5 * a.mass * va.norm2() + 0.5 * b.mass * vb.norm2();
    assert!(
        ke_after <= ke_before + 1e-9,
        "the contact created energy: {ke_before} J became {ke_after} J"
    );
    let lost = ke_before - ke_after;
    let heat = c.heat_a + c.heat_b;
    assert!(
        (heat - lost).abs() / lost.max(1e-30) < 1e-9,
        "{lost} J left the motion and {heat} J arrived as heat"
    );
    // Split so both sides rise by the same temperature: equal capacities here,
    // so equal shares.
    assert!((c.heat_a - c.heat_b).abs() < 1e-9);
}

/// The friction couple conserves angular momentum about any point.
///
/// The obvious spelling of it does not, and only stops doing so once the two
/// are interpenetrating rather than just touching — which is the only state the
/// engine ever actually sees. Measured at a 1.3% leak per contact before the
/// contact point was made a single point shared by both sides.
#[test]
fn the_friction_couple_conserves_angular_momentum() {
    let m = Material::STEEL;
    let mut a = side(0.0, 3.0, 125.0, 1.0, &m);
    let mut b = side(1.9, -3.0, 125.0, 1.0, &m);
    // Sliding as well as closing, or there is no friction to test.
    a.velocity = v3(3.0, 2.0, 0.0);
    b.velocity = v3(-3.0, -1.0, 0.0);
    let c = contact(&a, &b).expect("an approaching overlap is a contact");
    assert!(c.friction.norm() > 0.0, "no friction acted, so nothing is under test");

    let total = c.normal + c.friction;
    let l_before = a.pos.cross(a.velocity.scale(a.mass)) + b.pos.cross(b.velocity.scale(b.mass));
    let va = a.velocity - total.scale(1.0 / a.mass);
    let vb = b.velocity + total.scale(1.0 / b.mass);
    let l_after = a.pos.cross(va.scale(a.mass))
        + b.pos.cross(vb.scale(b.mass))
        + c.spin_a
        + c.spin_b;
    let scale = l_before.norm().max(l_after.norm()).max(1e-30);
    assert!(
        (l_after - l_before).norm() / scale < 1e-12,
        "angular momentum went from {l_before:?} to {l_after:?}"
    );
}

/// Friction cannot exceed the Coulomb limit, and stops the slide when it can
/// afford to.
#[test]
fn friction_is_bounded_by_the_normal_impulse() {
    let m = Material::STEEL;
    let mu = phys::neighbourhood::friction(&Surface::of(&m), &Surface::of(&m));
    // Sliding far faster than it is closing: friction saturates.
    let mut a = side(0.0, 0.2, 125.0, 1.0, &m);
    let mut b = side(1.9, -0.2, 125.0, 1.0, &m);
    a.velocity = v3(0.2, 50.0, 0.0);
    b.velocity = v3(-0.2, -50.0, 0.0);
    let c = contact(&a, &b).expect("an approaching overlap is a contact");
    let ratio = c.friction.norm() / c.normal.norm();
    assert!(
        ratio <= mu + 1e-12,
        "friction reached {ratio} of the normal impulse, above the limit of {mu}"
    );
    assert!(ratio > mu * 0.999, "friction should have saturated, and reached {ratio}");
}

/// Two promoted things with materials, overlapping and closing, in a node whose
/// own solver is quiet enough that what happens is the contact.
///
/// `Tier::Galactic` at ten metres, deliberately. A tier is a physics regime and
/// not a size, so a small node can carry the collisionless gravity solver — and
/// a hundred kilograms two metres apart pull on each other at 10^-9 m/s^2,
/// which leaves the contact as the only thing in the measurement.
///
/// **`Tier::Continuum` is where a metre-scale solid belongs and it cannot be
/// used yet**, which is `PLAY.md` §3.3 rather than anything about this test:
/// "`solvers::for_tier(Continuum)` is `Hydro`. A building, a wolf and a boulder
/// are all `Continuum`, and none of them is a fluid." SPH reads `Matter` through
/// a gas equation of state, and applied to condensed matter it answers with
/// pressures nothing can hold — measured, 1.0x10^8 Pa for this box and
/// 1.7x10^9 Pa for a bucket of water, against a tensile strength of 4.5x10^7 Pa
/// for green wood. The box bursts from its own equation of state before
/// anything touches it. §3.3's state-aware dispatch is the fix and is a Phase 1
/// item that has not been built.
///
/// The other alternative — the same test at galactic *distances* — silently
/// placed both children at the same point, because 50 m added to 2x10^20 m is
/// below what an `f64` can represent.
fn colliding_pair(approaching: bool) -> (phys::engine::World, phys::ids::NodeIdx, phys::ids::NodeIdx) {
    use phys::engine::{default_spec, World};
    use phys::morph::Program;
    use phys::sampler::{MassSpectrum, Profile, SampleSpec};
    use phys::state::{BodyKind, Composition, Matter};
    use phys::tree::Tree;

    let matter = Matter::neutral(1000.0, 10.0, 290.0, Composition::primordial());
    let spec = SampleSpec::new(8, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let mut w = World::new(Tree::new(0xC07AC7, matter, Tier::Galactic, spec), 1.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    let a = w.tree.promote(root, 0, default_spec(tier.finer()));
    let b = w.tree.promote(root, 1, default_spec(tier.finer()));
    // A structure is what has a material, and a material is what has a surface.
    w.emplace(a, Program::Tower, 1.0, None);
    w.emplace(b, Program::Tower, 1.0, None);
    // Sized by hand. A promoted child that comes out larger than the parent it
    // was promoted from is a separate defect and not what this measures.
    w.tree.nodes[a.get()].matter.radius = 1.0;
    w.tree.nodes[b.get()].matter.radius = 1.0;
    let at = w.tree.nodes[a.get()].motion.offset;
    w.tree.nodes[b.get()].motion.offset = at + v3(1.9, 0.0, 0.0);
    let s = if approaching { 1.0 } else { -1.0 };
    w.tree.nodes[a.get()].motion.velocity = v3(3.0 * s, 0.0, 1.0);
    w.tree.nodes[b.get()].motion.velocity = v3(-3.0 * s, 0.0, 0.0);
    w.tree.pin(a);
    w.tree.pin(b);
    (w, a, b)
}

/// The Phase 1 headline: two promoted things collide and rebound at a
/// restitution derived from what they are made of.
#[test]
fn two_promoted_things_collide_and_rebound() {
    let (mut w, a, b) = colliding_pair(true);
    let (ma, mb) = (
        w.tree.nodes[a.get()].matter.mass,
        w.tree.nodes[b.get()].matter.mass,
    );
    let closing_before = (w.tree.nodes[b.get()].motion.velocity
        - w.tree.nodes[a.get()].motion.velocity)
        .x;
    assert!(closing_before < 0.0, "the pair should start out approaching");

    w.advance_node(w.tree.root, 1.0);

    assert_eq!(w.stats.contacts_resolved, 1, "the overlap was never resolved");
    let closing_after = (w.tree.nodes[b.get()].motion.velocity
        - w.tree.nodes[a.get()].motion.velocity)
        .x;
    assert!(
        closing_after > 0.0,
        "they were closing at {closing_before} m/s and are still closing at \
         {closing_after} m/s"
    );
    assert!(
        closing_after < -closing_before,
        "they separated at {closing_after} m/s having approached at \
         {closing_before} m/s, which is more than was put in"
    );

    // And the rebound is the one the materials say, not a number picked to make
    // this pass. Both are reinforced frame; the separation ratio is the
    // restitution at the speed they met at.
    let surface = Surface::of(&Material::REINFORCED_FRAME);
    let expected = restitution(&surface, &surface, closing_before);
    let measured = closing_after / -closing_before;
    assert!(
        (measured - expected).abs() < 1e-3,
        "they rebounded at {measured} of the approach where the materials say \
         {expected}"
    );
    let _ = (ma, mb);
}

/// And the contact conserves what a contact has to conserve, measured against
/// the same pair passing without touching.
///
/// The control is the same geometry with both velocities reversed, so the
/// solver sees identical masses at identical positions and its own residual —
/// which is not zero — appears in both runs. Whatever separates them is the
/// contact.
#[test]
fn a_collision_in_a_running_world_conserves_momentum_and_spin() {
    let mut momenta = Vec::new();
    let mut angular = Vec::new();
    for approaching in [true, false] {
        let (mut w, a, b) = colliding_pair(approaching);
        let total = |w: &phys::engine::World| {
            let (na, nb) = (&w.tree.nodes[a.get()], &w.tree.nodes[b.get()]);
            let pa = na.motion.velocity.scale(na.matter.mass);
            let pb = nb.motion.velocity.scale(nb.matter.mass);
            (
                pa + pb,
                na.motion.offset.cross(pa)
                    + na.matter.spin
                    + nb.motion.offset.cross(pb)
                    + nb.matter.spin,
            )
        };
        let (p0, l0) = total(&w);
        w.advance_node(w.tree.root, 1.0);
        let (p1, l1) = total(&w);
        assert_eq!(
            w.stats.contacts_resolved,
            u64::from(approaching),
            "a pair moving apart must not be resolved as a contact"
        );
        momenta.push((p1 - p0).norm() / p0.norm());
        angular.push((l1 - l0).norm() / l0.norm().max(1e-300));
    }
    // The contact's own contribution to either is nothing: the run that had one
    // drifts no more than the run that did not.
    assert!(
        momenta[0] < momenta[1] * 2.0 + 1e-12,
        "the collision moved the pair's momentum by {} against {} for the same \
         pair passing untouched",
        momenta[0],
        momenta[1]
    );
    assert!(
        angular[0] < angular[1] * 2.0 + 1e-12,
        "the collision moved the pair's angular momentum by {} against {} for \
         the same pair passing untouched — the friction couple is unbalanced",
        angular[0],
        angular[1]
    );
}

/// Two bodies of the same node are not in contact, they are in it together.
///
/// The hazard is real and was measured rather than assumed — and it is not
/// universal, which is why it is worth pinning. A materialised `Wall` has 527
/// overlapping pairs among 55 members, because a wall is courses of blocks
/// packed against each other; a `Tower`, a `Tree` and a `Settlement` have none,
/// because their members are long and thin and sit a member-length apart.
/// Resolving a wall's courses as collisions would blow it apart on the frame it
/// was materialised.
///
/// The rule is broader than that case, though. A node's own bodies are already
/// coupled by whatever the node is — the structure solver for members of a
/// structure, the tier's solver for parcels of a continuum — so contact is for
/// what that coupling does not reach: a promoted child against another, or
/// against the bodies of the node it is sitting in.
#[test]
fn the_bodies_of_one_node_do_not_collide_with_each_other() {
    use phys::engine::World;
    use phys::morph::Program;
    use phys::sampler::{MassSpectrum, Profile, SampleSpec};
    use phys::state::{BodyKind, Composition, Matter};
    use phys::tree::Tree;

    let matter = Matter::neutral(1.0e5, 10.0, 290.0, Composition::primordial());
    let spec = SampleSpec::new(64, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let mut w = World::new(Tree::new(0x571FF, matter, Tier::Galactic, spec), 1.0);
    let root = w.tree.root;
    // A material for every body in it, which is the condition contact needs.
    // A wall, specifically: it is the program whose members actually overlap.
    w.emplace(root, Program::Wall, 1.0e4, None);
    w.tree.refine(root);
    let n = w.tree.nodes[root.get()].bodies.len();
    assert!(n > 4, "the structure should have members to test, and has {n}");

    let nb = w.tree.neighbourhood(root);
    let pairs = nb.pairs(0.0).expect("zero is within any reach");
    let overlaps = pairs.len();
    assert!(
        overlaps > 0,
        "a wall's courses should overlap each other and none of them do — this \
         test would then prove nothing"
    );

    // And they have to be *moving into* each other, or nothing is under test.
    // A materialised structure sits still, so every pair of its members is
    // refused by `contact` for having no closing speed at all — which is how
    // the first version of this test passed with the rule it exists to check
    // deleted outright.
    {
        let bodies = &mut w.tree.nodes[root.get()].bodies;
        for &(i, j) in pairs.iter() {
            let toward = (bodies[j].pos - bodies[i].pos).unit();
            bodies[i].vel = toward.scale(5.0);
            bodies[j].vel = toward.scale(-5.0);
        }
    }

    for _ in 0..3 {
        w.advance_node(root, 1.0);
    }
    assert_eq!(
        w.stats.contacts_resolved, 0,
        "{overlaps} overlapping pairs inside one structure were resolved as \
         collisions; a tree would come apart on the frame it was materialised"
    );
}

// ---------------------------------------------------------------------------
// A ball in a box
// ---------------------------------------------------------------------------

/// A wooden ball loose inside a wooden box, in deep space.
///
/// This is the contact half asked the question it is actually for: not one
/// impulse in isolation but a thing rattling around inside another thing, over
/// a thousand frames and eight collisions, with nothing else acting on either.
/// Deep space is not decoration — with no planetary gravity the momentum and
/// angular momentum of the whole assembly are *exactly* conserved quantities,
/// so any leak in the contact path has nowhere to hide.
///
/// # The box is a shell of panels, because a node's contents are spheres
///
/// Six spheres cannot enclose a volume: the first version of this used one per
/// face and the ball left through a corner on its second bounce, which is not a
/// defect in anything but the arrangement. So the box is a 4x4 grid of panels
/// per face — 96 of them, each a sphere wide enough that the grid has no holes.
/// That is the honest representation of a container in an engine whose contact
/// primitive is sphere-to-sphere, and it is worth knowing that it is what a box
/// costs.
///
/// # What this does not show, and it matters
///
/// **The box is not rigid.** Nothing in the contact path holds a structure's
/// members together: each panel takes its own impulse onto its own velocity and
/// drifts off with it, so the ball does not bounce off a two-tonne box, it
/// bounces off one 500 kg panel. Measured over 1500 frames, the panels start
/// 3.00 to 3.90 m from the centre and end 3.18 to 5.05 — the box grows by about
/// a third while being knocked around inside. `solvers::structure` has the
/// machinery for a rigid response and `drop_fragments` already uses it, through
/// `Mechanism::PointImpulse` and `damage`; contact does not call it. That gap
/// is in `docs/BACKLOG.md`.
///
/// So what is asserted here is conservation and response — the box does take up
/// the ball's momentum, and off-centre hits do spin it — and not rigidity,
/// which the engine does not yet have.
#[test]
fn a_ball_loose_in_a_box_conserves_momentum_and_angular_momentum() {
    use phys::engine::World;
    use phys::math::Vec3;
    use phys::sampler::{MassSpectrum, Profile, SampleSpec};
    use phys::state::{BodyKind, Composition, Matter};
    use phys::topology::{Material, Topology};
    use phys::tree::Tree;

    const L: f64 = 3.0; // half-width of the box, m
    const K: usize = 4; // panels per edge
    const PANEL_MASS: f64 = 500.0;
    const BALL_MASS: f64 = 10.0;
    const R_BALL: f64 = 0.4;
    const DT: f64 = 0.01;
    const FRAMES: usize = 1500;

    // A cubic shell of panel centres, and the radius that leaves no gap: a
    // square grid at `step` is covered by discs of `step / sqrt(2)`.
    let step = 2.0 * L / K as f64;
    let at = |i: usize| -L + step * (i as f64 + 0.5);
    let mut centres = Vec::new();
    for axis in 0..3 {
        for sign in [1.0f64, -1.0] {
            for i in 0..K {
                for j in 0..K {
                    centres.push(match axis {
                        0 => v3(sign * L, at(i), at(j)),
                        1 => v3(at(i), sign * L, at(j)),
                        _ => v3(at(i), at(j), sign * L),
                    });
                }
            }
        }
    }
    let r_panel = step * 0.75;
    let np = centres.len();

    // `Tier::Galactic` is the collisionless-gravity regime, and a tier is a
    // physics regime rather than a size. Forty-eight tonnes spread over six
    // metres pull on each other at about 10^-8 m/s^2, which is what deep space
    // is.
    //
    // Continuum is where this box belongs and is unusable until `PLAY.md` §3.3
    // lands: it dispatches to SPH on size alone, and SPH reads a solid through
    // a gas equation of state. Measured, this box's `pressure()` at Continuum is
    // 1.0x10^8 Pa — twice green wood's tensile strength — so it bursts before
    // the ball has moved. See `colliding_pair` above.
    let total = PANEL_MASS * np as f64 + BALL_MASS;
    let spec = SampleSpec::new(np + 1, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let matter = Matter::neutral(total, 2.0 * L, 290.0, Composition::primordial());
    let mut w = World::new(Tree::new(0xB0FFED, matter, Tier::Galactic, spec), 1.0);
    let root = w.tree.root;
    w.tree.refine(root);
    {
        let n = &mut w.tree.nodes[root.get()];
        for (i, c) in centres.iter().enumerate() {
            let b = &mut n.bodies[i];
            b.pos = *c;
            b.vel = Vec3::ZERO;
            b.spin = Vec3::ZERO;
            b.radius = r_panel;
            b.mass = PANEL_MASS;
        }
        let b = &mut n.bodies[np];
        b.pos = Vec3::ZERO;
        b.vel = Vec3::ZERO;
        b.spin = Vec3::ZERO;
        b.radius = R_BALL;
        b.mass = BALL_MASS;
        // Wood, for every body in it. A surface is what a material is for.
        n.topology = Some(Topology { material: Material::GREEN_WOOD, ..Default::default() });
    }
    // A small spec deliberately: the ball needs contents only because
    // `advance_node` will not advance the motion of a node that has none, and
    // `default_spec` at this tier would give it twenty thousand bodies and
    // fifty times the runtime.
    let ball_spec = SampleSpec::new(8, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    let ball = w.tree.promote(root, np, ball_spec);
    w.tree.refine(ball);
    {
        let n = &mut w.tree.nodes[ball.get()];
        n.matter.radius = R_BALL;
        n.matter.mass = BALL_MASS;
        n.motion.offset = Vec3::ZERO;
        // Off every axis, so it reaches several faces and strikes them
        // off-centre. A head-on bounce would never test the friction couple.
        n.motion.velocity = v3(5.0, 1.7, 0.9);
        n.matter.spin = Vec3::ZERO;
        n.topology = Some(Topology { material: Material::GREEN_WOOD, ..Default::default() });
    }
    w.tree.pin(root);
    w.tree.pin(ball);

    // Momentum and angular momentum of everything: the panels as bodies of the
    // node, the ball as the node it was promoted into. The ball's stand-in body
    // is skipped, or it would be counted twice.
    let totals = |w: &World| {
        let (mut p, mut l) = (Vec3::ZERO, Vec3::ZERO);
        let (mut box_p, mut box_l) = (Vec3::ZERO, Vec3::ZERO);
        for (i, b) in w.tree.nodes[root.get()].bodies.iter().enumerate() {
            if i == np {
                continue;
            }
            let q = b.vel.scale(b.mass);
            let m = b.pos.cross(q) + b.spin;
            p += q;
            box_p += q;
            l += m;
            box_l += m;
        }
        let n = &w.tree.nodes[ball.get()];
        let q = n.motion.velocity.scale(n.matter.mass);
        let m = n.motion.offset.cross(q) + n.matter.spin;
        p += q;
        l += m;
        (p, l, box_p, box_l)
    };

    let (p0, l0, _, _) = totals(&w);
    // A fixed denominator, and it has to be fixed. Summing the angular momentum
    // actually present collapses to nothing on the first frame — the ball is at
    // the origin, where `r x p` is zero, and the panels are at rest — so any
    // last-bit noise divided by it reads as a catastrophic leak. The
    // characteristic scale of this assembly is what the ball can carry about
    // the box's own half-width, and it does not move.
    let l_scale = p0.norm() * L;
    let ke0 = 0.5 * BALL_MASS * w.tree.nodes[ball.get()].motion.velocity.norm2();
    let mut worst_p: f64 = 0.0;
    let mut worst_l: f64 = 0.0;
    let mut escaped = false;
    for _ in 0..FRAMES {
        w.advance_node(root, DT);
        w.advance_node(ball, DT);
        let (p, l, _, _) = totals(&w);
        worst_p = worst_p.max((p - p0).norm() / p0.norm());
        worst_l = worst_l.max((l - l0).norm() / l_scale);
        let at = w.tree.nodes[ball.get()].motion.offset;
        escaped |= at.x.abs() > L || at.y.abs() > L || at.z.abs() > L;
    }
    let (_, _, box_p, box_l) = totals(&w);
    let hits = w.stats.contacts_resolved;

    // It has to have actually bounced around, or the rest proves nothing about
    // contact: a ball that never touched anything conserves everything.
    assert!(
        hits >= 5,
        "the ball struck the box {hits} times in {FRAMES} frames, which is too \
         few for this to be measuring collisions at all"
    );
    assert!(!escaped, "the ball left the box, so it stopped being a box");

    // The whole point. Nothing acts on this assembly from outside.
    assert!(
        worst_p < 1e-5,
        "{hits} collisions moved the assembly's momentum by {worst_p} of itself"
    );
    assert!(
        worst_l < 1e-5,
        "{hits} collisions moved the assembly's angular momentum by {worst_l} \
         of what the assembly can carry — the friction couple is unbalanced"
    );

    // And the box responded rather than absorbing the ball silently: it took up
    // the momentum, and the off-centre strikes spun it.
    assert!(
        box_p.norm() > 0.5 * p0.norm(),
        "the box ended with {} of the {} kg m/s the ball arrived with",
        box_p.norm(),
        p0.norm()
    );
    assert!(
        box_l.norm() > 0.0,
        "the box was struck off-centre eight times and never started turning"
    );

    // An off-centre strike turns what it hits. Asserted on the panels rather
    // than on the ball because the ball is always the second of each pair, and
    // a couple taken about one side's centre instead of the shared contact
    // point still conserves — it just puts all of the spin on one side. Zeroing
    // the ball's share leaves this at exactly zero and every conservation
    // assertion above still passing, which is how it would be missed.
    let panel_spin: f64 = w.tree.nodes[root.get()].bodies[..np]
        .iter()
        .map(|b| b.spin.norm())
        .sum();
    assert!(
        panel_spin > 0.0,
        "the box was struck off-centre {hits} times and no panel is turning"
    );

    // Wood does not bounce, and this is the derived restitution doing it rather
    // than a damping term anywhere. The whole assembly's kinetic energy is what
    // is measured, not the ball's: the ball slows partly by handing momentum to
    // a panel, which is not a loss. Measured, over eight impacts: 3.0% of the
    // energy is left. With restitution forced to 1 the same run leaves 23.5%,
    // which is what sets the bound between them.
    //
    // The restitution itself is pinned by `restitution_falls_with_the_speed_of_
    // the_impact` and by `two_promoted_things_collide_and_rebound`, which checks
    // the measured rebound against what the materials say. This is the
    // system-level consequence of it.
    let mut ke = 0.5 * BALL_MASS * w.tree.nodes[ball.get()].motion.velocity.norm2();
    for b in &w.tree.nodes[root.get()].bodies[..np] {
        ke += 0.5 * b.mass * b.vel.norm2();
    }
    assert!(
        ke < ke0 * 0.1,
        "the assembly still has {ke} J of the {ke0} J it started with, after \
         {hits} collisions between two pieces of green wood"
    );
}
