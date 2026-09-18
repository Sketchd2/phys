//! Convex collision geometry: the hull bake and the narrow phase.
//!
//! Every case here has an answer that can be written down without running the
//! engine — a point-to-segment distance, a face-to-face separation, a solid
//! sphere's inertia — because a narrow phase checked only against itself is a
//! narrow phase that agrees with its own mistakes.

use phys::math::{v3, Vec3};
use phys::shape::{closest, Hull, Sphere};

const TOL: f64 = 1e-12;

fn near(a: f64, b: f64, tol: f64, what: &str) {
    assert!(
        (a - b).abs() <= tol,
        "{what}: expected {b}, measured {a}, off by {}",
        (a - b).abs()
    );
}

// ---------------------------------------------------------------------------
// The bake
// ---------------------------------------------------------------------------

#[test]
fn a_hull_of_one_sphere_is_that_sphere() {
    let h = Hull::sphere(v3(1.0, 2.0, 3.0), 0.5);
    assert_eq!(h.len(), 1);
    near(h.bound(), 0.5, TOL, "bound of a lone sphere");
    assert_eq!(h.centre(), v3(1.0, 2.0, 3.0));
    // Support in any direction is the centre plus the radius along it.
    let d = v3(1.0, 1.0, 1.0);
    let want = v3(1.0, 2.0, 3.0) + d.unit().scale(0.5);
    assert!((h.support(d) - want).norm() <= TOL, "support of a lone sphere");
}

#[test]
fn the_bake_drops_a_sphere_swallowed_by_another() {
    // A pea inside a beachball contributes nothing to the surface and must not
    // cost a comparison on every support query.
    let h = Hull::of_spheres([
        Sphere::new(Vec3::ZERO, 10.0),
        Sphere::new(v3(1.0, 0.0, 0.0), 0.1),
        Sphere::new(v3(0.0, 2.0, 0.0), 0.5),
    ]);
    assert_eq!(h.len(), 1, "only the beachball reaches the surface");
    near(h.bound(), 10.0, TOL, "bound after culling");
}

#[test]
fn the_bake_keeps_two_identical_spheres_once() {
    // Both are contained in the other, so a naive test drops both and leaves an
    // empty hull — a thing that has stopped colliding with anything.
    let h = Hull::of_spheres([
        Sphere::new(v3(1.0, 1.0, 1.0), 2.0),
        Sphere::new(v3(1.0, 1.0, 1.0), 2.0),
    ]);
    assert_eq!(h.len(), 1, "duplicates must leave exactly one behind");
    assert!(!h.is_empty());
}

#[test]
fn a_hull_of_two_spheres_is_exactly_a_capsule() {
    let h = Hull::capsule(v3(-2.0, 0.0, 0.0), v3(2.0, 0.0, 0.0), 0.5);
    assert_eq!(h.len(), 2);
    // Along the axis the surface is at the end cap.
    near(h.support(v3(1.0, 0.0, 0.0)).x, 2.5, TOL, "capsule end cap");
    // Across it, at the cylinder wall, and at the same height anywhere along.
    near(h.support(v3(0.0, 1.0, 0.0)).y, 0.5, TOL, "capsule wall");
}

// ---------------------------------------------------------------------------
// The narrow phase, against distances that can be written down
// ---------------------------------------------------------------------------

#[test]
fn two_spheres_reproduce_the_sphere_path_exactly() {
    // The result the old code gave, which the hull path must not change for the
    // shape that already worked.
    let (ca, ra) = (v3(0.0, 0.0, 0.0), 1.0);
    let (cb, rb) = (v3(3.0, 4.0, 0.0), 0.5);
    let c = closest(&Hull::sphere(ca, ra), &Hull::sphere(cb, rb)).expect("two spheres");
    near(c.gap, 5.0 - 1.0 - 0.5, TOL, "sphere-sphere gap");
    let n = (cb - ca).unit();
    assert!((c.normal - n).norm() <= TOL, "sphere-sphere normal is the line of centres");
    assert!((c.on_a - (ca + n.scale(ra))).norm() <= TOL, "witness on a");
    assert!((c.on_b - (cb - n.scale(rb))).norm() <= TOL, "witness on b");
    assert!(!c.degenerate);
}

#[test]
fn a_sphere_against_a_capsule_is_the_point_to_segment_distance() {
    // The capsule runs along x; the sphere sits above its middle, where the
    // closest feature is the *shaft* and not either end. A sphere-only engine
    // measures to the midpoint bead instead, which is the defect being caught.
    let cap = Hull::capsule(v3(-5.0, 0.0, 0.0), v3(5.0, 0.0, 0.0), 0.25);
    let ball = Hull::sphere(v3(1.3, 2.0, 0.0), 0.4);
    let c = closest(&cap, &ball).expect("capsule and sphere");
    near(c.gap, 2.0 - 0.25 - 0.4, TOL, "point-to-shaft gap");
    assert!((c.normal - v3(0.0, 1.0, 0.0)).norm() <= 1e-9, "normal is perpendicular to the shaft");

    // Past the end, the closest feature becomes the cap and the distance is to
    // the end sphere.
    let beyond = Hull::sphere(v3(9.0, 0.0, 0.0), 0.4);
    let c2 = closest(&cap, &beyond).expect("capsule and sphere beyond the end");
    near(c2.gap, 4.0 - 0.25 - 0.4, TOL, "point-to-cap gap");
}

#[test]
fn two_parallel_capsules_measure_across_the_shafts() {
    let a = Hull::capsule(v3(-3.0, 0.0, 0.0), v3(3.0, 0.0, 0.0), 0.3);
    let b = Hull::capsule(v3(-3.0, 1.5, 0.0), v3(3.0, 1.5, 0.0), 0.2);
    let c = closest(&a, &b).expect("two capsules");
    near(c.gap, 1.5 - 0.3 - 0.2, TOL, "parallel capsule gap");
}

#[test]
fn a_slab_of_panels_presents_a_face_and_not_a_row_of_beads() {
    // This is the ball-in-box defect in miniature, and it is the reason the
    // hull exists. Sixteen panels in a 4x4 grid on the x = 3 plane make one
    // wall. A ball approaching along x must stop at the wall's *face*,
    // wherever across the wall it happens to be aimed — not at whichever bead
    // is nearest, which is what a sphere-only path measures.
    let r = 1.125;
    let at = |i: usize| -2.25 + 1.5 * i as f64;
    let mut panels = Vec::new();
    for i in 0..4 {
        for j in 0..4 {
            panels.push(Sphere::new(v3(3.0, at(i), at(j)), r));
        }
    }
    let wall = Hull::of_spheres(panels.clone());

    // Aimed dead at a panel centre, and aimed at the seam between four of
    // them. Both must give the same face distance.
    let face = 3.0 - r; // the wall's inner surface
    for (label, y, z) in [("at a panel centre", at(1), at(2)), ("at a seam", 0.0, 0.0)] {
        let ball = Hull::sphere(v3(0.0, y, z), 0.4);
        let c = closest(&ball, &wall).expect("ball and wall");
        near(c.gap, face - 0.4, 1e-9, &format!("gap {label}"));
        assert!(
            (c.normal - v3(1.0, 0.0, 0.0)).norm() <= 1e-9,
            "normal {label} should be the wall's own, measured {:?}",
            c.normal
        );
    }

    // And the defect it replaces: against the nearest single panel alone, the
    // seam shot reads a materially larger gap, because it measures to a bead.
    let nearest_bead = panels
        .iter()
        .map(|p| (p.centre - v3(0.0, 0.0, 0.0)).norm() - p.radius - 0.4)
        .fold(f64::INFINITY, f64::min);
    assert!(
        nearest_bead > face - 0.4 + 0.05,
        "the bead path should read a gap at least 5 cm larger at the seam, or \
         this test is not measuring the difference it claims: bead {nearest_bead}, \
         face {}",
        face - 0.4
    );
}

#[test]
fn the_narrow_phase_is_symmetric_to_the_bit() {
    // `contact` applies one impulse with both signs, so a normal that differed
    // by a rounding step depending on which side asked would put a different
    // force on each. Both sides must agree exactly.
    let a = Hull::capsule(v3(-1.0, 0.3, 0.2), v3(2.0, -0.4, 0.9), 0.35);
    let b = Hull::of_spheres([
        Sphere::new(v3(1.0, 3.0, 0.0), 0.5),
        Sphere::new(v3(2.2, 3.4, 0.7), 0.25),
        Sphere::new(v3(0.4, 4.1, -0.3), 0.6),
    ]);
    let ab = closest(&a, &b).expect("ab");
    let ba = closest(&b, &a).expect("ba");
    assert_eq!(ab.gap, ba.gap, "gap must be identical from either side");
    assert_eq!(ab.normal, ba.normal.scale(-1.0), "normals must be exact opposites");
    assert_eq!(ab.on_a, ba.on_b, "witness points must swap exactly");
    assert_eq!(ab.on_b, ba.on_a, "witness points must swap exactly");
}

#[test]
fn a_hull_translated_measures_the_same_shape_somewhere_else() {
    // Expressing a hull in another node's frame is a translation, and a narrow
    // phase whose answer drifted with position would break the moment two
    // nodes far from the origin touched.
    let a = Hull::capsule(v3(0.0, 0.0, 0.0), v3(1.0, 0.0, 0.0), 0.2);
    let b = Hull::sphere(v3(0.5, 2.0, 0.0), 0.3);
    let here = closest(&a, &b).expect("here");

    let by = v3(1e6, -3e5, 7e4);
    let there = closest(&a.translated(by), &b.translated(by)).expect("there");
    near(there.gap, here.gap, 1e-9, "gap after translation");
    assert!((there.normal - here.normal).norm() <= 1e-9, "normal after translation");
}

#[test]
fn an_overlap_is_reported_as_a_negative_gap_with_a_measured_normal() {
    // Surfaces interpenetrating is the ordinary contact case and must not be
    // degenerate — the cores are still a clear radius apart.
    let a = Hull::sphere(v3(0.0, 0.0, 0.0), 1.0);
    let b = Hull::sphere(v3(1.5, 0.0, 0.0), 1.0);
    let c = closest(&a, &b).expect("overlapping spheres");
    near(c.gap, -0.5, TOL, "overlap depth");
    assert!(!c.degenerate, "a surface overlap is not a degenerate query");
    assert!((c.normal - v3(1.0, 0.0, 0.0)).norm() <= TOL);
}

#[test]
fn cores_inside_each_other_are_flagged_rather_than_guessed_at() {
    // A ball whose centre is inside the wall's own core. There is no measured
    // normal left, and reporting one anyway is how a contact pushes the wrong
    // way.
    let wall = Hull::of_spheres([
        Sphere::new(v3(-1.0, 0.0, 0.0), 0.5),
        Sphere::new(v3(1.0, 0.0, 0.0), 0.5),
        Sphere::new(v3(0.0, 1.0, 0.0), 0.5),
    ]);
    let ball = Hull::sphere(v3(0.0, 0.3, 0.0), 0.1);
    let c = closest(&wall, &ball).expect("ball inside the core");
    assert!(c.degenerate, "a core interpenetration must say so");
}

#[test]
fn an_empty_hull_collides_with_nothing() {
    let a = Hull::empty();
    let b = Hull::sphere(Vec3::ZERO, 1.0);
    assert!(closest(&a, &b).is_none());
    assert!(closest(&b, &a).is_none());
}

// ---------------------------------------------------------------------------
// Mass properties, against the textbook
// ---------------------------------------------------------------------------

#[test]
fn a_lone_sphere_has_the_solid_spheres_inertia() {
    let h = Hull::sphere(v3(4.0, 5.0, 6.0), 2.0);
    let mp = h.mass_properties(&[3.0]).expect("mass properties");
    near(mp.mass, 3.0, TOL, "mass");
    assert!((mp.centre_of_mass - v3(4.0, 5.0, 6.0)).norm() <= TOL, "centre of mass");
    let want = 0.4 * 3.0 * 4.0; // (2/5) m r^2
    for k in 0..3 {
        near(mp.inertia.0[k][k], want, 1e-9, "diagonal inertia");
    }
    near(mp.inertia.0[0][1], 0.0, 1e-12, "off-diagonal inertia");
}

#[test]
fn two_masses_on_an_axis_obey_the_parallel_axis_theorem() {
    // Point-like spheres a distance d apart: about the centre of mass the
    // inertia across the axis is `2 m (d/2)^2` and along it is nothing.
    let d = 4.0;
    let m = 2.5;
    let h = Hull::of_spheres([
        Sphere::new(v3(-d / 2.0, 0.0, 0.0), 0.0),
        Sphere::new(v3(d / 2.0, 0.0, 0.0), 0.0),
    ]);
    let mp = h.mass_properties(&[m, m]).expect("mass properties");
    near(mp.mass, 2.0 * m, TOL, "total mass");
    assert!(mp.centre_of_mass.norm() <= TOL, "centre of mass is the midpoint");
    near(mp.inertia.0[0][0], 0.0, 1e-12, "no inertia about the axis through both");
    let across = 2.0 * m * (d / 2.0) * (d / 2.0);
    near(mp.inertia.0[1][1], across, 1e-9, "inertia across the axis");
    near(mp.inertia.0[2][2], across, 1e-9, "inertia across the axis");
}

#[test]
fn mass_properties_refuse_a_list_that_does_not_match() {
    let h = Hull::capsule(Vec3::ZERO, v3(1.0, 0.0, 0.0), 0.1);
    assert!(h.mass_properties(&[1.0]).is_none(), "two spheres, one mass");
    assert!(h.mass_properties(&[0.0, 0.0]).is_none(), "no mass at all");
}

/// Out of reach, the query answers from the bounds and says less than it knows.
///
/// `PLAY.md` §7 item 8. `closest_of` skips any piece that cannot be the nearest
/// pair, which is exact — a hull's bound contains it — and when *nothing* can
/// touch it answers from the two bounding spheres without descending at all.
///
/// The gap it then returns is a **lower bound** on the true separation rather
/// than the separation itself, which is worth pinning: every caller in the
/// engine uses it to decide whether a pair is in contact, and a lower bound is
/// sound for that. A caller that wanted the true distance between two far-apart
/// structures would need the descent, and would have to ask for it.
#[test]
fn pieces_out_of_reach_are_not_tested_and_the_gap_is_a_lower_bound() {
    use phys::math::v3;
    use phys::shape::{closest, closest_of, Hull};

    // A long capsule and a small sphere well off the end of it. The capsule's
    // bounding sphere is much fatter than the capsule, so the two answers
    // differ by a knowable amount.
    let capsule = Hull::capsule(v3(-5.0, 0.0, 0.0), v3(5.0, 0.0, 0.0), 0.1);
    let ball = Hull::sphere(v3(0.0, 40.0, 0.0), 0.5);

    let exact = closest(&capsule, &ball).expect("two non-empty hulls");
    let bounded = closest_of(std::slice::from_ref(&capsule), std::slice::from_ref(&ball))
        .expect("two non-empty sides");
    println!(
        "  capsule to ball: exact gap {:.4} m, from the bounds {:.4} m",
        exact.gap, bounded.gap
    );
    assert!(exact.gap > 0.0, "they are not touching, or this measures nothing");
    assert!(
        bounded.gap > 0.0,
        "a lower bound on a positive gap must still be positive, or contact \
         would be reported where there is none"
    );
    assert!(
        bounded.gap <= exact.gap + 1e-9,
        "the bound must not over-state the separation: {} against {}",
        bounded.gap,
        exact.gap
    );

    // And within reach it is the exact answer again, piece by piece.
    let near = Hull::sphere(v3(0.0, 0.7, 0.0), 0.5);
    let exact = closest(&capsule, &near).expect("two non-empty hulls");
    let through = closest_of(std::slice::from_ref(&capsule), std::slice::from_ref(&near))
        .expect("two non-empty sides");
    assert!(
        (through.gap - exact.gap).abs() < 1e-12,
        "within reach the answer must be the descent's: {} against {}",
        through.gap,
        exact.gap
    );

    // Many pieces, one of which is the answer: the rest are skipped and the
    // answer is unchanged.
    let many: Vec<Hull> = (0..64)
        .map(|i| Hull::capsule(v3(i as f64 * 3.0, 0.0, 0.0), v3(i as f64 * 3.0, 1.0, 0.0), 0.1))
        .collect();
    let one = [Hull::sphere(v3(0.0, 0.5, 0.8), 0.2)];
    let best = closest_of(&many, &one).expect("non-empty");
    let brute = many
        .iter()
        .filter_map(|h| closest(h, &one[0]))
        .map(|c| c.gap)
        .fold(f64::MAX, f64::min);
    println!("  64 pieces: skipping gave {:.6} m, testing all gave {:.6} m", best.gap, brute);
    assert!(
        (best.gap - brute).abs() < 1e-12,
        "skipping a piece that cannot win must not change the answer"
    );
}
