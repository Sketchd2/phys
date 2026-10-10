//! Water that moves atom by atom: the forces are minus the energy's gradient
//! (through the law, the induced dipoles, the virtual site and the surface), the
//! bonds hold when they are asked to, and a short run conserves its energy.

use phys::flex::{bisector_transfer, intramolecular, Flex, FlexForces, FlexSites, Surface};
use phys::liquid::{Polarisable, SiteSite};
use phys::math::{Quat, Vec3};

const AMU: f64 = 1822.888486;

fn law() -> (SiteSite, Vec<f64>, String) {
    let text = std::fs::read_to_string("law-water-mp2-cloud-final.txt").expect("the cloud law");
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    (law, alpha, text)
}

fn small_liquid() -> Flex {
    let (_, _, text) = law();
    let sites = FlexSites::from_text(&text, [0, 1, 1]);
    let surface = Surface::from_text(&std::fs::read_to_string("surface-water.txt").expect("the surface"));
    let r = 0.9789 / 0.529177210903;
    let th = 104.08f64.to_radians();
    let template = [Vec3::ZERO, Vec3 { x: r, y: 0.0, z: 0.0 }, Vec3 { x: r * th.cos(), y: r * th.sin(), z: 0.0 }];
    let mut pos = Vec::new();
    // 27 molecules 6.6 bohr apart in a box of 19.8: neighbours inside the cut-off, which the
    // first version of this test (eight molecules 11 bohr apart, cut off at 10.5) did not have -
    // its forces were the molecules' own and it tested nothing between them.
    let side = 19.8;
    let mut k = 12345u64;
    let mut rnd = || {
        k = k.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((k >> 33) as f64) / ((1u64 << 31) as f64)
    };
    for ix in 0..3 {
        for iy in 0..3 {
            for iz in 0..3 {
                let centre = Vec3 { x: 3.3 + 6.6 * ix as f64, y: 3.3 + 6.6 * iy as f64, z: 3.3 + 6.6 * iz as f64 };
                let q = Quat { w: rnd() * 2.0 - 1.0, v: Vec3 { x: rnd() * 2.0 - 1.0, y: rnd() * 2.0 - 1.0, z: rnd() * 2.0 - 1.0 } };
                let n = (q.w * q.w + q.v.norm2()).sqrt();
                let q = Quat { w: q.w / n, v: q.v.scale(1.0 / n) };
                let mut m = [Vec3::ZERO; 3];
                for a in 0..3 {
                    m[a] = centre + q.rotate(template[a]) + Vec3 { x: rnd() - 0.5, y: rnd() - 0.5, z: rnd() - 0.5 }.scale(0.1);
                }
                pos.push(m);
            }
        }
    }
    let n = pos.len();
    Flex {
        sites,
        masses: [15.9949146 * AMU, 1.00782503 * AMU, 1.00782503 * AMU],
        cell: [side; 3],
        pos,
        vel: vec![[Vec3::ZERO; 3]; n],
        r_on: 7.6,
        r_cut: 9.6,
        surface,
        bond: None,
        geo: None,
    }
}

fn energy(l: &Flex, law: &SiteSite, alpha: &[f64]) -> f64 {
    let pol = Polarisable { law, alpha: alpha.to_vec() };
    let f = l.forces(&pol, None);
    f.inter + f.intra
}

#[test]
fn the_forces_are_minus_the_gradient_of_the_energy() {
    let (sitesite, alpha, _) = law();
    let mut l = small_liquid();
    let pol = Polarisable { law: &sitesite, alpha: alpha.clone() };
    let f = l.forces(&pol, None);
    let mut worst: f64 = 0.0;
    let mut largest: f64 = 0.0;
    for (i, a, axis) in [(0usize, 0usize, 0usize), (0, 1, 1), (0, 2, 2), (3, 0, 1), (5, 1, 0), (7, 2, 2), (2, 0, 2), (6, 1, 1)] {
        let h = 1e-4;
        let comp = |v: Vec3, k: usize| [v.x, v.y, v.z][k];
        let mut shifted = |s: f64| {
            let old = l.pos[i][a];
            let mut v = old;
            match axis {
                0 => v.x += s,
                1 => v.y += s,
                _ => v.z += s,
            }
            l.pos[i][a] = v;
            let e = energy(&l, &sitesite, &alpha);
            l.pos[i][a] = old;
            e
        };
        let numeric = -(shifted(h) - shifted(-h)) / (2.0 * h);
        let analytic = comp(f.force[i][a], axis);
        println!("  molecule {i} atom {a} axis {axis}: force {analytic:+.8e}, minus the gradient {numeric:+.8e}");
        worst = worst.max((analytic - numeric).abs());
        largest = largest.max(analytic.abs());
    }
    assert!(largest > 1e-4, "no force to speak of: {largest:e}");
    assert!(worst < 2e-6 * largest.max(1.0) + 2e-8, "force error {worst:e} against forces up to {largest:e}");
    // The net force on the whole system is nothing.
    let mut net = Vec3::ZERO;
    for m in &f.force {
        for a in m {
            net += *a;
        }
    }
    assert!(net.norm() < 1e-9, "net force {:e}", net.norm());
}

#[test]
fn the_virtual_site_hands_its_force_to_the_atoms_consistently() {
    let pos = [Vec3 { x: 0.1, y: -0.2, z: 0.0 }, Vec3 { x: 1.8, y: 0.3, z: 0.1 }, Vec3 { x: -0.5, y: 1.7, z: -0.2 }];
    let f = Vec3 { x: 0.3, y: -0.7, z: 0.2 };
    let t = bisector_transfer(&pos, 0.3, f);
    let sum = t[0] + t[1] + t[2];
    assert!((sum - f).norm() < 1e-14, "{sum:?} against {f:?}");
    let site = |p: &[Vec3; 3]| {
        let u = (p[1] - p[0]).unit();
        let v = (p[2] - p[0]).unit();
        p[0] + (u + v).unit().scale(0.3)
    };
    for a in 0..3 {
        for axis in 0..3 {
            let h = 1e-6;
            let bump = |s: f64| {
                let mut p = pos;
                match axis {
                    0 => p[a].x += s,
                    1 => p[a].y += s,
                    _ => p[a].z += s,
                }
                f.dot(site(&p))
            };
            let numeric = (bump(h) - bump(-h)) / (2.0 * h);
            let analytic = [t[a].x, t[a].y, t[a].z][axis];
            assert!((numeric - analytic).abs() < 1e-8, "atom {a} axis {axis}: {numeric} against {analytic}");
        }
    }
}

#[test]
fn the_surface_has_its_minimum_at_zero_and_its_forces_are_its_gradient() {
    let surface = Surface::from_text(&std::fs::read_to_string("surface-water.txt").expect("the surface"));
    let pos = [Vec3::ZERO, Vec3 { x: 1.85, y: 0.0, z: 0.0 }, Vec3 { x: -0.45, y: 1.78, z: 0.0 }];
    let (e, f) = intramolecular(&surface, &pos);
    assert!(e >= -1e-9 && e < 0.05, "{e}");
    for a in 0..3 {
        for axis in 0..3 {
            let h = 1e-5;
            let at = |s: f64| {
                let mut p = pos;
                match axis {
                    0 => p[a].x += s,
                    1 => p[a].y += s,
                    _ => p[a].z += s,
                }
                intramolecular(&surface, &p).0
            };
            let numeric = -(at(h) - at(-h)) / (2.0 * h);
            let analytic = [f[a].x, f[a].y, f[a].z][axis];
            assert!((numeric - analytic).abs() < 1e-8, "atom {a} axis {axis}: {numeric} against {analytic}");
        }
    }
}

#[test]
fn bonds_hold_and_a_short_run_conserves_energy() {
    let (sitesite, alpha, _) = law();
    let pol = Polarisable { law: &sitesite, alpha: alpha.clone() };
    for bond in [None, Some(0.9789 / 0.529177210903)] {
        let mut l = small_liquid();
        l.bond = bond;
        if let Some(len) = bond {
            for m in l.pos.iter_mut() {
                for h in [1usize, 2] {
                    let d = (m[h] - m[0]).unit();
                    m[h] = m[0] + d.scale(len);
                }
            }
        }
        let h = if bond.is_some() { 40.0 } else { 20.0 };
        let mut f = l.forces(&pol, None);
        let total = |l: &Flex, f: &FlexForces| f.inter + f.intra + l.kinetic().1;
        let e0 = total(&l, &f);
        let mut worst_bond: f64 = 0.0;
        let mut spread: f64 = 0.0;
        for _ in 0..150 {
            f = l.step(&pol, f, h, None);
            spread = spread.max((total(&l, &f) - e0).abs());
            if let Some(len) = bond {
                for m in &l.pos {
                    worst_bond = worst_bond.max(((m[1] - m[0]).norm() - len).abs() / len).max(((m[2] - m[0]).norm() - len).abs() / len);
                }
            }
        }
        println!("  bond {bond:?}: energy wandered by {spread:.3e} hartree over 150 steps of {h}, worst bond error {worst_bond:.2e}");
        assert!(spread < 2e-3, "energy not conserved: {spread:e}");
        if bond.is_some() {
            assert!(worst_bond < 1e-9, "a bond moved by {worst_bond:e}");
        }
    }
}

/// E8c pair data: each molecule of a pair is drawn separately from the surface,
/// so the two bonds differ (an asymmetric stretch happens), the draw stays in the
/// fitted region, favours low energy, and leaves the unflexed draw alone.
#[test]
fn drawn_shapes_are_asymmetric_in_region_and_boltzmann() {
    use phys::pairs::FlexDraw;
    use phys::rng::{Purpose, Stream};
    let draw = FlexDraw::from_file("surface-water.txt", 1500.0);
    let mut s = Stream::at(7, 0, 0, Purpose::Positions);
    let n = 4000;
    let (mut asym, mut sum_r, mut sum_th, mut lo_r, mut hi_r) = (0.0f64, 0.0, 0.0, f64::MAX, 0.0f64);
    for _ in 0..n {
        let (r1, r2, th) = draw.shape(&mut s);
        asym += (r1 - r2).abs() / n as f64;
        sum_r += (r1 + r2) / 2.0 / n as f64;
        sum_th += th / n as f64;
        lo_r = lo_r.min(r1.min(r2));
        hi_r = hi_r.max(r1.max(r2));
    }
    println!("mean |r1-r2| {asym:.4} bohr, mean r {sum_r:.4} bohr, mean theta {:.2} deg, r range {lo_r:.3}..{hi_r:.3}", sum_th.to_degrees());
    assert!(asym > 0.03, "the bonds must differ: {asym}");
    assert!(lo_r >= 0.86 / 0.529177210903 && hi_r <= 1.16 / 0.529177210903);
    assert!((sum_r * 0.529177210903 - 0.97).abs() < 0.04 && (sum_th.to_degrees() - 104.0).abs() < 4.0);
    // Colder is narrower: the sampling temperature does what it says.
    let cold = FlexDraw::from_file("surface-water.txt", 300.0);
    let spread = |d: &FlexDraw| {
        let mut s = Stream::at(8, 0, 0, Purpose::Positions);
        (0..2000).map(|_| { let (a, b, _) = d.shape(&mut s); (a - b).abs() }).sum::<f64>() / 2000.0
    };
    assert!(spread(&cold) < 0.6 * spread(&draw));
}

/// A cloud law has sites placed from the atoms, and the biased pair draw must
/// place them: without them the pair is unneutral, its energy is thousands of
/// kcal/mol at any separation, and the draw never accepts anything. It once
/// spun for 37 minutes on pair 0. Also: the same draw with `--flex` terminates
/// and its molecules differ in bond length.
#[test]
fn the_biased_draw_places_a_cloud_laws_sites() {
    use phys::pairs::{Bias, FlexDraw, Monomer};
    if !std::path::Path::new("grow-water.state").exists() {
        return; // the state is not tracked (gitignored): nothing to draw from
    }
    let mut m = Monomer::load("water", false);
    let b = Bias::from_file("law-water-mp2-cloud-slopes.txt", 2.0);
    let seps: Vec<f64> = (0..20).map(|k| m.biased_pair(k, &b).2).collect();
    assert!(seps.iter().all(|s| (4.2..=8.0).contains(s)));
    assert!((seps[0] - 4.60).abs() < 0.01 && (seps[1] - 6.16).abs() < 0.01, "the draw is the one the earlier data was taken with: {seps:?}");
    m.flex = Some(FlexDraw::from_file("surface-water.txt", 1500.0));
    let (a, _, _) = m.biased_pair(3, &b);
    assert!(((a[1] - a[0]).norm() - (a[2] - a[0]).norm()).abs() > 1e-3);
}

/// With the bonds held (`--bonds-at`) both are the held length, only the angle varies.
#[test]
fn held_bonds_draw_only_the_angle() {
    use phys::pairs::FlexDraw;
    use phys::rng::{Purpose, Stream};
    let mut d = FlexDraw::from_file("surface-water.txt", 1500.0);
    let r = 0.97894 / 0.529177210903;
    d.held_bond = Some(r);
    let mut s = Stream::at(11, 0, 0, Purpose::Positions);
    let shapes: Vec<_> = (0..500).map(|_| d.shape(&mut s)).collect();
    assert!(shapes.iter().all(|&(a, b, _)| a == r && b == r));
    let (lo, hi) = shapes.iter().fold((9.0f64, 0.0f64), |(l, h), &(_, _, t)| (l.min(t), h.max(t)));
    assert!(hi.to_degrees() - lo.to_degrees() > 10.0, "the angle must vary: {lo} {hi}");
}

/// A cloud whose parameters are quadratics in the angle, with slopes large enough to matter.
fn geo_cloud(l: &Flex, text: &str, sitesite: &SiteSite) -> phys::geolaw::GeoCloud {
    let bis = SiteSite::bisector_from_text(text).expect("a bisector").1;
    let x0 = phys::geolaw::x_of_law(sitesite, bis);
    let mut g = phys::geolaw::GeoCloud::constant([sitesite.charge[0], sitesite.charge[1]], x0, 104.08f64.to_radians());
    let _ = l;
    // ln widths, bisector distance, populations: each given a slope and a curvature.
    let c1 = [0.15, -0.10, 0.20, 0.08, 0.12, -0.10];
    let c2 = [0.05, 0.04, -0.06, 0.02, -0.05, 0.03];
    for p in 0..6 {
        g.coef[p][1] = c1[p];
        g.coef[p][2] = c2[p];
    }
    g
}

#[test]
fn forces_are_minus_the_gradient_when_the_cloud_follows_the_angle() {
    let (sitesite, alpha, text) = law();
    let mut l = small_liquid();
    // Make the molecules differ in angle, so the parameters differ between them.
    for (i, m) in l.pos.iter_mut().enumerate() {
        let shift = ((i as f64 * 0.7).sin()) * 0.12;
        let (o, h1, h2) = (m[0], m[1], m[2]);
        let u = (h2 - o).unit();
        let v = (h1 - o).unit();
        let n = u.cross(v).unit();
        let rot = Quat { w: (shift / 2.0).cos(), v: n.scale((shift / 2.0).sin()) };
        m[2] = o + rot.rotate(h2 - o);
        let _ = (h1, v);
    }
    l.geo = Some(geo_cloud(&l, &text, &sitesite));
    let pol = Polarisable { law: &sitesite, alpha: alpha.clone() };
    let f = l.forces(&pol, None);
    let mut worst: f64 = 0.0;
    let mut largest: f64 = 0.0;
    let mut lone: f64 = 0.0;
    for (i, a, axis) in [(0usize, 0usize, 0usize), (0, 1, 1), (0, 2, 2), (3, 0, 1), (5, 1, 0), (7, 2, 2), (2, 0, 2), (6, 1, 1), (1, 2, 0), (4, 1, 2)] {
        let h = 1e-4;
        let mut shifted = |s: f64| {
            let old = l.pos[i][a];
            let mut v = old;
            match axis {
                0 => v.x += s,
                1 => v.y += s,
                _ => v.z += s,
            }
            l.pos[i][a] = v;
            let e = energy(&l, &sitesite, &alpha);
            l.pos[i][a] = old;
            e
        };
        let numeric = -(shifted(h) - shifted(-h)) / (2.0 * h);
        let analytic = [f.force[i][a].x, f.force[i][a].y, f.force[i][a].z][axis];
        println!("  geo: molecule {i} atom {a} axis {axis}: force {analytic:+.8e}, minus the gradient {numeric:+.8e}");
        worst = worst.max((analytic - numeric).abs());
        largest = largest.max(analytic.abs());
        lone = lone.max(numeric.abs());
    }
    assert!(largest > 1e-4, "no force to speak of: {largest:e}");
    assert!(worst < 2e-6 * largest.max(1.0) + 2e-8, "force error {worst:e} against forces up to {largest:e}");
    // And the same cloud with the angle's part of the force left out would have been wrong:
    // the slope terms are not small.
    let mut g0 = l.geo.clone().unwrap();
    for p in 0..6 {
        g0.coef[p][1] = 0.0;
        g0.coef[p][2] = 0.0;
    }
    let mut flat = small_liquid();
    flat.pos = l.pos.clone();
    flat.geo = Some(g0);
    let ff = flat.forces(&pol, None);
    println!("  energies with the angle's dependence {:.9} and without {:.9}", f.inter, ff.inter);
    let mut diff: f64 = 0.0;
    for (a, b) in f.force.iter().zip(&ff.force) {
        for k in 0..3 {
            diff = diff.max((a[k] - b[k]).norm());
        }
    }
    println!("  the angle's dependence changes the largest force by {diff:.3e} (forces up to {largest:.3e}; {lone:.3e})");
    assert!(diff > 10.0 * worst, "the test cannot tell the angle's terms from round-off: {diff:e} against {worst:e}");
    let mut net = Vec3::ZERO;
    for m in &f.force {
        for a in m {
            net += *a;
        }
    }
    assert!(net.norm() < 1e-9, "net force {:e}", net.norm());
}
