//! The induced-dipole model (`src/induction.rs`): its error function, its limits,
//! and that its forces are the derivative of its energy.

use phys::induction::{erf, molecular_polarisability, solve, Charge, Cluster, Link, PolSite};
use phys::liquid::switch;
use phys::math::Vec3;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3 { x, y, z }
}

#[test]
fn the_error_function_matches_known_values() {
    for (x, want) in [(0.0, 0.0), (0.5, 0.5204998778130465), (1.0, 0.8427007929497149), (2.0, 0.9953222650189527), (3.0, 0.9999779095030014), (4.5, 0.999999999803384), (-1.0, -0.8427007929497149)] {
        assert!((erf(x) - want).abs() < 1e-14, "erf({x}) = {} against {want}", erf(x));
    }
}

/// Far from a point charge a polarisable site is an ordinary polarisable atom:
/// `E = -1/2 alpha q^2 / r^4`, and a dipole `alpha q / r^2` along the line.
#[test]
fn far_from_a_charge_the_site_is_a_polarisable_atom() {
    let (alpha, q, r) = (5.0, 1.0, 10.0);
    let cl = Cluster { pol: vec![vec![PolSite { pos: v(0.0, 0.0, 0.0), alpha }], vec![]], charges: vec![vec![], vec![Charge { pos: v(r, 0.0, 0.0), q }]] };
    let links = [Link { i: 0, j: 1, shift: Vec3::ZERO, weight: 1.0, dweight: 0.0, d: v(-r, 0.0, 0.0) }];
    let out = solve(&cl, &links, None, 1e-12);
    let want = -0.5 * alpha * q * q / r.powi(4);
    println!("  energy {:.6e} against {want:.6e}; dipole {:?}", out.energy, out.dipoles[0][0]);
    assert!((out.energy / want - 1.0).abs() < 1e-8);
    assert!((out.dipoles[0][0].x.abs() - alpha * q / (r * r)).abs() < 1e-10);
}

/// A molecule's polarisability from its atoms: one site is its own; two sites
/// far apart add; two close ones in line, coupled head to tail along their
/// axis, respond more than the sum along it and less across it.
#[test]
fn a_molecules_polarisability_comes_from_its_atoms_and_their_coupling() {
    let one = molecular_polarisability(&[PolSite { pos: Vec3::ZERO, alpha: 6.0 }]);
    assert!((one[0][0] - 6.0).abs() < 1e-12 && one[0][1].abs() < 1e-12);
    let far = molecular_polarisability(&[PolSite { pos: Vec3::ZERO, alpha: 4.0 }, PolSite { pos: v(0.0, 0.0, 400.0), alpha: 2.0 }]);
    assert!((far[0][0] - 6.0).abs() < 1e-6 && (far[2][2] - 6.0).abs() < 1e-3, "{far:?}");
    let near = molecular_polarisability(&[PolSite { pos: Vec3::ZERO, alpha: 4.0 }, PolSite { pos: v(0.0, 0.0, 3.0), alpha: 4.0 }]);
    println!("  two sites 3 bohr apart, each 4: along {:.4}, across {:.4} (sum 8)", near[2][2], near[0][0]);
    assert!(near[2][2] > 8.0 && near[0][0] < 8.0);
}

/// Three molecules of two polarisable sites and two charges each, close enough
/// that everything interacts. `weights` switches the pair terms on the centres'
/// separation, as the liquid does.
fn build(offsets: &[Vec3; 3], weights: bool) -> (Cluster, Vec<Link>) {
    let base = [v(0.0, 0.0, 0.0), v(7.0, 1.0, -1.0), v(2.0, 8.0, 2.5)];
    let shapes = [[v(0.9, 0.2, 0.1), v(-0.8, 0.5, -0.2)], [v(-0.5, 0.9, 0.3), v(0.7, -0.4, 0.6)], [v(0.2, -0.8, 0.5), v(-0.3, 0.4, -0.9)]];
    let mut cl = Cluster::default();
    let mut coms = Vec::new();
    for m in 0..3 {
        let c = base[m] + offsets[m];
        coms.push(c);
        cl.pol.push(vec![PolSite { pos: c + shapes[m][0], alpha: 5.0 + m as f64 }, PolSite { pos: c + shapes[m][1], alpha: 2.0 }]);
        cl.charges.push(vec![Charge { pos: c + shapes[m][0].scale(1.1), q: 0.6 - 0.2 * m as f64 }, Charge { pos: c + shapes[m][1].scale(0.9), q: -0.5 + 0.1 * m as f64 }]);
    }
    let mut links = Vec::new();
    for i in 0..3 {
        for j in i + 1..3 {
            let d = coms[i] - coms[j];
            let (w, dw) = if weights { switch(d.norm(), 6.0, 14.0) } else { (1.0, 0.0) };
            links.push(Link { i, j, shift: Vec3::ZERO, weight: w, dweight: dw, d });
        }
    }
    (cl, links)
}

fn unit(axis: usize, s: f64) -> Vec3 {
    [Vec3 { x: s, y: 0.0, z: 0.0 }, Vec3 { x: 0.0, y: s, z: 0.0 }, Vec3 { x: 0.0, y: 0.0, z: s }][axis]
}

/// With the weights fixed, each site's force is minus the derivative of the
/// energy (the dipoles minimised again each time) by central differences; with
/// the weights a function of the separation, moving a whole molecule rigidly
/// changes the energy by the forces on its sites and the weights' terms.
#[test]
fn forces_are_the_derivative_of_the_energy() {
    let zero = [Vec3::ZERO; 3];
    let (cl, links) = build(&zero, false);
    let base = solve(&cl, &links, None, 1e-13);
    let h = 1e-5;
    let mut worst = 0.0f64;
    for (m, which, is_pol) in [(0usize, 0usize, true), (1, 1, true), (2, 0, false), (1, 1, false)] {
        for axis in 0..3 {
            let shifted = |s: f64| {
                let mut c = cl.clone();
                if is_pol {
                    c.pol[m][which].pos += unit(axis, s);
                } else {
                    c.charges[m][which].pos += unit(axis, s);
                }
                solve(&c, &links, None, 1e-13).energy
            };
            let numeric = -(shifted(h) - shifted(-h)) / (2.0 * h);
            let f = if is_pol { base.force_pol[m][which] } else { base.force_charge[m][which] };
            let analytic = [f.x, f.y, f.z][axis];
            worst = worst.max((numeric - analytic).abs() / (1e-9 + analytic.abs().max(numeric.abs())));
            assert!((numeric - analytic).abs() < 1e-8 * (1.0 + analytic.abs()), "molecule {m} {} {which} axis {axis}: numeric {numeric:.9e} analytic {analytic:.9e}", if is_pol { "dipole site" } else { "charge" });
        }
    }
    println!("  forces on single sites agree with differences of the energy; worst relative {worst:.1e}");
    let (cl2, links2) = build(&zero, true);
    let base2 = solve(&cl2, &links2, None, 1e-13);
    for m in 0..3 {
        for axis in 0..3 {
            let moved = |s: f64| {
                let mut o = zero;
                o[m] = unit(axis, s);
                let (c, l) = build(&o, true);
                solve(&c, &l, None, 1e-13).energy
            };
            let numeric = -(moved(h) - moved(-h)) / (2.0 * h);
            let mut f = Vec3::ZERO;
            for x in &base2.force_pol[m] {
                f += *x;
            }
            for x in &base2.force_charge[m] {
                f += *x;
            }
            for (k, l) in links2.iter().enumerate() {
                if l.i == m {
                    f += base2.link_force[k];
                }
                if l.j == m {
                    f -= base2.link_force[k];
                }
            }
            let analytic = [f.x, f.y, f.z][axis];
            assert!((numeric - analytic).abs() < 1e-8 * (1.0 + analytic.abs()), "molecule {m} axis {axis}: numeric {numeric:.9e} analytic {analytic:.9e}");
        }
    }
    println!("  the weights' terms close the rigid-body force too");
}

/// A warm start reaches the same dipoles in fewer iterations.
#[test]
fn a_warm_start_is_the_same_answer_sooner() {
    let zero = [Vec3::ZERO; 3];
    let (cl, links) = build(&zero, true);
    let cold = solve(&cl, &links, None, 1e-10);
    let moved = [v(0.01, 0.0, 0.0), Vec3::ZERO, v(0.0, 0.01, 0.0)];
    let (cl2, links2) = build(&moved, true);
    let from_cold = solve(&cl2, &links2, None, 1e-10);
    let from_warm = solve(&cl2, &links2, Some(&cold.dipoles), 1e-10);
    println!("  iterations from zero {}, from the last step's dipoles {}", from_cold.iterations, from_warm.iterations);
    assert!(from_warm.iterations < from_cold.iterations);
    assert!((from_warm.energy - from_cold.energy).abs() < 1e-12);
}
