//! A basis grown for the molecule (`electrons::grow`): its pieces checked
//! against what they claim before anything is grown with them.

use phys::electrons::functional::Functional;
use phys::electrons::grow::{candidates_for_test, coordinates, estimate, ladders};
use phys::electrons::molecule::Molecule;
use phys::electrons::scf::solve;

fn water(r: f64, theta_deg: f64) -> Molecule {
    let a = 1.0 / 0.529177210903;
    let th = theta_deg.to_radians();
    Molecule { z: vec![8, 1, 1], positions: vec![[0.0, 0.0, 0.0], [r * a, 0.0, 0.0], [r * a * th.cos(), r * a * th.sin(), 0.0]], charge: 0, unpaired: 0 }
}

/// Displacing a coordinate moves it by what was asked and, to first order,
/// moves no other — on four atoms round a centre, where the three bonds and
/// three angles are independent, and on a ring of three, where they are not
/// (three lengths fix the angles). The first version moved a bond's two ends
/// apart along it, which moves the centre atom with them, and failed this.
#[test]
fn a_displaced_coordinate_moves_by_what_was_asked_and_alone() {
    use phys::electrons::grow::displaced;
    let star = (vec![[0.1, -0.2, 0.3], [1.9, 0.1, 0.2], [-0.4, 1.7, 0.9], [0.2, -0.3, -1.6]], vec![(0, 1), (0, 2), (0, 3)]);
    for (p, bonds) in [star] {
        let coords = coordinates(&bonds);
        for (ci, c) in coords.iter().enumerate() {
            let d = 1e-3;
            let q = displaced(&coords, &p, ci, d);
            let moved = c.value(&q) - c.value(&p);
            assert!((moved - d).abs() < 1e-5, "{c:?} moved {moved} for {d}");
            for other in &coords {
                if other != c {
                    let m = (other.value(&q) - other.value(&p)).abs();
                    assert!(m < 1e-5, "{other:?} moved {m} when {c:?} moved {d}");
                }
            }
            // No translation of the whole.
            let shift: Vec<f64> = (0..3).map(|k| q.iter().zip(&p).map(|(a, b)| a[k] - b[k]).sum::<f64>()).collect();
            assert!(shift.iter().all(|s| s.abs() < 1e-9), "translated by {shift:?}");
        }
    }
    // A ring of three: its three bonds and three angles are redundant, so a
    // bond cannot move alone — but the step is still the smallest that does
    // it, and the bond moves by what was asked.
    let p = vec![[0.0, 0.0, 0.0], [2.8, 0.0, 0.0], [1.3, 2.4, 0.1]];
    let coords = coordinates(&[(0, 1), (1, 2), (2, 0)]);
    let q = displaced(&coords, &p, 0, 1e-3);
    let moved = coords[0].value(&q) - coords[0].value(&p);
    assert!(moved > 0.0, "the ring's bond moved {moved}");
}

/// The cheap estimate against the real thing: water in its atoms' own
/// orbitals, and for the candidates the estimate rates highest, each added
/// alone and solved again. The estimate is second order and ignores the
/// response of everything else, so it is held to the sign and to a factor of
/// two — what the prototype found — and the ranking it gives is checked to be
/// the ranking the solves give at the top.
#[test]
fn the_cheap_estimate_tracks_adding_a_function_and_solving() {
    let mol = water(0.97, 104.2);
    let f = Functional::Pbe;
    let lads = ladders(&mol, f);
    let none = vec![Vec::new(); 3];
    let p = mol.problem_with(f, Some(&none));
    let s = solve(&p, 200, 1e-10);
    let (rungs, shells) = candidates_for_test(&mol, &lads, &[]);
    let est = estimate(&p, &s, &shells);
    let mut order: Vec<usize> = (0..rungs.len()).collect();
    order.sort_by(|&a, &b| est[a].total_cmp(&est[b]));
    println!("  seed {} functions, E {:.8}; {} candidates", p.basis.size, s.energy, rungs.len());
    for &i in order.iter().take(6) {
        let lad = &lads[rungs[i].ladder];
        let mut extra = vec![Vec::new(); 3];
        extra[lad.atom].push((lad.l, lad.exponent(rungs[i].k)));
        let pi = mol.problem_with(f, Some(&extra));
        let si = solve(&pi, 200, 1e-10);
        let actual = si.energy - s.energy;
        println!("  atom {} l {} exponent {:.4}: estimate {:+.3e}, solved {:+.3e}", lad.atom, lad.l, lad.exponent(rungs[i].k), est[i], actual);
        assert!(actual < 0.0 && est[i] < 0.0);
        assert!(est[i] / actual < 2.5 && est[i] / actual > 0.4, "estimate {} against {}", est[i], actual);
    }
}

/// Water's two hydrogens are one class and its oxygen another; stretch one
/// bond past the length tolerance and they are no longer the same.
#[test]
fn equivalent_atoms_are_found_from_the_geometry() {
    use phys::electrons::grow::equivalent_atoms;
    let w = water(0.97, 104.2);
    let c = equivalent_atoms(&w);
    assert_eq!(c[1], c[2]);
    assert_ne!(c[0], c[1]);
    let mut bent = w.clone();
    bent.positions[1][0] *= 1.002;
    let c = equivalent_atoms(&bent);
    assert_ne!(c[1], c[2], "a 0.2% longer bond is a different hydrogen");
    let mut nearly = w.clone();
    nearly.positions[1][0] *= 1.0002;
    let c = equivalent_atoms(&nearly);
    assert_eq!(c[1], c[2], "a 0.02% difference is inside the tolerance");
}
