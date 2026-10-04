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

/// The estimate does not depend on which way the molecule faces: water as it
/// is and turned 30 degrees, for a compact and a diffuse candidate on its
/// oxygen. The diffuse one moved by 0.5% when the estimate cut its candidates'
/// remainders at 1% of their norm.
#[test]
fn the_estimate_does_not_depend_on_which_way_the_molecule_faces() {
    use phys::electrons::basis::Shell;
    let base = water(0.97724, 104.354).positions;
    let f = Functional::Pbe;
    let extra = vec![vec![(3usize, 0.54896)], vec![(0usize, 0.32365), (0, 2.09499)], vec![(0usize, 0.32365), (0, 2.09499)]];
    let mut out = Vec::new();
    for phi in [0.0f64, 30f64.to_radians()] {
        let (c, s) = (phi.cos(), phi.sin());
        let p: Vec<[f64; 3]> = base.iter().map(|v| [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]]).collect();
        let m = Molecule { z: vec![8, 1, 1], positions: p.clone(), charge: 0, unpaired: 0 };
        let prob = m.problem_with(f, Some(&extra));
        let sol = solve(&prob, 200, 1e-10);
        out.push(estimate(&prob, &sol, &[Shell::primitive(p[0], 2, 0.1), Shell::primitive(p[0], 2, 0.3)]));
    }
    for k in 0..2 {
        let rel = (out[0][k] / out[1][k] - 1.0).abs();
        println!("  candidate {k}: {:.10e} against {:.10e}, {rel:.1e}", out[0][k], out[1][k]);
        assert!(rel < 1e-4, "candidate {k} changed by {rel:e} when turned");
    }
}

/// What a run writes after every round reads back as what was written, to
/// the bit — otherwise a run carried on after a restart is a different run.
#[test]
fn a_saved_round_reads_back_exactly() {
    use phys::electrons::grow::{Ladder, Resume, Rung};
    let r = Resume {
        ladders: vec![Ladder { atom: 1, l: 2, first: 0.1 / 3.0, ratio: 2.570_000_000_000_1, lo: -2, hi: 7 }],
        chosen: vec![Rung { ladder: 0, k: -1 }, Rung { ladder: 0, k: 6 }],
        positions: vec![[0.1, -1.0 / 7.0, 1e-300], [2.0f64.sqrt(), 0.0, -0.0]],
        history: vec![vec![1.8314, 104.5f64.to_radians()], vec![1.0 / 3.0, 2.0]],
    };
    let back = Resume::from_text(&r.to_text()).expect("it parses");
    assert_eq!(back, r);
    assert!(Resume::from_text("ladder 1 2 three").is_none(), "a broken line is refused, not skipped");
}

/// A symmetry carrying one coordinate to another: each of methane's C-H bonds
/// onto every other with every distance kept, and nothing carrying a bond onto
/// an angle or water's O-H onto a different kind of bond.
#[test]
fn a_symmetry_carries_equivalent_coordinates_onto_each_other() {
    use phys::electrons::grow::{symmetry_taking, Coordinate};
    let a = 1.0 / 0.529177210903;
    let (r, t) = (1.0953 * a, 1.0 / 3f64.sqrt());
    let methane = Molecule { z: vec![6, 1, 1, 1, 1], positions: vec![[0.0; 3], [r * t, r * t, r * t], [-r * t, -r * t, r * t], [-r * t, r * t, -r * t], [r * t, -r * t, -r * t]], charge: 0, unpaired: 0 };
    let dist = |m: &Molecule, i: usize, j: usize| (0..3).map(|k| (m.positions[i][k] - m.positions[j][k]).powi(2)).sum::<f64>().sqrt();
    for to in 2..5 {
        let perm = symmetry_taking(&methane, &Coordinate::Bond(0, 1), &Coordinate::Bond(0, to)).expect("C-H onto C-H");
        assert_eq!((perm[0], perm[1]), (0, to));
        for i in 0..5 {
            for j in 0..5 {
                assert!((dist(&methane, i, j) - dist(&methane, perm[i], perm[j])).abs() < 1e-9, "the permutation keeps distances");
            }
        }
    }
    let w = water(0.97, 104.2);
    assert!(symmetry_taking(&w, &Coordinate::Bond(0, 1), &Coordinate::Bond(0, 2)).is_some(), "water's two O-H bonds are equivalent");
    assert!(symmetry_taking(&w, &Coordinate::Bond(0, 1), &Coordinate::Angle(1, 0, 2)).is_none(), "a bond is not an angle");
    let bent = water(0.97, 104.2);
    let mut lopsided = bent.clone();
    lopsided.positions[2][0] *= 1.05;
    assert!(symmetry_taking(&lopsided, &Coordinate::Bond(0, 1), &Coordinate::Bond(0, 2)).is_none(), "bonds of different lengths are not equivalent");
}

/// Probing one coordinate per symmetry class against probing every
/// coordinate, on water at its seed basis: the predictions agree to the
/// molecule's own residual asymmetry and the same functions are picked.
/// Measured on methane before it was adopted: 2.3e-3 of a tolerance at
/// worst, against predictions up to 8.2, the first eight groups picked the
/// same, and 208 s against 40.
#[test]
fn probing_by_symmetry_picks_what_probing_everything_picks() {
    use phys::electrons::grow::{equivalent_atoms, extras, predict_moves};
    let mol = water(0.97, 104.2);
    let bonds = vec![(0, 1), (0, 2)];
    let coords = coordinates(&bonds);
    let f = Functional::Pbe;
    let lads = ladders(&mol, f);
    let extra = extras(&mol, &lads, &[]);
    let p0 = mol.problem_with(f, Some(&extra));
    let s0 = solve(&p0, 200, 1e-10);
    let values: Vec<f64> = coords.iter().map(|c| c.value(&mol.positions)).collect();
    let (rungs, shells) = candidates_for_test(&mol, &lads, &[]);
    let full = predict_moves(&mol, &coords, &values, f, &extra, &s0, &lads, &rungs, &shells, false, &|_| {});
    let sym = predict_moves(&mol, &coords, &values, f, &extra, &s0, &lads, &rungs, &shells, true, &|_| {});
    let tols: Vec<f64> = coords.iter().map(|c| c.tolerance(&mol.positions)).collect();
    let mut worst: f64 = 0.0;
    for ri in 0..rungs.len() {
        for ci in 0..coords.len() {
            worst = worst.max((full[ri][ci] - sym[ri][ci]).abs() / tols[ci]);
        }
    }
    let score = |row: &Vec<f64>| row.iter().zip(&tols).map(|(p, t)| p.abs() / t).fold(0.0, f64::max);
    let class = equivalent_atoms(&mol);
    let groups = |pred: &Vec<Vec<f64>>| {
        let mut o: Vec<usize> = (0..rungs.len()).collect();
        o.sort_by(|&a, &b| score(&pred[b]).total_cmp(&score(&pred[a])).then(rungs[a].cmp(&rungs[b])));
        let mut g: Vec<(usize, usize, i32)> = Vec::new();
        for i in o {
            let l = &lads[rungs[i].ladder];
            let key = (class[l.atom], l.l, rungs[i].k);
            if !g.contains(&key) {
                g.push(key);
            }
        }
        g
    };
    let (gf, gs) = (groups(&full), groups(&sym));
    println!("  worst difference {worst:.2e} tolerances; first groups {:?}", &gf[..6]);
    assert!(worst < 1e-2, "probing by symmetry moved a prediction by {worst:.2e} tolerances");
    assert_eq!(gf[..6], gs[..6], "probing by symmetry picked differently");
}
