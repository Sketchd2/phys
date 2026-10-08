//! Boys localisation: the density does not change, the orbitals become what a
//! chemist draws (two bonds and two lone pairs for water), and the centroids
//! are the molecule's own.

use phys::electrons::functional::Functional;
use phys::electrons::grow::Resume;
use phys::electrons::hf::hartree_fock;
use phys::electrons::integrals::dipole;
use phys::electrons::localise::boys;
use phys::electrons::molecule::Molecule;

#[test]
fn water_localises_to_two_bonds_and_two_lone_pairs() {
    let state = Resume::from_text(&std::fs::read_to_string("grow-water.state").expect("a growth state")).expect("readable");
    let mol = Molecule { z: vec![8, 1, 1], positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let p = mol.problem(Functional::Pbe);
    let hf = hartree_fock(&p, 100, 1e-10);
    let n = p.basis.size;
    let dip = dipole(&p.basis, [0.0; 3]);
    // Valence orbitals: all but the oxygen 1s.
    let loc = boys(&hf.c, hf.m, n, 1..hf.occupied, &dip);
    // The density of the localised set is the canonical set's.
    let density = |c: &[f64]| -> Vec<f64> {
        let mut d = vec![0.0; n * n];
        for a in 0..n {
            for b in 0..n {
                d[a * n + b] = (1..hf.occupied).map(|i| c[a * hf.m + i] * c[b * hf.m + i]).sum();
            }
        }
        d
    };
    let (before, after) = (density(&hf.c), density(&loc.c));
    let worst = before.iter().zip(&after).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
    assert!(worst < 1e-10, "the density moved by {worst:.2e}");
    // The functional increased.
    let spread = |cs: &[[f64; 3]]| -> f64 { cs.iter().map(|c| c.iter().map(|x| x * x).sum::<f64>()).sum() };
    let canonical: Vec<[f64; 3]> = (1..hf.occupied)
        .map(|i| {
            let mut c = [0.0; 3];
            for k in 0..3 {
                for a in 0..n {
                    for b in 0..n {
                        c[k] += hf.c[a * hf.m + i] * hf.c[b * hf.m + i] * dip[k].get(a, b);
                    }
                }
            }
            c
        })
        .collect();
    assert!(spread(&loc.centroids) > spread(&canonical) + 0.1, "{} against {}", spread(&loc.centroids), spread(&canonical));
    // Where they are: distances of each centroid from the oxygen, and from the nearer hydrogen.
    let o = state.positions[0];
    let mut out = Vec::new();
    for c in &loc.centroids {
        let to_o = ((c[0] - o[0]).powi(2) + (c[1] - o[1]).powi(2) + (c[2] - o[2]).powi(2)).sqrt();
        let to_h = (1..3).map(|h| ((c[0] - state.positions[h][0]).powi(2) + (c[1] - state.positions[h][1]).powi(2) + (c[2] - state.positions[h][2]).powi(2)).sqrt()).fold(f64::INFINITY, f64::min);
        out.push((to_o, to_h));
        println!("  centroid {:?}: {:.3} bohr from O, {:.3} from the nearest H", c, to_o, to_h);
    }
    out.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    // Two bonds (near a hydrogen, about half the bond length from it) and two
    // lone pairs (far from both, close to the oxygen).
    assert!(out[0].1 < 1.6 && out[1].1 < 1.6, "{out:?}");
    assert!(out[2].1 > 2.0 && out[3].1 > 2.0 && out[2].0 < 0.9 && out[3].0 < 0.9, "{out:?}");
    assert!((out[2].0 - out[3].0).abs() < 1e-6 && (out[0].1 - out[1].1).abs() < 1e-3, "the two of each kind are equivalent: {out:?}");
}
