//! A molecule's density shared among its atoms (`electrons::partition`, E6).

use phys::electrons::functional::Functional;
use phys::electrons::molecule::Molecule;
use phys::electrons::partition::{atom_pair_c6, partition, FreeAtom};
use phys::electrons::scf::solve;
use phys::electrons::vdw::{c6, Z_AB_DF1};

fn water() -> Molecule {
    let a = 1.0 / 0.529177210903;
    let th = 104.52f64.to_radians() / 2.0;
    Molecule { z: vec![8, 1, 1], positions: vec![[0.0; 3], [0.9572 * a * th.sin(), 0.0, 0.9572 * a * th.cos()], [-0.9572 * a * th.sin(), 0.0, 0.9572 * a * th.cos()]], charge: 0, unpaired: 0 }
}

/// A free atom's tabulated density holds its electrons.
#[test]
fn a_free_atom_holds_its_electrons() {
    for z in [1u32, 8] {
        let fa = FreeAtom::new(z, Functional::Pbe);
        // 4 pi int r^2 rho dr, in ln r: 4 pi int r^3 rho d(ln r).
        let (r0, r1, n) = (1e-5f64, 30.0f64, 20000usize);
        let h = (r1 / r0).ln() / n as f64;
        let mut total = 0.0;
        for i in 0..n {
            let r = r0 * ((i as f64 + 0.5) * h).exp();
            total += 4.0 * std::f64::consts::PI * r.powi(3) * fa.density(r) * h;
        }
        println!("  Z = {z}: the tabulated density holds {total:.6} electrons");
        assert!((total - z as f64).abs() < 2e-3, "Z = {z} holds {total}");
    }
}

/// Water's density shared among its atoms: charges that add to nothing, an
/// oxygen that takes electrons from two equal hydrogens, atom-pair `C6`s that
/// add to the molecule's, and the dipole of the density itself.
#[test]
fn water_is_shared_among_its_atoms() {
    let mol = water();
    let p = mol.problem(Functional::Pbe);
    let s = solve(&p, 200, 1e-10);
    let free = vec![FreeAtom::new(1, Functional::Pbe), FreeAtom::new(8, Functional::Pbe)];
    let part = partition(&mol, &p, &s, &free, 50, 12, Z_AB_DF1);
    let total: f64 = part.charge.iter().sum();
    println!("  charges O {:+.4} H {:+.4} {:+.4} (sum {total:+.1e}); volumes {:.2} {:.2} {:.2} bohr^3", part.charge[0], part.charge[1], part.charge[2], part.volume[0], part.volume[1], part.volume[2]);
    assert!(total.abs() < 1e-4, "the charges add to {total}");
    assert!(part.charge[0] < 0.0 && part.charge[1] > 0.0, "oxygen takes electrons from hydrogen");
    assert!((part.charge[1] - part.charge[2]).abs() < 1e-6, "the two hydrogens are alike");
    let pairs = atom_pair_c6(&part, &part);
    let sum: f64 = pairs.iter().flatten().sum();
    let whole = c6(&part.sites, &part.sites);
    println!("  C6: O-O {:.3}, O-H {:.3}, H-H {:.3}; sum {sum:.6} against the molecule's {whole:.6}", pairs[0][0], pairs[0][1], pairs[1][1]);
    assert!((sum - whole).abs() < 1e-10 * whole, "the atom pairs add to {sum}, the molecule to {whole}");
    let mu = (part.dipole[0].powi(2) + part.dipole[1].powi(2) + part.dipole[2].powi(2)).sqrt() * 2.541746;
    let from_charges = {
        let mut m = [0.0f64; 3];
        for (q, c) in part.charge.iter().zip(&mol.positions) {
            for k in 0..3 {
                m[k] += q * c[k];
            }
        }
        (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt() * 2.541746
    };
    println!("  dipole of the density {mu:.3} D against 1.85 measured; of the Hirshfeld charges alone {from_charges:.3} D");
}
