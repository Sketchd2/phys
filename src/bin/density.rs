//! `phys-density`: one molecule's Hartree-Fock electron density on its own
//! quadrature grid, for deriving a model of it (`docs/PLAY.md` E8f).
//!
//! ```sh
//! phys-density water > density-water.txt      # grow-water.state in the working directory
//! ```
//!
//! The molecule at the shape its growth relaxed it to, in the per-element
//! basis. Header lines give the atoms (`atom Z x y z`, bohr) and each occupied
//! orbital's energy, and `centre x y z` for each Boys-localised orbital above
//! the 1s cores (a bond or a lone pair); then one line a grid point: `x y z weight rho rho_no1s
//! rho_valence`. `rho_no1s` leaves out the lowest orbital for each atom
//! heavier than helium, and `rho_valence` also the next for each nitrogen,
//! oxygen and fluorine (their 2s-like orbital), as `phys-overlap` counts them.

use phys::electrons::functional::Functional;
use phys::electrons::grow::Resume;
use phys::electrons::molecule::Molecule;

fn main() {
    let name = std::env::args().nth(1).unwrap_or_else(|| "water".into());
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a growth state")).expect("readable");
    let z = phys::pairs::atoms_of(&name);
    let mol = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let p = mol.problem(Functional::Pbe);
    let hf = phys::electrons::hf::hartree_fock(&p, 200, 1e-9);
    assert!(hf.converged, "the molecule did not converge");
    for (zz, r) in z.iter().zip(&state.positions) {
        println!("atom {zz} {:.10} {:.10} {:.10}", r[0], r[1], r[2]);
    }
    for o in 0..hf.occupied {
        println!("orbital {o} {:.6}", hf.levels[o]);
    }
    let heavy = z.iter().filter(|&&q| q > 2).count();
    // Where the molecule's own bonds and lone pairs are: the centroids of its
    // occupied orbitals above the 1s cores, Boys-localised (`centre x y z`).
    let dip = phys::electrons::integrals::dipole(&p.basis, [0.0; 3]);
    let local = phys::electrons::localise::boys(&hf.c, hf.m, p.basis.size, heavy..hf.occupied, &dip);
    for c in &local.centroids {
        println!("centre {:.10} {:.10} {:.10}", c[0], c[1], c[2]);
    }
    let deep = heavy + z.iter().filter(|&&q| (7..=9).contains(&q)).count();
    let sites: Vec<([f64; 3], f64)> = p.nuclei.iter().zip(&p.sizes).map(|(nu, sz)| (nu.1, *sz)).collect();
    let grid = phys::electrons::grid::molecular_pruned(&sites, p.radial, p.theta, p.prune);
    let nb = p.basis.size;
    let mut phi = vec![0.0f64; nb];
    for (g, point) in grid.points.iter().enumerate() {
        phys::electrons::values::at(&p.basis, *point, &mut phi, None);
        let mut rho = [0.0f64; 3];
        for o in 0..hf.occupied {
            let v: f64 = (0..nb).map(|x| hf.c[x * hf.m + o] * phi[x]).sum();
            let d = 2.0 * v * v;
            rho[0] += d;
            if o >= heavy {
                rho[1] += d;
            }
            if o >= deep {
                rho[2] += d;
            }
        }
        println!("{:.8} {:.8} {:.8} {:.8e} {:.8e} {:.8e} {:.8e}", point[0], point[1], point[2], grid.weights[g], rho[0], rho[1], rho[2]);
    }
}
