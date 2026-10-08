//! The S22 dimers at Hartree-Fock + RI-MP2, counterpoise corrected, against
//! CCSD(T)/CBS — `docs/PLAY.md` Phase 6, E8a, stages 1 and 2.
//!
//! ```sh
//! phys-s22-gpu [Water_dimer Ammonia_dimer ...]
//! ```
//!
//! The geometries and reference energies are the S22 set's (Jurecka, Sponer,
//! Cerny and Hobza 2006; the reference energies Takatani et al. 2010): external
//! benchmark data, here only to say how good the engine's own method is, never
//! an input to it. The GPU's Hartree-Fock passes are used for any dimer whose
//! whitened integrals fit in the card's memory and the CPU's for any that do
//! not (the 6 GB card holds up to about 600 functions in double precision);
//! `PHYS_NO_GPU_FOCK=1` leaves everything on the CPU.

use phys::electrons::functional::Functional;
use phys::electrons::hf::{frozen_core, hartree_fock, mp2};
use phys::electrons::molecule::Molecule;
use std::time::Instant;

const S22: &[(&str, &[u32], &[[f64; 3]], usize, f64)] = &[
    ("Ammonia_dimer", &[7, 1, 1, 1, 7, 1, 1, 1], &[[-1.578718, -0.046611, 0.0], [-2.158621, 0.136396, -0.809565], [-2.158621, 0.136396, 0.809565], [-0.849471, 0.658193, 0.0], [1.578718, 0.046611, 0.0], [2.158621, -0.136396, -0.809565], [0.849471, -0.658193, 0.0], [2.158621, -0.136396, 0.809565]], 4, -3.1708),
    ("Water_dimer", &[8, 1, 1, 8, 1, 1], &[[-1.551007, -0.11452, 0.0], [-1.934259, 0.762503, 0.0], [-0.599677, 0.040712, 0.0], [1.350625, 0.111469, 0.0], [1.680398, -0.373741, -0.758561], [1.680398, -0.373741, 0.758561]], 3, -5.0203),
    ("Methane_dimer", &[6, 1, 1, 1, 1, 6, 1, 1, 1, 1], &[[0.0, -0.00014, 1.859161], [-0.888551, 0.51306, 1.494685], [0.888551, 0.51306, 1.494685], [0.0, -1.026339, 1.494868], [0.0, 8.9e-05, 2.948284], [0.0, 0.00014, -1.859161], [0.0, -8.9e-05, -2.948284], [-0.888551, -0.51306, -1.494685], [0.888551, -0.51306, -1.494685], [0.0, 1.026339, -1.494868]], 5, -0.5304),
    ("Ethene-ethyne_complex", &[6, 6, 1, 1, 1, 1, 6, 6, 1, 1], &[[0.0, -0.667578, -2.124659], [0.0, 0.667578, -2.124659], [0.923621, -1.232253, -2.126185], [-0.923621, -1.232253, -2.126185], [-0.923621, 1.232253, -2.126185], [0.923621, 1.232253, -2.126185], [0.0, 0.0, 2.900503], [0.0, 0.0, 1.69324], [0.0, 0.0, 0.627352], [0.0, 0.0, 3.963929]], 6, -1.5105),
    ("Ethene_dimer", &[6, 6, 1, 1, 1, 1, 6, 6, 1, 1, 1, 1], &[[-0.471925, -0.471925, -1.859111], [0.471925, 0.471925, -1.859111], [-0.872422, -0.872422, -0.936125], [0.872422, 0.872422, -0.936125], [-0.870464, -0.870464, -2.783308], [0.870464, 0.870464, -2.783308], [-0.471925, 0.471925, 1.859111], [0.471925, -0.471925, 1.859111], [-0.872422, 0.872422, 0.936125], [0.872422, -0.872422, 0.936125], [-0.870464, 0.870464, 2.783308], [0.870464, -0.870464, 2.783308]], 6, -1.4989),
    ("Formic_acid_dimer", &[6, 8, 8, 1, 1, 6, 8, 8, 1, 1], &[[-1.888896, -0.179692, 0.0], [-1.49328, 1.073689, 0.0], [-1.170435, -1.16659, 0.0], [-2.979488, -0.258829, 0.0], [-0.498833, 1.107195, 0.0], [1.888896, 0.179692, 0.0], [1.49328, -1.073689, 0.0], [1.170435, 1.16659, 0.0], [2.979488, 0.258829, 0.0], [0.498833, -1.107195, 0.0]], 5, -18.7989),
];

fn main() {
    if std::env::var_os("PHYS_NO_GPU_FOCK").is_none() {
        match phys_gpu::install_fock() {
            Ok(name) => println!("Hartree-Fock passes on {name} (where the integrals fit)"),
            Err(e) => eprintln!("no GPU for Hartree-Fock: {e}"),
        }
        match phys_gpu::install_product() {
            Ok(name) => println!("MP2 products on {name}"),
            Err(e) => eprintln!("no GPU for MP2 products: {e}"),
        }
    }
    let wanted: Vec<String> = std::env::args().skip(1).collect();
    let a = 1.0 / 0.529177210903;
    for (name, z, pos, na, reference) in S22 {
        if !wanted.is_empty() && !wanted.iter().any(|w| w == name) {
            continue;
        }
        let t0 = Instant::now();
        let p: Vec<[f64; 3]> = pos.iter().map(|p| [p[0] * a, p[1] * a, p[2] * a]).collect();
        let mol = Molecule { z: z.to_vec(), positions: p, charge: 0, unpaired: 0 };
        let base = mol.problem(Functional::Pbe);
        let nat = z.len();
        let ne_a: f64 = z[..*na].iter().map(|&x| x as f64).sum();
        let ne_b: f64 = z[*na..].iter().map(|&x| x as f64).sum();
        let mut e_hf = [0.0; 3];
        let mut e_mp2 = [0.0; 3];
        let slots: Vec<(Vec<usize>, f64)> = vec![((0..nat).collect(), ne_a + ne_b), ((0..*na).collect(), ne_a), ((*na..nat).collect(), ne_b)];
        for (slot, (keep, ne)) in slots.into_iter().enumerate() {
            let p = base.with_ghosts(&keep, ne);
            let hf = hartree_fock(&p, 150, 1e-9);
            assert!(hf.converged, "{name}: Hartree-Fock did not converge");
            let c = mp2(&p, &hf, frozen_core(z, &keep));
            e_hf[slot] = hf.energy;
            e_mp2[slot] = c.correlation;
        }
        let k = 627.509474;
        let hf_int = (e_hf[0] - e_hf[1] - e_hf[2]) * k;
        let c_int = (e_mp2[0] - e_mp2[1] - e_mp2[2]) * k;
        println!("S22 {name}: {} functions; CP HF {hf_int:.3} + MP2 correlation {c_int:.3} = {:.3} kcal/mol against CCSD(T)/CBS {reference:.3} ({:+.1}%); {:.0} s", base.basis.size, hf_int + c_int, 100.0 * ((hf_int + c_int) / reference - 1.0), t0.elapsed().as_secs_f64());
    }
}
