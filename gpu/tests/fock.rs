//! Hartree-Fock and MP2 with the GPU's resident integrals against the CPU's,
//! on a water dimer: the same energies to round-off, and how much faster.

use phys::electrons::functional::Functional;
use phys::electrons::hf::{frozen_core, hartree_fock, mp2, set_fock_engine};
use phys::electrons::molecule::Molecule;
use std::time::Instant;

fn dimer() -> Molecule {
    let a = 1.0 / 0.529177210903;
    let pos: Vec<[f64; 3]> = [[-1.551007, -0.114520, 0.0], [-1.934259, 0.762503, 0.0], [-0.599677, 0.040712, 0.0], [1.350625, 0.111469, 0.0], [1.680398, -0.373741, -0.758561], [1.680398, -0.373741, 0.758561]].iter().map(|p| [p[0] * a, p[1] * a, p[2] * a]).collect();
    Molecule { z: vec![8, 1, 1, 8, 1, 1], positions: pos, charge: 0, unpaired: 0 }
}

#[test]
fn the_gpu_fock_engine_gives_the_cpus_energies() {
    let mol = dimer();
    let p = mol.problem(Functional::Pbe);
    let all: Vec<usize> = (0..6).collect();
    let fc = frozen_core(&mol.z, &all);
    // The CPU first, with nothing installed.
    set_fock_engine(None);
    let t = Instant::now();
    let hf_cpu = hartree_fock(&p, 100, 1e-10);
    let mp_cpu = mp2(&p, &hf_cpu, fc);
    let t_cpu = t.elapsed().as_secs_f64();
    let engine = phys_gpu::fock::GpuFock::new().expect("this test needs a GPU with double precision");
    let name = engine.name_string();
    set_fock_engine(Some(std::sync::Arc::new(engine)));
    let t = Instant::now();
    let hf_gpu = hartree_fock(&p, 100, 1e-10);
    let t_load = t.elapsed().as_secs_f64();
    let mp_gpu = mp2(&p, &hf_gpu, fc);
    // A second solve, the table already resident.
    let t = Instant::now();
    let hf_again = hartree_fock(&p, 100, 1e-10);
    let t_second = t.elapsed().as_secs_f64();
    set_fock_engine(None);
    println!("  CPU: Hartree-Fock {:.10} ({} iterations), MP2 {:.10}, {t_cpu:.1} s together", hf_cpu.energy, hf_cpu.iterations, mp_cpu.correlation);
    println!("  {name}: Hartree-Fock {:.10} ({} iterations, {t_load:.1} s with the integrals loaded), MP2 {:.10}; a second solve {t_second:.1} s", hf_gpu.energy, hf_gpu.iterations, mp_gpu.correlation);
    println!("  differences: Hartree-Fock {:.2e}, MP2 {:.2e} hartree; the second solve {:.2e}", hf_gpu.energy - hf_cpu.energy, mp_gpu.correlation - mp_cpu.correlation, hf_again.energy - hf_cpu.energy);
    assert!((hf_gpu.energy - hf_cpu.energy).abs() < 1e-8, "Hartree-Fock differs by {:.2e}", hf_gpu.energy - hf_cpu.energy);
    assert!((mp_gpu.correlation - mp_cpu.correlation).abs() < 1e-8, "MP2 differs by {:.2e}", mp_gpu.correlation - mp_cpu.correlation);
    assert!((hf_again.energy - hf_cpu.energy).abs() < 1e-8);
}
