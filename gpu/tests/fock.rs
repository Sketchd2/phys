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

/// The same table kept whole on the card, split between the card and memory,
/// and split three ways with the last on disk: the same bits out, and the
/// CPU's to round-off. Equality is asserted, not closeness — the segments are
/// visited in the same order whatever the tiers are.
#[test]
fn the_tiers_give_the_same_bits() {
    use phys::electrons::scf::SpillTo;
    use phys_gpu::fock::{Budget, GpuFock};
    let mol = dimer();
    let p = mol.problem(Functional::Pbe);
    let all: Vec<usize> = (0..6).collect();
    let fc = frozen_core(&mol.z, &all);
    set_fock_engine(None);
    let hf_cpu = hartree_fock(&p, 100, 1e-10);
    let dir = std::env::temp_dir().join("phys-fock-tiers-test");
    let base = Budget { vram: 0, ram: Some(0), disk: Some(SpillTo { dir: dir.clone(), cap_bytes: 20_000_000_000 }), segment_bytes: 96 << 20, x_bytes: 100 << 20 };
    let whole = Budget { vram: 8 << 30, ..base.clone() };
    let two = Budget { vram: 700 << 20, ram: Some(2 << 30), ..base.clone() };
    let three = Budget { vram: 700 << 20, ram: Some(300 << 20), ..base.clone() };
    let mut results = Vec::new();
    let mut plans = Vec::new();
    for (label, b) in [("whole", whole), ("card+memory", two), ("card+memory+disk", three)] {
        let engine = std::sync::Arc::new(GpuFock::with_budget(b).expect("a GPU with double precision"));
        set_fock_engine(Some(engine.clone()));
        let t = Instant::now();
        let hf = hartree_fock(&p, 100, 1e-10);
        let mp = mp2(&p, &hf, fc);
        println!("  {label}: HF {:.12}, MP2 {:.12}, {:.1} s, {:?}", hf.energy, mp.correlation, t.elapsed().as_secs_f64(), engine.current_plan());
        plans.push(engine.current_plan().expect("the table was kept"));
        results.push((hf.energy, mp.correlation));
        set_fock_engine(None);
    }
    let (pw, p2, p3) = (&plans[0], &plans[1], &plans[2]);
    assert_eq!(pw.on_card, pw.segments);
    assert!(p2.on_card > 0 && p2.in_memory > 0 && p2.on_disk == 0, "{p2:?}");
    assert!(p3.on_card > 0 && p3.in_memory > 0 && p3.on_disk > 0, "{p3:?}");
    // Nowhere to put the rest: the engine says no and the CPU does it.
    let engine = std::sync::Arc::new(GpuFock::with_budget(Budget { disk: None, ..plan_budget(&base) }).expect("a GPU"));
    set_fock_engine(Some(engine.clone()));
    let hf_declined = hartree_fock(&p, 100, 1e-10);
    set_fock_engine(None);
    assert!(engine.current_plan().is_none(), "kept a table it had no room for");
    assert_eq!(hf_declined.energy, hf_cpu.energy);
    assert!((results[0].0 - hf_cpu.energy).abs() < 1e-8);
    assert_eq!(results[0], results[1], "the card and memory differ from the card alone");
    assert_eq!(results[0], results[2], "the file differs from the card alone");
    let _ = std::fs::remove_dir_all(&dir);
}

fn plan_budget(base: &phys_gpu::fock::Budget) -> phys_gpu::fock::Budget {
    phys_gpu::fock::Budget { vram: 700 << 20, ram: Some(300 << 20), ..base.clone() }
}
