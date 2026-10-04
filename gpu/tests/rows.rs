//! The GPU's rows against the CPU's, on water's density. The owner's
//! precision verdict waits for E8 (PLAY.md E7); this checks only that the GPU
//! computes the same sums — a wrong port is out by orders of magnitude — and
//! reports both, with their times.
//!
//! Measured on the RTX 2060: E_nl 4.4e-5 relative from the CPU's, the core
//! rows 1.2e-4. It is single precision's own: each kernel value is good to
//! about 3e-7 absolute (a Rust copy of the shader's kernel against the f64
//! table), and a row's terms add to ten-odd in magnitude while the row comes
//! to a few hundredths, because the kernel integrates to zero. Kahan
//! summation does not move it (switched off: 2e-12 relative), and nor did
//! carrying each position as two singles.

use phys::electrons::functional::Functional;
use phys::electrons::molecule::Molecule;
use phys::electrons::scf::{solve, Batches};
use phys::electrons::vdw::{density_and_gradient_on, kernel_table, q0, CpuRows, RowEngine, RowPoint, FLOOR, Z_AB_DF1};
use std::time::Instant;

#[test]
fn the_gpu_computes_the_cpus_rows() {
    let gpu = match phys_gpu::GpuRows::new() {
        Ok(g) => g,
        Err(e) => panic!("no GPU to test: {e}"),
    };
    println!("  {}", gpu.name());
    let a = 1.0 / 0.529177210903;
    let th = 104.52f64.to_radians() / 2.0;
    let water = Molecule { z: vec![8, 1, 1], positions: vec![[0.0; 3], [0.9572 * a * th.sin(), 0.0, 0.9572 * a * th.cos()], [-0.9572 * a * th.sin(), 0.0, 0.9572 * a * th.cos()]], charge: 0, unpaired: 0 };
    let p = water.problem(Functional::Pbe);
    let s = solve(&p, 200, 1e-10);
    let mut d = s.density_alpha.clone();
    for k in 0..d.a.len() {
        d.a[k] += s.density_beta.a[k];
    }
    let atoms: Vec<([f64; 3], f64)> = p.nuclei.iter().zip(&p.sizes).map(|((_, q), r)| (*q, *r)).collect();
    let grid = phys::electrons::grid::molecular_pruned(&atoms, 50, 12, false);
    let b = Batches::new(&p.basis, &grid);
    let dens = density_and_gradient_on(&p.basis, &b, &d, grid.points.len());
    let points: Vec<RowPoint> = grid.points.iter().zip(&grid.weights).zip(&dens).filter_map(|((r, w), (n, g))| {
        let wn = w * n;
        (*n > 0.0 && wn.abs() >= FLOOR).then(|| RowPoint { r: *r, wn, q: q0(*n, g[0] * g[0] + g[1] * g[1] + g[2] * g[2], Z_AB_DF1) })
    }).collect();
    let table = kernel_table();
    let t = Instant::now();
    let cpu = CpuRows.rows(&points, table, true);
    let t_cpu = t.elapsed().as_secs_f64();
    let _ = gpu.rows(&points[..256.min(points.len())], table, true);
    let t = Instant::now();
    let on_gpu = gpu.rows(&points, table, true);
    let t_gpu = t.elapsed().as_secs_f64();
    let energy = |rows: &[phys::electrons::vdw::RowSums]| 0.5 * points.iter().zip(rows).map(|(p, r)| p.wn * r.a).sum::<f64>();
    let (e_cpu, e_gpu) = (energy(&cpu), energy(&on_gpu));
    let worst = |f: &dyn Fn(&phys::electrons::vdw::RowSums) -> f64| {
        let scale = cpu.iter().map(|r| f(r).abs()).fold(0.0, f64::max);
        cpu.iter().zip(&on_gpu).map(|(x, y)| (f(x) - f(y)).abs()).fold(0.0, f64::max) / scale
    };
    let (wa, wb, wc) = (worst(&|r| r.a), worst(&|r| r.b), worst(&|r| r.c[0].abs() + r.c[1].abs() + r.c[2].abs()));
    println!("  {} points: CPU {t_cpu:.2} s, GPU {t_gpu:.3} s ({:.0}x)", points.len(), t_cpu / t_gpu);
    println!("  E_nl {e_cpu:.12} on the CPU, {e_gpu:.12} on the GPU, relative {:.2e}", (e_gpu - e_cpu) / e_cpu);
    println!("  worst row against the largest: a {wa:.2e}, b {wb:.2e}, c {wc:.2e}");
    assert!(((e_gpu - e_cpu) / e_cpu).abs() < 1e-4, "the GPU's energy is {e_gpu} against {e_cpu}");
    assert!(wa < 1e-3 && wb < 1e-3 && wc < 1e-3, "a row differs by more than a single-precision port can: {wa:.2e} {wb:.2e} {wc:.2e}");
}
