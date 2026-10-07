//! The GPU's single-precision products against the CPU's double: how far apart
//! they are, over inner dimensions as long as the correlated methods use.

use phys::electrons::hf::ProductEngine;
use phys::electrons::linalg::product_nt;
use phys_gpu::product::GpuProduct;

fn fill(len: usize, seed: &mut u64) -> Vec<f64> {
    (0..len)
        .map(|_| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((*seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        })
        .collect()
}

/// Sizes that are not multiples of the tile, so that the edges are tested too,
/// and an inner dimension of several thousand, as the exchange matrix and the
/// MP2 pair integrals have.
#[test]
fn the_gpu_product_agrees_with_the_cpu_to_single_precision() {
    let gpu = match GpuProduct::new() {
        Ok(g) => g,
        Err(e) => panic!("this test needs a GPU: {e}"),
    };
    let mut seed = 12345u64;
    for (m, n, kd) in [(5, 7, 3), (70, 130, 17), (300, 211, 1000), (257, 129, 5622), (1254, 513, 4000)] {
        let a = fill(m * kd, &mut seed);
        let b = fill(n * kd, &mut seed);
        let mut cpu = vec![0.0; m * n];
        let mut out = vec![0.0; m * n];
        product_nt(&a, &b, m, n, kd, &mut cpu);
        gpu.product_nt(&a, &b, m, n, kd, &mut out);
        let scale = (cpu.iter().map(|x| x * x).sum::<f64>() / cpu.len() as f64).sqrt();
        let worst = cpu.iter().zip(&out).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);
        let rms = (cpu.iter().zip(&out).map(|(x, y)| (x - y).powi(2)).sum::<f64>() / cpu.len() as f64).sqrt();
        println!("  {m} x {n} x {kd}: rms {:.2e} and worst {:.2e} of {:.2e} (rms of the product)", rms, worst, scale);
        assert!(rms < 2e-6 * scale, "{m} x {n} x {kd}: rms error {:.2e} against a product of rms {:.2e}", rms, scale);
        assert!(worst < 2e-5 * scale.max(1e-300) * 10.0, "{m} x {n} x {kd}: worst error {:.2e}", worst);
    }
}
