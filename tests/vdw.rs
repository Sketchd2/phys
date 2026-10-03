//! The van der Waals density functional's kernel (`electrons::vdw`), against
//! what it must be whatever the numerics: its large-separation limit, and a
//! uniform electron gas getting no non-local energy at all.

use phys::electrons::vdw::{converged_quadrature, phi_asymptotic, KernelQuadrature, KernelTable, ASYMPTOTE_C, ASYMPTOTIC_FROM};
use std::f64::consts::PI;

/// The quadrature has converged: the stored kernel's (1024 points to
/// `a_max = 256`) against four times the points over twice the range.
#[test]
fn the_kernel_quadrature_has_converged() {
    let stored = converged_quadrature();
    let check = KernelQuadrature::new(4096, 512.0);
    let mut worst: f64 = 0.0;
    for &(d1, d2) in &[(0.1, 0.1), (0.5, 1.0), (1.0, 1.0), (2.0, 3.0), (4.0, 4.0), (1.0, 6.0), (0.3, 11.0)] {
        let (a, b) = (stored.phi(d1, d2), check.phi(d1, d2));
        worst = worst.max((a - b).abs());
        println!("  phi({d1}, {d2}) = {a:+.10}  against {b:+.10}");
    }
    assert!(worst < 1e-8, "the stored kernel's quadrature is {worst:.2e} from a finer one");
}

/// Far apart, two pieces of density attract as `R^-6`: where both arguments
/// are at least `ASYMPTOTIC_FROM`, `phi` is `-C / (d1^2 d2^2 (d1^2 + d2^2))`
/// with `C = 12 (4 pi / 9)^3`; and where one is small it is not, which is why
/// the switch tests the smaller.
#[test]
fn the_kernel_falls_off_as_the_sixth_power() {
    let q = converged_quadrature();
    // Dion's constant, written out here rather than taken from the module:
    // a wrong `gamma` rescales the kernel's arguments, which leaves the zero
    // integral at zero and moves the module's own asymptote with it, so only
    // a value stated independently can catch it (checked: `gamma` 1% high
    // passed every other test).
    let c = 12.0 * (4.0 * PI / 9.0f64).powi(3);
    assert!((ASYMPTOTE_C - c).abs() < 1e-12 * c, "the module's C is {ASYMPTOTE_C}, Dion's {c}");
    for &(d1, d2) in &[(12.0, 12.0), (12.0, 18.0), (14.0, 14.0)] {
        let dion = -c / (d1 * d1 * d2 * d2 * (d1 * d1 + d2 * d2));
        let ratio = q.phi(d1, d2) / dion;
        println!("  phi({d1}, {d2}) / asymptote = {ratio:.6}");
        assert!((ratio - 1.0).abs() < 2e-4, "({d1}, {d2}): {ratio}");
        assert!(d1.min(d2) >= ASYMPTOTIC_FROM);
    }
    let lopsided = q.phi(23.76, 0.24) / phi_asymptotic(23.76, 0.24);
    println!("  phi(23.76, 0.24) / asymptote = {lopsided:.4}: one small argument, nothing like it");
    assert!(lopsided < 0.1);
    println!("  C = {ASYMPTOTE_C:.6}");
}

/// The table the double sum reads, against the quadrature it was built from,
/// at points between its nodes.
#[test]
fn the_kernel_table_interpolates_to_its_quadrature() {
    let q = converged_quadrature();
    let table = KernelTable::build(128, 64.0, &q);
    let mut worst: f64 = 0.0;
    let mut seed: u64 = 1;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..60 {
        let (a, b) = (next(), next());
        let (d1, d2) = (0.05 + 14.0 * a * a, 0.05 + 14.0 * b * b);
        worst = worst.max((table.phi(d1, d2) - q.phi(d1, d2)).abs());
    }
    println!("  128 x 128 to d = 64: worst |error| {worst:.2e} over 60 points");
    // Measured: 2.3e-6 at 64 nodes, 2.8e-7 at 128, 3.6e-8 at 256 — third order.
    assert!(worst < 5e-7, "the table is {worst:.2e} from its quadrature");
}

/// In a uniform gas `q0` is the same everywhere, and the non-local energy per
/// electron is `n / (2 q0^3) ∫ 4 pi D^2 phi(D, D) dD` — which has to vanish,
/// because the uniform gas's correlation is already all in LDA. The kernel is
/// built to integrate to zero; this checks the quadrature kept that.
#[test]
fn a_uniform_gas_gets_no_nonlocal_energy() {
    let q = converged_quadrature();
    // Composite Gauss-Legendre in D on [0, cut], the asymptote beyond.
    let cut = 24.0;
    let (x, w) = phys::electrons::grid::gauss_legendre(16);
    let panels = 96;
    let (mut total, mut size) = (0.0, 0.0);
    for p in 0..panels {
        let (lo, hi) = (cut * p as f64 / panels as f64, cut * (p + 1) as f64 / panels as f64);
        for (xi, wi) in x.iter().zip(&w) {
            let d = 0.5 * (lo + hi) + 0.5 * (hi - lo) * xi;
            let f = 4.0 * PI * d * d * q.phi(d, d) * 0.5 * (hi - lo) * wi;
            total += f;
            size += f.abs();
        }
    }
    // ∫_cut^∞ 4 pi D^2 (-C / (2 D^6)) dD = -2 pi C / (3 cut^3).
    let tail = -2.0 * PI * ASYMPTOTE_C / (3.0 * cut * cut * cut);
    total += tail;
    println!("  ∫ 4 pi D^2 phi(D, D) dD = {total:+.3e}, against ∫|...| = {size:.4}, tail {tail:+.3e}");
    assert!(total.abs() < 1e-3 * size, "the uniform gas gets {total:.3e} against a scale of {size:.3}");
}
