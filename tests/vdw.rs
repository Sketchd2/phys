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

/// `C6` from the non-local functional, against the published vdW-DF values
/// (Vydrov and Van Voorhis, arXiv:1004.4850, Table II): neon and water, with
/// vdW-DF1's `Z_ab` and with vdW-DF2's, which differ by 2.2x.
///
/// The published values were computed on densities from a range-separated
/// hybrid with the right long-range exchange; these are on PBE's, whose tails
/// are too diffuse, and `C6` reads the tails. So the engine comes out high,
/// by about the same for both `Z_ab` (Ne 3.4% and 3.1%, water 11.4% and
/// 7.7%), which is the density and not the functional. What the test can
/// catch is everything else: `Z_ab = 0` or of the wrong sign puts `C6` out by
/// a factor of 1e6 to 1e14, because `q0` falls to nothing in the tails.
#[test]
fn c6_coefficients_match_the_published_ones_for_both_z_ab() {
    use phys::electrons::functional::Functional;
    use phys::electrons::molecule::Molecule;
    use phys::electrons::scf::{solve, Batches};
    use phys::electrons::vdw::{c6, density_on, sites, Z_AB_DF1};
    let a = 1.0 / 0.529177210903;
    let th = 104.52f64.to_radians() / 2.0;
    let water = Molecule { z: vec![8, 1, 1], positions: vec![[0.0; 3], [0.9572 * a * th.sin(), 0.0, 0.9572 * a * th.cos()], [-0.9572 * a * th.sin(), 0.0, 0.9572 * a * th.cos()]], charge: 0, unpaired: 0 };
    let neon = Molecule { z: vec![10], positions: vec![[0.0; 3]], charge: 0, unpaired: 0 };
    // (name, molecule, published vdW-DF-04, published vdW-DF-10 = vdW-DF2)
    for (name, mol, df1, df2) in [("Ne", &neon, 9.45, 3.07), ("H2O", &water, 46.96, 17.17)] {
        let p = mol.problem(Functional::Pbe);
        let s = solve(&p, 200, 1e-10);
        let mut d = s.density_alpha.clone();
        for k in 0..d.a.len() {
            d.a[k] += s.density_beta.a[k];
        }
        let atoms: Vec<([f64; 3], f64)> = p.nuclei.iter().zip(&p.sizes).map(|((_, q), r)| (*q, *r)).collect();
        // C6 is converged to four digits at this grid: 52.307, 52.309, 52.308
        // for water at (50, 12), (75, 18), (100, 24).
        let grid = phys::electrons::grid::molecular_pruned(&atoms, 50, 12, false);
        let b = Batches::new(&p.basis, &grid);
        let dens = density_on(&p.basis, &b, &d, grid.points.len());
        for (z, published) in [(Z_AB_DF1, df1), (-1.887, df2)] {
            let ss = sites(&grid, &dens, z, 0.0);
            let got = c6(&ss, &ss);
            let ratio = got / published;
            println!("  {name:<4} Z_ab {z:+.4}: C6 {got:.3} against published {published:.2} ({:+.1}%)", (ratio - 1.0) * 100.0);
            assert!((1.0..1.15).contains(&ratio), "{name} at Z_ab {z}: {got} against {published}");
        }
    }
}

/// The double sum through the tabulated kernel, between two neon atoms far
/// apart, against the same sum through the kernel's asymptote at every pair's
/// true distance: they must agree, which checks the table's use, its edge and
/// its switch to the asymptote together. What is left between that and
/// `-C6 / R^6` is the atoms' size, which falls as `R^-2` (averaging
/// `|r - r'|^-6` over two clouds adds about `5 <x^2> / R^2`); checking that it
/// does is what tells it from an error. Measured: 1.10782 at R = 20 bohr and
/// 1.05250 at 28, where the size alone predicts 1.0552.
#[test]
fn the_double_sum_approaches_c6_over_r6_as_the_size_of_the_atoms_says() {
    use phys::electrons::functional::Functional;
    use phys::electrons::molecule::Molecule;
    use phys::electrons::scf::{solve, Batches};
    use phys::electrons::vdw::{c6, density_on, nonlocal_energy, sites, Site, Z_AB_DF1};
    let neon = Molecule { z: vec![10], positions: vec![[0.0; 3]], charge: 0, unpaired: 0 };
    let p = neon.problem(Functional::Pbe);
    let s = solve(&p, 200, 1e-10);
    let mut d = s.density_alpha.clone();
    for k in 0..d.a.len() {
        d.a[k] += s.density_beta.a[k];
    }
    let atoms: Vec<([f64; 3], f64)> = p.nuclei.iter().zip(&p.sizes).map(|((_, q), r)| (*q, *r)).collect();
    let grid = phys::electrons::grid::molecular_pruned(&atoms, 50, 12, false);
    let b = Batches::new(&p.basis, &grid);
    let one = sites(&grid, &density_on(&p.basis, &b, &d, grid.points.len()), Z_AB_DF1, 0.0);
    let c = c6(&one, &one);
    let table = KernelTable::build(128, 64.0, &converged_quadrature());
    let alone = nonlocal_energy(&one, &table);
    let mut excess = Vec::new();
    for r in [20.0f64, 28.0] {
        let other: Vec<Site> = one.iter().map(|x| Site { r: [x.r[0], x.r[1], x.r[2] + r], ..*x }).collect();
        let both: Vec<Site> = one.iter().chain(&other).cloned().collect();
        let cross = nonlocal_energy(&both, &table) - 2.0 * alone;
        let mut asym = 0.0;
        for x in &one {
            for y in &other {
                let dd = [x.r[0] - y.r[0], x.r[1] - y.r[1], x.r[2] - y.r[2]];
                let rr = (dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2]).sqrt();
                asym += x.wn * y.wn * phi_asymptotic(x.q * rr, y.q * rr);
            }
        }
        let over = cross / (-c / r.powi(6));
        println!("  R = {r}: cross {cross:.6e}, asymptote at true distances {asym:.6e} ({:.7}), against -C6/R^6 {over:.5}", cross / asym);
        assert!((cross / asym - 1.0).abs() < 1e-5, "at R = {r} the table's sum is {} of the asymptote's", cross / asym);
        excess.push(over - 1.0);
    }
    let scaling = excess[1] / excess[0] / (20.0f64 / 28.0).powi(2);
    println!("  the excess over C6/R^6 scales as R^-2 to within {:.1}%", (scaling - 1.0) * 100.0);
    assert!((scaling - 1.0).abs() < 0.08, "the excess does not fall as R^-2: {scaling}");
}

/// The semilocal partners of the non-local correlation: each is LDA
/// correlation with an exchange enhancement `F(s)`, so its energy density less
/// LDA's is the uniform gas's exchange times `F(s) - 1`. The forms are written
/// out here independently of `functional.rs`: PBE's with `kappa = 0.804`,
/// revPBE's with `kappa = 1.245`, and refit PW86's
/// `(1 + 1.851 s^2 + 17.33 s^4 + 0.163 s^6)^(1/15)`.
#[test]
fn the_exchange_partners_are_the_published_forms() {
    use phys::electrons::functional::{evaluate, Functional};
    let mu = 0.2195149727645171;
    let pbe = |k: f64| move |s2: f64| 1.0 + k - k / (1.0 + mu * s2 / k);
    let rpw86 = |s2: f64| (1.0 + 1.851 * s2 + 17.33 * s2 * s2 + 0.163 * s2 * s2 * s2).powf(1.0 / 15.0);
    let forms: [(Functional, &dyn Fn(f64) -> f64); 3] = [(Functional::PbeXLdaC, &pbe(0.804)), (Functional::RevPbeXLdaC, &pbe(1.245)), (Functional::Rpw86XLdaC, &rpw86)];
    for &(n, s) in &[(0.3, 0.2), (0.05, 1.0), (0.002, 2.5)] {
        let kf = (3.0 * PI * PI * n).powf(1.0 / 3.0);
        let g2 = (2.0 * kf * n * s).powi(2);
        let ex_unif = -0.75 * (3.0 / PI).powf(1.0 / 3.0) * n.powf(4.0 / 3.0);
        // Spin halves, each carrying a quarter of the squared gradient.
        let (lda, _) = evaluate(Functional::Lda, n / 2.0, n / 2.0, g2 / 4.0, g2 / 4.0, g2 / 4.0);
        for (f, form) in forms.iter() {
            let (e, _) = evaluate(*f, n / 2.0, n / 2.0, g2 / 4.0, g2 / 4.0, g2 / 4.0);
            let want = ex_unif * (form(s * s) - 1.0);
            println!("  {f:?} n {n} s {s}: {:.10e} against {want:.10e}", e - lda);
            assert!((e - lda - want).abs() < 1e-12 * want.abs().max(1e-12), "{f:?} at n {n}, s {s}");
        }
    }
}

/// The non-local potential is the derivative of the non-local energy: a small
/// change of the density matrix moves the energy by `sum V_mn dD_mn`, checked
/// by central differences on water. And the energy computed row by row, for
/// the derivatives, is the energy computed pair by pair.
#[test]
fn the_nonlocal_potential_is_the_derivative_of_its_energy() {
    use phys::electrons::functional::Functional;
    use phys::electrons::molecule::Molecule;
    use phys::electrons::scf::{solve, Batches};
    use phys::electrons::vdw::{density_and_gradient_on, nonlocal, nonlocal_energy, nonlocal_matrix, sites, Z_AB_DF1};
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
    let grid = phys::electrons::grid::molecular_pruned(&atoms, 30, 8, false);
    let b = Batches::new(&p.basis, &grid);
    let table = KernelTable::build(128, 64.0, &converged_quadrature());
    let energy_of = |dm: &phys::electrons::linalg::Matrix| {
        let dens = density_and_gradient_on(&p.basis, &b, dm, grid.points.len());
        nonlocal(&grid, &dens, Z_AB_DF1, 0.0, &table)
    };
    let dens = density_and_gradient_on(&p.basis, &b, &d, grid.points.len());
    let nl = nonlocal(&grid, &dens, Z_AB_DF1, 0.0, &table);
    let by_pairs = nonlocal_energy(&sites(&grid, &dens.iter().map(|(n, g)| (*n, g[0] * g[0] + g[1] * g[1] + g[2] * g[2])).collect::<Vec<_>>(), Z_AB_DF1, 0.0), &table);
    println!("  E_nl {:.12} by rows, {by_pairs:.12} by pairs", nl.energy);
    assert!((nl.energy - by_pairs).abs() < 1e-11 * by_pairs.abs());
    let v = nonlocal_matrix(&p.basis, &b, &dens, &nl);
    // A physical change, water's LDA density less its PBE one. A random
    // change spread over every element was tried first and made the density
    // negative in the far tails, where points are then dropped, so points
    // came and went with the step and the difference did not converge (its
    // error changed sign between steps of 1e-2 and 3e-2). With this one it
    // falls as h^2: 8.7e-5, 8.2e-6, 9.0e-7, 1.1e-7, 2.0e-8 at h = 1 to 0.01.
    let lda = {
        let mut q = p.clone();
        q.functional = Functional::Lda;
        solve(&q, 200, 1e-10)
    };
    let mut dd = lda.density_alpha.clone();
    for k in 0..dd.a.len() {
        dd.a[k] += lda.density_beta.a[k] - d.a[k];
    }
    let predicted: f64 = v.a.iter().zip(&dd.a).map(|(x, y)| x * y).sum();
    let h = 1e-2;
    let shifted = |sign: f64| {
        let mut m = d.clone();
        for k in 0..m.a.len() {
            m.a[k] += sign * h * dd.a[k];
        }
        energy_of(&m).energy
    };
    let measured = (shifted(1.0) - shifted(-1.0)) / (2.0 * h);
    println!("  dE along LDA less PBE: {measured:.10e} by differences, {predicted:.10e} from the potential");
    assert!((measured - predicted).abs() < 1e-6 * predicted.abs(), "{measured} against {predicted}");
}

/// A counterpoise partner: the ghost atoms keep their functions and lose
/// their nuclei, their electrons and their share of the starting guess, so
/// the guess holds exactly the partner's electrons.
#[test]
fn a_ghosted_partner_starts_with_its_own_electrons() {
    use phys::electrons::functional::Functional;
    use phys::electrons::molecule::Molecule;
    let a = 1.0 / 0.529177210903;
    let pos = [[-1.551007, -0.114520, 0.0], [-1.934259, 0.762503, 0.0], [-0.599677, 0.040712, 0.0], [1.350625, 0.111469, 0.0], [1.680398, -0.373741, -0.758561], [1.680398, -0.373741, 0.758561]];
    let dimer = Molecule { z: vec![8, 1, 1, 8, 1, 1], positions: pos.iter().map(|p| [p[0] * a, p[1] * a, p[2] * a]).collect(), charge: 0, unpaired: 0 };
    let base = dimer.problem(Functional::Pbe);
    let partner = base.with_ghosts(&[3, 4, 5], 10.0);
    assert_eq!(partner.basis.size, base.basis.size, "the ghosts keep their functions");
    assert_eq!(partner.nuclei.iter().map(|n| n.0).collect::<Vec<_>>(), vec![0.0, 0.0, 0.0, 8.0, 1.0, 1.0]);
    assert_eq!((partner.alpha, partner.beta), (5.0, 5.0));
    // The guess's electron count: trace(D S).
    let (s, _, _) = phys::electrons::integrals::one_electron(&partner.basis, &partner.nuclei);
    let (da, db) = partner.guess.as_ref().expect("a guess");
    let count: f64 = (0..s.n * s.n).map(|k| (da.a[k] + db.a[k]) * s.a[k]).sum();
    let whole: f64 = {
        let (da, db) = base.guess.as_ref().unwrap();
        (0..s.n * s.n).map(|k| (da.a[k] + db.a[k]) * s.a[k]).sum()
    };
    println!("  the guess holds {count:.6} electrons for the partner, {whole:.6} for the dimer");
    assert!((count - 10.0).abs() < 1e-6 * whole, "the partner's guess holds {count} electrons");
}
