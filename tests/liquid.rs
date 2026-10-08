//! Rigid molecules in a periodic box (`liquid`), against what the mechanics
//! must do whatever the law. The law here is artificial — Lennard-Jones on
//! oxygen and fixed charges, roughly SPC water — and only exercises the
//! machinery; the laws the boiling points use come from the electronic
//! structure.

use phys::liquid::{switch, Forces, Kind, Liquid, SiteLaw, AMU, K_B};
use phys::math::{Quat, Vec3};

const A: f64 = 1.0 / 0.529177210903;

struct Spc;
impl SiteLaw for Spc {
    fn site_pair(&self, ta: usize, tb: usize, r: f64) -> (f64, f64) {
        let q = [-0.82, 0.41];
        let (mut u, mut du) = (q[ta] * q[tb] / r, -q[ta] * q[tb] / (r * r));
        if ta == 0 && tb == 0 {
            let (sigma, eps) = (3.166 * A, 0.1553 / 627.5095);
            let s6 = (sigma / r).powi(6);
            u += 4.0 * eps * (s6 * s6 - s6);
            du += 4.0 * eps * (-12.0 * s6 * s6 + 6.0 * s6) / r;
        }
        (u, du)
    }
}

fn water_kind() -> Kind {
    let th = 109.47f64.to_radians() / 2.0;
    let r = 1.0 * A;
    let pos = [[0.0, 0.0, 0.0], [r * th.sin(), 0.0, r * th.cos()], [-r * th.sin(), 0.0, r * th.cos()]];
    Kind::from_atoms(&pos, &[15.999 * AMU, 1.008 * AMU, 1.008 * AMU], &[0, 1, 1])
}

/// A cubic lattice of `side^3` molecules, oriented and moving by a seeded draw.
fn lattice(side: usize, spacing: f64, temperature: f64) -> Liquid {
    let kind = water_kind();
    let mut s = Liquid::noise(11);
    let (mut com, mut orientation, mut momentum, mut spin) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for i in 0..side {
        for j in 0..side {
            for k in 0..side {
                com.push(Vec3 { x: i as f64 * spacing, y: j as f64 * spacing, z: k as f64 * spacing });
                let axis = s.normal3().unit();
                orientation.push(Quat::from_axis_angle(axis, s.uniform() * 6.283));
                momentum.push(s.normal3().scale((kind.mass * K_B * temperature).sqrt()));
                spin.push(Vec3 {
                    x: s.normal() * (kind.inertia[0] * K_B * temperature).sqrt(),
                    y: s.normal() * (kind.inertia[1] * K_B * temperature).sqrt(),
                    z: s.normal() * (kind.inertia[2] * K_B * temperature).sqrt(),
                });
            }
        }
    }
    let n = momentum.len() as f64;
    let mean = momentum.iter().fold(Vec3::ZERO, |a, p| a + *p).scale(1.0 / n);
    for p in momentum.iter_mut() {
        *p -= mean;
    }
    let edge = side as f64 * spacing;
    Liquid { kind, cell: [edge; 3], com, orientation, momentum, spin, r_on: 0.45 * edge - 2.0 * A, r_cut: 0.45 * edge }
}

/// A kind is centred on its mass and turned onto its principal axes: its
/// sites keep their distances, and its inertia tensor there is diagonal.
#[test]
fn a_kind_sits_on_its_principal_axes() {
    let k = water_kind();
    let d = |a: usize, b: usize| (k.sites[a] - k.sites[b]).norm();
    assert!((d(0, 1) - 1.0 * A).abs() < 1e-12 && (d(0, 2) - 1.0 * A).abs() < 1e-12, "the shape is kept");
    let masses = [15.999 * AMU, 1.008 * AMU, 1.008 * AMU];
    let com = k.sites.iter().zip(&masses).fold(Vec3::ZERO, |a, (s, m)| a + s.scale(*m));
    assert!(com.norm() < 1e-9, "centred on its mass");
    let offdiag = |i: usize, j: usize| k.sites.iter().zip(&masses).map(|(s, m)| -m * [s.x, s.y, s.z][i] * [s.x, s.y, s.z][j]).sum::<f64>();
    let worst = offdiag(0, 1).abs().max(offdiag(0, 2).abs()).max(offdiag(1, 2).abs());
    println!("  inertia {:.4e} {:.4e} {:.4e}; largest off-diagonal {worst:.1e}", k.inertia[0], k.inertia[1], k.inertia[2]);
    assert!(worst < 1e-9 * k.inertia[2], "the inertia tensor is diagonal on the principal axes");
}

/// A free asymmetric rotor keeps its energy and its angular momentum in the
/// box frame — which tests the splitting into rotations about each axis.
#[test]
fn a_free_rotor_keeps_its_energy_and_angular_momentum() {
    let kind = water_kind();
    let spin = Vec3 { x: 3.0, y: -2.0, z: 1.5 };
    let mut l = Liquid { kind, cell: [1e6; 3], com: vec![Vec3::ZERO], orientation: vec![Quat::IDENTITY], momentum: vec![Vec3::ZERO], spin: vec![spin], r_on: 1.0, r_cut: 2.0 };
    let lab = |l: &Liquid| l.orientation[0].rotate(l.spin[0]);
    let (e0, l0) = (l.kinetic().1, lab(&l));
    let mut f = l.forces(&Spc);
    for _ in 0..20_000 {
        f = l.step(&Spc, f, 5.0, None);
    }
    let (e1, l1) = (l.kinetic().1, lab(&l));
    println!("  rotational energy {e0:.10e} -> {e1:.10e}; angular momentum moved {:.1e}", (l1 - l0).norm() / l0.norm());
    assert!((e1 - e0).abs() < 1e-6 * e0, "the free rotor's energy drifted");
    assert!((l1 - l0).norm() < 1e-9 * l0.norm(), "the free rotor's angular momentum moved");
}

/// Without a thermostat a periodic liquid keeps its energy (no drift past the
/// integrator's own oscillation) and its total momentum.
#[test]
fn a_liquid_without_a_bath_keeps_its_energy_and_momentum() {
    let mut l = lattice(4, 3.1 * A, 300.0);
    let total = |l: &Liquid, f: &Forces| {
        let (t, r) = l.kinetic();
        t + r + f.energy
    };
    let mut f = l.forces(&Spc);
    let e0 = total(&l, &f);
    let mut energies = Vec::new();
    for _ in 0..2000 {
        f = l.step(&Spc, f, 20.0, None);
        energies.push(total(&l, &f));
    }
    let p = l.momentum.iter().fold(Vec3::ZERO, |a, x| a + *x).norm();
    let (lo, hi) = energies.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), e| (a.min(*e), b.max(*e)));
    let first = energies[..200].iter().sum::<f64>() / 200.0;
    let last = energies[energies.len() - 200..].iter().sum::<f64>() / 200.0;
    let (kt, kr) = l.kinetic();
    println!("  64 molecules, 2000 steps of 20 a.u.: energy {e0:.8} -> mean {last:.8}, spread {:.2e}, drift {:.2e}; kinetic {:.2e}; momentum {p:.1e}", hi - lo, last - first, kt + kr);
    assert!((last - first).abs() < 1e-3 * (kt + kr), "the energy drifted by {:.2e}", last - first);
    assert!(p < 1e-9, "total momentum {p:.2e}");
}

/// A bath brings both kinds of motion to its temperature: three halves kT a
/// molecule in translation and in rotation alike.
#[test]
fn a_bath_brings_translation_and_rotation_to_its_temperature() {
    let mut l = lattice(4, 3.1 * A, 100.0);
    let mut stream = Liquid::noise(3);
    let target = 350.0;
    let mut f = l.forces(&Spc);
    let (mut tr, mut rot, mut count) = (0.0, 0.0, 0.0);
    for step in 0..6000 {
        f = l.step(&Spc, f, 20.0, Some((target, 1.0 / 2000.0, &mut stream)));
        if step >= 2000 {
            let (t, r) = l.kinetic();
            tr += t;
            rot += r;
            count += 1.0;
        }
    }
    let n = l.com.len() as f64;
    let (t_tr, t_rot) = (2.0 * tr / count / (3.0 * n * K_B), 2.0 * rot / count / (3.0 * n * K_B));
    println!("  target {target} K: translation {t_tr:.1} K, rotation {t_rot:.1} K");
    assert!((t_tr / target - 1.0).abs() < 0.05 && (t_rot / target - 1.0).abs() < 0.05, "translation {t_tr:.1} K, rotation {t_rot:.1} K");
}

/// The switch is one inside, zero outside, smooth between.
#[test]
fn the_switch_is_smooth() {
    let (on, cut) = (10.0, 12.0);
    assert_eq!(switch(9.0, on, cut), (1.0, 0.0));
    assert_eq!(switch(12.5, on, cut), (0.0, 0.0));
    let h = 1e-6;
    for r in [10.3, 11.0, 11.7] {
        let fd = (switch(r + h, on, cut).0 - switch(r - h, on, cut).0) / (2.0 * h);
        assert!((fd - switch(r, on, cut).1).abs() < 1e-8);
    }
}

/// The site-site law's slope is the derivative of its energy, and Tang and
/// Toennies' damping goes from nothing at contact to one far out.
#[test]
fn the_site_site_law_is_consistent() {
    use phys::liquid::{tang_toennies, SiteSite};
    let law = SiteSite { charge: vec![-0.7, 0.35], pair: vec![[40.0, 1.9, 15.0, 300.0], [3.0, 2.1, 4.0, 60.0], [3.0, 2.1, 4.0, 60.0], [0.5, 2.4, 1.0, 15.0]], sigma: Vec::new(), damp: Vec::new() };
    for (ta, tb) in [(0, 0), (0, 1), (1, 1)] {
        for r in [2.0, 3.5, 5.0, 9.0] {
            let h = 1e-6;
            let fd = (law.site_pair(ta, tb, r + h).0 - law.site_pair(ta, tb, r - h).0) / (2.0 * h);
            let an = law.site_pair(ta, tb, r).1;
            assert!((fd - an).abs() < 1e-7 * an.abs().max(1e-6), "({ta},{tb}) at {r}: {an} against {fd}");
        }
    }
    assert!(tang_toennies(6, 1e-3).0 < 1e-15 && (tang_toennies(6, 60.0).0 - 1.0).abs() < 1e-15);
    let h = 1e-6;
    for x in [0.5, 3.0, 8.0] {
        let fd = (tang_toennies(8, x + h).0 - tang_toennies(8, x - h).0) / (2.0 * h);
        assert!((fd - tang_toennies(8, x).1).abs() < 1e-8);
    }
}

/// The density profile counts every molecule once, wherever the slab sits.
#[test]
fn a_density_profile_counts_every_molecule() {
    let mut l = lattice(4, 3.1 * A, 300.0);
    l.cell[2] *= 3.0;
    for r in l.com.iter_mut() {
        r.z += 0.9 * l.cell[2];
    }
    let bins = 60;
    let prof = l.density_profile(bins);
    let vol = l.cell[0] * l.cell[1] * l.cell[2] / bins as f64;
    let total: f64 = prof.iter().map(|d| d * vol).sum();
    assert!((total - l.com.len() as f64).abs() < 1e-9, "{total} molecules counted");
    let middle: f64 = prof[bins / 3..2 * bins / 3].iter().map(|d| d * vol).sum();
    println!("  {} of {} molecules in the middle third", middle, l.com.len());
    assert!(middle > 0.9 * l.com.len() as f64, "the profile is centred on the slab");
}

/// The fitter recovers a law from energies that law produced: from
/// perturbed numbers it comes back to a residual at round-off and the
/// original charges.
#[test]
fn the_fitter_recovers_a_law_from_its_own_energies() {
    use phys::liquid::{fit_site_site, pair_energy, PairEnergy, SiteSite};
    let truth = SiteSite { charge: vec![-0.74, 0.37], pair: vec![[60.0, 2.0, 18.0, 350.0], [4.0, 2.2, 5.0, 70.0], [4.0, 2.2, 5.0, 70.0], [0.6, 2.5, 1.2, 16.0]], sigma: Vec::new(), damp: Vec::new() };
    let kind = water_kind();
    let mut s = Liquid::noise(5);
    let mut data = Vec::new();
    while data.len() < 120 {
        let dir = s.normal3().unit();
        let dist = s.range(5.0, 12.0);
        let (qa, qb) = (Quat::from_axis_angle(s.normal3().unit(), s.uniform() * 6.283), Quat::from_axis_angle(s.normal3().unit(), s.uniform() * 6.283));
        let a: Vec<(Vec3, usize)> = kind.sites.iter().zip(&kind.types).map(|(p, t)| (qa.rotate(*p), *t)).collect();
        let b: Vec<(Vec3, usize)> = kind.sites.iter().zip(&kind.types).map(|(p, t)| (qb.rotate(*p) + dir.scale(dist), *t)).collect();
        let closest = a.iter().flat_map(|(p, _)| b.iter().map(move |(q, _)| (*p - *q).norm())).fold(f64::INFINITY, f64::min);
        if closest < 2.8 {
            continue;
        }
        let mut pe = PairEnergy { a, b, energy: 0.0 };
        pe.energy = pair_energy(&truth, &pe);
        data.push(pe);
    }
    let mut start = truth.clone();
    start.charge = vec![-0.6, 0.3];
    for p in start.pair.iter_mut() {
        p[0] *= 1.5;
        p[1] *= 0.9;
        p[2] *= 1.3;
        p[3] *= 0.7;
    }
    let fit = fit_site_site(&data, &[1, 2], &start, 400.0, 1e-3, 3000);
    let spread = (data.iter().map(|d| d.energy * d.energy).sum::<f64>() / data.len() as f64).sqrt();
    println!("  {} iterations: rms {:.2e} against energies of rms {spread:.2e}; charges {:.5} {:.5}", fit.iterations, fit.rms, fit.law.charge[0], fit.law.charge[1]);
    assert!(fit.rms < 1e-6 * spread, "the fit left {:.2e}", fit.rms);
    assert!((fit.law.charge[0] + 0.74).abs() < 1e-3 && (fit.law.charge[1] - 0.37).abs() < 1e-3, "the charges came back as {:?}", fit.law.charge);
}

/// A slab is built at the density asked for, in the middle of its box, with
/// the vapour's regions empty to begin with.
#[test]
fn a_slab_is_built_at_its_density() {
    use phys::liquid::coexisting_densities;
    // Water's liquid density, molecules per bohr^3 (33.4 per nm^3).
    let density = 33.4e-3 * 0.529177210903f64.powi(3);
    let side = 30.0;
    let slab = Liquid::slab(water_kind(), 216, density, side, 3.0 * side, 300.0, 7, 10.0, 12.0);
    assert_eq!(slab.com.len(), 216);
    let prof = slab.density_profile(90);
    let (liquid, vapour) = coexisting_densities(&prof, 0.3, 0.2);
    println!("  liquid {:.3e} against {density:.3e} asked; vapour {vapour:.1e}", liquid);
    assert!((liquid / density - 1.0).abs() < 0.15, "the slab's middle is at {liquid:.3e}");
    assert_eq!(vapour, 0.0, "the vapour starts empty");
    let p = slab.momentum.iter().fold(Vec3::ZERO, |a, x| a + *x).norm();
    assert!(p < 1e-9, "the slab is at rest overall");
}

/// A fitted law survives being written and read back, to the bit.
#[test]
fn a_law_reads_back_exactly() {
    use phys::liquid::SiteSite;
    let law = SiteSite { charge: vec![-0.7123456789012345, 0.35617283945061725], pair: vec![[40.1, 1.93, 15.2, 301.0], [3.01, 2.11, 4.02, 60.3], [3.01, 2.11, 4.02, 60.3], [0.51, 2.42, 1.03, 15.4]], sigma: Vec::new(), damp: Vec::new() };
    assert_eq!(SiteSite::from_text(&law.to_text()), Some(law));
}

/// The molecular virial is minus the energy's derivative with respect to a
/// uniform scaling of the box and every centre in it (the arms fixed): that is
/// what makes `(2 K + W) / 3 V` the pressure. Checked by central differences
/// on 64 molecules at liquid density, where the cut-off's switch is crossed.
#[test]
fn the_virial_is_minus_the_energy_slope_under_scaling() {
    let l = lattice(4, 3.1 * A, 300.0);
    let f = l.forces(&Spc);
    let at = |lambda: f64| {
        let mut m = l.clone();
        for r in m.com.iter_mut() {
            *r = r.scale(lambda);
        }
        for c in m.cell.iter_mut() {
            *c *= lambda;
        }
        m.forces(&Spc).energy
    };
    let eps = 1e-5;
    let minus_slope = -(at(1.0 + eps) - at(1.0 - eps)) / (2.0 * eps);
    println!("  virial {:.8e}, minus dU/dlambda {minus_slope:.8e}; pressure {:.1} bar", f.virial, l.pressure(&f) * phys::liquid::BAR_PER_HARTREE_PER_BOHR3);
    assert!((f.virial - minus_slope).abs() < 1e-6 * f.virial.abs().max(1e-3), "{} against {minus_slope}", f.virial);
}

/// A massless site rides along: it does not move the centre of mass or change
/// the inertia, it sits where it was put (here on the bisector at the distance
/// asked), and it is carried into the principal frame with the atoms.
#[test]
fn a_massless_site_rides_along() {
    let th = 104.52f64.to_radians() / 2.0;
    let r = 0.9572 * A;
    let pos = [[0.0, 0.0, 0.0], [r * th.sin(), 0.0, r * th.cos()], [-r * th.sin(), 0.0, r * th.cos()]];
    let plain = Kind::of_molecule(&[8, 1, 1], &pos, &[0, 1, 1], None);
    let site = Kind::of_molecule(&[8, 1, 1], &pos, &[0, 1, 1], Some((2, 0.267)));
    assert_eq!(site.sites.len(), 4);
    for k in 0..3 {
        assert!((plain.inertia[k] - site.inertia[k]).abs() < 1e-9 * plain.inertia[k], "the inertia moved");
    }
    assert!((plain.mass - site.mass).abs() < 1e-9 * plain.mass);
    // The site is 0.267 bohr from the oxygen, and on the line from it
    // through the hydrogens' midpoint.
    let o = site.sites[0];
    let m = site.sites[3];
    assert!(((m - o).norm() - 0.267).abs() < 1e-9, "{}", (m - o).norm());
    let mid = (site.sites[1] + site.sites[2]).scale(0.5) - o;
    assert!((mid.unit().dot((m - o).unit()) - 1.0).abs() < 1e-12, "the site is off the bisector");
    assert_eq!(site.types, vec![0, 1, 1, 2]);
}

/// The off-atom site's line in a law's text.
#[test]
fn a_laws_bisector_line_is_read() {
    use phys::liquid::SiteSite;
    let text = "charges 1e0 5e-1 -2e0\npair 0 0 1e0 1e0 1e0 1e0\nbisector 2 2.67e-1\n";
    assert_eq!(SiteSite::bisector_from_text(text), Some((2, 0.267)));
    assert_eq!(SiteSite::bisector_from_text("charges 1e0\n"), None);
}

/// A law's text with its bisector line still reads as a pair law: the loader
/// once rejected the line as unknown, and a fitted law could not be run.
#[test]
fn a_law_with_a_bisector_line_still_reads() {
    use phys::liquid::SiteSite;
    let text = "charges 1e0 -1e0\npair 0 0 1e0 1e0 1e0 1e0\npair 0 1 1e0 1e0 1e0 1e0\npair 1 1 1e0 1e0 1e0 1e0\nbisector 1 2.67e-1\n";
    let law = SiteSite::from_text(text).expect("a law with a bisector line");
    assert_eq!(law.charge, vec![1.0, -1.0]);
    assert_eq!(SiteSite::bisector_from_text(text), Some((1, 0.267)));
}

/// A fit can hold numbers where they were started: a charge-only site keeps
/// zero repulsion and dispersion, a dispersion coefficient given from outside
/// stays what it was given, and everything else is still recovered. Held
/// parameters that moved would make the bisector law's "derived C6" a fitted one.
#[test]
fn a_fit_holds_what_it_is_told_to_hold() {
    use phys::liquid::{fit_site_site_held, pair_energy, Held, PairEnergy, SiteSite};
    // Types: 0 and 1 are atoms (1 and 2 of them), 2 is a charge-only site.
    let zero = [0.0, 1.0, 0.0, 0.0];
    let truth = SiteSite {
        charge: vec![0.9, 0.6, -2.1],
        pair: vec![[60.0, 2.0, 18.0, 350.0], [4.0, 2.2, 5.0, 70.0], zero, [4.0, 2.2, 5.0, 70.0], [0.6, 2.5, 1.2, 16.0], zero, zero, zero, zero],
        sigma: Vec::new(),
        damp: Vec::new(),
    };
    let kind = water_kind();
    let mut s = Liquid::noise(9);
    let mut data = Vec::new();
    let with_site = |m: Vec<(Vec3, usize)>| -> Vec<(Vec3, usize)> {
        let o = m[0].0;
        let bis = ((m[1].0 - o).unit() + (m[2].0 - o).unit()).unit();
        let mut out = m;
        out.push((o + bis.scale(0.3), 2));
        out
    };
    while data.len() < 150 {
        let dir = s.normal3().unit();
        let dist = s.range(5.0, 12.0);
        let (qa, qb) = (Quat::from_axis_angle(s.normal3().unit(), s.uniform() * 6.283), Quat::from_axis_angle(s.normal3().unit(), s.uniform() * 6.283));
        let a: Vec<(Vec3, usize)> = kind.sites.iter().zip(&kind.types).map(|(p, t)| (qa.rotate(*p), *t)).collect();
        let b: Vec<(Vec3, usize)> = kind.sites.iter().zip(&kind.types).map(|(p, t)| (qb.rotate(*p) + dir.scale(dist), *t)).collect();
        let closest = a.iter().flat_map(|(p, _)| b.iter().map(move |(q, _)| (*p - *q).norm())).fold(f64::INFINITY, f64::min);
        if closest < 2.8 {
            continue;
        }
        let mut pe = PairEnergy { a: with_site(a), b: with_site(b), energy: 0.0 };
        pe.energy = pair_energy(&truth, &pe);
        data.push(pe);
    }
    let mut start = truth.clone();
    start.charge = vec![0.5, 0.3, -0.9];
    for k in [0usize, 1, 3, 4] {
        for j in [0usize, 1, 3] {
            start.pair[k][j] *= if j == 1 { 0.9 } else { 1.4 };
        }
    }
    // The held C6 of the atom pair (0, 1) is the truth's, wherever else the start is.
    let held_c6 = truth.pair[1][2];
    let fit = fit_site_site_held(&data, &[1, 2, 1], &start, 400.0, 1e-3, 3000, &Held { pairs: &[(0, 2), (1, 2), (2, 2)], dispersion: &[(0, 1)], sigma: false, charges: false, damp: false, repulsion: false, no_dispersion: false });
    let spread = (data.iter().map(|d| d.energy * d.energy).sum::<f64>() / data.len() as f64).sqrt();
    println!("  {} iterations: rms {:.2e} against {spread:.2e}; charges {:?}; held C6 {:e}", fit.iterations, fit.rms, fit.law.charge, fit.law.pair[1][2]);
    assert!((fit.law.pair[1][2] / held_c6 - 1.0).abs() < 1e-12, "the held C6 moved to {}", fit.law.pair[1][2]);
    for (a, b) in [(0usize, 2usize), (1, 2), (2, 2)] {
        assert!(fit.law.pair[a * 3 + b][0] < 1e-250 && fit.law.pair[a * 3 + b][2] < 1e-250, "the charge-only site's pair {a}-{b} gained {:?}", fit.law.pair[a * 3 + b]);
    }
    assert!(fit.rms < 1e-4 * spread, "the fit left {:.2e}", fit.rms);
}

/// A law with induced dipoles for the tests: charges on the types, a
/// polarisability on each.
fn polarisable_law() -> (phys::liquid::SiteSite, Vec<f64>) {
    let law = phys::liquid::SiteSite { charge: vec![-0.7, 0.35], pair: vec![[40.0, 1.9, 15.0, 300.0], [3.0, 2.1, 4.0, 60.0], [3.0, 2.1, 4.0, 60.0], [0.5, 2.4, 1.0, 15.0]], sigma: Vec::new(), damp: Vec::new() };
    (law, vec![5.5, 2.0])
}

/// With induced dipoles, the force on a molecule and the torque on it are minus
/// the derivative of the energy with respect to moving it and turning it,
/// the dipoles minimised again each time: a small box (nearest neighbours at 3 A, the cut-off
/// at 4 A), so that the switch's region and the images are all in play.
#[test]
fn induced_forces_and_torques_are_the_derivative_of_the_energy() {
    use phys::liquid::Polarisable;
    let (base, alpha) = polarisable_law();
    let law = Polarisable { law: &base, alpha };
    let l = lattice(3, 3.0 * A, 300.0);
    let f = l.forces(&law);
    let plain = l.forces(&base);
    println!("  27 molecules: pair-law energy {:.6e}, with induction {:.6e}", plain.energy, f.energy);
    assert!(f.energy < plain.energy, "induction lowers the energy");
    assert!(f.dipoles.iter().any(|m| m.iter().any(|d| d.norm() > 1e-4)), "the dipoles are there");
    let h = 1e-5;
    let axes = [Vec3 { x: 1.0, y: 0.0, z: 0.0 }, Vec3 { x: 0.0, y: 1.0, z: 0.0 }, Vec3 { x: 0.0, y: 0.0, z: 1.0 }];
    let mut worst = 0.0f64;
    for m in [0usize, 13] {
        for (k, axis) in axes.iter().enumerate() {
            let energy_moved = |s: f64| {
                let mut c = l.clone();
                c.com[m] += axis.scale(s);
                c.forces(&law).energy
            };
            let numeric = -(energy_moved(h) - energy_moved(-h)) / (2.0 * h);
            let analytic = [f.force[m].x, f.force[m].y, f.force[m].z][k];
            worst = worst.max((numeric - analytic).abs() / (1e-6 + analytic.abs()));
            assert!((numeric - analytic).abs() < 1e-8 + 1e-6 * analytic.abs(), "molecule {m} force axis {k}: numeric {numeric:.9e} analytic {analytic:.9e}");
            let energy_turned = |s: f64| {
                let mut c = l.clone();
                c.orientation[m] = Quat::from_axis_angle(*axis, s).then(c.orientation[m]);
                c.forces(&law).energy
            };
            let numeric = -(energy_turned(h) - energy_turned(-h)) / (2.0 * h);
            let analytic = [f.torque[m].x, f.torque[m].y, f.torque[m].z][k];
            worst = worst.max((numeric - analytic).abs() / (1e-6 + analytic.abs()));
            assert!((numeric - analytic).abs() < 1e-8 + 1e-6 * analytic.abs(), "molecule {m} torque axis {k}: numeric {numeric:.9e} analytic {analytic:.9e}");
        }
    }
    println!("  forces and torques agree with differences of the energy; worst relative {worst:.1e}");
}

/// The molecular virial with induction is still minus the energy's slope under a
/// uniform scaling of the box and the centres.
#[test]
fn the_virial_with_induction_is_minus_the_energy_slope_under_scaling() {
    use phys::liquid::Polarisable;
    let (base, alpha) = polarisable_law();
    let law = Polarisable { law: &base, alpha };
    let l = lattice(3, 3.4 * A, 300.0);
    let f = l.forces(&law);
    let at = |lambda: f64| {
        let mut m = l.clone();
        for r in m.com.iter_mut() {
            *r = r.scale(lambda);
        }
        for c in m.cell.iter_mut() {
            *c *= lambda;
        }
        m.r_on *= 1.0;
        m.forces(&law).energy
    };
    let eps = 1e-5;
    let minus_slope = -(at(1.0 + eps) - at(1.0 - eps)) / (2.0 * eps);
    println!("  virial {:.8e}, minus dU/dlambda {minus_slope:.8e}", f.virial);
    assert!((f.virial - minus_slope).abs() < 1e-6 * f.virial.abs().max(1e-3), "{} against {minus_slope}", f.virial);
}

/// A polarisable liquid without a bath keeps its energy as the plain one does.
#[test]
fn a_polarisable_liquid_keeps_its_energy() {
    use phys::liquid::Polarisable;
    let (base, alpha) = polarisable_law();
    let law = Polarisable { law: &base, alpha };
    let mut l = lattice(3, 3.4 * A, 300.0);
    let total = |l: &Liquid, f: &Forces| {
        let (t, r) = l.kinetic();
        t + r + f.energy
    };
    let mut f = l.forces(&law);
    let e0 = total(&l, &f);
    let mut energies = Vec::new();
    for _ in 0..600 {
        f = l.step(&law, f, 10.0, None);
        energies.push(total(&l, &f));
    }
    let (lo, hi) = energies.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), e| (a.min(*e), b.max(*e)));
    let (kt, kr) = l.kinetic();
    println!("  27 polarisable molecules, 600 steps of 10 a.u.: energy {e0:.8} spread {:.2e} against kinetic {:.2e}", hi - lo, kt + kr);
    assert!(hi - lo < 2e-3 * (kt + kr), "the energy wandered by {:.2e}", hi - lo);
}

/// A law's text carries the polarisabilities of its types, and a text without
/// them gives none.
#[test]
fn a_laws_polarisabilities_are_read_from_its_text() {
    use phys::liquid::SiteSite;
    let text = "charges 1e0 5e-1 -2e0\npair 0 0 1e0 1e0 1e0 1e0\nbisector 2 2.67e-1\nalpha 0 5.5e0\nalpha 1 2e0\n";
    assert_eq!(SiteSite::alpha_from_text(text, 3), vec![5.5, 2.0, 0.0]);
    assert_eq!(SiteSite::alpha_from_text("charges 1e0\n", 2), vec![0.0, 0.0]);
    assert!(SiteSite::from_text(text).is_some(), "an alpha line does not stop the law reading");
}

/// Enough molecules that the induced-dipole solve runs across threads (more
/// than four hundred links): the answer is the same twice running to the last
/// bit, and one molecule's force and torque are still the derivative of the
/// energy.
#[test]
fn the_threaded_induction_solve_is_repeatable_and_correct() {
    use phys::liquid::Polarisable;
    let (base, alpha) = polarisable_law();
    let law = Polarisable { law: &base, alpha };
    let l = lattice(4, 3.0 * A, 300.0);
    let f1 = l.forces(&law);
    let f2 = l.forces(&law);
    assert_eq!(f1.energy.to_bits(), f2.energy.to_bits(), "two evaluations of one configuration differ");
    assert!(f1.force.iter().zip(&f2.force).all(|(a, b)| a.x.to_bits() == b.x.to_bits() && a.y.to_bits() == b.y.to_bits() && a.z.to_bits() == b.z.to_bits()));
    let h = 1e-5;
    let m = 21usize;
    for (k, axis) in [Vec3 { x: 1.0, y: 0.0, z: 0.0 }, Vec3 { x: 0.0, y: 0.0, z: 1.0 }].iter().enumerate() {
        let energy_moved = |s: f64| {
            let mut c = l.clone();
            c.com[m] += axis.scale(s);
            c.forces(&law).energy
        };
        let numeric = -(energy_moved(h) - energy_moved(-h)) / (2.0 * h);
        let analytic = if k == 0 { f1.force[m].x } else { f1.force[m].z };
        println!("  64 molecules, molecule {m}: force numeric {numeric:.9e} analytic {analytic:.9e}");
        assert!((numeric - analytic).abs() < 1e-7 + 1e-5 * analytic.abs(), "numeric {numeric:.9e} analytic {analytic:.9e}");
    }
}

#[test]
fn a_smeared_charge_is_a_point_charge_far_away_and_softer_near() {
    use phys::liquid::{SiteLaw, SiteSite};
    let zero = [0.0, 1.0, 0.0, 0.0];
    let point = SiteSite { charge: vec![0.6, -0.6], pair: vec![zero; 4], sigma: Vec::new(), damp: Vec::new() };
    let smeared = SiteSite { sigma: vec![0.7, 0.4], ..point.clone() };
    // Far apart the two agree; close, the smeared pair attracts less (a charge
    // inside another's cloud sees less of it); the derivative is the energy's.
    let (far_p, _) = point.site_pair(0, 1, 12.0);
    let (far_s, _) = smeared.site_pair(0, 1, 12.0);
    assert!((far_p - far_s).abs() < 1e-12 * far_p.abs().max(1.0), "{far_p} {far_s}");
    let (near_p, _) = point.site_pair(0, 1, 0.5);
    let (near_s, _) = smeared.site_pair(0, 1, 0.5);
    assert!(near_s.abs() < 0.8 * near_p.abs(), "{near_s} {near_p}");
    for r in [0.8, 1.5, 3.0, 6.0] {
        let h = 1e-6;
        let (_, d) = smeared.site_pair(0, 1, r);
        let fd = (smeared.site_pair(0, 1, r + h).0 - smeared.site_pair(0, 1, r - h).0) / (2.0 * h);
        assert!((d - fd).abs() < 1e-7 * d.abs().max(1.0), "r {r}: {d} against {fd}");
    }
}
