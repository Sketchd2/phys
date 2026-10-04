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
    let law = SiteSite { charge: vec![-0.7, 0.35], pair: vec![[40.0, 1.9, 15.0, 300.0], [3.0, 2.1, 4.0, 60.0], [3.0, 2.1, 4.0, 60.0], [0.5, 2.4, 1.0, 15.0]] };
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
