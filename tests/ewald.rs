//! The Ewald sum (`src/ewald.rs`, `PLAY.md` E8d): what makes it right is that the
//! answer does not depend on how the interaction was split, that its forces and
//! torques are the derivative of its energy, and that its virial is the
//! derivative of the energy under rigid scaling of the box. Each is checked on a
//! small liquid of polarisable water in random orientations.

use phys::electrons::grow::{equivalent_atoms, Resume};
use phys::electrons::molecule::Molecule;
use phys::ewald::Ewald;
use phys::liquid::{extra_sites_from_text, EwaldLaw, Kind, Liquid, Polarisable, SiteLaw, SiteSite, AMU};
use phys::math::{Quat, Vec3};
use phys::rng::{Purpose, Stream};

const SIDE: f64 = 23.5;
const N: usize = 64;

struct Setup {
    liquid: Liquid,
    law: SiteSite,
    alpha: Vec<f64>,
}

fn setup() -> Setup {
    let text = std::fs::read_to_string("law-water-mp2-cloud-final.txt").expect("the cloud law");
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    let bisector = SiteSite::bisector_from_text(&text);
    let state = Resume::from_text(&std::fs::read_to_string("grow-water.state").expect("a state")).expect("readable");
    let z = vec![8u32, 1, 1];
    let mono = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let types = equivalent_atoms(&mono);
    let frame = extra_sites_from_text(&text);
    let kind = Kind::of_molecule_sites(&z, &state.positions, &types, bisector, &frame);
    let molecule_g = kind.mass / AMU * 1.66053906660e-24;
    let per_bohr3 = 1.0 / molecule_g * 0.529177210903e-8f64.powi(3);
    let _ = per_bohr3;
    let mut liquid = Liquid::slab(kind, N, N as f64 / SIDE.powi(3), SIDE, SIDE, 298.0, 7, 9.5, 11.5);
    // Random orientations and a jitter of the lattice, seeded.
    let mut s = Stream::at(0xe1a1d, 0, 0, Purpose::Positions);
    for i in 0..N {
        let (u1, u2, u3) = (s.uniform(), s.uniform(), s.uniform());
        let tau = std::f64::consts::TAU;
        let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
        liquid.orientation[i] = Quat { w: a * (tau * u2).sin(), v: Vec3 { x: a * (tau * u2).cos(), y: b * (tau * u3).sin(), z: b * (tau * u3).cos() } };
        liquid.com[i] += Vec3 { x: s.uniform() - 0.5, y: s.uniform() - 0.5, z: s.uniform() - 0.5 }.scale(0.8);
    }
    Setup { liquid, law, alpha }
}

fn ewald_law<'a>(inner: &'a dyn SiteLaw, sp: f64, cell: [f64; 3]) -> EwaldLaw<'a> {
    EwaldLaw { inner, ewald: Ewald::new(sp, cell, 1e-7) }
}

/// The energy does not depend on the splitting width, and the truncated sum is
/// visibly a different thing (so this can fail).
#[test]
fn the_energy_does_not_depend_on_how_it_was_split() {
    let s = setup();
    let polar = Polarisable { law: &s.law, alpha: s.alpha.clone() };
    let cell = s.liquid.cell;
    let truncated = s.liquid.forces(&polar).energy;
    let e: Vec<f64> = [2.8, 2.9, 3.0].iter().map(|&sp| s.liquid.forces(&ewald_law(&polar, sp, cell)).energy).collect();
    let spread = e.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - e.iter().cloned().fold(f64::INFINITY, f64::min);
    println!("Ewald energies {e:?}, spread {spread:.3e}; truncated {truncated:.6}");
    assert!(spread < 3e-4, "the energy moved by {spread:e} hartree with the splitting width");
    assert!((truncated - e[1]).abs() > 20.0 * spread.max(1e-5), "the truncated energy {truncated} is not distinguishable from Ewald {}", e[1]);
}

/// Force and torque are minus the derivative of the energy (the dipoles
/// re-minimised at each point).
#[test]
fn forces_and_torques_are_the_energys_derivative() {
    let s = setup();
    let polar = Polarisable { law: &s.law, alpha: s.alpha.clone() };
    let law = ewald_law(&polar, 3.0, s.liquid.cell);
    let f0 = s.liquid.forces(&law);
    for &i in &[2usize, 41] {
        for axis in 0..3 {
            let d = 2e-3;
            let e_pos = |dx: f64| -> f64 {
                let mut l = s.liquid.clone();
                let mut v = Vec3::ZERO;
                match axis {
                    0 => v.x = dx,
                    1 => v.y = dx,
                    _ => v.z = dx,
                }
                l.com[i] += v;
                l.forces(&law).energy
            };
            let fd = -(e_pos(d) - e_pos(-d)) / (2.0 * d);
            let fa = [f0.force[i].x, f0.force[i].y, f0.force[i].z][axis];
            assert!((fa - fd).abs() < 2e-6, "molecule {i} axis {axis}: force {fa:e} against {fd:e}");
            let e_rot = |th: f64| -> f64 {
                let mut l = s.liquid.clone();
                let e = [Vec3 { x: 1.0, y: 0.0, z: 0.0 }, Vec3 { x: 0.0, y: 1.0, z: 0.0 }, Vec3 { x: 0.0, y: 0.0, z: 1.0 }][axis];
                l.orientation[i] = l.orientation[i].then(Quat::from_axis_angle(e, th)).unit();
                l.forces(&law).energy
            };
            let td = -(e_rot(d) - e_rot(-d)) / (2.0 * d);
            let tb = s.liquid.orientation[i].conjugate().rotate(f0.torque[i]);
            let ta = [tb.x, tb.y, tb.z][axis];
            assert!((ta - td).abs() < 2e-6, "molecule {i} axis {axis}: torque {ta:e} against {td:e}");
        }
    }
}

/// The molecular virial is minus the derivative of the energy under rigid
/// scaling of the box and the centres.
#[test]
fn the_virial_is_the_energys_derivative_under_scaling() {
    let s = setup();
    let polar = Polarisable { law: &s.law, alpha: s.alpha.clone() };
    let law = ewald_law(&polar, 3.0, s.liquid.cell);
    let w = s.liquid.forces(&law).virial;
    let e_at = |lam: f64| -> f64 {
        let mut l = s.liquid.clone();
        l.cell = [s.liquid.cell[0] * lam, s.liquid.cell[1] * lam, s.liquid.cell[2] * lam];
        for c in l.com.iter_mut() {
            *c = c.scale(lam);
        }
        let scaled = ewald_law(&polar, 3.0, l.cell);
        l.forces(&scaled).energy
    };
    let eps = 2e-4;
    let w_fd = -(e_at(1.0 + eps) - e_at(1.0 - eps)) / ((1.0 + eps).ln() - (1.0 - eps).ln());
    assert!((w - w_fd).abs() < 1e-4 * w.abs().max(1.0), "virial {w:e} against {w_fd:e}");
}

/// The table-based `erf_fast` the Ewald paths use agrees with the series to rounding, everywhere
/// including the sign, the small arguments and the saturated tail.
#[test]
fn the_fast_erf_agrees_with_the_series() {
    let mut worst: f64 = 0.0;
    let mut at = 0.0;
    for i in -6000..=6000 {
        let x = i as f64 / 1000.0 + 1.37e-4;
        let d = (phys::induction::erf_fast(x) - phys::induction::erf(x)).abs();
        // The series sums ~100 terms of e^(x^2) size beyond x = 4 and carries ~1e-14 of its own there.
        let limit = if x.abs() < 4.0 { 2e-15 } else { 3e-14 };
        if d / limit > worst {
            worst = d / limit;
            at = x;
        }
    }
    assert!(worst < 1.0, "erf_fast differs from erf by {worst} of its allowance at {at}");
}

// ---- a liquid of two kinds of molecule (docs/SOLUTIONS.md 4.1) ---------------------------------

/// The law with every site type given a second copy (types `T..2T`), the copies' own charges and
/// polarisabilities scaled, and every pair entry between any two copies the original's.
fn doubled(law: &SiteSite, alpha: &[f64], charge_scale: f64, alpha_scale: f64) -> (SiteSite, Vec<f64>) {
    let t = law.charge.len();
    let pick = |v: &[f64], a: usize, b: usize| v.get((a % t) * t + (b % t)).copied().unwrap_or(0.0);
    let mut out = law.clone();
    out.charge = (0..2 * t).map(|a| law.charge[a % t] * if a >= t { charge_scale } else { 1.0 }).collect();
    out.sigma = if law.sigma.is_empty() { Vec::new() } else { (0..2 * t).map(|a| law.sigma[a % t]).collect() };
    out.pair = (0..2 * t * 2 * t).map(|k| law.pair[((k / (2 * t)) % t) * t + (k % (2 * t)) % t]).collect();
    out.damp = if law.damp.is_empty() { Vec::new() } else { (0..2 * t * 2 * t).map(|k| pick(&law.damp, k / (2 * t), k % (2 * t))).collect() };
    out.overlap = if law.overlap.is_empty() { Vec::new() } else { (0..2 * t * 2 * t).map(|k| pick(&law.overlap, k / (2 * t), k % (2 * t))).collect() };
    let a2: Vec<f64> = (0..2 * t).map(|a| alpha[a % t] * if a >= t { alpha_scale } else { 1.0 }).collect();
    (out, a2)
}

/// A kind with its site types moved up by `shift`.
fn shifted(kind: &Kind, shift: usize) -> Kind {
    let mut k = kind.clone();
    k.types = k.types.iter().map(|t| t + shift).collect();
    k
}

fn mixed_setup(charge_scale: f64, alpha_scale: f64) -> (Setup, Liquid, SiteSite, Vec<f64>) {
    let pure = setup();
    let t = pure.law.charge.len();
    let (law2, alpha2) = doubled(&pure.law, &pure.alpha, charge_scale, alpha_scale);
    let (k0, k1) = (pure.liquid.kind.clone(), shifted(&pure.liquid.kind, t));
    let n = pure.liquid.com.len();
    let mut mix = Liquid::mixture(vec![(k0, n / 2), (k1, n - n / 2)], n as f64 / SIDE.powi(3), SIDE, SIDE, 298.0, 7, 9.5, 11.5);
    // The same geometry as the pure liquid, so only the labels differ.
    mix.com = pure.liquid.com.clone();
    mix.orientation = pure.liquid.orientation.clone();
    (pure, mix, law2, alpha2)
}

/// Two labels for one kind of molecule, with every unlike pair given the like pair's entries, are one
/// liquid: the mixture's energy, forces, torques and virial are the pure liquid's. A species
/// bookkeeping slip (the wrong charges, the wrong polarisable sites, a constant taken off the wrong
/// molecule) would move them.
#[test]
fn two_labels_of_one_molecule_are_one_liquid() {
    let (pure, mix, law2, alpha2) = mixed_setup(1.0, 1.0);
    assert!(!mix.homogeneous(), "the mixture has two species");
    let polar = Polarisable { law: &pure.law, alpha: pure.alpha.clone() };
    let polar2 = Polarisable { law: &law2, alpha: alpha2 };
    let (a, b) = (pure.liquid.forces(&ewald_law(&polar, 3.0, pure.liquid.cell)), mix.forces(&ewald_law(&polar2, 3.0, mix.cell)));
    let rel = |x: f64, y: f64| (x - y).abs() / y.abs().max(1e-12);
    println!("energy {:.10} against {:.10}; virial {:.10} against {:.10}", b.energy, a.energy, b.virial, a.virial);
    assert!(rel(b.energy, a.energy) < 1e-9, "energy");
    assert!(rel(b.virial, a.virial) < 1e-8, "virial");
    let worst = (0..a.force.len()).map(|i| (a.force[i] - b.force[i]).norm().max((a.torque[i] - b.torque[i]).norm())).fold(0.0, f64::max);
    let scale = a.force.iter().map(|f| f.norm()).fold(0.0, f64::max);
    println!("worst force or torque difference {worst:.2e} against the largest force {scale:.2e}");
    assert!(worst < 1e-8 * scale, "forces or torques moved");
    // And it is not trivially so: the same mixture with the second kind's charges scaled is different.
    let (_, mix3, law3, alpha3) = mixed_setup(0.9, 1.0);
    let polar3 = Polarisable { law: &law3, alpha: alpha3 };
    let c = mix3.forces(&ewald_law(&polar3, 3.0, mix3.cell));
    assert!(rel(c.energy, a.energy) > 1e-3, "a changed charge on one kind must change the energy");
}

/// In a mixture whose kinds differ (charges and polarisabilities), force and torque are still minus
/// the derivative of the energy, for a molecule of each kind.
#[test]
fn a_mixtures_forces_are_the_energys_derivative() {
    let (_, mix, law2, alpha2) = mixed_setup(0.9, 1.15);
    let polar = Polarisable { law: &law2, alpha: alpha2 };
    let law = ewald_law(&polar, 3.0, mix.cell);
    let f0 = mix.forces(&law);
    let first0 = (0..mix.com.len()).find(|&i| mix.species_of(i) == 0).expect("a molecule of the first kind");
    let first1 = (0..mix.com.len()).find(|&i| mix.species_of(i) == 1).expect("a molecule of the second kind");
    for &i in &[first0, first1] {
        for axis in 0..3 {
            let d = 2e-3;
            let e_pos = |dx: f64| -> f64 {
                let mut l = mix.clone();
                let mut v = Vec3::ZERO;
                match axis {
                    0 => v.x = dx,
                    1 => v.y = dx,
                    _ => v.z = dx,
                }
                l.com[i] += v;
                l.forces(&law).energy
            };
            let fd = -(e_pos(d) - e_pos(-d)) / (2.0 * d);
            let fa = [f0.force[i].x, f0.force[i].y, f0.force[i].z][axis];
            assert!((fa - fd).abs() < 2e-6, "molecule {i} (kind {}) axis {axis}: force {fa:e} against {fd:e}", mix.species_of(i));
            let e_rot = |th: f64| -> f64 {
                let mut l = mix.clone();
                let e = [Vec3 { x: 1.0, y: 0.0, z: 0.0 }, Vec3 { x: 0.0, y: 1.0, z: 0.0 }, Vec3 { x: 0.0, y: 0.0, z: 1.0 }][axis];
                l.orientation[i] = l.orientation[i].then(Quat::from_axis_angle(e, th)).unit();
                l.forces(&law).energy
            };
            let td = -(e_rot(d) - e_rot(-d)) / (2.0 * d);
            let tb = mix.orientation[i].conjugate().rotate(f0.torque[i]);
            let ta = [tb.x, tb.y, tb.z][axis];
            assert!((ta - td).abs() < 2e-6, "molecule {i} (kind {}) axis {axis}: torque {ta:e} against {td:e}", mix.species_of(i));
        }
    }
}
