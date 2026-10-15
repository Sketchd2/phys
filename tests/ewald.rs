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
