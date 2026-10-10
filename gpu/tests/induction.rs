//! The GPU's induced-dipole solve against the reference: the same energy,
//! dipoles, forces and virial for a liquid of 216 molecules, to far inside the
//! solve's tolerance, and the same bits twice.

use phys::flex::{Flex, FlexSites, Surface};
use phys::liquid::{Polarisable, SiteSite};
use phys::math::{Quat, Vec3};

const AMU: f64 = 1822.888486;

fn liquid(side_count: usize) -> (Flex, SiteSite, Vec<f64>) {
    let text = std::fs::read_to_string("../law-water-mp2-cloud-final.txt").expect("the cloud law");
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    let sites = FlexSites::from_text(&text, [0, 1, 1]);
    let surface = Surface::from_text(&std::fs::read_to_string("../surface-water.txt").expect("the surface"));
    let r = 0.9789 / 0.529177210903;
    let th = 104.08f64.to_radians();
    let template = [Vec3::ZERO, Vec3 { x: r, y: 0.0, z: 0.0 }, Vec3 { x: r * th.cos(), y: r * th.sin(), z: 0.0 }];
    let spacing = 5.87;
    let side = spacing * side_count as f64;
    let mut k = 987654321u64;
    let mut rnd = || {
        k = k.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((k >> 33) as f64) / ((1u64 << 31) as f64)
    };
    let mut pos = Vec::new();
    for ix in 0..side_count {
        for iy in 0..side_count {
            for iz in 0..side_count {
                let centre = Vec3 { x: (ix as f64 + 0.5) * spacing, y: (iy as f64 + 0.5) * spacing, z: (iz as f64 + 0.5) * spacing };
                let q = Quat { w: rnd() * 2.0 - 1.0, v: Vec3 { x: rnd() * 2.0 - 1.0, y: rnd() * 2.0 - 1.0, z: rnd() * 2.0 - 1.0 } };
                let n = (q.w * q.w + q.v.norm2()).sqrt();
                let q = Quat { w: q.w / n, v: q.v.scale(1.0 / n) };
                let mut m = [Vec3::ZERO; 3];
                for a in 0..3 {
                    m[a] = centre + q.rotate(template[a]) + Vec3 { x: rnd() - 0.5, y: rnd() - 0.5, z: rnd() - 0.5 }.scale(0.2);
                }
                pos.push(m);
            }
        }
    }
    let n = pos.len();
    let flex = Flex {
        sites,
        masses: [15.9949146 * AMU, 1.00782503 * AMU, 1.00782503 * AMU],
        cell: [side; 3],
        pos,
        vel: vec![[Vec3::ZERO; 3]; n],
        r_on: (0.49 * side).min(16.0) - 2.0,
        r_cut: (0.49 * side).min(16.0),
        surface,
        bond: None,
        geo: None,
    };
    (flex, law, alpha)
}

#[test]
fn the_gpu_solve_is_the_reference_and_is_repeatable() {
    let (flex, law, alpha) = liquid(6);
    let pol = Polarisable { law: &law, alpha: alpha.clone() };
    phys::induction::set_induction_engine(None);
    let reference = flex.forces(&pol, None);
    let name = match phys_gpu::install_induction() {
        Ok(n) => n,
        Err(e) => {
            eprintln!("no GPU here ({e}); skipped");
            return;
        }
    };
    println!("{name}");
    let gpu = flex.forces(&pol, None);
    let again = flex.forces(&pol, None);
    let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-12);
    println!("energy {:.12} against {:.12} (relative {:.2e}); virial {:.9} against {:.9} (relative {:.2e})", gpu.inter, reference.inter, rel(gpu.inter, reference.inter), gpu.virial, reference.virial, rel(gpu.virial, reference.virial));
    assert!(rel(gpu.inter, reference.inter) < 1e-9, "energy");
    assert!(rel(gpu.virial, reference.virial) < 1e-8, "virial");
    let mut worst: f64 = 0.0;
    let mut scale: f64 = 0.0;
    let mut worst_dipole: f64 = 0.0;
    for (g, r) in gpu.force.iter().zip(&reference.force) {
        for a in 0..3 {
            worst = worst.max((g[a] - r[a]).norm());
            scale = scale.max(r[a].norm());
        }
    }
    for (g, r) in gpu.dipoles.iter().zip(&reference.dipoles) {
        for (a, b) in g.iter().zip(r) {
            worst_dipole = worst_dipole.max((*a - *b).norm());
        }
    }
    println!("largest force difference {worst:.2e} against a largest force {scale:.2e}; dipoles differ by at most {worst_dipole:.2e}");
    assert!(worst < 1e-8 * scale.max(1e-3), "forces");
    assert!(worst_dipole < 1e-7, "dipoles");
    // The same bits twice.
    assert_eq!(gpu.inter.to_bits(), again.inter.to_bits(), "energy repeats");
    assert_eq!(gpu.virial.to_bits(), again.virial.to_bits(), "virial repeats");
    for (g, r) in gpu.force.iter().zip(&again.force) {
        for a in 0..3 {
            assert_eq!((g[a].x.to_bits(), g[a].y.to_bits(), g[a].z.to_bits()), (r[a].x.to_bits(), r[a].y.to_bits(), r[a].z.to_bits()));
        }
    }
    phys::induction::set_induction_engine(None);
}

/// The whole force routine of the flexible liquid with both engines against none, for a law
/// with fixed charges and for a cloud that follows each molecule's angle.
#[test]
fn the_gpu_pair_loop_and_solve_give_the_cpu_forces_with_and_without_a_geometry_following_cloud() {
    for follow in [false, true] {
        let (mut flex, law, alpha) = liquid(6);
        if follow {
            let text = std::fs::read_to_string("../law-water-mp2-cloud-final.txt").unwrap();
            let bis = SiteSite::bisector_from_text(&text).unwrap().1;
            let mut g = phys::geolaw::GeoCloud::constant([law.charge[0], law.charge[1]], phys::geolaw::x_of_law(&law, bis), 104.08f64.to_radians());
            let c1 = [0.15, -0.10, 0.20, 0.08, 0.12, -0.10];
            let c2 = [0.05, 0.04, -0.06, 0.02, -0.05, 0.03];
            for p in 0..6 {
                g.coef[p][1] = c1[p];
                g.coef[p][2] = c2[p];
            }
            flex.geo = Some(g);
            // Angles that differ from molecule to molecule.
            for (i, m) in flex.pos.iter_mut().enumerate() {
                let shift = ((i as f64 * 0.7).sin()) * 0.12;
                let (o, h2) = (m[0], m[2]);
                let n = (m[2] - o).unit().cross((m[1] - o).unit()).unit();
                let rot = Quat { w: (shift / 2.0).cos(), v: n.scale((shift / 2.0).sin()) };
                m[2] = o + rot.rotate(h2 - o);
            }
        }
        let pol = Polarisable { law: &law, alpha: alpha.clone() };
        phys::induction::set_induction_engine(None);
        phys::liquid::set_pair_engine(None);
        let cpu = flex.forces(&pol, None);
        phys_gpu::install_induction().expect("a GPU");
        phys_gpu::install_pairs().expect("a GPU");
        let gpu = flex.forces(&pol, None);
        let again = flex.forces(&pol, None);
        phys::induction::set_induction_engine(None);
        phys::liquid::set_pair_engine(None);
        let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-12);
        let mut worst: f64 = 0.0;
        let mut scale: f64 = 0.0;
        for (g, c) in gpu.force.iter().zip(&cpu.force) {
            for a in 0..3 {
                worst = worst.max((g[a] - c[a]).norm());
                scale = scale.max(c[a].norm());
            }
        }
        println!("follow {follow}: energy {:.12} against {:.12} (relative {:.2e}), virial relative {:.2e}, largest force difference {worst:.2e} against {scale:.2e}", gpu.inter, cpu.inter, rel(gpu.inter, cpu.inter), rel(gpu.virial, cpu.virial));
        assert!(rel(gpu.inter, cpu.inter) < 1e-10, "energy");
        assert!(rel(gpu.virial, cpu.virial) < 1e-9, "virial");
        assert!(worst < 1e-9 * scale, "forces: {worst:e} against {scale:e}");
        assert_eq!(gpu.inter.to_bits(), again.inter.to_bits(), "repeatable");
        for (g, r) in gpu.force.iter().zip(&again.force) {
            for a in 0..3 {
                assert_eq!((g[a].x.to_bits(), g[a].y.to_bits(), g[a].z.to_bits()), (r[a].x.to_bits(), r[a].y.to_bits(), r[a].z.to_bits()));
            }
        }
    }
}
