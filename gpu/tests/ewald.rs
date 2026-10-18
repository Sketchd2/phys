//! The GPU's reciprocal sums against the reference's: the same structure factors, fields,
//! response, forces and arm sums for 216 molecules' worth of charges and dipoles, to far
//! inside any tolerance a force evaluation needs, and the same bits twice.

use phys::ewald::{begin_reciprocal, set_reciprocal_engine, CpuReciprocal, Ewald, ReciprocalSession, C};
use phys::math::Vec3;

fn system() -> (Ewald, Vec<Vec3>, Vec<f64>, Vec<Vec3>, Vec<Vec3>, Vec<Vec3>, Vec<Vec3>) {
    let cell = [35.2; 3];
    let ew = Ewald::new(4.2, cell, 1e-5);
    let mut k = 12345u64;
    let mut rnd = move || {
        k = k.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((k >> 11) as f64) / ((1u64 << 53) as f64)
    };
    let mut cpos = Vec::new();
    let mut cq = Vec::new();
    let mut carm = Vec::new();
    let mut dpos = Vec::new();
    let mut darm = Vec::new();
    let mut mu = Vec::new();
    for _ in 0..216 {
        let com = Vec3 { x: rnd() * cell[0], y: rnd() * cell[1], z: rnd() * cell[2] };
        for s in 0..3 {
            let arm = Vec3 { x: rnd() * 3.0 - 1.5, y: rnd() * 3.0 - 1.5, z: rnd() * 3.0 - 1.5 };
            cpos.push(com + arm);
            carm.push(arm);
            cq.push(if s == 0 { -1.1 } else { 0.55 });
            let arm = Vec3 { x: rnd() * 3.0 - 1.5, y: rnd() * 3.0 - 1.5, z: rnd() * 3.0 - 1.5 };
            dpos.push(com + arm);
            darm.push(arm);
            mu.push(Vec3 { x: rnd() * 0.6 - 0.3, y: rnd() * 0.6 - 0.3, z: rnd() * 0.6 - 0.3 });
        }
    }
    (ew, cpos, cq, carm, dpos, darm, mu)
}

fn worst_c(a: &[C], b: &[C]) -> f64 {
    let scale = b.iter().map(|c| c.re.abs().max(c.im.abs())).fold(0.0, f64::max);
    a.iter().zip(b).map(|(x, y)| (x.re - y.re).abs().max((x.im - y.im).abs())).fold(0.0, f64::max) / scale
}

fn worst_v(a: &[Vec3], b: &[Vec3]) -> f64 {
    let scale = b.iter().map(|v| v.x.abs().max(v.y.abs()).max(v.z.abs())).fold(0.0, f64::max);
    a.iter().zip(b).map(|(x, y)| (x.x - y.x).abs().max((x.y - y.y).abs()).max((x.z - y.z).abs())).fold(0.0, f64::max) / scale
}

#[test]
fn the_gpu_sums_are_the_reference_and_are_repeatable() {
    let (ew, cpos, cq, carm, dpos, darm, mu) = system();
    set_reciprocal_engine(None);
    let mut cpu = CpuReciprocal::new(&ew, &cpos, &cq, &carm, &dpos, &darm);
    let sq0 = cpu.charge_sum();
    let e00 = cpu.field_at_dipoles();
    let r0 = cpu.dipole_response(&mu);
    let f0 = cpu.finish(&mu);

    let name = match phys_gpu::install_ewald() {
        Ok(n) => n,
        Err(e) => {
            eprintln!("no GPU here ({e}); skipped");
            return;
        }
    };
    println!("{name}, {} wave vectors", ew.nk());
    let mut run = || {
        let mut s = begin_reciprocal(&ew, &cpos, &cq, &carm, &dpos, &darm);
        let sq = s.charge_sum();
        let e0 = s.field_at_dipoles();
        let r = s.dipole_response(&mu);
        let f = s.finish(&mu);
        (sq, e0, r, f)
    };
    let (sq, e0, r, f) = run();
    let again = run();
    let mut report = |what: &str, v: f64| {
        println!("{what:<28} {v:.2e}");
        assert!(v < 1e-12, "{what}");
    };
    report("charges' structure factor", worst_c(&sq, &sq0));
    report("field at the dipoles", worst_v(&e0, &e00));
    report("dipole response", worst_v(&r, &r0));
    report("dipoles' structure factor", worst_c(&f.d, &f0.d));
    report("charge arm sum", worst_c(&f.tc, &f0.tc));
    report("dipole arm sum", worst_c(&f.td, &f0.td));
    report("force on the charges", worst_v(&f.grad_charge, &f0.grad_charge));
    report("force on the dipoles' sites", worst_v(&f.grad_dipole_pos, &f0.grad_dipole_pos));
    // The scalar loops are the reference's own, given the card's sums.
    let stot = |sq: &[C], d: &[C]| -> Vec<C> { sq.iter().zip(d).map(|(a, b)| C { re: a.re + b.re, im: a.im + b.im }).collect() };
    let (v0, v1) = (ew.virial_from_sums(&stot(&sq0, &f0.d), &f0.tc, &f0.td, &f0.d), ew.virial_from_sums(&stot(&sq, &f.d), &f.tc, &f.td, &f.d));
    println!("virial {v1:.12e} against {v0:.12e}");
    assert!((v1 - v0).abs() / v0.abs() < 1e-11, "virial");
    let (e1, e0e) = (ew.energy(&sq), ew.energy(&sq0));
    assert!((e1 - e0e).abs() / e0e.abs() < 1e-11, "energy");
    // Two calls, the same bits.
    let same = sq.iter().zip(&again.0).all(|(a, b)| a.re.to_bits() == b.re.to_bits() && a.im.to_bits() == b.im.to_bits())
        && r.iter().zip(&again.2).all(|(a, b)| a.x.to_bits() == b.x.to_bits() && a.y.to_bits() == b.y.to_bits() && a.z.to_bits() == b.z.to_bits())
        && f.grad_charge.iter().zip(&again.3.grad_charge).all(|(a, b)| a.x.to_bits() == b.x.to_bits());
    assert!(same, "two calls gave different bits");
    set_reciprocal_engine(None);
}

#[test]
fn where_a_submission_goes() {
    let (ew, cpos, cq, carm, dpos, darm, mu) = system();
    let gpu = match phys_gpu::GpuEwald::new() {
        Ok(g) => std::sync::Arc::new(g),
        Err(e) => {
            eprintln!("no GPU here ({e}); skipped");
            return;
        }
    };
    set_reciprocal_engine(Some(gpu.clone()));
    let mut s = begin_reciprocal(&ew, &cpos, &cq, &carm, &dpos, &darm);
    let _ = s.charge_sum();
    let _ = s.dipole_response(&mu);
    let (nk, nd, nc) = (ew.nk(), dpos.len(), cpos.len());
    let reps = 50;
    println!("an empty-ish submission (k_stot over {} groups): {:.3} ms", nk.div_ceil(64), gpu.bench(&[("k_stot", nk.div_ceil(64))], reps));
    println!("k_dsum ({nk} workgroups):                    {:.3} ms", gpu.bench(&[("k_dsum", nk)], reps));
    println!("k_resp ({nd} workgroups):                    {:.3} ms", gpu.bench(&[("k_resp", nd)], reps));
    println!("k_sq   ({nk} workgroups):                    {:.3} ms", gpu.bench(&[("k_sq", nk)], reps));
    println!("k_gch  ({nc} workgroups):                    {:.3} ms", gpu.bench(&[("k_gch", nc)], reps));
    println!("k_arm_c, k_arm_d, k_gdpos:                   {:.3} ms", gpu.bench(&[("k_arm_c", nk), ("k_arm_d", nk), ("k_gdpos", nd)], reps));
    set_reciprocal_engine(None);
}
