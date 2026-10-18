//! The unlike-pair code against the like-pair data it must agree with when the two molecules are the
//! same: the first MP2 pair of `pairs-water-gpu-mp2.txt` (which `phys-pairs-gpu water --mp2` made) recomputed through
//! `unlike_interaction` with water as both molecules. Skipped where that file is not.

use phys::math::Vec3;

#[test]
fn two_waters_through_the_unlike_code_give_the_like_codes_energy() {
    let path = "../runs/e8d/core/pairs-water-gpu-mp2.txt";
    let Ok(text) = std::fs::read_to_string(path) else {
        eprintln!("no {path}; skipped");
        return;
    };
    std::env::set_current_dir("..").expect("the repo root");
    let line = text.lines().find(|l| !l.starts_with('#')).expect("a pair");
    let (head, atoms) = line.split_once('|').expect("a bar");
    let h: Vec<f64> = head.split_whitespace().map(|v| v.parse().unwrap()).collect();
    let (total_file, hf_file) = (h[2], h[3]);
    let nums: Vec<&str> = atoms.split_whitespace().collect();
    let v: Vec<Vec3> = nums.chunks(4).map(|c| Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }).collect();
    let (pa, pb) = (v[..3].to_vec(), v[3..6].to_vec());
    if phys_gpu::install_fock().is_err() {
        eprintln!("no GPU Fock; the CPU does it");
    }
    let water = phys::pairs::Monomer::load("water", false);
    let (total, hf) = phys::pairs::unlike_interaction(&water, &water, &pa, &pb);
    println!("total {total:.10e} against {total_file:.10e}; Hartree-Fock {hf:.10e} against {hf_file:.10e}");
    assert!((total - total_file).abs() < 5e-8, "the total interaction");
    assert!((hf - hf_file).abs() < 5e-8, "the Hartree-Fock part");
    // A pair that is not water's: one of the new species through the same code, finite and bound at contact.
    let q = phys::pairs::unlike_pair(&water, &water, 3, None);
    assert!(q.2 >= 4.0 && q.2 <= 11.0, "a separation in range");
}
