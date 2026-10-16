//! A law's energy for each pair in a pair file, next to the reference's
//! (`docs/PLAY.md` E8d): one line per pair, `index O-O(A) E_reference E_law` in
//! hartree, for fitting a term to what a law leaves over.
//!
//! ```sh
//! phys-lawenergy law.txt -gpu-wall-mp2        # pairs-water<tag>.txt in the working directory
//! ```

use phys::liquid::{extra_sites_from_text, pair_energy, with_bisector_site, with_frame_sites, PairEnergy, Polarisable, SiteSite};
use phys::math::Vec3;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let text = std::fs::read_to_string(&args[0]).unwrap_or_else(|_| panic!("no {}", args[0]));
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    let bis = SiteSite::bisector_from_text(&text).map(|b| b.1);
    let frame = extra_sites_from_text(&text);
    let mut distinct: Vec<usize> = frame.iter().map(|s| s.ty).collect();
    distinct.sort();
    distinct.dedup();
    let atom_types = law.charge.len() - bis.is_some() as usize - distinct.len();
    let polar = Polarisable { law: &law, alpha };
    let pairs = std::fs::read_to_string(format!("pairs-water{}.txt", args[1])).unwrap_or_else(|_| panic!("no pairs-water{}.txt", args[1]));
    for line in pairs.lines().filter(|l| !l.starts_with('#') && l.contains('|')) {
        let (head, tail) = line.split_once('|').unwrap();
        let h: Vec<&str> = head.split_whitespace().collect();
        let (k, e_ref): (&str, f64) = (h[0], h[2].parse().unwrap());
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<(Vec3, usize)> = t.chunks(4).map(|c| (Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }, c[3].parse().unwrap())).collect();
        let n = atoms.len() / 2;
        let place = |m: &[(Vec3, usize)]| with_frame_sites(&with_bisector_site(m, atom_types, bis.unwrap_or(0.0)), &frame);
        let e_law = pair_energy(&polar, &PairEnergy { a: place(&atoms[..n]), b: place(&atoms[n..]), energy: 0.0 });
        println!("{k} {:.4} {:.10e} {:.10e}", (atoms[0].0 - atoms[n].0).norm() * 0.529177210903, e_ref, e_law);
    }
}
