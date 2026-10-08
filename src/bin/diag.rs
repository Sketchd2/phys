//! Where a pair law is wrong — `PLAY.md` E8b.
//!
//! ```sh
//! cargo run --release --bin phys-diag -- water law-total.txt [law-hf.txt] -mp2,-liq-gpu-mp2,-trimer-mp2
//! ```
//!
//! Reads the same `pairs-<name><tag>.txt` as `phys-fit` (column 3 the total
//! interaction energy, column 4 Hartree-Fock's), applies the law(s) with the
//! induction and bisector site they carry, and prints the error by the
//! separation of the two first atoms (the oxygens) — mean signed and RMS, in
//! kcal/mol, with the reference's mean — so that a shell the liquid lives in
//! and a fit's average do not hide one another. A second law is compared with
//! the Hartree-Fock column.

use phys::liquid::{extra_sites_from_text, pair_energy, with_bisector_site, with_frame_sites, ExtraSite, PairEnergy, Polarisable, SiteSite};
use phys::math::Vec3;

const KCAL: f64 = 627.509474;
const BOHR: f64 = 0.529177210903;

fn load(path: &str, types_atoms: usize) -> (SiteSite, Vec<f64>, Option<f64>, Vec<ExtraSite>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| panic!("no {path}"));
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    let _ = types_atoms;
    (law, alpha, SiteSite::bisector_from_text(&text).map(|b| b.1), extra_sites_from_text(&text))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args[0].clone();
    let tags = args.last().expect("tags").clone();
    let laws: Vec<&String> = args[1..args.len() - 1].iter().collect();
    let text: String = tags.split(',').map(|t| std::fs::read_to_string(format!("pairs-{name}{t}.txt")).unwrap_or_else(|_| panic!("no pairs-{name}{t}.txt"))).collect::<Vec<_>>().join("\n");
    let mut rows: Vec<(f64, f64, f64, Vec<(Vec3, usize)>)> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()) {
        let (head, tail) = line.split_once('|').expect("a |");
        let h: Vec<f64> = head.split_whitespace().map(|x| x.parse().expect("number")).collect();
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<(Vec3, usize)> = t.chunks(4).map(|c| (Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }, c[3].parse().unwrap())).collect();
        let n = atoms.len() / 2;
        let roo = (atoms[0].0 - atoms[n].0).norm();
        rows.push((roo, h[2], h[3], atoms));
    }
    let n_atoms = rows[0].3.len() / 2;
    let edges = [0.0, 2.5, 2.8, 3.1, 3.4, 3.8, 4.4, 5.5, 99.0];
    for (li, lp) in laws.iter().enumerate() {
        let (law, alpha, bis, frame) = load(lp, n_atoms);
        let types = law.charge.len();
        let atom_types = types - bis.is_some() as usize - !frame.is_empty() as usize;
        println!("law {lp} against the {} column", if li == 0 { "total (MP2)" } else { "Hartree-Fock" });
        println!("  O-O (A)        n   ref mean   mean err   rms err   (kcal/mol)");
        for w in edges.windows(2) {
            let sel: Vec<&(f64, f64, f64, Vec<(Vec3, usize)>)> = rows.iter().filter(|r| r.0 * BOHR >= w[0] && r.0 * BOHR < w[1]).collect();
            if sel.is_empty() {
                continue;
            }
            let (mut s, mut s2, mut r) = (0.0, 0.0, 0.0);
            for row in &sel {
                let (a, b) = (&row.3[..n_atoms], &row.3[n_atoms..]);
                let (a, b) = match bis {
                    Some(d) => (with_frame_sites(&with_bisector_site(a, atom_types, d), &frame), with_frame_sites(&with_bisector_site(b, atom_types, d), &frame)),
                    None => (a.to_vec(), b.to_vec()),
                };
                let reference = if li == 0 { row.1 } else { row.2 };
                let pe = PairEnergy { a, b, energy: reference };
                let e = pair_energy(&Polarisable { law: &law, alpha: alpha.clone() }, &pe);
                let d = (e - reference) * KCAL;
                s += d;
                s2 += d * d;
                r += reference * KCAL;
            }
            let n = sel.len() as f64;
            println!("  {:>4.1} - {:>4.1}  {:>4}  {:>9.3}  {:>9.3}  {:>9.3}", w[0], w[1].min(9.9), sel.len(), r / n, s / n, (s2 / n).sqrt());
        }
    }
}
