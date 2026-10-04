//! Route A's law from its pair energies — `PLAY.md` E8.
//!
//! ```sh
//! cargo run --release --bin phys-fit -- water [temperature for the weights, K] [tag]
//! ```
//!
//! Reads `pairs-<name><tag>.txt` (written by `phys-pairs`, or by
//! `phys-pairs-gpu` with tag `-gpu`), holds every fifth pair
//! out of the fit, fits the site-site law to the rest for each exchange
//! partner the pairs carry — PBE exchange (nothing fitted in it) and revPBE
//! (published vdW-DF1) — and reports the residual on the pairs it fitted and,
//! the one that is the law's error by the owner's condition, on the pairs it
//! never saw. Writes `law-<name><tag>-pbe.txt` and `law-<name><tag>-revpbe.txt`.

use phys::liquid::{fit_site_site, pair_energy, PairEnergy, SiteSite};
use phys::math::Vec3;

const KCAL: f64 = 627.509474;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let temperature: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(500.0);
    let tag = args.get(2).cloned().unwrap_or_default();
    let text = std::fs::read_to_string(format!("pairs-{name}{tag}.txt")).expect("a pairs file");
    // Each line: index R E_pbe E_rev | x y z type ... (first molecule, then second).
    let mut rows: Vec<(f64, f64, Vec<(Vec3, usize)>)> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()) {
        let (head, tail) = line.split_once('|').expect("a | in each line");
        let h: Vec<f64> = head.split_whitespace().map(|x| x.parse().expect("a number")).collect();
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<(Vec3, usize)> = t.chunks(4).map(|c| (Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }, c[3].parse().unwrap())).collect();
        rows.push((h[2], h[3], atoms));
    }
    let n_atoms = rows[0].2.len() / 2;
    let types = rows[0].2.iter().map(|(_, t)| *t).max().unwrap() + 1;
    let mut multiplicity = vec![0usize; types];
    for (_, t) in &rows[0].2[..n_atoms] {
        multiplicity[*t] += 1;
    }
    println!("{name}: {} pairs, {n_atoms} atoms a molecule, {types} site types {multiplicity:?}", rows.len());
    for (label, column) in [("pbe", 0usize), ("revpbe", 1)] {
        let data: Vec<PairEnergy> = rows.iter().map(|(ep, er, atoms)| PairEnergy { a: atoms[..n_atoms].to_vec(), b: atoms[n_atoms..].to_vec(), energy: if column == 0 { *ep } else { *er } }).collect();
        let (train, held): (Vec<(usize, &PairEnergy)>, Vec<(usize, &PairEnergy)>) = data.iter().enumerate().partition(|(i, _)| i % 5 != 4);
        let train: Vec<PairEnergy> = train.into_iter().map(|(_, d)| d.clone()).collect();
        let held: Vec<PairEnergy> = held.into_iter().map(|(_, d)| d.clone()).collect();
        // A neutral, plain start: the fit is what decides.
        let start = SiteSite { charge: vec![0.0; types], pair: vec![[10.0, 2.0, 10.0, 200.0]; types * types] };
        let fit = fit_site_site(&train, &multiplicity, &start, temperature, 1e-3, 5000);
        let rms = |set: &[PairEnergy]| (set.iter().map(|d| (pair_energy(&fit.law, d) - d.energy).powi(2)).sum::<f64>() / set.len() as f64).sqrt();
        let low = |set: &[PairEnergy]| {
            let e_min = set.iter().map(|d| d.energy).fold(f64::INFINITY, f64::min);
            let near: Vec<&PairEnergy> = set.iter().filter(|d| d.energy < e_min + 3.0 / KCAL).collect();
            ((near.iter().map(|d| (pair_energy(&fit.law, d) - d.energy).powi(2)).sum::<f64>() / near.len().max(1) as f64).sqrt(), near.len())
        };
        let (lt, nt) = low(&train);
        let (lh, nh) = low(&held);
        println!("  {label}: {} iterations; RMS {:.3} kcal/mol fitted, {:.3} unseen; within 3 kcal/mol of the lowest: {:.3} over {nt} fitted, {:.3} over {nh} unseen", fit.iterations, rms(&train) * KCAL, rms(&held) * KCAL, lt * KCAL, lh * KCAL);
        println!("    charges {:?}", fit.law.charge);
        std::fs::write(format!("law-{name}{tag}-{label}.txt"), fit.law.to_text()).expect("the law written");
    }
}
