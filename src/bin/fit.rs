//! Route A's law from its pair energies — `PLAY.md` E8.
//!
//! ```sh
//! cargo run --release --bin phys-fit -- water [temperature for the weights, K] [tag] [--bisector] [--c6 0-0=17.48,0-1=6.331,1-1=2.376]
//! ```
//!
//! `--bisector` adds one more charge-only site to every molecule, on the
//! bisector of the first atom's two bonds (the lone pairs' charge, TIP4P-like):
//! its charge follows from neutrality, its pairs carry no repulsion or
//! dispersion, and its distance from the first atom is found by a search over
//! the fits (the weighted residual decides). The law is written with its
//! `bisector <type> <distance>` line, which `phys-bulk` reads. `--c6` holds
//! those dispersion coefficients at the values given (type pair `a-b`, in
//! hartree bohr^6) instead of fitting them: for E6's, the ones the engine
//! derived from the non-local kernel.
//!
//! Reads `pairs-<name><tag>.txt` (written by `phys-pairs`, or by
//! `phys-pairs-gpu` with tag `-gpu`), holds every fifth pair
//! out of the fit, fits the site-site law to the rest for each exchange
//! partner the pairs carry — PBE exchange (nothing fitted in it) and revPBE
//! (published vdW-DF1) — and reports the residual on the pairs it fitted and,
//! the one that is the law's error by the owner's condition, on the pairs it
//! never saw. Writes `law-<name><tag>-pbe.txt` and `law-<name><tag>-revpbe.txt`.

use phys::liquid::{fit_site_site_held, pair_energy, Fitted, Held, PairEnergy, SiteSite};
use phys::math::Vec3;

const KCAL: f64 = 627.509474;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let temperature: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(500.0);
    // One tag, or several joined by commas (`-gpu,-liq-gpu`): every file's
    // pairs are fitted together, held-out ones taken from the whole.
    let tag_arg = args.get(2).filter(|a| !a.starts_with("--")).cloned().unwrap_or_default();
    let tag = tag_arg.replace(',', "+");
    let bisector = args.iter().any(|a| a == "--bisector");
    let c6_held: Vec<((usize, usize), f64)> = args
        .iter()
        .position(|a| a == "--c6")
        .and_then(|i| args.get(i + 1))
        .map(|list| {
            list.split(',')
                .map(|item| {
                    let (k, v) = item.split_once('=').expect("a-b=value");
                    let (a, b) = k.split_once('-').expect("a-b");
                    ((a.parse().expect("a type"), b.parse().expect("a type")), v.parse().expect("a number"))
                })
                .collect()
        })
        .unwrap_or_default();
    let suffix = format!("{}{}", if bisector { "-bis" } else { "" }, if c6_held.is_empty() { "" } else { "-c6" });
    let text: String = tag_arg.split(',').map(|t| std::fs::read_to_string(format!("pairs-{name}{t}.txt")).unwrap_or_else(|_| panic!("no pairs-{name}{t}.txt"))).collect::<Vec<_>>().join("\n");
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
        let mut bisector_distance: Option<f64> = None;
        let (train, held): (Vec<(usize, &PairEnergy)>, Vec<(usize, &PairEnergy)>) = data.iter().enumerate().partition(|(i, _)| i % 5 != 4);
        let train: Vec<PairEnergy> = train.into_iter().map(|(_, d)| d.clone()).collect();
        let held: Vec<PairEnergy> = held.into_iter().map(|(_, d)| d.clone()).collect();
        // A plain start, the fit deciding every number. The first charge is
        // not zero: the electrostatic energy goes as the product of two
        // charges, so with every charge at zero it is stationary and the fit
        // cannot move them — the first run left them at 0.0 and fitted a law
        // with no hydrogen bond in it. (The last charge follows from
        // neutrality.)
        let fit = if bisector {
            let (fit, d) = fit_with_bisector(&train, &multiplicity, types, n_atoms, temperature, &c6_held);
            println!("    bisector site {d:.4} bohr from the first atom (weighted residual {:.3e})", fit.weighted_rms);
            bisector_distance = Some(d);
            fit
        } else {
            let mut charge = vec![0.0; types];
            charge[0] = -0.1;
            let mut start = SiteSite { charge, pair: vec![[10.0, 2.0, 10.0, 200.0]; types * types] };
            for &((a, b), v) in &c6_held {
                start.pair[a * types + b][2] = v;
                start.pair[b * types + a][2] = v;
            }
            let dispersion: Vec<(usize, usize)> = c6_held.iter().map(|&(p, _)| p).collect();
            fit_site_site_held(&train, &multiplicity, &start, temperature, 1e-3, 5000, &Held { pairs: &[], dispersion: &dispersion })
        };
        let with_site = |d: &PairEnergy| -> PairEnergy { match bisector_distance { Some(x) => add_bisector(d, n_atoms, types, x), None => d.clone() } };
        let train: Vec<PairEnergy> = train.iter().map(&with_site).collect();
        let held: Vec<PairEnergy> = held.iter().map(&with_site).collect();
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
        let mut text = fit.law.to_text();
        if let Some(d) = bisector_distance {
            text += &format!("bisector {types} {d:e}\n");
        }
        std::fs::write(format!("law-{name}{tag}-{label}{suffix}.txt"), text).expect("the law written");
    }
}

/// `d` with the extra site added to each molecule: type `types` (one past the
/// atoms' types, which are `0..types`), `distance` bohr from the first atom
/// along the bisector of its bonds to the second and third.
fn add_bisector(d: &PairEnergy, n_atoms: usize, types: usize, distance: f64) -> PairEnergy {
    let extend = |m: &[(Vec3, usize)]| -> Vec<(Vec3, usize)> {
        let o = m[0].0;
        let (u, v) = ((m[1].0 - o).unit(), (m[2].0 - o).unit());
        let mut out = m.to_vec();
        out.push((o + (u + v).unit().scale(distance), types));
        out
    };
    let _ = n_atoms;
    PairEnergy { a: extend(&d.a), b: extend(&d.b), energy: d.energy }
}

/// Fit with one bisector site per molecule, searching its distance. `types`
/// counts the atoms' site types; the new site is type `types`. For each
/// distance the fit is started from several charges (the energy goes as a
/// product of charges and has more than one basin) and the best kept; the
/// distance is the coarse-grid minimum of the weighted residual, refined.
fn fit_with_bisector(train: &[PairEnergy], multiplicity: &[usize], types: usize, n_atoms: usize, temperature: f64, c6_held: &[((usize, usize), f64)]) -> (Fitted, f64) {
    let mut mult = multiplicity.to_vec();
    mult.push(1);
    let all = types + 1;
    let site_pairs: Vec<(usize, usize)> = (0..all).map(|a| (a, types)).collect();
    let dispersion: Vec<(usize, usize)> = c6_held.iter().map(|&(p, _)| p).collect();
    let best_at = |d: f64| -> Fitted {
        let data: Vec<PairEnergy> = train.iter().map(|p| add_bisector(p, n_atoms, types, d)).collect();
        let mut best: Option<Fitted> = None;
        for q_first in [-0.6, -0.2, 0.3, 0.8, 1.3] {
            let mut charge = vec![0.0; all];
            charge[0] = q_first;
            charge[1] = 0.3;
            let mut start = SiteSite { charge, pair: vec![[10.0, 2.0, 10.0, 200.0]; all * all] };
            for &(a, b) in &site_pairs {
                start.pair[a * all + b] = [0.0, 1.0, 0.0, 0.0];
                start.pair[b * all + a] = [0.0, 1.0, 0.0, 0.0];
            }
            for &((a, b), v) in c6_held {
                start.pair[a * all + b][2] = v;
                start.pair[b * all + a][2] = v;
            }
            let fit = fit_site_site_held(&data, &mult, &start, temperature, 1e-3, 5000, &Held { pairs: &site_pairs, dispersion: &dispersion });
            if fit.weighted_rms.is_finite() && best.as_ref().map_or(true, |b| fit.weighted_rms < b.weighted_rms) {
                best = Some(fit);
            }
        }
        best.expect("a fit that is finite")
    };
    let mut scored: Vec<(f64, f64)> = Vec::new();
    for k in 0..9 {
        let d = 0.10 + 0.05 * k as f64;
        let f = best_at(d);
        println!("      d {d:.3}: weighted residual {:.4e}", f.weighted_rms);
        scored.push((f.weighted_rms, d));
    }
    scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let centre = scored[0].1;
    let mut best = (scored[0].0, centre);
    for step in [-0.025, -0.0125, 0.0125, 0.025] {
        let d = centre + step;
        let f = best_at(d);
        println!("      d {d:.4}: weighted residual {:.4e}", f.weighted_rms);
        if f.weighted_rms < best.0 {
            best = (f.weighted_rms, d);
        }
    }
    (best_at(best.1), best.1)
}
