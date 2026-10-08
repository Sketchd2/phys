//! Route A's law from its pair energies — `PLAY.md` E8.
//!
//! ```sh
//! cargo run --release --bin phys-fit -- water [temperature for the weights, K] [tag] [--bisector] [--c6 0-0=17.48,0-1=6.331,1-1=2.376]
//! ```
//!
//! `--alpha polar-water.txt` gives each atom class's polarisability (as
//! `phys-polar` derives it from the molecule's own electrons): every pair's
//! energy in the fit then includes the induction of the two molecules, in the
//! charges being fitted, and the law is written with its `alpha` lines.
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

use phys::liquid::{extra_sites_from_text, fit_site_site_polarised, pair_energy, with_bisector_site, with_frame_sites, ExtraSite, Fitted, Held, PairEnergy, Polarisable, SiteSite};
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
    let use_sigma = args.iter().any(|a| a == "--sigma");
    // `--start law.txt`: begin the atom pairs' numbers from this law's
    // rather than from a plain guess, so that terms added to it (`--lp`) are
    // perturbations of a fit that already works.
    let start_law: Option<SiteSite> = args.iter().position(|a| a == "--start").and_then(|i| args.get(i + 1)).map(|f| SiteSite::from_text(&std::fs::read_to_string(f).unwrap_or_else(|_| panic!("no {f}"))).expect("a readable law"));
    // `--overlap`: the repulsion is that of the charge clouds' overlap (the
    // density's own sizes, one strength per pair) and not a Born-Mayer per atom pair.
    let use_overlap = args.iter().any(|a| a == "--overlap");
    // `--no-c8`: no C8 at all.
    let no_c8 = args.iter().any(|a| a == "--no-c8");
    // `--esp esp-water.txt`: the charges and the bisector site from the
    // molecule's own density (`phys-esp`), held while the rest is fitted.
    let esp: Option<(Vec<f64>, f64, Vec<f64>)> = args.iter().position(|a| a == "--esp").and_then(|i| args.get(i + 1)).map(|f| {
        let t = std::fs::read_to_string(f).unwrap_or_else(|_| panic!("no {f}"));
        let law = SiteSite::from_text(&t).expect("charges in the esp file");
        (law.charge, SiteSite::bisector_from_text(&t).expect("a bisector line").1, law.sigma)
    });
    let alpha_file = args.iter().position(|a| a == "--alpha").and_then(|i| args.get(i + 1)).cloned();
    let suffix = format!("{}{}{}{}{}", if bisector { "-bis" } else { "" }, if c6_held.is_empty() { "" } else { "-c6" }, if esp.is_some() { "-esp" } else { "" }, if alpha_file.is_some() { "-ind" } else { "" }, if use_sigma { "-sig" } else { "" });
    // Each atom type's polarisability, from `phys-polar`'s file, if asked.
    let alpha_atoms: Option<Vec<f64>> = alpha_file.as_ref().map(|f| {
        let t = std::fs::read_to_string(f).unwrap_or_else(|_| panic!("no {f}"));
        SiteSite::alpha_from_text(&t, 16)
    });
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
    // `--lp esp-sites.txt`: sites (from the localised orbitals, `phys-esp
    // --sites`) that carry a repulsion of their own and no charge: the
    // anisotropy an atom-atom Born-Mayer lacks. One more site type, after the
    // bisector site's.
    // `--charge-sites file`: sites that carry charge only, of the types the
    // file gives (the density's own core-and-cloud model, `phys-esp --cloud`).
    let charge_mode = args.iter().any(|a| a == "--charge-sites");
    let charge_sites: Vec<ExtraSite> = args.iter().position(|a| a == "--charge-sites").and_then(|i| args.get(i + 1)).map(|f| extra_sites_from_text(&std::fs::read_to_string(f).unwrap_or_else(|_| panic!("no {f}")))).unwrap_or_default();
    let lp_sites: Vec<ExtraSite> = if charge_mode { charge_sites } else { args
        .iter()
        .position(|a| a == "--lp")
        .and_then(|i| args.get(i + 1))
        .map(|f| extra_sites_from_text(&std::fs::read_to_string(f).unwrap_or_else(|_| panic!("no {f}"))).into_iter().map(|s| ExtraSite { ty: types + 1, frame: s.frame }).collect())
        .unwrap_or_default() };
    let extra_classes: usize = { let mut t: Vec<usize> = lp_sites.iter().map(|s| s.ty).collect(); t.sort(); t.dedup(); t.len() };
    // The polarisability of every site type the fit sees: the atoms', and zero
    // for the bisector site.
    let alpha_full: Option<Vec<f64>> = alpha_atoms.as_ref().map(|a| {
        let mut v: Vec<f64> = a[..types].to_vec();
        if bisector {
            v.push(0.0);
        }
        for _ in 0..extra_classes {
            v.push(0.0);
        }
        println!("  induced dipoles on, polarisabilities by type {v:?} bohr^3");
        v
    });
    // `--two-stage`: the Hartree-Fock column first, with no dispersion, for the
    // repulsion; then the total energy with that repulsion held, for the
    // dispersion and its own damping (needs `--bisector --esp`).
    let two_stage = args.iter().any(|a| a == "--two-stage");
    assert!(!two_stage || (bisector && esp.is_some()), "--two-stage needs --bisector and --esp");
    let mut hf_law: Option<SiteSite> = None;
    let order: [(&str, usize); 2] = if two_stage { [("revpbe", 1), ("pbe", 0)] } else { [("pbe", 0), ("revpbe", 1)] };
    for (label, column) in order {
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
            let (fit, d) = fit_with_bisector(&train, &multiplicity, types, temperature, &c6_held, alpha_full.as_deref(), use_sigma, esp.as_ref(), &lp_sites, charge_mode, no_c8, use_overlap, start_law.as_ref(), match (two_stage, &hf_law) { (false, _) => Stage::Single, (true, None) => Stage::Hf, (true, Some(l)) => Stage::Correlation(l) });
            println!("    bisector site {d:.4} bohr from the first atom (weighted residual {:.3e})", fit.weighted_rms);
            bisector_distance = Some(d);
            fit
        } else {
            let mut charge = vec![0.0; types];
            charge[0] = -0.1;
            let mut start = SiteSite { charge, pair: vec![[10.0, 2.0, 10.0, 200.0]; types * types], sigma: Vec::new(), damp: Vec::new(), overlap: Vec::new() };
            for &((a, b), v) in &c6_held {
                start.pair[a * types + b][2] = v;
                start.pair[b * types + a][2] = v;
            }
            let dispersion: Vec<(usize, usize)> = c6_held.iter().map(|&(p, _)| p).collect();
            fit_site_site_polarised(&train, &multiplicity, &start, temperature, 1e-3, 5000, &Held { pairs: &[], dispersion: &dispersion, sigma: use_sigma, charges: esp.is_some(), damp: false, repulsion: false, no_dispersion: false, no_dispersion_pairs: &[], no_c8, overlap: use_overlap, no_born_mayer: use_overlap }, alpha_full.as_deref())
        };
        if two_stage && column == 1 {
            hf_law = Some(fit.law.clone());
        }
        let with_site = |d: &PairEnergy| -> PairEnergy { match bisector_distance { Some(x) => add_bisector(d, types, x, &lp_sites), None => d.clone() } };
        let train: Vec<PairEnergy> = train.iter().map(&with_site).collect();
        let held: Vec<PairEnergy> = held.iter().map(&with_site).collect();
        let eval = |d: &PairEnergy| -> f64 {
            match &alpha_full {
                Some(a) => pair_energy(&Polarisable { law: &fit.law, alpha: a.clone() }, d),
                None => pair_energy(&fit.law, d),
            }
        };
        let rms = |set: &[PairEnergy]| (set.iter().map(|d| (eval(d) - d.energy).powi(2)).sum::<f64>() / set.len() as f64).sqrt();
        let low = |set: &[PairEnergy]| {
            let e_min = set.iter().map(|d| d.energy).fold(f64::INFINITY, f64::min);
            let near: Vec<&PairEnergy> = set.iter().filter(|d| d.energy < e_min + 3.0 / KCAL).collect();
            ((near.iter().map(|d| (eval(d) - d.energy).powi(2)).sum::<f64>() / near.len().max(1) as f64).sqrt(), near.len())
        };
        let (lt, nt) = low(&train);
        let (lh, nh) = low(&held);
        println!("  {label}: {} iterations; RMS {:.3} kcal/mol fitted, {:.3} unseen; within 3 kcal/mol of the lowest: {:.3} over {nt} fitted, {:.3} over {nh} unseen", fit.iterations, rms(&train) * KCAL, rms(&held) * KCAL, lt * KCAL, lh * KCAL);
        println!("    charges {:?}", fit.law.charge);
        let mut text = fit.law.to_text();
        for s in &lp_sites {
            text += &format!("site {} {:e} {:e} {:e}\n", s.ty, s.frame[0], s.frame[1], s.frame[2]);
        }
        if let Some(d) = bisector_distance {
            text += &format!("bisector {types} {d:e}\n");
        }
        if let Some(a) = &alpha_atoms {
            for (t, v) in a[..types].iter().enumerate() {
                text += &format!("alpha {t} {v:e}\n");
            }
        }
        std::fs::write(format!("law-{name}{tag}-{label}{suffix}.txt"), text).expect("the law written");
    }
}

/// Which energies a fit is of: all of them at once, or the two stages of
/// `--two-stage`.
#[derive(Clone, Copy)]
enum Stage<'a> {
    Single,
    /// Hartree-Fock's energies, dispersion held at zero.
    Hf,
    /// The total energy around this Hartree-Fock law's repulsion.
    Correlation(&'a SiteSite),
}

/// `d` with the extra site added to each molecule: type `types` (one past the
/// atoms' types, which are `0..types`), `distance` bohr from the first atom
/// along the bisector of its bonds to the second and third.
fn add_bisector(d: &PairEnergy, types: usize, distance: f64, extra: &[ExtraSite]) -> PairEnergy {
    PairEnergy { a: with_frame_sites(&with_bisector_site(&d.a, types, distance), extra), b: with_frame_sites(&with_bisector_site(&d.b, types, distance), extra), energy: d.energy }
}

/// Fit with one bisector site per molecule, searching its distance. `types`
/// counts the atoms' site types; the new site is type `types`. For each
/// distance the fit is started from several charges (the energy goes as a
/// product of charges and has more than one basin) and the best kept; the
/// distance is the coarse-grid minimum of the weighted residual, refined.
fn fit_with_bisector(train: &[PairEnergy], multiplicity: &[usize], types: usize, temperature: f64, c6_held: &[((usize, usize), f64)], alpha: Option<&[f64]>, use_sigma: bool, esp: Option<&(Vec<f64>, f64, Vec<f64>)>, lp: &[ExtraSite], charge_only: bool, no_c8: bool, use_overlap: bool, start_law: Option<&SiteSite>, stage: Stage) -> (Fitted, f64) {
    let mut mult = multiplicity.to_vec();
    mult.push(1);
    let mut classes: Vec<usize> = lp.iter().map(|s| s.ty).collect();
    classes.sort();
    classes.dedup();
    for t in &classes {
        mult.push(lp.iter().filter(|s| s.ty == *t).count());
    }
    let all = types + 1 + classes.len();
    let site_pairs: Vec<(usize, usize)> = if charge_only { (types..all).flat_map(|t| (0..=t).map(move |a| (a, t))).collect() } else { (0..all).map(|a| (a, types)).collect() };
    // The repelling sites' pairs with the atoms and with each other are fitted
    // for repulsion and held at zero for dispersion.
    let lp_pairs: Vec<(usize, usize)> = if lp.is_empty() || charge_only { Vec::new() } else { (0..types).chain(std::iter::once(types + 1)).map(|a| (a, types + 1)).collect() };
    let dispersion: Vec<(usize, usize)> = c6_held.iter().map(|&(p, _)| p).collect();
    let best_at = |d: f64| -> Fitted {
        let data: Vec<PairEnergy> = train.iter().map(|p| add_bisector(p, types, d, lp)).collect();
        let mut best: Option<Fitted> = None;
        let starts: Vec<f64> = if esp.is_some() { vec![0.0] } else { vec![-0.6, -0.2, 0.3, 0.8, 1.3] };
        for q_first in starts {
            let mut charge = vec![0.0; all];
            charge[0] = q_first;
            charge[1] = 0.3;
            let mut sigma_start = Vec::new();
            if let Some((q, _, sg)) = esp {
                charge = q.clone();
                sigma_start = sg.clone();
                while charge.len() < all {
                    charge.push(0.0);
                    sigma_start.push(0.5);
                }
            }
            let mut start = SiteSite { charge, pair: vec![[10.0, 2.0, 10.0, 200.0]; all * all], sigma: sigma_start, damp: Vec::new(), overlap: Vec::new() };
            if let Stage::Correlation(base) = stage {
                start = base.clone();
                // The Hartree-Fock stage left no dispersion, and a number at
                // zero (its logarithm far below any step) cannot be moved:
                // start from the engine's C6 where it has one, and a plain C8.
                for &(a, b) in &site_pairs {
                    let _ = (a, b);
                }
                for ta in 0..types {
                    for tb in 0..types {
                        let v = c6_held.iter().find(|((x, y), _)| (*x, *y) == (ta, tb) || (*y, *x) == (ta, tb)).map(|(_, v)| *v).unwrap_or(10.0);
                        start.pair[ta * all + tb][2] = v;
                        start.pair[ta * all + tb][3] = 200.0;
                    }
                }
            }
            for &(a, b) in &site_pairs {
                start.pair[a * all + b] = [0.0, 1.0, 0.0, 0.0];
                start.pair[b * all + a] = [0.0, 1.0, 0.0, 0.0];
            }
            if no_c8 {
                for p in start.pair.iter_mut() {
                    p[3] = 0.0;
                }
            }
            if use_overlap {
                // No Born-Mayer of its own; B stays as the dispersion's damping.
                for p in start.pair.iter_mut() {
                    p[0] = 1e-300;
                    p[1] = 2.0;
                }
            }
            if let Some(law) = start_law {
                let lt = law.charge.len();
                for a in 0..types {
                    for b in 0..types {
                        start.pair[a * all + b] = law.pair[a * lt + b];
                    }
                }
            }
            // A repelling site begins nearly absent when there is a law to
            // perturb, and at a plain guess when there is not.
            let lp_start = if start_law.is_some() { [1e-4, 2.0, 0.0, 0.0] } else { [5.0, 1.5, 0.0, 0.0] };
            for &(a, b) in &lp_pairs {
                start.pair[a * all + b] = lp_start;
                start.pair[b * all + a] = lp_start;
            }
            if matches!(stage, Stage::Hf) {
                for p in start.pair.iter_mut() {
                    p[2] = 0.0;
                    p[3] = 0.0;
                }
            }
            if matches!(stage, Stage::Single) {
                for &((a, b), v) in c6_held {
                    start.pair[a * all + b][2] = v;
                    start.pair[b * all + a][2] = v;
                }
            }
            let held = Held { pairs: &site_pairs, dispersion: if matches!(stage, Stage::Correlation(_)) { &[] } else { &dispersion }, sigma: use_sigma, charges: esp.is_some(), damp: matches!(stage, Stage::Correlation(_)), repulsion: matches!(stage, Stage::Correlation(_)), no_dispersion: matches!(stage, Stage::Hf), no_dispersion_pairs: &lp_pairs, no_c8, overlap: use_overlap, no_born_mayer: use_overlap };
            let fit = fit_site_site_polarised(&data, &mult, &start, temperature, 1e-3, 5000, &held, alpha);
            if fit.weighted_rms.is_finite() && best.as_ref().map_or(true, |b| fit.weighted_rms < b.weighted_rms) {
                best = Some(fit);
            }
        }
        best.expect("a fit that is finite")
    };
    if let Some((_, d, _)) = esp {
        return (best_at(*d), *d);
    }
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
