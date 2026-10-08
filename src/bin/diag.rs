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

fn distinct_types(frame: &[ExtraSite]) -> usize {
    let mut t: Vec<usize> = frame.iter().map(|s| s.ty).collect();
    t.sort();
    t.dedup();
    t.len()
}

fn load(path: &str, types_atoms: usize) -> (SiteSite, Vec<f64>, Option<f64>, Vec<ExtraSite>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| panic!("no {path}"));
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    let _ = types_atoms;
    (law, alpha, SiteSite::bisector_from_text(&text).map(|b| b.1), extra_sites_from_text(&text))
}

/// `phys-diag water law.txt --clusters cluster-mp2-water-3-c0.txt,...`: the
/// three-body energy of each trimer (its cluster interaction less the sum of its
/// pairs') as the reference gives it and as the law does, which has only its
/// induction to make one.
fn clusters(law_path: &str, files: &str) {
    use phys::induction::{solve, Charge, Cluster, Link, PolSite};
    let (law, alpha, bis, frame) = load(law_path, 3);
    let types = law.charge.len();
    let atom_types = types - bis.is_some() as usize - distinct_types(&frame);
    let mut total_ref = 0.0;
    let mut total_law = 0.0;
    let mut n = 0.0;
    println!("trimers against {law_path}: three-body energy (kcal/mol), reference then the law's, and each cluster's interaction");
    for f in files.split(',') {
        let text = std::fs::read_to_string(f).unwrap_or_else(|_| panic!("no {f}"));
        let mut mols: Vec<Vec<(Vec3, usize)>> = Vec::new();
        let mut pairs_sum = 0.0;
        let mut cluster = 0.0;
        for l in text.lines() {
            let w: Vec<&str> = l.split_whitespace().collect();
            if w.get(1) == Some(&"mol") {
                let v: Vec<f64> = w[3..].iter().map(|x| x.parse().expect("number")).collect();
                let atoms: Vec<(Vec3, usize)> = v.chunks(3).zip([0usize, 1, 1]).map(|(c, t)| (Vec3 { x: c[0], y: c[1], z: c[2] }, t)).collect();
                let atoms = match bis {
                    Some(d) => with_bisector_site(&atoms, atom_types, d),
                    None => atoms,
                };
                mols.push(with_frame_sites(&atoms, &frame));
            } else if w.first() == Some(&"pair") {
                pairs_sum += w[5].parse::<f64>().expect("number") + w[6].parse::<f64>().expect("number");
            } else if w.first() == Some(&"cluster") {
                cluster = w[2].parse::<f64>().expect("number") + w[3].parse::<f64>().expect("number");
            }
        }
        if mols.len() != 3 {
            continue;
        }
        let mut cl = Cluster::default();
        for m in &mols {
            cl.pol.push(m.iter().filter(|(_, t)| alpha[*t] > 0.0).map(|(p, t)| PolSite { pos: *p, alpha: alpha[*t] }).collect());
            cl.charges.push(m.iter().filter(|(_, t)| law.charge[*t] != 0.0).map(|(p, t)| Charge { pos: *p, q: law.charge[*t], sigma: law.sigma.get(*t).copied().unwrap_or(0.0) }).collect());
        }
        let mut links = Vec::new();
        for i in 0..3 {
            for j in i + 1..3 {
                links.push(Link { i, j, shift: Vec3::ZERO, weight: 1.0, dweight: 0.0, d: Vec3::ZERO });
            }
        }
        let whole = solve(&cl, &links, None, 1e-12).energy;
        let mut pairwise = 0.0;
        for i in 0..3 {
            for j in i + 1..3 {
                let pe = PairEnergy { a: mols[i].clone(), b: mols[j].clone(), energy: 0.0 };
                pairwise += phys::liquid::induction_pair_energy(&alpha, &law.charge, &law.sigma, &pe);
            }
        }
        let law3 = (whole - pairwise) * KCAL;
        let ref3 = (cluster - pairs_sum) * KCAL;
        println!("  {f}: reference {ref3:+.3}, law {law3:+.3}; cluster {:.3} kcal/mol", cluster * KCAL);
        total_ref += ref3;
        total_law += law3;
        n += 1.0;
    }
    println!("  mean over {n}: reference {:+.3}, law {:+.3}", total_ref / n, total_law / n);
}

/// `phys-diag water law.txt --es es-water.txt tags`: the law's electrostatics
/// alone (its charges and widths, no repulsion, dispersion or induction)
/// against the exact electrostatic energy of the two molecules' own densities
/// (`phys-es-gpu`), pair by pair, by distance and by how nearly the contact is
/// a hydrogen bond.
fn electrostatics(law_path: &str, es_file: &str, tags: &str) {
    let (law, _alpha, bis, frame) = load(law_path, 3);
    let types = law.charge.len();
    let atom_types = types - bis.is_some() as usize - distinct_types(&frame);
    let mut only_es = law.clone();
    for p in only_es.pair.iter_mut() {
        *p = [0.0, 1.0, 0.0, 0.0];
    }
    let es: std::collections::HashMap<usize, f64> = std::fs::read_to_string(es_file)
        .unwrap_or_else(|_| panic!("no {es_file}"))
        .lines()
        .filter_map(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            Some((w.first()?.parse().ok()?, w.get(2)?.parse().ok()?))
        })
        .collect();
    let text: String = tags.split(',').map(|t| std::fs::read_to_string(format!("pairs-water{t}.txt")).unwrap_or_default()).collect::<Vec<_>>().join("\n");
    // (O-O in A, smallest H-O...O angle, reference, law)
    let mut rows: Vec<(f64, f64, f64, f64)> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#') && l.contains('|')) {
        let (head, tail) = line.split_once('|').expect("a |");
        let k: usize = head.split_whitespace().next().and_then(|x| x.parse().ok()).expect("an index");
        let Some(&reference) = es.get(&k) else { continue };
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<(Vec3, usize)> = t.chunks(4).map(|c| (Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }, c[3].parse().unwrap())).collect();
        let n = atoms.len() / 2;
        let (a, b) = (&atoms[..n], &atoms[n..]);
        let angle = |o: Vec3, h: Vec3, other: Vec3| -> f64 { (h - o).unit().dot((other - o).unit()).clamp(-1.0, 1.0).acos().to_degrees() };
        let theta = [angle(a[0].0, a[1].0, b[0].0), angle(a[0].0, a[2].0, b[0].0), angle(b[0].0, b[1].0, a[0].0), angle(b[0].0, b[2].0, a[0].0)].iter().cloned().fold(f64::INFINITY, f64::min);
        let place = |m: &[(Vec3, usize)]| -> Vec<(Vec3, usize)> {
            let m = match bis {
                Some(d) => with_bisector_site(m, atom_types, d),
                None => m.to_vec(),
            };
            with_frame_sites(&m, &frame)
        };
        let pe = PairEnergy { a: place(a), b: place(b), energy: 0.0 };
        rows.push(((a[0].0 - b[0].0).norm() * BOHR, theta, reference * KCAL, pair_energy(&only_es, &pe) * KCAL));
    }
    println!("electrostatics of {law_path} against the densities' own, {} pairs (kcal/mol)", rows.len());
    println!("  by O-O (A):      n   exact mean   mean err   rms err");
    for w in [0.0, 2.5, 2.8, 3.1, 3.4, 3.8, 4.7].windows(2) {
        let sel: Vec<&(f64, f64, f64, f64)> = rows.iter().filter(|r| r.0 >= w[0] && r.0 < w[1]).collect();
        if sel.is_empty() {
            continue;
        }
        let n = sel.len() as f64;
        println!("  {:>4.1} - {:>4.1}   {:>4}  {:>9.3}  {:>9.3}  {:>9.3}", w[0], w[1], sel.len(), sel.iter().map(|r| r.2).sum::<f64>() / n, sel.iter().map(|r| r.3 - r.2).sum::<f64>() / n, (sel.iter().map(|r| (r.3 - r.2).powi(2)).sum::<f64>() / n).sqrt());
    }
    println!("  inside 3.3 A by the smallest H-O...O angle:");
    for (lo, hi) in [(0.0, 15.0), (15.0, 30.0), (30.0, 50.0), (50.0, 80.0), (80.0, 181.0)] {
        let sel: Vec<&(f64, f64, f64, f64)> = rows.iter().filter(|r| r.0 < 3.3 && r.1 >= lo && r.1 < hi).collect();
        if sel.is_empty() {
            continue;
        }
        let n = sel.len() as f64;
        println!("    {lo:>3.0}-{hi:>3.0} deg {:>4}  {:>9.3}  {:>9.3}  {:>9.3}", sel.len(), sel.iter().map(|r| r.2).sum::<f64>() / n, sel.iter().map(|r| r.3 - r.2).sum::<f64>() / n, (sel.iter().map(|r| (r.3 - r.2).powi(2)).sum::<f64>() / n).sqrt());
    }
}

/// `phys-diag water cloud-law.txt --fit-es es-water.txt tags`: the core-and-cloud
/// model's widths, populations and bisector distance fitted to the exact
/// electrostatic energies themselves (the densities' first-order interaction,
/// `phys-es-gpu`), every other pair held out. Nothing here is a supermolecular
/// energy: the reference is the monomers' own densities, and the cores stay the
/// elements'. Writes `<law>-fit.txt`.
fn fit_electrostatics(law_path: &str, es_file: &str, tags: &str) {
    let text0 = std::fs::read_to_string(law_path).unwrap_or_else(|_| panic!("no {law_path}"));
    let law0 = SiteSite::from_text(&text0).expect("a readable law");
    let frame0 = extra_sites_from_text(&text0);
    let bis0 = SiteSite::bisector_from_text(&text0).expect("a bisector line");
    // Types: 0, 1 cores; 2 the bisector cloud; 3 the oxygen cloud; 4 the hydrogen clouds.
    assert_eq!(law0.charge.len(), 5, "a core-and-cloud law (phys-esp --cloud)");
    let es: std::collections::HashMap<usize, f64> = std::fs::read_to_string(es_file)
        .unwrap_or_else(|_| panic!("no {es_file}"))
        .lines()
        .filter_map(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            Some((w.first()?.parse().ok()?, w.get(2)?.parse().ok()?))
        })
        .collect();
    let text: String = tags.split(',').map(|t| std::fs::read_to_string(format!("pairs-water{t}.txt")).unwrap_or_default()).collect::<Vec<_>>().join("\n");
    let mut rows: Vec<(usize, Vec<(Vec3, usize)>, Vec<(Vec3, usize)>, f64)> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#') && l.contains('|')) {
        let (head, tail) = line.split_once('|').expect("a |");
        let k: usize = head.split_whitespace().next().and_then(|x| x.parse().ok()).expect("an index");
        let Some(&reference) = es.get(&k) else { continue };
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<(Vec3, usize)> = t.chunks(4).map(|c| (Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }, c[3].parse().unwrap())).collect();
        let n = atoms.len() / 2;
        rows.push((k, atoms[..n].to_vec(), atoms[n..].to_vec(), reference * KCAL));
    }
    // A model from parameters: ln widths of the three clouds, the bisector
    // distance, the bisector and oxygen cloud populations (the hydrogens' follow
    // from the eight valence electrons).
    let r_oh = frame0.iter().find(|s| s.ty == 4).map(|s| s.frame[0]).expect("hydrogen clouds");
    let build = |x: &[f64]| -> (SiteSite, f64) {
        let (n_m, n_o) = (x[4], x[5]);
        let n_h = (8.0 - n_m - n_o) / 2.0;
        let mut l = law0.clone();
        l.charge = vec![law0.charge[0], law0.charge[1], -n_m, -n_o, -n_h];
        l.sigma = vec![0.0, 0.0, x[0].exp(), x[1].exp(), x[2].exp()];
        for p in l.pair.iter_mut() {
            *p = [0.0, 1.0, 0.0, 0.0];
        }
        (l, x[3].abs().max(1e-3))
    };
    let energy = |x: &[f64], r: &(usize, Vec<(Vec3, usize)>, Vec<(Vec3, usize)>, f64)| -> f64 {
        let (l, d) = build(x);
        let place = |m: &[(Vec3, usize)]| with_frame_sites(&with_bisector_site(m, 2, d), &frame0);
        let pe = PairEnergy { a: place(&r.1), b: place(&r.2), energy: 0.0 };
        pair_energy(&l, &pe) * KCAL
    };
    let (train, test): (Vec<_>, Vec<_>) = rows.iter().enumerate().partition(|(i, _)| i % 2 == 0);
    let train: Vec<_> = train.into_iter().map(|(_, r)| r).collect();
    let test: Vec<_> = test.into_iter().map(|(_, r)| r).collect();
    // Weights: the contact region counts most, as it does in a liquid.
    let weight = |r: &(usize, Vec<(Vec3, usize)>, Vec<(Vec3, usize)>, f64)| -> f64 { 1.0 / (1.0 + ((r.1[0].0 - r.2[0].0).norm() * BOHR - 2.7).max(0.0).powi(2)) };
    // `--dipole <au>`: the molecule's dipole held to this value (the MP2 one from
    // the field derivative of its energy), because the densities the reference
    // energies come from are Hartree-Fock's and a few per cent too polar.
    let dipole_target: Option<f64> = std::env::args().position(|a| a == "--dipole").and_then(|i| std::env::args().nth(i + 1)).and_then(|v| v.parse().ok());
    let dipole_of = |x: &[f64]| -> f64 {
        let (l, d) = build(x);
        let m = with_frame_sites(&with_bisector_site(&rows[0].1, 2, d), &frame0);
        let mut mu = Vec3::ZERO;
        for (p, t) in &m {
            mu = mu + p.scale(l.charge[*t]);
        }
        mu.norm()
    };
    let cost = |x: &[f64]| -> f64 {
        let base = train.iter().map(|r| weight(r) * (energy(x, r) - r.3).powi(2)).sum::<f64>() / train.len() as f64;
        match dipole_target {
            Some(t) => base + 3000.0 * ((dipole_of(x) - t) / t).powi(2),
            None => base,
        }
    };
    let x0: Vec<f64> = vec![law0.sigma[2].ln(), law0.sigma[3].ln(), law0.sigma[4].ln(), bis0.1, -law0.charge[2], -law0.charge[3]];
    let report = |label: &str, x: &[f64]| {
        let stat = |set: &[&(usize, Vec<(Vec3, usize)>, Vec<(Vec3, usize)>, f64)]| -> (f64, f64) {
            let n = set.len() as f64;
            ((set.iter().map(|r| energy(x, r) - r.3).sum::<f64>()) / n, (set.iter().map(|r| (energy(x, r) - r.3).powi(2)).sum::<f64>() / n).sqrt())
        };
        println!("  {label}: train mean/rms {:+.3}/{:.3}, held-out {:+.3}/{:.3} kcal/mol", stat(&train).0, stat(&train).1, stat(&test).0, stat(&test).1);
        for w in [0.0, 2.8, 3.1, 3.4, 4.7].windows(2) {
            let sel: Vec<&(usize, Vec<(Vec3, usize)>, Vec<(Vec3, usize)>, f64)> = rows.iter().filter(|r| { let d = (r.1[0].0 - r.2[0].0).norm() * BOHR; d >= w[0] && d < w[1] }).collect();
            if !sel.is_empty() {
                let (m, s) = stat(&sel);
                println!("      O-O {:.1}-{:.1} A: {} pairs, mean {m:+.3}, rms {s:.3}", w[0], w[1], sel.len());
            }
        }
    };
    report("before", &x0);
    let nn = x0.len();
    let mut simplex: Vec<Vec<f64>> = vec![x0.clone()];
    for i in 0..nn {
        let mut v = x0.clone();
        v[i] += if i < 3 { 0.2 } else if i == 3 { 0.1 } else { 0.3 };
        simplex.push(v);
    }
    let mut vals: Vec<f64> = simplex.iter().map(|v| cost(v)).collect();
    for _ in 0..3000 {
        let mut idx: Vec<usize> = (0..=nn).collect();
        idx.sort_by(|&a, &b| vals[a].partial_cmp(&vals[b]).unwrap());
        simplex = idx.iter().map(|&i| simplex[i].clone()).collect();
        vals = idx.iter().map(|&i| vals[i]).collect();
        let centroid: Vec<f64> = (0..nn).map(|kk| simplex[..nn].iter().map(|v| v[kk]).sum::<f64>() / nn as f64).collect();
        let worst = simplex[nn].clone();
        let along = |t: f64| -> Vec<f64> { (0..nn).map(|kk| centroid[kk] + t * (worst[kk] - centroid[kk])).collect() };
        let refl = along(-1.0);
        let fr = cost(&refl);
        if fr < vals[0] {
            let ex = along(-2.0);
            let fe = cost(&ex);
            if fe < fr {
                simplex[nn] = ex;
                vals[nn] = fe;
            } else {
                simplex[nn] = refl;
                vals[nn] = fr;
            }
        } else if fr < vals[nn - 1] {
            simplex[nn] = refl;
            vals[nn] = fr;
        } else {
            let con = along(0.5);
            let fc = cost(&con);
            if fc < vals[nn] {
                simplex[nn] = con;
                vals[nn] = fc;
            } else {
                for i in 1..=nn {
                    for kk in 0..nn {
                        simplex[i][kk] = simplex[0][kk] + 0.5 * (simplex[i][kk] - simplex[0][kk]);
                    }
                    vals[i] = cost(&simplex[i]);
                }
            }
        }
    }
    let best = simplex[0].clone();
    report("after ", &best);
    println!("  dipole {:.4} au ({:.3} D)", dipole_of(&best), dipole_of(&best) / 0.393430307);
    let (l, d) = build(&best);
    let out = format!(
        "charges {}\nsigma {}\nbisector 2 {:e}\nsite 3 0e0 0e0 0e0\nsite 4 {r_oh:e} 0e0 0e0\nsite 4 0e0 {r_oh:e} 0e0\n",
        l.charge.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "),
        l.sigma.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "),
        d
    );
    println!("  charges {:?}, widths {:?}, bisector {d:.4} bohr", l.charge, l.sigma);
    std::fs::write(law_path.replace(".txt", "-fit.txt"), out).expect("written");
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if let Some(i) = args.iter().position(|a| a == "--fit-es") {
        fit_electrostatics(&args[1], &args[i + 1], args.last().expect("tags"));
        return;
    }
    if let Some(i) = args.iter().position(|a| a == "--es") {
        electrostatics(&args[1], &args[i + 1], args.last().expect("tags"));
        return;
    }
    if let Some(i) = args.iter().position(|a| a == "--clusters") {
        clusters(&args[1], &args[i + 1]);
        return;
    }
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
        let atom_types = types - bis.is_some() as usize - distinct_types(&frame);
        println!("law {lp} against the {} column", if li == 0 { "total (MP2)" } else { "Hartree-Fock" });
        println!("  O-O (A)        n   ref mean   mean err   rms err   (kcal/mol)");
        let mut by_angle: Vec<(f64, f64, f64, f64, usize)> = Vec::new();
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
                if row.0 * BOHR < 3.3 {
                    // The smallest angle H-O...O over both molecules' hydrogens: how
                    // nearly the contact is a hydrogen bond.
                    let (o1, o2) = (row.3[0].0, row.3[n_atoms].0);
                    let angle = |o: Vec3, h: Vec3, other: Vec3| -> f64 { (h - o).unit().dot((other - o).unit()).clamp(-1.0, 1.0).acos().to_degrees() };
                    let theta = [angle(o1, row.3[1].0, o2), angle(o1, row.3[2].0, o2), angle(o2, row.3[n_atoms + 1].0, o1), angle(o2, row.3[n_atoms + 2].0, o1)].iter().cloned().fold(f64::INFINITY, f64::min);
                    by_angle.push((row.0 * BOHR, theta, reference * KCAL, d, 0));
                }
                s += d;
                s2 += d * d;
                r += reference * KCAL;
            }
            let n = sel.len() as f64;
            println!("  {:>4.1} - {:>4.1}  {:>4}  {:>9.3}  {:>9.3}  {:>9.3}", w[0], w[1].min(9.9), sel.len(), r / n, s / n, (s2 / n).sqrt());
        }
        println!("  inside 3.3 A by the smallest H-O...O angle (0 is a straight hydrogen bond):");
        for (lo, hi) in [(0.0, 15.0), (15.0, 30.0), (30.0, 50.0), (50.0, 80.0), (80.0, 181.0)] {
            let sel: Vec<&(f64, f64, f64, f64, usize)> = by_angle.iter().filter(|e| e.1 >= lo && e.1 < hi).collect();
            if sel.is_empty() {
                continue;
            }
            let n = sel.len() as f64;
            println!("    {lo:>3.0}-{hi:>3.0} deg {:>4}  {:>9.3}  {:>9.3}  {:>9.3}", sel.len(), sel.iter().map(|e| e.2).sum::<f64>() / n, sel.iter().map(|e| e.3).sum::<f64>() / n, (sel.iter().map(|e| e.3 * e.3).sum::<f64>() / n).sqrt());
        }
    }
}
