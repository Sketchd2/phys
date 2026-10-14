//! The many-body energy of clusters cut from a liquid snapshot, as a pair law
//! with induced dipoles gives it (PLAY.md E8d).
//!
//! ```sh
//! phys-manybody water bulk-water-298.snap law.txt [--max 10] [--centres 0,20,45 | --every 4]
//! ```
//!
//! For a cluster of `n` molecules (the centre and its `n - 1` nearest), the law's
//! energy with every induced dipole solved together, less the sum over the pairs
//! of the law's energy with each pair solved alone, is what the law has that is
//! not additive: only its induction makes any. It is split into the sum over
//! triples of each triple's non-additive part (three-body) and what is left
//! (four-body and up), and set beside the reference's by `phys-cluster-gpu`.

use phys::induction::{solve, Charge, Cluster, Link, PolSite};
use phys::liquid::{extra_sites_from_text, pair_energy, with_bisector_site, with_frame_sites, PairEnergy, SiteSite};
use phys::pairs::Snapshot;
use phys::math::Vec3;

const KCAL: f64 = 627.509474;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let snapshot = args.get(1).expect("a snapshot");
    let law_path = args.get(2).expect("a law file");
    let flag = |f: &str| args.iter().position(|a| a == f).and_then(|i| args.get(i + 1)).cloned();
    let max: usize = flag("--max").and_then(|s| s.parse().ok()).unwrap_or(10);
    let text = std::fs::read_to_string(law_path).unwrap_or_else(|_| panic!("no {law_path}"));
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    let bis = SiteSite::bisector_from_text(&text).map(|b| b.1);
    let frame = extra_sites_from_text(&text);
    let mut distinct: Vec<usize> = frame.iter().map(|s| s.ty).collect();
    distinct.sort();
    distinct.dedup();
    let atom_types = law.charge.len() - bis.is_some() as usize - distinct.len();
    let snap = Snapshot::load(snapshot, 3);
    let centres: Vec<usize> = match flag("--centres") {
        Some(c) => c.split(',').filter_map(|s| s.parse().ok()).collect(),
        None => {
            let every: usize = flag("--every").and_then(|s| s.parse().ok()).unwrap_or(4);
            (0..snap.molecules()).step_by(every).collect()
        }
    };
    let sites = |m: &[Vec3]| -> Vec<(Vec3, usize)> {
        let atoms: Vec<(Vec3, usize)> = m.iter().cloned().zip([0usize, 1, 1]).collect();
        let atoms = match bis {
            Some(d) => with_bisector_site(&atoms, atom_types, d),
            None => atoms,
        };
        with_frame_sites(&atoms, &frame)
    };
    let induction = |mols: &[&Vec<(Vec3, usize)>]| -> f64 {
        let mut cl = Cluster::default();
        for m in mols {
            cl.pol.push(m.iter().filter(|(_, t)| alpha[*t] > 0.0).map(|(p, t)| PolSite { pos: *p, alpha: alpha[*t] }).collect());
            cl.charges.push(m.iter().filter(|(_, t)| law.charge[*t] != 0.0).map(|(p, t)| Charge { pos: *p, q: law.charge[*t], sigma: law.sigma.get(*t).copied().unwrap_or(0.0) }).collect());
        }
        let mut links = Vec::new();
        for i in 0..mols.len() {
            for j in i + 1..mols.len() {
                links.push(Link { i, j, shift: Vec3::ZERO, weight: 1.0, dweight: 0.0, d: Vec3::ZERO });
            }
        }
        solve(&cl, &links, None, 1e-12).energy
    };
    // Per size: the many-body energy, its three-body part, the pair-sum energy.
    let mut rows: Vec<Vec<(f64, f64, f64)>> = vec![Vec::new(); max + 1];
    for &c in &centres {
        let mols: Vec<Vec<(Vec3, usize)>> = snap.cluster(c, max).iter().map(|(_, m)| sites(m)).collect();
        // The induction of every pair and every triple of the ten, once.
        let n = mols.len();
        let mut pair_ind = vec![vec![0.0; n]; n];
        let mut pair_other = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in i + 1..n {
                pair_ind[i][j] = induction(&[&mols[i], &mols[j]]);
                pair_other[i][j] = pair_energy(&law, &PairEnergy { a: mols[i].clone(), b: mols[j].clone(), energy: 0.0 });
            }
        }
        for size in 3..=max.min(n) {
            let refs: Vec<&Vec<(Vec3, usize)>> = mols[..size].iter().collect();
            let full = induction(&refs);
            let mut pairs = 0.0;
            let mut other = 0.0;
            for i in 0..size {
                for j in i + 1..size {
                    pairs += pair_ind[i][j];
                    other += pair_other[i][j];
                }
            }
            let mut three = 0.0;
            for i in 0..size {
                for j in i + 1..size {
                    for k in j + 1..size {
                        let e = induction(&[&mols[i], &mols[j], &mols[k]]);
                        three += e - pair_ind[i][j] - pair_ind[i][k] - pair_ind[j][k];
                    }
                }
            }
            rows[size].push(((full - pairs) * KCAL, three * KCAL, (other + pairs) * KCAL));
        }
    }
    println!("{} clusters per size about centres of {snapshot}; law {law_path}; kcal/mol", centres.len());
    println!("  n   many-body (mean +- s.e.)   sd     three-body    four-body up   pair-sum interaction   many-body / pair-sum");
    for size in 3..=max {
        let r = &rows[size];
        if r.is_empty() {
            continue;
        }
        let k = r.len() as f64;
        let mean = |f: &dyn Fn(&(f64, f64, f64)) -> f64| r.iter().map(f).sum::<f64>() / k;
        let m = mean(&|x| x.0);
        let sd = (r.iter().map(|x| (x.0 - m).powi(2)).sum::<f64>() / (k - 1.0).max(1.0)).sqrt();
        let t = mean(&|x| x.1);
        let p = mean(&|x| x.2);
        println!("{size:3}   {m:+9.3} +- {:5.3}   {sd:6.3}   {t:+9.3}   {:+9.3}   {p:+12.2}   {:+7.3}", sd / k.sqrt(), m - t, m / p);
    }
}
