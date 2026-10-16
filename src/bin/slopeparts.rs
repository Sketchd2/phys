//! The law's slope along the line of centres on the pairs the engine has MP2
//! slopes for, split into its terms (`docs/PLAY.md` E8d): which term carries the
//! force in a shell, and so which could be the one that is off.
//!
//! ```sh
//! phys-slopeparts law.txt deriv-water.txt pairs-tag      # pairs-water<tag>.txt in the working directory
//! ```
//!
//! kcal/mol per bohr, by oxygen-oxygen separation: MP2's slope, the law's, and the law's
//! permanent electrostatics, pair induction, repulsion and dispersion.

use phys::liquid::{extra_sites_from_text, pair_energy, with_bisector_site, with_frame_sites, PairEnergy, Polarisable, SiteSite};
use phys::math::Vec3;

const KCAL: f64 = 627.509474;
const BOHR: f64 = 0.529177210903;

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
    let types = law.charge.len();
    let atom_types = types - bis.is_some() as usize - distinct.len();
    let geometry = std::fs::read_to_string(format!("pairs-water{}.txt", args[2])).unwrap_or_else(|_| panic!("no pairs-water{}.txt", args[2]));
    let mut by_index: std::collections::HashMap<usize, Vec<(Vec3, usize)>> = std::collections::HashMap::new();
    for line in geometry.lines().filter(|l| !l.starts_with('#') && l.contains('|')) {
        let (head, tail) = line.split_once('|').unwrap();
        let k: usize = head.split_whitespace().next().and_then(|x| x.parse().ok()).unwrap();
        let t: Vec<&str> = tail.split_whitespace().collect();
        by_index.insert(k, t.chunks(4).map(|c| (Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }, c[3].parse().unwrap())).collect());
    }
    let mass = |t: usize| if t == 0 { 15.999 } else { 1.008 };
    let mut no_disp = law.clone();
    no_disp.pair.iter_mut().for_each(|p| {
        p[2] = 0.0;
        p[3] = 0.0;
    });
    let mut no_rep = law.clone();
    no_rep.pair.iter_mut().for_each(|p| p[0] = 0.0);
    no_rep.overlap.iter_mut().for_each(|k| *k = 0.0);
    let polar = Polarisable { law: &law, alpha: alpha.clone() };
    // rows: O-O (A), MP2, total, electrostatics, induction, repulsion, dispersion
    let mut rows: Vec<[f64; 7]> = Vec::new();
    for line in std::fs::read_to_string(&args[1]).unwrap().lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        if w.len() < 4 {
            continue;
        }
        let (k, delta, ep, em): (usize, f64, f64, f64) = (w[0].parse().unwrap(), w[1].parse().unwrap(), w[2].parse().unwrap(), w[3].parse().unwrap());
        let Some(atoms) = by_index.get(&k) else { continue };
        let n = atoms.len() / 2;
        let (a, b) = (&atoms[..n], &atoms[n..]);
        let com = |m: &[(Vec3, usize)]| -> Vec3 {
            let total: f64 = m.iter().map(|(_, t)| mass(*t)).sum();
            m.iter().fold(Vec3::ZERO, |acc, (p, t)| acc + p.scale(mass(*t) / total))
        };
        let axis = (com(b) - com(a)).unit();
        let place = |m: &[(Vec3, usize)]| with_frame_sites(&with_bisector_site(m, atom_types, bis.unwrap_or(0.0)), &frame);
        let shifted = |sgn: f64| -> Vec<(Vec3, usize)> { b.iter().map(|(p, t)| (*p + axis.scale(sgn * delta), *t)).collect() };
        let slope = |f: &dyn Fn(&PairEnergy) -> f64| -> f64 {
            let e = |bb: Vec<(Vec3, usize)>| f(&PairEnergy { a: place(a), b: place(&bb), energy: 0.0 });
            (e(shifted(1.0)) - e(shifted(-1.0))) / (2.0 * delta) * KCAL
        };
        let full = slope(&|d| pair_energy(&polar, d));
        let base = slope(&|d| pair_energy(&law, d));
        let nd = slope(&|d| pair_energy(&no_disp, d));
        let nr = slope(&|d| pair_energy(&no_rep, d));
        let disp = base - nd;
        let rep = base - nr;
        let ind = full - base;
        let es = base - disp - rep;
        let o_o = (a[0].0 - b[0].0).norm() * BOHR;
        if std::env::var("VERBOSE").is_ok() && o_o > 3.0 && o_o < 3.9 {
            // the shortest H...O between the molecules, and the H-bond angle there
            let mut best = (f64::INFINITY, 0.0);
            for (ha, oa) in [(a, b), (b, a)] {
                for h in 1..3 {
                    let d = (ha[h].0 - oa[0].0).norm() * BOHR;
                    if d < best.0 {
                        let (u, v) = ((ha[0].0 - ha[h].0).unit(), (oa[0].0 - ha[h].0).unit());
                        best = (d, u.dot(v).clamp(-1.0, 1.0).acos().to_degrees());
                    }
                }
            }
            println!("  pair {k:4} O-O {o_o:.2}  shortest H..O {:.2} A at O-H..O {:.0} deg | MP2 slope {:+.3} law {:+.3} err {:+.3} | es {:+.2} ind {:+.2} rep {:+.2} disp {:+.2}", best.0, best.1, (ep - em) / (2.0 * delta) * KCAL, full, full - (ep - em) / (2.0 * delta) * KCAL, es, ind, rep, disp);
        }
        rows.push([o_o, (ep - em) / (2.0 * delta) * KCAL, full, es, ind, rep, disp]);
    }
    println!("slopes along the line of centres, kcal/mol per bohr, mean over the pairs in each O-O shell: {} pairs", rows.len());
    println!("  O-O (A)      n     MP2    law   law-MP2 | electrostatics  induction  repulsion  dispersion");
    for w in [0.0, 2.8, 3.1, 3.4, 3.8, 4.7].windows(2) {
        let sel: Vec<&[f64; 7]> = rows.iter().filter(|r| r[0] >= w[0] && r[0] < w[1]).collect();
        if sel.is_empty() {
            continue;
        }
        let m = |k: usize| sel.iter().map(|r| r[k]).sum::<f64>() / sel.len() as f64;
        println!("  {:.1}-{:.1}  {:>4}  {:>7.3} {:>7.3} {:>7.3} | {:>10.3} {:>10.3} {:>10.3} {:>10.3}", w[0], w[1], sel.len(), m(1), m(2), m(2) - m(1), m(3), m(4), m(5), m(6));
    }
}
