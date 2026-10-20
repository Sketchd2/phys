use phys::geolaw::{correction, GeoCloud};
use phys::liquid::{extra_sites_from_text, pair_energy, with_bisector_site, with_frame_sites, PairEnergy, Polarisable, SiteSite};
use phys::math::Vec3;

// A law file against a pairs file: by the nearest-atom distance, the MP2 energy, the law's, and the error.
// LAW, PAIRS, optionally GEO (a file with geo lines) and ALPHA (not needed: read from the law).
#[test]
fn wallcheck() {
    let law_file = std::env::var("LAW").unwrap();
    let pairs_file = std::env::var("PAIRS").unwrap();
    let text = std::fs::read_to_string(&law_file).unwrap();
    let law = SiteSite::from_text(&text).unwrap();
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    let bis = SiteSite::bisector_from_text(&text).unwrap();
    let frame = extra_sites_from_text(&text);
    let geo = GeoCloud::from_text(&text);
    let place = |m: &[(Vec3, usize)]| with_frame_sites(&with_bisector_site(m, bis.0, bis.1), &frame);
    let mut rows: Vec<(f64, f64, f64)> = Vec::new();
    for l in std::fs::read_to_string(&pairs_file).unwrap().lines() {
        if l.starts_with('#') || !l.contains('|') {
            continue;
        }
        let (head, tail) = l.split_once('|').unwrap();
        let e: f64 = head.split_whitespace().nth(2).unwrap().parse().unwrap();
        let v: Vec<f64> = tail.split_whitespace().map(|x| x.parse().unwrap()).collect();
        let at = |i: usize| (Vec3 { x: v[4 * i], y: v[4 * i + 1], z: v[4 * i + 2] }, v[4 * i + 3] as usize);
        let a: Vec<(Vec3, usize)> = (0..3).map(at).collect();
        let b: Vec<(Vec3, usize)> = (3..6).map(at).collect();
        let pe = PairEnergy { a: place(&a), b: place(&b), energy: 0.0 };
        let mut m = pair_energy(&Polarisable { law: &law, alpha: alpha.clone() }, &pe);
        if let Some(g) = &geo {
            m += correction(g, &frame, &alpha, &a, &b);
        }
        let near = a.iter().flat_map(|x| b.iter().map(move |y| (x.0 - y.0).norm())).fold(f64::INFINITY, f64::min) * 0.529177210903;
        rows.push((near, e * 627.509474, m * 627.509474));
    }
    println!("{} pairs of {pairs_file} against {law_file}", rows.len());
    for (lo, hi) in [(0.6, 0.8), (0.8, 1.0), (1.0, 1.3), (1.3, 1.6), (1.6, 1.9), (1.9, 2.3), (2.3, 9.0)] {
        let sel: Vec<&(f64, f64, f64)> = rows.iter().filter(|r| r.0 >= lo && r.0 < hi).collect();
        if sel.is_empty() {
            continue;
        }
        let n = sel.len() as f64;
        let mean_ref = sel.iter().map(|r| r.1).sum::<f64>() / n;
        let mean_law = sel.iter().map(|r| r.2).sum::<f64>() / n;
        let rms = (sel.iter().map(|r| (r.2 - r.1).powi(2)).sum::<f64>() / n).sqrt();
        let rel = (sel.iter().map(|r| ((r.2 - r.1) / r.1.abs().max(5.0)).powi(2)).sum::<f64>() / n).sqrt();
        println!("  nearest atoms {lo:.1}-{hi:.1} A: {:>3} pairs, MP2 mean {mean_ref:>8.1}, law mean {mean_law:>8.1}, rms error {rms:>7.1} kcal/mol, relative {:.0}%", sel.len(), rel * 100.0);
    }
}
