//! The core-and-cloud model's six numbers as quadratics in the molecule's angle,
//! fitted to the exact electrostatic energies of pairs whose molecules each have
//! their own angle — `docs/PLAY.md` E8c, treatment (2).
//!
//! ```sh
//! cargo run --release --bin phys-geofit -- water cloud-fit.txt es-water.txt -gpu-bias-bend-mp2 \
//!     [--dipoles props.txt] [--theta0 104.078] [--ridge 0.5] [--hold-c0] [--out cloud-geo.txt]
//! ```
//!
//! `cloud-fit.txt` is the cloud at one angle (`phys-diag --fit-es` writes it): the
//! starting point and the reference (`theta0`). `es-water.txt` is
//! `phys-es-gpu`'s first-order electrostatic energies of the pairs in
//! `pairs-water<tag>.txt`. `--dipoles` is `phys-monomer-props`'s table (the
//! `scripts/e8c-monomer-props.sh` output): the molecule's MP2 dipole at the
//! angles where it was computed, which the model's dipole is held to, because the
//! densities the electrostatic energies come from are Hartree-Fock's and a few
//! per cent too polar. Only lines whose two bond lengths are equal and equal to the
//! cloud's hydrogen distance are used.
//!
//! Writes the cloud at `theta0` as an ordinary law's text (`charges`, `sigma`,
//! `bisector`, `site` lines, which `phys-fit --esp` reads) followed by the `geo`
//! lines that make it follow the angle (`phys::geolaw`).
//!
//! The fit is Levenberg-Marquardt on the residuals (weighted kcal/mol), half the
//! pairs held out, with a ridge on the slopes and curvatures so that parameters
//! the data cannot tell apart stay near the reference.

use phys::geolaw::{pair_electrostatics, placed, x_of_law, CloudParams, GeoCloud, NPAR};
use phys::liquid::{extra_sites_from_text, ExtraSite, SiteSite};
use phys::math::Vec3;

const KCAL: f64 = 627.509474;
const BOHR: f64 = 0.529177210903;

type Row = (Vec<(Vec3, usize)>, Vec<(Vec3, usize)>, f64);

fn molecule_dipole(geo: &GeoCloud, frame: &[ExtraSite], r_oh: f64, theta: f64) -> f64 {
    let atoms = vec![(Vec3::ZERO, 0usize), (Vec3 { x: r_oh, y: 0.0, z: 0.0 }, 1), (Vec3 { x: r_oh * theta.cos(), y: r_oh * theta.sin(), z: 0.0 }, 1)];
    let p = geo.params(theta);
    let sites = placed(&atoms, p.d, frame);
    let mu = sites.iter().fold(Vec3::ZERO, |acc, (x, t)| acc + x.scale(p.charge[*t]));
    mu.norm()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let flag = |n: &str| args.iter().position(|a| a == n).and_then(|i| args.get(i + 1)).cloned();
    let name = args[0].clone();
    let cloud_file = args[1].clone();
    let es_file = args[2].clone();
    let tags = args[3].clone();
    let theta0 = flag("--theta0").and_then(|v| v.parse::<f64>().ok()).unwrap_or(104.078).to_radians();
    let ridge: f64 = flag("--ridge").and_then(|v| v.parse().ok()).unwrap_or(0.5);
    let out = flag("--out").unwrap_or_else(|| "cloud-geo.txt".into());
    // `--hold-c0`: the cloud at the reference angle stays as it was given (a cloud already
    // fitted and tried), and only its slopes and curvatures in the angle are fitted: the
    // continuation from the nearest fitted geometry that the plan asks for.
    let hold_c0 = args.iter().any(|a| a == "--hold-c0");

    let text0 = std::fs::read_to_string(&cloud_file).unwrap_or_else(|_| panic!("no {cloud_file}"));
    let law0 = SiteSite::from_text(&text0).expect("a readable cloud law");
    assert_eq!(law0.charge.len(), 5, "a core-and-cloud law (phys-esp --cloud)");
    let bis0 = SiteSite::bisector_from_text(&text0).expect("a bisector line").1;
    let frame = extra_sites_from_text(&text0);
    let r_oh = frame.iter().find(|s| s.ty == 4).map(|s| s.frame[0]).expect("hydrogen clouds");
    let x0 = x_of_law(&law0, bis0);

    // The pairs and their exact electrostatics.
    let es: std::collections::HashMap<usize, f64> = std::fs::read_to_string(&es_file)
        .unwrap_or_else(|_| panic!("no {es_file}"))
        .lines()
        .filter_map(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            Some((w.first()?.parse().ok()?, w.get(2)?.parse().ok()?))
        })
        .collect();
    let text: String = tags.split(',').map(|t| std::fs::read_to_string(format!("pairs-{name}{t}.txt")).unwrap_or_default()).collect::<Vec<_>>().join("\n");
    let mut rows: Vec<Row> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#') && l.contains('|')) {
        let (head, tail) = line.split_once('|').expect("a |");
        let k: usize = head.split_whitespace().next().and_then(|x| x.parse().ok()).expect("an index");
        let Some(&reference) = es.get(&k) else { continue };
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<(Vec3, usize)> = t.chunks(4).map(|c| (Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }, c[3].parse().unwrap())).collect();
        let n = atoms.len() / 2;
        rows.push((atoms[..n].to_vec(), atoms[n..].to_vec(), reference * KCAL));
    }
    assert!(!rows.is_empty(), "no pairs with an electrostatic energy");

    // The dipole table: angle (rad) and MP2 dipole (au) at the bond length of the cloud's hydrogens.
    let mut dipoles: Vec<(f64, f64)> = Vec::new();
    if let Some(f) = flag("--dipoles") {
        for l in std::fs::read_to_string(&f).unwrap_or_else(|_| panic!("no {f}")).lines() {
            let (head, rest) = l.split_once('|').expect("a |");
            let h: Vec<f64> = head.split_whitespace().skip(1).filter_map(|x| x.parse().ok()).collect();
            let m: Vec<f64> = rest.split('|').next().unwrap().split_whitespace().filter_map(|x| x.parse().ok()).collect();
            if h.len() == 3 && (h[0] - r_oh).abs() < 1e-6 && (h[1] - r_oh).abs() < 1e-6 && m.len() >= 2 {
                dipoles.push((h[2], m[1]));
            }
        }
        println!("  {} MP2 dipoles at the bond length {:.4} A", dipoles.len(), r_oh * BOHR);
    }

    let geo_of = |z: &[f64]| -> GeoCloud {
        let mut g = GeoCloud::constant([law0.charge[0], law0.charge[1]], x0, theta0);
        for p in 0..NPAR {
            g.coef[p] = [z[3 * p], z[3 * p + 1], z[3 * p + 2]];
        }
        g
    };
    let model = |g: &GeoCloud, a: &[(Vec3, usize)], b: &[(Vec3, usize)]| -> f64 {
        let pa = g.params(GeoCloud::angle(&[a[0].0, a[1].0, a[2].0]));
        let pb = g.params(GeoCloud::angle(&[b[0].0, b[1].0, b[2].0]));
        pair_electrostatics(&placed(a, pa.d, &frame), &pa, &placed(b, pb.d, &frame), &pb) * KCAL
    };
    let weight = |r: &Row| -> f64 { 1.0 / (1.0 + ((r.0[0].0 - r.1[0].0).norm() * BOHR - 2.7).max(0.0).powi(2)) };
    let (train, test): (Vec<&Row>, Vec<&Row>) = {
        let (mut t, mut h) = (Vec::new(), Vec::new());
        for (i, r) in rows.iter().enumerate() {
            if i % 2 == 0 { t.push(r) } else { h.push(r) }
        }
        (t, h)
    };
    // Scales for the ridge: a slope or curvature of one is a change by one unit of the parameter over `scale` radians.
    let unit = [0.3, 0.3, 0.3, 0.2, 0.3, 0.3];
    let residuals = |z: &[f64]| -> Vec<f64> {
        let g = geo_of(z);
        let mut r: Vec<f64> = train.iter().map(|row| weight(row).sqrt() * (model(&g, &row.0, &row.1) - row.2)).collect();
        for &(th, mu) in &dipoles {
            r.push(30.0 * (molecule_dipole(&g, &frame, r_oh, th) - mu) / mu);
        }
        for p in 0..NPAR {
            for k in 1..3 {
                r.push(ridge * z[3 * p + k] / unit[p] * 4.0);
            }
        }
        r
    };
    let stat = |z: &[f64], set: &[&Row]| -> (f64, f64) {
        let g = geo_of(z);
        let n = set.len() as f64;
        (set.iter().map(|r| model(&g, &r.0, &r.1) - r.2).sum::<f64>() / n, (set.iter().map(|r| (model(&g, &r.0, &r.1) - r.2).powi(2)).sum::<f64>() / n).sqrt())
    };
    let report = |label: &str, z: &[f64]| {
        let (tm, tr) = stat(z, &train);
        let (hm, hr) = stat(z, &test);
        println!("  {label}: train mean/rms {tm:+.3}/{tr:.3}, held-out {hm:+.3}/{hr:.3} kcal/mol");
    };
    let mut z = vec![0.0; 3 * NPAR];
    for p in 0..NPAR {
        z[3 * p] = x0[p];
    }
    report("before (the cloud at the reference angle, no dependence)", &z);

    // Levenberg-Marquardt with a finite-difference Jacobian.
    let n = z.len();
    let mut lambda = 1e-2;
    let mut r = residuals(&z);
    let mut cost: f64 = r.iter().map(|x| x * x).sum();
    for _iter in 0..200 {
        let mut jac = vec![vec![0.0; n]; r.len()];
        for j in 0..n {
            if hold_c0 && j % 3 == 0 {
                continue;
            }
            let h = 1e-6 * (1.0 + z[j].abs());
            let mut zp = z.clone();
            zp[j] += h;
            let rp = residuals(&zp);
            for i in 0..r.len() {
                jac[i][j] = (rp[i] - r[i]) / h;
            }
        }
        let mut jtj = vec![vec![0.0; n]; n];
        let mut jtr = vec![0.0; n];
        for i in 0..r.len() {
            for a in 0..n {
                jtr[a] += jac[i][a] * r[i];
                for b in 0..n {
                    jtj[a][b] += jac[i][a] * jac[i][b];
                }
            }
        }
        let mut improved = false;
        for _try in 0..12 {
            let mut m = jtj.clone();
            for a in 0..n {
                m[a][a] += lambda * (jtj[a][a].max(1e-12));
            }
            let step = solve(&m, &jtr);
            let zn: Vec<f64> = z.iter().zip(&step).map(|(a, b)| a - b).collect();
            let rn = residuals(&zn);
            let cn: f64 = rn.iter().map(|x| x * x).sum();
            if cn < cost {
                let gain = cost - cn;
                z = zn;
                r = rn;
                cost = cn;
                lambda = (lambda / 3.0).max(1e-9);
                improved = gain > 1e-12 * cost.max(1e-12);
                break;
            }
            lambda *= 4.0;
        }
        if !improved {
            break;
        }
    }
    report("after ", &z);
    let g = geo_of(&z);
    println!("  parameters (c0, c1, c2) of ln sigma_bis, ln sigma_O, ln sigma_H, d_bis, n_bis, n_O:");
    for p in 0..NPAR {
        println!("    {p}: {:+.5} {:+.5} {:+.5}", g.coef[p][0], g.coef[p][1], g.coef[p][2]);
    }
    println!("  the cloud over the angle (degrees: charges bisector / oxygen / hydrogen, widths, d, dipole):");
    for deg in [88.0f64, 92.0, 96.0, 100.0, 104.078, 108.0, 112.0, 116.0, 120.0] {
        let th = deg.to_radians();
        let p: CloudParams = g.params(th);
        let mu = molecule_dipole(&g, &frame, r_oh, th);
        let table = dipoles.iter().min_by(|a, b| (a.0 - th).abs().partial_cmp(&(b.0 - th).abs()).unwrap()).filter(|d| (d.0 - th).abs() < 0.01).map(|d| format!(" (MP2 {:.4})", d.1)).unwrap_or_default();
        println!("    {deg:>7.3}: {:+.3} {:+.3} {:+.3}  {:.3} {:.3} {:.3}  {:.4}  {:.4} au{table}", p.charge[2], p.charge[3], p.charge[4], p.sigma[2], p.sigma[3], p.sigma[4], p.d, mu);
    }
    let base = g.params(theta0);
    let mut law = law0.clone();
    law.charge = base.charge.to_vec();
    law.sigma = base.sigma.to_vec();
    let text = format!(
        "charges {}\nsigma {}\nbisector 2 {:e}\nsite 3 0e0 0e0 0e0\nsite 4 {r_oh:e} 0e0 0e0\nsite 4 0e0 {r_oh:e} 0e0\n{}",
        law.charge.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "),
        law.sigma.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>().join(" "),
        base.d,
        g.to_text()
    );
    std::fs::write(&out, text).expect("the cloud written");
    println!("  written to {out}");
}

/// Solve `m x = b` for a small dense system by Gaussian elimination with pivoting.
fn solve(m: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut a: Vec<Vec<f64>> = m.to_vec();
    let mut x = b.to_vec();
    for c in 0..n {
        let piv = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs())).unwrap();
        a.swap(piv, c);
        x.swap(piv, c);
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            for k in c..n {
                a[r][k] -= f * a[c][k];
            }
            x[r] -= f * x[c];
        }
    }
    for r in (0..n).rev() {
        let mut s = x[r];
        for k in r + 1..n {
            s -= a[r][k] * x[k];
        }
        x[r] = s / a[r][r];
    }
    x
}
