//! The pair virial of one liquid snapshot, term by term and shell by shell
//! (PLAY.md E8d).
//!
//! ```sh
//! phys-virial law.txt bulk-water-298.snap [temperature K] [r_cut bohr] [r_on bohr]
//! ```
//!
//! For every pair of molecules within the cut-off the law's energy is taken
//! as the molecules' centres are moved along their line (orientations fixed),
//! which for rigid molecules is minus the force on the centres, and
//! `-d dU/ds` summed over pairs, over three times the volume, is the pair
//! virial's share of the pressure. The terms are separated by switching them
//! off in turn - dispersion, repulsion, induction (a pair's, not the liquid's
//! many-body solve), and the permanent electrostatics are what is left - and
//! the pressure is the translational kinetic term plus these plus the
//! dispersion beyond the cut-off from the law's own C6. It is one snapshot, so
//! it says where the pressure comes from and how it would answer a change in
//! one term; it is not the run's average.

use phys::liquid::{extra_sites_from_text, pair_energy, switch, with_bisector_site, with_frame_sites, PairEnergy, Polarisable, SiteSite, BAR_PER_HARTREE_PER_BOHR3, K_B};
use phys::math::Vec3;

const BOHR: f64 = 0.529177210903;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let law_path = args.first().expect("a law file");
    let snap_path = args.get(1).expect("a snapshot");
    let temperature: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(298.0);
    let r_cut: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(16.0);
    let r_on: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(r_cut - 2.0);
    let text = std::fs::read_to_string(law_path).unwrap_or_else(|_| panic!("no {law_path}"));
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    let bis = SiteSite::bisector_from_text(&text).map(|b| b.1);
    let frame = extra_sites_from_text(&text);
    let mut distinct: Vec<usize> = frame.iter().map(|s| s.ty).collect();
    distinct.sort();
    distinct.dedup();
    let types = law.charge.len();
    let atom_types = types - bis.is_some() as usize - distinct.len();

    let snap = std::fs::read_to_string(snap_path).unwrap_or_else(|_| panic!("no {snap_path}"));
    let mut cell = [0.0; 3];
    let mut mols: Vec<Vec<Vec3>> = Vec::new();
    for line in snap.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w.first().copied() {
            Some("box") => (0..3).for_each(|k| cell[k] = w[1 + k].parse().unwrap()),
            Some("mol") => {
                let v: Vec<f64> = w[2..].iter().map(|x| x.parse().unwrap()).collect();
                mols.push(v.chunks(3).map(|c| Vec3 { x: c[0], y: c[1], z: c[2] }).collect());
            }
            _ => {}
        }
    }
    let n = mols.len();
    let volume = cell[0] * cell[1] * cell[2];
    let mass = |k: usize| if k == 0 { 15.999 } else { 1.008 };
    let com = |m: &[Vec3]| -> Vec3 {
        let total: f64 = (0..m.len()).map(mass).sum();
        m.iter().enumerate().fold(Vec3::ZERO, |acc, (k, p)| acc + p.scale(mass(k) / total))
    };
    let image = |mut d: Vec3| {
        d.x -= cell[0] * (d.x / cell[0]).round();
        d.y -= cell[1] * (d.y / cell[1]).round();
        d.z -= cell[2] * (d.z / cell[2]).round();
        d
    };
    let sites = |m: &[Vec3]| -> Vec<(Vec3, usize)> {
        let atoms: Vec<(Vec3, usize)> = m.iter().cloned().zip([0usize, 1, 1]).collect();
        let atoms = match bis {
            Some(d) => with_bisector_site(&atoms, atom_types, d),
            None => atoms,
        };
        with_frame_sites(&atoms, &frame)
    };

    // The law with terms off.
    let mut no_disp = law.clone();
    no_disp.pair.iter_mut().for_each(|p| {
        p[2] = 0.0;
        p[3] = 0.0;
    });
    let mut no_rep = law.clone();
    no_rep.pair.iter_mut().for_each(|p| p[0] = 0.0);
    no_rep.overlap.iter_mut().for_each(|k| *k = 0.0);
    let energy = |l: &SiteSite, ind: bool, a: &[(Vec3, usize)], b: &[(Vec3, usize)]| -> f64 {
        let pe = PairEnergy { a: a.to_vec(), b: b.to_vec(), energy: 0.0 };
        if ind {
            pair_energy(&Polarisable { law: l, alpha: alpha.clone() }, &pe)
        } else {
            pair_energy(l, &pe)
        }
    };

    {
        use phys::induction::{molecular_polarisability, PolSite};
        let m = sites(&mols[0]);
        let pol: Vec<PolSite> = m.iter().filter(|(_, t)| alpha[*t] > 0.0).map(|(p, t)| PolSite { pos: *p, alpha: alpha[*t] }).collect();
        let t = molecular_polarisability(&pol);
        let sum_sites: f64 = pol.iter().map(|p| p.alpha).sum();
        println!("molecule's polarisability from the law's sites: mean {:.3} bohr^3 (tensor diagonal {:.3} {:.3} {:.3}); sites' sum {:.3}", (t[0][0] + t[1][1] + t[2][2]) / 3.0, t[0][0], t[1][1], t[2][2], sum_sites);
    }
    let delta = 0.02;
    let edges = [0.0, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.5];
    // [total, dispersion, repulsion, induction, electrostatics] per shell.
    let mut w = vec![[0.0f64; 5]; edges.len() - 1];
    let mut counts = vec![0usize; edges.len() - 1];
    let mut energy_sum = 0.0f64;
    let mut ind_sum = 0.0f64;
    let coms: Vec<Vec3> = mols.iter().map(|m| com(m)).collect();
    for i in 0..n {
        for j in i + 1..n {
            let d = image(coms[j] - coms[i]);
            let dist = d.norm();
            if dist >= r_cut {
                continue;
            }
            let shift = d - (coms[j] - coms[i]);
            let a = sites(&mols[i]);
            let b0: Vec<Vec3> = mols[j].iter().map(|p| *p + shift).collect();
            let axis = d.unit();
            let at = |s: f64| -> Vec<(Vec3, usize)> { sites(&b0.iter().map(|p| *p + axis.scale(s)).collect::<Vec<_>>()) };
            let (bp, bm, b) = (at(delta), at(-delta), at(0.0));
            let (sw, dsw) = switch(dist, r_on, r_cut);
            let term = |l: &SiteSite, ind: bool| -> f64 {
                let u0 = energy(l, ind, &a, &b);
                let du = (energy(l, ind, &a, &bp) - energy(l, ind, &a, &bm)) / (2.0 * delta);
                // d(S u)/ds, then -d times it: the pair's share of the virial.
                -dist * (sw * du + dsw * u0)
            };
            let (e_full, e_base) = (energy(&law, true, &a, &b), energy(&law, false, &a, &b));
            energy_sum += sw * e_full;
            ind_sum += sw * (e_full - e_base);
            let full = term(&law, true);
            let base = term(&law, false);
            let nd = term(&no_disp, false);
            let nr = term(&no_rep, false);
            let disp = base - nd;
            let rep = base - nr;
            let ind = full - base;
            let coul = base - disp - rep;
            let shell = edges.windows(2).position(|e| dist >= e[0] && dist < e[1]).unwrap_or(edges.len() - 2);
            for (k, v) in [full, disp, rep, ind, coul].into_iter().enumerate() {
                w[shell][k] += v;
            }
            counts[shell] += 1;
        }
    }

    // The induction the liquid actually solves: every molecule's dipoles together,
    // links switched, as `Liquid::add_induction` does. Its energy at a scaling of
    // the molecules' centres (each molecule rigid) gives the many-body virial.
    let many_body = |scale: f64| -> f64 {
        use phys::induction::{solve, Charge, Cluster, Link, PolSite};
        let shifted: Vec<Vec<Vec3>> = mols.iter().zip(&coms).map(|(m, c)| m.iter().map(|p| *p + c.scale(scale - 1.0)).collect()).collect();
        let all: Vec<Vec<(Vec3, usize)>> = shifted.iter().map(|m| sites(m)).collect();
        let mut cl = Cluster::default();
        for m in &all {
            cl.pol.push(m.iter().filter(|(_, t)| alpha[*t] > 0.0).map(|(p, t)| PolSite { pos: *p, alpha: alpha[*t] }).collect());
            cl.charges.push(m.iter().filter(|(_, t)| law.charge[*t] != 0.0).map(|(p, t)| Charge { pos: *p, q: law.charge[*t], sigma: law.sigma.get(*t).copied().unwrap_or(0.0) }).collect());
        }
        let mut links = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                // The links are those of the unscaled box; their separations scale.
                let raw0 = coms[i] - coms[j];
                let d0 = image(raw0);
                if d0.norm() >= r_cut {
                    continue;
                }
                let raw = raw0.scale(scale);
                let d = d0.scale(scale);
                let (wt, dw) = switch(d.norm(), r_on, r_cut);
                links.push(Link { i, j, shift: raw - d, weight: wt, dweight: dw, d });
            }
        }
        solve(&cl, &links, None, 1e-11).energy
    };
    let (e_minus, e_zero, e_plus) = (many_body(0.999), many_body(1.0), many_body(1.001));
    let ind_virial_energy = -(e_plus - e_minus) / (2.0 * 0.001);
    println!("induction, many-body solve: {:.3} kcal/mol per molecule (pair sum: {:.3}); non-additive {:.3}", e_zero / n as f64 * 627.509474, ind_sum / n as f64 * 627.509474, (e_zero - ind_sum) / n as f64 * 627.509474);
    let to_bar = BAR_PER_HARTREE_PER_BOHR3 / (3.0 * volume);
    println!("{n} molecules, box {:.2} bohr, {temperature} K, cut-off {r_on}-{r_cut} bohr, law {law_path}", cell[0]);
    println!("pair virial by centre-to-centre shell, in bar: total  dispersion  repulsion  pair-induction  permanent-electrostatics");
    let mut sum = [0.0f64; 5];
    for (s, e) in edges.windows(2).enumerate() {
        println!("  {:5.2}-{:5.2} A  n={:5}  {:9.0} {:9.0} {:9.0} {:9.0} {:9.0}", e[0] * BOHR, e[1] * BOHR, counts[s], w[s][0] * to_bar, w[s][1] * to_bar, w[s][2] * to_bar, w[s][3] * to_bar, w[s][4] * to_bar);
        for k in 0..5 {
            sum[k] += w[s][k];
        }
    }
    let kinetic = n as f64 * K_B * temperature / volume * BAR_PER_HARTREE_PER_BOHR3;
    println!("  all pairs               {:9.0} {:9.0} {:9.0} {:9.0} {:9.0}", sum[0] * to_bar, sum[1] * to_bar, sum[2] * to_bar, sum[3] * to_bar, sum[4] * to_bar);
    // The dispersion beyond the cut-off, from the law's own C6.
    let mut c6 = 0.0;
    for a in [0usize, 1, 1] {
        for b in [0usize, 1, 1] {
            c6 += law.pair[a * types + b][2];
        }
    }
    let rho = n as f64 / volume;
    let (mut integral, mut r, h) = (0.0, r_on, 0.001);
    while r < 2000.0 {
        let (s, _) = switch(r + h / 2.0, r_on, r_cut);
        integral += (1.0 - s) / (r + h / 2.0).powi(4) * h;
        r += h;
    }
    let tail = -2.0 * std::f64::consts::PI * rho * rho * c6 * integral * BAR_PER_HARTREE_PER_BOHR3;
    let tail_placeholder = tail;
    println!("pair energy per molecule (switched, pair induction): {:.3} kcal/mol", energy_sum / n as f64 * 627.509474);
    println!("translational kinetic term N kT / V      {kinetic:9.0} bar");
    println!("pair virial (all terms, pair induction)  {:9.0} bar", sum[0] * to_bar);
    println!("pressure of the truncated pair model     {:9.0} bar   (many-body induction not included)", kinetic + sum[0] * to_bar);
    // W = -dE/dln(lambda); the pair part of induction is already in `sum`, so the
    // many-body virial's extra is the solve's less the pairs'.
    let w_ind_total = ind_virial_energy;
    let w_ind_pairs = sum[3];
    let p_mb = (w_ind_total - w_ind_pairs) * to_bar;
    println!("many-body induction virial beyond the pairs' {p_mb:9.0} bar  (solve {:.0}, pairs {:.0})", w_ind_total * to_bar, w_ind_pairs * to_bar);
    println!("PRESSURE (kinetic + pairs + many-body induction){:9.0} bar   <- compare the run", kinetic + sum[0] * to_bar + p_mb);
    println!("                           with the tail beyond the cut-off{:9.0} bar", kinetic + sum[0] * to_bar + p_mb + tail_placeholder);
    println!("dispersion beyond the cut-off (C6 {c6:.2})  {tail:9.0} bar  ->  with the tail {:9.0} bar", kinetic + sum[0] * to_bar + tail);
}
