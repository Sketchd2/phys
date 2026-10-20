p = r"E:\Projects\phys\src\bin\fit.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
marker = "    let n_atoms = rows[0].2.len() / 2;\n"
assert marker in s
block = '''    // `--geo-correct cloud-geo.txt`: molecules whose cloud follows their angle
    // (`phys::geolaw`, written by `phys-geofit`). This fit's model has the cloud
    // at the reference angle, so the part of every energy that the angles change
    // - electrostatics and induction with each molecule's own cloud, less the
    // same with the reference cloud - is taken off the target, and the terms
    // fitted here (repulsion, dispersion) are what is left for them. Pairs,
    // slopes and turns alike.
    if let Some(gf) = args.iter().position(|a| a == "--geo-correct").and_then(|i| args.get(i + 1)) {
        let gt = std::fs::read_to_string(gf).unwrap_or_else(|_| panic!("no {gf}"));
        let geo = phys::geolaw::GeoCloud::from_text(&gt).expect("geo lines in the --geo-correct file");
        let frame = extra_sites_from_text(&gt);
        let alpha = alpha_atoms.clone().expect("--geo-correct needs --alpha (the induction is corrected too)");
        let corr = |a: &[(Vec3, usize)], b: &[(Vec3, usize)]| phys::geolaw::correction(&geo, &frame, &alpha, a, b);
        let mut worst: f64 = 0.0;
        for (ep, er, atoms) in rows.iter_mut() {
            let n = atoms.len() / 2;
            let c = corr(&atoms[..n], &atoms[n..]);
            worst = worst.max(c.abs() * KCAL);
            *ep -= c;
            *er -= c;
        }
        for sl in raw_slopes.iter_mut() {
            let (cp, cm) = (corr(&sl.plus.a, &sl.plus.b), corr(&sl.minus.a, &sl.minus.b));
            sl.plus.energy -= cp;
            sl.minus.energy -= cm;
            sl.target = sl.plus.energy - sl.minus.energy;
        }
        println!("  angle-following cloud taken off {} pairs and {} slopes and turns (largest correction {worst:.3} kcal/mol)", rows.len(), raw_slopes.len());
    }
'''
s = s.replace(marker, block + marker, 1)
open(p, "w", newline="").write(s.replace("\n", nl))
print("patched")
