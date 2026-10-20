p = r"E:\Projects\phys\src\flex.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")

# 1. the links, built once, and the engine call
old = '        let job = |idx: &mut dyn Iterator<Item = usize>| {\n            idx.map(|i| {\n                let mut fsite = vec![Vec3::ZERO; ns];\n                let (mut gq, mut gs)'
assert old in s
new = '''        // The pairs of molecules inside the cut-off, for the induction and for a pair-loop engine.
        let pair_engine = crate::liquid::pair_engine().filter(|_| law.site_site().is_some());
        let mut links: Vec<Link> = Vec::new();
        if law.induction().is_some() || pair_engine.is_some() {
            for i in 0..n {
                for j in i + 1..n {
                    let raw = centres[i] - centres[j];
                    let d = self.image(raw);
                    let dist = d.norm();
                    if dist >= self.r_cut {
                        continue;
                    }
                    let (wgt, dw) = switch(dist, self.r_on, self.r_cut);
                    links.push(Link { i, j, shift: raw - d, weight: wgt, dweight: dw, d });
                }
            }
        }
        let engine_result = pair_engine.as_ref().and_then(|e| {
            let tab = law.site_site()?;
            let rel: Vec<Vec<Vec3>> = (0..n).map(|i| sites[i].iter().map(|x| *x - centres[i]).collect()).collect();
            let (qs, ss): (Vec<Vec<f64>>, Vec<Vec<f64>>) = match &geo_params {
                Some(gp) => ((0..n).map(|i| types.iter().map(|&t| gp[i].charge[t]).collect()).collect(), (0..n).map(|i| types.iter().map(|&t| gp[i].sigma[t]).collect()).collect()),
                None => ((0..n).map(|_| types.iter().map(|&t| tab.charge[t]).collect()).collect(), (0..n).map(|_| types.iter().map(|&t| tab.sigma.get(t).copied().unwrap_or(0.0)).collect()).collect()),
            };
            e.pairs(&crate::liquid::PairProblem { table: tab, types, rel: &rel, charge: &qs, sigma: &ss, links: &links })
        });
        let job = |idx: &mut dyn Iterator<Item = usize>| {
            idx.map(|i| {
                let mut fsite = vec![Vec3::ZERO; ns];
                let (mut gq, mut gs)'''
s = s.replace(old, new, 1)

# 2. use the engine's answer, else the CPU loop
old2 = '''        for part in crate::electrons::scf::parallel_interleaved(n, &job) {
            for (i, e, fsite, fc, w, q, sg) in part {
                energy += e;
                virial += 0.5 * w;'''
assert old2 in s
new2 = '''        let parts = match engine_result {
            Some(r) => {
                energy = r.energy;
                virial = r.virial;
                vec![(0..n).map(|i| (i, 0.0, r.fsite[i].clone(), r.fcentre[i], 0.0, r.gq[i].clone(), r.gs[i].clone())).collect::<Vec<_>>()]
            }
            None => crate::electrons::scf::parallel_interleaved(n, &job),
        };
        for part in parts {
            for (i, e, fsite, fc, w, q, sg) in part {
                energy += e;
                virial += 0.5 * w;'''
s = s.replace(old2, new2, 1)

# 3. the induction block uses the links built above
old3 = '''            let mut links = Vec::new();
            for i in 0..n {
                for j in i + 1..n {
                    let raw = centres[i] - centres[j];
                    let d = self.image(raw);
                    let dist = d.norm();
                    if dist >= self.r_cut {
                        continue;
                    }
                    let (wgt, dw) = switch(dist, self.r_on, self.r_cut);
                    links.push(Link { i, j, shift: raw - d, weight: wgt, dweight: dw, d });
                }
            }
            let out = solve(&cl, &links, warm, 1e-9);'''
assert old3 in s
s = s.replace(old3, "            let out = solve(&cl, &links, warm, 1e-9);", 1)
open(p, "w", newline="").write(s.replace("\n", nl))
print("flex.rs patched")
