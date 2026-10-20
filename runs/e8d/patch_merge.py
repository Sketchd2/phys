p = r'E:\Projects\phys\src\liquid.rs'
s = open(p, encoding='utf-8').read()

old = """        let (mut cpos, mut cq, mut ccom, mut carm) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for i in 0..n {
            for &k in &charged {
                cpos.push(self.com[i] + arms[i][k]);
                cq.push(q_of(k));
                ccom.push(self.com[i]);
                carm.push(arms[i][k]);
            }
        }
        let t_a = std::time::Instant::now();"""
new = """        // In reciprocal space every charge is a point, so charges at one place (a core and its cloud
        // on the same atom) are one charge there: the structure factor, and with it the energy, the
        // forces on the molecule and the virial, are unchanged, and there are fewer sites to sum.
        // The self term and the same-molecule terms below keep the original charges.
        let mut groups: Vec<(usize, f64)> = Vec::new(); // (a member site for the arm, the summed charge)
        for &k in &charged {
            match groups.iter_mut().find(|(k0, _)| (arms[0][*k0] - arms[0][k]).norm() < 1e-9) {
                Some(g) => g.1 += q_of(k),
                None => groups.push((k, q_of(k))),
            }
        }
        let (mut cpos, mut cq, mut ccom, mut carm) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for i in 0..n {
            for &(k, q) in &groups {
                cpos.push(self.com[i] + arms[i][k]);
                cq.push(q);
                ccom.push(self.com[i]);
                carm.push(arms[i][k]);
            }
        }
        let cq_original: Vec<f64> = (0..n).flat_map(|_| charged.iter().map(|&k| q_of(k))).collect();
        let t_a = std::time::Instant::now();"""
assert old in s
s = s.replace(old, new, 1)

old = "forces.energy += ew.energy(&sq) - ew.charge_self(&cq) - intra * n as f64;"
new = "forces.energy += ew.energy(&sq) - ew.charge_self(&cq_original) - intra * n as f64;"
assert old in s
s = s.replace(old, new, 1)

# the final application loop: one force per merged site
old = """        let mut idx = 0;
        for i in 0..n {
            for &k in &charged {
                let f = fch_rec[idx];
                forces.force[i] += f;
                forces.torque[i] += arms[i][k].cross(f);
                idx += 1;
            }
        }"""
new = """        let mut idx = 0;
        for i in 0..n {
            for &(k, _) in &groups {
                let f = fch_rec[idx];
                forces.force[i] += f;
                forces.torque[i] += arms[i][k].cross(f);
                idx += 1;
            }
        }"""
assert old in s, 'final loop'
s = s.replace(old, new, 1)
open(p, 'w', encoding='utf-8').write(s)
print('ok')
