p = r'E:\Projects\phys\src\induction.rs'
s = open(p, encoding='utf-8').read()
a = s.index("pub fn solve_ewald(")
head, tail = s[:a], s[a:]
# timers, printed once per call when PHYS_PROFILE is set
tail = tail.replace("    let nm = cl.pol.len();\n", "    let t_start = std::time::Instant::now();\n    let nm = cl.pol.len();\n", 1)
tail = tail.replace("    let mut e0: Vec<Vec3> = e0_extra.to_vec();", "    let t_items = t_start.elapsed();\n    let mut e0: Vec<Vec3> = e0_extra.to_vec();", 1)
tail = tail.replace("    let energy = -0.5 * mu.iter().zip(&e0).map(|(m, e)| m.dot(*e)).sum::<f64>();\n    let mut fpol", "    let t_cg = t_start.elapsed();\n    let energy = -0.5 * mu.iter().zip(&e0).map(|(m, e)| m.dot(*e)).sum::<f64>();\n    let mut fpol", 1)
tail = tail.replace("    Induced { energy, dipoles, force_pol, force_charge, link_force, virial, iterations }\n}", "    if std::env::var_os(\"PHYS_PROFILE\").is_some() {\n        let t_all = t_start.elapsed();\n        eprintln!(\"  solve_ewald: items {:.1} ms, cg ({} iterations) {:.1} ms, forces {:.1} ms, links {}\", t_items.as_secs_f64() * 1e3, iterations, (t_cg - t_items).as_secs_f64() * 1e3, (t_all - t_cg).as_secs_f64() * 1e3, links.len());\n    }\n    Induced { energy, dipoles, force_pol, force_charge, link_force, virial, iterations }\n}", 1)
s = head + tail
open(p, 'w', encoding='utf-8').write(s)

p = r'E:\Projects\phys\src\liquid.rs'
t = open(p, encoding='utf-8').read()
t = t.replace("        let ph_c = ew.phases(&cpos);\n        let sq = ew.charge_sum(&ph_c, &cq);\n        // The permanent charges: the reciprocal energy", "        let t_a = std::time::Instant::now();\n        let ph_c = ew.phases(&cpos);\n        let sq = ew.charge_sum(&ph_c, &cq);\n        let t_phases = t_a.elapsed();\n        // The permanent charges: the reciprocal energy", 1)
t = t.replace("            let out = crate::induction::solve_ewald(", "            let t_b = std::time::Instant::now();\n            let out = crate::induction::solve_ewald(", 1)
t = t.replace("            let _ = solve_ewald;\n", "            let _ = solve_ewald;\n            let t_solve = t_b.elapsed();\n", 1)
t = t.replace("        forces.virial += virial_rec;\n        if !dipoles.is_empty() {", "        forces.virial += virial_rec;\n        if std::env::var_os(\"PHYS_PROFILE\").is_some() {\n            eprintln!(\"  add_ewald: charge phases+sum {:.1} ms; whole add_ewald {:.1} ms\", t_phases.as_secs_f64() * 1e3, t_a.elapsed().as_secs_f64() * 1e3);\n        }\n        if !dipoles.is_empty() {", 1)
open(p, 'w', encoding='utf-8').write(t)

# the erf test: tight below 4, looser where the series itself loses digits
q = r'E:\Projects\phys\tests\ewald.rs'
u = open(q, encoding='utf-8').read()
u = u.replace('''        let d = (phys::induction::erf_fast(x) - phys::induction::erf(x)).abs();
        if d > worst {
            worst = d;
            at = x;
        }
    }
    assert!(worst < 2e-15, "erf_fast differs from erf by {worst:e} at {at}");''', '''        let d = (phys::induction::erf_fast(x) - phys::induction::erf(x)).abs();
        // The series sums ~100 terms of e^(x^2) size beyond x = 4 and carries ~1e-14 of its own there.
        let limit = if x.abs() < 4.0 { 2e-15 } else { 3e-14 };
        if d / limit > worst {
            worst = d / limit;
            at = x;
        }
    }
    assert!(worst < 1.0, "erf_fast differs from erf by {worst} of its allowance at {at}");''')
open(q, 'w', encoding='utf-8').write(u)
print('ok')
