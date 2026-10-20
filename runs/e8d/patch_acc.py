p = r'E:\Projects\phys\src\bin\ewald.rs'
s = open(p, encoding='utf-8').read()
s = s.replace("Ewald::new(sp, cell, 1e-7)", "Ewald::new(sp, cell, acc)")
s = s.replace("Ewald::new(sp, l.cell, 1e-7)", "Ewald::new(sp, l.cell, acc)")
s = s.replace("    let check = args.iter().any(|a| a == \"check\");", "    let check = args.iter().any(|a| a == \"check\");\n    let acc: f64 = std::env::var(\"PHYS_EWALD_ACC\").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-7);")
open(p, 'w', encoding='utf-8').write(s)

p = r'E:\Projects\phys\src\bin\bulk.rs'
s = open(p, encoding='utf-8').read()
s = s.replace("ewald: Ewald::new(sp, l.cell, 1e-7) };", "ewald: Ewald::new(sp, l.cell, std::env::var(\"PHYS_EWALD_ACC\").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-7)) };")
open(p, 'w', encoding='utf-8').write(s)
print('ok')
