p = r'E:\Projects\phys\src\bin\bulk.rs'
s = open(p, encoding='utf-8').read()

# import
old = "use phys::liquid::{Kind, Liquid, Polarisable, SiteLaw, SiteSite, AMU, BAR_PER_HARTREE_PER_BOHR3, K_B};"
new = "use phys::ewald::Ewald;\nuse phys::liquid::{EwaldLaw, Kind, Liquid, Polarisable, SiteLaw, SiteSite, AMU, BAR_PER_HARTREE_PER_BOHR3, K_B};"
assert old in s
s = s.replace(old, new, 1)

# wrap the law after the liquid is made
old = "    println!(\"{name} at {temperature} K: {n} molecules, {density_g} g/cm3, box {side:.2} bohr ({:.2} A), cut-off {r_cut:.1} bohr\", side * 0.529177);"
new = old + """
    // `PHYS_EWALD=<s' in bohr>` carries the long-range part of the charges and
    // the induced dipoles in reciprocal space (ewald.rs, PLAY.md E8d); the pair
    // terms keep only their short-range Coulomb remainder. Unset, the run is
    // the truncated one it has always been.
    let ewald_law;
    let law_ref: &dyn SiteLaw = match std::env::var("PHYS_EWALD").ok().and_then(|v| v.parse::<f64>().ok()) {
        Some(sp) => {
            ewald_law = EwaldLaw { inner: law_ref, ewald: Ewald::new(sp, l.cell, 1e-7) };
            println!("  Ewald sum on: s' = {sp} bohr, {} reciprocal vectors", ewald_law.ewald.nk());
            &ewald_law
        }
        None => law_ref,
    };"""
assert old in s
s = s.replace(old, new, 1)

# tail at the end: after the pressure line
old = "    println!(\"  temperatures: translation {:.1} K, rotation {:.1} K (bath {temperature} K)\", mean(2), mean(3));"
new = old + """
    // The dispersion beyond the cut-off, from the law's own C6 (isotropic, g = 1):
    // P = -2 pi rho^2 C6 I, I = integral from r_on of (1 - S) / r^4; and the same for U/N.
    {
        let types = law.charge.len();
        let mut c6 = 0.0;
        for ta in &l.kind.types {
            for tb in &l.kind.types {
                c6 += law.pair[ta * types + tb][2];
            }
        }
        let rho = n as f64 / (l.cell[0] * l.cell[1] * l.cell[2]);
        let (mut integral, mut r, h) = (0.0, l.r_on, 0.001);
        while r < 2000.0 {
            let (sw, _) = phys::liquid::switch(r + h / 2.0, l.r_on, l.r_cut);
            integral += (1.0 - sw) / (r + h / 2.0).powi(4) * h;
            r += h;
        }
        let p_tail = -2.0 * std::f64::consts::PI * rho * rho * c6 * integral * BAR_PER_HARTREE_PER_BOHR3;
        let u_tail = -2.0 * std::f64::consts::PI * rho * c6 * integral * KCAL;
        println!("  dispersion beyond the cut-off (C6 {c6:.2}): {p_tail:.0} bar and {u_tail:.3} kcal/mol a molecule, so pressure with it {:.0} bar, U/N {:.3}", mean(1) + p_tail, mean(0) + u_tail);
    }"""
assert old in s
s = s.replace(old, new, 1)
open(p, 'w', encoding='utf-8').write(s)
print('patched')
