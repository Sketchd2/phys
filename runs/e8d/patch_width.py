p = r'E:\Projects\phys\src\induction.rs'
s = open(p, encoding='utf-8').read()
old = """/// The width of a smeared dipole of polarisability `alpha`.
pub fn sigma_of(alpha: f64) -> f64 {
    ((2.0 / PI).sqrt() * alpha / 3.0).cbrt()
}"""
new = """/// The factor on every smeared dipole's width, 1 unless `PHYS_IND_WIDTH` says
/// otherwise (`PLAY.md` E8d). The width `sigma_of` gives is the one under which a
/// lone smeared dipole has the polarisability it is given; a larger one weakens
/// the coupling of dipoles and the field of charges on them at contact, which is
/// where Pauli repulsion quenches real induction and a classical dipole does not
/// (exchange-induction). The molecule's polarisability at large separation is
/// what the sites' `alpha` are fitted to (`phys-polar --target`), so scaling the
/// width changes only the short range. Read once.
pub fn dipole_width_scale() -> f64 {
    static SCALE: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *SCALE.get_or_init(|| std::env::var("PHYS_IND_WIDTH").ok().and_then(|v| v.parse().ok()).unwrap_or(1.0))
}

/// The width of a smeared dipole of polarisability `alpha`.
pub fn sigma_of(alpha: f64) -> f64 {
    ((2.0 / PI).sqrt() * alpha / 3.0).cbrt() * dipole_width_scale()
}"""
assert old in s
open(p, 'w', encoding='utf-8').write(s.replace(old, new, 1))

p = r'E:\Projects\phys\src\bin\polar.rs'
s = open(p, encoding='utf-8').read()
old = """    let (hf, mp2) = polarisabilities(&p, frozen, step);
    println!("  19 fields in {:.0} s", t.elapsed().as_secs_f64());"""
new = """    // `--target xx,yy,zz,xy` (bohr^3): the MP2 tensor already measured, to refit the
    // atoms' polarisabilities under another dipole width (`PHYS_IND_WIDTH`) without
    // the fields being run again.
    let target: Option<Vec<f64>> = args.iter().position(|a| a == "--target").and_then(|i| args.get(i + 1)).map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect());
    let (hf, mp2) = match &target {
        Some(v) if v.len() == 4 => {
            let m = [[v[0], v[3], 0.0], [v[3], v[1], 0.0], [0.0, 0.0, v[2]]];
            (m, m)
        }
        _ => polarisabilities(&p, frozen, step),
    };
    println!("  19 fields in {:.0} s (dipole width scale {})", t.elapsed().as_secs_f64(), phys::induction::dipole_width_scale());"""
assert old in s
open(p, 'w', encoding='utf-8').write(s.replace(old, new, 1))
print('patched')
