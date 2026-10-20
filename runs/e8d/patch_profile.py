p = r'E:\Projects\phys\src\bin\ewald.rs'
s = open(p, encoding='utf-8').read()
old = "    if !check {\n        return;\n    }"
new = """    // `profile`: where a force evaluation spends its time, warm-started as in a run.
    if args.iter().any(|a| a == "profile") {
        let sp = widths[0];
        let plain = EwaldLaw { inner: &law, ewald: Ewald::new(sp, cell, acc) };
        let full = EwaldLaw { inner: &polarised, ewald: Ewald::new(sp, cell, acc) };
        let time = |label: &str, f: &dyn Fn(Option<&[Vec<Vec3>]>) -> Vec<Vec<Vec3>>| {
            let mut warm = f(None);
            let t = std::time::Instant::now();
            let reps = 6;
            for _ in 0..reps {
                warm = f(Some(&warm));
            }
            println!("  {label:46} {:7.1} ms", t.elapsed().as_secs_f64() * 1000.0 / reps as f64);
        };
        println!("profile at s' = {sp}, accuracy {acc} ({} k), per force evaluation:", full.ewald.nk());
        time("pair loop only (no induction, no Ewald)", &|_| {
            liquid.forces(&law);
            Vec::new()
        });
        time("pair loop + charge reciprocal sum (Ewald, no induction)", &|_| {
            liquid.forces(&plain);
            Vec::new()
        });
        time("truncated with induction (as phys-bulk)", &|w| liquid.forces_from(&polarised, w).dipoles);
        time("Ewald with induction (the production run)", &|w| liquid.forces_from(&full, w).dipoles);
        return;
    }
    if !check {
        return;
    }"""
assert old in s
open(p, 'w', encoding='utf-8').write(s.replace(old, new, 1))
print('ok')
