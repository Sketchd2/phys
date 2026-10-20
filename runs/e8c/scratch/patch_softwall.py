import re
# liquid.rs: the field and the weights
p = r"E:\Projects\phys\src\liquid.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
old = "    pub no_born_mayer: bool,\n"
assert old in s
s = s.replace(old, old + """    /// A soft weight for the repulsive wall (hartree; 0 is off): a pair more than 25
    /// kcal/mol above the lowest weighs `0.05 / (1 + (E / soft_wall)^2)` and not the
    /// Boltzmann weight's floor, so that a wall pair counts for its *relative* error
    /// (for E well above `soft_wall` the weighted residual is about `0.22 soft_wall`
    /// times the fractional miss) and a pair of 2000 kcal/mol does not outweigh the
    /// whole well. The floor alone weighs every wall pair the same, so the few
    /// highest decide.
    pub soft_wall: f64,
""", 1)
old_w = "    let weights: Vec<f64> = data.iter().map(|d| (-(d.energy - e_min) / kt).exp().max(floor)).collect();\n"
assert old_w in s
new_w = """    let weights: Vec<f64> = data
        .iter()
        .map(|d| {
            let above = d.energy - e_min;
            if held.soft_wall > 0.0 && above > 25.0 / 627.509474 {
                0.05 / (1.0 + (d.energy / held.soft_wall).powi(2))
            } else {
                (-above / kt).exp().max(floor)
            }
        })
        .collect();
"""
s = s.replace(old_w, new_w, 1)
open(p, "w", newline="").write(s.replace("\n", nl))

# fit.rs: the flag and the two struct literals
p = r"E:\Projects\phys\src\bin\fit.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
s = s.replace("no_born_mayer: use_overlap && !with_bm, slopes: &[] }", "no_born_mayer: use_overlap && !with_bm, slopes: &[], soft_wall }", 1)
s = s.replace("no_born_mayer: use_overlap && !with_bm, slopes: &slopes }", "no_born_mayer: use_overlap && !with_bm, slopes: &slopes, soft_wall }", 1)
flag = '''    // `--soft-wall E0` (kcal/mol): the repulsive wall weighed by its relative error (see
    // `Held::soft_wall`), for pair sets that reach the wall (`phys-pairs-gpu --wall`).
    let soft_wall: f64 = args.iter().position(|a| a == "--soft-wall").and_then(|i| args.get(i + 1)).and_then(|v| v.parse::<f64>().ok()).map(|k| k / KCAL).unwrap_or(0.0);
'''
marker = "    let use_overlap = args.iter().any(|a| a == \"--overlap\");\n"
assert marker in s
s = s.replace(marker, marker + flag, 1)
open(p, "w", newline="").write(s.replace("\n", nl))

# the test's literal
p = r"E:\Projects\phys\tests\liquid.rs"
s = open(p, newline="").read()
s = s.replace("no_born_mayer: false, slopes: &[] }", "no_born_mayer: false, slopes: &[], soft_wall: 0.0 }", 1)
open(p, "w", newline="").write(s)
print("patched")
