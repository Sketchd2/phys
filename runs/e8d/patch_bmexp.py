import re
# liquid.rs: a Held field and its use
p = r'E:\Projects\phys\src\liquid.rs'
s = open(p, encoding='utf-8').read()
old = "    pub no_born_mayer: bool,\n"
assert old in s
s = s.replace(old, old + """    /// If positive, the Born-Mayer exponent of every pair is held at this value
    /// (bohr^-1) and only the amplitudes are fitted: the exponent of exchange
    /// repulsion is set by the density's tail, `2 sqrt(2 IP)` (`PLAY.md` E8d).
    pub bm_exponent: f64,
""", 1)
old = "    let base = pack(start);\n    let is_held"
assert old in s
s = s.replace(old, """    let mut base = pack(start);
    if held.bm_exponent > 0.0 {
        for k in 0..pairs.len() {
            base[types - 1 + 4 * k + 1] = held.bm_exponent;
        }
    }
    let is_held""", 1)
old = "            if (held.repulsion && j < 2) || (held.no_born_mayer && j < 2)"
assert old in s
s = s.replace(old, "            if (j == 1 && held.bm_exponent > 0.0) || (held.repulsion && j < 2) || (held.no_born_mayer && j < 2)", 1)
open(p, 'w', encoding='utf-8').write(s)

# fit.rs: read the flag and pass it
p = r'E:\Projects\phys\src\bin\fit.rs'
s = open(p, encoding='utf-8').read()
n = s.count("no_born_mayer: use_overlap && !with_bm,")
print('Held literals in fit.rs:', n)
s = s.replace("no_born_mayer: use_overlap && !with_bm,", "no_born_mayer: use_overlap && !with_bm, bm_exponent: bm_exponent_arg(),")
old = "fn soft_wall_arg() -> f64 {"
assert old in s
s = s.replace(old, """/// `--bm-exponent B` (bohr^-1): hold every Born-Mayer exponent at B and fit only the amplitudes.
fn bm_exponent_arg() -> f64 {
    let args: Vec<String> = std::env::args().collect();
    args.iter().position(|a| a == "--bm-exponent").and_then(|i| args.get(i + 1)).and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0)
}

""" + old, 1)
open(p, 'w', encoding='utf-8').write(s)
print('ok')
