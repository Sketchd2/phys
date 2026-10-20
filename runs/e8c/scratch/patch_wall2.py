p = r"E:\Projects\phys\src\pairs.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
old = '    let out = format!("pairs-{name}{tag}{}{}{}.txt", if mono.bias.is_some() { "-bias" } else { "" }, match &mono.flex'
assert old in s
new = '''    if let Some(i) = args.iter().position(|a| a == "--wall") {
        let lo: f64 = args.get(i + 1).and_then(|v| v.parse().ok()).expect("--wall LO HI (angstrom)");
        let hi: f64 = args.get(i + 2).and_then(|v| v.parse().ok()).expect("--wall LO HI (angstrom)");
        mono.wall = Some((lo, hi));
    }
    let out = format!("pairs-{name}{tag}{}{}{}{}.txt", if mono.bias.is_some() { "-bias" } else { "" }, if mono.wall.is_some() { "-wall" } else { "" }, match &mono.flex'''
s = s.replace(old, new, 1)
open(p, "w", newline="").write(s.replace("\n", nl))
print("ok")
