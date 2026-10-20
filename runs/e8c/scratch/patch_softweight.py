p = r"E:\Projects\phys\src\liquid.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
s = s.replace("    pub soft_wall: f64,\n", "    pub soft_wall: f64,\n    /// The soft weight's size at zero energy (0 means 0.05): the smaller it is, the less the\n    /// wall pulls on the well.\n    pub soft_weight: f64,\n", 1)
s = s.replace("                0.05 / (1.0 + (d.energy / held.soft_wall).powi(2))", "                (if held.soft_weight > 0.0 { held.soft_weight } else { 0.05 }) / (1.0 + (d.energy / held.soft_wall).powi(2))", 1)
open(p, "w", newline="").write(s.replace("\n", nl))

p = r"E:\Projects\phys\src\bin\fit.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
s = s.replace("soft_wall: soft_wall_arg() }", "soft_wall: soft_wall_arg(), soft_weight: soft_weight_arg() }")
s = s.replace("\nfn main() {", """
/// `--soft-weight W0`: the soft wall weight's size (see `Held::soft_weight`).
fn soft_weight_arg() -> f64 {
    let args: Vec<String> = std::env::args().collect();
    args.iter().position(|a| a == "--soft-weight").and_then(|i| args.get(i + 1)).and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0)
}

fn main() {""", 1)
open(p, "w", newline="").write(s.replace("\n", nl))

p = r"E:\Projects\phys\tests\liquid.rs"
s = open(p, newline="").read()
s = s.replace("soft_wall: 0.0 }", "soft_wall: 0.0, soft_weight: 0.0 }", 1)
open(p, "w", newline="").write(s)
print("ok")
