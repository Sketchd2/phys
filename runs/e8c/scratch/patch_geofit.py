p = r"E:\Projects\phys\src\bin\geofit.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
# the flag
s = s.replace('    let out = flag("--out").unwrap_or_else(|| "cloud-geo.txt".into());\n',
              '    let out = flag("--out").unwrap_or_else(|| "cloud-geo.txt".into());\n    // `--hold-c0`: the cloud at the reference angle stays as it was given (a cloud already\n    // fitted and tried), and only its slopes and curvatures in the angle are fitted: the\n    // continuation from the nearest fitted geometry that the plan asks for.\n    let hold_c0 = args.iter().any(|a| a == "--hold-c0");\n', 1)
old = "        for j in 0..n {\n            let h = 1e-6 * (1.0 + z[j].abs());"
assert old in s
s = s.replace(old, "        for j in 0..n {\n            if hold_c0 && j % 3 == 0 {\n                continue;\n            }\n            let h = 1e-6 * (1.0 + z[j].abs());", 1)
s = s.replace("//!     [--dipoles props.txt] [--theta0 104.078] [--ridge 0.5] [--out cloud-geo.txt]\n",
              "//!     [--dipoles props.txt] [--theta0 104.078] [--ridge 0.5] [--hold-c0] [--out cloud-geo.txt]\n", 1)
open(p, "w", newline="").write(s.replace("\n", nl))
print("patched")
