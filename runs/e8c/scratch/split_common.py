import re
p = r"E:\Projects\phys\gpu\src\induction.wgsl.tpl"
s = open(p, encoding="utf8").read()
a = s.index("const SQRT_PI")
b = s.index("// (ta, tbp, da, dbp) of the reference's `dd_terms`.")
common = s[a:b]
s = s[:a] + "@COMMON@\n\n" + s[b:]
open(p, "w", encoding="utf8", newline="\n").write(s)
hdr = """// Shared by the GPU kernels that need the reference's special functions in
// double precision: WGSL has no f64 exp, so `gexp` is built here; `gerf` is the
// reference's own series. Included by gen_induction.py in each template at
// `@COMMON@`; edit it here.

"""
open(r"E:\Projects\phys\gpu\src\common.wgsl.tpl", "w", encoding="utf8", newline="\n").write(hdr + common)
g = r"E:\Projects\phys\gpu\gen_induction.py"
t = open(g, encoding="utf8").read()
t = t.replace('''text = open(os.path.join(here, "src", "induction.wgsl.tpl"), encoding="utf8").read()
text = text.replace("@INV_ODD@", inv_odd).replace("@HORNER@", horner.rstrip("\\n"))
open(os.path.join(here, "src", "induction.wgsl"), "w", encoding="utf8", newline="\\n").write(text)
print("wrote src/induction.wgsl,", len(text.splitlines()), "lines")''', '''common = open(os.path.join(here, "src", "common.wgsl.tpl"), encoding="utf8").read()
for name in ("induction", "sitepair"):
    path = os.path.join(here, "src", name + ".wgsl.tpl")
    if not os.path.exists(path):
        continue
    text = open(path, encoding="utf8").read().replace("@COMMON@", common)
    text = text.replace("@INV_ODD@", inv_odd).replace("@HORNER@", horner.rstrip("\\n"))
    open(os.path.join(here, "src", name + ".wgsl"), "w", encoding="utf8", newline="\\n").write(text)
    print("wrote src/%s.wgsl, %d lines" % (name, len(text.splitlines())))''')
t = t.replace("Fill the tables of src/induction.wgsl.tpl and write src/induction.wgsl.", "Fill the tables of src/*.wgsl.tpl (with src/common.wgsl.tpl included at @COMMON@) and write src/*.wgsl.")
open(g, "w", encoding="utf8", newline="\n").write(t)
print("split")
