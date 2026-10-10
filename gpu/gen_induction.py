"""Fill the tables of src/induction.wgsl.tpl and write src/induction.wgsl.

    python gen_induction.py

The reciprocals 1/(2n+1) of erf's series and the coefficients 1/k! of the
exponential's polynomial are computed here, not typed. Run it after editing the
template; the generated file is committed so a build needs no Python.
"""
import math
import os

here = os.path.dirname(os.path.abspath(__file__))


def lit(x):
    s = "%.17g" % x
    if "." not in s and "e" not in s:
        s += ".0"
    return s + "lf"


inv_odd = ", ".join(lit(1.0 / (2 * n + 1)) for n in range(130))
c = [1.0 / math.factorial(k) for k in range(15)]
horner = "var poly = %s;\n" % lit(c[14]) + "".join("    poly = poly * r + %s;\n" % lit(c[k]) for k in range(13, -1, -1))

text = open(os.path.join(here, "src", "induction.wgsl.tpl"), encoding="utf8").read()
text = text.replace("@INV_ODD@", inv_odd).replace("@HORNER@", horner.rstrip("\n"))
open(os.path.join(here, "src", "induction.wgsl"), "w", encoding="utf8", newline="\n").write(text)
print("wrote src/induction.wgsl,", len(text.splitlines()), "lines")
