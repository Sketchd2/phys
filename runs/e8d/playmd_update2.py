p = r'E:\Projects\phys\docs\PLAY.md'
s = open(p, encoding='utf-8', newline='').read()
marker = "     - **Open, the owner's.** (1) What puts nine neighbours inside 4.15 A"
assert marker in s and 'hole-free law' not in s
bullet = """     - **Looking for a better law (the owner: "the best solution, even a novel
       one"), and what it showed.** Tried against the engine's own reference
       before any liquid was run, mostly negative: widening the induced dipoles
       (with the site polarisabilities refitted to the same MP2 tensor, so the
       gas-phase value is kept) makes the 300 Hartree-Fock triples fit worse,
       0.115 -> 0.134 -> 0.161 -> 0.281 kcal/mol at width scale 1.0, 1.2, 1.4,
       1.6; holding the Born-Mayer exponent at the density tail's 2 sqrt(2 IP),
       about 2 per bohr (`phys-fit --bm-exponent`), costs pair accuracy and leaves
       the shell error; zero-charge lone-pair repulsion sites are driven to zero
       strength by the fit; a population-weighted fit leaves the 3.1-3.4 A energy
       bias at -0.2 kcal/mol. The 3.1-3.8 A error all laws share is on bent H..O
       contacts (95-125 degrees) and is a repulsion the isotropic form lacks
       (`phys-slopeparts`). **What worked**: the no-wall law has the best Ewald
       liquid and holes inside 1 A (-577 kcal/mol at 0.6-0.8 A); a steep Born-Mayer
       core on O-H alone (exponent 7, amplitude fitted to what the law leaves on the
       wall pairs, the rest held; `runs/e8d/core/make_core.py`) removes every hole
       (+33.5, +20.6, +10.2, +2.6 kcal/mol at 0.6-1.4 A; deepest dimer -4.63 at
       2.92 A) and leaves the liquid range unchanged. Its Ewald liquid: -88 +- 97
       bar (-204 with the tail), U/N -10.21, heat of vaporisation 10.8, g(OO) first
       peak 2.57 at 2.85 A (measured 2.57-2.75 at 2.80), 7.56 neighbours inside
       3.95 A, q 0.57, 3.61 hydrogen bonds: the best combination so far, still with
       no first-shell gap. The Ewald force evaluation was also made 1.6 times
       cheaper on a quiet machine (117 -> 72 ms: analytic reciprocal virial,
       coincident charges merged), with identical results.
"""
s = s.replace(marker, bullet + marker, 1)
open(p, 'w', encoding='utf-8', newline='').write(s)
print('ok')
