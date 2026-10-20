play = r'E:\Projects\phys\docs\PLAY.md'
s = open(play, encoding='utf-8', newline='').read()

old = """*Detail, tables and the run log: `docs/E8D.md`. What follows is what was
     measured and what it supports; the Ewald dynamics were still running when
     this was written and their numbers are in `docs/E8D.md`.*"""
new = """*Detail, tables and the run log: `docs/E8D.md`. What follows is what was
     measured and what it supports; the Ewald dynamics of the second law (the
     no-wall one) were still running when this was written and their numbers are
     in `docs/E8D.md`.*"""
assert old in s
s = s.replace(old, new, 1)

a = s.index("U/N is 0.6 kcal/mol less bound. **Ewald dynamics**")
b = s.index("     - **Open, the owner's.**")
new_block = """U/N is 0.6 kcal/mol less bound. **Ewald dynamics**, E8c wall-fit law, 216
       molecules, 1 g/cm3, 298 K, s' 4.2, 30 ps (24 sampled): **pressure -701
       bar, -916 with the dispersion tail** (truncated -1,872 +- 73; the 5-block
       standard error printed is 20 bar but successive samples scatter by about
       1,000 bar, so take +- 200), U/N -10.43 (truncated -11.18), heat of
       vaporisation 11.0 (11.8), g(OO) first peak 2.68 at 2.85 A (2.83), first
       minimum 0.86 at 4.15 A (0.88 at 4.15), **9.07 neighbours inside it
       (8.99)**, tetrahedral q 0.541 (0.558; 0.547 with a 17.5-bohr cut-off),
       hydrogen bonds 3.45 a molecule (3.43; 3.50). **So the truncation cost
       the pressure about 1.2 kbar and the binding energy 0.75 kcal/mol a
       molecule, and cost the structure nothing: the 9-neighbour, q 0.54
       liquid is the law's, not the cut-off's** - which is what the pressure
       budget and the literature's collapse of pair-based ab initio water say.
       What is left, -0.7 to -0.9 kbar, is within the sensitivity found above (an
       8% polarisability, a 0.04 kcal/mol/bohr slope in the 3.8-4.7 A shell).
"""
s = s[:a] + new_block + s[b:]

a = s.index("     - **Open, the owner's.**")
b = s.index("   - **E8z, remove what was not chosen.**")
new_open = """     - **Open, the owner's.** (1) What puts nine neighbours inside 4.15 A
       when water has four to five: not the truncation (above) and not a
       three-body repulsion (no such term is in the Hartree-Fock triples); the
       pair law's contact and mid-range force and the induction's strength at
       contacts are what the budget says it is sensitive to, and the
       reference (MP2 in a basis 5% short in the dimer) may itself sit denser than
       water; (2) the polarisability / induction at liquid contacts (the triples
       say weaker, the pressure a little weaker, exchange-induction quenching is
       the physics a point-dipole model lacks) - a law-level addition the plan
       does not name; (3) the Born-Mayer second repulsion shape as the E8c refit
       (closes the holes; its Ewald liquid is untried); (4) a long-range sum and
       the dispersion tail (-117 to -214 bar) as part of every liquid run, now
       that they exist; (5) a second molecule, since only water has pair data.
       Not measured: Ewald dynamics for more than one law, runs longer than 30
       ps, the effect of the tin-foil boundary, any molecule but water.
"""
s = s[:a] + new_open + s[b:]
open(play, 'w', encoding='utf-8', newline='').write(s)

notes = r'E:\Projects\phys\runs\many\E8d-notes.md'
t = open(notes, encoding='utf-8').read()
t += """
### 7.6 Result: Ewald dynamics of the E8c wall-fit law (finished 05:4x AEST, 16 Oct)

216 molecules, 1.0 g/cm3, 298 K, s' = 4.2 bohr, k accuracy 1e-5 (1535 vectors; 1e-5 changes
the static pressure by 5 bar against 1e-7), r_cut 16, 30 ps with 24 ps sampled
(`runs/e8d/ewald/opt1`), against the truncated runs of the same law:

| | Ewald, 30 ps | truncated 16 bohr, 80 ps (E8c) | truncated 17.5 bohr, 40 ps |
|---|---|---|---|
| pressure, bar | **-701** (-916 with the dispersion tail of -214) | -1,872 +- 73 | -1,667 +- 67 |
| U/N, kcal/mol | -10.427 +- 0.013 (-10.519 with the tail) | -11.181 | -11.119 |
| heat of vaporisation | 11.02 (11.11 with the tail) | 11.77 | 11.71 |
| g(OO) first peak | 2.68 at 2.85 A | 2.83 at 2.85 | 2.79 at 2.85 |
| first minimum | 0.86 at 4.15 A | 0.88 at 4.15 | 0.88 at 4.15 |
| neighbours inside the minimum | **9.07** | 8.99 | 9.01 |
| tetrahedral q (one frame) | **0.541** | 0.558 | 0.547 |
| hydrogen bonds a molecule (one frame) | 3.45 | 3.43 | 3.50 |

The block standard error printed for the Ewald pressure (20 bar) is not believable: the
samples run -1,934 ... +211 with a typical successive scatter of ~1,000 bar; +- 200 is the
honest figure for 24 ps. **The truncation took 1.2 kbar off the pressure and 0.75 kcal/mol a
molecule off the binding, and took nothing off the structure.** The collapsed 9-neighbour
liquid is the law's.
"""
open(notes, 'w', encoding='utf-8').write(t)
print('updated')
