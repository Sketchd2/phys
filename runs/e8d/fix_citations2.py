p = r'E:\Projects\phys\runs\many\E8d-notes.md'
s = open(p, encoding='utf-8').read()
a = s.index("4. **Quantum nuclei**")
b = s.index("## 5. Hypotheses")
new = """4. **Quantum nuclei** shift the density by ~0.3 % for MB-pol (1.004 → 1.001, Reddy et al.); they
   broaden g(r) and weaken H-bonds (Wang, Ceriotti, Markland 2014). By the liquid's compressibility
   0.3 % is ~70 bar: my inference from the density shift, not a statement of those papers.
5. **MP2 for the liquid**: MP2/TZ-class baselines over-structure water (first peak ~11 %
   high); results depend strongly on basis convergence (2026 reassessment: "MP2 with a
   triple-ζ basis demonstrates poor performance indicative of insufficient convergence"; an
   MP2/RPA isobaric-isothermal Monte Carlo study by Del Ben, VandeVondele and Hutter reports RPA
   0.994 g/cm³ and an over-structured MP2). The engine's basis is derived and ~5 % short in the dimer.
6. **Long-range dispersion/cutoff** artefacts bias pressure and density (the study of Ewald-summed
   dispersion, arXiv:0704.0474, finds them secondary to Coulomb in water). The dispersion tail from the
   law's own C6 is −117 to −214 bar (section 2d corrects the first, sharp-cut-off estimate).

"""
s = s[:a] + new + s[b:]
open(p, 'w', encoding='utf-8').write(s)

q = r'E:\Projects\phys\runs\e8d\make_e8d_doc.py'
t = open(q, encoding='utf-8').read()
t = t.replace("hydrogen bond and broaden g(r) but broaden g(r): Wang,", "hydrogen bond and broaden g(r): Wang,")
open(q, 'w', encoding='utf-8').write(t)
print('ok')
