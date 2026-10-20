p = r'E:\Projects\phys\runs\many\E8d-notes.md'
s = open(p, encoding='utf-8').read()
fixes = [
    ("Babin, Leforestier, Paesani, HBB2-pol,\n   2013", "Babin, Medders, Paesani, HBB2-pol,\n   2012"),
    ("Babin, Leforestier, Paesani", "Babin, Medders, Paesani"),
    ("6. **Long-range dispersion/cutoff** artefacts bias pressure and density (Isele-Holder,\n   Mitchell, Ismail 2012; in water secondary to Coulomb). We compute it: -350 ... -525 bar.",
     "6. **Long-range dispersion/cutoff** artefacts bias pressure and density (the study of Ewald-summed\n   dispersion, arXiv:0704.0474, finds them secondary to Coulomb in water). We compute the\n   dispersion tail from the law's own C6: -117 to -214 bar (section 2d corrects the first estimate)."),
    ("Del Ben et al.: RPA gives 0.994 g/cm3, MP2 over-structured", "an MP2/RPA isobaric-isothermal Monte Carlo study (Del Ben, VandeVondele, Hutter) reports RPA 0.994 g/cm3 and an over-structured MP2"),
    ("they broaden\n   g(r) and weaken H-bonds (Wang, Ceriotti, Markland 2014) but do not move the pressure by kilobars.",
     "they broaden\n   g(r) and weaken H-bonds (Wang, Ceriotti, Markland 2014). The density shift above is 0.3 %, which by\n   the liquid's compressibility is ~70 bar: my inference, not a statement of those papers."),
]
for old, new in fixes:
    if old in s:
        s = s.replace(old, new)
        print('fixed:', old[:50].replace('\n', ' '))
    else:
        print('NOT FOUND:', old[:60].replace('\n', ' '))
open(p, 'w', encoding='utf-8').write(s)

q = r'E:\Projects\phys\runs\e8d\make_e8d_doc.py'
t = open(q, encoding='utf-8').read()
t = t.replace("Babin, Leforestier, Paesani (HBB2-pol), arXiv:1210.7022", "Babin, Medders, Paesani (HBB2-pol), arXiv:1210.7022")
t = t.replace("VandeVondele, Hutter, J. Chem. Phys. 143, 054506 (2015) (RPA 0.994 g/cm3).\nLong-range dispersion and cut-off artefacts in pressure and density: Isele-Holder,\nMitchell, Ismail, arXiv:0704.0474 (Ewald summation of dispersion).",
              "VandeVondele, Hutter (MP2 and RPA isobaric-isothermal Monte Carlo; RPA 0.994 g/cm3).\nEwald summation of long-range dispersion, and cut-off artefacts in pressure and density: arXiv:0704.0474.")
t = t.replace("Ceriotti, Markland, arXiv:1405.7366 . MP2 with a triple-zeta", "Ceriotti, Markland, arXiv:1405.7366 (the pressure size quoted in the notes is my inference from the density shift). MP2 with a triple-zeta")
t = t.replace("do not move the pressure by kilobars", "broaden g(r)")
open(q, 'w', encoding='utf-8').write(t)
print('doc script updated')
