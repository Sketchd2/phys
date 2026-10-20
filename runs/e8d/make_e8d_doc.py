notes = open(r'E:\Projects\phys\runs\many\E8d-notes.md', encoding='utf-8').read()
head = """# E8d: where the liquid's pressure comes from

The detailed record behind the E8d entry in `PLAY.md` section 7: the many-body
probe, the pressure cross-check, the literature, the designs tried and the runs,
as they were made (2026-10-14 to 16). Sections 1-6 are the probe, the pressure
budget, the refits and the cut-off study; section 7 is the Ewald sum. Numbers are
as measured at the time they were written; where a run was still going the text
says so, and the table at the end of section 7 is the place the finished numbers
go. Everything is reproducible from the tools named here (`phys-cluster-gpu
--manybody`, `phys-manybody`, `phys-virial`, `phys-ewald`, `phys-bulk` with
`PHYS_EWALD`, `phys-fit --bm`) and the scripts under `runs/many`, `runs/virial`
and `runs/e8d`; the data are not committed (they are 60 MB and regenerable).

"""
sources = """

## Sources

Pair-only ab initio water collapses (coordination 11.9 for the CC-pol pair
potential alone) and three-body terms restore tetrahedral coordination:
Bukowski, Szalewicz, Groenenboom, van der Avoird, Science 315, 1249 (2007),
https://pubmed.ncbi.nlm.nih.gov/17332406/ . The neglect of short-range
three-body exchange-repulsion and charge-transfer terms leads to an incorrect
liquid structure: Babin, Medders, Paesani (HBB2-pol), arXiv:1210.7022 .
MB-pol: three-body about 20% of the interaction energy; classical against
quantum density 1.004 against 1.001 g/cm3 at 298 K: Reddy et al., arXiv:1609.02884 .
Beyond-two-body errors of PBE come from close approach of monomers that are not
hydrogen-bonded to each other, BLYP's from hydrogen-bond cooperativity: Gillan,
Alfe, Bartok, Csanyi, arXiv:1303.0751 . Nuclear quantum effects weaken the
hydrogen bond and broaden g(r): Wang,
Ceriotti, Markland, arXiv:1405.7366 (the pressure size quoted in the notes is my inference from the density shift). MP2 with a triple-zeta basis over-structures
liquid water and depends on basis convergence: arXiv:2603.20451; Del Ben,
VandeVondele, Hutter (MP2 and RPA isobaric-isothermal Monte Carlo; RPA 0.994 g/cm3).
Ewald summation of long-range dispersion, and cut-off artefacts in pressure and density: arXiv:0704.0474.
"""
open(r'E:\Projects\phys\docs\E8D.md', 'w', encoding='utf-8').write(head + notes + sources)
print('written', len(head + notes + sources))
