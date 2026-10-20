p = r'E:\Projects\phys\docs\PLAY.md'
s = open(p, encoding='utf-8', newline='').read()
marker = "     - **Open, the owner's.** (1) What puts nine neighbours inside 4.15 A"
assert marker in s and 'ball-and-stick' not in s
bullet = """     - **The owner's simpler method: ball-and-stick molecules, atom-atom forces
       by brute force** (`phys-sticks`, `runs/e8d/sticks/`; docs/E8D.md section 8).
       Springs from the Hessian of the engine's monomer surface (it reproduces the
       engine's harmonic frequencies, 1626, 3847, 3973 cm^-1), every intermolecular
       atom pair summed, charges by Ewald, Langevin dynamics, about 16 ms a step
       for 64 molecules (roughly 7 times cheaper a picosecond than the rigid liquid
       with induction, though the flexible hydrogens need a step half as long).
       The atom-atom law is point charges, an exponential repulsion and a damped
       C6 per atom-pair type, fitted to the engine's 410 MP2 pair energies, with
       no clouds, no extra site and no induction. Where it falls short, so far: the
       pair error is 0.72-0.74 kcal/mol held out against 0.24-0.30 for the cloud +
       induction law; with the dispersion free the fit takes the molecule-to-
       molecule C6 to about 160 hartree bohr^6 (about 50 physical) and one atom-pair
       C6 negative, so dispersion is standing in for what pair additivity and point
       charges cannot do; a 64-molecule 3 ps test collapses to 12.8 neighbours
       inside the first minimum with P about -2.8 kbar. The 216-molecule 20 ps
       runs for three variants (gas-phase charge, fitted charge, dispersion held
       at the engine's derived C6) were still running when this was written.
"""
s = s.replace(marker, bullet + marker, 1)
open(p, 'w', encoding='utf-8', newline='').write(s)
print('ok')
