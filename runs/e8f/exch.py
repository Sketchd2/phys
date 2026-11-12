#!/usr/bin/env python3
"""Exchange repulsion (Heitler-London less exact electrostatics) against the
density overlap, species by species. Reads overlap-*.txt, columns
k sep E_total E_hf E_es S dT_tf dT_w dE_x E_hl."""
import glob, math
import numpy as np
KCAL = 627.509474
allr = {}
for f in sorted(glob.glob('overlap-*.txt')):
    sp = f[len('overlap-'):-4]; rows = []
    for l in open(f):
        v = l.split()
        if len(v) < 10: continue
        et, eh, es, s = map(float, v[2:6]); hl = float(v[9])
        rows.append((s, (hl - es) * KCAL, (hl - es) / s, (eh - hl) * KCAL, es * KCAL, (et - eh) * KCAL))
    allr[sp] = sorted(rows)
    print(f"{sp}:     S       exchange   K = exch/S   HF beyond   E_es    E_corr   [kcal/mol]")
    for r in sorted(rows):
        if r[0] > 2e-5: print(f"   {r[0]:.3e} {r[1]:+8.3f}  {r[2]:7.2f}   {r[3]:+8.3f} {r[4]:+8.3f} {r[5]:+8.3f}")
print("\nspecies            n(S>1e-4)  K: mean  spread   log-log slope")
for sp, rows in allr.items():
    a = np.array([r for r in rows if r[0] > 1e-4])
    if len(a) < 2: continue
    g = np.polyfit(np.log(a[:, 0]), np.log(a[:, 1]), 1)[0] if len(a) > 2 else float('nan')
    print(f"{sp:18s} {len(a):3d}      {a[:,2].mean():6.2f}  {a[:,2].std()/a[:,2].mean()*100:4.1f}%   {g:.3f}")
