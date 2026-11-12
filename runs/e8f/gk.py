#!/usr/bin/env python3
"""The repulsion against candidate functionals of the two monomer densities.
Reads overlap-*.txt (k sep E_total E_hf E_es S dT_tf dT_w dE_x). Target:
E_hf - E_es, less the law's induction where ind/ind-<species>.txt has it."""
import glob, os
import numpy as np
KCAL = 627.509474
for f in sorted(glob.glob('overlap-*.txt')):
    sp = f[len('overlap-'):-4]
    ind = {}
    if os.path.exists(f'ind/ind-{sp}.txt'):
        ind = {int(l.split()[0]): float(l.split()[1]) for l in open(f'ind/ind-{sp}.txt')}
    print(f"{sp}{' (law induction removed)' if ind else ''}:   S      target   | K=target/S | TF+X   ratio | TF+W/9+X ratio | parts: TF  W  X   [kcal/mol]")
    rows = []
    for l in open(f):
        v = l.split(); k = int(v[0]); et, eh, es, s, tf, w, x = map(float, v[2:9])
        t = eh - es - ind.get(k, 0.0)
        rows.append((s, t * KCAL, t / s, (tf + x) * KCAL, (tf + x) / t, (tf + w / 9 + x) * KCAL, (tf + w / 9 + x) / t, tf * KCAL, w * KCAL, x * KCAL))
    for r in sorted(rows):
        if r[0] < 2e-5: continue
        print(f"   {r[0]:.3e} {r[1]:+8.3f} | {r[2]:6.2f} | {r[3]:+8.3f} {r[4]:5.2f} | {r[5]:+8.3f} {r[6]:5.2f} | {r[7]:+8.3f} {r[8]:+8.3f} {r[9]:+8.3f}")
