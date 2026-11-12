#!/usr/bin/env python3
"""The engine's pair energy, part by part, against the monomer-record model.
usage: parts.py PAIRS OVERLAP "a-b:C6,..."   (C6 by atom-type pair, hartree bohr^6)
E_corr = E_total - E_hf against -sum C6/r^6 (undamped), and
E_hf - E_es against the density overlap."""
import sys, math
import numpy as np
KCAL = 627.509474
c6 = {}
for a in sys.argv[3].split(','):
    k, v = a.split(':'); i, j = map(int, k.split('-')); c6[(i, j)] = c6[(j, i)] = float(v)
geo = {}
for l in open(sys.argv[1]):
    if l.startswith('#') or '|' not in l: continue
    h, t = l.split('|'); t = t.split(); at = [(np.array(list(map(float, t[i:i+3]))), int(t[i+3])) for i in range(0, len(t), 4)]
    geo[int(h.split()[0])] = at
print("   k  near(A)   S        E_es    E_hf-E_es  ratio/S   E_corr   -C6/r^6  ratio   (kcal/mol)")
R = []
for l in open(sys.argv[2]):
    k, sep, et, eh, es, s = l.split(); k = int(k); et, eh, es, s = map(float, (et, eh, es, s))
    at = geo[k]; n = len(at) // 2
    d6 = 0.0; near = 1e9
    for pa, ta in at[:n]:
        for pb, tb in at[n:]:
            r = np.linalg.norm(pa - pb); near = min(near, r); d6 -= c6[(ta, tb)] / r ** 6
    R.append((near * 0.529177, s, es * KCAL, (eh - es) * KCAL, (eh - es) / s, (et - eh) * KCAL, d6 * KCAL))
for r in sorted(R):
    print(f"      {r[0]:6.2f} {r[1]:.3e} {r[2]:+8.3f} {r[3]:+8.3f} {r[4]:8.3f} {r[5]:+8.3f} {r[6]:+8.3f} {r[5]/r[6]:6.2f}")
