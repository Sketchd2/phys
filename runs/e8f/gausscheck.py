#!/usr/bin/env python3
"""The Gaussian model's overlap between two molecules against the exact one.
usage: gausscheck.py NAME PAIRS OVERLAP [column of the exact overlap, default 7 = valence]"""
import sys, numpy as np
name, pf, of = sys.argv[1:4]; col = int(sys.argv[4]) if len(sys.argv) > 4 else 7
Z = {'water': [8, 1, 1], 'methane': [6, 1, 1, 1, 1], 'ammonia': [7, 1, 1, 1]}[name]
M = np.load(f'gauss-{name}.npy')
geo = {}
for l in open(pf):
    if l.startswith('#') or '|' not in l: continue
    h, t = l.split('|'); t = t.split(); geo[int(h.split()[0])] = [np.array(list(map(float, t[i:i+3]))) for i in range(0, len(t), 4)]
def sites(at): return [(p, s, a) for p, z in zip(at, Z) for zz, s, a in M if zz == z]
ex, mo = [], []
for l in open(of):
    v = l.split()
    if len(v) != 9: continue
    g = geo[int(v[0])]; n = len(Z); A, B = sites(g[:n]), sites(g[n:]); S = 0.0
    for pa, sa, qa in A:
        for pb, sb, qb in B:
            s2 = sa * sa + sb * sb; S += qa * qb * (np.pi * s2) ** -1.5 * np.exp(-np.sum((pa - pb) ** 2) / s2)
    ex.append(float(v[col])); mo.append(S)
ex, mo = np.array(ex), np.array(mo); m = ex > 3e-5; r = mo[m] / ex[m]
print(f"{name} {of}: {m.sum()} pairs, model/exact overlap: mean {r.mean():.3f}, spread {r.std()/r.mean()*100:.1f}%, range {r.min():.2f}-{r.max():.2f}")
