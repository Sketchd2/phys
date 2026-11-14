#!/usr/bin/env python3
"""The `charges`, `sigma` and `site` lines of a law from a Gaussian model of a molecule's own
density (gaussfit.py --col 5 --moments: the density without the 1s cores; each nucleus a point of
charge Z less its core).  Every (kind of site, width) is a site type; a site is placed in the frame
of the first three atoms as the simulator's `site` lines are.  Gaussians under 1e-4 electrons are
left out.   usage: build3.py NAME > cloud.txt   (run from runs/e8f after gaussfit)"""
import sys, numpy as np
name = sys.argv[1]
atoms = [(int(l.split()[1]), np.array(list(map(float, l.split()[2:5])))) for l in open(f'density-{name}.txt') if l.startswith('atom')]
model = np.load(f'gauss-{name}.npy')
zs = []
for z, _ in atoms:
    if z not in zs: zs.append(z)
assert atoms[0][0] == zs[0]
o, p1, p2 = atoms[0][1], atoms[1][1], atoms[2][1]
u, v = (p1 - o) / np.linalg.norm(p1 - o), (p2 - o) / np.linalg.norm(p2 - o); n = np.cross(u, v); n /= np.linalg.norm(n)
F = np.array([u, v, n]).T
types = {}; sites = []
for m in model:
    if abs(m[4]) < 1e-4: continue
    key = (round(m[3], 6), round(m[4], 10))       # sites alike by symmetry share width and population
    t = types.setdefault(key, len(zs) + len(types))
    sites.append((t, np.linalg.solve(F, m[:3] - o)))
q = [float(z - 2 if z > 2 else z) for z in zs] + [-k[1] for k in types]
s = [0.0] * len(zs) + [k[0] for k in types]
print('charges ' + ' '.join(f'{x:.12e}' for x in q))
print('sigma ' + ' '.join(f'{x:.12e}' for x in s))
for t, f in sites: print(f'site {t} {f[0]:.12e} {f[1]:.12e} {f[2]:.12e}')
