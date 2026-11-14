#!/usr/bin/env python3
"""A law from a molecule's own record and one constant: the charge clouds
(phys-cloud / phys-cloudfit: `charges`, `sigma`, `site` lines), the atom-pair
C6 (phys-c6), the atomic polarisabilities (phys-polar), and a repulsion
K q_a q_b (pi (s_a^2 + s_b^2))^(-3/2) exp(-r^2 / (s_a^2 + s_b^2)) between every
two clouds. The dispersion's damping exponent is 2.0 as in the stored water law.
usage: build2.py CLOUD.txt "0-0=..,0-1=..,1-1=.." POLAR.txt K [core "a-b=A:B,..."] > law"""
import sys, math
cloud = [l for l in open(sys.argv[1]).read().splitlines() if l.strip()]
c6 = {tuple(sorted(map(int, k.split('-')))): float(v) for k, v in (x.split('=') for x in sys.argv[2].replace('--c6', '').split(','))}
alpha = [l for l in open(sys.argv[3]).read().splitlines() if l.startswith('alpha')]
K = float(sys.argv[4])
core = {}
if len(sys.argv) > 6:
    for x in sys.argv[6].split(','):
        k, v = x.split('='); core[tuple(sorted(map(int, k.split('-'))))] = tuple(map(float, v.split(':')))
q = [float(x) for x in next(l for l in cloud if l.startswith('charges')).split()[1:]]
s = [float(x) for x in next(l for l in cloud if l.startswith('sigma')).split()[1:]]
n = len(q); atoms = sum(1 for x in s if x == 0)
print(next(l for l in cloud if l.startswith('charges')))
for a in range(n):
    for b in range(a, n):
        A, B = core.get((a, b), (0.0 if (a, b) in c6 else 1e-300, 5.0 if (a, b) in c6 else 2.0))
        print(f"pair {a} {b} {A:.12e} {B} {c6.get((a, b), 1e-300):.12e} 1.0000000000000237e-300")
for a in range(n):
    for b in range(a, n):
        if s[a] > 0 and s[b] > 0 and K != 0:
            print(f"overlap {a} {b} {K*q[a]*q[b]*(math.pi*(s[a]**2+s[b]**2))**-1.5:.12e}")
for l in cloud:
    if l.split()[0] in ('sigma', 'site', 'bisector'): print(l)
for l in alpha: print(l)
for a in range(atoms):
    for b in range(a, atoms): print(f"damp {a} {b} 2.0")
