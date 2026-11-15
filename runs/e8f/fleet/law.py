#!/usr/bin/env python3
"""A species' law from its own record and its pairs' energies (PLAY.md E8f, HANDOFF.md).
usage, from runs/e8f/species/NAME after the stages grow, mono, pairs and full:

    law.py NAME [cp|uncorrected]

The law is the molecule's charge clouds (cloud-NAME.txt, phys-cloud), its atom-pair C6 (c6.log,
phys-c6), its polarisability (polar-NAME.txt; where the split among atoms did not finish, the
molecule's MP2 mean shared equally among its heavy atoms), and a repulsion
K q_a q_b (pi (s_a^2 + s_b^2))^(-3/2) exp(-r^2 / (s_a^2 + s_b^2)) between every two clouds.
Two numbers are fitted to the pairs' total energies (counterpoise corrected unless
`uncorrected` is asked for): K, and the exponent of the dispersion's damping, scanned from
1.0 to 2.4. A pair counts 1 if it is bound and exp(-E / e) if not, e the depth of the
species' deepest pair, so the fit is of what a liquid sees.

Writes law.txt, and law.ok with the two numbers, the errors by depth, and the second virial
coefficient at two temperatures: a number of hundreds of digits there means the law has a hole
at contact, and law.ok then begins HOLE and the law is not to be run as a liquid."""
import sys, os, math, subprocess
import numpy as np
H = 627.509474
name = sys.argv[1]; which = sys.argv[2] if len(sys.argv) > 2 else 'cp'
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')); B = ROOT + '/target/release'

cloud = [l for l in open(f'cloud-{name}.txt').read().splitlines() if l.strip()]
q = [float(x) for x in next(l for l in cloud if l.startswith('charges')).split()[1:]]
s = [float(x) for x in next(l for l in cloud if l.startswith('sigma')).split()[1:]]
atoms = sum(1 for x in s if x == 0)
c6 = {}
for l in open('c6.log'):
    if l.startswith('--c6'):
        for x in l.split()[1].split(','):
            k, v = x.split('='); c6[tuple(sorted(map(int, k.split('-'))))] = float(v)
classes = eval(next(l for l in open('c6.log') if 'atom classes' in l).split('classes')[1].strip())
z = [int(x) for x in open(f'grow-{name}.atoms').read().split()]
if os.path.exists(f'polar-{name}.txt') and any(l.startswith('alpha') for l in open(f'polar-{name}.txt')):
    alpha = [l.strip() for l in open(f'polar-{name}.txt') if l.startswith('alpha')]; alpha_note = 'split among atoms by phys-polar'
else:
    mean = float(next(l for l in open('polar.log') if 'MP2: mean' in l).split('mean')[1].split()[0])
    heavy = [i for i, zz in enumerate(z) if zz > 1]; share = mean / len(heavy)
    alpha = [f'alpha {c} {share}' for c in sorted({classes[i] for i in heavy})]; alpha_note = f'the MP2 mean {mean:.4f} shared among {len(heavy)} heavy atoms'

def law_text(K, damp):
    out = [next(l for l in cloud if l.startswith('charges'))]; n = len(q)
    for a in range(n):
        for b in range(a, n):
            has = (a, b) in c6
            out.append(f"pair {a} {b} {0.0 if has else 1e-300:.12e} {5.0 if has else 2.0} {c6.get((a, b), 1e-300):.12e} 1.0000000000000237e-300")
    for a in range(n):
        for b in range(a, n):
            if s[a] > 0 and s[b] > 0 and K != 0:
                out.append(f"overlap {a} {b} {K * q[a] * q[b] * (math.pi * (s[a] ** 2 + s[b] ** 2)) ** -1.5:.12e}")
    out += [l for l in cloud if l.split()[0] in ('sigma', 'site', 'bisector')] + alpha
    out += [f"damp {a} {b} {damp}" for a in range(atoms) for b in range(a, atoms)]
    return '\n'.join(out) + '\n'

# the reference: the pairs given the pair's own solve, as a pairs file phys-lawenergy reads
col = 12 if which == 'cp' else 10
ref = {}
for l in open(f'full-{name}.txt'):
    t = l.split()
    if len(t) > col: ref[int(t[0])] = float(t[col])
with open('pairs-water-ref.txt', 'w') as f:     # phys-lawenergy's file name, whatever the molecule
    for l in open(f'pairs-{name}-sel.txt'):
        if l.startswith('#') or '|' not in l: continue
        h, tail = l.split('|'); k = int(h.split()[0])
        if k in ref: f.write(f"{k} {h.split()[1]} {ref[k]:.10e} 0 |{tail.rstrip()}\n")
def energies(K, damp):
    open('law-tmp.txt', 'w').write(law_text(K, damp))
    r = subprocess.run([B + '/phys-lawenergy', 'law-tmp.txt', '-ref'], capture_output=True, text=True).stdout
    return np.array([[float(x) for x in l.split()[:4]] for l in r.splitlines() if l.strip()])
best = None
for damp in (1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.4):
    e0, e1 = energies(0, damp), energies(1, damp)
    r, b, sl = e0[:, 2] * H, e0[:, 3] * H, (e1[:, 3] - e0[:, 3]) * H
    eps = max(-r.min(), 0.05); w = np.where(r < 0, 1.0, np.exp(-np.maximum(r, 0) / eps))
    K = (w * (r - b) * sl).sum() / (w * sl * sl).sum(); d = b + K * sl - r
    rms = math.sqrt((w * d * d).sum() / w.sum())
    if best is None or rms < best[0]: best = (rms, damp, K, r, d, w)
rms, damp, K, r, d, w = best
open('law.txt', 'w').write(law_text(K, damp)); os.remove('law-tmp.txt')
lines = [f"K {K:.4f} damping {damp} weighted rms {rms:.3f} kcal/mol over {len(r)} pairs ({which}); deepest pair {r.min():.3f}; polarisability: {alpha_note}"]
bound = r < 0
if bound.any(): lines.append(f"bound pairs {bound.sum()}: mean error {d[bound].mean():+.3f}, rms {d[bound].std():.3f}; deepest: reference {r.min():.3f}, law {(r + d)[r.argmin()]:.3f}")
if (~bound).any(): lines.append(f"repulsive pairs {(~bound).sum()}: reference {', '.join(f'{x:.1f}' for x in r[~bound])}; law {', '.join(f'{x:.1f}' for x in (r + d)[~bound])}")
# the hole check: the second virial coefficient at two temperatures
b2 = subprocess.run([B + '/phys-b2', name, 'law.txt', '300', '600', '--orientations', '1500'], capture_output=True, text=True).stdout.split('\n')
vals = [l.split() for l in b2 if l.strip()]
hole = (not vals) or any((not math.isfinite(float(v[2]))) or abs(float(v[2])) > 1e7 for v in vals)
lines.append('second virial coefficient (cm3/mol): ' + '; '.join(f"{v[0]} K: {float(v[2]):.4g}" for v in vals))
open('law.ok', 'w').write(('HOLE ' if hole else '') + '\n'.join(lines) + '\n')
print(open('law.ok').read())
