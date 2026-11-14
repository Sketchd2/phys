#!/usr/bin/env python3
"""A law whose repulsion is one constant times the overlap of the molecule's own
charge clouds: K_ab = K q_a q_b (pi (s_a^2 + s_b^2))^(-3/2) for every pair of
cloud types, the clouds being the Gaussians the monomer's potential was fitted
with (electrons each: their charges sum to the valence count).
usage: build.py LAW_IN K [--no-core] > law_out"""
import sys, math
src = open(sys.argv[1]).read().splitlines(); K = float(sys.argv[2])
q = [float(x) for x in next(l for l in src if l.startswith('charges')).split()[1:]]
s = [float(x) for x in next(l for l in src if l.startswith('sigma')).split()[1:]]
out = []
for l in src:
    w = l.split()
    if w[0] == 'overlap': continue
    if w[0] == 'pair' and '--no-core' in sys.argv: w[3] = '1.0000000000000237e-300' if float(w[3]) != 0 else w[3]; l = ' '.join(w)
    out.append(l)
    if w[0] == 'pair' and w[1:3] == [str(len(q)-1)]*2:
        for a in range(len(q)):
            for b in range(a, len(q)):
                if s[a] > 0 and s[b] > 0:
                    out.append(f"overlap {a} {b} {K*q[a]*q[b]*(math.pi*(s[a]**2+s[b]**2))**-1.5:.12e}")
print('\n'.join(out))
