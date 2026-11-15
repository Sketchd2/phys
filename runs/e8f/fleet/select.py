#!/usr/bin/env python3
"""Choose the pairs a liquid's neighbours look like from a file of drawn geometries.
usage: select.py PAIRS-GEOM NATOMS COUNT > selected
A pair is kept if its closest atom-atom distance is between 1.6 and 4.5 A (inside 1.6 two atoms
have collided; beyond 4.5 nothing overlaps), and COUNT are taken spread evenly over that
distance, nearest first, so the wall, the contact shell and the first approach are all there."""
import sys, math
n = int(sys.argv[2]); want = int(sys.argv[3]); rows = []
for l in open(sys.argv[1]):
    if l.startswith('#') or '|' not in l: continue
    t = l.split('|')[1].split(); p = [tuple(map(float, t[i:i + 3])) for i in range(0, len(t), 4)]
    d = min(math.dist(a, b) for a in p[:n] for b in p[n:]) * 0.529177210903
    if 1.6 <= d <= 4.5: rows.append((d, l.rstrip('\n')))
rows.sort()
if len(rows) > want:
    rows = [rows[round(i * (len(rows) - 1) / (want - 1))] for i in range(want)]
for d, l in rows: print(l)
sys.stderr.write(f"{len(rows)} pairs, closest contacts {rows[0][0]:.2f} to {rows[-1][0]:.2f} A\n" if rows else "no pairs in range\n")
