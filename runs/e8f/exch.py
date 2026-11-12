#!/usr/bin/env python3
"""Exchange repulsion (Heitler-London less exact electrostatics) against three
density overlaps, species by species. Reads overlap-*.txt, columns
k sep E_total E_hf E_es S S_no1s S_valence E_hl."""
import glob, math, sys
import numpy as np
KCAL = 627.509474
allr = {}
for f in sorted(glob.glob('overlap-*.txt')):
    sp = f[len('overlap-'):-4]; rows = []
    for l in open(f):
        v = l.split()
        if len(v) != 9: continue
        et, eh, es, s, s1, sv, hl = map(float, v[2:9])
        rows.append((s, s1, sv, (hl - es) * KCAL, (eh - hl) * KCAL, es * KCAL, (et - eh) * KCAL))
    if rows: allr[sp] = np.array(sorted(rows))
names = ['all electrons', 'without 1s', 'valence (no 1s, no N/O/F 2s)']
SMIN = 5e-5
for col in range(3):
    print(f"\nexchange / overlap, {names[col]}  [hartree bohr^3]; pairs with whole overlap above {SMIN:g}")
    lx, ly, who = [], [], []
    for sp, a in allr.items():
        m = a[:, 0] > SMIN; k = a[m, 3] / KCAL / a[m, col]
        if m.sum() < 2: continue
        g = np.polyfit(np.log(a[m, col]), np.log(a[m, 3]), 1)[0]
        print(f"  {sp:14s} n {m.sum():3d}  K mean {k.mean():7.2f}  spread {k.std()/k.mean()*100:5.1f}%  range {k.min():6.2f}-{k.max():6.2f}  log-log slope {g:.3f}")
        lx += list(np.log(a[m, col])); ly += list(np.log(a[m, 3] / KCAL)); who += [sp] * m.sum()
    if len(lx) > 4:
        lx, ly = np.array(lx), np.array(ly)
        g, lk = np.polyfit(lx, ly, 1); res = ly - g * lx - lk
        print(f"  one law for all: E = {math.exp(lk):.3f} S^{g:.3f}; rms error {np.std(res)*100:.1f}% of the exchange; by species: " + ", ".join(f"{sp} {np.mean(res[np.array(who)==sp])*100:+.1f}%" for sp in allr if (np.array(who)==sp).any()))
        res1 = ly - lx - np.mean(ly - lx)
        print(f"  with the power held at 1: K = {math.exp(np.mean(ly-lx)):.3f}; rms error {np.std(res1)*100:.1f}%; by species: " + ", ".join(f"{sp} {np.mean(res1[np.array(who)==sp])*100:+.1f}%" for sp in allr if (np.array(who)==sp).any()))
if '-v' in sys.argv:
    for sp, a in allr.items():
        print(f"\n{sp}:  S_all  S_no1s  S_val  exchange  HF-beyond  E_es  E_corr [kcal/mol]")
        for r in a: print("  " + " ".join(f"{x:.3e}" for x in r[:3]) + " " + " ".join(f"{x:+8.3f}" for x in r[3:]))
