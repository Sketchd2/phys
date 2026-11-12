#!/usr/bin/env python3
"""Is the repulsion one constant times the density overlap, for every species?
Reads overlap-*.txt (k sep E_total E_hf E_es S). The Hartree-Fock interaction
less the exact electrostatics is exchange repulsion plus Hartree-Fock's
induction (negative), so the ratio (E_hf - E_es) / S is a lower bound on the
repulsion's constant and is cleanest where the overlap is large."""
import glob, math, sys
import numpy as np
KCAL = 627.509474
rows = {}
for f in sorted(glob.glob('overlap-*.txt')):
    sp = f[len('overlap-'):-4]
    for l in open(f):
        k, sep, et, eh, es, s = l.split()
        rows.setdefault(sp, []).append((float(sep), float(et), float(eh), float(es), float(s)))
print("species            n   bands of overlap S: median (E_hf - E_es)/S [hartree bohr^3] (count)")
edges = [1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 1e-1]
allx, ally, lab = [], [], []
for sp, r in rows.items():
    a = np.array(r); rep = a[:, 2] - a[:, 3]; S = a[:, 4]
    out = []
    for lo, hi in zip(edges, edges[1:]):
        m = (S >= lo) & (S < hi)
        out.append(f"{np.median(rep[m]/S[m]):6.2f} ({m.sum():2d})" if m.sum() else "      -     ")
    print(f"{sp:16s} {len(r):3d}  " + "  ".join(out))
    m = (S > 1e-4) & (rep > 0)
    allx += list(np.log(S[m])); ally += list(np.log(rep[m])); lab += [sp] * m.sum()
print("bands:              " + "  ".join(f"{lo:.0e}-{hi:.0e}" for lo, hi in zip(edges, edges[1:])))
x, y = np.array(allx), np.array(ally)
if len(x) > 3:
    g, lk = np.polyfit(x, y, 1)
    print(f"\nall species, S > 1e-4: E_rep = K S^g with g = {g:.3f}, K = {math.exp(lk):.3f}; rms of ln(ratio) {np.std(y - g*x - lk):.3f}")
    k1 = math.exp(np.mean(y - x))
    print(f"with g = 1 held: K = {k1:.3f}, rms of ln(ratio) {np.std(y - x):.3f}")
    for sp in rows:
        m = np.array([l == sp for l in lab])
        if m.sum(): print(f"  {sp:16s} mean ln(E_rep / (K S^g)) {np.mean((y - g*x - lk)[m]):+.3f}  (g = 1: {np.mean((y - x)[m]) - math.log(k1):+.3f})  n {m.sum()}")
