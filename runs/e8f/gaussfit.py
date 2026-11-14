#!/usr/bin/env python3
"""A molecule's density as Gaussians on its atoms and on its own bond and
lone-pair centres (the Boys centroids phys-density prints), by linear least
squares from the density itself, then checked against exact pair overlaps.

usage: gaussfit.py NAME "w1,w2,..." [--power p] [--col c] [--atoms-only] [--check PAIRS OVERLAP [excol]]...
  widths in bohr; fit weight 1/rho^p (default 1); col 6 = valence density, 5 = without 1s, 4 = whole
A Gaussian of width s is (pi s^2)^(-3/2) exp(-r^2/s^2), the simulator's cloud.
Sites alike by symmetry (the same sorted distances to the atoms) share amplitudes."""
import sys, numpy as np
a = sys.argv[1:]; name = a[0]; widths = [float(x) for x in a[1].split(',')]
power = float(a[a.index('--power') + 1]) if '--power' in a else 1.0
col = int(a[a.index('--col') + 1]) if '--col' in a else 6
atoms = []; centres = []; pts = []
for l in open(f'density-{name}.txt'):
    v = l.split()
    if v[0] == 'atom': atoms.append((int(v[1]), np.array(list(map(float, v[2:5])))))
    elif v[0] == 'centre': centres.append(np.array(list(map(float, v[1:4]))))
    elif v[0] != 'orbital': pts.append(list(map(float, v)))
P = np.array(pts); xyz, w, rho = P[:, :3], P[:, 3], P[:, col]
ref = np.array([r for _, r in atoms])
sites = [(r, ('atom', z)) for z, r in atoms]
if '--atoms-only' not in a:
    ext = [float(x) for x in a[a.index('--extend') + 1].split(',')] if '--extend' in a else []
    for c in centres:
        hv = min((r for z, r in atoms if z > 2), key=lambda r: np.linalg.norm(c - r))
        for f in [1.0] + ext:      # --extend f,...: more centres on the line from the heavy atom through each centroid
            p = hv + f * (c - hv)
            key = tuple(sorted((z, round(float(np.linalg.norm(p - r)), 1)) for z, r in atoms))
            sites.append((p, ('centre', key)))
classes = []
for _, k in sites:
    if k not in classes: classes.append(k)
cols = [(k, s) for k in classes for s in widths]
G = np.zeros((len(P), len(cols)))
for j, (k, s) in enumerate(cols):
    for r, kk in sites:
        if kk == k: G[:, j] += (np.pi * s * s) ** -1.5 * np.exp(-np.sum((xyz - r) ** 2, axis=1) / (s * s))
rmax = float(a[a.index("--rmax") + 1]) if "--rmax" in a else 1e9
keep = (rho > 1e-7) & (rho < rmax)   # --rmax: fit only the tail, where molecules meet
wt = np.sqrt(w[keep] / np.maximum(rho[keep], 1e-7) ** power)
if '--moments' in a:
    # hold the model to the density's own charge, dipole and second moments (rows of great weight):
    # a pair's electrostatic energy at range is these, and a least-squares fit of the density alone
    # leaves them a few hundredths of an electron out
    f = [np.ones(len(P))] + [xyz[:, i] for i in range(3)] + [xyz[:, i] * xyz[:, j] for i in range(3) for j in range(i, 3)]
    rows = []; rhs = []
    for fi in f:
        row = np.zeros(len(cols))
        for j, (k, s) in enumerate(cols):
            for r, kk in sites:
                if kk == k:
                    # the moment of a normalised Gaussian at r of width s: 1, r_i, r_i r_j + delta_ij s^2 / 2
                    row[j] += 0.0
        rows.append(row); rhs.append(np.sum(w * rho * fi))
    idx = 0
    def gm(r, s):
        out = [1.0] + list(r) + [r[i] * r[j] + (s * s / 2 if i == j else 0.0) for i in range(3) for j in range(i, 3)]
        return np.array(out)
    M = np.zeros((10, len(cols)))
    for j, (k, s) in enumerate(cols):
        for r, kk in sites:
            if kk == k: M[:, j] += gm(r, s)
    big = 1e3
    c = np.linalg.lstsq(np.vstack([G[keep] * wt[:, None], big * M]), np.concatenate([rho[keep] * wt, big * np.array(rhs)]), rcond=None)[0]
else:
    c = np.linalg.lstsq(G[keep] * wt[:, None], rho[keep] * wt, rcond=None)[0]
fit = G @ c
amp = {ks: x for ks, x in zip(cols, c)}
print(f"{name}, {len(classes)} kinds of site x {len(widths)} widths: {np.sum(w*rho):.4f} electrons, model {np.sum(w*fit):.4f}; integral |model - rho| = {np.sum(w*np.abs(fit-rho)):.4f}")
model = [(r, s, amp[(k, s)]) for r, k in sites for s in widths]        # (reference position, width, electrons)
np.save(f'gauss-{name}.npy', np.array([list(r) + [s, q] for r, s, q in model]))
def place(at):
    """the model's sites carried onto a molecule's atoms `at` (rigid fit of the reference atoms)"""
    at = np.array(at); ca, cr = at.mean(0), ref.mean(0)
    U, _, Vt = np.linalg.svd((ref - cr).T @ (at - ca)); d = np.sign(np.linalg.det(U @ Vt))
    R = (U @ np.diag([1, 1, d]) @ Vt)
    return [((r - cr) @ R + ca, s, q) for r, s, q in model]
i = 0
while '--check' in a[i:]:
    i = a.index('--check', i) + 1; pf, of = a[i], a[i + 1]
    excol = int(a[i + 2]) if len(a) > i + 2 and a[i + 2].isdigit() else {6: 7, 5: 6, 4: 5}[col]
    geo = {}
    for l in open(pf):
        if l.startswith('#') or '|' not in l: continue
        h, t = l.split('|'); t = t.split(); geo[int(h.split()[0])] = [np.array(list(map(float, t[q:q+3]))) for q in range(0, len(t), 4)]
    ex, mo = [], []; n = len(atoms)
    for l in open(of):
        v = l.split()
        if len(v) != 9: continue
        g = geo[int(v[0])]; A, B = place(g[:n]), place(g[n:]); S = 0.0
        for pa, sa, qa in A:
            for pb, sb, qb in B:
                s2 = sa * sa + sb * sb; S += qa * qb * (np.pi * s2) ** -1.5 * np.exp(-np.sum((pa - pb) ** 2) / s2)
        ex.append(float(v[excol])); mo.append(S)
    ex, mo = np.array(ex), np.array(mo); m = ex > 3e-5; r = mo[m] / ex[m]
    print(f"   {of}: {m.sum()} pairs, model/exact overlap mean {r.mean():.3f}, spread {r.std()/r.mean()*100:.1f}%, range {r.min():.2f}-{r.max():.2f}")
