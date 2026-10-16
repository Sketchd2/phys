"""The no-wall law (good liquid-range pressure, holes inside 1 A) plus a steep Born-Mayer core fitted only to what
the law leaves over on the wall pairs (PLAY.md E8d). Everything else in the law is held; the dispersion's damping
exponent is written out as 2.0, the value it had implicitly, so a steep core exponent cannot change it."""
import subprocess, sys, math, numpy as np
from scipy.optimize import nnls

HART = 627.509474
BOHR = 0.529177210903
base = r'E:\Projects\phys\runs\e8d\core'
law_in = r'E:\Projects\phys\runs\e8c\avg\law-final.txt'
exe = r'E:\Projects\phys\target\dev9\release\phys-lawenergy.exe'
sets = ['-gpu-mp2', '-gpu-bias-mp2', '-gpu-wall-mp2', '-liq-gpu-mp2']

rows = []   # (E_ref, E_law, distances by type pair (OO, OH, HH) as arrays of r in bohr)
for tag in sets:
    out = subprocess.run([exe, law_in, tag], capture_output=True, text=True, cwd=base).stdout.splitlines()
    law_e = {}
    for l in out:
        w = l.split()
        law_e[w[0]] = (float(w[2]), float(w[3]))
    for l in open(base + r'\pairs-water%s.txt' % tag):
        if l.startswith('#') or '|' not in l:
            continue
        head, tail = l.split('|')
        k = head.split()[0]
        t = tail.split()
        atoms = [(np.array([float(t[4 * i]), float(t[4 * i + 1]), float(t[4 * i + 2])]), int(t[4 * i + 3])) for i in range(6)]
        d = {0: [], 1: [], 2: []}
        for a in range(3):
            for b in range(3, 6):
                r = np.linalg.norm(atoms[a][0] - atoms[b][0])
                d[atoms[a][1] + atoms[b][1]].append(r)
        e_ref, e_law = law_e[k]
        rows.append((e_ref, e_law, {kk: np.array(v) for kk, v in d.items()}))
print(len(rows), 'pairs')
e_ref = np.array([r[0] for r in rows]) * HART
e_law = np.array([r[1] for r in rows]) * HART
res = e_ref - e_law
print('no-wall law: rms error on pairs below 10 kcal/mol %.3f; on pairs above 30 kcal/mol mean %.1f (reference minus law)' % (
    np.sqrt(np.mean(res[e_ref < 10] ** 2)), np.mean(res[e_ref > 30])))
w = 1.0 / (np.abs(e_ref) + 5.0)
best = None
for B in [4.0, 5.0, 6.0, 7.0, 8.0]:
    phi = np.array([[np.exp(-B * r[2][k]).sum() for k in (0, 1, 2)] for r in rows])
    # scale columns for conditioning
    scale = np.maximum(phi.max(axis=0), 1e-300)
    sol, _ = nnls((phi / scale * w[:, None]), res * w)
    A = sol / scale
    fit = phi @ A
    err = res - fit
    liquid = e_ref < 10
    print('B %.1f: A (hartree... kcal) %s ; rms residual on pairs below 10 kcal/mol %.3f (core adds up to %.2f kcal there), on pairs above 30: rel %.1f%%' % (
        B, np.array2string(A, precision=3), np.sqrt(np.mean(err[liquid] ** 2)), np.max(np.abs(fit[liquid])), 100 * np.mean(np.abs(err[e_ref > 30]) / e_ref[e_ref > 30])))
    cost = np.sum((err * w) ** 2)
    if best is None or cost < best[0]:
        best = (cost, B, A.copy())
print('chosen B', best[1])
B = float(sys.argv[1]) if len(sys.argv) > 1 else best[1]
phi = np.array([[np.exp(-B * r[2][k]).sum() for k in (0, 1, 2)] for r in rows])
scale = np.maximum(phi.max(axis=0), 1e-300)
sol, _ = nnls((phi / scale * w[:, None]), res * w)
A = sol / scale / HART      # to hartree
lines = open(law_in).read().splitlines()
out = []
names = {(0, 0): 0, (0, 1): 1, (1, 1): 2}
for l in lines:
    p = l.split()
    if p and p[0] == 'pair' and (int(p[1]), int(p[2])) in names:
        p[3] = repr(float(A[names[(int(p[1]), int(p[2]))]]))
        p[4] = repr(B)
        l = ' '.join(p)
    out.append(l)
out += ['damp 0 0 2.0', 'damp 0 1 2.0', 'damp 1 1 2.0']
name = r'E:\Projects\phys\runs\e8d\core\law-core-B%g.txt' % B
open(name, 'w').write('\n'.join(out) + '\n')
print('written', name, 'A (hartree)', A)
