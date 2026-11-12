#!/usr/bin/env python3
"""A law fitted mechanism by mechanism to the engine's own exact parts.

Per pair the engine gives (overlap-water-*.txt, ind/parts-water-*.txt):
  E_es   exact electrostatics of the two Hartree-Fock densities
  E_hl   Heitler-London: E_es + exchange repulsion
  E_hf   Hartree-Fock: + induction and whatever else responds (charge transfer)
  E_tot  MP2: + correlation (dispersion, and the change correlation makes to the electrostatics)
and the model's own cloud electrostatics and induced-dipole energy on the same geometry.

Terms, each against its own target, sites = atoms:
  exchange        sum A_ab exp(-B_ab r)                 target E_hl - E_es
  short-range     -sum A'_ab exp(-B'_ab r)              target (E_hf - E_hl) - model induction
    attraction (what Hartree-Fock has beyond the antisymmetrised molecules that point dipoles do not give)
  the rest        cloud electrostatics - sum f6(B_ab r) C6_ab / r^6, nothing fitted:
                  target E_tot - (E_hl - E_es) - (E_hf - E_hl)  =  E_es + E_corr
Half the pairs (odd index) are held out of every fit.
"""
import sys, math
import numpy as np
KCAL = 627.509474
C6 = {(0, 0): 17.48, (0, 1): 6.33, (1, 0): 6.33, (1, 1): 2.38}       # phys-c6, water (PLAY.md E8b)
sets = [('bias', '../e8c/avg/pairs-water-gpu-bias-mp2.txt'), ('wall', '../e8c/wall-rigid/pairs-water-gpu-wall-mp2.txt')]
rows = []
for tag, pf in sets:
    geo = {}
    for l in open(pf):
        if l.startswith('#') or '|' not in l: continue
        h, t = l.split('|'); t = t.split()
        geo[int(h.split()[0])] = [(np.array(list(map(float, t[i:i+3]))), int(t[i+3])) for i in range(0, len(t), 4)]
    parts = {int(l.split()[0]): (float(l.split()[1]), float(l.split()[2])) for l in open(f'ind/parts-water-{tag}.txt')}
    try: lines = open(f'overlap-water-{tag}.txt').readlines()
    except FileNotFoundError: continue
    for l in lines:
        v = l.split()
        if len(v) != 9: continue
        k = int(v[0]); et, eh, es, s, s1, sv, hl = map(float, v[2:9])
        at = geo[k]; n = len(at) // 2
        d = {}   # distances by type pair
        for pa, ta in at[:n]:
            for pb, tb in at[n:]:
                d.setdefault((min(ta, tb), max(ta, tb)), []).append(np.linalg.norm(pa - pb))
        rows.append(dict(tag=tag, k=k, d={q: np.array(x) for q, x in d.items()}, et=et*KCAL, eh=eh*KCAL, es=es*KCAL, hl=hl*KCAL,
                         es_cloud=parts[k][0]*KCAL, ind=parts[k][1]*KCAL, near=min(min(x) for x in d.values())*0.529177))
types = sorted({q for r in rows for q in r['d']})
print(f"{len(rows)} pairs ({sum(r['tag']=='bias' for r in rows)} liquid-drawn, {sum(r['tag']=='wall' for r in rows)} wall); site pairs {types}")
train = [r for r in rows if r['k'] % 2 == 0]; test = [r for r in rows if r['k'] % 2 == 1]

def design(rs, B):
    return np.array([[np.exp(-B[i] * r['d'][q]).sum() for i, q in enumerate(types)] for r in rs])
def fit_bm(target, wfun, B0, sign=1.0):
    """amplitudes linear (non-negative by clipping and refitting), exponents by a coordinate search."""
    y = np.array([target(r) for r in train]) * sign; w = np.array([wfun(r) for r in train])
    def solve(B):
        X = design(train, B); keep = list(range(len(types)))
        while True:
            A = np.zeros(len(types))
            sol = np.linalg.lstsq(X[:, keep] * w[:, None], y * w, rcond=None)[0]
            A[keep] = sol
            if (sol >= 0).all() or len(keep) == 1: break
            keep = [k for k, a in zip(keep, sol) if a >= 0]
        A = np.maximum(A, 0); return A, float(np.sum((w * (X @ A - y)) ** 2))
    B = np.array(B0, float); A, best = solve(B); step = 0.4
    while step > 0.005:
        moved = False
        for i in range(len(B)):
            for s in (+step, -step):
                Bn = B.copy(); Bn[i] = max(0.3, Bn[i] + s); An, c = solve(Bn)
                if c < best - 1e-12: B, A, best, moved = Bn, An, c, True
        if not moved: step /= 2
    return A * sign, B
def model(rs, A, B): return design(rs, B) @ A
def tt6(x):
    s = sum(x ** k / math.factorial(k) for k in range(7)); return 1 - np.exp(-x) * s
def report(name, target, pred_train, pred_test, band=True):
    for lab, rs, p in (('fitted', train, pred_train), ('held out', test, pred_test)):
        y = np.array([target(r) for r in rs]); e = p - y
        liq = np.array([r['tag'] == 'bias' for r in rs])
        s = f"  {name:22s} {lab:8s} liquid-drawn rms {np.sqrt(np.mean(e[liq]**2)):.3f} mean {np.mean(e[liq]):+.3f} (of rms {np.sqrt(np.mean(y[liq]**2)):.2f})"
        if (~liq).any():
            rel = e[~liq] / np.maximum(np.abs(y[~liq]), 1.0)
            s += f"; wall rms {np.sqrt(np.mean(e[~liq]**2)):.2f}, relative {np.sqrt(np.mean(rel**2))*100:.0f}% (of rms {np.sqrt(np.mean(y[~liq]**2)):.1f})"
        print(s + "  kcal/mol")

# exchange
t_ex = lambda r: r['hl'] - r['es']
w_ex = lambda r: 1.0 / (1.0 + abs(t_ex(r)) / 5.0)
A_ex, B_ex = fit_bm(t_ex, w_ex, [2.0] * len(types))
print("exchange:   " + "  ".join(f"{q}: A {a:.4g} B {b:.3f}" for q, a, b in zip(types, A_ex / KCAL, B_ex)) + "   (hartree, 1/bohr)")
report('exchange', t_ex, model(train, A_ex, B_ex), model(test, A_ex, B_ex))
# short-range attraction
t_ct = lambda r: (r['eh'] - r['hl']) - r['ind']
w_ct = lambda r: 1.0 / (1.0 + abs(t_ct(r)) / 2.0)
A_ct, B_ct = fit_bm(t_ct, w_ct, [1.5] * len(types), sign=-1.0)
print("attraction: " + "  ".join(f"{q}: A {a:.4g} B {b:.3f}" for q, a, b in zip(types, A_ct / KCAL, B_ct)))
report('short-range attraction', t_ct, model(train, A_ct, B_ct), model(test, A_ct, B_ct))
print(f"  (the model's induction alone against Hartree-Fock's: mean ratio {np.mean([r['ind']/(r['eh']-r['hl']) for r in rows if abs(r['eh']-r['hl'])>0.1]):.2f})")
# the rest: nothing fitted
def rest(rs):
    out = []
    for r in rs:
        e = r['es_cloud']
        for i, q in enumerate(types):
            d = r['d'][q]; e -= KCAL * np.sum(tt6(B_ex[i] * d) * C6[q] / d ** 6)
        out.append(e)
    return np.array(out)
t_rest = lambda r: r['et'] - (r['hl'] - r['es']) - (r['eh'] - r['hl'])
report('clouds + dispersion', t_rest, rest(train), rest(test))
# the whole law against MP2
tot = lambda rs: model(rs, A_ex, B_ex) + model(rs, A_ct, B_ct) + np.array([r['ind'] for r in rs]) + rest(rs)
report('the whole law vs MP2', lambda r: r['et'], tot(train), tot(test))
if '-v' in sys.argv:
    print("\n  near(A)  MP2    law   | exch: exact model | beyond-HL: exact  ind+attr | rest: exact model")
    p = dict(zip([id(r) for r in rows], zip(model(rows, A_ex, B_ex), model(rows, A_ct, B_ct), rest(rows))))
    for r in sorted(rows, key=lambda r: r['near']):
        ex, ct, rs_ = p[id(r)]
        print(f"  {r['tag'][0]} {r['near']:.2f} {r['et']:+8.2f} {ex+ct+r['ind']+rs_:+8.2f} | {t_ex(r):+8.2f} {ex:+8.2f} | {r['eh']-r['hl']:+7.2f} {ct+r['ind']:+7.2f} | {t_rest(r):+8.2f} {rs_:+8.2f}")
