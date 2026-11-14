#!/usr/bin/env python3
"""The break-down tests: frozen orbitals against the ghost-basis rows, the uncorrected pair
energy against the counterpoise one, and what raising the table cutoff does. Run in runs/e8f/break."""
import numpy as np, glob, os
H = 627.509474
def rows(f):
    d = {}
    for l in open(f):
        t = l.split()
        try: d[int(t[0])] = [float(x) for x in t[1:]]
        except ValueError: pass
    return d
def timing(f):
    p = f[:-4] + '.time'
    if not os.path.exists(p): return 'still running'
    t = open(p).read().split(); return f"{int(t[0]) / max(1, int(t[1])):.0f} s a pair, peak {t[2]} MB"
print('1. frozen orbitals against the ghost-basis rows')
for sp, ghost in (('water', ['../overlap-water-bias.txt']), ('ammonia', ['../overlap-ammonia.txt']), ('methane', ['../overlap-methane.txt'])):
    f = f'r-{sp}-frozen-1e-14.txt'
    if not os.path.exists(f): continue
    fr = rows(f); gh = {}
    for g in ghost: gh.update(rows(g))
    k = [i for i in fr if i in gh and len(gh[i]) >= 8 and len(fr[i]) >= 8]
    if not k: print(f'   {sp}: no pairs in both yet'); continue
    a = np.array([gh[i][:8] for i in k]); b = np.array([fr[i][:8] for i in k])
    es = (b[:, 3] - a[:, 3]) * H; exg = (a[:, 7] - a[:, 3]) * H; exf = (b[:, 7] - b[:, 3]) * H; sv = b[:, 6] / a[:, 6] - 1
    K = lambda S, x: (S * H * x).sum() / ((S * H) ** 2).sum()
    print(f"   {sp}: {len(k)} pairs. electrostatics rms {es.std():.4f} mean {es.mean():+.4f} kcal/mol; exchange rms {(exf-exg).std():.4f} mean {(exf-exg).mean():+.4f} of {np.sqrt((exg**2).mean()):.2f}; valence overlap {sv.mean():+.4f} +- {sv.std():.4f}; K {K(a[:,6],exg):.3f} -> {K(b[:,6],exf):.3f}; {timing(f)}")
print('2. the uncorrected pair energy against the counterpoise one, and the cutoff')
for sp in ('water', 'ammonia', 'methane', 'methanol'):
    base = None
    for c in ('1e-14', '1e-10', '1e-8', '1e-6', '1e-5', '1e-4', '1e-3', '1e-2'):
        f = f'r-{sp}-full-{c}.txt'
        if not os.path.exists(f): continue
        r = {i: v for i, v in rows(f).items() if len(v) >= 10}
        if not r: print(f'   {sp} cutoff {c}: no row (the solve did not converge, or still running)'); continue
        if base is None:
            base = r; cp = np.array([v[1] for v in r.values()]) * H; un = np.array([v[9] for v in r.values()]) * H
            print(f"   {sp}, {len(r)} pairs, counterpoise / uncorrected kcal/mol: " + ', '.join(f'{a:+.3f}/{b:+.3f}' for a, b in zip(cp, un)) + f"; mean difference {np.mean(un - cp):+.3f}; {timing(f)}")
        else:
            both = [i for i in r if i in base]
            d = max(abs(r[i][9] - base[i][9]) for i in both) * H; dh = max(abs(r[i][7] - base[i][7]) for i in both) * H
            print(f"   {sp} cutoff {c}: {len(both)} pairs, largest change uncorrected total {d:.1e}, Heitler-London {dh:.1e} kcal/mol; {timing(f)}")
print('3. methanol')
for f in sorted(glob.glob('r-methanol-*.txt')):
    r = rows(f); print(f"   {f}: {len(r)} pairs; {timing(f)}; " + ' '.join(f"pair {i}: E_hl {v[7]*H:+.4f}" for i, v in r.items()))
