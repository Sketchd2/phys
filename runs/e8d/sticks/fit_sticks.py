"""Ball-and-stick water, derived from the engine's own reference (PLAY.md E8d, the owner's idea).

Springs from the engine's monomer surface (surface-water.txt: HF + RI-MP2 over r1, r2, theta);
atom-atom pair law from the engine's MP2 pair energies (point charges on the atoms, an
exponential repulsion and a damped C6 per atom-pair type, nothing else: no clouds, no
induction, no extra sites). Writes params-<variant>.txt for phys-sticks.
"""
import math, sys, json
import numpy as np
from scipy.optimize import minimize, least_squares

HART = 627.509474
BOHR = 0.529177210903

# ---- the monomer ------------------------------------------------------------------------
terms = []
ref = None
for l in open(r'E:\Projects\phys\surface-water.txt'):
    w = l.split()
    if not w or w[0].startswith('#'):
        continue
    if w[0] == 'ref':
        ref = [float(x) for x in w[1:4]]
    elif w[0] == 'term':
        terms.append((int(w[1]), int(w[2]), int(w[3]), float(w[4])))


def surface(r1, r2, th):
    d1, d2, dt = (r1 - ref[0]) / 0.1, (r2 - ref[0]) / 0.1, (th - ref[1]) / 0.2
    return ref[2] + sum(c * d1 ** a * d2 ** b * dt ** k for a, b, k, c in terms)


# 'ref r (bohr), theta (rad), energy': ref = [r_ref, theta_ref, e_ref]
res = minimize(lambda x: surface(x[0], x[0], x[1]), [ref[0], ref[1]], method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-14})
r_min, th_min = res.x
print('surface minimum: r_OH %.5f A, theta %.3f deg' % (r_min * BOHR, math.degrees(th_min)))


def hessian(x0, h=2e-3):
    n = 3
    H = np.zeros((n, n))
    f = lambda v: surface(v[0], v[1], v[2])
    for i in range(n):
        for j in range(n):
            e = np.zeros(n); g = np.zeros(n)
            e[i] = h; g[j] = h
            H[i, j] = (f(x0 + e + g) - f(x0 + e - g) - f(x0 - e + g) + f(x0 - e - g)) / (4 * h * h)
    return H


x0 = np.array([r_min, r_min, th_min])
H = hessian(x0)
print('Hessian in (r1, r2, theta), hartree/bohr^2, /bohr/rad, /rad^2:')
print(H)
# frequencies as a check: harmonic normal modes with water masses
m_o, m_h = 15.999 * 1822.888486, 1.008 * 1822.888486


def freqs(H, r, th):
    # Wilson G matrix for (r1, r2, theta) of a bent XY2
    g_rr = 1 / m_h + 1 / m_o
    g_r1r2 = math.cos(th) / m_o
    g_rt = -math.sin(th) / (m_o * r)
    g_tt = 2 / (m_h * r * r) + 2 / (m_o * r * r) * (1 - math.cos(th))
    G = np.array([[g_rr, g_r1r2, g_rt], [g_r1r2, g_rr, g_rt], [g_rt, g_rt, g_tt]])
    w2 = np.linalg.eigvals(G @ H).real
    cm = 219474.63
    return sorted(math.sqrt(abs(x)) * cm * np.sign(x) for x in w2)


print('harmonic frequencies (cm-1) from the surface Hessian:', ['%.0f' % f for f in freqs(H, r_min, th_min)])

# the shape the pair data were computed at
R0 = 1.8499284924411517
TH0 = math.radians(104.078)

# ---- the pairs --------------------------------------------------------------------------
base = r'E:\Projects\phys\runs\e8c\avg\refit1'
files = ['pairs-water-gpu-mp2.txt', 'pairs-water-gpu-bias-mp2.txt', 'pairs-water-gpu-wall-mp2.txt']
pairs = []
for fn in files:
    for l in open(base + '\\' + fn):
        if l.startswith('#') or '|' not in l:
            continue
        head, tail = l.split('|')
        h = head.split()
        e = float(h[2])
        t = tail.split()
        atoms = [(np.array([float(t[4 * k]), float(t[4 * k + 1]), float(t[4 * k + 2])]), int(t[4 * k + 3])) for k in range(6)]
        pairs.append((e, atoms))
print(len(pairs), 'pairs; E range %.2f .. %.1f kcal/mol' % (min(p[0] for p in pairs) * HART, max(p[0] for p in pairs) * HART))

E = np.array([p[0] for p in pairs])
# atoms: 0-2 molecule A, 3-5 molecule B; types 0 = O, 1 = H
pos = np.array([[a[0] for a in p[1]] for p in pairs])   # (n, 6, 3)
typ = np.array([[a[1] for a in p[1]] for p in pairs])   # (n, 6)
n = len(pairs)
# all A-B atom pairs
ia, ib = np.meshgrid(range(3), range(3, 6), indexing='ij')
ia, ib = ia.ravel(), ib.ravel()
r = np.linalg.norm(pos[:, ia, :] - pos[:, ib, :], axis=2)       # (n, 9)
ta, tb = typ[:, ia], typ[:, ib]
kind = (ta + tb)  # 0 = OO, 1 = OH, 2 = HH
tt = np.zeros((n, 9), dtype=int)
tt[(ta == 0) & (tb == 0)] = 0
tt[((ta == 0) & (tb == 1)) | ((ta == 1) & (tb == 0))] = 1
tt[(ta == 1) & (tb == 1)] = 2


def tang_toennies6(x):
    s = 1.0; term = 1.0
    for k in range(1, 7):
        term *= x / k
        s += term
    return 1 - np.exp(-x) * s


def model(params, qh):
    # params: logA[3], B[3], C6[3]
    logA, B, C6 = params[0:3], params[3:6], params[6:9]
    qo = -2 * qh
    q = np.where(typ == 0, qo, qh)
    qq = q[:, ia] * q[:, ib]
    coul = (qq / r).sum(axis=1)
    A = np.exp(logA)[tt]; Bk = B[tt]; C = C6[tt]
    rep = (A * np.exp(-Bk * r)).sum(axis=1)
    disp = (tang_toennies6(Bk * r) * C / r ** 6).sum(axis=1)
    return coul + rep - disp


# weights: the low-energy pairs a liquid is made of, as the project does, and the wall by relative error
low = E.min()
w = 1.0 / np.sqrt(1.0 + (np.maximum(E - low, 0) * HART / 4.0) ** 2)
held = np.arange(n) % 5 == 0
train = ~held


def fit(variant, qh_fixed=None, c6_fixed=None):
    def unpack(p):
        qh = qh_fixed if qh_fixed is not None else p[-1]
        pp = np.concatenate([p[:6], np.array(c6_fixed)]) if c6_fixed is not None else p[:9]
        return pp, qh

    def resid(p):
        pp, qh = unpack(p)
        return (model(pp, qh) - E) * w * HART

    # start: A ~ 1, B ~ 1.5..2.5, C6 ~ 5
    p0 = np.array([math.log(5.0), math.log(1.0), math.log(1.0), 1.4, 1.6, 1.8] + ([] if c6_fixed is not None else [20.0, 5.0, 2.0]) + ([] if qh_fixed is not None else [0.35]))
    best = None
    rng = np.random.default_rng(1)
    for trial in range(6):
        start = p0 + rng.normal(0, 0.3, size=p0.shape) * (trial > 0)
        try:
            sol = least_squares(lambda p: resid(p)[train], start, method='trf', max_nfev=4000)
        except Exception as ex:
            continue
        if best is None or sol.cost < best.cost:
            best = sol
    pp, qh = unpack(best.x)
    pred = model(pp, qh)
    err = (pred - E) * HART
    near = (E - low) * HART < 3.0
    print('%s: q_H %.4f  A %s  B %s  C6 %s' % (variant, qh, np.round(np.exp(pp[0:3]), 3), np.round(pp[3:6], 3), np.round(pp[6:9], 3)))
    print('   rms within 3 kcal/mol of the lowest: fitted %.3f (n=%d)  held out %.3f (n=%d); all pairs below 20 kcal/mol: %.3f' % (
        np.sqrt(np.mean(err[near & train] ** 2)), (near & train).sum(), np.sqrt(np.mean(err[near & held] ** 2)), (near & held).sum(),
        np.sqrt(np.mean(err[(E * HART < 20)] ** 2))))
    # the deepest dimer and the long-range sums
    cp = 'params-%s.txt' % variant
    with open(r'E:\Projects\phys\runs\e8d\sticks\\' + cp, 'w') as f:
        f.write('qh %r\n' % qh)
        f.write('A %s\n' % ' '.join(repr(float(x)) for x in np.exp(pp[0:3])))
        f.write('B %s\n' % ' '.join(repr(float(x)) for x in pp[3:6]))
        f.write('C6 %s\n' % ' '.join(repr(float(x)) for x in pp[6:9]))
        f.write('r0 %r\n' % R0)
        f.write('theta0 %r\n' % TH0)
        f.write('hessian %s\n' % ' '.join(repr(float(x)) for x in H.ravel()))
    return best


# the dipole the monomer's own MP2 gives (1.873 D from phys-esp in the E8c averaged-shape pipeline)
mu_au = 1.873 * 0.393430
qh_gas = mu_au / (2 * R0 * math.cos(TH0 / 2))
print('gas-phase charge from the monomer dipole: q_H = %.4f e' % qh_gas)
E6 = [17.48, 6.331, 2.376]
if len(sys.argv) > 1 and sys.argv[1] == 'e6':
    fit('e6disp', qh_gas, E6)
else:
    fit('gascharge', qh_gas)
    fit('freecharge', None)
