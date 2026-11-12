#!/usr/bin/env python3
"""The liquid's excess free energy from a path of bulk states, and the vapour
pressure it implies.

usage: path.py DIR MOLAR_MASS T_HOT RHO_REF "T:B2,..." [--no-tail]
       (B2 in cm3/mol, from phys-b2 with the same law)

Every state is a phys-bulk run in DIR/s-<T>-<g/cm3>/.  The path:
  (1) the isotherm T_HOT (supercritical) from the dilute gas to RHO_REF,
        F_ex/NkT = integral_0^rho (Z - 1) d rho / rho ,
      with the integrand's value at rho = 0 pinned to B2;
  (2) the isochore RHO_REF from T_HOT down,  d(F_ex/NkT)/d beta = U/N;
  (3) at a temperature that has runs at several densities, along it to the
      density rho_0 where the liquid's pressure is the vapour's,
        d(F_ex/N) = P_ex d rho / rho^2 .
      With one density only, the liquid is taken to zero pressure to first
      order (the Gibbs correction cancels P/rho), which is good while
      kappa P^2 / 2 rho is small and is printed as "1st order".
Then  mu_ex = F_ex/N + P_ex/rho_0,  and the vapour by its second virial
coefficient:  ln rho_V + 2 B2 rho_V = ln rho_0 + mu_ex/kT,
              p = rho_V kT (1 + B2 rho_V).
U and P carry phys-bulk's dispersion tail beyond the cut-off unless --no-tail
(a truncated law compared with a slab that has no tail either).
"""
import sys, re, glob, os, math
import numpy as np
R = 1.987204259e-3       # kcal/mol/K
RBAR = 83.14462618       # bar cm3/mol/K
BARCM3 = R / RBAR        # kcal/mol per bar cm3/mol
notail = '--no-tail' in sys.argv
argv = [a for a in sys.argv if a != '--no-tail']
d, M, Thot, rhoR = argv[1], float(argv[2]), float(argv[3]), float(argv[4])
B2 = {float(a.split(':')[0]): float(a.split(':')[1]) for a in argv[5].split(',')}
rng = np.random.default_rng(1)

def blocks(a, nb=8):
    n = len(a) // nb
    m = [np.mean(a[i*n:(i+1)*n]) for i in range(nb)]
    return float(np.mean(a)), float(np.std(m, ddof=1) / math.sqrt(nb))

st = {}
for f in glob.glob(os.path.join(d, 's-*', 'liquid.log')):
    dd = os.path.dirname(f)
    if not os.path.exists(os.path.join(dd, 'done')): continue
    _, T, rho = os.path.basename(dd).split('-'); T, rho = float(T), float(rho)
    tail = re.search(r'cut-off \(C6 [\d.]+\): (-?[\d.]+) bar and (-?[\d.]+) kcal', open(f).read())
    tp, tu = (0.0, 0.0) if notail else (float(tail.group(1)), float(tail.group(2)))
    a = np.loadtxt([n for n in glob.glob(os.path.join(dd, 'bulk-*-*.txt')) if not n.endswith('gOO.txt')][0])
    # the excess pressure sample by sample: the kinetic part at the sample's own translational temperature removed
    pex, dpex = blocks(a[:, 2] - rho / M * RBAR * a[:, 3])
    u, du = blocks(a[:, 1])
    st[(T, rho)] = dict(U=u + tu, dU=du, Pex=pex + tp, dP=dpex)

# (1) the isotherm
iso = sorted((rho, v) for (T, rho), v in st.items() if T == Thot)
rhos = np.array([r for r, _ in iso]); vm = M / rhos
y = np.array([v['Pex'] for _, v in iso]) * vm / (RBAR * Thot) / rhos
dy = np.array([v['dP'] for _, v in iso]) * vm / (RBAR * Thot) / rhos
x = np.concatenate([[0.0], rhos]); yy = np.concatenate([[B2[Thot] / M], y])
w = np.concatenate([[0.0], 1 / np.maximum(dy, 1e-3 * abs(y).max())]); w[0] = 10 * w[1:].max()
def iso_F(yv, deg=5):
    ci = np.polyint(np.polyfit(x, yv, deg, w=w)); return np.polyval(ci, rhoR) - np.polyval(ci, 0.0)
print(f"isotherm {Thot:g} K, {len(rhos)} states, B2/M = {B2[Thot]/M:.4f} cm3/g at rho = 0:")
for deg in (3, 4, 5, 6):
    c = np.polyfit(x, yy, deg, w=w)
    chi = np.sqrt(np.mean(((np.polyval(c, rhos) - y) / dy) ** 2))
    print(f"  degree {deg}: F_ex/NkT = {iso_F(yy, deg):+.4f}, residual {chi:.2f} sigma")
F_hot = iso_F(yy); dF_hot = np.std([iso_F(np.concatenate([[yy[0]], y + rng.normal(0, dy)])) for _ in range(400)])
print(f"  taken (degree 5): F_ex/NkT({Thot:g} K, {rhoR}) = {F_hot:+.4f} +- {dF_hot:.4f}")

# (2) the isochore
ch = sorted((T, v) for (T, rho), v in st.items() if abs(rho - rhoR) < 1e-9)
Ts = np.array([T for T, _ in ch]); beta = 1 / (R * Ts); U = np.array([v['U'] for _, v in ch]); dU = np.maximum(np.array([v['dU'] for _, v in ch]), 1e-4)
deg = min(4, len(Ts) - 2); bh = 1 / (R * Thot)
def chore(Uv, b):
    ci = np.polyint(np.polyfit(beta, Uv, deg, w=1 / dU)); return np.polyval(ci, b) - np.polyval(ci, bh)
print(f"isochore {rhoR} g/cm3, {len(Ts)} states; U/N fit residuals (kcal/mol): " + " ".join(f"{a:+.4f}" for a in np.polyval(np.polyfit(beta, U, deg, w=1/dU), beta) - U))

# (3) and the vapour
print("\n   T    rho_0   F_ex/N    mu_ex    +-     p_sat(atm)  rho_V(g/cm3)   [kcal/mol]")
out = []
for T, v in ch:
    if T not in B2 or T == Thot: continue
    b = 1 / (R * T)
    bF = F_hot + chore(U, b)
    dbF = math.hypot(dF_hot, np.std([chore(U + rng.normal(0, dU), b) for _ in range(200)]))
    leg = sorted((rho, s) for (TT, rho), s in st.items() if TT == T)
    note = ""
    if len(leg) >= 3:
        lr = np.array([r for r, _ in leg]); lp = np.array([s['Pex'] for _, s in leg]); ldp = np.array([s['dP'] for _, s in leg])
        def solve(lpv, p_target):
            c = np.polyfit(lr, lpv, 2, w=1 / ldp)
            tot = np.polyadd(c, [RBAR * T / M, -p_target])      # P_ex(rho) + rho R T / M - p
            r0 = [r.real for r in np.roots(tot) if abs(r.imag) < 1e-9 and lr.min() * 0.8 < r.real < lr.max() * 1.2]
            if not r0: return None
            r0 = min(r0, key=lambda r: abs(r - rhoR))
            g = np.linspace(rhoR, r0, 400); dFl = np.trapezoid(np.polyval(c, g) / g ** 2, g) * M * BARCM3   # kcal/mol
            return r0, dFl, np.polyval(c, r0)
        s0 = solve(lp, 0.0)
        if s0 is None: continue
        rho0, dFl, pex0 = s0
        sp = [solve(lp + rng.normal(0, ldp), 0.0) for _ in range(200)]
        sp = [(a[1] + a[2] * M / a[0] * BARCM3) for a in sp if a]
        dbF = math.hypot(dbF, np.std(sp) * b)
    else:
        rho0, dFl, pex0 = rhoR, 0.0, -rhoR * RBAR * T / M          # first order: liquid at zero pressure
        note = f" 1st order (liquid at {v['Pex'] + rhoR*RBAR*T/M:+.0f} bar)"
    F = bF / b + dFl
    mu = F + pex0 * M / rho0 * BARCM3
    rn = rho0 / M; rv = rn * math.exp(mu * b)
    try:
        for _ in range(100): rv = rn * math.exp(mu * b - 2 * B2[T] * rv)
    except OverflowError:
        continue
    if abs(B2[T] * rv) > 0.35: note += " (B2 rho_V %.2f: virial series doubtful)" % (B2[T]*rv)        # the virial series is no use this near the critical point
    p = rv * (1 + B2[T] * rv) * RBAR * T / 1.01325
    print(f"{T:6.1f} {rho0:.4f} {F:+.4f} {mu:+.4f} {dbF/b:.4f}  {p:9.4f} x/{math.exp(dbF):.3f}  {rv*M:.5f}{note}")
    out.append((T, p))
if len(out) >= 2:
    lt = np.array([1 / T for T, _ in out]); lp = np.log([p for _, p in out])
    k = np.polyfit(lt, lp, min(2, len(out) - 1))
    roots = [1 / r.real for r in np.roots(k) if abs(r.imag) < 1e-12 and r.real > 0 and min(Ts) * 0.8 < 1 / r.real < max(Ts)]
    print("ln p against 1/T crosses one atmosphere at", ", ".join(f"{r:.1f} K" for r in roots))
