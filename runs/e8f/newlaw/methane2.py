#!/usr/bin/env python3
"""Methane with the one-constant law, below 148 K. At 0.4226 g/cm3 the law's liquid is under
tension that the box does not hold below about 140 K (the pressure is flat at -400 to -670 bar
from 0.36 to 0.44 g/cm3: it has opened a cavity), so the 0.4226 isochore is not a path there.
From path.py's free energy of the 148 K liquid at zero pressure: along 148 K to 0.48 g/cm3,
down the 0.48 isochore to 134 and 122 K, and back to each temperature's zero-pressure density.
usage: methane2.py DIR F148 RHO148 "T:B2,..."   (F148 = path.py's F_ex/N at 148 K, kcal/mol)"""
import sys, re, glob, os, math, numpy as np
d, F148, r148 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]); B2 = {float(a.split(':')[0]): float(a.split(':')[1]) for a in sys.argv[4].split(',')}
M = 16.043; R = 1.987204e-3; RBAR = 83.14462; BARCM3 = 2.390057e-5; RI = 0.48
st = {}
for f in glob.glob(os.path.join(d, 's-*', 'liquid.log')):
    dd = os.path.dirname(f)
    if not os.path.exists(dd + '/done'): continue
    _, T, rho = os.path.basename(dd).split('-'); T, rho = float(T), float(rho)
    t = re.search(r'cut-off \(C6 [\d.]+\): (-?[\d.]+) bar and (-?[\d.]+) kcal', open(f).read())
    a = np.loadtxt([n for n in glob.glob(dd + '/bulk-*-*.txt') if not n.endswith('gOO.txt')][0])
    st[(T, rho)] = (a[:, 1].mean() + float(t.group(2)), (a[:, 2] - rho / M * RBAR * a[:, 3]).mean() + float(t.group(1)))
def leg(T, lo=0.455):
    pts = sorted((rho, v[1]) for (TT, rho), v in st.items() if TT == T and rho > lo)
    return np.polyfit([p[0] for p in pts], [p[1] for p in pts], min(2, len(pts) - 1)), [p[0] for p in pts]
def dF(c, a, b):       # integral of P_ex / rho^2 from a to b, kcal/mol
    g = np.linspace(a, b, 400); return np.trapezoid(np.polyval(c, g) / g ** 2, g) * M * BARCM3
def rho0(c, T, near):
    tot = np.polyadd(c, [RBAR * T / M, 0.0]); r = [x.real for x in np.roots(tot) if abs(x.imag) < 1e-9 and 0.3 < x.real < 0.6]
    return min(r, key=lambda x: abs(x - near))
c148, used = leg(148, 0.42); print('148 K leg densities', used)
F = {148.0: F148 + dF(c148, r148, RI)}
Ts = [148.0, 134.0, 122.0]; beta = np.array([1 / (R * T) for T in Ts]); U = np.array([st[(T, RI)][0] for T in Ts])
ci = np.polyint(np.polyfit(beta, U, 2))
for T, b in zip(Ts[1:], beta[1:]): F[T] = (F[148.0] * beta[0] + np.polyval(ci, b) - np.polyval(ci, beta[0])) / b
print('   T    rho_0   F_ex/N    mu_ex   p_sat(atm)')
out = []
for T in Ts:
    c, used = leg(T, 0.42 if T == 148 else 0.455); r0 = rho0(c, T, 0.47); Fl = F[T] + dF(c, RI, r0); mu = Fl + np.polyval(c, r0) * M / r0 * BARCM3
    b = 1 / (R * T); rn = r0 / M; rv = rn * math.exp(mu * b)
    for _ in range(100): rv = rn * math.exp(mu * b - 2 * B2[T] * rv)
    p = rv * (1 + B2[T] * rv) * RBAR * T / 1.01325; out.append((T, p))
    print(f"{T:6.1f} {r0:.4f} {Fl:+.4f} {mu:+.4f} {p:9.4f}   (leg densities {used})")
k = np.polyfit([1 / T for T, _ in out], np.log([p for _, p in out]), 1)
print('ln p against 1/T (straight line) crosses one atmosphere at %.1f K; heat of vaporisation %.2f kcal/mol' % (-k[0] / k[1], -k[0] * R))
