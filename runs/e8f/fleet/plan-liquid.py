#!/usr/bin/env python3
"""Choose a species' liquid states from its own law (HANDOFF.md D3): no measured number goes in.
usage, from runs/e8f/species/NAME with law.txt in place:   plan-liquid.py NAME

The law's second virial coefficient is scanned in temperature. Two numbers come out of it:
the Boyle temperature T_B, where it changes sign, and its largest value B_max, the molecule's
size as the law has it. For the two laws already run to a boiling point, the liquid's molar
volume at the boiling point is 1.11 B_max (water) and 1.08 B_max (methane), and the boiling
point 0.29 T_B and 0.25 T_B; those set the ladder, not the answer.

  the isochore   at rho_ref = M / (1.05 B_max), a little denser than the liquid so that no state
                 on it is under tension, twelve temperatures from T_B down to 0.20 T_B
  the isotherm   at T_B, ten densities from rho_ref down to the dilute gas
  the legs       at 0.22, 0.25, 0.28 and 0.32 T_B, three lower densities each

Writes liquid/states.txt (one "T density ps" a line), liquid/b2.arg (the "T:B2,..." list
path.py takes, at the ladder's temperatures) and plan.ok."""
import sys, os, subprocess, math
name = sys.argv[1]
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')); B = ROOT + '/target/release'
MASS = {1: 1.008, 6: 12.011, 7: 14.007, 8: 15.999, 17: 35.45}
M = sum(MASS[int(z)] for z in open(f'grow-{name}.atoms').read().split())
def b2(temps, n=1500):
    out = subprocess.run([B + '/phys-b2', name, 'law.txt'] + [f'{t:g}' for t in temps] + ['--orientations', str(n)], capture_output=True, text=True).stdout
    return {float(l.split()[0]): float(l.split()[2]) for l in out.splitlines() if l.strip()}
grid = [100 * 1.3 ** i for i in range(21)]                 # 100 K to 19,000 K
scan = b2(grid); ts = sorted(scan)
cross = [(a, b) for a, b in zip(ts, ts[1:]) if scan[a] < 0 <= scan[b]]
if not cross: sys.exit('the second virial coefficient does not change sign between 100 and 19000 K')
a, b = cross[0]; tb = a + (b - a) * (0 - scan[a]) / (scan[b] - scan[a])       # the Boyle temperature
bmax = max(scan.values())
rho = M / (1.05 * bmax)
r = lambda x, s=3: float(f'{x:.{s}g}')
iso_t = [r(tb * f) for f in (1, 0.83, 0.7, 0.59, 0.5, 0.43, 0.37, 0.32, 0.28, 0.25, 0.22, 0.2)]
leg_t = [r(tb * f) for f in (0.22, 0.25, 0.28, 0.32)]
dipole = 0.0
for l in open('cloud.log'):
    if 'dipole' in l and ' au' in l:
        dipole = float(l.split('au')[0].split()[-1]); break
ewald = dipole > 0.05          # a polar molecule's charges need the long-range sum
ps = 12 if ewald else 20
states = [(t, r(rho)) for t in iso_t] + [(iso_t[0], r(rho * f)) for f in (0.88, 0.76, 0.64, 0.52, 0.4, 0.28, 0.18, 0.1, 0.05, 0.02)]
states += [(t, r(rho * f)) for t in leg_t for f in (0.95, 0.9, 0.85)]
os.makedirs('liquid', exist_ok=True)
open('liquid/states.txt', 'w').write(''.join(f'{t:g} {d:g} {ps}\n' for t, d in states))
vals = b2(sorted(set(iso_t)), 8000)
open('liquid/b2.arg', 'w').write(','.join(f'{t:g}:{vals[t]:.1f}' for t in sorted(vals)))
open('liquid/setup.txt', 'w').write(f"M {M:.3f}\nThot {iso_t[0]:g}\nrho_ref {r(rho):g}\newald {'4.2' if ewald else 'off'}\nmolecules {64 if ewald else 108}\n")
open('plan.ok', 'w').write(f"Boyle temperature {tb:.0f} K; B_max {bmax:.2f} cm3/mol; molar mass {M:.3f}; rho_ref {r(rho):g} g/cm3; dipole {dipole:.3f} au, {'Ewald, 64 molecules, 12 ps' if ewald else 'cut-off, 108 molecules, 20 ps'}; {len(states)} states\n")
print(open('plan.ok').read())
