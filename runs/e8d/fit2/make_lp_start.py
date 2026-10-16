"""A start law for the lone-pair repulsion sites: the Born-Mayer law with one more site type (5), neutral, with
non-zero overlap strengths so the fit can move them (the fit works in logarithms, so a start of zero is stuck)."""
import sys, glob
sigma_lp = sys.argv[1]
src = glob.glob(r'E:\Projects\phys\runs\e8d\fit\bm_s30w0.01\law-water-*-pbe-bis-esp-ind.txt')[0]
lines = open(src).read().splitlines()
out = []
types = 5
for l in lines:
    w = l.split()
    if not w:
        continue
    if w[0] == 'charges':
        out.append(l + ' 0e0')
    elif w[0] == 'sigma':
        out.append(l + ' ' + sigma_lp)
    else:
        out.append(l)
# pair rows involving the new type: no atom-type repulsion or dispersion
for a in range(types + 1):
    out.append('pair %d 5 1e-300 2e0 1e-300 1e-300' % a if a < 5 else 'pair 5 5 1e-300 2e0 1e-300 1e-300')
# overlap strengths for the new type with every cloud type, non-zero starts
for a in (2, 3, 4, 5):
    out.append('overlap %d 5 0.05' % a)
open(r'E:\Projects\phys\runs\e8d\fit2\law0-lp-%s.txt' % sigma_lp, 'w').write('\n'.join(out) + '\n')
print('written law0-lp-%s.txt' % sigma_lp, len(out), 'lines')
