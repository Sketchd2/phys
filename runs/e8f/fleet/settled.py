#!/usr/bin/env python3
"""Has a basis growth settled?  usage: settled.py grow-NAME.state [tol] [band]
The owner's rule: stop when the values a round reports have not changed for 5 rounds, or have
oscillated inside a band for 10. A round's values are its `hist` line (the molecule's relaxed
bond lengths and angles, bohr and radians). "Not changed" is every value within `tol` of the
round before, five rounds running; "oscillating" is every value's whole range over the last
ten rounds within `band`. Prints SETTLED (exit 0) or GROWING (exit 1) and the measure."""
import sys
tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3
band = float(sys.argv[3]) if len(sys.argv) > 3 else 5e-3
h = [[float(x) for x in l.split()[1:]] for l in open(sys.argv[1]) if l.startswith('hist')]
h = [r for r in h if len(r) == len(h[-1])] if h else h
step = lambda a, b: max(abs(x - y) for x, y in zip(a, b)) if a else 0.0
last5 = max((step(h[i], h[i - 1]) for i in range(len(h) - 5, len(h))), default=9.0) if len(h) >= 6 else 9.0
rng10 = max(max(c) - min(c) for c in zip(*h[-10:])) if len(h) >= 10 and h[-1] else 9.0
ok = last5 < tol or rng10 < band
print(('SETTLED' if ok else 'GROWING') + f' rounds {len(h)}; largest change over the last 5 rounds {last5:.2e} (limit {tol:g}); range over the last 10 {rng10:.2e} (limit {band:g})')
sys.exit(0 if ok else 1)
