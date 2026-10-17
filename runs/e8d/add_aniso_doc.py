p='docs/PLAY.md'
s=open(p,encoding='utf-8',newline='').read()
k="     - **A second molecule: methane.**"
add="""     - **A direction-dependent repulsion was tested before being built, and the
       data does not support it** (the owner's idea: an antenna-style strength map
       on the repulsion; tests allowed, a law needs the owner's say). On 47
       liquid-sampled pairs (O-O 2.7-4.3 A) the hole-free law's energy error is
       0.18 kcal/mol rms against 2.12 of pair energy, and correlates 0.00 with the
       O-H..O angle and 0.02 with the side the hydrogen approaches (`aniso_test.py`).
       The slope error that drives the shells, over 32 pairs at 3.0-3.9 A, is
       +0.17 (bent, 90-115 deg, n 7), +0.11 (115-135, n 10) and +0.12 (135-180,
       n 15) kcal/mol per bohr with spreads of 0.13-0.17: a nearly constant radial
       offset, not an angular pattern (correlation with angle -0.20). Nothing was
       fitted or added. What is left as candidates for the missing shell gap:
       induction strength, liquid effects beyond three bodies, and nuclear quantum
       effects (being investigated).
"""
assert k in s
s=s.replace(k,add+k,1)
open(p,'w',encoding='utf-8',newline='').write(s)
n='runs/many/E8d-notes.md'
t=open(n,encoding='utf-8',newline='').read()
t+="""

## 14. A direction-dependent repulsion: tested, not supported

Idea (the owner's): give each atom's repulsion an antenna-style strength map, a size that varies with direction
(an anisotropic atom-atom potential), as the cure for the shared 3.1-3.8 A over-attraction on bent H..O contacts.
Test before any build, on the hole-free law (`runs/e8d/aniso_test.py`, `phys-slopeparts` with `VERBOSE=1`):

- Energies, 47 liquid-sampled pairs with O-O 2.7-4.3 A: rms error 0.18 kcal/mol against 2.12 of pair energy;
  correlation with the O-H..O angle 0.00, with the approach side to the acceptor's bisector 0.02, with H..O
  distance -0.20. At O-O 3.2-3.9 A the mean error is +0.085 +- 0.144 (90-115 deg, n 7), +0.018 +- 0.120
  (115-135, n 9), -0.052 +- 0.056 (135-180, n 4).
- Slopes along the line of centres, 32 pairs at 3.0-3.9 A: error +0.171 +- 0.131 (90-115 deg, n 7),
  +0.105 +- 0.154 (115-135, n 10), +0.119 +- 0.169 (135-180, n 15) kcal/mol per bohr; correlation with angle -0.20.

The error is a nearly constant radial offset, not an angular pattern, so a direction-dependent term fitted here would
mostly act as a radial adjustment, which the Born-Mayer and refit searches already exhausted. Not built.
Candidates left: induction strength, many-body effects beyond three bodies, nuclear quantum effects.
"""
open(n,'w',encoding='utf-8',newline='').write(t)
