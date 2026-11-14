#!/bin/bash
# scan.sh DIR CLOUD C6 POLAR tags [core] : the weighted fit of K at each dispersion damping exponent
for D in 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.4; do echo -n "damp $D  "; DAMP=$D python3 fitk4.py "$@" | sed -n 1,2p | tr '\n' ' ' | sed 's/deepest pair [-0-9.]* kcal\/mol; //; s/(unweighted[^)]*)//; s/ *weighted  */ /'; echo; done
