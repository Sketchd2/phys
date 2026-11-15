#!/bin/bash
# D1 for methane and ammonia (scale 1): water's two higher-angular-momentum steps (A, then A+E), the
# exponents scaled from oxygen's by the heavy atom's (C 0.70, N 0.85 of O's). Diffuse functions are left
# out: on water they barely moved the bracket.
B=../../../target/release/phys-shortfall
CA="6:4:1.4,1:3:1.6"; CE="6:4:0.55,1:3:0.6,6:3:0.18,1:2:0.2"
NA="7:4:1.7,1:3:1.6"; NE="7:4:0.68,1:3:0.6,7:3:0.21,1:2:0.2"
for v in "methane A $CA" "methane AE $CA,$CE" "ammonia A $NA" "ammonia AE $NA,$NE"; do set -- $v
  [ -s $1-more-$2.txt ] && grep -q '^ *[0-9]' $1-more-$2.txt || $B $1 1.0 --more $3 > $1-more-$2.txt 2>&1 || exit 1; done
echo finished
