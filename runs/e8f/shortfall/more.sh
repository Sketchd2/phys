#!/bin/bash
# What added functions do to the water dimer's bracket (scale 1): higher angular momentum, diffuse functions, both.
B=../../../target/release/phys-shortfall
A="8:4:2.0,1:3:1.6"; D="8:0:0.06,8:1:0.05,8:2:0.12,1:0:0.03,1:1:0.08"; E="8:4:0.8,1:3:0.6,8:3:0.25,1:2:0.2"
for v in "A $A" "D $D" "AD $A,$D" "ADE $A,$D,$E"; do set -- $v; [ -s water-more-$1.txt ] || $B water 1.0 --more $2 > water-more-$1.txt 2>&1; done
echo finished
