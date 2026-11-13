#!/bin/bash
# Density legs: at 1.0 g/cm3 this law's liquid is under 0.5 to 2.7 kbar between 300 and 400 K,
# too far for a first-order correction, so each temperature is walked down in density from its
# cooled 1.0 g/cm3 state until the pressure changes sign (64 molecules, Ewald, 12 ps each).
export PHYS_EWALD=4.2 PHYS_EASE=0.15
B=$(pwd)/../../../target/release/phys-bulk
for T in 373 350 400 320; do
  prev=$(pwd)/cool/s-$T-1.0/bulk-water-$T.snap
  for D in 0.96 0.92 0.88 0.84; do
    d=cool/s-$T-$D; mkdir -p $d
    if [ ! -f $d/done ]; then
      ( cd $d && cp ../../law.txt ../../grow-water.state . && PHYS_START_SNAP=$prev $B water law.txt $T 64 $D 12 > liquid.log 2>&1 && touch done )
    fi
    prev=$(pwd)/$d/bulk-water-$T.snap
  done
done
echo finished3
