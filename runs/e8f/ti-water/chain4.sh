#!/bin/bash
# Density legs at 400 and 440 K for the stored water law (64 molecules, Ewald, 12 ps each),
# each started from the state before it.
export PHYS_EWALD=4.2 PHYS_EASE=0.15
B=$(pwd)/../../../target/release/phys-bulk
for T in 400 440; do
  prev=$(pwd)/cool/s-$T-1.0/bulk-water-$T.snap
  for D in 0.96 0.92 0.88 0.84 0.8 0.76; do
    d=cool/s-$T-$D; mkdir -p $d
    if [ ! -f $d/done ]; then
      ( cd $d && cp ../../law.txt ../../../../../grow-water.state . && PHYS_START_SNAP=$prev $B water law.txt $T 64 $D 12 > liquid.log 2>&1 && touch done )
    fi
    prev=$(pwd)/$d/bulk-water-$T.snap
  done
done
echo finished4
