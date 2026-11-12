#!/bin/bash
# Water's free-energy path with the stored hole-free law (runs/e8d/core/law-core-B5.txt):
# up the 1.0 g/cm3 isochore from 300 K to 900 K, then down the 900 K isotherm to the
# dilute gas, each state started from the one before (64 molecules, Ewald, 12 ps).
export PHYS_EWALD=4.2 PHYS_EASE=0.15
B=$(pwd)/../../../target/release/phys-bulk
prev=$(pwd)/start.snap
run() {  # T density
  d=s-$1-$2; mkdir -p $d
  if [ ! -f $d/done ]; then
    ( cd $d && cp ../law.txt ../grow-water.state . && PHYS_START_SNAP=$prev $B water law.txt $1 64 $2 12 > liquid.log 2>&1 && touch done )
  fi
  prev=$(pwd)/$d/bulk-water-$1.snap
}
for T in 300 320 350 373 400 440 490 550 630 740 900; do run $T 1.0; done
for D in 0.88 0.76 0.64 0.52 0.4 0.3 0.2 0.12 0.06 0.02; do run 900 $D; done
echo finished
