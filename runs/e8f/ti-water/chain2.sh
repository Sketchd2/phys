#!/bin/bash
# The low-temperature end of the isochore again, cooling from the equilibrated 550 K liquid and
# twice as long (24 ps): the first pass went up from a 300 K liquid that had had 1 ps since its
# lattice and was still drifting (U/N -8.73 -> -8.98 over its 8 sampled ps).
export PHYS_EWALD=4.2 PHYS_EASE=0.15
B=$(pwd)/../../../target/release/phys-bulk
prev=$(pwd)/s-550-1.0/bulk-water-550.snap
mkdir -p cool
for T in 490 440 400 373 350 320 300; do
  d=cool/s-$T-1.0; mkdir -p $d
  if [ ! -f $d/done ]; then
    ( cd $d && cp ../../law.txt ../../grow-water.state . && PHYS_START_SNAP=$prev $B water law.txt $T 64 1.0 24 > liquid.log 2>&1 && touch done )
  fi
  prev=$(pwd)/$d/bulk-water-$T.snap
done
echo finished2
