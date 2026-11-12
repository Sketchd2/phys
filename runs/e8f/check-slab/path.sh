#!/bin/bash
# As ../ti-methane/path.sh, with the pair law cut off as phys-coexist cuts it
# (switched 12 -> 14 bohr on the centres), one state at a time.
export PHYS_RCUT=14
B=$(pwd)/../../../target/release/phys-bulk
while read T D; do
  d=s-$T-$D; mkdir -p $d; [ -f $d/done ] && continue
  ( cd $d && cp ../law.txt ../grow-methane.state . && $B methane law.txt $T 108 $D 20 > liquid.log 2>&1 && touch done )
done < states.txt
echo finished
