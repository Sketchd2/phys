#!/bin/bash
# liquid-state.sh NAME T RHO PS : one liquid state of NAME's law, from runs/e8f/species/NAME, in
# liquid/s-T-RHO. Lattice start; Ewald or cut-off and the molecule count as plan-liquid.py chose
# (liquid/setup.txt). Finished when the directory holds "done"; safe to start again.
NAME=$1; T=$2; RHO=$3; PS=$4; ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
D=liquid/s-$T-$RHO; mkdir -p $D; [ -f $D/done ] && exit 0
E=$(awk '$1=="ewald"{print $2}' liquid/setup.txt); N=$(awk '$1=="molecules"{print $2}' liquid/setup.txt)
[ "$E" != off ] && export PHYS_EWALD=$E
cp law.txt grow-$NAME.state grow-$NAME.atoms $D/ && cd $D && $ROOT/target/release/phys-bulk $NAME law.txt $T $N $RHO $PS > liquid.log 2>&1 && touch done
