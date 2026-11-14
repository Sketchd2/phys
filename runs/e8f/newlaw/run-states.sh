#!/bin/bash
# run-states.sh DIR LIST : run every state of LIST (lines "T density ps start-snapshot", the
# snapshot relative to runs/e8f) with DIR/law.txt, 64 molecules, Ewald, each in DIR/s-T-density.
# A finished state holds a file "done" and is skipped, so this can be run again after a stop.
# usage, from runs/e8f/newlaw:  ./run-states.sh water water/part1.txt
cd "$(dirname "$0")"; DIR=$1; LIST=$(pwd)/$2; ROOT=$(pwd)/../../..; E8F=$(pwd)/..
NAME=$(basename $DIR); NMOL=${NMOL:-64}
export PHYS_EWALD=${PHYS_EWALD:-4.2} PHYS_EASE=0.15
while read T D PS SNAP; do
  [ -z "$T" ] && continue
  d=$DIR/s-$T-$D; mkdir -p $d
  if [ ! -f $d/done ]; then
    ( cd $d && cp ../law.txt $ROOT/grow-$NAME.state . && PHYS_START_SNAP=$E8F/$SNAP $ROOT/target/release/phys-bulk $NAME law.txt $T $NMOL $D $PS > liquid.log 2>&1 && touch done )
  fi
done < $LIST
echo "finished $2"
