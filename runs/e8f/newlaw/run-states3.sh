#!/bin/bash
# run-states2.sh DIR LIST NMOL : run every state of LIST (lines "T density ps") with DIR/law.txt,
# NMOL molecules, from a lattice, plain cut-off (no Ewald), each in DIR/s-T-density.
# A finished state holds a file "done" and is skipped, so this can be run again after a stop.
# usage, from runs/e8f/newlaw:  ./run-states2.sh methane methane/part1.txt 108
cd "$(dirname "$0")"; DIR=$1; LIST=$(pwd)/$2; NMOL=$3; ROOT=$(pwd)/../../..
NAME=${NAME:-$(basename $DIR)}
while read T D PS; do
  [ -z "$T" ] && continue
  d=$DIR/s-$T-$D; mkdir -p $d
  if [ ! -f $d/done ]; then
    ( cd $d && cp ../law.txt $ROOT/grow-$NAME.state . && $ROOT/target/release/phys-bulk $NAME law.txt $T $NMOL $D $PS > liquid.log 2>&1 && touch done )
  fi
done < $LIST
echo "finished $2"
