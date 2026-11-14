#!/bin/bash
# grow.sh NAME "ATOMS;BONDS" : grow NAME's basis until settled.py says it has settled, at most
# MAXH hours (default 20). Run in runs/e8f/species/NAME; safe to start again (phys-grow resumes
# from grow-NAME.state). Exit 0 with the state file in place either way; the last line of
# grow.status says SETTLED, or STOPPED AT THE TIME LIMIT with where it had got to.
NAME=$1; SPEC=$2; ROOT=$(cd "$(dirname "$0")/../../.." && pwd); F=$ROOT/runs/e8f/fleet
export PHYS_SPILL_DIR=${PHYS_SPILL_DIR:-/tmp/spill}; mkdir -p $PHYS_SPILL_DIR
if [ -f grow-$NAME.state ] && python3 $F/settled.py grow-$NAME.state > grow.status; then exit 0; fi
if [ -n "$SPEC" ]; then $ROOT/target/release/phys-grow "$SPEC" $NAME >> grow.log 2>&1 & else $ROOT/target/release/phys-grow $NAME >> grow.log 2>&1 & fi
P=$!; END=$(( $(date +%s) + ${MAXH:-20} * 3600 ))
while kill -0 $P 2>/dev/null; do
  sleep 30
  if [ -f grow-$NAME.state ] && python3 $F/settled.py grow-$NAME.state > grow.status; then kill $P; wait $P 2>/dev/null; exit 0; fi
  if [ $(date +%s) -gt $END ]; then kill $P; wait $P 2>/dev/null; echo "STOPPED AT THE TIME LIMIT; $(cat grow.status)" > grow.status; exit 0; fi
done
python3 $F/settled.py grow-$NAME.state > grow.status; exit 0
