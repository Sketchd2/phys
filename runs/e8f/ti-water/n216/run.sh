#!/bin/sh
# Re-runnable: skips a state that has `done`.
cd "$(dirname "$0")"
ROOT=$(pwd)
BIN=$ROOT/../../../../target/release/phys-bulk
mkdir -p s-900-1.0 s-400-1.0
if [ ! -f s-900-1.0/done ]; then
  cp law.txt grow-water.state s-900-1.0/
  (cd s-900-1.0 && PHYS_EWALD=4.2 $BIN water law.txt 900 216 1.0 12 > liquid.log 2>&1 && touch done) || exit 1
fi
if [ ! -f s-400-1.0/done ]; then
  cp law.txt grow-water.state s-400-1.0/
  (cd s-400-1.0 && PHYS_EWALD=4.2 PHYS_EASE=0.15 PHYS_START_SNAP=$ROOT/s-900-1.0/bulk-water-900.snap $BIN water law.txt 400 216 1.0 24 > liquid.log 2>&1 && touch done) || exit 1
fi
