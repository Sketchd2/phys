#!/usr/bin/env bash
# E8c option 2, the stages that follow the bend pair data (the same ones
# water-shape.sh runs on the rigid shape): exact electrostatic energies of the
# two molecules' densities, their slopes along the line between the centres,
# and the turns - all of the biased bend pairs, in runs/e8c/bend. Waits for
# bend/.complete; resumable and retrying like the other queue scripts.
#
#   scripts/e8c-bend-stages.sh
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; D=$ROOT/bend
LOCK=$D/stages.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "already running"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
mkdir -p "$ROOT/gpubin-stages"
for b in phys-es-gpu; do [ -f "$ROOT/gpubin-stages/$b.exe" ] || cp "$HERE/target/dev5/release/$b.exe" "$ROOT/gpubin-stages/"; done
G=$ROOT/gpubin-stages
cd "$D"
qlog() { echo "[$(date -u '+%m-%d %H:%M')] $*" >> stages.log; }
while [ ! -f .complete ]; do sleep 60; done
P=pairs-water-gpu-bias-bend-mp2.txt
fails=0
run_stage() { # name outfile command...
  name=$1; out=$2; shift 2
  while [ ! -f ".done-$name" ]; do
    qlog "$name: start"
    if "$@" > "$name.log" 2>&1; then touch ".done-$name"; qlog "$name: done"; else fails=$((fails+1)); qlog "$name: stopped ($fails)"; [ $fails -ge 8 ] && { echo "$name failed" > STAGES_FAILED; exit 1; }; sleep 60; fi
  done
}
run_stage electrostatics es-water.txt "$G/phys-es-gpu.exe" water $P --max-oo 4.4 --every 3
run_stage slopes deriv-water.txt "$G/phys-es-gpu.exe" water $P --deriv --max-oo 4.6 --every 3
run_stage turns rotate-water.txt "$G/phys-es-gpu.exe" water $P --deriv --rotate --max-oo 4.0 --every 4
touch .stages-complete; qlog "all stages complete"
