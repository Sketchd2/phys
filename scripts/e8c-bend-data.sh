#!/usr/bin/env bash
# E8c, the data for option 2 (bonds held, bend free): pairs whose molecules each
# have both O-H bonds at BOND angstrom (the ground-state average, 0.97894) and
# the H-O-H angle drawn from the monomer's surface at that bond length, so the
# law of option 2 is fitted to the one coordinate that moves in it.
#
#   scripts/e8c-bend-data.sh      (idempotent and resumable)
#
# runs/e8c/bend/: pairs-water-gpu-bend-mp2.txt (uniform) and
# pairs-water-gpu-bias-bend-mp2.txt (biased by the relaxed law).
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; D=$ROOT/bend; mkdir -p "$D" "$ROOT/gpubin-bend"
LOCK=$D/queue.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "already running"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
UNIFORM=${UNIFORM:-60}; BIASED=${BIASED:-250}; T=${T:-1500}; MAX_TRIES=${MAX_TRIES:-5}
qlog() { echo "[$(date -u '+%m-%d %H:%M')] $*" | tee -a "$D/queue.log"; }
[ -f "$ROOT/gpubin-bend/phys-pairs-gpu.exe" ] || cp "$HERE/target/dev5/release/phys-pairs-gpu.exe" "$ROOT/gpubin-bend/"
for f in grow-water.state; do [ -f "$D/$f" ] || cp "$HERE/$f" "$D/"; done
cp "$HERE/surface-water.txt" "$D/surface-water.txt"
[ -f "$D/law0.txt" ] || cp "$ROOT/law-bias.txt" "$D/law0.txt"
P=$ROOT/gpubin-bend/phys-pairs-gpu.exe
cd "$D"
[ -f .complete ] && { echo done; exit 0; }
tries=0
while [ ! -f .complete ]; do
  before=$(cat pairs-water-gpu*-bend-mp2.txt 2>/dev/null | wc -l)
  qlog "start (try $((tries+1)))"
  { "$P" water $UNIFORM --mp2 --flex surface-water.txt $T --bonds-at ${BOND:-0.97894} >> uniform.log 2>&1 \
    && "$P" water $BIASED --mp2 --bias law0.txt 2.0 --flex surface-water.txt $T --bonds-at ${BOND:-0.97894} >> bias.log 2>&1; } \
    && { touch .complete; qlog "complete"; break; }
  after=$(cat pairs-water-gpu*-bend-mp2.txt 2>/dev/null | wc -l)
  if [ "$after" = "$before" ]; then tries=$((tries+1)); else tries=0; fi
  qlog "stopped; no-progress tries $tries/$MAX_TRIES"
  [ $tries -ge $MAX_TRIES ] && { qlog "giving up"; echo failed > FAILED; exit 1; }
  sleep 60
done
