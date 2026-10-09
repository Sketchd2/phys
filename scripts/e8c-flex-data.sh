#!/usr/bin/env bash
# E8c, the data for water that moves (options 2 and 3): pairs whose molecules
# are each drawn from the monomer's energy surface, so every internal
# coordinate varies, an asymmetric stretch included. The same pairs are used
# for both flexible treatments; the law that follows the geometry is fitted to
# them afterwards (PLAY.md E8c).
#
#   scripts/e8c-flex-data.sh      (idempotent and resumable, like e8c-queue.sh)
#
# runs/e8c/flex/: `pairs-water-gpu-flex-mp2.txt` (uniform draw, UNIFORM pairs) and
# `pairs-water-gpu-bias-flex-mp2.txt` (drawn where the relaxed law puts a liquid's
# neighbours, BIASED pairs). T is the sampling temperature of the shapes.
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; D=$ROOT/flex; mkdir -p "$D" "$ROOT/gpubin-flex"
LOCK=$D/queue.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "already running"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
UNIFORM=${UNIFORM:-100}; BIASED=${BIASED:-300}; T=${T:-1500}; MAX_TRIES=${MAX_TRIES:-5}
qlog() { echo "[$(date -u '+%m-%d %H:%M')] $*" | tee -a "$D/queue.log"; }
[ -f "$ROOT/gpubin-flex/phys-pairs-gpu.exe" ] || cp "$HERE/target/dev5/release/phys-pairs-gpu.exe" "$ROOT/gpubin-flex/"
for f in grow-water.state; do [ -f "$D/$f" ] || cp "$HERE/$f" "$D/"; done
cp "$HERE/surface-water.txt" "$D/surface-water.txt"
[ -f "$D/law0.txt" ] || cp "$ROOT/law-bias.txt" "$D/law0.txt"
P=$ROOT/gpubin-flex/phys-pairs-gpu.exe
cd "$D"
[ -f .complete ] && { echo done; exit 0; }
tries=0
while [ ! -f .complete ]; do
  before=$(cat pairs-water-gpu*-flex-mp2.txt 2>/dev/null | wc -l)
  qlog "start (try $((tries+1)))"
  { "$P" water $UNIFORM --mp2 --flex surface-water.txt $T >> uniform.log 2>&1 \
    && "$P" water $BIASED --mp2 --bias law0.txt 2.0 --flex surface-water.txt $T >> bias.log 2>&1; } \
    && { touch .complete; qlog "complete"; break; }
  after=$(cat pairs-water-gpu*-flex-mp2.txt 2>/dev/null | wc -l)
  if [ "$after" = "$before" ]; then tries=$((tries+1)); else tries=0; fi
  qlog "stopped; no-progress tries $tries/$MAX_TRIES"
  [ $tries -ge $MAX_TRIES ] && { qlog "giving up"; echo failed > FAILED; exit 1; }
  sleep 60
done
