#!/usr/bin/env bash
# E8c, more pair data for water that moves, drawn wider (the owner, 2026-10-11:
# "collect more first, and if the GPU is idle run the wider data too"): the same
# draw as e8c-flex-data.sh at a higher sampling temperature (T, default 2700 K,
# the O-H stretch's zero-point-equivalent), run ONLY WHEN THE GPU IS IDLE: before
# every chunk of CHUNK pairs it waits until no other phys-*-gpu program is
# running, so it fills the gaps (a liquid, a fit) and never slows another stage
# by more than one chunk. Idempotent and resumable.
#
#   scripts/e8c-flexwide-data.sh
#
# runs/e8c/flexwide/: pairs-water-gpu-flex-mp2.txt (UNIFORM) and
# pairs-water-gpu-bias-flex-mp2.txt (BIASED).
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; D=$ROOT/flexwide; mkdir -p "$D" "$ROOT/gpubin-flexwide"
LOCK=$D/queue.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "already running"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
UNIFORM=${UNIFORM:-100}; BIASED=${BIASED:-300}; T=${T:-2700}; CHUNK=${CHUNK:-5}; MAX_TRIES=${MAX_TRIES:-5}
qlog() { echo "[$(date -u '+%m-%d %H:%M')] $*" >> "$D/queue.log"; }
[ -f "$ROOT/gpubin-flexwide/phys-pairs-gpu.exe" ] || cp "$HERE/target/dev5/release/phys-pairs-gpu.exe" "$ROOT/gpubin-flexwide/"
[ -f "$D/grow-water.state" ] || cp "$HERE/grow-water.state" "$D/"
cp "$HERE/surface-water.txt" "$D/surface-water.txt"
[ -f "$D/law0.txt" ] || cp "$ROOT/law-bias.txt" "$D/law0.txt"
P=$ROOT/gpubin-flexwide/phys-pairs-gpu.exe
cd "$D"
[ -f .complete ] && { echo done; exit 0; }
busy() { tasklist 2>/dev/null | grep -qiE "phys-(pairs|es|cluster|recheck)-gpu"; }
done_n() { [ -f "$1" ] && grep -vc '^#' "$1" || echo 0; }
fails=0
for stage in uniform bias; do
  if [ $stage = uniform ]; then total=$UNIFORM; file=pairs-water-gpu-flex-mp2.txt; extra=""; else total=$BIASED; file=pairs-water-gpu-bias-flex-mp2.txt; extra="--bias law0.txt 2.0"; fi
  while [ "$(done_n $file)" -lt "$total" ]; do
    while busy; do sleep 60; done
    have=$(done_n $file); target=$((have + CHUNK)); [ $target -gt $total ] && target=$total
    qlog "$stage: $have/$total, running to $target"
    "$P" water $target --mp2 $extra --flex surface-water.txt $T >> $stage.log 2>&1
    if [ "$(done_n $file)" -le "$have" ]; then fails=$((fails+1)); qlog "no progress ($fails/$MAX_TRIES)"; [ $fails -ge $MAX_TRIES ] && { echo failed > FAILED; qlog "giving up"; exit 1; }; sleep 60; else fails=0; fi
  done
done
touch .complete; qlog complete
