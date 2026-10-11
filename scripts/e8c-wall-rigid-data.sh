#!/usr/bin/env bash
# E8c, the repulsive wall for the RIGID law (option 1): the same close-contact MP2 pairs as
# e8c-wall-data.sh, but with every molecule at the averaged shape the rigid law was
# built on (runs/e8c/avg/grow-water.state), no angle drawn. Waits until no other GPU
# program is running; resumable.
#
#   scripts/e8c-wall-rigid-data.sh
#
# runs/e8c/wall-rigid/pairs-water-gpu-wall-mp2.txt
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; D=$ROOT/wall-rigid; mkdir -p "$D" "$ROOT/gpubin-wall-rigid"
LOCK=$D/queue.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "already running"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
COUNT=${COUNT:-150}; LO=${WALL_LO:-1.0}; HI=${WALL_HI:-2.3}; T=${T:-1500}; CHUNK=${CHUNK:-5}; MAX_TRIES=${MAX_TRIES:-5}
qlog() { echo "[$(date -u '+%m-%d %H:%M')] $*" >> "$D/queue.log"; }
[ -f "$ROOT/gpubin-wall-rigid/phys-pairs-gpu.exe" ] || cp "$HERE/target/dev5/release/phys-pairs-gpu.exe" "$ROOT/gpubin-wall-rigid/"
[ -f "$D/grow-water.state" ] || cp "$ROOT/avg/grow-water.state" "$D/"
cp "$HERE/surface-water.txt" "$D/surface-water.txt"
P=$ROOT/gpubin-wall-rigid/phys-pairs-gpu.exe
cd "$D"
FILE=pairs-water-gpu-wall-mp2.txt
[ -f .complete ] && { echo done; exit 0; }
busy() { tasklist 2>/dev/null | grep -qiE "phys-(pairs|es|cluster|recheck|flexbulk)-gpu"; }
done_n() { [ -f "$FILE" ] && grep -vc '^#' "$FILE" || echo 0; }
fails=0
while [ "$(done_n)" -lt "$COUNT" ]; do
  while busy; do sleep 60; done
  have=$(done_n); target=$((have + CHUNK)); [ $target -gt $COUNT ] && target=$COUNT
  qlog "$have/$COUNT, running to $target"
  "$P" water $target --mp2 --wall $LO $HI >> pairs.log 2>&1
  if [ "$(done_n)" -le "$have" ]; then fails=$((fails+1)); qlog "no progress ($fails/$MAX_TRIES)"; [ $fails -ge $MAX_TRIES ] && { echo failed > FAILED; qlog "giving up"; exit 1; }; sleep 60; else fails=0; fi
done
touch .complete; qlog complete
