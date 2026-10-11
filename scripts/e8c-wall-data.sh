#!/usr/bin/env bash
# E8c, the repulsive wall (PLAY.md E8c, scheduled by the owner 2026-10-11): MP2 pair
# energies at close contact - the nearest atoms of the two molecules WALL_LO to
# WALL_HI angstrom apart - which the uniform and biased draws never reach and
# which the laws' wall must be fitted to. Angles drawn from the monomer's surface
# at the held bond length (so the average angle and the bend are both covered).
# Waits until no other GPU program is running; resumable.
#
#   scripts/e8c-wall-data.sh
#
# runs/e8c/wall/pairs-water-gpu-wall-bend-mp2.txt
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; D=$ROOT/wall; mkdir -p "$D" "$ROOT/gpubin-wall"
LOCK=$D/queue.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "already running"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
COUNT=${COUNT:-150}; LO=${WALL_LO:-1.0}; HI=${WALL_HI:-2.3}; T=${T:-1500}; CHUNK=${CHUNK:-5}; MAX_TRIES=${MAX_TRIES:-5}
qlog() { echo "[$(date -u '+%m-%d %H:%M')] $*" >> "$D/queue.log"; }
[ -f "$ROOT/gpubin-wall/phys-pairs-gpu.exe" ] || cp "$HERE/target/dev5/release/phys-pairs-gpu.exe" "$ROOT/gpubin-wall/"
for f in grow-water.state; do [ -f "$D/$f" ] || cp "$HERE/$f" "$D/"; done
cp "$HERE/surface-water.txt" "$D/surface-water.txt"
P=$ROOT/gpubin-wall/phys-pairs-gpu.exe
cd "$D"
FILE=pairs-water-gpu-wall-bend-mp2.txt
[ -f .complete ] && { echo done; exit 0; }
busy() { tasklist 2>/dev/null | grep -qiE "phys-(pairs|es|cluster|recheck|flexbulk)-gpu"; }
done_n() { [ -f "$FILE" ] && grep -vc '^#' "$FILE" || echo 0; }
fails=0
while [ "$(done_n)" -lt "$COUNT" ]; do
  while busy; do sleep 60; done
  have=$(done_n); target=$((have + CHUNK)); [ $target -gt $COUNT ] && target=$COUNT
  qlog "$have/$COUNT, running to $target"
  "$P" water $target --mp2 --wall $LO $HI --flex surface-water.txt $T --bonds-at ${BOND:-0.97894} >> pairs.log 2>&1
  if [ "$(done_n)" -le "$have" ]; then fails=$((fails+1)); qlog "no progress ($fails/$MAX_TRIES)"; [ $fails -ge $MAX_TRIES ] && { echo failed > FAILED; qlog "giving up"; exit 1; }; sleep 60; else fails=0; fi
done
touch .complete; qlog complete
