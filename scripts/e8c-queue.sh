#!/usr/bin/env bash
# The E8c data queue: every shape water-shape.sh is run at, one after another,
# surviving a lost session, a crashed step or a reboot-and-relaunch.
#
#   scripts/e8c-queue.sh LANE       (LANE a or b; idempotent: a second copy of a lane exits at once)
#
# Everything lives in runs/e8c/ inside the repository, outside any session's
# scratch directory: private copies of the binaries (so a rebuild cannot hit a
# running .exe), one directory a shape, and queue-$LANE.log / status.txt. A shape is
# finished when its directory holds `.complete`. A step that fails is retried
# (the shape script resumes from its files); a shape that fails MAX_TRIES times
# in a row without gaining a file stops the queue and says so in FAILED.
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; mkdir -p "$ROOT/bin"
LANE=${1:?lane a or b}
LOCK=$ROOT/queue-$LANE.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "queue already running (pid $(cat "$LOCK"))"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
qlog() { echo "[$(date -u '+%m-%d %H:%M')] $*" | tee -a "$ROOT/queue-$LANE.log"; }
MAX_TRIES=${MAX_TRIES:-5}

# private binaries, taken once
for b in phys-polar phys-esp phys-fit phys-diag phys-bulk; do [ -f "$ROOT/bin/$b.exe" ] || cp "$HERE/target/dev2/release/$b.exe" "$ROOT/bin/"; done
mkdir -p "$ROOT/gpubin"; for b in phys-pairs-gpu phys-es-gpu; do [ -f "$ROOT/gpubin/$b.exe" ] || cp "$HERE/target/dev4/release/$b.exe" "$ROOT/gpubin/"; done
export PHYS_BIN=$ROOT/bin PHYS_GPU_BIN=$ROOT/gpubin
[ -f "$ROOT/grow-water.state" ] || cp "$HERE/grow-water.state" "$ROOT/grow-water.state"
[ -f "$ROOT/law-bias.txt" ] || cp "$HERE/law-water-mp2-cloud-slopes.txt" "$ROOT/law-bias.txt"

# name  r_OH(A)  theta(deg)  skip-liquid.
# Lane a is the plan's option (1) and then the stretches; lane b is the bend
# (2). Two lanes share the one GPU, so each runs slower than alone.
case $LANE in
a) SHAPES="
avg      0.97894 104.078 0
";;
b) SHAPES="
bend-lo  0.97894  98.5   1
bend-hi  0.97894 109.5   1
";;
esac
echo "$SHAPES" | while read -r name r th skip; do
  [ -z "${name:-}" ] && continue
  d=$ROOT/$name
  [ -f "$d/.complete" ] && continue
  tries=0
  while [ ! -f "$d/.complete" ]; do
    before=$(ls "$d" 2>/dev/null | wc -l); before_sz=$(du -sk "$d" 2>/dev/null | cut -f1)
    qlog "$name: start (try $((tries+1)))"; echo "$name running" > "$ROOT/status-$LANE.txt"
    # a shape's stages write their own pipeline.log; the queue only decides retry
    SKIP_LIQUID=$skip bash scripts/water-shape.sh "$d" "$r" "$th" "$ROOT/law-bias.txt" >> "$ROOT/$name.out" 2>&1 \
      && { touch "$d/.complete"; qlog "$name: complete"; break; }
    after=$(ls "$d" | wc -l); after_sz=$(du -sk "$d" | cut -f1)
    if [ "$after" = "$before" ] && [ "$after_sz" = "$before_sz" ]; then tries=$((tries+1)); else tries=0; fi
    qlog "$name: stopped (exit nonzero), no-progress tries $tries/$MAX_TRIES"
    if [ $tries -ge $MAX_TRIES ]; then qlog "$name: giving up"; echo "$name failed; see $ROOT/$name.out" > "$ROOT/FAILED-$LANE"; echo "$name FAILED" > "$ROOT/status-$LANE.txt"; exit 1; fi
    sleep 60
  done
done
rc=$?; [ -f "$ROOT/FAILED-$LANE" ] || { qlog "queue finished"; echo "all shapes complete" > "$ROOT/status-$LANE.txt"; }
