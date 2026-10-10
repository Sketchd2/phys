#!/usr/bin/env bash
# E8c option 2, from the bend pair data to its liquid: fit the cloud's
# dependence on the angle, take it off the energies and fit the rest, assemble the
# law, run the liquid on the GPU. Waits for bend/.stages-complete; resumable (a
# step whose output exists is not repeated).
#
#   scripts/e8c-bend-fit.sh
#
# In runs/e8c/bend: cloud-geo.txt (phys-geofit), law-fit.txt (phys-fit),
# law-geo.txt (the two together), liquid.log and the files phys-flexbulk writes.
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; D=$ROOT/bend
LOCK=$D/fit.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "already running"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
mkdir -p "$ROOT/bin-geo"
for b in phys-geofit phys-fit phys-flexbulk phys-diag; do [ -f "$ROOT/bin-geo/$b.exe" ] || cp "$HERE/target/dev5/release/$b.exe" "$ROOT/bin-geo/" 2>/dev/null; done
[ -f "$ROOT/bin-geo/phys-flexbulk-gpu.exe" ] || cp "$HERE/target/dev5/release/phys-flexbulk-gpu.exe" "$ROOT/bin-geo/"
B=$ROOT/bin-geo
cd "$D"
qlog() { echo "[$(date -u '+%m-%d %H:%M')] $*" >> fit.log; }
while [ ! -f .stages-complete ]; do sleep 60; done
cp "$HERE/surface-water.txt" . 2>/dev/null
cp "$ROOT/avg/polar-water.txt" . 2>/dev/null
cp "$ROOT/avg/esp-water-cloud-mp2-fit.txt" cloud-avg.txt
cp "$ROOT/monomer-bend/props.txt" dipoles.txt
step() { name=$1; shift; [ -f ".done-fit-$name" ] && return 0; qlog "$name: start"; if "$@"; then touch ".done-fit-$name"; qlog "$name: done"; else qlog "$name: FAILED"; echo "$name failed" > FIT_FAILED; exit 1; fi; }

# 1. the cloud as quadratics in the angle, fitted to the exact electrostatics, the dipole held to MP2's
step geofit bash -c "\"$B/phys-geofit.exe\" water cloud-avg.txt es-water.txt -gpu-bias-bend-mp2 --dipoles dipoles.txt --theta0 104.078 --out cloud-geo.txt > geofit.log 2>&1"

# 2. what the angles change taken off every energy, slope and turn; the rest fitted
step fit bash -c "\"$B/phys-fit.exe\" water 1000 -gpu-bend-mp2,-gpu-bias-bend-mp2 --bisector --alpha polar-water.txt --esp cloud-geo.txt --charge-sites cloud-geo.txt --overlap --no-c8 --slopes deriv-water.txt --slope-pairs -gpu-bias-bend-mp2 --rotations rotate-water.txt --slope-weight 3 --geo-correct cloud-geo.txt > fit.out 2>&1"
cp law-water-gpu-bend-mp2+-gpu-bias-bend-mp2-pbe-bis-esp-ind.txt law-fit.txt
# 3. the law the liquid reads: the fit's, with the cloud's geo lines
[ -f law-geo.txt ] || { cp law-fit.txt law-geo.txt; grep -E "^geo" cloud-geo.txt >> law-geo.txt; }

# 4. the option 2 liquid on the GPU: bonds held at the average, 100 ps
step liquid bash -c "\"$B/phys-flexbulk-gpu.exe\" law-geo.txt 298 216 1.0 100 --bond 0.97894 --tag flex2 > liquid.log 2>&1"
touch .fit-complete; qlog "option 2 complete"
