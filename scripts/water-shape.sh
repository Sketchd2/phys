#!/usr/bin/env bash
# One rigid shape of water, from its monomer to its liquid (PLAY.md E8c).
#
#   scripts/water-shape.sh DIR R_OH_ANGSTROM THETA_DEGREES [BIAS_LAW]
#
# In DIR (made if absent): the shape's state file, its polarisabilities, its
# cloud electrostatics fitted to the exact first-order energies, MP2 pairs (a
# uniform draw and one biased by BIAS_LAW), the slopes and turns of their
# energies, the fit, and a 100 ps liquid. Every step is resumable: a file
# already complete is not made again, an interrupted one carries on.
# Sizes are the leaner ones the comparison uses, the same for each shape.
set -u
DIR=$1; R=$2; THETA=$3; BIAS=${4:-$(pwd)/law-water-mp2-cloud-slopes.txt}
HERE=$(pwd)
BIN=${PHYS_BIN:-$HERE/target/dev2/release}
GPUBIN=${PHYS_GPU_BIN:-$HERE/target/dev4/release}
RANDOM_PAIRS=${RANDOM_PAIRS:-60}
BIAS_PAIRS=${BIAS_PAIRS:-200}
mkdir -p "$DIR" && cd "$DIR" || exit 1
log() { echo "[$(date -u +%H:%M)] $*" | tee -a pipeline.log; }

# 1. the shape: the engine's growth state with the atoms moved
if [ ! -f grow-water.state ]; then
  python - "$HERE/grow-water.state" "$R" "$THETA" > grow-water.state <<'EOF'
import sys, math
src, r, th = sys.argv[1], float(sys.argv[2]) / 0.529177210903, math.radians(float(sys.argv[3]))
pos = [(0.0, 0.0, 0.0), (r, 0.0, 0.0), (r * math.cos(th), r * math.sin(th), 0.0)]
k = 0
for line in open(src):
    if line.startswith("pos "):
        x, y, z = pos[k]
        print("pos %r %r %r" % (x, y, z))
        k += 1
    else:
        print(line, end="")
EOF
fi
cp "$BIAS" law0.txt

# 2. what the molecule's electrons say of it
[ -f polar-water.txt ] || { log "polarisabilities"; "$BIN/phys-polar.exe" water > polar.log 2>&1; }
[ -f esp-water-cloud-mp2.txt ] || { log "cloud potential fit"; "$BIN/phys-esp.exe" water --gauss --cloud --mp2-dipole > esp.log 2>&1; }
MU=$(grep -o "MP2 from the field [0-9.]* au" esp.log | grep -o "[0-9.]*" | head -1)

# 3. pairs: uniform, then biased (the pair program keeps its own count)
log "uniform pairs"; "$GPUBIN/phys-pairs-gpu.exe" water $RANDOM_PAIRS --mp2 > pairs-uniform.log 2>&1
log "biased pairs"; "$GPUBIN/phys-pairs-gpu.exe" water $BIAS_PAIRS --mp2 --bias law0.txt 2.0 > pairs-bias.log 2>&1

# 4. exact electrostatics, slopes, turns
log "electrostatics"; "$GPUBIN/phys-es-gpu.exe" water pairs-water-gpu-bias-mp2.txt --max-oo 4.4 --every 3 > es.log 2>&1
log "cloud refit"; "$BIN/phys-diag.exe" water esp-water-cloud-mp2.txt --fit-es es-water.txt --dipole "$MU" -gpu-bias-mp2 > cloudfit.log 2>&1
log "slopes"; "$GPUBIN/phys-es-gpu.exe" water pairs-water-gpu-bias-mp2.txt --deriv --max-oo 4.6 --every 3 > slopes.log 2>&1
log "turns"; "$GPUBIN/phys-es-gpu.exe" water pairs-water-gpu-bias-mp2.txt --deriv --rotate --max-oo 4.0 --every 4 > turns.log 2>&1

# 5. the law and its liquid
log "fit"
"$BIN/phys-fit.exe" water 1000 -gpu-mp2,-gpu-bias-mp2 --bisector --alpha polar-water.txt \
  --esp esp-water-cloud-mp2-fit.txt --charge-sites esp-water-cloud-mp2-fit.txt --overlap --no-c8 \
  --slopes deriv-water.txt --slope-pairs -gpu-bias-mp2 --rotations rotate-water.txt --slope-weight 3 > fit.log 2>&1
cp law-water-gpu-mp2+-gpu-bias-mp2-pbe-bis-esp-ind.txt law-final.txt
log "liquid"; "$BIN/phys-bulk.exe" water law-final.txt 298 216 1.0 100 > liquid.log 2>&1
log "done"
