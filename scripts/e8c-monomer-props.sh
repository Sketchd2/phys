#!/usr/bin/env bash
# E8c, the monomer-derived terms of the hybrid law, first stage (PLAY.md E8c,
# the owner 2026-10-11): the MP2 dipole and polarisability tensor of water at
# a set of geometries drawn from the flexible pair data, on the CPU, while the
# GPU is busy. Resumable; the result is runs/e8c/monomer/props.txt, one line
# per geometry: id r1 r2 theta (bohr, bohr, rad) | dipole HF, MP2 (au) | MP2
# polarisability mean and diagonal (bohr^3).
#
#   scripts/e8c-monomer-props.sh [N_GEOMETRIES] [PARALLEL]
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"
ROOT=$HERE/runs/e8c; D=$ROOT/monomer; mkdir -p "$D"
N=${1:-60}; PAR=${2:-2}; BIN=$ROOT/bin
LOCK=$D/queue.lock
if [ -f "$LOCK" ] && tasklist //FI "PID eq $(cat "$LOCK")" 2>/dev/null | grep -q "$(cat "$LOCK")"; then echo "already running"; exit 0; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
[ -f "$D/geometries.txt" ] || python - "$N" "$ROOT" > "$D/geometries.txt" <<'PY'
import sys, math, random
n, root = int(sys.argv[1]), sys.argv[2]
A = 0.529177210903
mols = []
for f in ("flex/pairs-water-gpu-flex-mp2.txt", "flex/pairs-water-gpu-bias-flex-mp2.txt"):
    for l in open(f"{root}/{f}"):
        if l.startswith("#") or not l.strip(): continue
        v = [float(x) for x in l.split("|")[1].split()]
        for m in (0, 1):
            p = [(v[4*(3*m+i)], v[4*(3*m+i)+1], v[4*(3*m+i)+2]) for i in range(3)]
            d1 = [p[1][k]-p[0][k] for k in range(3)]; d2 = [p[2][k]-p[0][k] for k in range(3)]
            r1 = math.sqrt(sum(x*x for x in d1)); r2 = math.sqrt(sum(x*x for x in d2))
            th = math.acos(sum(a*b for a, b in zip(d1, d2)) / (r1*r2))
            mols.append((r1, r2, th))
random.seed(20261011)
# the averaged and the relaxed shapes first, then a spread over the draw, r1 >= r2 (the molecule is symmetric under exchange)
chosen = [(0.97894/A, 0.97894/A, math.radians(104.078)), (0.96154/A, 0.96154/A, math.radians(104.145))]
pool = [(max(a, b), min(a, b), t) for a, b, t in mols]
chosen += random.sample(pool, n - len(chosen))
for i, (r1, r2, th) in enumerate(chosen):
    print(i, repr(r1), repr(r2), repr(th))
PY
one() {
  id=$1; r1=$2; r2=$3; th=$4
  g=$D/g$id; mkdir -p "$g"; cd "$g"
  grep -q "^$id " "$D/props.txt" 2>/dev/null && exit 0
  python - "$HERE/grow-water.state" "$r1" "$r2" "$th" > grow-water.state <<'PY'
import sys, math
src, r1, r2, th = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
pos = [(0.0, 0.0, 0.0), (r1, 0.0, 0.0), (r2*math.cos(th), r2*math.sin(th), 0.0)]
k = 0
for line in open(src):
    if line.startswith("pos "):
        x, y, z = pos[k]; print("pos %r %r %r" % (x, y, z)); k += 1
    else: print(line, end="")
PY
  # phys-polar's last step, the fit of atomic polarisabilities by symmetry class, does not
  # end for a molecule whose two hydrogens differ (it ran 40 minutes on r1 2.01, r2 1.70
  # bohr); the tensor this reads is printed before it, so the step is cut off, not waited for.
  # The same for phys-esp, whose dipole is its first line.
  timeout 420 "$BIN/phys-polar.exe" water > polar.log 2>&1
  timeout 420 "$BIN/phys-esp.exe" water --gauss --cloud --mp2-dipole > esp.log 2>&1
  hf=$(grep -o "Hartree-Fock from the density [0-9.]*" esp.log | grep -o "[0-9.]*$" | head -1)
  mp=$(grep -o "MP2 from the field [0-9.]*" esp.log | grep -o "[0-9.]*$" | head -1)
  pm=$(grep "MP2: mean" polar.log | sed 's/.*mean \([0-9.]*\) bohr^3, diagonal \([0-9. -]*\),.*/\1 \2/')
  [ -n "$hf" ] && [ -n "$mp" ] && [ -n "$pm" ] || exit 1
  echo "$id $r1 $r2 $th | $hf $mp | $pm" >> "$D/props.txt"
}
export -f one; export D BIN HERE
stuck=0
while true; do
  have=$(wc -l < "$D/props.txt" 2>/dev/null || echo 0)
  todo=$(while read -r id r1 r2 th; do grep -q "^$id " "$D/props.txt" 2>/dev/null || echo "$id $r1 $r2 $th"; done < "$D/geometries.txt")
  [ -z "$todo" ] && break
  echo "$todo" | xargs -P "$PAR" -L1 bash -c 'one $0 $1 $2 $3' || sleep 60
  if [ "$(wc -l < "$D/props.txt" 2>/dev/null || echo 0)" = "$have" ]; then stuck=$((stuck+1)); [ $stuck -ge 3 ] && { echo "no progress" > "$D/FAILED"; exit 1; }; else stuck=0; fi
done
touch "$D/.complete"
