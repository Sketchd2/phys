#!/bin/bash
# job.sh SPECIES MODE CUTOFF MAXSEP MAX : MAX close pairs of SPECIES (separation under MAXSEP bohr,
# always the same ones) through phys-overlap in MODE (ghost | frozen | full), with shell pairs
# below CUTOFF left out of the table. Writes r-SPECIES-MODE-CUTOFF.txt (the tool's rows) and
# appends "seconds pairs peak_MB" to the matching .time file. Safe to start again: rows already
# written are skipped. Run from runs/e8f/break.
S=$1; MODE=$2; CUT=$3; SEP=$4; MAX=$5; ROOT=$(cd ../../.. && pwd)
case $S in
  water) P=$ROOT/runs/e8c/avg/pairs-water-gpu-bias-mp2.txt;;
  ammonia) P=$ROOT/runs/e8e/pairs/pairs-ammonia-gpu-mp2.txt;;
  methane) P=$ROOT/runs/e8d/methane/pairs-methane-gpu-mp2.txt;;
  methanol) P=$ROOT/runs/e8e/pairs-methanol/pairs-methanol-gpu-mp2.txt;;
esac
case $MODE in ghost) F="";; frozen) F="--frozen";; full) F="--frozen --full";; esac
cp -n $ROOT/grow-$S.state . 2>/dev/null
export PHYS_PAIR_CUTOFF=$CUT PHYS_SPILL_DIR=/tmp/spill; mkdir -p /tmp/spill
O=r-$S-$MODE-$CUT.txt; T0=$(date +%s)
$ROOT/target/release/phys-overlap $S $P --max-sep $SEP --max $MAX $F --out $O > log-$S-$MODE-$CUT.txt 2>&1 &
PID=$!; PEAK=0
while kill -0 $PID 2>/dev/null; do   # the process's own high-water mark, sampled while it lives
  M=$(awk '/VmHWM/ {print $2}' /proc/$PID/status 2>/dev/null); [ -n "$M" ] && PEAK=$M; sleep 2
done
wait $PID || exit 1
echo "$(( $(date +%s) - T0 )) $(grep -c . $O) $(( PEAK / 1024 ))" >> r-$S-$MODE-$CUT.time
