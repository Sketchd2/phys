#!/bin/bash
# E8e S1: the exact first-order electrostatic energy of each pair already computed (phys-es-gpu), for the
# cloud refit (phys-cloudfit). Resumable (it skips what es-<name>.txt holds); repeated until the pair queue is done.
cd /e/Projects/phys/runs/e8e/pairs
log() { echo "[$(date '+%m-%d %H:%M')] $*" >> es.log; }
while true; do
  for m in ammonia methanol ethanol methane; do
    for f in pairs-$m-gpu-mp2.txt pairs-$m-gpu-wall-mp2.txt; do
      [ -f "$f" ] || continue
      [ -f grow-$m.state ] || continue
      log "$m: $f"
      ./phys-es-gpu.exe $m $f --max-oo 9 >> es-run.log 2>&1
    done
  done
  [ -f all-done ] && break
  sleep 1800
done
touch es-done
