#!/bin/bash
# The exact first-order electrostatic energy of each pair already computed (phys-es-gpu), for the cloud refit. Each
# pairs file gets its OWN output directory, because the energies are keyed by pair index and the uniform and wall files
# both number from 0 (the first version wrote both to es-<name>.txt and silently skipped the wall file).
#   uniform pairs -> ./es-<name>.txt        wall pairs -> ./wall/es-<name>.txt
cd /e/Projects/phys/runs/e8e/pairs
log() { echo "[$(date '+%m-%d %H:%M')] $*" >> es.log; }
for m in "$@"; do
  [ -f grow-$m.state ] || continue
  log "$m: uniform"; ./phys-es-gpu.exe $m pairs-$m-gpu-mp2.txt --max-oo 9 >> es-run.log 2>&1
  if [ -f pairs-$m-gpu-wall-mp2.txt ]; then
    mkdir -p wall; cp -n grow-$m.state wall/; cp -f pairs-$m-gpu-wall-mp2.txt wall/
    log "$m: wall"; (cd wall && ../phys-es-gpu.exe $m pairs-$m-gpu-wall-mp2.txt --max-oo 9 >> ../es-run-wall.log 2>&1)
  fi
done
touch es2-done
