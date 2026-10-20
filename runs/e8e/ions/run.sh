#!/bin/bash
# E8e S5 data that needs nothing unsettled: MP2 pairs of an ion with a water molecule, uniform draw
# and close contact. Water's shape is the repo's grown one. Resumable.
cd /e/Projects/phys/runs/e8e/ions
log() { echo "[$(date '+%m-%d %H:%M')] $*" >> queue.log; }
try() { for i in 1 2 3 4 5; do "$@" >> run.log 2>&1 && return 0; log "failed (try $i): $*"; sleep 30; done; return 1; }
for ion in sodium chloride; do
  log "$ion: starting"
  try ./phys-pairs-unlike-gpu.exe water $ion 60
  try ./phys-pairs-unlike-gpu.exe water $ion 40 --wall 1.0 2.3
  log "$ion: done"
done
touch all-done
