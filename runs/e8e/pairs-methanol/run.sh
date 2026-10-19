#!/bin/bash
# E8e S1: methanol's MP2 pair data at the shape methanol growth has reached after round 12 (a snapshot of its
# state, taken 2026-10-19 10:05). Provisional: the shape moved by under 0.01 % in length and 0.02 degrees in
# angle over the last several rounds (against the 0.1 % and 0.1 degrees of E5b), and the pairs use the per-element
# basis, so only the shape matters. If the owner has the growth go on and the shape ever moves past the
# tolerance, these pairs are the ones to redo. Resumable.
cd /e/Projects/phys/runs/e8e/pairs-methanol
log() { echo "[$(date '+%m-%d %H:%M')] $*" >> queue.log; }
try() { for i in 1 2 3 4 5; do "$@" >> run.log 2>&1 && return 0; log "failed (try $i): $*"; sleep 30; done; return 1; }
log "methanol: starting pair data on the round-12 shape"
try ./phys-pairs-gpu.exe methanol 150 --mp2
try ./phys-pairs-unlike-gpu.exe water methanol 40
try ./phys-pairs-gpu.exe methanol 100 --mp2 --wall 1.0 2.3
try ./phys-pairs-unlike-gpu.exe water methanol 40 --wall 1.0 2.3
log "methanol: done"
touch all-done
