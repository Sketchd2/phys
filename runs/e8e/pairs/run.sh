#!/bin/bash
# E8e: MP2 pair data that depends on nothing still undecided (docs/SOLUTIONS.md 9.3, S1 and S2):
#   like pairs of each new molecule, uniform draw and wall draw; unlike pairs with water.
# Each molecule waits for its basis growth to settle (runs/e8e/grow/done-<name>), so no pair is
# computed at a shape growth still moves. Resumable: the programs skip the indices they hold.
cd /e/Projects/phys/runs/e8e/pairs
G=/e/Projects/phys/runs/e8e/grow
log() { echo "[$(date '+%m-%d %H:%M')] $*" >> queue.log; }
ready() { [ -f "$G/done-$1" ]; }
state() { [ -f grow-$1.state ] || cp "$G/grow-$1.state" .; }
water() { [ -f grow-water.state ] || cp /e/Projects/phys/grow-water.state .; }
try() { # try <command...> up to 5 times
  for i in 1 2 3 4 5; do "$@" >> run.log 2>&1 && return 0; log "failed (try $i): $*"; sleep 30; done; return 1; }
water
for m in ammonia methanol ethanol methane; do
  while ! ready $m; do sleep 120; done
  state $m
  log "$m: grown, starting pair data"
  try ./phys-pairs-gpu.exe $m 150 --mp2
  try ./phys-pairs-unlike-gpu.exe water $m 40
  try ./phys-pairs-gpu.exe $m 100 --mp2 --wall 1.0 2.3
  try ./phys-pairs-unlike-gpu.exe water $m 40 --wall 1.0 2.3
  log "$m: done"
done
touch all-done
