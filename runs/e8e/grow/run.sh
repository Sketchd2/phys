#!/bin/bash
# E8e: grow each molecule's basis to its stopping rule (PLAY.md E5b), one after another, resumable.
# Independent of every law and decision in SOLUTIONS.md: the owner asked (2026-10-18) that the CPU
# and GPU not sit idle while work waits on earlier stages.
cd /e/Projects/phys/runs/e8e/grow
for n in ammonia methanol ethanol methane; do
  [ -f done-$n ] && continue
  ./phys-grow.exe $n >> grow-$n.log 2>&1 && touch done-$n
done
touch all-done
