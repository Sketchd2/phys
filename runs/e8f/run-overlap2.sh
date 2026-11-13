#!/bin/bash
# Exchange repulsion, three density overlaps and exact electrostatics on stored MP2 pairs of
# ammonia and methane, close pairs only, coarse passes first. Each output file is appended to
# and the program resumes by pair index, so this can simply be run again after a stop.
cp -n ../../grow-ammonia.state ../../grow-methane.state . 2>/dev/null
B=../../target/release/phys-overlap
for every in 8 4 2 1; do
  $B ammonia ../e8e/pairs/pairs-ammonia-gpu-mp2.txt --max-sep 8 --every $every --out overlap-ammonia.txt
  $B methane ../e8d/methane/pairs-methane-gpu-mp2.txt --max-sep 9 --every $every --out overlap-methane.txt
  $B ammonia ../e8e/pairs/pairs-ammonia-gpu-wall-mp2.txt --every $((every*2)) --out overlap-ammonia-wall.txt
done
echo finished
