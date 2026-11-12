#!/bin/sh
# Density overlap, electron-gas terms and exact electrostatics on stored MP2
# pairs, four species, close pairs only, coarse passes first (resumes by index).
B=../../target/release/phys-overlap
for every in 8 4 2 1; do
  $B methane ../e8d/methane/pairs-methane-gpu-mp2.txt --max-sep 9 --every $every --out overlap-methane.txt
  $B water ../e8c/avg/pairs-water-gpu-bias-mp2.txt --max-sep 7 --every $((every*2)) --out overlap-water-bias.txt
  $B ammonia ../e8e/pairs/pairs-ammonia-gpu-mp2.txt --max-sep 8 --every $every --out overlap-ammonia.txt
  $B methanol ../e8e/pairs-methanol/pairs-methanol-gpu-mp2.txt --max-sep 10 --every $((every/4+1)) --out overlap-methanol.txt
  $B ammonia ../e8e/pairs/pairs-ammonia-gpu-wall-mp2.txt --every $((every*2)) --out overlap-ammonia-wall.txt
  $B water ../e8c/wall-rigid/pairs-water-gpu-wall-mp2.txt --every $((every*2)) --out overlap-water-wall.txt
done
echo finished
