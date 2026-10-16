#!/bin/bash
# Density scan of the hole-free law under Ewald (the fast build): 0.95 and 1.05 g/cm3, 30 ps each.
cd /e/Projects/phys/runs/e8d
export PHYS_EWALD=4.2 PHYS_EWALD_ACC=1e-5
for D in 0.95 1.05; do
  mkdir -p ewald/core_B7_d$D && cd ewald/core_B7_d$D
  cp ../../core/law-core-B7.txt law.txt; cp /e/Projects/phys/runs/e8c/avg/grow-water.state .
  [ -f done ] || { E:/Projects/phys/target/dev11/release/phys-bulk.exe water law.txt 298 216 $D 30 > liquid.log 2>&1 && touch done; }
  cd ../..
done
