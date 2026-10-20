#!/bin/bash
# Liquids (1.0 g/cm3, 80 ps) of refit variants, one at a time, after the density scan.
cd /e/Projects/phys/runs/e8d
true
for V in bm_s30w0.01 bm_s100w0.02 bm_s30w002; do
  mkdir -p liq/$V && cd liq/$V
  L=$(ls ../../fit/$V/law-water-*-pbe-bis-esp-ind.txt | head -1)
  cp "$L" law.txt; cp ../../fit/$V/grow-water.state .
  [ -f done ] || { E:/Projects/phys/target/release/phys-bulk.exe water law.txt 298 216 1.0 80 > liquid.log 2>&1 && touch done; }
  cd ../..
done
