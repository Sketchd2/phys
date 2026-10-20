#!/bin/bash
# Screening liquids (1.0 g/cm3, 40 ps): the Born-Mayer refits with weaker induction.
cd /e/Projects/phys/runs/e8d
for V in alpha_0.80 alpha_0.70 alpha_0.90; do
  mkdir -p liq/$V && cd liq/$V
  L=$(ls ../../fit/$V/law-water-*-pbe-bis-esp-ind.txt | head -1)
  cp "$L" law.txt; cp ../../fit/$V/grow-water.state .
  [ -f done ] || { E:/Projects/phys/target/release/phys-bulk.exe water law.txt 298 216 1.0 40 > liquid.log 2>&1 && touch done; }
  cd ../..
done
