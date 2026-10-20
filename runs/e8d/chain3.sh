#!/bin/bash
# After the running alpha_0.80 liquid: cut-off study, then the other induction scales.
cd /e/Projects/phys/runs/e8d
while tasklist | grep -q phys-bulk; do sleep 30; done
touch liq/alpha_0.80/done
for RC in 17.5 12; do
  mkdir -p rc/$RC && cd rc/$RC
  cp /e/Projects/phys/runs/e8c/avg/refit1/law-opt1.txt law.txt; cp /e/Projects/phys/runs/e8c/avg/refit1/grow-water.state .
  [ -f done ] || { PHYS_RCUT=$RC E:/Projects/phys/target/dev5/release/phys-bulk.exe water law.txt 298 216 1.0 40 > liquid.log 2>&1 && touch done; }
  cd ../..
done
for V in alpha_0.70 alpha_0.90; do
  mkdir -p liq/$V && cd liq/$V
  L=$(ls ../../fit/$V/law-water-*-pbe-bis-esp-ind.txt | head -1)
  cp "$L" law.txt; cp ../../fit/$V/grow-water.state .
  [ -f done ] || { E:/Projects/phys/target/release/phys-bulk.exe water law.txt 298 216 1.0 40 > liquid.log 2>&1 && touch done; }
  cd ../..
done
