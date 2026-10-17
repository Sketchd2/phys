#!/bin/bash
cd /e/Projects/phys/runs/e8d
export PHYS_EWALD=4.2 PHYS_EWALD_ACC=1e-5
for f in 0.9 0.8; do
  d=ewald/ind$f; mkdir -p $d; cd $d
  cp ../../indscale/law-ind$f.txt law.txt; cp /e/Projects/phys/runs/e8c/avg/grow-water.state .
  [ -f done ] || { E:/Projects/phys/target/dev12/release/phys-bulk.exe water law.txt 298 216 1.0 30 > liquid.log 2>&1 && touch done; }
  cd ../..
done
touch ewald/indscale-all-done
