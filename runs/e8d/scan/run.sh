#!/bin/bash
# Density scan of one law: phys-bulk at several densities, 60 ps each. Usage: run.sh LAWFILE TAG
cd /e/Projects/phys/runs/e8d/scan
LAW=$1; TAG=$2
for D in 0.95 1.05 1.10; do
  mkdir -p $TAG-$D && cd $TAG-$D
  cp /e/Projects/phys/runs/e8c/avg/refit1/grow-water.state . 2>/dev/null
  [ -f done ] || { E:/Projects/phys/target/release/phys-bulk.exe water $LAW 298 216 $D 60 > liquid.log 2>&1 && touch done; }
  cd ..
done
