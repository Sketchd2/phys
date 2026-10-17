#!/bin/bash
cd /e/Projects/phys/runs/e8d
export PHYS_EWALD=4.2 PHYS_EWALD_ACC=1e-5
mkdir -p ewald/core_B7_long60 && cd ewald/core_B7_long60
cp ../../core/law-core-B7.txt law.txt; cp /e/Projects/phys/runs/e8c/avg/grow-water.state .
[ -f done ] || { E:/Projects/phys/target/dev12/release/phys-bulk.exe water law.txt 298 216 1.0 60 > liquid.log 2>&1 && touch done; }
