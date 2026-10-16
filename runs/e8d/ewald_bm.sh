#!/bin/bash
# Ewald dynamics (30 ps) of the Born-Mayer refits, after the no-wall run has finished.
cd /e/Projects/phys/runs/e8d
until [ -f ewald/nowall/done ]; do sleep 60; done
export PHYS_EWALD=4.2
export PHYS_EWALD_ACC=1e-5
run() { # name law state
  mkdir -p ewald/$1 && cd ewald/$1
  cp "$2" law.txt; cp "$3" grow-water.state
  [ -f done ] || { E:/Projects/phys/target/dev7/release/phys-bulk.exe water law.txt 298 216 1.0 30 > liquid.log 2>&1 && touch done; }
  cd ../..
}
S=/e/Projects/phys/runs/e8c/avg/refit1/grow-water.state
run bm_s30w001 "$(ls /e/Projects/phys/runs/e8d/fit/bm_s30w0.01/law-water-*-pbe-bis-esp-ind.txt | head -1)" $S
run bm_g1_s100w002 "$(ls /e/Projects/phys/runs/e8d/fit2/f_s100_w0.02/law-water-*-pbe-bis-esp-ind.txt | head -1)" $S
