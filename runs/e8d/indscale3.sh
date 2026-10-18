#!/bin/bash
cd /e/Projects/phys/runs/e8d
export PHYS_EWALD=4.2 PHYS_EWALD_ACC=1e-5
run() { # name law density width
  d=ewald/$1; mkdir -p $d; cd $d
  cp ../../indscale/$2 law.txt; cp /e/Projects/phys/runs/e8c/avg/grow-water.state .
  [ -f done ] || { PHYS_IND_WIDTH=$4 E:/Projects/phys/target/dev12/release/phys-bulk.exe water law.txt 298 216 $3 30 > liquid.log 2>&1 && touch done; }
  cd ../..
}
run ind1.2 law-ind1.2.txt 1.0 1.0
run ind1.1_w1.2 law-ind1.1.txt 1.0 1.2
touch ewald/indscale3-done
