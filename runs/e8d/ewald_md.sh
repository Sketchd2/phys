#!/bin/bash
# Ewald dynamics (PHYS_EWALD): the E8c wall-fit law and the no-wall law at 1.0 g/cm3, 298 K.
# usage: ewald_md.sh PS   (picoseconds of production + settling)
cd /e/Projects/phys/runs/e8d
PS=${1:-40}
export PHYS_EWALD=4.2
export PHYS_EWALD_ACC=1e-5
run() { # name law state
  mkdir -p ewald/$1 && cd ewald/$1
  cp "$2" law.txt; cp "$3" grow-water.state
  [ -f done ] || { E:/Projects/phys/target/dev7/release/phys-bulk.exe water law.txt 298 216 1.0 $PS > liquid.log 2>&1 && touch done; }
  cd ../..
}
run opt1 /e/Projects/phys/runs/e8c/avg/refit1/law-opt1.txt /e/Projects/phys/runs/e8c/avg/refit1/grow-water.state
run nowall /e/Projects/phys/runs/e8c/avg/law-final.txt /e/Projects/phys/runs/e8c/avg/grow-water.state
