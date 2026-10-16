#!/bin/bash
# Is the structure at 1.0 g/cm3 the equilibrium one? Start 1.0 from the compressed tetrahedral 0.95 configuration, and 0.95 from the expanded dense 1.0 one.
cd /e/Projects/phys/runs/e8d
export PHYS_EWALD=4.2 PHYS_EWALD_ACC=1e-5
run() { # name density snapshot
  mkdir -p ewald/$1 && cd ewald/$1
  cp ../../core/law-core-B7.txt law.txt; cp /e/Projects/phys/runs/e8c/avg/grow-water.state .
  [ -f done ] || { PHYS_START_SNAP=$3 E:/Projects/phys/target/dev12/release/phys-bulk.exe water law.txt 298 216 $2 30 > liquid.log 2>&1 && touch done; }
  cd ../..
}
run core_B7_up1.0_from0.95 1.0 /e/Projects/phys/runs/e8d/ewald/core_B7_d0.95/bulk-water-298.snap
run core_B7_down0.95_from1.0 0.95 /e/Projects/phys/runs/e8d/ewald/core_B7/bulk-water-298.snap
