#!/bin/bash
# The ball-and-stick liquid: 216 flexible waters at 1.0 g/cm3, 298 K, 20 ps, three pair laws.
cd /e/Projects/phys/runs/e8d/sticks
for V in e6disp gascharge freecharge; do
  mkdir -p run-$V && cd run-$V
  [ -f done ] || { E:/Projects/phys/target/dev8/release/phys-sticks.exe ../params-$V.txt 298 216 1.0 20 > out.log 2>&1 && touch done; }
  cd ..
done
