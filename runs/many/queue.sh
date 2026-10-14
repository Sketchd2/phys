#!/bin/bash
cd /e/Projects/phys
. runs/many/env.sh
B=target/release/phys-cluster-gpu.exe
run() { # centre size
  $B water bulk-water-298.snap --manybody --method hf --centre $1 --size $2 --out runs/many/mb-hf-c$1.txt >> runs/many/mb-hf-c$1.log 2>&1
}
for c in 0 20 45 70 95 120 145 170 195 210; do run $c 6; done
run 0 10
