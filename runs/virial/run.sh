#!/bin/bash
cd /e/Projects/phys
. runs/many/env.sh
B=target/release/phys-cluster-gpu.exe
for L in 0.99 1.01; do
  $B water bulk-water-298.snap --manybody --method hf --centre 0 --size 5 --scale $L --out runs/virial/vir-c0-n5-L$L.txt >> runs/virial/vir-c0-n5-L$L.log 2>&1
done
