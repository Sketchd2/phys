#!/bin/bash
cd /e/Projects/phys
for rc in 16 20 24 28 32; do
  echo "== r_cut $rc (r_on $((rc-2)))"
  target/release/phys-virial.exe runs/e8c/avg/refit1/law-opt1.txt runs/e8d/g1/wall-liquid-x8.snap 298 $rc $((rc-2)) | grep -E "non-additive|pair virial \(all|many-body induction virial|PRESSURE|with the tail"
done
