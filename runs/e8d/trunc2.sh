#!/bin/bash
cd /e/Projects/phys
for rc in 16 20 24 28 32; do
  echo "== no-wall law, r_cut $rc"
  target/release/phys-virial.exe runs/e8c/avg/law-final.txt runs/e8d/g1/nowall-liquid-x8.snap 298 $rc $((rc-2)) | grep -E "non-additive|PRESSURE|with the tail" | sed 's/^ *//'
done
