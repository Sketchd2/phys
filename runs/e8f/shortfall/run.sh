#!/bin/bash
# The basis bracket (counterpoise and uncorrected) along each S22 dimer's approach line.
B=../../../target/release/phys-shortfall
[ -s water.txt ] || $B water 0.9 0.95 1.0 1.05 1.1 1.2 1.35 1.6 > water.txt 2>&1
[ -s methane.txt ] || $B methane 0.9 1.0 1.1 1.25 1.5 > methane.txt 2>&1
[ -s ammonia.txt ] || $B ammonia 0.9 1.0 1.1 1.25 1.5 > ammonia.txt 2>&1
echo finished
