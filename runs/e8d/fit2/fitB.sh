#!/bin/bash
# fitB.sh NAME B SLOPEW WALLW : Born-Mayer exponent held at B (bohr^-1); 0 = free
NAME=$1; B=$2; SW=$3; WW=$4
mkdir -p $NAME; cd $NAME
cp ../{pairs-water-gpu-bias-mp2.txt,pairs-water-gpu-mp2.txt,pairs-water-gpu-wall-mp2.txt,pairs-water-liq-gpu-mp2.txt,pairs-water-slp.txt,deriv-train.txt,rotate-water.txt,polar-water.txt,esp-water-cloud-mp2-fit.txt,grow-water.state,law0.txt} .
cp deriv-train.txt deriv-water.txt
EXTRA=""; [ "$B" != "0" ] && EXTRA="--bm-exponent $B"
E:/Projects/phys/target/dev9/release/phys-fit.exe water 1000 -gpu-mp2,-gpu-bias-mp2,-gpu-wall-mp2,-liq-gpu-mp2 --bisector --alpha polar-water.txt --esp esp-water-cloud-mp2-fit.txt --charge-sites esp-water-cloud-mp2-fit.txt --overlap --no-c8 --slopes deriv-water.txt --slope-pairs -slp --rotations rotate-water.txt --slope-weight $SW --soft-wall 30 --soft-weight $WW --bm $EXTRA > fit.log 2>&1
L=$(ls law-water-*-pbe-bis-esp-ind.txt | head -1)
echo "== $NAME  B=$B slopeW=$SW wallW=$WW"
grep -E "^ +pbe:" fit.log | sed 's/within 3 kcal\/mol of the lowest: //; s/^ */  /'
cp ../deriv-test.txt dtest.txt
echo "  held-out liquid slopes (law - MP2):"
E:/Projects/phys/target/dev9/release/phys-diag.exe water "$L" --slopes dtest.txt -slp 2>&1 | sed -n 3,7p
echo "  E8c in-sample slopes (inner, first shell):"
cp ../../fit/v0/deriv-water.txt d8c.txt 2>/dev/null; cp ../pairs-water-gpu-bias-mp2.txt pairs-water-b8.txt
E:/Projects/phys/target/dev9/release/phys-diag.exe water "$L" --slopes d8c.txt -b8 2>&1 | sed -n 3,5p
grep -E "^pair (0 0|0 1|1 1)" "$L" | awk '{printf "  pair %s %s: B=%.3f C6=%.3f\n",$2,$3,$5,$6}'
awk '/^pair/ && $2<=1 && $3<=1 {c[$2" "$3]=$6} END {printf "  C6 mol-mol = %.2f\n", c["0 0"]+4*c["0 1"]+4*c["1 1"]}' "$L"
