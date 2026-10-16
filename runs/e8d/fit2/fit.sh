#!/bin/bash
# fit.sh NAME SLOPEW WALLW [extra]  -> a Born-Mayer + overlap fit that includes the liquid pairs and their slopes
NAME=$1; SW=$2; WW=$3; shift 3
mkdir -p $NAME; cd $NAME
cp ../{pairs-water-gpu-bias-mp2.txt,pairs-water-gpu-mp2.txt,pairs-water-gpu-wall-mp2.txt,pairs-water-liq-gpu-mp2.txt,pairs-water-slp.txt,deriv-train.txt,rotate-water.txt,polar-water.txt,esp-water-cloud-mp2-fit.txt,grow-water.state,law0.txt} .
cp deriv-train.txt deriv-water.txt
E:/Projects/phys/target/release/phys-fit.exe water 1000 -gpu-mp2,-gpu-bias-mp2,-gpu-wall-mp2,-liq-gpu-mp2 --bisector --alpha polar-water.txt --esp esp-water-cloud-mp2-fit.txt --charge-sites esp-water-cloud-mp2-fit.txt --overlap --no-c8 --slopes deriv-water.txt --slope-pairs -slp --rotations rotate-water.txt --slope-weight $SW --soft-wall 30 --soft-weight $WW --bm "$@" > fit.log 2>&1
L=$(ls law-water-*-pbe-bis-esp-ind.txt | head -1)
echo "== $NAME  slopeW=$SW wallW=$WW"
grep -E "^ +pbe:" fit.log | sed 's/within 3 kcal\/mol of the lowest: //; s/^ */  /'
echo "  held-out G1 slopes (odd pairs, never fitted), law - MP2:"
cp ../deriv-test.txt dtest.txt; cp ../pairs-water-slp.txt pairs-water-slp.txt
E:/Projects/phys/target/release/phys-diag.exe water "$L" --slopes dtest.txt -slp 2>&1 | sed -n 3,7p
awk '/^pair/ && $2<=1 && $3<=1 {c[$2" "$3]=$6} END {printf "  C6 mol-mol = %.2f\n", c["0 0"]+4*c["0 1"]+4*c["1 1"]}' "$L"
