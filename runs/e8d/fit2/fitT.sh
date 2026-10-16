#!/bin/bash
# fitT.sh NAME TEMP SLOPEW WALLW [extra phys-fit flags]: weights exp(-(E-Emin)/kT) at TEMP kelvin (high = flat)
NAME=$1; T=$2; SW=$3; WW=$4; shift 4
mkdir -p $NAME; cd $NAME
cp ../{pairs-water-gpu-bias-mp2.txt,pairs-water-gpu-mp2.txt,pairs-water-gpu-wall-mp2.txt,pairs-water-liq-gpu-mp2.txt,pairs-water-slp.txt,deriv-train.txt,rotate-water.txt,polar-water.txt,esp-water-cloud-mp2-fit.txt,grow-water.state,law0.txt} .
cp deriv-train.txt deriv-water.txt
E:/Projects/phys/target/dev9/release/phys-fit.exe water $T -gpu-mp2,-gpu-bias-mp2,-gpu-wall-mp2,-liq-gpu-mp2 --bisector --alpha polar-water.txt --esp esp-water-cloud-mp2-fit.txt --charge-sites esp-water-cloud-mp2-fit.txt --overlap --no-c8 --slopes deriv-water.txt --slope-pairs -slp --rotations rotate-water.txt --slope-weight $SW --soft-wall 30 --soft-weight $WW --bm "$@" > fit.log 2>&1
L=$(ls law-water-*-pbe-bis-esp-ind.txt | head -1)
echo "== $NAME  T=$T slopeW=$SW wallW=$WW $*"
grep -E "^ +pbe:" fit.log | sed 's/within 3 kcal\/mol of the lowest: //; s/^ */  /'
cp ../deriv-test.txt dtest.txt
echo "  held-out liquid slopes (law - MP2):"
E:/Projects/phys/target/dev9/release/phys-diag.exe water "$L" --slopes dtest.txt -slp 2>&1 | sed -n 3,7p
echo "  energy error on the 100 liquid pairs (mean err, kcal/mol; 2.8-3.1, 3.1-3.4, 3.4-3.8):"
cp ../pairs-water-liq-gpu-mp2.txt .
E:/Projects/phys/target/dev9/release/phys-diag.exe water "$L" -liq-gpu-mp2 2>&1 | sed -n 4,6p
awk '/^pair/ && $2<=1 && $3<=1 {c[$2" "$3]=$6} END {printf "  C6 mol-mol = %.2f\n", c["0 0"]+4*c["0 1"]+4*c["1 1"]}' "$L"
