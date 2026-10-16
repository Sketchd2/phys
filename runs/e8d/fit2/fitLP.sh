#!/bin/bash
# fitLP.sh NAME SIGMA_LP TEMP SLOPEW WALLW : Born-Mayer + overlap with repulsion-only lone-pair sites of width SIGMA_LP
NAME=$1; SG=$2; T=$3; SW=$4; WW=$5
mkdir -p $NAME; cd $NAME
cp ../{pairs-water-gpu-bias-mp2.txt,pairs-water-gpu-mp2.txt,pairs-water-gpu-wall-mp2.txt,pairs-water-liq-gpu-mp2.txt,pairs-water-slp.txt,deriv-train.txt,rotate-water.txt,polar-water.txt,grow-water.state,law0.txt} .
cp ../esp-lp-$SG.txt esp-water-cloud-mp2-fit.txt
cp deriv-train.txt deriv-water.txt; cp ../law0-lp-$SG.txt law0.txt
E:/Projects/phys/target/dev9/release/phys-fit.exe water $T -gpu-mp2,-gpu-bias-mp2,-gpu-wall-mp2,-liq-gpu-mp2 --bisector --alpha polar-water.txt --esp esp-water-cloud-mp2-fit.txt --charge-sites esp-water-cloud-mp2-fit.txt --overlap --no-c8 --slopes deriv-water.txt --slope-pairs -slp --rotations rotate-water.txt --slope-weight $SW --soft-wall 30 --soft-weight $WW --bm --start law0.txt > fit.log 2>&1
L=$(ls law-water-*-pbe-bis-esp-ind.txt 2>/dev/null | head -1)
echo "== $NAME  sigma_lp=$SG T=$T slopeW=$SW wallW=$WW"
if [ -z "$L" ]; then tail -5 fit.log; exit; fi
grep -E "^ +pbe:" fit.log | sed 's/within 3 kcal\/mol of the lowest: //; s/^ */  /'
cp ../deriv-test.txt dtest.txt; cp ../pairs-water-liq-gpu-mp2.txt .
echo "  held-out liquid slopes (law - MP2):"
E:/Projects/phys/target/dev9/release/phys-diag.exe water "$L" --slopes dtest.txt -slp 2>&1 | sed -n 3,7p
echo "  energy error on the 100 liquid pairs (2.8-3.1, 3.1-3.4, 3.4-3.8):"
E:/Projects/phys/target/dev9/release/phys-diag.exe water "$L" -liq-gpu-mp2 2>&1 | sed -n 4,6p
grep -E "^overlap" "$L" | tr '\n' ';'; echo
