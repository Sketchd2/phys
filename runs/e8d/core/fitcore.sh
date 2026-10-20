#!/bin/bash
# fitcore.sh NAME B WALLW : the no-wall law with a steep Born-Mayer core (exponent B) fitted to everything with the rest held
NAME=$1; B=$2; WW=$3
mkdir -p $NAME; cd $NAME
cp ../{pairs-water-gpu-bias-mp2.txt,pairs-water-gpu-mp2.txt,pairs-water-gpu-wall-mp2.txt,pairs-water-liq-gpu-mp2.txt,pairs-water-slp.txt,deriv-train.txt,deriv-test.txt,rotate-water.txt,polar-water.txt,esp-water-cloud-mp2-fit.txt,grow-water.state,law-start.txt} .
cp deriv-train.txt deriv-water.txt
E:/Projects/phys/target/dev9/release/phys-fit.exe water 1000 -gpu-mp2,-gpu-bias-mp2,-gpu-wall-mp2,-liq-gpu-mp2 --bisector --alpha polar-water.txt --esp esp-water-cloud-mp2-fit.txt --charge-sites esp-water-cloud-mp2-fit.txt --no-c8 --c6 0-0=0.0,0-1=2.7244184049839726,1-1=3.764045228833158 --slopes deriv-water.txt --slope-pairs -slp --rotations rotate-water.txt --slope-weight 30 --soft-wall 30 --soft-weight $WW --bm-exponent $B --start law-start.txt > fit.log 2>&1
L=$(ls law-water-*-pbe-bis-esp-ind.txt 2>/dev/null | head -1)
echo "== $NAME  core exponent $B wallW=$WW"
if [ -z "$L" ]; then tail -5 fit.log; exit; fi
grep -E "^ +pbe:" fit.log | sed 's/within 3 kcal\/mol of the lowest: //; s/^ */  /'
grep -E "^pair (0 0|0 1|1 1)" "$L" | awk '{printf "  core %s %s: A=%.3e B=%.2f C6=%.3f\n",$2,$3,$4,$5,$6}'
echo "  liquid-range check, held-out liquid slopes and energy errors (law - MP2):"
E:/Projects/phys/target/dev9/release/phys-diag.exe water "$L" --slopes deriv-test.txt -slp 2>&1 | sed -n 3,7p
cp ../pairs-water-liq-gpu-mp2.txt .
E:/Projects/phys/target/dev9/release/phys-diag.exe water "$L" -liq-gpu-mp2 2>&1 | sed -n 4,6p
