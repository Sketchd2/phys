#!/bin/bash
# variant.sh NAME SETS SLOPEW SOFTW [extra flags]   -> fits in NAME/, scores slopes by shell
NAME=$1; SETS=$2; SW=$3; WW=$4; shift 4
mkdir -p $NAME; cd $NAME
cp ../v0/{pairs-water-gpu-bias-mp2.txt,pairs-water-gpu-mp2.txt,pairs-water-gpu-wall-mp2.txt,deriv-water.txt,rotate-water.txt,polar-water.txt,esp-water-cloud-mp2-fit.txt,grow-water.state,law0.txt} .
WALLFLAGS=""
case "$SETS" in *wall*) WALLFLAGS="--soft-wall 30 --soft-weight $WW";; esac
E:/Projects/phys/target/release/phys-fit.exe water 1000 $SETS --bisector --alpha polar-water.txt --esp esp-water-cloud-mp2-fit.txt --charge-sites esp-water-cloud-mp2-fit.txt --overlap --no-c8 --slopes deriv-water.txt --slope-pairs -gpu-bias-mp2 --rotations rotate-water.txt --slope-weight $SW $WALLFLAGS "$@" > fit.log 2>&1
L=$(ls law-water-*-pbe-bis-esp-ind.txt | head -1)
echo "== $NAME  sets=$SETS slopeW=$SW wallW=$WW  $*"
grep -E "^ +pbe:" fit.log | sed 's/^ */  /'
E:/Projects/phys/target/release/phys-diag.exe water "$L" --slopes deriv-water.txt -gpu-bias-mp2 2>&1 | sed -n 3,7p
awk '/^pair/ {if ($2<=1 && $3<=1) c+=$6} END {printf "  dispersion (O-O, O-H, H-H summed pair table) %.2f\n", c}' "$L"
