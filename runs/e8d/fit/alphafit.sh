#!/bin/bash
# alphafit.sh F : refit (bm, slopeW 30, wallW 0.01) with polarisabilities x F
F=$1; NAME=alpha_$F
mkdir -p $NAME; cd $NAME
cp ../v0/{pairs-water-gpu-bias-mp2.txt,pairs-water-gpu-mp2.txt,pairs-water-gpu-wall-mp2.txt,deriv-water.txt,rotate-water.txt,esp-water-cloud-mp2-fit.txt,grow-water.state,law0.txt} .
python - $F <<'P'
import sys
f=float(sys.argv[1]); out=[]
for l in open('../v0/polar-water.txt'):
    w=l.split()
    out.append('alpha %s %r\n'%(w[1],float(w[2])*f) if w and w[0]=='alpha' else l)
open('polar-water.txt','w').writelines(out)
P
E:/Projects/phys/target/release/phys-fit.exe water 1000 -gpu-mp2,-gpu-bias-mp2,-gpu-wall-mp2 --bisector --alpha polar-water.txt --esp esp-water-cloud-mp2-fit.txt --charge-sites esp-water-cloud-mp2-fit.txt --overlap --no-c8 --slopes deriv-water.txt --slope-pairs -gpu-bias-mp2 --rotations rotate-water.txt --slope-weight 30 --soft-wall 30 --soft-weight 0.01 --bm > fit.log 2>&1
L=$(ls law-water-*-pbe-bis-esp-ind.txt | head -1)
echo "== alpha x $F"; grep -E "^ +pbe:" fit.log | sed 's/within 3 kcal\/mol of the lowest: //; s/^ */  /'
E:/Projects/phys/target/release/phys-diag.exe water "$L" --slopes deriv-water.txt -gpu-bias-mp2 2>&1 | sed -n 3,7p | awk '{printf "   %s slope err %s (n=%s)\n",$1,$4,$2}'
E:/Projects/phys/target/release/phys-virial.exe "$L" /e/Projects/phys/runs/e8d/g1/wall-liquid.snap 298 | grep -E "non-additive|PRESSURE|polaris"
