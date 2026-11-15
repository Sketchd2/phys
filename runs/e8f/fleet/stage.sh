#!/bin/bash
# stage.sh STAGE NAME : one stage of one species' pipeline, run from runs/e8f/species/NAME.
# Every stage is safe to start again (it skips what it has finished) and ends by writing
# STAGE.ok, which is what the next stage and the orchestrator look for. The species' line in
# fleet/species.tsv gives its atoms and bonds.
#   grow     the basis and relaxed shape (grow-NAME.state, grow-NAME.atoms), by the owner's
#            stopping rule (settled.py), at most MAXH hours
#   mono     the single molecule's record: density and orbital energies, charge clouds,
#            dispersion coefficients, polarisability
#   pairs    600 drawn pair geometries, COUNT (default 24) chosen across the contact range, and
#            each one's electrostatics, density overlaps and exchange from frozen orbitals
#   full     the first NFULL (default 8) chosen pairs again with the pair's own Hartree-Fock and
#            MP2, uncorrected and counterpoise corrected (only sensible for small molecules)
#   law      the species' law from its record and the full stage's energies (law.py): law.txt,
#            and law.ok, which begins HOLE if the law must not be run as a liquid
#   plan     the liquid states, chosen from the law's own second virial coefficient
#            (plan-liquid.py): liquid/states.txt, one "T density ps" a line, each of which is
#            then a job of its own, `../../fleet/liquid-state.sh NAME T DENSITY PS`
#   boil     the boiling point from the finished states (path.py): boil.ok
STAGE=$1; NAME=$2; ROOT=$(cd "$(dirname "$0")/../../.." && pwd); F=$ROOT/runs/e8f/fleet; B=$ROOT/target/release
SPEC=$(awk -F'\t' -v n=$NAME '$1==n {print $3}' $F/species.tsv); [ -z "$SPEC" ] && { echo "no species $NAME"; exit 1; }
export PHYS_SPILL_DIR=${PHYS_SPILL_DIR:-/tmp/spill}; mkdir -p $PHYS_SPILL_DIR
[ -f $STAGE.ok ] && exit 0
need() { [ -f $1.ok ] || { echo "stage $1 has not finished for $NAME"; exit 1; }; }
case $STAGE in
grow)
  # a species already grown at the repository root starts from that state
  [ -f grow-$NAME.state ] || cp $ROOT/grow-$NAME.state . 2>/dev/null
  $F/grow.sh $NAME "$SPEC" || exit 1
  [ -f grow-$NAME.atoms ] || $B/phys-grow "$SPEC" $NAME --atoms-only > /dev/null 2>&1
  [ -f grow-$NAME.atoms ] || { echo "no atoms file"; exit 1; }
  cat grow.status > grow.ok ;;
mono)
  need grow
  [ -s density-$NAME.txt ] || $B/phys-density $NAME > density-$NAME.txt 2> density.log || exit 1
  [ -s cloud-$NAME.txt ] || $B/phys-cloud $NAME > cloud.log 2>&1 || exit 1
  [ -s c6.log ] && grep -q "molecule-molecule" c6.log || $B/phys-c6 $NAME > c6.log 2>&1 || exit 1
  # the polarisability tensor is printed before the split among atoms, which has hung (ammonia):
  # forty minutes, then the tensor alone is kept
  if ! grep -q "MP2: mean" polar.log 2>/dev/null; then timeout 2400 $B/phys-polar $NAME > polar.log 2>&1; fi
  grep -q "MP2: mean" polar.log || { echo "no polarisability"; exit 1; }
  { echo "orbital energies: $(grep '^orbital' density-$NAME.txt | awk '{printf "%s ", $3}')"; grep "molecule-molecule" c6.log; grep "MP2: mean" polar.log; } > mono.ok ;;
pairs)
  need mono
  N=$(wc -w < grow-$NAME.atoms)
  [ -s pairs-$NAME-geom.txt ] && [ $(grep -vc '^#' pairs-$NAME-geom.txt) -ge 600 ] || $B/phys-pairs $NAME 600 --geometry --from 250 > geom.log 2>&1 || exit 1
  [ -s pairs-$NAME-sel.txt ] || python3 $F/select.py pairs-$NAME-geom.txt $N ${COUNT:-24} > pairs-$NAME-sel.txt 2> select.log || exit 1
  $B/phys-overlap $NAME pairs-$NAME-sel.txt --frozen --out frozen-$NAME.txt >> frozen.log 2>&1 || exit 1
  echo "$(grep -c . frozen-$NAME.txt) pairs; $(cat select.log)" > pairs.ok ;;
full)
  need pairs
  head -${NFULL:-8} pairs-$NAME-sel.txt > pairs-$NAME-full.txt
  $B/phys-overlap $NAME pairs-$NAME-full.txt --frozen --full --cp --out full-$NAME.txt >> full.log 2>&1 || exit 1
  echo "$(grep -c . full-$NAME.txt) pairs" > full.ok ;;
law)
  need full
  python3 $F/law.py $NAME ${REF:-cp} > law.log 2>&1 || exit 1 ;;
plan)
  need law
  grep -q '^HOLE' law.ok && { echo "the law has a hole; no liquid"; exit 1; }
  python3 $F/plan-liquid.py $NAME > plan.log 2>&1 || exit 1 ;;
boil)
  need plan
  WANT=$(grep -c . liquid/states.txt); HAVE=$(ls liquid/s-*/done 2>/dev/null | wc -l)
  [ $HAVE -ge $((WANT * 9 / 10)) ] || { echo "$HAVE of $WANT liquid states are finished"; exit 1; }
  M=$(awk '$1=="M"{print $2}' liquid/setup.txt); TH=$(awk '$1=="Thot"{print $2}' liquid/setup.txt); RR=$(awk '$1=="rho_ref"{print $2}' liquid/setup.txt)
  python3 $ROOT/runs/e8f/path.py liquid $M $TH $RR "$(cat liquid/b2.arg)" > boil.log 2>&1 || exit 1
  grep "crosses one atmosphere" boil.log | sed 's/.*atmosphere at //' > boil.ok
  [ -s boil.ok ] || { echo "no crossing found" > boil.ok; } ;;
*) echo "no stage $STAGE"; exit 1 ;;
esac
