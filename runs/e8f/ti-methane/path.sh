#!/bin/sh
# The free-energy path for methane's liquid, by states that all exist: down a
# supercritical isotherm from the dilute gas to the liquid's density, then
# down the isochore to the liquid's temperature. One phys-bulk run a state.
# usage: path.sh N PS  (reads states from states.txt: "T density")
N=$1; PS=$2
B=$(pwd)/../../../target/release/phys-bulk
run() {
  d=s-$1-$2; mkdir -p $d; [ -f $d/done ] && return
  ( cd $d && cp ../law.txt ../../grow-methane.state . && $B methane law.txt $1 $N $2 $PS > liquid.log 2>&1 && touch done )
}
while read T D; do
  run $T $D &
  if [ $(jobs -r | wc -l) -ge 2 ]; then wait -n 2>/dev/null || wait; fi
done < states.txt
wait
echo finished
