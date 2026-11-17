#!/bin/bash
# pairs-part.sh NAME K N : part K (1-based) of the pairs stage's frozen-orbital pairs, N pairs a part, run from
# runs/e8f/species/NAME after its selection exists. One process per part: a whole 24-pair run grew until the
# machine killed it (dimethyl ether, 4 pairs in), and it does not resume. Writes frozen-part-K.txt; done
# when it holds its N lines (fewer for the last part). assemble: pairs-part.sh NAME all N writes pairs.ok.
NAME=$1; K=$2; N=$3; ROOT=$(cd "$(dirname "$0")/../../.." && pwd); B=$ROOT/target/release
SEL=pairs-$NAME-sel.txt; [ -s $SEL ] || { echo "no selection"; exit 1; }
TOTAL=$(grep -c . $SEL); PARTS=$(( (TOTAL + N - 1) / N ))
if [ "$K" = all ]; then
  for k in $(seq 1 $PARTS); do [ -f frozen-part-$k.done ] || { echo "part $k not done"; exit 1; }; done
  for k in $(seq 1 $PARTS); do cat frozen-part-$k.txt; done > frozen-$NAME.txt
  echo "$(grep -c . frozen-$NAME.txt) pairs; $(cat select.log)" > pairs.ok; exit 0
fi
[ -f frozen-part-$K.done ] && exit 0
WANT=$(( K*N > TOTAL ? TOTAL - (K-1)*N : N ))
sed -n "$(( (K-1)*N + 1 )),$(( K*N ))p" $SEL > sel-part-$K.txt
export PHYS_SPILL_DIR=/tmp/spill PHYS_SPILL_MAX_GB=0; mkdir -p /tmp/spill; rm -rf /tmp/spill/*
$B/phys-overlap $NAME sel-part-$K.txt --frozen --out frozen-part-$K.txt > frozen-part-$K.log 2>&1 || exit 1
[ $(grep -c . frozen-part-$K.txt) -eq $WANT ] && touch frozen-part-$K.done
