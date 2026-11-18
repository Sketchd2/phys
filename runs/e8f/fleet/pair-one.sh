#!/bin/bash
# pair-one.sh NAME K J : the frozen stage for pair J alone of pairs-part K (sel-part-K.txt), for a part that was
# killed on memory even at 3 pairs; `pair-one.sh NAME K all` joins frozen-part-K-*.txt into frozen-part-K.txt.
NAME=$1; K=$2; J=$3; ROOT=$(cd "$(dirname "$0")/../../.." && pwd); B=$ROOT/target/release
[ -s sel-part-$K.txt ] || sed -n "$(( (K-1)*3 + 1 )),$(( K*3 ))p" pairs-$NAME-sel.txt > sel-part-$K.txt  # 3-pair parts, as pairs-part.sh cuts them
[ -s sel-part-$K.txt ] || { echo "no sel-part-$K.txt"; exit 1; }; T=$(grep -c . sel-part-$K.txt)
if [ "$J" = all ]; then
  for j in $(seq 1 $T); do [ -f frozen-part-$K-$j.done ] || { echo "pair $j not done"; exit 1; }; done
  for j in $(seq 1 $T); do cat frozen-part-$K-$j.txt; done > frozen-part-$K.txt; touch frozen-part-$K.done; exit 0
fi
[ -f frozen-part-$K-$J.done ] && exit 0
sed -n "${J}p" sel-part-$K.txt > sel-part-$K-$J.txt
export PHYS_SPILL_DIR=/tmp/spill PHYS_SPILL_MAX_GB=0; mkdir -p /tmp/spill; rm -rf /tmp/spill/*
$B/phys-overlap $NAME sel-part-$K-$J.txt --frozen --out frozen-part-$K-$J.txt > frozen-part-$K-$J.log 2>&1 || exit 1
[ $(grep -c . frozen-part-$K-$J.txt) -eq 1 ] && touch frozen-part-$K-$J.done
