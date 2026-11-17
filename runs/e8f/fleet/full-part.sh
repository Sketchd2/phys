#!/bin/bash
# full-part.sh NAME K N : stage full (counterpoise, --full --cp) for pairs (K-1)*N+1..K*N of the first NFULL (8)
# selected pairs, in its own process, spill off; `full-part.sh NAME all N` assembles full-NAME.txt and full.ok.
# Whole full runs were killed on memory (chloroform after 2 of 8 pairs), as the pairs stage was.
NAME=$1; K=$2; N=$3; ROOT=$(cd "$(dirname "$0")/../../.." && pwd); B=$ROOT/target/release
[ -s pairs-$NAME-sel.txt ] || { echo "no selection"; exit 1; }
head -${NFULL:-8} pairs-$NAME-sel.txt > pairs-$NAME-full.txt
TOTAL=$(grep -c . pairs-$NAME-full.txt); PARTS=$(( (TOTAL + N - 1) / N ))
if [ "$K" = all ]; then
  for k in $(seq 1 $PARTS); do [ -f full-part-$k.done ] || { echo "part $k not done"; exit 1; }; done
  for k in $(seq 1 $PARTS); do cat full-part-$k.txt; done > full-$NAME.txt
  echo "$(grep -c . full-$NAME.txt) pairs" > full.ok; exit 0
fi
[ -f full-part-$K.done ] && exit 0
WANT=$(( K*N > TOTAL ? TOTAL - (K-1)*N : N ))
sed -n "$(( (K-1)*N + 1 )),$(( K*N ))p" pairs-$NAME-full.txt > full-sel-part-$K.txt
export PHYS_SPILL_DIR=/tmp/spill PHYS_SPILL_MAX_GB=0; mkdir -p /tmp/spill; rm -rf /tmp/spill/*
$B/phys-overlap $NAME full-sel-part-$K.txt --frozen --full --cp --out full-part-$K.txt > full-part-$K.log 2>&1 || exit 1
[ $(grep -c . full-part-$K.txt) -eq $WANT ] && touch full-part-$K.done
