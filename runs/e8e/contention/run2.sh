#!/bin/bash
# GPU contention, second attempt (the first stalled: the ion queue it meant to time had already finished,
# and it left the ammonia and growth jobs suspended for 53 minutes). Now: a sodium job extends the ion
# data (pairs 60-119, useful in itself) as the second GPU job, every stage has a time limit, and every
# process is resumed on exit however the script ends.
cd /e/Projects/phys/runs/e8e/contention
PS="powershell -NoProfile -ExecutionPolicy Bypass -File proc.ps1"
pid() { tasklist | grep -i "$1" | head -1 | awk '{print $2}'; }
GROW=$(pid phys-grow); AMM=$(pid phys-pairs-gpu)
AL=/e/Projects/phys/runs/e8e/pairs/run.log; NL=/e/Projects/phys/runs/e8e/ions/run2.log
count() { [ -f "$1" ] && grep -c "^pair" "$1" || echo 0; }
times() { grep "^pair" "$1" | tail -n "$2" | sed -E 's/.*, ([0-9]+) s \(.*/\1/;t;s/.*, ([0-9]+) s$/\1/' | tr '\n' ' '; }
NA=""
cleanup() { for p in $GROW $AMM $NA; do [ -n "$p" ] && $PS resume $p > /dev/null 2>&1; done; echo "cleanup: resumed" >> result2.txt; }
trap cleanup EXIT
echo "pids grow $GROW ammonia $AMM; started $(date)" > result2.txt
susp() { for p in "$@"; do [ -n "$p" ] && $PS suspend $p > /dev/null 2>&1; done; }
res() { for p in "$@"; do [ -n "$p" ] && $PS resume $p > /dev/null 2>&1; done; }
# wait for k+1 new lines in a log (the first straddles the change), at most 2400 s
stage() { local name=$1 log=$2 k=$3 log2=$4 k2=$5 limit=2400 t=0
  local c0=$(count $log) d0=$(count $log2)
  while [ $(count $log) -lt $((c0 + 1 + k)) ] && [ $t -lt $limit ]; do sleep 10; t=$((t+10)); done
  local line="$name: ammonia/first log $(times $log $k) s per pair"
  [ -n "$log2" ] && line="$line | sodium $(times $log2 $k2) s per pair"
  echo "$line (waited $t s)" >> result2.txt; }
# 1 ammonia alone
susp $GROW; stage "1 ammonia alone" $AL 4
# 2 ammonia + growth
res $GROW; stage "2 ammonia + growth" $AL 4
# 3 ammonia + sodium (growth suspended)
susp $GROW
cd /e/Projects/phys/runs/e8e/ions; ./phys-pairs-unlike-gpu.exe water sodium 120 > run2.log 2>&1 &
cd /e/Projects/phys/runs/e8e/contention
sleep 20; NA=$(pid phys-pairs-unlike-gpu)
stage "3 ammonia + sodium job" $AL 4 $NL 10
# 4 sodium alone
susp $AMM; stage "4 sodium alone" $NL 10
# 5 all three
res $AMM $GROW; stage "5 all three" $AL 4 $NL 10
echo "finished $(date)" >> result2.txt
