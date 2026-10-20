#!/bin/bash
# Does running the GPU jobs together cost or gain? Time the ammonia pairs and the ion pairs alone and in
# company by suspending and resuming the other jobs' processes (nothing is killed). The first pair finished
# after each change is discarded, since it straddles it.
cd /e/Projects/phys/runs/e8e/contention
PS="powershell -NoProfile -ExecutionPolicy Bypass -File proc.ps1"
pid() { tasklist | grep -i "$1" | head -1 | awk '{print $2}'; }
GROW=$(pid phys-grow); AMM=$(pid phys-pairs-gpu); ION=$(pid phys-pairs-unlike-gpu)
echo "pids grow $GROW ammonia $AMM ions $ION" > result.txt
AL=/e/Projects/phys/runs/e8e/pairs/run.log; IL=/e/Projects/phys/runs/e8e/ions/run.log
count() { grep -c "^pair" "$1"; }
times() { grep "^pair" "$1" | tail -n "$2" | sed -E 's/.*, ([0-9]+) s \(.*/\1/;t;s/.*, ([0-9]+) s$/\1/' | tr '\n' ' '; }
suspend() { for p in "$@"; do [ -n "$p" ] && $PS suspend $p > /dev/null 2>&1; done; }
resume() { for p in "$@"; do [ -n "$p" ] && $PS resume $p > /dev/null 2>&1; done; }
stage() { # name  logfile-to-watch  pairs-to-measure  [log2 pairs2]
  local name=$1 log=$2 k=$3 log2=$4 k2=$5
  local c0=$(count $log); local d0=0; [ -n "$log2" ] && d0=$(count $log2)
  # wait for one straddling pair then k more
  while [ $(count $log) -lt $((c0 + 1 + k)) ]; do sleep 10; done
  local line="$name: $(times $log $k) s per pair (first log)"
  if [ -n "$log2" ]; then line="$line | $(times $log2 $k2) s per pair (second log)"; fi
  echo "$line" >> result.txt
}
# make sure all three are running to begin with
resume $GROW $AMM $ION
echo "stage 1: ions alone" >> result.txt; suspend $GROW $AMM; stage "ions alone" $IL 8
echo "stage 2: ammonia + ions" >> result.txt; resume $AMM; stage "ammonia with ions (grow suspended)" $AL 3 $IL 12
echo "stage 3: ammonia alone" >> result.txt; suspend $ION; stage "ammonia alone" $AL 3
echo "stage 4: ammonia + grow" >> result.txt; resume $GROW; stage "ammonia with grow (ions suspended)" $AL 3
echo "stage 5: all" >> result.txt; resume $ION; stage "all three" $AL 3 $IL 12
resume $GROW $AMM $ION
echo done >> result.txt
