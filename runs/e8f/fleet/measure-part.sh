#!/bin/bash
# measure-part.sh NAME K N : run part K of the pairs stage exactly as pairs-part.sh does (spill off), into
# measure-part-K.* so nothing finished is touched, sampling every 30 s: the process's resident and peak memory
# (VmRSS, VmHWM from /proc), the machine's available memory, /tmp/spill's size, free disk, and the pair it is on.
# Asked by the owner: what kills a pairs run that should fit (the engine rebuilds what it cannot store).
NAME=$1; K=$2; N=$3; ROOT=$(cd "$(dirname "$0")/../../.." && pwd); B=$ROOT/target/release
SEL=pairs-$NAME-sel.txt; sed -n "$(( (K-1)*N + 1 )),$(( K*N ))p" $SEL > measure-sel-$K.txt
export PHYS_SPILL_DIR=/tmp/spill PHYS_SPILL_MAX_GB=0; mkdir -p /tmp/spill; rm -rf /tmp/spill/*
OUT=measure-part-$K; : > $OUT.mem
echo "# t_s rss_GB peak_GB avail_GB spill_GB disk_free_GB pairs_done" >> $OUT.mem
$B/phys-overlap $NAME measure-sel-$K.txt --frozen --out $OUT.txt > $OUT.log 2>&1 & P=$!; T0=$(date +%s)
while kill -0 $P 2>/dev/null; do
  r=$(awk '/VmRSS/{print $2}' /proc/$P/status 2>/dev/null); h=$(awk '/VmHWM/{print $2}' /proc/$P/status 2>/dev/null)
  a=$(awk '/MemAvailable/{print $2}' /proc/meminfo); s=$(du -sk /tmp/spill 2>/dev/null | cut -f1); d=$(df -k --output=avail /tmp | tail -1)
  echo "$(( $(date +%s) - T0 )) $(awk -v r=${r:-0} -v h=${h:-0} -v a=$a -v s=${s:-0} -v d=$d 'BEGIN{printf "%.2f %.2f %.2f %.2f %.1f", r/1048576, h/1048576, a/1048576, s/1048576, d/1048576}') $(grep -c '^pair' $OUT.log)" >> $OUT.mem
  sleep 30; done
wait $P; RC=$?; echo "# exit $RC after $(( $(date +%s) - T0 )) s; dmesg: $(dmesg 2>/dev/null | grep -i 'killed process' | tail -1)" >> $OUT.mem
echo "# peak resident $(awk '!/^#/{if($3>m)m=$3}END{print m}' $OUT.mem) GB; lowest available $(awk '!/^#/{if(n==""||$4<n)n=$4}END{print n}' $OUT.mem) GB" >> $OUT.mem
exit 0
