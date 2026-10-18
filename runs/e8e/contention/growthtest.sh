#!/bin/bash
# How many cores does methanol growth use when the ammonia pair job is suspended? Bounded: 15 minutes, and the
# pair job is resumed from an exit trap however this ends.
cd /e/Projects/phys/runs/e8e/contention
PS="powershell -NoProfile -ExecutionPolicy Bypass -File proc.ps1"
AMM=$(tasklist | grep -i phys-pairs-gpu | head -1 | awk '{print $2}')
cpu() { powershell -NoProfile -Command "(Get-Process phys-grow).CPU" | tr -d '\r'; }
trap '$PS resume $AMM >/dev/null 2>&1; echo "resumed $(date)" >> growthtest.txt' EXIT
echo "start $(date) ammonia pid $AMM" > growthtest.txt
a=$(cpu); sleep 300; b=$(cpu)
echo "with the pair job running: $(awk -v a=$a -v b=$b 'BEGIN{printf "%.2f", (b-a)/300}') cores" >> growthtest.txt
$PS suspend $AMM >/dev/null 2>&1
a=$(cpu); sleep 600; b=$(cpu)
echo "with the pair job suspended: $(awk -v a=$a -v b=$b 'BEGIN{printf "%.2f", (b-a)/600}') cores" >> growthtest.txt
