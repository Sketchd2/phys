#!/usr/bin/env bash
# stop each job right after its next pair is written; both resume with the same launch commands
cd /e/Projects/phys/runs/e8c
n(){ grep -vc '^#' "$1"; }
A=avg/pairs-water-gpu-bias-mp2.txt; F=flex/pairs-water-gpu-bias-flex-mp2.txt
a0=$(n $A); f0=$(n $F)
stopa=0; stopf=0
while [ $stopa = 0 ] || [ $stopf = 0 ]; do
  if [ $stopa = 0 ] && [ $(n $A) -gt $a0 ]; then
    powershell -NoProfile -Command "Get-CimInstance Win32_Process | ? { \$_.CommandLine -match 'e8c-queue.sh a|water-shape.sh|runs.e8c.gpubin.phys-pairs' } | % { Stop-Process -Id \$_.ProcessId -Force }"; stopa=1; echo "avg stopped at $(n $A)"; fi
  if [ $stopf = 0 ] && [ $(n $F) -gt $f0 ]; then
    powershell -NoProfile -Command "Get-CimInstance Win32_Process | ? { \$_.CommandLine -match 'e8c-flex-data.sh|gpubin-flex' } | % { Stop-Process -Id \$_.ProcessId -Force }"; stopf=1; echo "flex stopped at $(n $F)"; fi
  sleep 3
done
rm -f queue-a.lock flex/queue.lock; echo paused
