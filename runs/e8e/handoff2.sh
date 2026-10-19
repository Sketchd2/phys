#!/bin/bash
# After the uniform ammonia pairs' electrostatics (done: 60 pairs), the wall pairs' (es2.sh, in pairs/wall/): let it reach
# 60 pairs - the close-contact pairs are what constrain the lone-pair cloud, which the first fit emptied - then stop it
# and start the methanol worker (memory: the two do not fit together with growth). Bounded at 5 hours.
cd /e/Projects/phys/runs/e8e
log() { echo "[$(date '+%m-%d %H:%M')] $*" >> alerts.log; }
start=$(date +%s)
n() { [ -f pairs/wall/es-ammonia.txt ] && wc -l < pairs/wall/es-ammonia.txt || echo 0; }
while [ "$(n)" -lt 60 ] && [ $(( $(date +%s) - start )) -lt 18000 ]; do sleep 60; done
log "HANDOFF2: ammonia wall electrostatics at $(n) pairs; stopping it and starting the methanol worker"
powershell -NoProfile -Command "Get-CimInstance Win32_Process | Where-Object { \$_.Name -eq 'bash.exe' -and \$_.CommandLine -match 'es2.sh' -and \$_.CommandLine -notmatch 'snapshot' -and \$_.CommandLine -notmatch 'handoff' } | ForEach-Object { Stop-Process -Id \$_.ProcessId -Force }; Stop-Process -Name phys-es-gpu -Force -ErrorAction SilentlyContinue"
sleep 5
cd pairs-methanol
export PHYS_QUEUE_TOKEN=e8e-methanol PHYS_SPILL_DIR=/h/phys-spill PHYS_SPILL_MAX_GB=40
nohup sh ./worker-loop.sh 127.0.0.1:7878 desk ./phys-worker-gpu.exe >> worker.out 2>&1 &
sleep 120
log "HANDOFF2 done: $(tail -n 1 worker.out | cut -c1-100)"
