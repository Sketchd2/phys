#!/bin/bash
# The ammonia exact-electrostatics job (the input to its cloud fit) and the methanol queue's worker cannot both
# fit in memory with basis growth running: the worker, which learned that a methanol task takes 7.5 GB, was told by
# the server that the 2 GB it found free was too little (exit 4) and stopped. So: let the electrostatics job reach
# 100 pairs (enough for the cloud fit), stop it, and start the worker. Bounded: gives up waiting after 6 hours and
# starts the worker anyway.
cd /e/Projects/phys/runs/e8e
log() { echo "[$(date '+%m-%d %H:%M')] $*" >> alerts.log; }
start=$(date +%s)
while [ "$(wc -l < pairs/es-ammonia.txt)" -lt 60 ] && [ $(( $(date +%s) - start )) -lt 21600 ]; do sleep 60; done
n=$(wc -l < pairs/es-ammonia.txt)
log "HANDOFF: ammonia electrostatics at $n pairs; stopping it and starting the methanol queue's worker"
powershell -NoProfile -Command "Get-CimInstance Win32_Process | Where-Object { \$_.Name -eq 'bash.exe' -and \$_.CommandLine -eq '\"C:\Program Files\Git\usr\bin\bash.exe\" ./es.sh' } | ForEach-Object { Stop-Process -Id \$_.ProcessId -Force }; Stop-Process -Name phys-es-gpu -Force -ErrorAction SilentlyContinue"
sleep 5
cd pairs-methanol
export PHYS_QUEUE_TOKEN=e8e-methanol PHYS_SPILL_DIR=/h/phys-spill PHYS_SPILL_MAX_GB=40
nohup sh ./worker-loop.sh 127.0.0.1:7878 desk ./phys-worker-gpu.exe >> worker.out 2>&1 &
sleep 120
log "HANDOFF done: $(tail -n 1 worker.out | cut -c1-100)"
