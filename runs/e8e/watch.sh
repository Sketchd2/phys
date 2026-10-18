#!/bin/bash
# Stall watchdog for the long jobs in jobs.txt (docs/SOLUTIONS.md 9.3). Every ten minutes it checks, for each job:
#   DEAD       its process is gone and its done-file is absent (twice running, to ride out a restart between steps);
#   STALL      its process is there but none of its logs has been written for the allowed minutes;
#   NOT RUNNING its process is there but has used under 5 s of CPU in the last half hour (suspended, or hung),
#              unless it is waiting on the GPU, which still shows in its log;
# and it warns when free memory is under 1 GB or the disk under 5 GB. Everything goes to alerts.log with a
# time; the file ALERT exists while any alert is unresolved, and ok.log is touched every round so that a dead
# watchdog is itself visible (its mtime is older than ten minutes). Nothing here stops or changes a job.
cd /e/Projects/phys/runs/e8e
declare -A seen cpu_hist
log() { echo "[$(date '+%m-%d %H:%M')] $*" >> alerts.log; }
log "watchdog started"
while true; do
  bad=0
  while IFS='|' read -r name image logs minutes donefile; do
    name=$(echo $name | xargs); [[ "$name" == \#* || -z "$name" ]] && continue
    image=$(echo $image | xargs); logs=$(echo $logs | xargs); minutes=$(echo $minutes | xargs); donefile=$(echo $donefile | xargs)
    [ -f "$donefile" ] && { seen[$name,dead]=0; continue; }
    cpu=$(powershell -NoProfile -Command "(Get-Process -Name ${image%.exe} -ErrorAction SilentlyContinue | Measure-Object CPU -Sum).Sum" 2>/dev/null | tr -d '\r')
    if [ -z "$cpu" ]; then
      seen[$name,dead]=$(( ${seen[$name,dead]:-0} + 1 ))
      if [ ${seen[$name,dead]} -ge 2 ]; then log "DEAD: $name ($image is not running and $donefile does not exist)"; bad=1; fi
      continue
    fi
    seen[$name,dead]=0
    # silence: the newest write among the logs
    newest=0
    for f in ${logs//,/ }; do [ -f "$f" ] && m=$(stat -c %Y "$f") && [ "$m" -gt "$newest" ] && newest=$m; done
    silent=$(( ( $(date +%s) - newest ) / 60 ))
    if [ "$newest" -gt 0 ] && [ "$silent" -gt "$minutes" ]; then log "STALL: $name has written nothing for $silent minutes (allowed $minutes)"; bad=1; fi
    # cpu used over the last three rounds
    hist="${cpu_hist[$name]:-}"; hist="$hist $cpu"; hist=$(echo $hist | awk '{n=NF; s=(n>4)?n-3:1; for(i=s;i<=n;i++) printf "%s ", $i}')
    cpu_hist[$name]="$hist"
    set -- $hist
    if [ $# -ge 4 ]; then used=$(awk -v a="$1" -v b="${@: -1}" 'BEGIN{print b-a}'); if awk -v u="$used" 'BEGIN{exit !(u<5)}'; then log "NOT RUNNING: $name used $used s of CPU in half an hour (suspended or hung)"; bad=1; fi; fi
  done < jobs.txt
  free_gb=$(powershell -NoProfile -Command "[math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB,1)" | tr -d '\r')
  disk_gb=$(df -BG /e 2>/dev/null | awk 'NR==2{gsub("G","",$4); print $4}')
  awk -v f="$free_gb" 'BEGIN{exit !(f<1.0)}' && { log "MEMORY: only $free_gb GB free"; bad=1; }
  [ -n "$disk_gb" ] && [ "$disk_gb" -lt 5 ] && { log "DISK: only $disk_gb GB free on E:"; bad=1; }
  if [ $bad = 1 ]; then touch ALERT; else rm -f ALERT; fi
  touch ok.log
  sleep 600
done
