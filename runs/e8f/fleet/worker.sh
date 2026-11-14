#!/bin/bash
# worker.sh W : one step of worker number W. Call it again and again; each call returns within
# about nine minutes and prints one status line as its LAST line:
#   BUILDING            the programs are being compiled; call again
#   RUNNING <job>       a job is running; call again
#   FINISHED <job>      a job ended (results committed and pushed); call again for the next
#   FAILED <job>        a job ended badly (its log tail is committed); call again for the next
#   IDLE                nothing in this worker's queue yet; call again
#   ALL DONE            the queue says STOP; do not call again
# All the logic is here. The queue is runs/e8f/fleet/queue/W.txt on branch $MAIN, re-read at
# every call, one job a line: "<job id> <directory under runs/e8f> <command...>". A job is run
# from its directory with nohup, its output in fleet/log/<job id>.log. A job must be safe to
# start again (the run scripts skip what is finished), because a reclaimed machine starts over.
# Results go to branch claude/e8f-w-W and nowhere else.
W=$1; MAIN=${MAIN:-claude/laughing-keller-sagqde}; BR=${BR:-claude/e8f-w-$W}
ROOT=$(cd "$(dirname "$0")/../../.." && pwd); F=$ROOT/runs/e8f/fleet; cd "$ROOT"
mkdir -p $F/log $F/done $F/failed /tmp/fleet; PID=/tmp/fleet/pid; CUR=/tmp/fleet/job
push() {  # commit everything under runs/e8f and push to this worker's branch, four tries
  git add -A runs/e8f >/dev/null 2>&1; git commit -qm "fleet worker $W: $1" >/dev/null 2>&1
  for t in 2 4 8 16; do git push -q origin HEAD:$BR 2>/dev/null && return 0; sleep $t; done; return 1
}
wait_for() {  # wait up to 8 minutes for process $1 to end; 0 if it ended
  for i in $(seq 1 96); do kill -0 $1 2>/dev/null || return 0; sleep 5; done; return 1
}
# 0. on a fresh machine, take back what this worker had already pushed (finished jobs, partial runs)
if [ ! -f /tmp/fleet/synced ]; then
  git fetch -q origin $BR 2>/dev/null && git checkout -q FETCH_HEAD -- runs/e8f 2>/dev/null
  touch /tmp/fleet/synced
fi
# 1. the programs
if [ ! -x target/release/phys-bulk ] || [ ! -f /tmp/fleet/built ]; then
  if [ ! -f /tmp/fleet/build.pid ] || ! kill -0 $(cat /tmp/fleet/build.pid) 2>/dev/null; then
    if [ -f /tmp/fleet/build.ok ]; then touch /tmp/fleet/built; else
      nohup bash -c 'cargo build --release > /tmp/fleet/build.log 2>&1 && touch /tmp/fleet/build.ok' > /dev/null 2>&1 &
      echo $! > /tmp/fleet/build.pid
    fi
  fi
  if [ ! -f /tmp/fleet/built ]; then wait_for $(cat /tmp/fleet/build.pid); [ -f /tmp/fleet/build.ok ] && touch /tmp/fleet/built; fi
  if [ ! -f /tmp/fleet/built ]; then echo BUILDING; exit 0; fi
fi
# 2. a job already running on this machine
if [ -f $PID ] && [ -f $CUR ]; then
  J=$(cat $CUR)
  if wait_for $(cat $PID); then
    rm -f $PID $CUR
    if [ -f /tmp/fleet/ok-$J ]; then touch $F/done/$J; push "$J done"; echo "FINISHED $J"
    else tail -30 $F/log/$J.log > $F/failed/$J; push "$J failed"; echo "FAILED $J"; fi
  else push "$J in progress"; echo "RUNNING $J"; fi
  exit 0
fi
# 3. the next job of this worker's queue (the queue is read fresh from the main branch)
git fetch -q origin $MAIN 2>/dev/null; Q=$(git show FETCH_HEAD:runs/e8f/fleet/queue/$W.txt 2>/dev/null)
while read -r J D CMD; do
  [ -z "$J" ] && continue
  [ "$J" = STOP ] && { echo "ALL DONE"; exit 0; }
  [ -f $F/done/$J ] || [ -f $F/failed/$J ] && continue
  # inputs the job needs come from the main branch as it is now
  git checkout -q FETCH_HEAD -- runs/e8f/$D 2>/dev/null; mkdir -p runs/e8f/$D
  echo $J > $CUR; rm -f /tmp/fleet/ok-$J
  ( cd runs/e8f/$D && nohup bash -c "$CMD && touch /tmp/fleet/ok-$J" > $F/log/$J.log 2>&1 & echo $! > $PID )
  sleep 2
  if wait_for $(cat $PID); then
    rm -f $PID $CUR
    if [ -f /tmp/fleet/ok-$J ]; then touch $F/done/$J; push "$J done"; echo "FINISHED $J"
    else tail -30 $F/log/$J.log > $F/failed/$J; push "$J failed"; echo "FAILED $J"; fi
  else push "$J started"; echo "RUNNING $J"; fi
  exit 0
done <<< "$Q"
sleep 240; echo IDLE
