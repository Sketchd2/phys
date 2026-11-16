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
W=$1; MAIN=${MAIN:-claude/e8f-orchestrator}; BR=${BR:-claude/e8f-w-$W}
ROOT=$(cd "$(dirname "$0")/../../.." && pwd); F=$ROOT/runs/e8f/fleet; cd "$ROOT"
mkdir -p $F/log $F/done $F/failed /tmp/fleet; PID=/tmp/fleet/pid; CUR=/tmp/fleet/job
# This script as the main branch has it now: a worker cannot be told anything once it is running
# (it reads no messages), so a correction reaches it only by being fetched here.
git fetch -q origin $MAIN 2>/dev/null
if [ -z "$FLEET_UPDATED" ] && git show FETCH_HEAD:runs/e8f/fleet/worker.sh > /tmp/fleet/worker.new 2>/dev/null && ! cmp -s /tmp/fleet/worker.new "$0"; then
  cp /tmp/fleet/worker.new "$0"; FLEET_UPDATED=1 exec bash "$0" "$@"
fi
# One call at a time on a machine: two calls at once each found no job running and started the same
# job twice, and the second copy's output overwrote the first's after the first had finished.
exec 8>/tmp/fleet/call.lock; flock -n 8 || { echo "RUNNING (another call of this worker is active)"; exit 0; }
# The local branch carries the name it is pushed under: a push of HEAD to another name was
# refused for every worker of the first fleet, silently, because its errors were thrown away.
git checkout -q -B $BR 2>/dev/null
PUSHED=""
push() {  # commit everything under runs/e8f and push to this worker's branch, four tries
  git add -A runs/e8f >/dev/null 2>&1; git commit -qm "fleet worker $W: $1" >/dev/null 2>&1
  for t in 2 4 8 16; do git push -q -u origin $BR 2>/tmp/fleet/push.err && return 0
    # a replacement worker starts from the main branch, not from the history it replaces, so its
    # first push is not a fast-forward: record the old history as a parent, keeping this tree (which
    # step 0 already filled from it), and try again
    if grep -q 'non-fast-forward\|fetch first\|rejected' /tmp/fleet/push.err; then
      git fetch -q origin $BR 2>/dev/null && git merge -q -s ours --allow-unrelated-histories --no-edit FETCH_HEAD >/dev/null 2>&1; fi
    sleep $t; done
  PUSHED=" (PUSH FAILED: $(tail -1 /tmp/fleet/push.err | cut -c1-120))"; return 1
}
wait_for() {  # wait up to 5 minutes for process $1 to end; 0 if it ended
  for i in $(seq 1 60); do kill -0 $1 2>/dev/null || return 0; sleep 5; done; return 1
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
      nohup bash -c 'cargo build --release > /tmp/fleet/build.log 2>&1 && touch /tmp/fleet/build.ok' > /dev/null 2>&1 8>&- &
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
    if [ "$(cat /tmp/fleet/rc-$J 2>/dev/null)" = 0 ]; then touch $F/done/$J; push "$J done"; echo "FINISHED $J$PUSHED"; exit 0
    elif [ -f /tmp/fleet/rc-$J ]; then tail -30 $F/log/$J.log > $F/failed/$J; push "$J failed"; echo "FAILED $J$PUSHED"; exit 0; fi
    # no exit status at all: the job was killed from outside (a session stopped at its limit takes
    # its processes with it), which is not the job failing. It is started again below.
  else push "$J in progress"; echo "RUNNING $J$PUSHED"; exit 0; fi
fi
# 3. the next job of this worker's queue (the queue is read fresh from the main branch)
git fetch -q origin $MAIN 2>/dev/null; Q=$(git show FETCH_HEAD:runs/e8f/fleet/queue/$W.txt 2>/dev/null)
while read -r J D CMD; do
  [ -z "$J" ] && continue
  [ "$J" = STOP ] && { echo "ALL DONE"; exit 0; }
  [ -f $F/done/$J ] || [ -f $F/failed/$J ] && continue
  # inputs the job needs come from the main branch as it is now
  # except the files this worker's own commits wrote: the main branch learns of those only at the
  # orchestrator's next merge, and taking its older copy cut methane's twelve full pairs to two
  # between one stage and the next. A file the main branch has since removed (set aside to be
  # redone) goes, unless this worker wrote it: a checkout adds files but never takes one away.
  OWN=$(git log --format= --name-only --grep "^fleet worker $W:" HEAD -- runs/e8f/$D 2>/dev/null | sort -u)
  git ls-tree -r --name-only FETCH_HEAD -- runs/e8f/$D 2>/dev/null | grep -vxF -f <(echo "$OWN") | while read -r f; do
    git checkout -q FETCH_HEAD -- "$f" 2>/dev/null; done; mkdir -p runs/e8f/$D
  git ls-files -co runs/e8f/$D 2>/dev/null | while read -r f; do
    git cat-file -e FETCH_HEAD:"$f" 2>/dev/null || grep -qxF "$f" <<< "$OWN" || rm -f "$f"; done
  echo $J > $CUR; rm -f /tmp/fleet/rc-$J
  ( cd runs/e8f/$D && nohup bash -c "( $CMD ); echo \$? > /tmp/fleet/rc-$J" > $F/log/$J.log 2>&1 8>&- & echo $! > $PID )
  sleep 2
  if wait_for $(cat $PID); then
    rm -f $PID $CUR
    if [ "$(cat /tmp/fleet/rc-$J 2>/dev/null)" = 0 ]; then touch $F/done/$J; push "$J done"; echo "FINISHED $J$PUSHED"
    elif [ -f /tmp/fleet/rc-$J ]; then tail -30 $F/log/$J.log > $F/failed/$J; push "$J failed"; echo "FAILED $J$PUSHED"
    else echo "RUNNING $J (restarting)"; fi
  else push "$J started"; echo "RUNNING $J$PUSHED"; fi
  exit 0
done <<< "$Q"
sleep 240; echo IDLE
