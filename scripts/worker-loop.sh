#!/bin/sh
# Keep a queue worker going on a Linux machine (a Pi, a rack server).
#
#   scripts/worker-loop.sh host:port [name] [class-binary]
#
# Runs `phys-worker` from the current directory, which must hold the same
# grow-<name>.state files (and any snapshot the plan names) as the server. A
# worker exits 10 after ten tasks so that each stretch runs in a fresh process
# (a long-lived one slows down: see PLAY.md E8); this starts the next one. It
# stops when the server says nothing is left (exit 0) and waits a minute and
# tries again if the server cannot be reached (exit 3), so a rebooted
# server or a network blip does not need anyone to restart the fleet.
#
# Memory: a water pair needs several GB for its tables. If PHYS_SPILL_DIR is
# set it spills to that directory (cap PHYS_SPILL_MAX_GB, default 40); on an SD
# card that is slower and wears the card, so give a Pi an external disk or none.

server=${1:?usage: worker-loop.sh host:port [name] [binary]}
name=${2:-$(hostname)}
bin=${3:-./phys-worker}
unreachable=0
while true; do
    "$bin" "$server" --name "$name" --jobs 10
    case $? in
        0) exit 0 ;;
        10) unreachable=0 ;;
        *) unreachable=$((unreachable + 1))
           [ "$unreachable" -ge 60 ] && exit 3
           sleep 60 ;;
    esac
done
