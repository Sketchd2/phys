#!/bin/bash
# Restart the water runs if nothing is running: each chain skips the states already marked done.
cd "$(dirname "$0")"
if pgrep -f phys-bulk > /dev/null; then echo running; exit 0; fi
if grep -q '^finished3$' chain.log; then echo complete; exit 0; fi
nohup bash -c './chain.sh >> chain.log 2>&1; ./chain2.sh >> chain.log 2>&1; ./chain3.sh >> chain.log 2>&1' > /dev/null 2>&1 &
echo restarted
