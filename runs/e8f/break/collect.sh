#!/bin/bash
# Copy every result row file from the fleet workers' branches into this directory.
cd "$(dirname "$0")/../../.."
for w in "$@"; do git fetch -q origin claude/e8f-w-$w 2>/dev/null || continue
  for f in $(git ls-tree --name-only FETCH_HEAD runs/e8f/break/ | grep -E "/r-.*\.(txt|time)$"); do git show FETCH_HEAD:$f > $f 2>/dev/null; done
  echo "w$w: done $(git ls-tree --name-only FETCH_HEAD runs/e8f/fleet/done/ 2>/dev/null | wc -l) failed $(git ls-tree --name-only FETCH_HEAD runs/e8f/fleet/failed/ 2>/dev/null | wc -l) | $(git log -1 --format='%cr %s' FETCH_HEAD)"
done
