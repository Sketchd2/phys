#!/bin/bash
# merge.sh N... : copy into this tree only the files worker N's own commits wrote ("fleet worker N:"),
# as they stand at the tip of its branch, except what fleet/merge-skip.txt ("N path-prefix" a line) excludes. A worker branch also carries stale copies of everything
# that was on its source branch when it was made, so checking out whole directories undoes newer work.
cd "$(git rev-parse --show-toplevel)"
for n in "$@"; do R=origin/claude/e8f-w-$n
  git log --format= --name-only --grep "^fleet worker $n:" $R -- runs/e8f/species runs/e8f/shortfall runs/e8f/fleet/done runs/e8f/fleet/failed | sort -u |
  while read -r f; do skip=; while read -r m pre; do [ "$m" = "$n" ] && [[ "$f" == $pre* ]] && skip=1; done < runs/e8f/fleet/merge-skip.txt
    # a liquid state is taken only from the worker whose queue runs it: every branch carries stale
    # in-progress copies of other workers' states (with their "done"), and taking those undid seven
    # finished ethanol states on 2026-11-26.
    if [[ "$f" =~ /liquid/s-([^/-]+)-([^/]+)/ ]]; then
      grep -q "liquid-state.sh [^ ]* ${BASH_REMATCH[1]} ${BASH_REMATCH[2]} " runs/e8f/fleet/queue/$n.txt 2>/dev/null || skip=1; fi
    [ -z "$skip" ] && [ -n "$f" ] && git cat-file -e $R:"$f" 2>/dev/null && { mkdir -p "$(dirname "$f")"; git show $R:"$f" > "$f"; }; done
done
