# Fleet state (the orchestrator's own notes; newest first)

## 2026-11-15 18:10 UTC, wake 4

- FOUND: worker.sh's job-start checkout overwrote a worker's own newer files with main's (methane full 12 -> 2 pairs,
  law HOLE from 2 pairs). Fixed: a worker keeps files its own commits wrote. merge-skip.txt excludes 41's stale ammonia
  and its truncated methane full files. Methane lawb/planb queued on 41 (restores 12 pairs from its commit 36db96c).
- Ammonia on 42: 24 pairs done, full 8/12. Methanol full 6/12 on 40. 17-21 still growing.
- Created 44 (DME, formic acid grows, then dichloromethane) and 45 (ethanol, acetaldehyde, then isobutane). Count 11.
- Window allowed, resets 22:10 UTC. Orchestrator cost 3.78 USD.
- Next: check methane law (no HOLE) and plan -> spread its liquid states over idle workers (41 is idle after planb).

## 2026-11-15 17:10 UTC, wake 3

- FOUND: wake 2's merge undid the ammonia/methane moves (whole-directory checkout of w-30, an old-main branch).
  Fixed: merge with `bash runs/e8f/fleet/merge.sh N...` ONLY. worker.sh now prunes files main removed.
- Ammonia redo on 42 (new ids ...b). Methane on 41: grow/mono/pairs done, full12 running, then law, plan2.
- Methanol full12: 4 of 12 pairs (40). 17-21 still in grow (rounds 5-9).
- Created 42 (ammonia, propane, ethylamine, DME, formic acid grows) and 43 (acetonitrile, ethanol, acetaldehyde grows). Count 9.
- Usage: cost_usd now appears for the orchestrator (2.44); worker costs not yet summed. Window reset 17:10.
- Next: rebalance unstarted grow lines onto new workers as count grows; when 17-21 finish law, queue their plan
  (job ids <name>-plan) and liquid states; when methane/methanol plans exist, spread liquid states.

## 2026-11-15 16:08 UTC, wake 2

- No limit stop at wakes 1-2: created workers 40 and 41 (count 7 with 17-21). Sessions in state.json.
- 17-21 still in grow (ethane r5, ethylene r7, methylamine r8, hydrazine r9, h2o2 r7). 17 unpushed 51 min but running.
- 30: water dilute third try still running.
- D1 extra queued at end of 40 (shortfall/more2.sh).
- Next wake: if no limit stop, +2 workers (42, 43) with pending grows: propane, dimethyl_ether, ethanol, acetaldehyde,
  formic_acid, acetonitrile, ethylamine (two each per queue, smallest first). Window resets 17:10 UTC: check for stops.

## 2026-11-15 15:20 UTC, wake 1 (launch)

- Branch claude/e8f-orchestrator made; worker.sh's MAIN default now points at it.
- Merged w-10..15, 17..30. Water: D3 met (33 of 34 states, boil 325.6 K vs 327.5 K target, within 3 K).
  Methanol: grow, mono, pairs done; full had 2 of 8 pairs (worker 22, stopped).
- Running: 17 ethane, 18 ethylene, 19 methylamine, 20 hydrazine, 21 hydrogen_peroxide (all grow, old queues,
  then mono/pairs/full NFULL=8/law and STOP); 30 water dilute state (old queue).
- Stopped: 22-29. Their species' grow is partial (state merged here). Waiting for workers:
  propane, dimethyl_ether, ethanol, acetaldehyde, formic_acid, acetonitrile, ethylamine (grow continues).
- Queue 40 (not yet created): methanol full12/law/plan, then water grow/mono/pairs/full12.
- Queue 41 (not yet created): ammonia pairs24/full12/law/plan, then methane grow..plan.
- Usage: status allowed, resets 17:10 UTC. No cost field in session metadata.
- Next wake: if no limit stop, create workers 40 and 41 (count 7). Then queue the pending grows on further workers
  as the count grows (42..), smallest first. When 17-21 finish (law), queue their plan stage.
- Later: once water's full is done, move water's hand law/plan/boil/liquid to species/water/hand/ and queue water law+plan.
  D1 extra: phys-shortfall --more on methane and ammonia dimers still to queue.
