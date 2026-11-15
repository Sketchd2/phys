# Fleet state (the orchestrator's own notes; newest first)

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
