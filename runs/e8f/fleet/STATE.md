# Fleet state (the orchestrator's own notes; newest first)

## WAKE CADENCE: every 62 min, self-scheduled. At the start of each wake call send_later (+62 min) with the same message; backstop task every 6 h re-arms if the chain breaks.

## 2026-11-16 00:15 UTC, wake 10

- Worker 19 (replacement) never pushed: archived, queue 19 STOP. Methylamine on new worker 56 (grow3..plan2, NFULL 12).
- Worker 45 idle again at 23:08 -> second resume. If it stops a third time, replace it (same number 45).
- Ammonia liquid 29/34; ammonia-boil at end of 46 (if it fails for want of states, requeue as ammonia-boil2).
- Hydrazine pairs done (20; old queue -> full 8, law). Acetonitrile full (43), formic acid pairs (49), ethylene full (18).
- Created 56 (methylamine), 57 (acetic acid), 58 (chloroform). Count 23 (19 retired). Next wake 01:16 UTC.

## 2026-11-15 23:11 UTC, wake 9

- Water stage law HOLE too: 4 of 5 stage laws HOLE (methane, methanol, H2O2, water); ammonia not HOLE, suspect.
- D1 extra finished (40): shortfall/D1-summary.md. Ammonia's bracket 1.8% wide, middle -1.2%; methane 10% wide,
  middle -3.7%. merge.sh now merges runs/e8f/shortfall.
- Replacement 19 still refused: shallow clone needs --allow-unrelated-histories (added). Check w-19 next wake;
  if still nothing, archive it and requeue methylamine on a new number.
- Acetonitrile at full (43), formic acid at pairs (49), hydrazine at pairs (20), ethylene at full (18).
- 40 -> isopropanol grow. Created 54 (acetone), 55 (ethylene glycol). Count 21. Next wake 00:12 UTC.

## 2026-11-15 22:14 UTC, check after the window reset (owner asked)

- No worker stopped at the limit; new window allowed, resets 03:10 UTC. 52, 53 building.
- Worker 45 had ended its turn at 21:39 (not a limit stop): resume sent.
- Worker 19 replacement: every push since 21:15 refused (non-ff); the fixed worker.sh reaches it from its next call
  (its 22:05 call started before the fix). merge -s ours tested here against w-19: works.

## 2026-11-15 22:08 UTC, wake 8

- Laws so far: methane HOLE, methanol HOLE, hydrogen_peroxide HOLE (8 pairs), ammonia not HOLE but suspect.
- Replacement 19 had every push refused (non-ff); worker.sh push() now merges -s ours and retries. Check next wake
  that w-19 pushes again. Methylamine full set to NFULL=12.
- Formic acid grown (49) -> mono..plan on 49. Hydrazine grown (20, old queue continues). Ammonia liquid 20/34.
- Water law/plan on 40 after the D1 extra (still running, 46 min).
- Created 52 (butane grow), 53 (propanol grow). Count 19. Five-hour window passed with no limit stop.

## 2026-11-15 21:09 UTC, wake 7

- PATTERN FOR THE OWNER: two of three stage laws so far begin HOLE (methane 12 pairs, methanol 12 pairs; B2 -inf),
  and ammonia's is not HOLE but has attractive "repulsive" pairs. law.py gives laws with an attractive hole at
  contact; HANDOFF names a core (newlaw/core.py) as the remedy. Not the pipeline's to add (owner decision); by the
  brief each such species stops at law and is recorded. Data for D1/D2 is unaffected.
- Water stages grow..full done (40); water law + plan queued on 40 (marker moves hand law aside on the worker).
  D3 result (325.6 K) kept at species/water/hand/.
- Ammonia liquid 16/34. Methanol: law HOLE. Acetonitrile at pairs (43). H2O2 at full (21), ethylene at pairs (18).
- Worker 19 replaced (new session, my branch); old archived. Created 50 (dichloromethane) and 51 (isobutane). Count 17.

## 2026-11-15 20:08 UTC, wake 6

- Grown: acetonitrile (SETTLED r9) -> 43 runs mono, pairs, full (NFULL 8, 3 heavy), law, plan. Ethylene grow ended
  unsettled (phys-grow exited, r8); 18 on pairs. H2O2 grown, pairs done (21). Both on old queues -> full (8), law.
- Ammonia liquid 7/34 (46, 47). Methanol full on 40 still. 19 resumed (pushed 54 min ago).
- Created 48 (ethylamine grow) and 49 (formic acid grow). Count 15. Window allowed, resets 22:10. Orchestrator 5.78 USD.
- When 17-21 reach law (STOP on old queue): queue their plan + liquid on new workers (not on 17-21, old main).

## 2026-11-15 19:09 UTC, wake 5

- Methane: law from 12 pairs is HOLE -> stops (recorded). Ammonia: law K 7.00 damping 1.0, suspect (see DECISIONS);
  plan 34 states, spread on 46/47. Queue ammonia-boil once >= 31 of 34 are done (NOT on 41: 41 holds stale ammonia files).
- Water: 34/34 states done (dilute rerun finished on 30); boil.ok 325.6 K stands (from 33). Re-run boil? not needed.
- Worker 19 idle since 17:33 -> resume sent 19:08. Check next wake it pushed.
- 41 now: acetaldehyde grow, butane grow. 45: ethanol, isobutane. Created 46, 47. Count 13.
- Methanol full 9/12 on 40, then law, plan (on 40), then water stages.

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
