# Fleet state (the orchestrator's own notes; newest first)

## WORKER COUNT: 37 from 19:05 UTC (blocks to 13:20 and 18:20 clean at 31, 35). Block ends 23:20 UTC. Every species is now started: add workers only to split real work (e.g. liquid states); if the limit hits: back to 31 and hold. The orchestrator's own context (409k) is the main cost: keep wakes minimal. Each wake: sample orchestrator + a few workers' cost (get_session) into usage.json; any worker rising >1 USD between wakes is looping: fix the cause. Full snapshot of every session once a block.

## WAKE CADENCE: every 62 min, self-scheduled. At the start of each wake call send_later (+62 min) with the same message; backstop task every 6 h re-arms if the chain breaks.

## 2026-11-16 22:16 UTC, wake 30

- DME parts done 1, 4, 6, 8; 2 on 44, 3 moved to 95, 5 on 93, 7 on 94. Chloroform (79), ethylene glycol (55) mono done,
  now on pairs: if killed like DME, switch to pairs-part.sh parts (selection exists once pairs starts).
- Resumed 76, 81, 82, 85, 88 (they stopped ~10 min after the last resume). Block ends 23:20. Next 23:18.

## 2026-11-16 21:12 UTC, wake 29

- Quiet hour. DME parts 4, 6, 8 running (93-95); 2 on 44. Acetonitrile 32/34.
- 7 idle (69, 76, 81, 82, 84, 85, 88), resumed. POLICY CHANGE: fresh workers stop too (1-3 h), so retiring to a new
  number gains nothing and costs a rebuild: just resume every idle worker each wake. Next 22:14.

## 2026-11-16 20:08 UTC, wake 28

- Acetic acid grown -> mono (4 h polar)..plan on 57.
- DME pairs killed after 4 of 24 pairs (memory grows over one process; no resume). Now 8 parts of 3 via
  fleet/pairs-part.sh: part 1 seeded from the 4 done; 2-3 on 44, 4-5 on 93, 6-7 on 94, 8 on 95 (new; count 40).
  WHEN ALL frozen-part-K.done EXIST: run here `../../fleet/pairs-part.sh dimethyl_ether all 3` in species/dimethyl_ether
  (writes frozen + pairs.ok), commit, then queue full..plan (watch full for the same memory kill).
- Resumed 42, 76, 82, 88 (first stops). Acetonitrile 32/34. Next 21:11.

## 2026-11-16 19:05 UTC, wake 27

- Block 13:20-18:20 PASSED CLEAN at 35. Orchestrator +16.66 USD over it (vs 7.42 the block before): my context is 409k
  tokens and every call re-reads it. Full worker snapshot skipped (stale mid-turn; ~70k tokens), logged.
- 35 -> 37: 89 chlorobenzene, 90 octane (the last unstarted species). 77, 80 (second stop) -> 91 acetone, 92 butanol.
- Acetonitrile 32/34. Next 20:07.

## 2026-11-16 18:02 UTC, wake 26

- Quiet hour: no stage finished. DME at pairs2 (44); chloroform, ethylene glycol on mono (4 h polar).
- 72, 73 (second stop) retired to 87 (hexane), 88 (cyclohexane). 43's queue was done: takes acetonitrile 462/0.627
  (from 83) and phenol. Resumed 77, 80. Count 35. Unstarted: chlorobenzene, octane.
- Block ends 18:20 UTC; judge it at the 19:04 wake (full snapshot, then +4 if clean and there is work).

## 2026-11-16 16:58 UTC, wake 25

- ACETONITRILE BOILING POINT: 414.9 K from 31 of 34 states (run here, stage.sh boil; expt 354.8 K). The other 3 states
  (43, 83) are running: when 34/34, rerun the boil HERE (worker copies may be stale), compare.
- Cause of both ammonia's and acetonitrile's truncated states: a worker committed the whole runs/e8f tree, including
  other workers' half-finished copies it took from main at job start; merging later workers laid those over finished
  ones. worker.sh now commits only files written since its job started (stamp; tested). merge-skip lines added for
  the 7 acetonitrile states. A boil should be run here, not on a worker.
- Chloroform (79) and ethylene glycol (55) grown -> mono (4 h polar) .. plan queued. DME mono2 running (44).
- 46 -> aniline grow. 49, 74 (second stop) retired to 86 (benzene), 85 (ethyl acetate). Resumed 73. Count 35.
- Unstarted: phenol, chlorobenzene, octane. Block ends 18:20 UTC: next wake (18:01) is before it; judge at 19:0x.

## 2026-11-16 15:55 UTC, wake 24

- Merge list: see state.json merge_workers (61 retired -> 84, ethanol grow).
- CARBON TETRACHLORIDE STOPS AT PAIRS: phys-overlap killed twice (18 GB fit table, 15 GB worker; no-spill also killed).
  species/carbon_tetrachloride/STOPPED.txt. Owner's decision (bigger machine or smaller basis); not retried.
- DME mono failed the same way CCl4's first did (polar timeout): mono2 with 4 h polar + pairs2..plan2 on 44.
- ACETONITRILE 31/34 -> acetonitrile-boil queued on 46. Expect a boiling point next wake.
- 65 -> heptane grow; 69 -> glycerol grow (its states were done). Resumed 72, 74, 49 (first time). 69 had been blocked
  once by the auto-mode classifier; resumed.
- Unstarted: aniline, phenol, chlorobenzene, octane. Orchestrator 50.83 USD (+6.0; heavy wakes). Count 35. Next 16:57.

## 2026-11-16 14:52 UTC, wake 23

- Merge list now (retired workers dropped after their final merge): 17 40 41 42 43 44 46 47 48 49 50 52 53 55 57 59 61 62
  65 69 72 73 74 75 76 77 78 79 80 81 82 83.
- Retired 51, 54, 56, 58, 60, 63, 64, 68: idle at every wake, and each resume restarted the job on a fresh machine.
  Their remaining jobs moved to new workers 76-83 (ids -wNN); old queues STOP, sessions archived. Fresh workers
  (72-75) loop fine. Count stays 35. Resumed 61 (first time idle).
- CCl4 pairs killed: 18.1 GB fit table, spill disk full, 15 GB machine. Retry pairs3 with PHYS_SPILL_MAX_GB=0 (rebuild,
  no spill). If killed again, CCl4 stops at pairs (machine-size limit) for the owner.
- Acetonitrile 27/34 (boil at >=31). DME mono running (44). Next wake 15:54.

## 2026-11-16 13:48 UTC, wake 22

- Block 08:20-13:20 PASSED CLEAN at 31. Full snapshot (usage.json): orchestrator +7.42 USD over the block; workers'
  cost_usd only refreshes when a turn ends, so long-turn workers read stale; the ones that refreshed rose ~0.25 each
  -> workers ~8 USD est. Window ~15 USD vs ceiling est. 40-50. Count 31 -> 35: 72 hexane, 73 cyclohexane,
  74 ethyl acetate, 75 toluene (grows).
- Dimethyl ether grown -> mono..plan on 44. Dichloromethane law HOLE (8th). Pyridine grow on 50 (its queue was done).
- Acetonitrile 23/34. CCl4 pairs2 running on 65 (mono2 passed with the 4 h polar). Resumed 56, 58, 60, 63, 64, 68.
- Unstarted: glycerol, heptane, aniline, phenol, chlorobenzene, octane. Next wake 14:51.

## 2026-11-16 12:44 UTC, wake 21

- AMMONIA BOILING POINT: 577.2 K from all 34 liquid states (ammonia-boil4 on 46; reruns complete). Experiment 239.8 K.
  First species through the whole pipeline. (Provisional 578.2 K from 30 states removed.)
- My merged copy had two truncated ammonia states again: worker 47 committed stale copies 16 h ago and merged after 46.
  merge-skip.txt: "47 runs/e8f/species/ammonia/". Merged tree now 34/34 complete.
- 46's queue was empty: 6 unstarted acetonitrile states moved to it from 43/68/69 (ids -w46). Acetonitrile 18/34.
- 8 idle again (51, 54, 56, 58, 60, 63, 64, 68): resumed. The reminder line is not enough; if it continues, retire repeat stoppers to new numbers.
- CCl4 mono2 (4 h polar) running on 65. Orchestrator 43.90 USD. Block ends ~13:20; next wake 13:46 judges it.

## 2026-11-16 11:40 UTC, wake 20

- ammonia-boil3 failed: path.py found 4 more states marked done with truncated logs (1080/0.455, 238/0.605, 347/0.64,
  401/0.712; from before the liquid-state lock). Run here on the 30 complete states it gives 578.2 K (provisional,
  species/ammonia/boil-provisional.txt). The 4 set aside and rerun on 46 (-fourth, clear-once markers), then ammonia-boil4.
- CCl4 mono failed: phys-polar hit the 40-min timeout (510 functions). mono2 runs it with 4 h first; pairs2..plan2 follow on 65.
- Acetonitrile liquid 18/34. Acetaldehyde at full (41). Resumed 58, 60, 64 (3 idle, down from 7: reminder helps).
- Orchestrator 41.87 USD (+1.04). Holding 31 to ~13:20. Next wake 12:42.

## 2026-11-16 10:36 UTC, wake 19

- Merged. Seven workers (51, 54, 56, 58, 60, 64, 68) had ended their turns with jobs still running (Haiku stops after a long
  call). Resumed all seven; worker.sh now prints a keep-looping reminder each call (reaches them by self-update).
- Acetonitrile liquid 18/34. Ammonia 1080/0.0142-third still running (46). CCl4 mono running (65).
- Orchestrator 40.83 USD (+0.68). Holding 31 to ~13:20. Next wake 11:38.

## 2026-11-16 09:33 UTC, wake 18

- Merged. Haiku workers 51, 54, 58, 60 had ended their turns while their grows were still running ("continuing the loop"); resumed.
- Carbon tetrachloride grown -> mono..plan queued on 65. Acetonitrile liquid 15/34. Ammonia 1080/0.0142-third still running (46).
- Usage sample: orchestrator 40.15 USD (+2.8 since 09:02); 46 and 41 unchanged; no looper. Holding 31. Next wake 10:35.

## 2026-11-16 08:33 UTC, wake 17

- Fleet-wide stop ~07:33-07:42 WAS THE FIVE-HOUR LIMIT (owner: 17:39 AEST = 07:39 UTC); workers misreported it
  as a permission-check failure. Count cut 35 -> 31 (66, 67, 70, 71 archived, queues STOP). My 07:39 wake was cut by a restart; chain re-armed at 08:30 (next 09:32). All 32 workers resumed at 08:30.
- Next wake: check they all pushed again; resume any still idle. +2 only after two clean wakes; hold at 31 until a window passes clean.
- Acetonitrile liquid 11/34. Ammonia 1080/0.0142-third still unfinished (46).

## 2026-11-16 06:37 UTC, wake 16

- Formic acid law HOLE. Tally: HOLE 7 (methane, methanol, H2O2, water, ethylene, hydrazine, formic acid); not 2
  (ammonia, acetonitrile). Acetonitrile liquid 6/34 (43, 68, 69). Ammonia 1080/0.0142-third running on 46.
- All resumed workers pushing again. 49 -> benzene. Created 70 (pyridine), 71 (ethyl acetate). Count 35.
- Unstarted: glycerol, heptane, toluene, aniline, phenol, chlorobenzene, octane. Next wake 07:38.

## 2026-11-16 05:36 UTC, wake 15

- The per-machine call lock (wake 14) was held indefinitely; workers looped on "another call active" and 8 stopped
  (41, 46, 51, 56, 58, 60, 66, 67). Lock removed; all 8 resumed. Two workers committed by hand while stuck (41, 46);
  harmless (their own files). The per-state lock in liquid-state.sh stays.
- Laws: hydrazine HOLE (8 pairs). Acetonitrile NOT HOLE (K 7.74) -> 34 states on 43/68/69; queue acetonitrile-boil
  once >= 31 done. Ammonia: -third reruns on 46, then boil3.
- Created 68, 69 (acetonitrile states). Count 33. Next wake 06:34.

## 2026-11-16 04:30 UTC, wake 14

- ammonia-boil2 failed: the reruns of 271/0.712 and 1080/0.0142 again "done" with a 5-6 line log. Probable double
  run. liquid-state.sh: per-state flock + done only if log has "  snapshot in". worker.sh: per-machine call lock
  (fd 8 closed for children; tested). -third reruns + ammonia-boil3 on 46. If -third truncates again, the cause is
  phys-bulk itself at those states: then boil without them (path.py needs every done state complete) - owner note.
- Acetaldehyde grown (settled r11) -> mono..plan on 41. Others progressing.
- Created 66 (hexane), 67 (cyclohexane). Count 31 (cap lifted). Next wake 05:30.

## 2026-11-16 03:26 UTC, wake 13

- Second window passed, no limit stop; new window resets 08:10. Orchestrator cost 17.77 USD.
- 5-heavy grows (47 pentane, 59 cyclopentane, 60 butanol, 62, 63) show no push for 1-2 h: checked 59, it is
  RUNNING; round 0 of a large basis commits nothing until it ends. Not quiet.
- Progress: ethane pairs (17), methylamine pairs (56), dichloromethane pairs (50), formic acid full (49),
  acetonitrile full (43), hydrazine full (20). Ammonia 1080/0.0142 rerun still going (46).
- Created 64 (butanone), 65 (CCl4). Count 29 (max 30). Next wake 04:27.

## 2026-11-16 02:21 UTC, wake 12

- 45 stopped a third time: archived, ethanol on new 61. Policy now: retire a twice-resumed stopper to a NEW number.
- Grown: dichloromethane (-> mono..plan on 50), methylamine (56 continues), ethane (unsettled r7; 17 continues).
- Ammonia: 271/0.712 rerun done? 1080/0.0142 rerunning on 46, then boil2.
- Created 61 (ethanol), 62 (diethyl ether), 63 (THF). Count 27 of 30 max. Next wake 03:23 (window resets 03:10).

## 2026-11-16 01:19 UTC, wake 11

- Ethylene law HOLE (8 pairs, old queue). Five HOLE of six stage laws; ammonia the exception.
- ammonia-boil failed: two states truncated but marked done; rerun on 46 (…-again), then ammonia-boil2.
- 47 finished its states -> pentane grow. 18 and 21 finished their old queues. 20: hydrazine full.
- Formic acid at full (49), acetonitrile full (43). Created 59 (cyclopentane), 60 (butanol). Count 25. Next wake 02:19.

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
