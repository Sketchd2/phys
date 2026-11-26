# Fleet state (the orchestrator's own notes; newest first)

## (lifted 04:04 UTC 11-18: weekly limit reset; workers now report five_hour allowed) WEEKLY LIMIT WARNING (15:19 UTC 11-17): seven_day status allowed_warning, resets 2026-11-18 02:00 UTC. Hold at 38, add NO workers (replacements only) until it resets; if it hits, everything pauses until 02:00 and the backstop re-arms.

## WORKER COUNT: 37 from 19:05 UTC (blocks to 13:20 and 18:20 clean at 31, 35). Block ends 23:20 UTC. Every species is now started: add workers only to split real work (e.g. liquid states); if the limit hits: back to 31 and hold. The orchestrator's own context (409k) is the main cost: keep wakes minimal. Each wake: sample orchestrator + a few workers' cost (get_session) into usage.json; any worker rising >1 USD between wakes is looping: fix the cause. Full snapshot of every session once a block.

## WAKE CADENCE: every 62 min, self-scheduled. At the start of each wake call send_later (+62 min) with the same message; backstop task every 6 h re-arms if the chain breaks.

## 2026-11-18 13:39 UTC, wake 67

- CHLOROFORM LAW ON 110 WAS INVALID: it used a stale 2-pair full-chloroform.txt that 110 had itself committed at its
  part-4 start, and a worker keeps its OWN copy over main's. Set aside (law.*.2pairs-invalid); full reassembled HERE
  (8 pairs); 110 merge-skipped for chloroform/; law + plan on 132, which never touched chloroform.
  RULE: run law/plan on a worker that never wrote that species' files.
- Propane mono running (124); ethylamine parts 1 and 4 killed at 3 pairs (2,3,5-8 running).
- Acetone single pairs p4-1, p4-3, p5-1 and part 8 killed; isopropanol p3-1, p4-3 killed. Blocked species: no new work.
- Session check: 142 idle again: its queue was only isopropanol (blocked), so NOT replaced (STOP). 104 (THF) resumed
  (first). 91, 110 archived. Next 14:41.

## 2026-11-18 12:35 UTC, wake 66

- CHLOROFORM FULL ASSEMBLED HERE (4 parts x 2 = 8 pairs, full.ok). Law needs phys-lawenergy (not built here): law,
  plan on 110 (workers check out main's species dir at job start, so assembling here works).
- Propane grown (124): mono (12 h polar) + pairs on 124.
- Ethylamine whole pairs killed: 8 three-pair parts over idle 48, 53, 140, 141 (no new workers).
- BENZENE pairs killed (whole run): STOPPED for the owner, with CCl4, ethylene glycol, acetone, isopropanol, propanol.
  86 archived. Propanol parts 5-8 also killed (all 8 parts dead).
- Acetone single pairs p2-2, p2-3 FINISHED (p1-2, p1-3 killed): single pairs fit sometimes. Isopropanol p1-3, p4-1
  killed. Keep the queued ones running.
- Session check: 130 idle again -> fresh 145 (butanol); 142 resumed (first). 40 and 43 not on page 1 (running last
  wake). Next 13:37.

## 2026-11-18 11:31 UTC, wake 65

- BENZENE POLARISABILITY DONE (2nd 12 h try): mono.ok; pairs running on 86 (expect a memory kill).
- ACETONE / ISOPROPANOL / PROPANOL PAIRS BLOCKED ON MEMORY (STOPPED.txt in each): single pairs now killed too (a 46.5
  GB table rebuilt). Acetone 3-pair part 3 DID finish (on 133), so it is pair- or machine-dependent. Queued parts run
  out; no further splitting. Owner's call (bigger machine or smaller basis), with CCl4 and ethylene glycol.
- Propanol 3-pair parts 1-4 killed; 5-8 still queued.
- Session check: 138, 123 idle again -> fresh 143 (ethyl acetate), 144 (ethanol); 130 resumed (first). Others RUNNING.
  Next 12:33.

## 2026-11-18 10:27 UTC, wake 64

- 3-PAIR PARTS NO LONGER FIT for acetone and isopropanol: acetone parts 1,2,4,5 and isopropanol 2,3,5,8 killed (no
  part of either has finished). Those parts now run as single pairs (pair-one.sh) over each species' 4 workers;
  remaining 3-pair parts (acetone 3,6,7,8; isopropanol 6,7) still running. If they die too: singles. Join each part
  HERE with `pair-one.sh NAME K all`, then `pairs-part.sh NAME all 3`. Expect the same for propanol (139-141, 53).
- Acetic acid full part 1 killed even as ONE pair (28.7 GB table rebuilt): one retry on 139 (fresh machine). If it
  dies again, acetic acid full has 7 of 8 pairs: owner's call whether law may use 7.
- Session check: 135 idle again -> fresh 142 (its jobs); 91 (new jobs), 138, 123 idle (first) -> resumed. Next 11:29.

## 2026-11-18 09:23 UTC, wake 63

- Isopropanol parts 1 and 4 killed even at 3 pairs: single pairs (pair-one.sh, which now cuts its own sel-part file)
  over 40, 135, 136, 137; join each HERE with `pair-one.sh isopropanol K all` then `pairs-part.sh isopropanol all 3`.
  Killed partial frozen-part-1/4 removed and skipped.
- Propanol whole pairs killed: 8 three-pair parts over 53 + new 139-141. Count +3 (real work).
- Session check: 135 idle (first) -> resumed with new jobs; 53 resumed. All others RUNNING. Next 10:26.

## 2026-11-18 08:20 UTC, wake 62

- ACETALDEHYDE LAW HAS A HOLE (K 3.57, 8 pairs): no liquid. Tenth HOLE species. Its plan correctly refused.
- Propanol mono done (pairs next on 53). Acetic acid full parts 2, 3, 7 of 8 done; chloroform 1-3 of 4.
- Session check: 128 (ethyl acetate) idle again after its resume -> fresh 138. 41 and 109 finished their queues:
  archived. All others RUNNING. Next 09:22.

## 2026-11-18 07:17 UTC, wake 61

- Acetone, isopropanol mono done; their whole pairs runs were killed (47 GB table). Each split into 8 three-pair
  parts over 4 workers: acetone 91 + new 132-134, isopropanol 40 + new 135-137 (owner-only part skips set).
  Assemble HERE with `pairs-part.sh NAME all 3`; a part killed at 3 pairs -> pair-one.sh. Count +6 (real work).
- Session check: 128 (ethyl acetate) idle -> resumed (first). All others RUNNING. Next 08:19.

## 2026-11-18 06:13 UTC, wake 60

- Session check: every active worker RUNNING. No job finished this hour.
- Benzene polarisability on 86 restarted at 00:54 after 12.5 h with no exit status (killed from outside; the session
  reports exit 137 = SIGKILL, likely memory). Second attempt 5 h in. If it restarts again, stop benzene for the
  owner (like CCl4 / ethylene glycol): 840 functions may not fit phys-polar on a 15 GB worker. Next 07:16.

## 2026-11-18 05:07 UTC, wake 59

- Session check (one call, limit 100): every active worker RUNNING. No job finished this hour. Next 06:10.

## 2026-11-18 04:04 UTC, wake 58

- Session check, both pages: 127, 118, 113 stopped again within 30 min of the 03:00 resume -> fresh 129 (acetic acid
  full p4, p8), 130 (butanol), 131 (isobutane). Everything else RUNNING. No job finished this hour.
- Weekly limit has reset (five_hour allowed): hold lifted; add workers only to split real work. Next 05:06.

## 2026-11-18 03:00 UTC, wake 57

- SESSION CHECK (skipped at wakes 52-56, against the 20:41 rule): 127 (acetic acid full p4/p8) idle since 21:27,
  118 (butanol) since 21:49, 113 (isobutane) since 01:29: resumed. 121 (ethyl acetate) idle again after its resume
  -> fresh 128. All of 40-112, 120, 122-126 RUNNING. Do the session check EVERY wake (list_sessions children, 2 pages).
- No job finished this hour. Next 04:02.

## 2026-11-18 01:57 UTC, wake 56

- Quiet: no job finished since 00:53, no worker idle. Weekly limit resets 02:00: check rate_limit_info next wake
  and lift the hold if it reset. Next 02:59.

## 2026-11-18 00:53 UTC, wake 55

- Full parts done: acetic acid 2, 3 (of 8); chloroform 3 (of 4; p4 on 110). 111 idle -> archived.
  Benzene polarisability still running on 86. Weekly limit resets 02:00. Next 01:56.

## 2026-11-17 23:50 UTC, wake 54

- HEPTANE, TOLUENE, CYCLOHEXANE: grow.sh's 20 h limit fell inside round 0, so no grow state exists, yet grow.ok was
  written (heptane's mono panicked: "a growth state: NotFound"). grow.ok removed; regrow from scratch with MAXH=48 on
  65, 97, 101. Their old grow.ok is merge-skipped: REMOVE those three skip lines (merge-skip.txt, last grow.ok
  entries) once the regrow's grow.status says SETTLED or STOPPED with rounds, else the new grow.ok never arrives.
- Grown with a state (time limit, as planned): acetone (91), cyclopentane (59), isopropanol (40) -> mono + pairs.
- Chloroform full part 2 done (79): 79 archived (count 41 incl. 124-127). Parts 3, 4 on 111, 110.

## 2026-11-17 22:46 UTC, wake 53

- Grown: phenol (43), aniline (46), pyridine (50), diethyl ether (62), heptane (65): mono + pairs queued on each.
  29 of 44 grown. Ethylamine mono done (pairs next on 48). Next 23:49; benzene polar 12 h limit ~00:10.

## 2026-11-17 21:42 UTC, wake 52

- Grown: propanol (53), hexane (87: it was alive), chlorobenzene (89). mono (12 h polar) + pairs queued on each.
  24 of 44 grown. Others unchanged. Next 22:45.

## 2026-11-17 20:50 UTC (owner asked: full check of workers 1-127)

- Every session checked by status. Active 35 + 4 new (124-127): all RUNNING. Every retired worker is archived and
  every unfinished job of a retired worker has been moved (42..117 -> 118..127); 1-10 hold only pre-fleet test jobs
  for species already settled.
- Watch: 120 octane grow (log only "starting", 7 h), 87 hexane (no log change 6 h after "added 8"), 86 benzene
  polarisability (8 h; prints only at the end; 12 h limit ends ~00:10). Whether these are long silent steps or
  dead processes cannot be seen from here: if still silent at 00:10, replace.

## 2026-11-17 20:41 UTC (owner asked about propane)

- Propane grow: 42 stalled (16 h) -> 116 at 12:05; 116 carried on from its saved round 8, stopped 13:50, resumed
  15:19, restarted round 8 at 15:29 and went idle at 15:38. Missed for 5 h because wakes 47-51 read a quiet branch
  as a long job. -> fresh 124 (round 8 again; rounds 1-7 are kept in grow-propane.state).
- Session check found 3 more idle with unfinished queues: 117 pentane, 105 butanone, 115 (acetic acid full p4, p8)
  -> fresh 125, 126, 127; 121 resumed. Workers 40-91 not yet checked by status (next wake).
- IDLE CHECK FROM NOW: every wake, list_sessions(children) and resume/replace any worker whose session is IDLE with
  an unfinished queue; do not infer a long job from a quiet branch.

## 2026-11-17 20:38 UTC, wake 51

- Butane polarisability timed out at 4 h (884 functions; spill disk full): retry with 12 h, spill off, on 52 (as benzene
  on 86). Polar jobs for larger molecules get 12 h from now on (ethylamine's on 48 changed to match).
- Otherwise running: acetic acid full parts, chloroform full p2-4, DME full, 20 grows. Next 21:41.

## 2026-11-17 19:34 UTC, wake 50

- ETHANE BOILS AT 316.5 K (34 states; experiment 184.6 K). 11 of its 34 state logs had been laid over by stale
  truncated copies (108/109/110 committing other states' early logs): each restored from the worker that ran it,
  then the boil run HERE. All ethane liquid/ merges now skipped. Expect the same on any later species' states:
  restore by owner (fleet/done/<job> names the worker) before boiling.
- Ethylamine grown (48): mono, pairs on 48. Acetic acid full parts, chloroform full p2-4, DME full still running.

## 2026-11-17 18:30 UTC, wake 49

- Quiet hour: ethane 23 states; acetic acid full parts 1-4 running (93/94/109/115); chloroform full p2-4, DME full,
  benzene polarisability still running. Nothing to advance.
- Noted: octane's grow.log holds only "starting" lines across 102 and 120 (grow writes no state before its first step
  ends, so a worker whose machine is reclaimed starts over). If 120 stops again, keep the replacement alive rather
  than splitting: there is nothing to split.

## 2026-11-17 17:27 UTC, wake 48

- Acetic acid frozen-part-2.txt came back at 2 pairs AGAIN: 109 had committed a stale copy. Every worker that ran an
  acetic acid single pair (93 94 109 115) is now skipped for frozen-part-1/2/3, frozen-acetic_acid and pairs.ok.
  Reassembled HERE: 24 pairs. Full parts (1 pair each) on 93/94/109/115 still running.
- Ethane 23 complete states. 119 finished its ethane states: archived (count 37, weekly hold).
- Quiet workers are all on long single jobs (full parts, polarisability, fresh grows): no resumes this wake.

## 2026-11-17 16:23 UTC, wake 47

- ACETIC ACID PAIRS ASSEMBLED HERE: 24. (A merge had laid 115's killed 2-pair frozen-part-2.txt over the single-pair
  join: rejoined, and 115 is now skipped for part 2 and the assembled files.) Full as 8 single-pair parts
  (full-part.sh acetic_acid K 1) on 93/94/109/115, two each; assemble HERE with `full-part.sh acetic_acid all 1`,
  then law HERE, plan on a worker.
- Ethane 22 complete states. Chloroform full p2-4 and DME full still running (no pushes while they run).
- Resumed 93 94 109 115 (new jobs), 117, 120. Weekly-limit hold at 38 stands (resets 02:00).

## 2026-11-17 15:19 UTC, wake 46

- Block 09:20-14:20 clean (five_hour allowed). But the SEVEN-DAY limit shows allowed_warning (resets 02:00 UTC):
  hold at 38, replacements only. Orchestrator 139.79 USD total (+18 since 09:57), context 189k after compaction.
- Ethane 21 complete states. Acetic acid part 2 joined HERE from single pairs (3, each ~1 h on its own): pair-one.sh
  works where 3-pair parts are killed. Part 3 as single pairs on 93/94/109 too (115's 3-pair try still running; the
  first to finish wins: `pair-one.sh acetic_acid 3 all` or 115's frozen-part-3.done). Then pairs-part.sh ... all 3.
- 100 (ethyl acetate), 106 (glycerol), 107 (ethanol) stopped after a resume -> fresh 121, 122, 123. Resumed 105, 116.

## 2026-11-17 14:15 UTC, wake 45

- Ethane 18 complete states. Acetic acid part 2 KILLED even at 3 pairs (a 29.9 GB table): new fleet/pair-one.sh runs
  one pair alone; p2 pairs 1/2/3 on 93/94/109 (owner-only skips). Join HERE with `pair-one.sh acetic_acid 2 all`
  (in species/acetic_acid), then `pairs-part.sh acetic_acid all 3`. p3 on 115 may fail the same way.
- DME full on 112 has 1 pair in its log after 4 h (alive). Chloroform full parts 2-4 running.
- Resumed 41 46 48 86 104 117 (quiet >1 h; long full/polar jobs push nothing while running). Block 09:20-14:20 ends
  as this wake runs: judge at 15:18. Count 38.

## 2026-11-17 13:09 UTC, wake 44

- Ethane 13 complete states. CHLOROFORM FULL KILLED on memory after 2 of 8 pairs (like pairs): new fleet/full-part.sh
  (2-pair parts, spill off; `full-part.sh NAME all 2` assembles HERE). Part 1 = the 2 pairs already done; p2 on 79,
  p3 on 111, p4 on 110 (owner-only skips set). Then law HERE, plan on a worker. DME full still whole on 112 (3 h, no
  output yet): if it is killed, same parts (4 x 2).
- Butane grown: mono, pairs on 52. 102 (octane), 108 (ethane states), 114 (butanol) stopped again -> fresh 120, 119, 118.
  Resumed 52, 106, 107. Rate limit status allowed. Count 38. Block ends 14:20: judge at the 14:12 wake.

## 2026-11-17 12:03 UTC, wake 43

- Ethane 9 complete liquid states (21 logs of 34). Acetaldehyde pairs done (24): full..plan on 41. Acetic acid p1-3 on 115.
- Benzene polarisability killed on 86 (2-line log; its spill disk was full): retry queued on 86 (12 h, spill cleared and off).
- 42 (propane grow, 16 h quiet) and 47 (pentane, 12 h) -> fresh 116, 117; 17 idle (ALL DONE): archived. Count 38.
- Resumed 40 41 52 65 79 86 91 100 102 108 112 114 (quiet >1 h). Block ends 14:20: judge at the 14:09 wake. Next 13:06.

## 2026-11-17 10:59 UTC, wake 42

- Ethane 6/34 states. Acetic acid p4-8 done; p1-3 had failed on 57 9 h ago (old machine, no pairs-part.sh) ->
  fresh 115 (owner-only skips set); 57 archived. When 8/8: assemble HERE, then full..plan on a worker.
- Idle check skipped (context ~740k). Next 12:02.

## 2026-11-17 09:57 UTC, wake 41

- Block 04:20-09:20 PASSED CLEAN at 35-39. Orchestrator +26.84 USD (context 734k; auto-compaction at 784k is near:
  after it, this file is the whole memory). Holding 39.
- 95 (DME full), 98 (isobutane), 103 (butanol) stopped twice -> fresh 112, 113, 114. Resumed 100.
- Ethane 0/34 states yet; acetic acid 5/8 parts; chloroform full on 79. DME: when full.ok arrives, run law HERE
  (python3 ../../fleet/law.py dimethyl_ether cp via stage.sh law), then plan on a worker. Next 10:59.

## 2026-11-17 08:54 UTC, wake 40

- ETHANE planned: 34 liquid states, round-robin over 93, 94, 108, 109, 110, 111 (new). Boil HERE at >=31 complete
  (check every state's log has the tail line first).
- Chloroform pairs assembled here (24) -> full, law, plan on 79. DME full running on 95 (law here when full.ok).
- Acetic acid 5/8 parts. Count 39. Block ends 09:20; judge at 09:56.

## 2026-11-17 07:50 UTC, wake 39

- ETHANE LAW HAS NO HOLE (K 6.59, 8 pairs; computed by worker 17 on the OLD branch claude/laughing-keller-sagqde, so its law.py may predate main): third species with a liquid. Plan needs phys-b2 (not built here): on 93
  after acetic acid p7. Then queue its liquid states over several workers; boil HERE.
- DME: 95 assembled only 17 pairs: it (and 94) had committed stale partial part files before the commit fix. Each part
  is now taken only from the worker that ran it (merge-skip), restored, and assembled HERE: 24 pairs, pairs.ok.
  95 runs full; DME law and plan removed from 95 -> run law here (python) when full.ok arrives, plan on a worker.
  Same owner-only skips set for acetic acid's parts.
- Next 08:53.

## 2026-11-17 06:46 UTC, wake 38

- Chloroform p8 had failed on 79 for the same missing pairs-part.sh: requeued (p8b); worker.sh now refreshes scripts.
- 69 glycerol, 82 butanone, 84 ethanol, 99 THF (repeat stoppers) -> fresh 104 THF, 105 butanone, 106 glycerol,
  107 ethanol. 83 stopped (acetonitrile boil done; its last states not needed). Count 35.
- Acetic acid 5/8 parts. Resumed 95 (DME p2 + assemble + full..plan), 98, 103. Next 07:49.

## 2026-11-17 05:42 UTC, wake 37

- Worker 44 had no pairs-part.sh (an old machine keeps the fleet scripts it started with): DME p2/p3 and acetic acid
  p7/p8 failed there. worker.sh now refreshes runs/e8f/fleet/*.sh and *.py from main at every job start.
  DME p2, then assemble (pairs-part.sh dimethyl_ether all 3), full, law, plan, all on 95. 44 archived: count 36.
- Acetic acid parts 3/8 (p1-3 on 57, p7 on 93, p8 on 94). Chloroform 7/8 (p8 on 79): assemble here when done.
- Acetonitrile 1850/0.0148 (lowest density, 117k Ewald vectors) slow on 83; boil already done from 31. Next 06:45.

## 2026-11-17 04:38 UTC, wake 36

- Block 23:20-04:20 PASSED CLEAN at 37-40. Orchestrator +17.84 USD (context 647k; auto-compaction at 784k). Holding 37.
- Quiet: acetic acid parts 2/8, chloroform 7/8 (p8 on 79), DME 7/8 (p2 on 44). Fresh 97-103 running (99 stopped once).
- Resumed 69, 82, 99. Next 05:41.

## 2026-11-17 03:34 UTC, wake 35

- ETHYLENE GLYCOL STOPS AT PAIRS: parts 3 and 6 killed before their first pair (50.6 GB table): one pair does not fit
  a 15 GB worker. STOPPED.txt; owner's decision (with CCl4). Its other parts cancelled; 55 archived.
- Acetic acid part 4 done (parts fit): p1-3 on 57, p5-6 on 95, p7 on 93, p8 on 94. Chloroform 7/8 (p8 on 79). DME 7/8.
- Full instructions kept 82 and 84 going; 75, 76, 81, 85, 88, 90, 92 still stopped within ~10 min -> archived, jobs to
  FRESH workers 97 toluene, 98 isobutane, 99 THF, 100 ethyl acetate, 101 cyclohexane, 102 octane, 103 butanol.
  Rule now: a worker that stops within minutes of two resumes is archived and its job goes to a fresh session.
- Count 37. Block ends 04:20; judge at 04:37.

## 2026-11-17 02:30 UTC, wake 34

- Chloroform pairs killed after 21 of 24 (4.5 h, whole run) with NO spill attempt in its log: so a failed spill is not
  the only trigger; something builds up over many pairs (3-pair measurement saw none) or a heavy pair. Owner told.
  21 kept as parts 1-7; part 8 (pairs 22-24) on 79; then assemble here.
- 76, 81, 82, 84, 85, 88, 90, 92 stop ~10 min after each short resume. Sent the FULL worker instructions instead
  ("a long RUNNING is normal ... do not end your turn"). If they still stop, archive them and requeue their grows.
- EG and acetic acid parts started (0 done). DME 7/8 (p2 on 44). Next 03:33.

## 2026-11-17 01:26 UTC, wake 33

- Measurement finished: exit 0, 85 min, peak 11.62 GB, lowest available 3.55 GB, results identical to part 1. Done.
- Methylamine law HOLE (12 pairs; 9th HOLE). 78 (queue done) and 96 (measurement) archived: count 38.
- Ethylene glycol (50.6 GB table) and acetic acid pairs killed -> 8 parts of 3 each: EG on 55 (1-2), 93 (3-5), 94 (6-8);
  acetic acid on 57 (1-3), 95 (4-6), 44 (7-8, after DME p2). Assemble here with pairs-part.sh NAME all 3, then full..plan.
- Chloroform pairs running whole (916 functions, ~13 min a pair). Resumed 48, 76, 81, 82, 84, 85, 88, 90, 92. Next 02:29.

## 2026-11-17 00:22 UTC, wake 32

- Block 18:20-23:20 PASSED CLEAN at 37-40. Orchestrator +15.62 USD (context 537k). Holding 40 (no unstarted work).
- MEMORY MEASUREMENT (DME part 1, spill off; measure-part-1.mem): one pair at a time rises to 9-11.6 GB and falls back;
  peak 11.62 GB, 4.1 GB still available at the peak; NO growth from pair 1 to 3; spill 0, disk 29 GB free. So the
  killed whole runs were not a per-pair leak. They ran with spill ON, and each tried to write an 18-30 GB table onto
  ~29 GB of disk and failed; all 7 parts run with spill OFF survived. Likely trigger: the spill attempt (or a pair whose
  peak passes ~15 GB); not proven. Pair 3 still running; final line at the next wake.
- DME parts: all but 2 (on 44). Resumed 69, 76, 81, 82, 83, 85, 88, 90, 92. Next 01:25.

## 2026-11-16 23:19 UTC, wake 31

- Benzene grown -> mono (4 h polar)..plan on 86. DME parts done 1, 4, 6, 7, 8 (left 2 on 44, 3 on 95, 5 on 93).
- Measurement (96) first samples: 9.1 GB resident, peak 9.24, 6.0 GB available, no spill, 29 GB disk free, on pair 1.
- Resumed 75, 76, 81, 82, 84, 85, 88, 92. Block 18:20-23:20 ends now: judge at 00:22 (rate_limit status).

## 2026-11-16 23:01 UTC, owner request

- Owner asked whether memory-killed runs spill to disk: they try, but the 18-30 GB tables do not fit the ~19 GB free
  disk ("No space left"); the engine then rebuilds instead, and is still killed (DME after 4 pairs: growth across pairs).
  /tmp is on disk (ext4), not RAM. MEASURING: worker 96 runs fleet/measure-part.sh on DME part 1 (known results):
  species/dimethyl_ether/measure-part-1.mem has RSS/peak/available/spill/disk every 30 s, ~1.75 h. Report the peak and
  where it grows (per pair? during the rebuild?) to the owner; an engine fix is the owner's call.

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

### Wake 68 (14:41 UTC)
- Chloroform: law HOLE K 4.7059 over 8 pairs. Recorded as no liquid (the eleventh HOLE).
- Butane: polarisability timed out at 12 h. A 24 h retry is queued on 52.
- Ethylamine: pp2 is done. pp1, pp4 and pp5 were killed and have been split into single pairs on 48, 133 and 132. pp3, pp6, pp7 and pp8 are still running.
- Acetone: p2-1 is done. p4-2 and p5-2 were killed. Acetone stays blocked.
- Worker 143 (ethyl acetate grow) went idle for the first time and has been resumed.

### Wake 69 (15:44 UTC)
- Merge: the worker branches had not been fetched since the last wake, so the merge brought nothing in until I fetched them. merge.sh does not fetch.
- Ethylamine: pp7 is still running. pp3 was killed after 2 of its 3 pairs and has been split into single pairs on 140. The single pair p5-1 was killed too. A single pair cannot be split further, so if p5-1 stays killed, ethylamine will be one pair short.
- Propane: the polarisability step is printing results. The mono stage is still running on 124.
- Acetone p5-3 and isopropanol p5-1 and p5-2 were killed. Both species stay blocked.
- Usage limit hit at 15:13 (reset 15:20). The fleet is cut from 36 to 32: 40, 134, 136 and 137 are retired. 143 went idle twice, so it is replaced by 146 (ethyl acetate grow). The 21 workers that stopped at the limit have been resumed.

### Wake 70 (16:48 UTC)
- Propane: mono is done, and the pairs stage is next on 124.
- Ethylamine: pp6 and pp8 are done, and so are the single pairs p1-1, p1-2 and p5-2. The single pair p4-1 was killed on memory. So far two single pairs (p4-1 and p5-1) cannot run on a 15 GB machine. If ethylamine finishes 22 of 24 pairs, whether to go on is the owner's call.
- Workers 53 and 141 have nothing left in their queues, so they are archived. The fleet is down to 31.

### Wake 71 (17:51 UTC)
- Propane: the pairs stage was killed on memory. It is split into 3-pair parts pp1 to pp8 across 48, 124, 132 and 133, two parts each.
- Ethylamine: the single pairs p1-3, p4-2, p4-3 and p5-3 are done. Workers 132 and 133 had committed empty copies of other workers' pair files, and the merge laid those over the real results; parts 6 and 8 were hit the same way. Each part was restored from the worker that ran it, and merge-skip lines now keep every ethylamine part with its owner. I joined part 1 here (3 pairs). Done so far: parts 1, 2, 6 and 8 (12 pairs), plus 2 of 3 in each of parts 4 and 5 (p4-1 and p5-1 cannot run at 15 GB). Part 3 (single pairs) and part 7 are still on 140.
- Resumed 145 and 104, each idle for the first time since the limit resume.

### Wake 72 (18:54 UTC)
- Ethylamine: part 7 is done (3 pairs). Worker 48, which owns part 1, committed an empty copy of the joined part-1 file when its 3-pair job was killed, and the merge laid it over the join done here. It is restored, and every owner's joined-part file (frozen-part-K.) is now merge-skipped. The single-pair files stay with their owners. Done so far: parts 1, 2, 6, 7 and 8 (15 pairs) and 2 of 3 in each of parts 4 and 5. Part 3's single pairs are running on 140.
- Propane: pair parts 1 to 4 are running.
- 145 went idle a second time (18:01, ten minutes after its resume). The butanol grow moves to 147 under the two-strike rule.

### Wake 73 (19:58 UTC)
- Ethylamine: single pairs p3-1 and p3-2 are done; p3-3 is still on 140.
- Phenol and chlorobenzene: the polarisability run timed out at 12 h, as butane's did. Each gets a 24 h retry on its own worker (43, 89), with pairs after it.
- Propane: pp4 was killed on memory and is split into single pairs on 133.
- 104 went idle a second time (18:53, after the 17:51 resume). The tetrahydrofuran grow moves to 148.

### Wake 74 (21:02 UTC)
- Ethylamine: part 3 joined here. It is now 22 of 24 pairs: parts 1, 2, 3, 6, 7 and 8, plus 2 of 3 in each of parts 4 and 5. It waits for the owner's call on the 2 pairs that cannot run.
- Glycerol: grow is done. Mono (24 h) and pairs are queued on 122.
- Propane: pp1 and pp3 are done.
- 140 has nothing left in its queue and is archived. The fleet is down to 29.

### Wake 75 (22:05 UTC)
- Propane: pp2 is done (parts 1 to 3 are done). pp8 was killed on memory and is split into single pairs on 133, after p4's singles. pp5, pp6 and pp7 are still running.
- Cyclopentane: the polarisability run timed out at 12 h (the fourth such timeout). It gets a 24 h retry on 59.

### Wake 76 (23:09 UTC)
- Propane: pp7 is done. pp6 was killed on memory. The single pair p4-1 was killed too, so it is the first propane pair that cannot run at 15 GB.
  - Stale empty copies had replaced parts 1, 2, 3 and 7. Each was restored from its owner (48, 124, 132 and 132), and propane's merge is fenced like ethylamine's: each part file comes only from its owner, the split parts 4, 6 and 8 are joined here, and each single pair comes only from its assigned worker.
  - Part 6's singles go to 124 (p6-1, p6-2) and 132 (p6-3). Part 8's singles are spread: p8-1 on 133, p8-2 on 132, p8-3 on 124. pp5 is still on 48.
- Resumed 139 (acetic acid full p1 retry), idle for the first time since the limit resume. Acetic acid full parts 4, 5, 6 and 8 are still running on 129, 93 and 94.

### Wake 77 (00:13 UTC)
- Acetic acid: full part 5 is done (parts 2, 3, 5 and 7 now). Parts 1 (139), 4 and 8 (129) and 6 (94) are still running.
- Propane: pp5 is done (parts 1, 2, 3, 5 and 7, so 15 pairs) plus the single pairs p4-3 and p6-1. The single pairs p4-2 and p6-3 were killed on memory, as p4-1 was. Three propane pairs cannot run at 15 GB. Still pending: p6-2 and p8-1 to p8-3. If propane ends short of 24 pairs, whether to go on is the owner's call, as with ethylamine.
- 48 and 93 have nothing left in their queues and are archived. The fleet is down to 27.
- 139 went idle a second time (23:19, ten minutes after its resume). The acetic acid full p1 retry moves to 149. Resumed 147 (butanol grow), idle for the first time.

### Wake 78 (01:16 UTC)
- Ethanol and isobutane: grow is done. Mono (24 h) and pairs are queued on their own workers (144, 131).
- Acetic acid: full part 6 is done (parts 2, 3, 5, 6 and 7 now). Parts 4 and 8 are on 129, and part 1 is on 149.
- Propane: the single pair p8-1 is done. p6-2 and p8-2 were killed on memory, so five propane pairs now cannot run (p4-1, p4-2, p6-2, p6-3 and p8-2). That leaves 18 done and p8-3 pending on 124, so at most 19 of 24. Whether to go on is the owner's call.
- 94, 132 and 133 have nothing left in their queues and are archived. The fleet is down to 24.
- 147 went idle a second time (00:33, after the 00:14 resume). The butanol grow moves to 150.

### Wake 79 (02:20 UTC)
- Propane: the single pair p8-3 is done. Propane stops at 19 of 24 pairs and waits for the owner's call.
- 124 has nothing left in its queue and is archived. The fleet is down to 23.

### Wake 80 (03:23 UTC)
- No new results. All 23 workers are running.

### Wake 81 (04:27 UTC)
- No new results. Resumed 125 (pentane grow), idle for the first time since the limit resume.

### Wake 82 (05:31 UTC)
- No new results. Workers 122 to 150 are running. The older ones were running at 04:27.

### Wake 83 (06:35 UTC)
- No new results for a third hour. I checked the last push on every worker branch. 52 is running butane-mono-w52c (its "started" commit never reached the branch).
- Watch: 46 (aniline), 50 (pyridine) and 62 (diethyl ether) still report RUNNING their mono stage, started 15:56 yesterday. Each has a 12 h limit on phys-polar, and the rest of stage.sh mono has none. Nothing has been pushed past polar.log's header line. If they are still RUNNING at the next wake (about 16 h), look at them as possibly stuck.

### Wake 84 (07:39 UTC)
- No new results. 46, 50 and 62 are still on their mono stage at about 15.7 h, with nothing pushed past polar.log's header. worker.sh commits only changed files, so phys-polar itself is still running past its 12 h `timeout`. Precedent: butane-mono-w52b (12 h) started 20:51 and was reported failed at 14:11, about 17 h. So the timeout lands hours late, and these three should end around 09:00. If they have not failed or finished by the 08:42 wake plus one (about 09:45), they get replaced by fresh workers with the 24 h limit.
- 125 went idle a second time (06:11, after the 04:27 resume). The pentane grow moves to 151. Resumed 149 (acetic acid full p1), idle for the first time.

### Wake 85 (08:43 UTC)
- No new results. 46, 50 and 62 are at about 16.8 h on their monos. They get replaced at the next wake if they are still running.
- 149 went idle a second time (08:32, after the 07:40 resume). The acetic acid full p1 retry moves to 152.

### Wake 86 (09:46 UTC)
- Acetic acid: full part 4 is done (parts 2 to 7 now). Part 8 is on 129 and part 1 on 152.
- 46 (aniline), 50 (pyridine) and 62 (diethyl ether) were still RUNNING their mono at about 17.8 h, past the 12 h polar limit by more than butane's 17 h. They are replaced, as planned, by 153, 154 and 155: mono with the 24 h limit, then pairs. The old workers are stopped and archived.
- Worker 46's session had cost 17.66 USD (339 M cache-read tokens) by archive, against 0.9 to 1.3 USD for 50 and 62 over similar lifetimes. A long-lived session's context grows every call; replacing old workers is cheaper than keeping them.

### Wake 87 (10:49 UTC)
- No new results. 152 to 155 are RUNNING their jobs; their "started" commits have not reached their branches, as with 52 before, which is harmless. Workers 120 to 155 are running; the older ones were running at the last full check.

### Wake 88 (11:53 UTC)
- No new results. Workers 120 to 155 are running.

### Wake 89 (12:57 UTC)
- No new results. Full session check: every live worker from 59 to 155 is running (52 and 43 were running at the last check).

### Wake 90 (14:00 UTC)
- No new results. 151 (pentane grow) and 152 (acetic acid full p1) went idle and were resumed. Under the new rule in DECISIONS, long-job workers are resumed rather than replaced.

### Wake 91 (15:04 UTC)
- No new results. 151 and 152 pushed new "started" commits at 14:09 after their resume, so their jobs restarted (see the DECISIONS correction). Resuming saves nothing over replacing; idle time itself loses the job.
- 151 and 152 went idle again within 20 minutes of their resume (14:09, 14:17). Replaced by 156 (pentane grow) and 157 (acetic acid full p1).

### Wake 92 (16:09 UTC)
- No new results for the ninth hour. The newest workers (148 to 157) are running. Watch: butane-mono-w52c (24 h limit) has run about 25 h, so the 24 h limit also lands late.

### Wake 93 (17:12 UTC)
- No new results. Progress pushes since 16:57 from 65, 87, 97, 101, 126, 146, 148 and 150 (grows, plus hexane's mono). The fleet is working, on long jobs.

### Wake 94 (18:16 UTC)
- Butane: mono is done after about 27 h. MP2 mean polarisability is 53.4 bohr^3 (experiment about 55), and C6 is 1161 hartree bohr^6. The whole pairs stage (w52c) is already queued on 52. Propane's whole stage was killed on memory and butane is larger, so the 3-pair parts pp1 to pp8 are queued now on 52, 158, 159 and 160 (two each), with propane-style merge fencing so each part's files come only from its owner.

### Wake 95 (19:20 UTC)
- Acetic acid: full part 8 is done, so parts 2 to 8 are done (7 of 8). Part 1 is on 157. 129 has nothing left and is archived.
- Butane: the whole pairs stage (w52c) was killed on memory, as expected. pp2 (3 pairs) was killed too and is split into single pairs on 158. The singles are fenced to 158, and the part-2 joined file to nobody, since it is joined here.
- 155 (diethyl ether mono) went idle at 17:41 asking whether to commit runs/e8f/fleet/log/. It was resumed and told worker.sh does all committing (first idle). Its polar run will have restarted.

### Wake 96 (20:23 UTC)
- Dimethyl ether: full is done (8 pairs). Law and plan are queued on 126, which never wrote dimethyl ether files. 112 has nothing left and is archived.
- Butanone: grow is done. Mono (24 h) and pairs follow DME on 126.
- Octane: grow stopped at the 20 h limit inside round 0 with no grow state, as heptane, toluene and cyclohexane did (round 0 had about 26 h left). Its grow.ok is removed and merge-skipped from 120, and it regrows on 120 with MAXH=48.
- Butane: pp3 was killed on memory and is split into single pairs on 159, fenced as for pp2.
- 20:30: 157 (acetic acid full p1) went idle for the first time and was resumed. Its job restarts. 112 is archived.

## Wake 97 (21:27 UTC)
- DIMETHYL ETHER: law is a HOLE (K 4.3745, rms 0.548 kcal/mol over 8 pairs; B2 -inf at 300 and 600 K). No liquid. Its plan failing is the expected result of that.
- Hexane mono stopped at its old 12 h limit (queued before the 24 h rule) with no polarisability. Rerun on 87 with the 24 h limit (hexane-mono24-w87, then pairs).
- Ethyl acetate (146) and THF (148) grows are done. Mono (24 h) and pairs are queued on the same workers.
- Butane pairs parts 1 (52), 4 (160) and 6 (158) were killed on memory. Their pairs now run one at a time on the same workers, fenced as part 3's are.

## Wake 98 (22:31 UTC)
- Pyridine mono is done (MP2 mean polarisability 64.54 bohr^3). Pairs follow on 154.
- Butane pairs parts 5 (52) and 7 (159) were killed on memory. Their pairs now run one at a time on the same workers.
- 150 (butanol grow) went idle for the first time and was resumed. 157 went idle a second time (21:42, after the 20:30 resume), so it is archived and 161 (session_016vymdwokd2GcDNCGPT3p1s) takes acetic acid full part 1 (job acetic_acid-fullp1c-w161). Part 1 files are fenced from 157.
- 43, 52 and 59 are running. The fleet is 24 workers.

## Wake 99 (23:34 UTC)
- Pyridine pairs were killed as one run (an 81.6 GB table). Split into 8 three-pair parts: 1-3 on 154, 4-6 on new worker 162, 7-8 on new worker 163. Fenced per owner.
- Butane: single pairs p2-1 and p2-2 (158) were killed on memory (a 102.2 GB table at spill 0). Those two pairs cannot run on these machines; that is the owner's call, as for ethylamine and propane. Part 8 (160) was killed as a part and goes to single pairs on 160. The other queued singles still run.
- New workers 162 and 163 (pyridine parts). 159 went idle for the first time and was resumed. 155 went idle a second time (after its 17:41 resume), so it is archived and its diethyl ether mono (24 h) moves to new worker 164 (session_016vymdygFyqyAZ8PY9JGmCb). Diethyl ether files are fenced from 155. The fleet is 26 workers.

## Wake 100 (00:38 UTC)
- Ethanol mono is done (MP2 mean polarisability 33.71 bohr^3). Pairs follow on 144.
- Butane single pair p1-1 (52) was killed on memory, at a 102.8 GB table, like p2-1 and p2-2. All three singles killed so far hit about 102 GB. The remaining queued singles still run, because that is the only way to find out which butane pairs fit. A killed single is the owner's call.
- 161 (acetic acid full p1) went idle for the first time and was resumed. 159 and 150 went idle a second time, so both are archived. 165 (session_016vymdy7ahbzRd5yURUDs4p) takes butane singles p3-1..3 and p7-1..3. 166 (session_016vymdzj8A5Y1QTWd8YAKXA) takes the butanol grow with MAXH=48, resuming from its saved state (round 2). Each is fenced from its predecessor. The fleet stays at 26 workers.

## Wake 101 (01:42 UTC)
- Butane single p2-3 (158) was killed on memory, so all three pairs of butane part 2 are killed. That makes 4 of 24 butane pairs that cannot run here; the owner's call.
- Pyridine part 7 (163) was killed (an 80.8 GB table). Its pairs go to single pairs on 163, after part 8.
- 161 went idle a second time (00:57, 18 min after its resume), so it is archived and 167 (session_016vymdzkHqojhWroAoH9KgY) takes acetic acid full part 1 (acetic_acid-fullp1d-w167). This part has now restarted four times, on 157, 161 and 167. The fleet stays at 26 workers.

## Wake 102 (02:46 UTC)
- Butane singles p1-2 (52), p3-1 (165) and p4-1 (160) were killed on memory, at tables of about 102-103 GB. That is 7 of 24 butane pairs killed and none of its single pairs finished. The rest still run; the owner's call.
- Pyridine parts 1 (154), 4 (162) and 8 (163) were killed (81.6-81.8 GB tables). Their pairs go to single pairs on the same workers, after the parts still queued there.
- 164 (diethyl ether mono) went idle for the first time and was resumed. Its job restarts.
- 03:0x: the owner asked about running the out-of-memory pairs on their own Windows machine (large E: drive, GPU). This session cannot attach to it (no Remote Control folder offered). Added fleet/local-pair.ps1: one pair (default butane 2-1, the 102 GB table) with PHYS_SPILL_DIR on E:, timing and peak working set sampled every 30 s. Its result decides local against a rented high-memory machine; that choice is the owner's.

## Wake 103 (03:49 UTC)
- Butane singles p4-2 (160) and p6-1 (158) were killed on memory: 9 of 24 butane pairs killed, none finished.
- Pyridine part 5 (162) was killed. Its pairs go to single pairs on 162.
- The owner's local test (local-pair.ps1, butane 2-1, spill to E:) is running on their machine. No claude/e8f-local branch yet.

## Wake 104 (04:53 UTC)
- Butane singles p1-3 (52) and p3-2 (165) were killed on memory: 11 of 24 butane pairs killed, none finished. Part 1's three pairs are all killed.
- Pyridine single p7-1 (163) was killed: the first pyridine single pair, so pyridine's singles hit the same wall. The owner's call, with butane's.
- Pyridine part 2 (154) was killed as a part; its pairs go to single pairs on 154.
- No claude/e8f-local branch yet. The owner enabled Remote Control, but the phys folder on James_PC11 is listed as not available to cloud sessions, so this session still cannot read the local run.
- 164 went idle a second time (04:47, after its 02:46 resume), so it is archived and 168 (session_016vymdyge41BcXqmUn1XJdb) takes the diethyl ether mono (24 h). Diethyl ether files are fenced from 164. The fleet stays at 26 workers.

## Wake 105 (05:56 UTC)
- Butane singles p3-3 (165), p4-3 (160), p5-1 (52) and p6-2 (158) were killed on memory: 15 of 24 butane pairs killed, none finished. Parts 1 to 4 are wholly killed.
- Pyridine parts 3 (154) and 6 (162) were killed, so all 8 pyridine parts failed as parts. Their pairs go to single pairs on the same workers.
- No claude/e8f-local branch yet.

## Wake 106 (06:59 UTC)
- Pyridine singles p4-1 (162) and p7-2 (163) were killed on memory: 3 pyridine singles killed, none finished. Butane stands at 15 of 24 killed.
- No claude/e8f-local branch yet.

## Wake 107 (08:03 UTC)
- Butane singles p6-3 (158), p7-1 (165) and p8-1 (160) were killed on memory: 18 of 24 butane pairs killed, none finished. Parts 1-4 and 6 are wholly killed.
- Pyridine single p1-1 (154) was killed: 4 pyridine singles killed, none finished.
- 158 has nothing left and is archived. The fleet is 25 workers.
- No claude/e8f-local branch yet.
- 168 (diethyl ether mono) went idle for the first time and was resumed. Its job restarts.

## Wake 108 (09:07 UTC)
- Butane singles p5-2 (52) and p8-2 (160) were killed on memory: 20 of 24 butane pairs killed, none finished. Left: p5-3, p7-2, p7-3, p8-3.
- Pyridine singles p1-2 (154), p4-2 (162) and p7-3 (163) were killed: 7 pyridine singles killed, none finished. Part 7 is wholly killed.
- No claude/e8f-local branch yet.

## Wake 109 (10:11 UTC)
- Butane single p7-2 (165) was killed on memory: 21 of 24 butane pairs killed, none finished.
- Pyridine singles p1-3 (154) and p4-3 (162) were killed: 9 pyridine singles killed, none finished. Parts 1, 4 and 7 are wholly killed.
- No claude/e8f-local branch yet.
- 168 went idle a second time (09:13, after its 08:04 resume), so it is archived and 169 (session_016vymdzH1mszcuVNmaNqhSY) takes the diethyl ether mono (24 h). This mono has restarted on 155, 164, 168 and now 169. The fleet stays at 25 workers.

## Wake 110 (11:17 UTC)
- THE OWNER'S LOCAL TEST WORKED. Butane pair 2-1 (the 102 GB table that killed every fleet attempt) ran on the owner's Windows machine with PHYS_SPILL_DIR on E: and the cap at 300 GB: 6.60 h, peak working set 10.8 GB, spill 101.0 GB, one result line (overlap 1.9122e-3; electrostatic -3.447, exchange 11.533, HF beyond both -8.087 kcal/mol; 1768 functions). Built from a193b78b; nothing under src/ or Cargo.lock differs from the orchestrator branch, so it is the same calculation. Folded in as species/butane/frozen-part-2-1.txt and .done; local-2-1.* copied from claude/e8f-local. Fenced from 52 and 160.
- So the kills are not the table as such: with the table spilled, peak RAM is 10.8 GB, under a fleet worker's ~15 GB. What a run with spill off does instead (store half, rebuild the rest) is what is killed. Running the rest locally at 6.6 h a pair, or renting, is the owner's call.
- Butane singles p5-3 (52) and p8-3 (160) were killed on the fleet. Butane: 1 of 24 done (locally), 22 killed on the fleet, p7-3 still running on 165.
- Pyridine single p8-1 (163) was killed: 10 pyridine singles killed, none finished.
- ETHANOL PAIRS ARE DONE (24 pairs, closest contacts 1.60 to 4.43 A). Its full stage is split one pair a part: 1-3 on 144, 4-6 on 52, 7-8 on 160 (the queues that had emptied). Joined here when all 8 are done, then law and plan on a worker that never wrote ethanol files.

## Wake 111 (12:21 UTC)
- Butane single p7-3 (165) was killed on memory. Butane is finished on the fleet: 1 of 24 pairs done (2-1, on the owner's machine), 23 killed. What happens to the 23 is the owner's call (their machine at about 6.6 h a pair, or rented machines with enough disk to spill).
- Pyridine singles p2-1 (154) and p5-1 (162) were killed: 12 pyridine singles killed, none finished.
- 165 has nothing left and is archived. The fleet is 24 workers.

## Wake 112 (13:25 UTC)
- Pyridine singles p5-2 (162) and p8-2 (163) were killed on memory: 14 pyridine singles killed, none finished.

## Wake 113 (14:28 UTC)
- ACETIC ACID FULL IS DONE. Part 1 finished on 167 (after restarts on 157, 161 and 167), and all 8 parts are joined here: full.ok "8 pairs". Law and plan go to new worker 170, which never wrote acetic acid files. 167 has nothing left and is archived.
- Pyridine single p2-2 (154) was killed on memory: 15 pyridine singles killed, none finished.
- 170 is session_016vyme15UG6yXy1kCkTYMCc. The fleet stays at 24 workers.

## Wake 114 (15:32 UTC)
- ACETIC ACID: the law is a HOLE (K 4.3705, weighted rms 2.680 kcal/mol over 8 pairs; deepest pair -3.438 reference against -0.962 from the law; B2 -3.6e303 at 300 K). No liquid. Its plan failing is the expected result of that.
- Pyridine singles p2-3 (154), p5-3 (162) and p8-3 (163) were killed on memory: 18 of 24 pyridine pairs killed, none finished. Left: part 3 on 154 and part 6 on 162.
- 163 and 170 have nothing left and are archived. The fleet is 22 workers.

## Wake 115 (16:36 UTC)
- No new results or failures. Workers 131 to 169 are running; none idle.

## Wake 116 (17:39 UTC)
- Pyridine singles p3-1 (154) and p6-1 (162) were killed on memory: 20 of 24 pyridine pairs killed, none finished.

## Wake 117 (18:43 UTC)
- Ethanol full parts 1 (144) and 7 (160) are done: 2 of 8. Parts 2-3, 4-6 and 8 follow on 144, 52 and 160.

## Wake 118 (19:46 UTC)
- Ethanol full part 4 (52) is done: 3 of 8.
- Pyridine singles p3-2 (154) and p6-2 (162) were killed on memory: 22 of 24 pyridine pairs killed, none finished. Left: p3-3 and p6-3.

## Wake 119 (20:50 UTC)
- Pentane grow (156) stopped at the 20 h limit in round 2 with a grow state; accepted as the 2026-11-17 decision says. Mono (24 h) and pairs queued on 156.

## Wake 120 (2026-11-25 06:45 UTC) — restart after the weekly usage limit
- The weekly usage limit stopped the fleet and this session at about 21:55 UTC on Nov 20; it reset at 02:00 UTC on Nov 25. No wake ran in between; the backstop's 18 checks were all queued behind the stopped session. All running jobs died.
- Pyridine singles p3-3 (154) and p6-3 (162) had been killed on memory just before: all 24 pyridine pairs killed, none finished.
- Dropped and archived: 65 (heptane grow), 97 (toluene grow), 101 (cyclohexane grow), 120 (octane grow), 154 and 162 (empty). See DECISIONS.
- Resumed: 43 52 59 87 89 122 126 131 144 146 148 153 156 160 166 169 (16 workers).

## Wake 121 (07:42 UTC)
- All 16 resumed workers are running except 169, which went idle right after its resume (06:51) and was resumed again; the limit stop does not count as its first strike. Several workers skipped the "git pull" in the resume message as unexpected; harmless, worker.sh reads the queue it has.
- No new results yet (every job restarted at about 06:50).

## Wake 122 (08:46 UTC)
- No new results. 169 went idle again (07:51, 60 min after its second resume): archived, and 171 (session_016vyme1Pm6rk681Sj1sNa7a) takes the diethyl ether mono (24 h) and pairs. Diethyl ether files fenced from 169. The fleet stays at 16 workers.

## Wake 123 (09:50 UTC)
- No new results or failures. 171 has started; the workers checked are running.

## Wake 124 (10:53 UTC)
- No new results or failures. The workers checked are running.

## Wake 125 (11:57 UTC)
- No new results or failures. 166 (butanol grow) went idle at 10:59, the first since the restart, and was resumed; the grow restarts from its saved state.

## Wake 126 (13:00 UTC)
- No new results or failures. 166 went idle again (12:40, 41 min after its resume): archived, and 172 (session_016vymdxFLWSPeC9BqoATptS) takes the butanol grow (MAXH=48), resuming from its saved state (round 3). Butanol files fenced from 166. The fleet stays at 16 workers.

## Wake 127 (14:03 UTC)
- Ethanol full part 2 (144) is done: 4 of 8 (1, 2, 4, 7). Parts 3 (144), 5-6 (52) and 8 (160) remain.
- 171 (diethyl ether mono) went idle for the first time (13:33) and was resumed.

## Wake 128 (15:07 UTC)
- Ethanol full part 8 (160) is done: 5 of 8 (1, 2, 4, 7, 8). Part 6 moved from 52 (where it waited behind part 5) to 160 as ethanol-fullp6b-w160; part-6 files fenced from 52. Part 3 runs on 144, part 5 on 52.

## Wake 129 (16:12 UTC)
- Ethanol full part 5 (52) is done: 6 of 8. Parts 3 (144) and 6 (160) remain.
- 52 has nothing left and is archived. The fleet is 15 workers.

## Wake 130 (17:16 UTC)
- No new results or failures. 171 went idle a second time (16:39): archived, and 173 (session_016vymdzqKxfLE7DGFnKRiZg) takes the diethyl ether mono and pairs. That mono has now been on 155, 164, 168, 169, 171 and 173; it restarts each time. The fleet stays at 15 workers.

## Wake 131 (18:19 UTC)
- No new results or failures. 173 has started the diethyl ether mono; the workers checked are running.

## Wake 132 (19:22 UTC)
- No new results or failures. The workers checked are running. Ethanol full parts 3 (144) and 6 (160) are due around 21:00-22:00.

## Wake 133 (20:25 UTC)
- No new results or failures. The workers checked are running.

## Wake 134 (21:28 UTC)
- Ethanol full part 3 (144) is done: 7 of 8. Part 6 (160) remains; then join here, law and plan on a worker that never wrote ethanol files.
- 144 has nothing left and is archived. The fleet is 14 workers.
- Usage: worker 144 reported the seven-day limit at "allowed_warning" again (resets 2026-12-02 02:00 UTC), less than a day after the reset. If it is hit again the fleet stops until Dec 2. Fleet is not grown; finished workers are archived, not replaced.

## Wake 135 (22:31 UTC)
- No new results or failures. Ethanol full part 6 (160) still running. The workers checked are running.

## Wake 136 (23:35 UTC)
- ETHANOL FULL IS DONE. Part 6 finished on 160; all 8 parts joined here: full.ok "8 pairs". Law and plan go to new worker 174 (session_016vymdweQ46XzdvzdqaXraH), which never wrote ethanol files. 160 has nothing left and is archived. The fleet stays at 14 workers.

## Wake 137 (00:39 UTC, Nov 26)
- ETHANOL HAS A LIQUID. Law: K 9.2981, weighted rms 0.030 kcal/mol over 8 pairs (cp); B2 -97.77 cm3/mol at 300 K, +49.41 at 600 K. Plan: Boyle temperature 432 K, rho_ref 0.435 g/cm3, dipole 0.679 au, 34 states. The 34 liquid states are round-robin over 174 and new workers 175 (session_016vyme14roNGuagCnoeEcUD), 176 (session_016vymdwGZeSygi9BeuQHTDB) and 177 (session_016vymdz1tBYHn7mB65GJwNy), about 8-9 states each at about an hour a state (ethane's took 50 min). Boil HERE once at least 31 of 34 are complete (check every state's log has its last line).
- Isobutane mono is done (MP2 mean 53.28 bohr^3); pairs follow on 131. THF mono is done (MP2 mean 52.59 bohr^3); pairs follow on 148. Butane, an isomer of isobutane, needed 102 GB tables, so isobutane's pairs will probably be killed the same way.
- The fleet is 17 workers (3 added for the liquid states only; they are archived when their queues empty).

## Wake 138 (01:43 UTC)
- Ethanol liquid states: 4 of 34 done in under an hour each (432, 359, 302, 255 K).
- Isobutane and THF pairs killed at 102.6 and 102.8 GB tables; not split (see DECISIONS). They wait for a machine that can spill, with butane and pyridine. 131 and 148 archived. The fleet is 15 workers.

## Wake 139 (02:47 UTC)
- Ethanol liquid states: 8 of 34 done. The workers checked are running.

## Wake 140 (03:51 UTC)
- Ethanol liquid states: 12 of 34 done.
- 173 (diethyl ether mono) went idle for the first time (03:01) and was resumed; the mono restarts again.

### Wake 141, 2026-11-26 04:53 UTC
- Ethanol liquid states: 16 of 34 done (432 K at 0.226, 0.278, 0.331 and 0.383 merged this wake). 174-177 running.
- 173 (diethyl ether mono) went idle a second time at 04:52, so it was replaced: queue 173 STOP and fenced, archived; 178 created with the same two jobs (mono, then pairs).
- Why workers go idle: 173's transcript shows the worker writing a status sentence ("still RUNNING ... continuing the loop") and ending its turn (end_turn), usually in a turn opened by a task notification, instead of running worker.sh again. The calculation itself was not failing. worker-prompt.txt now says never to end a turn on a text reply; 178 got that line.
- Next wake trig_01SzQTGr269QTHFSfnBrq562 at 05:57.

### Wake 142, 2026-11-26 05:57 UTC
- The merge had undone seven finished ethanol states (108, 121, 186 and 95 K at 0.435; 432 K at 0.278, 0.331, 0.383): their "done" stayed but liquid.log came back mid-run. Cause: every worker branch carries stale in-progress copies of other workers' states, with a "done" from its source branch, and merge.sh took them in numeric order. Restored each from the worker that owns it; merge.sh now takes a liquid state only from the worker whose queue runs it. Re-merged: 16 of 34 done, every done log ends on its snapshot line.
- Ammonia, acetonitrile and ethane were boiled with every log tail checked, so they were not affected.
- All 16 workers running; no new done or failed jobs outside ethanol. Next wake trig_01GNupVe6ydMoesU7D4WxkkB at 07:00.

### Wake 143, 2026-11-26 07:01 UTC
- Ethanol liquid states: 18 of 34 done (432 K at 0.122 and 0.174 this wake), every done log ends on its snapshot line under the new merge fence. 174-177 running.
- Monos: butanone, ethyl acetate, glycerol and pentane reached their density runs. Phenol (43), cyclopentane (59) and chlorobenzene (89) are in 24 h polarisability runs whose logs have not changed, so their branches have nothing to push; their sessions are looping and report RUNNING. 178 (diethyl ether mono) is looping and running.
- Next wake trig_01DFLzhFWKiAHm3vRB9hQevN at 08:04.
