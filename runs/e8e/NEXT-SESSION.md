# Instructions for the next session (written 2026-10-20 ~17:30 AEST, before a weather shutdown)

Read `CLAUDE.md`, then `docs/SOLUTIONS.md` (§9.3 to §9.4 are the operational parts) and the E8e entry in
`docs/PLAY.md`. Branch: `claude/build-phase-6`. Last push: 1471e4d.

## State at shutdown
Nothing is running. Everything was stopped on purpose (expected power loss). Committed and pushed: all run data
and notes. **Not in git, on disk only:** `runs/e8e/RESTART.txt`, `handoff.sh`, `handoff2.sh` (they carry the queue
token the queue token (see RESTART.txt on disk), which must not be published), `gpu/tests/zz_f64.rs` (owner's scratch, leave alone), executables
(`*.exe`, ignored) and `runs/**/spill/` (ignored).

## First, bring the jobs back (in this order; `runs/e8e/RESTART.txt` has the exact commands)
1. Check the machine is sane: H:\ SSD present (`H:\phys-spill`, spill cap 40 GB), the executables still in their run
   dirs, nothing already running (`tasklist | grep phys-`).
2. **Growth** (`runs/e8e/grow/run.sh`): methanol resumes from `grow-methanol.state` (saved 14:47 on 20 Oct, round 14;
   the round in progress was lost). Ethanol and methane follow in the script. Ammonia is done (`done-ammonia`).
3. **Methanol queue server** in `runs/e8e/pairs-methanol`: `phys-queue.exe plan.txt --listen 0.0.0.0:7878 --lease 180
   --log queue`, with `PHYS_QUEUE_TOKEN=<token, see RESTART.txt on disk>`. It resumes from the pair files: 38 methanol-methanol pairs
   and 2 water-methanol pairs were done; 325 tasks in all. Check with
   `PHYS_QUEUE_TOKEN=<token, see RESTART.txt on disk> ./phys-queue.exe --status 127.0.0.1:7878`.
4. **GPU worker**: `nohup sh ./worker-loop.sh 127.0.0.1:7878 desk ./phys-worker-gpu.exe >> worker.out 2>&1 &` with
   `PHYS_QUEUE_TOKEN=<token, see RESTART.txt on disk> PHYS_SPILL_DIR=/h/phys-spill PHYS_SPILL_MAX_GB=40`. It needs 4 GB free; growth can
   squeeze it and the worker exits with `TooSmall` (exit 4). Restart it by hand when that happens (it happened once
   at 06:14 on 20 Oct). A pair takes ~23 min on the GPU; the queue needed ~11,000 min at that rate (~8 days).
5. **Watchdog**: exactly one instance, `nohup ./watch.sh &` in `runs/e8e` (jobs in `jobs.txt`; alerts in `alerts.log`,
   `ALERT` file when something is wrong, `ok.log` written every 10 min).

Do not run a second copy of anything: check `tasklist` first.

## Standing instructions from the owner (still in force)
- **Stall check-ins**: schedule a `send_later` about every 60 min while jobs in `runs/e8e/jobs.txt` are unfinished.
  Each one: read `alerts.log`, the `ALERT` file, `ok.log` age (over 15 min means restart the watchdog), report
  stalled/dead/finished jobs, restart a dead methanol worker, report whether any worker other than `desk` shows in
  `phys-queue --status`, then reschedule. Never wait without a bound (memory: no-unbounded-waits).
- **Keep CPU and GPU busy**: run independent pair work early; port to GPU and keep it if faster; GPU ports stay f64;
  speed-ups must not change numbers.
- **Don't trust earlier laws** (the methane one in particular); recorded in `docs/SOLUTIONS.md` S1.
- Don't wait for the owner on in-scope work. Do ask when the plan is silent (`CLAUDE.md`).
- Do not commit anything containing the queue token or other credentials. Check with a search before committing.

## Second machines
- A **Z440** was due to join the methanol queue to share the pair work; it was not available on 20 Oct. Steps are
  in `docs/CLUSTER.md` ("Methanol's MP2 pairs on more than one machine"). Server address on the LAN is
  192.168.0.22:7878, ZeroTier 10.33.0.26 (check the address is still right). Ask the owner whether it is ready.
- The owner also has an **HP ProLiant ML350p Gen 8** with one Xeon E5-2630 (6 cores, no GPU, no AVX2/FMA).
  Estimated ~5% extra throughput, only worth it if cheap. Before setting it up, check that the worker build does not
  assume AVX2/FMA, and ask for its RAM (16 GB or more keeps the 9.4 GB table in memory). The owner has not decided.

## Open owner decisions (do not decide for them)
- Stopping growth for methanol/ethanol/methane: the movement rule versus "nothing added".
- Pair budget and basis for methanol and ethanol (about 46 min/pair on CPU before the GPU worker).

## Work waiting
- **Ammonia cloud** was refit with both pair sets (`runs/e8e/cloud/ammonia/all/`): held-out 0.551 mean / 7.70 rms
  kcal/mol; close contacts below 2.8 A are the worst (rms 10.1); beyond 3.4 A rms is 0.9 or better; lone-pair
  population 0.323. Next: look at why contact is poor, then write the ammonia law (S1).
- Law-building chain still water-specific: `phys-diag` other modes, `phys-fit --bisector` and the [0,1,1] site types;
  the biased pair draw needs a first law. Wall-pair electrostatics for the other molecules.
- Then SOLUTIONS.md stages: cross terms and transferability (S2), chemical potentials (S3), onwards. S10 (calcium
  carbonate, rigid then flexible ions) comes after the E8c comparison.
- Possible GPU work: the Ewald part (21 ms against 8.8 ms flexible solve), the latency-bound pair loop.

## Things that went wrong, so they don't happen again
- Watchdog false alarms came from queue programs changing, idle servers and clock jumps; `watch.sh` handles them now.
- A contention test left jobs suspended for 53 min: bound every wait, resume what you suspend.
- The tool classifier blocks pushes of files holding the queue token; keep it out of anything committed.
- `pkill` does not exist here; kill by exact command line with PowerShell `Get-CimInstance Win32_Process`.
