# The fleet orchestrator's brief

You run the boiling-point data gathering for the 44 liquids of `fleet/species.tsv`, unattended,
until every species has gone as far through the pipeline as it can and the final table is
written. The owner is not watching. Read `runs/e8f/HANDOFF.md` once for what is established and
why; this file is what you do.

## What you are and are not

- You **never run a long job yourself**. Workers do. You edit queue files, read results,
  merge, and decide what each species does next.
- You **never message a worker to correct it**: a running worker reads no messages. You change
  `fleet/worker.sh` or a queue file on the main branch (workers fetch both at every step). The
  one message that works is the resume message below, to a worker that has stopped.
- You may add scripts and tools under `runs/e8f/`. You may **not** edit `docs/PLAY.md`,
  `docs/PHYSICS.md`, `src/chem/analyse.rs`, or change what existing engine code computes. Where
  HANDOFF.md marks an owner decision, take its stated default, write one line in
  `fleet/DECISIONS.md` saying you did, and carry on. Do not stop to ask.
- Do not open a pull request. Push only to your own branch (given at launch).
- If something is genuinely broken and you cannot repair it from here, write what and why at
  the top of `fleet/STATE.md`, publish the status page, and keep everything else moving.

## Every wake (you are woken hourly; a wake must finish in a few minutes)

1. `git pull` your branch. Read `fleet/STATE.md` (your own notes from last time) first.
2. Fetch the workers: `git fetch origin 'refs/heads/claude/e8f-w-*:refs/remotes/origin/claude/e8f-w-*'`.
3. For each worker branch, copy its new results into your tree (`git checkout origin/claude/e8f-w-N -- runs/e8f/species runs/e8f/fleet/done runs/e8f/fleet/failed`,
   one worker at a time; workers only ever add files, so later ones do not undo earlier ones).
4. Check each worker (see "Workers").
5. Advance each species (see "Pipeline"), writing new lines to queue files.
6. Read usage and set the child count (see "Usage").
7. Write `fleet/state.json` and `fleet/STATE.md`, run `python3 runs/e8f/fleet/status.py`, and
   publish `runs/e8f/fleet/status.html` with the Artifact tool, passing
   `url: https://claude.ai/artifact/QxtAopG5NkFZzHrCayyzSV` so the link stays the same.
8. Commit and push. End your turn. Do not wait or poll inside a wake.

## Workers

A worker is a Haiku session whose whole prompt is to run `bash runs/e8f/fleet/worker.sh N`
repeatedly. Create one with `create_session`: model `claude-haiku-5-5`, source the repository at
your branch, `outcome_branch: claude/e8f-w-N`, `permission_mode: auto`, tag `phys-fleet`, and
the prompt in `fleet/worker-prompt.txt` with N filled in. Set `MAIN` handling: the worker reads
its queue from the branch named in `worker.sh` (`MAIN`); if your branch is not that one, change
the default in `worker.sh` first.

- **Numbering:** never reuse a number for a different purpose. A worker's queue is
  `fleet/queue/N.txt`, one job a line (`<job id> <directory under runs/e8f> <command>`), ending
  `STOP` only when you want that worker to finish.
- **Stopped at the usage limit** (its session shows failed or idle and its queue is unfinished):
  send it exactly `Continue your loop: run bash runs/e8f/fleet/worker.sh N again and keep
  repeating it until its last line is ALL DONE.` That resumed all six in the test.
- **Quiet** (unfinished queue, no push for 45 minutes, and the resume message an hour ago did not
  help): create a replacement with the **same number**; it takes back what the old one pushed.
- **Failed jobs** (`fleet/failed/<job>` holds the log's tail): read it once. A job killed from
  outside is restarted by the worker itself and never appears here, so a failure is real. Requeue
  under a new id only if the cause is one you have removed; otherwise record it against the
  species and move on.
- A job id is never reused: to run something again, add a suffix.

## Pipeline (one species, in order; each stage is one job run in `species/NAME`)

```
../../fleet/stage.sh grow NAME      basis and shape, by the owner's stopping rule
../../fleet/stage.sh mono NAME      density, orbital energies, clouds, C6, polarisability
../../fleet/stage.sh pairs NAME     24 chosen pair geometries, frozen-orbital energies   (D2 data)
../../fleet/stage.sh full NAME      8 of them with the pair's own solve, both bracket ends (D1 data)
../../fleet/stage.sh law NAME       K and damping fitted to those; hole check
../../fleet/stage.sh plan NAME      liquid states chosen from the law itself            (D3)
../../fleet/liquid-state.sh NAME T DENSITY PS     one job per line of liquid/states.txt
../../fleet/stage.sh boil NAME      the boiling point
```

A stage is finished when `species/NAME/<stage>.ok` exists. Queue a stage only when the one
before it is finished. Liquid states are independent: spread them over all idle workers.

- **Size limits on these machines** (4 cores, 15 GB): `full` is for species of up to 3 heavy
  atoms (methanol took 11 minutes a pair and 9.4 GB). For 4 heavy atoms try 2 pairs first
  (`NFULL=2`) and continue only if a pair takes under an hour. Larger species stop after
  `pairs` and are recorded as "no pair solve on this hardware".
- **`grow`** may run for many hours on larger species and stops itself at 20 (`MAXH`). A species
  that stops at the time limit continues to `mono` with the basis it reached; say so in the table.
- **A law that begins HOLE** in `law.ok` does not get a liquid. Record it.
- **Order of work:** smallest species first (`species.tsv` is sorted that way). Keep every
  worker busy; never hold a worker idle waiting for a particular species.

## The three tests (HANDOFF.md section 4)

- **D1, the basis bracket:** the `full` stage gives both ends for every pair. Report per
  species the gap as a fraction of the binding. Also run `phys-shortfall` with `--more` on the
  methane and ammonia dimers as HANDOFF.md D1 describes (water is done).
- **D2, the repulsion constant against the orbital energy:** from each species'
  `frozen-NAME.txt` (exchange = column 9 minus column 5, valence overlap = column 8) fit K and
  compare with 5.435 / |highest orbital energy| from `mono.ok`. One table over all species.
- **D3, the planned liquid path:** water's 34 planned states are running or run in
  `species/water/liquid`. Its done-when is a boiling point within 3 K of 327.5 K. If it fails,
  say so in `STATE.md`, and still run the liquids (the data is wanted), flagging every boiling
  point as from an unvalidated plan.

## Usage

Start with **5** workers. Each wake:

1. `get_session` on yourself: read `external_metadata.rate_limit_info` (status, resetsAt).
2. Sum `external_metadata.usage.cost_usd` over yourself and every worker session; the difference
   from last wake (kept in `state.json`) is the hour's burn.
3. If the window closed since last wake (workers failed at the limit): record the total burned
   since the previous reset as the ceiling estimate, cut the child count by 10% (at least one),
   and send resume messages after `resetsAt`.
4. Otherwise, if the last two wakes had no limit stop: add 2 workers, or with a ceiling
   estimate, as many as keep the forecast (burn rate x time to reset) under 85% of it.
5. Never exceed 30. Trust `status` and `resetsAt` over your own forecast.

## Finishing

When no species can advance further, write `runs/e8f/RESULTS-fleet.md`: per species the stage
reached, K fitted and K predicted, the bracket gap, the law's errors, the boiling point against
`tests/liquids.rs`'s measured value, and what stopped it if it stopped. Put `STOP` on every
queue, publish the status page a last time, disable your hourly task, and stop.

## Where things stand at your launch (2026-11-15)

- Your branch is `claude/e8f-orchestrator`, made from `claude/laughing-keller-sagqde`. First
  act: in `fleet/worker.sh` change the default of `MAIN` to `claude/e8f-orchestrator`, commit,
  push. Every worker you create then reads its queue and the script from your branch.
- **Five workers are already running and count as your first five:** 17 (ethane), 18
  (ethylene), 19 (methylamine), 20 (hydrazine), 21 (hydrogen_peroxide). Each has a fixed queue
  of grow, mono, pairs, full, law on the *old* main branch, which you cannot edit; they run it
  to its STOP and finish. Merge their branches like any worker's. Number your own workers from 40.
- Workers 22 to 29 were stopped at launch to bring the count to five. Their species (methanol,
  propane, dimethyl_ether, ethanol, acetaldehyde, formic_acid, acetonitrile, ethylamine) are
  partly done on branches `claude/e8f-w-22` to `-29`: merge those, then queue each species'
  next unfinished stage on a worker of yours. Worker 30 holds water's last dilute state.
- The three checks before launch passed: the planned liquid path gives the one-constant water
  law 325.6 K against 327.5 K (D3's done-when is met); `law.py` ran end to end on ammonia; the
  test suite has the same failures as before this work (three the owner said to ignore, two
  that need an environment variable) plus `frames_stay_within_budget`, a wall-clock test that
  passed two runs of three on unchanged code.
- Ammonia's law from six pairs is poor at contact (6.5 kcal/mol in the reference, 77 in the
  law). Use `NFULL=12` for species of one or two heavy atoms so the fit has repulsive pairs.
- Water, methane and ammonia have not been through `mono`, `pairs` and `full` by the stage
  scripts on a worker (ammonia has, on the launching session's machine, and is committed).
  Queue water and methane like any other species.
