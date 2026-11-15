# Orchestrator decisions (one line each; the owner may reverse any)

- 2026-11-15 Session metadata carries no `usage.cost_usd`; the usage rule runs on `rate_limit_info` (status, resetsAt) alone, burn recorded as unknown.
- 2026-11-15 First wake: no wake history, so no workers added (count stays at the five running, 17-21). Workers 40 and 41 are queued and are created at the next wake with no limit stop (two wakes, this one and that, without a stop).
- 2026-11-15 Ammonia's six-pair law moved to species/ammonia/six-pairs/; ammonia rerun with 24 selected pairs and NFULL=12 (brief: NFULL=12 for one or two heavy atoms). A new selection must not share the old frozen/full files, which append and skip by index.
- 2026-11-15 Methane's hand-made law, plan and empty liquid moved to species/methane/hand/ (no state of it had run); methane goes through the stages like any species.
- 2026-11-15 Water's D3 artefacts (law.ok, plan.ok, boil.ok 325.6 K, liquid/) stay where they are; water's grow/mono/pairs/full run in species/water first, and the hand law is moved to species/water/hand/ before water's own law stage is queued.
- 2026-11-15 Methanol's full stage rerun with NFULL=12 (2 of 8 pairs done by worker 22; the first eight indices are the same pairs, so they are kept).
