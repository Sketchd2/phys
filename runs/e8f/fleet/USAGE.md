# What uses the five-hour window (measured 2026-11-16 09:02 UTC)

Session metadata carries cumulative usage (cost_usd and tokens). Snapshots go in usage.json; a window's burn is the
difference between snapshots at its two ends, summed over every session, archived ones included.

| Session | Cost to date | Rate | Note |
|---|---|---|---|
| orchestrator (Opus) | 37.39 USD | ~3.5 USD/h since 03:24 (17.77 -> 37.39 in 5.6 h) | context 510k tokens; every tool call re-reads it |
| typical worker (Haiku; 17, 40-44) | 0.9-2.9 USD over 15-40 h | ~0.06-0.10 USD/h | ~130k context |
| worker 46 | 15.41 USD | | 296 M cache reads: the stuck call lock (04:29-05:33), a call every 2 s |
| worker 66 | 10.23 USD | | same loop |

Estimate of the 03:10-07:39 window that hit the limit: orchestrator ~15 USD, the stuck-lock loop ~20-25 USD
(46, 66 and probably others caught in it), the remaining ~30 workers ~3 USD. So the worker count was not what
used the window; the stuck lock and the orchestrator were. Ceiling estimate: ~40-50 USD-equivalent a window.

Forecast for a clean window now (lock removed, calls every 8.5 min): orchestrator ~15-17 USD + 31 workers x
~0.4 USD = ~12 USD -> ~28 USD, well under the ceiling.

Rules from this:
- Stuck loop check every wake: a worker whose cost rises more than ~1 USD between wakes (normal ~0.1) is looping;
  look at its last events and fix the cause (it was the lock) rather than cutting the count.
- The orchestrator is the largest single consumer: keep wakes short (few tool calls, little output).
- Child count from measurement: target = (0.85 x ceiling - orchestrator's window burn) / per-worker window burn,
  approached at most +4 a wake, once a window has passed clean at the current count (owner's limit search).
