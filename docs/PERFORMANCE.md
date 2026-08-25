# Performance

Target: **20 updates per second** on a Ryzen 5 3600 (6 cores / 12 threads) with
an RTX 2060 (1920 CUDA cores, 6 GB, 336 GB/s, ~6.5 TFLOP/s FP32).

Everything in §1 is **measured** on the reference implementation. Everything in
§3 is **projected**, and the projection's assumptions are stated so they can be
checked rather than believed.

## 1. Measured: single core, release build

`cargo run --release --example bench`

### Materialisation and restriction

| n | prolong | per body | restrict | per body |
|---:|---:|---:|---:|---:|
| 1,000 | 0.74 ms | 0.74 µs | — | — |
| 10,000 | 8.1 ms | 0.81 µs | 0.60 ms | 0.060 µs |
| 100,000 | 102 ms | 1.03 µs | 12.1 ms | 0.121 µs |
| 500,000 | 524 ms | 1.05 µs | 62.3 ms | 0.125 µs |

Flat in *n*, as it must be — the sampler is a single pass plus a fixed number of
projection passes. Restriction is 8× cheaper than materialisation, which is the
right asymmetry: the engine coarsens far more often than it refines.

### Gravity — Barnes-Hut, one leapfrog step (two tree builds, two traversals)

| n | θ=0.5 | per body |
|---:|---:|---:|
| 1,000 | 7.8 ms | 7.8 µs |
| 10,000 | 233 ms | 23.3 µs |
| 50,000 | 1.69 s | 33.9 µs |

### Cost of the optional terms, n = 50,000

| Configuration | per body | relative |
|---|---:|---:|
| monopole, θ = 0.7 | 13.8 µs | 0.41× |
| θ = 0.5 (baseline) | 34.0 µs | 1.00× |
| θ = 0.5 + quadrupole | 51.3 µs | 1.51× |
| θ = 0.5 + retardation | 55.8 µs | 1.64× |
| θ = 0.5 + retardation + 1PN | 62.7 µs | 1.84× |
| θ = 0.3 | 121.4 µs | 3.57× |

Retardation — the term that makes gravity causal, and the whole reason the
scheduler's lookahead is physically justified — costs **64%**. That is the
single most important number in this table, because it is the price of the
engine's central premise, and it is affordable.

### Short-range solvers, constant density

| n | SPH | per body | MD | per body |
|---:|---:|---:|---:|---:|
| 1,000 | 2.6 ms | 2.6 µs | 4.1 ms | 4.1 µs |
| 10,000 | 33.9 ms | 3.4 µs | 60.9 ms | 6.1 µs |
| 50,000 – 100,000 | 285 ms | 5.7 µs | 1.34 s | 13.4 µs |

### Structures: analysis and dynamics

| Structure | Members | Static analysis | One dynamic substep | Substeps per 50 ms frame |
|---|---:|---:|---:|---:|
| Tree, determinate | 500 | 14 µs | 1.1 ms | 45 |
| Tree, determinate | 2,000 | 52 µs | 5.0 ms | 10 |
| Tree, determinate | 8,000 | 342 µs | 24 ms | 2 |
| Tower, braced | 256 | 12 ms | 2.4 ms | 21 |

The determinate static column is the O(n) reverse pass and costs essentially
nothing. The braced tower's 12 ms is a conjugate-gradient solve, and it is the
single most expensive structural operation in the engine.

The dynamic column exists at all because of the tree preconditioner. A slender
chain of `n` beam elements has a condition number growing like `n⁴`, so Jacobi
preconditioning needed **3,645 iterations and 692 ms** for one substep of a
2,000-member tree. The load path is a tree, so the stiffness matrix has a
perfect elimination ordering with no fill-in; factorising it exactly costs one
6×6 inverse per joint and conjugate gradient then finishes in **one** iteration.

| 2,000-member tree, one substep | Iterations | Time |
|---|---:|---:|
| Jacobi | 3,645 | 692 ms |
| Tree factorisation | 1 | 5.0 ms |

Two things had to be right for that. The forest must be a *maximum* spanning
forest by stiffness, not a breadth-first one: a stay skipping five joints along
a chain is a shortcut, so breadth-first discovery reaches the far joint through
the stay and drops a chain element instead — severing the load path to keep a
member a hundred times softer, which on a 200-element chain with ten stays took
the iteration count from 1 to 1,124. Prim's algorithm takes the stiffest
available edge at every step and cannot make that mistake; the same case then
runs in 31.

And the factorisation is declined above a threshold on the structure's degree of
static indeterminacy. Applying it costs about an order of magnitude more than a
diagonal division, so it has to remove about that much from the iteration count.
On a moment frame, where every bay between two floors is a closed loop, it
removes a factor of two and loses on the exchange — measured as a tower's static
analysis going from 12 ms to 23 ms before the threshold was added.

### Memory

| Structure | Bytes |
|---|---:|
| `Body` | 184 |
| `Aggregate` | 232 |
| `Node` | 576 |
| `Snapshot` (history) | 80 |

5.4 M bodies per GB. On a 6 GB card with 60% given to bodies: **19.6 M bodies
resident** — the hard ceiling on the working set, independent of time.

`Body` is 184 bytes because it carries an 8-species composition. The GPU layout
in `docs/GPU.md` splits this into a 32-byte hot record plus cold arrays, which
raises the resident ceiling to ~110 M.

## 2. What the optimisation history cost and bought

These were found by benchmarking, not by inspection, and each was a real bug in
something that looked fine:

| Change | Before | After |
|---|---:|---:|
| Removed a needless per-particle sort in the neighbour search | 111 µs/body | 62 µs/body |
| Stopped allocating `Vec<Body>` inside the projection loop | 3.9 µs/body | 1.0 µs/body |
| Shrank the octree cell from 152 to 56 bytes | 58 µs/body | 34 µs/body |
| Allocated only occupied octants | 3.9 cells/body | 1.5 cells/body |

The neighbour-search sort existed to guarantee determinism against HashMap
iteration order — which the query loop never used, because it visits the 27
cells in a fixed order and each cell's contents are already in body-index order.
It was 60% of the molecular dynamics runtime, defending against nothing.

### Known remaining headroom in the CPU path

| Optimisation | Standard? | Expected |
|---|---|---|
| Bucket leaves (8–16 bodies) instead of single-body leaves | yes, universal in tree codes | 2–4× |
| Grouped traversal — one interaction list per nearby group | yes (Gadget, PKDGRAV) | 3–5× |
| Multi-threading across 6 cores | — | 5–6× |
| SoA layout + SIMD in the force loop | yes | 2–3× |

Taken together these would put the CPU path at roughly **0.3–1 µs/body/step**
for gravity, which is production-code territory. The reference implementation is
deliberately single-threaded and scalar so it can serve as the bit-exact oracle
the GPU path is checked against.

## 3. Projected: the RTX 2060 budget

### 3.1 What one frame is worth

At 20 UPS a frame is 50 ms. Allocating 70% to physics gives **35 ms**, with the
rest for rendering, input and slack.

### 3.2 Why the naive bandwidth estimate is wrong

The tempting calculation is: 336 GB/s × 35 ms = 11.8 GB per frame; at 48 bytes
per body over ten passes that is 24 M bodies. **This is wrong**, and it is worth
saying why, because it is the number an optimistic design document would quote.

Tree traversal is not a streaming workload. It is pointer-chasing with divergent
control flow: each body walks a different path, warps diverge at every opening
test, and the cells it touches are effectively random in memory. The limit is
latency and occupancy, not bandwidth. Bandwidth bounds the *resident set*
(§1, 19.6 M bodies); it does not bound the *step rate*.

### 3.3 The defensible number

Burtscher and Pingali's GPU Barnes-Hut — still the reference implementation for
irregular tree codes on GPUs — reports **~21× over an optimised serial CPU
implementation**. Taking an optimised serial baseline of 1–3 µs/body/step
(consistent with §1 plus the headroom in §2), the GPU path lands at roughly
**0.05–0.15 µs/body/step** for long-range gravity.

In 35 ms that is:

| Solver | Interactions/body | Projected bodies per frame at 20 UPS |
|---|---:|---:|
| Gravity (Barnes-Hut, θ=0.5, retarded) | ~1,000 | **0.2 – 1 M** |
| SPH / hydrodynamics | ~50 | **2 – 6 M** |
| Molecular dynamics | ~100 | **1 – 4 M** |
| Statistical (nuclear tier) | ~1 | **50 – 200 M** |

So the honest working-set budget is roughly **10⁵–10⁷ bodies per frame**,
depending on what they are doing — not the 2.4 × 10⁷ the bandwidth argument
suggests, and nowhere near the 1.2 × 10⁶⁶ a galaxy contains.

The cost model in `budget::cost` carries `GPU_SPEEDUP = 60`, at the conservative
end of that range and relative to *this* implementation rather than an optimised
one. It is a single constant in one place, and the runtime calibrates around it
from measured frames anyway (`FrameBudget::observe_frame`), so an error in it
costs one frame of over-commitment, not a broken engine.

### 3.4 What that budget buys

The demo's galaxy-to-nucleus descent materialises 132,400 bodies across 24
live nodes — 20 MB. That is a *complete* chain from a 9.2 × 10⁸ M☉ disc to a
2.7 fm nucleus, and it fits in 1–2% of the projected per-frame budget.

The budget is not spent on depth. It is spent on **breadth**: how many regions
can be resolved simultaneously. One deep zoom is nearly free; ten thousand
simultaneously-resolved star systems is what costs.

## 4. Measured behaviour of the frame loop

From the demo, single-core, cost model told the truth about its hardware:

```
frame  wall (ms)  accepted  deferred   bodies     sim step       debt
    0      36.81         8        16   132400  1.000e-21 s     3.68e8
    1       9.21         8        16   132400  1.000e-21 s     3.68e8
    ...
    7       9.33         8        16   132400  1.000e-21 s     3.68e8

cost-model calibration after 8 frames: 0.905x (learned from measurement)
```

The first frame pays materialisation costs the model has not yet calibrated for;
it then settles. Sixteen tasks are deferred every frame and reported as detail
debt. That is the design working: an observer demanding femtometre resolution
across a galaxy is asking for more than the frame can buy, and the engine says
so instead of missing its deadline.

`tests/budget.rs::frames_stay_within_budget` asserts this end to end with a
deliberately absurd observer (10⁻⁹ rad resolution, priority 1000).

## 5. Where time actually goes, and the granularity limit

One 8,000-body gravity node costs ~90 ms on one CPU core — larger than an entire
frame. The scheduler cannot subdivide a single node's step, so **no node may be
allowed to exceed a fraction of the frame budget**. This is a real constraint on
the refinement policy (`engine::default_spec`), not an implementation detail:
node sizes are chosen so that a single step fits comfortably, which on the GPU
path means 4,000–8,000 bodies per node.

The alternative — sub-stepping a node across frames — is possible but would
break the causal-ordering guarantee unless the node's clock were tracked
mid-step. Capping node size is simpler and costs nothing.
