# phys — a multiscale galaxy engine

A physics engine that presents a galaxy consistent down to the subatomic level,
in real time, on a consumer GPU.

## The problem, stated honestly

A galaxy of 10⁹ stars contains about **1.2 × 10⁶⁶ baryons**. An RTX 2060 can
step, at 20 updates per second, somewhere between **10⁵ and 10⁷ bodies** per
frame depending on the interaction range (derivation in
[docs/PERFORMANCE.md](docs/PERFORMANCE.md)). The gap is around **10⁶⁰**.

No optimisation closes a gap like that. Neither will hardware: 10⁶⁰ is forty-five
orders of magnitude beyond what Moore's law has left. Any engine that promises
"a galaxy down to the atom" and then talks about SIMD is not being serious about
the arithmetic.

So this engine does not simulate a galaxy. It is **indistinguishable from one**
to any observer inside it, at any resolution they can actually achieve — and it
makes that a precise, tested claim rather than a slogan.

## The four ideas

**1. Lazy materialisation with exact conservation.** Detail is generated on
demand from a seeded distribution and destroyed when nobody is looking. The
generator is constrained so that coarsening the generated detail returns the
original bulk state *exactly*: energy, momentum, angular momentum, charge,
baryon and lepton number. Measured worst case over 126 configurations spanning
60 orders of magnitude in mass: **5.8 × 10⁻¹⁶**, which is machine epsilon. No
experiment performed at the coarse scale can detect the deception, because
there is nothing there to detect.

**2. Causality as a scheduling primitive.** Nothing propagates faster than
light, so a region at distance *d* has a guaranteed lookahead of *d/c*.
Conservative parallel discrete-event simulation normally struggles to find a
lookahead at all; relativity hands us one for free, and it is enormous exactly
where the distances are large. Two clouds a kiloparsec apart can be stepped
independently for three thousand years of simulated time. Two nucleons 2 fm
apart must be stepped in lockstep. The same rule produces both.

**3. Observation as commitment.** An unobserved quantity has no value; a
measured one is recorded permanently in a ledger and never re-sampled. At the
subatomic tier this is not an approximation of quantum mechanics — it *is*
quantum mechanics, which is why the trick survives all the way down instead of
breaking at the bottom.

**4. A frame budget that spends detail, not time.** Each frame gets 50 ms and
solves a knapsack over the available work. Frame rate is the invariant;
fidelity is the free variable. A user who zooms too far does not get a
slideshow, they get a slightly coarser world and a visible detail-debt readout.

## Run it

```sh
cargo run --release --bin phys-demo     # guided tour, galaxy to nucleus
cargo run --release --example bench     # measured cost of every hot path
cargo test --release                    # 46 tests
```

The demo descends 23 levels from a 9.2 × 10⁸ M☉ disc to a 2.7 fm nucleus,
verifies that the detail regenerates bit-for-bit after being thrown away,
observes retarded light, commits a measurement to the ledger, fires an impulse
that takes 25,000 years to arrive, and runs a frame loop under budget.

## Documentation

| | |
|---|---|
| [docs/DESIGN.md](docs/DESIGN.md) | The architecture, and which parts are standard versus new |
| [docs/PHYSICS.md](docs/PHYSICS.md) | The scale ladder, the solvers, and their validation |
| [docs/PERFORMANCE.md](docs/PERFORMANCE.md) | Measured costs, the RTX 2060 / Ryzen 5 3600 budget |
| [docs/GPU.md](docs/GPU.md) | Mapping the CPU reference onto the GPU |

## Layout

```
src/
  math.rs      deterministic vector/tensor arithmetic
  units.rs     SI constants and the seven-tier scale ladder
  rng.rs       address-derived randomness — the basis of regeneration
  ids.rs       path keys: identity that survives being forgotten
  coords.rs    nested frames, relativity, retarded time, honest precision
  state.rs     aggregates, bodies, the conserved tuple, restriction
  prolong.rs   constraint-projected materialisation
  tree.rs      the scale tree; refine, promote, coarsen
  causal.rs    light cones, history rings, the conservative scheduler
  observe.rs   observers, instruments, the commitment ledger
  budget.rs    the frame knapsack
  engine.rs    the orchestrator and the interaction API
  solvers/     gravity, hydro, molecular dynamics, nuclear, quantum
tests/         consistency, causality, solvers, statistics, budget, interaction
```

## Status

This is a complete, tested architecture with a CPU reference implementation of
every subsystem. The GPU compute path described in `docs/GPU.md` is designed but
not written; the cost model carries an explicit, stated assumption for it rather
than pretending to have measured it. `docs/PERFORMANCE.md` is careful about
which numbers are measurements and which are projections.
