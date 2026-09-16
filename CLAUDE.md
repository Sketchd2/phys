# Working on phys

`README.md` says what this project is and why. This says how to work on it.

## Whose decision it is

**Do not make assumptions on the owner's behalf that deviate from planned
decisions. When a question comes up that the plan does not answer, stop and
ask.**

The plan is `docs/PLAY.md` and the decisions in it are settled — follow them
without relitigating. Routine judgement *inside* a planned item is not a
deviation: choose, say what you chose, and carry on. What needs asking is the
question the plan did not foresee, and those are recognisable:

- **A law the plan names but does not specify.** D3 says exchange moves a
  conserved quantity "at a rate set by a transport coefficient derived from
  their `Mixture`s". No thermal conductivity exists anywhere in the codebase.
  Picking one is a `PHYSICS.md`-weight decision wearing an implementation
  detail's clothes.
- **An item that turns out to be wider than it reads.** "`PaceMode::Fixed(1.0)`
  as what a world is" looks like a changed default. It is two mechanisms, and
  retiring the second one needs a scheduler rule D1 explicitly says does not
  exist yet. How far to go was not the plan's to decide silently.
- **A scope boundary the plan draws without saying which side something is
  on.** Contact is specified for things with materials. Whether a rock — which
  has none — should collide is a question about the play space, not about code.

The cost of asking is one round trip. The cost of not asking is a decision
buried in a commit, defended by a doc comment, and found six weeks later by
somebody reading it as settled.

## The axioms

These are not preferences. Code that breaks one of them is wrong even if it
passes.

**Only physical law is axiom-side.** Nothing in the engine knows what salt is,
what a forest is, or what a neutron does to water. It knows electronegativity,
bond order, Stefan–Boltzmann and the second law, and everything else is
*derived*. When you find yourself about to write `if is_forest` or a table of
melting points, stop: the answer is a law one level down.

**Measure, never be told.** An actor can lie; an instrument cannot. Properties
are measured from state — water availability is the liquid phase fraction of
whatever a node holds, not a `water: bool`. `Interaction::Author` is the single
path that sets a value by hand, and it is audited precisely because it is the
exception.

**Derived, with shortcuts stored.** Deriving from first principles every frame
is not the goal; deriving *once* and storing the result is. A shortcut may hold
a *rule* as well as a constant: how a kind of thing answers its conditions is
derived once for the kind and then run for each individual of it — one
derivation for how an oak meets light, water and wind, and a different one for a
birch. A stored rule stays a function of its conditions and never a frozen
outcome, so a drought still reaches every tree in the region it touches. A
tabulated constant that was never derived is the thing to avoid, and so is a
shortcut that no regional event can move.

**Detail exists where something is happening.** Most of the tree is deleted
every frame. Anything regenerable must regenerate bit-identically; anything
*touched* is pinned and persisted instead. This is the invariant the whole
design rests on — see `tests/consistency.rs`.

**No special cases.** A new scenario the engine was never told about should
work because the laws compose. Salt dissolving, a season ending, and a town
half-built are all the same handful of mechanisms seen from different angles.

## The model

A `Node` is a region of space at a `Tier`, holding a `Matter` — mass, momentum,
composition, temperature, radius. That is what a node *is*; everything else is
optional and regenerable.

Two independent LOD axes, and conflating them is the classic mistake:

- **Materialising** (`refine`/`coarsen`) gives a node its *own* contents — a
  `Vec<Body>` produced by `sample` from `(matter, spec, world_seed, path_key,
  epoch)`. Deterministic, so ~300 bytes *are* the bodies.
- **Promoting** (`promote`) turns one of those bodies into a `Node` of its own,
  one level down the tree. `children` runs parallel to `bodies`.

A **tier is a physics regime, not a tree level.** Many refinements happen within
one tier. `solvers::for_tier` maps tier to solver; `Tier::containing(metres)`
maps a size to a tier.

**Scale transform:** `sample` (coarse → fine) and `summarise` (fine → coarse).
The guarantee is `summarise(sample(m)) == m` on the conserved set, to
`IDEMPOTENT_TOLERANCE`. Tree-level verbs are `refine`, `coarsen`, `promote`,
`reparent`.

**Structures** are the other kind of node: a `morph::Program` plus a genome,
which *generates* geometry rather than sampling it statistically. `plant` seeds
one that grows; `emplace` states one that is already there. Programs are
`Tree`, `Coral`, `Tower`, `Wall`, `Terrain`, `Settlement` — and a program's
output bodies are the next level's nodes, which is how a moon becomes a town
becomes a building.

**Time:** one instant, every scale. A node is either *solved* or *coasted*.
Local rate is the product up the chain of `(1/γ) · √(1+2Φ/c²) · bubble`.
Trajectory runs on coordinate time; interiors run on local time.

## Commands

```sh
cargo test                 # ~2 min, must be green
cargo test --features postgres --test postgres   # needs $PHYS_PG; separate
cargo run --release --bin phys-demo       # the ladder, galaxy to nucleus
cargo run --release --bin phys-rehome     # re-parenting, frame to frame
cargo run --release --bin phys-bubble     # admin time dilation
cargo run --release --bin phys-headless   # render from a byte stream only
```

`[profile.test]` is optimised **with `overflow-checks` and `debug-assertions`
forced back on**. Do not drop those: `opt-level` above 0 turns overflow checks
off by default, and this project depends on integer overflow panicking in tests.
`tests/profile_guard.rs` asserts it.

`step_frame(wall_us)` takes a **wall-clock budget in microseconds**, not a span.
How much world time a frame covers is the *pace*. To simulate a season, use
`World::pace_fixed(seconds_per_frame)`.

## Traps

Things that have actually cost time. Read before diagnosing anything.

**Measure before asserting a cause.** This is the single most expensive failure
in this project's history. One bug got three confident wrong diagnoses in a row
— each plausible, each disproved by a five-minute probe. Write the probe first.

**A test that cannot fail proves nothing.** Verify a new test against the bug it
claims to catch. Two tests in this repo passed against the very defect they were
written for, and one went tautological when a refactor rewrote both sides of a
comparison.

**Identity is issued, not derived — `EntityId`, not `PathKey`.** `PathKey` is
the *address*: it derives children and seeds the sampler, and it changes when a
node moves. `EntityId` is the *name*: issued once, never reused, unchanged by a
move, a coarsen or a reload. `mixtures` and `environments` are keyed by name, so
`reparent` does not touch them. Use `identify` on a write path (it issues),
`identity` to read (it does not).

**Only an *event* may name a node** — chemistry set, an environment authored, a
node pinned, an actor interacting. Never a scheduler-driven path: which nodes
the frame budget advances depends on a wall-clock allowance, so naming a node
for its clock made identity depend on machine speed, and `next_entity` is
persisted. `clocks` and `histories` are therefore keyed by *address*, and
`reparent` migrates them along with the identity index.

**Three things still move in `reparent`:** `identities` (a node discarded and
rebuilt recovers its name from its address and nothing else), plus `clocks` and
`histories` for the reason above. Ledger and audit entries stay keyed by
`PathKey` deliberately — a measurement was made *of a place*.

**Tier is cached, never revisited.** Set at promotion from the body's radius.
`plant`, `emplace` and growth all change a node's size without updating it, so a
node can be two tiers from what its radius says. Known; in the backlog.

**Wire format encodes enum *positions*.** Renaming a variant is safe; reordering
or inserting silently reinterprets old saves. Append only. `FORMAT_VERSION` in
`wire.rs`; `tests/persistence.rs` catches size changes but not reorderings.

**Temperature means two things.** Thermodynamic at Planetary and finer; a
velocity dispersion at Stellar and Galactic. Blackbody laws apply only to the
first — `evolve_matter` gates on tier for exactly this reason.

**A hand-set temperature no longer holds.** `evolve_matter` radiates it away
within a few frames. Drive thermal scenarios with light.

## Naming

A name says what the thing *is*, in words a stranger would recognise, and no two
things share one. `docs/NAMING.md` records the decisions and the reasoning —
including which uses of "bulk" are correct (centre-of-mass motion, bulk density)
and which were retired (the coarse-vs-bodies sense, now "matter").

A rename is not finished when the identifiers compile: grep the *old word* in
prose too, and expect some hits to be a different sense that must be left alone.

## Where truth lives

- **`docs/BACKLOG.md`** — the live list. Every entry carries what was *measured*
  and a trigger saying when it becomes worth doing. Read it before proposing
  architecture; most obvious gaps are already there with numbers.
- `docs/DESIGN.md` — decisions and their reasoning.
- `docs/PHYSICS.md` — what is modelled and how.
- **`docs/PLAY.md`** — the plan for turning the engine into an inhabited world,
  and the decisions behind it. Read it before proposing anything about actors,
  surfaces, contact, creatures, construction or sessions; those arguments have
  been had and written down.
- `docs/PERFORMANCE.md` — the budget arithmetic.
- `docs/VIEWING.md` — the plan for showing a test rather than asserting it.
- `docs/NAMING.md`, `docs/GPU.md`.

Doc comments carry the *why*, at length, including what was tried and rejected.
They are load-bearing; keep them that way.

## Current frontier

The engine models what happens *inside* a node very well and what happens
*between* nodes barely at all. Closing that is the whole of the near-term work,
and **`docs/PLAY.md` is the plan** — the decisions are made, the order is set,
and the reasoning for each is recorded there. The four things everything else
waits on:

1. **The adjacency relation** — nothing knows which nodes are next to each
   other. This single gap blocks contact, fire spread, flooding, heat
   conduction, mass diffusion, friction, and debris landing on anything but its
   own parent. Four backlog entries are one missing primitive. **Start here.**
2. A promoted child never feels a force — `motion.velocity` is written only at
   promotion, so a promoted node is ballistic forever, and two promoted things
   cannot affect each other at all.
3. An issued identity, so an object keeps its name across a move — and so the
   side tables stop needing `reparent` to move them.
4. Nodes cannot split, so contents that legitimately expand are tracked by a
   node claiming a volume they have left. The spread measurement it needs is the
   same one surface handoff and detached fragments need.

`PLAY.md` D11 also finds the largest standing axiom violation in the codebase:
**`morph::Program` is a species table.** Six variants, fourteen dispatch sites,
and seven per-variant columns including a tabulated per-species decay rate.
Five of those columns are properties of the *material* or the *measured
environment* rather than of a species, and belong there. Do not add a seventh
variant — that is what D11 exists to prevent.

Two decisions in `PLAY.md` change things already written down, so do not treat
the older text as current where they disagree: **the world runs at one second
per second at every tier** (observer-following pace becomes a single-player
tool, and slow motion becomes replay of a recording), and **`PathKey` stops
being an identity**.
