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
cargo check --target wasm32-unknown-unknown --lib   # 2 s; nothing else builds it
```

`[profile.test]` is optimised **with `overflow-checks` and `debug-assertions`
forced back on**. Do not drop those: `opt-level` above 0 turns overflow checks
off by default, and this project depends on integer overflow panicking in tests.
`tests/profile_guard.rs` asserts it.

`step_frame(wall_us)` takes a **wall-clock budget in microseconds**, not a span.
How much world time a frame covers is the *pace*. To simulate a season, use
`World::pace_fixed(seconds_per_frame)`.

A world runs at **one second per second**: `World::new` calls `pace_realtime`,
which is `PaceMode::Fixed` at `budget.target_us` — one over the update rate, so
0.05 s per frame at 20 ups, not 1.0. A hand-set fixed pace of `1.0` is twenty
times real time, and the way that was caught is worth keeping: `resolution_floor`
reported 5.3 m for air where §3.4 says 0.27, exactly twenty times too coarse.

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

**The spec travels with the tier, and `retier` is what reconciles them.** A
tier is derived from a radius by `tier_for`, and `spec_for` comes with it —
because a node materialising under a policy meant for another scale is how eight
thousand molecules once ended up inside a node the size of an atom. One rule,
two callers: `promote` (a body becoming a node) and `retier` (a node whose size
changed), the latter called from `plant`, `emplace`, growth, damage, severing
and an authored radius. **Deliberately not from `coarsen`**, where a radius is
being restored rather than changed. Do not add a third path that sets `tier`.

**`coarsen` has an idempotent early return.** A node that has not been disturbed
takes it and never reaches the branch that rewrites the matter. Two tests passed
against a defect because a round trip through the early return is unchanged
either way; if you are testing `coarsen`, disturb the node first.

**A solver may cover a third of the span it was handed.** `hydro` and `md` both
substep to their own stability limit, capped at `MAX_SUBSTEPS`, and report what
they actually covered in `SolveReport::dt_used` — the cap is a budget, not a
licence. So a node's clock can advance more slowly than the frame while nothing
is wrong, and a test that counts *frames* rather than world time will drift.
`a_ball_loose_in_a_box` needs 4800 frames for what used to take 1500, for
exactly this reason.

**A whole module never gets compiled.** `src/wasm.rs` is
`#![cfg(target_arch = "wasm32")]`, so neither `cargo test` nor `cargo check`
touches it, and a signature change can leave it broken for as long as nobody
looks. It has been. The check is in the command list above and takes two
seconds.

**`forgettable` requires no promoted children.** A node with any is `unreachable`
rather than forgettable, so a test that asserts on the ensemble-crossing counter
by building a ladder measures nothing. Give it one small node with a floor above
its own radius instead.

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

**Phase 1 of `docs/PLAY.md` §7 is done, and §7 has been reordered.** The engine
used to model what happens *inside* a node very well and what happens *between*
nodes barely at all; that is what Phase 1 closed.

**Read `PLAY.md` §2A before anything else.** It records Issue 1 — *the engine has
no representation for the shape of a solid, at any scale* — and D13 to D16
follow from it. Two phases were inserted ahead of Ground as a result, so the
numbering below the insertion has moved: **Things**, **Crossings**, then Ground,
Water, Bodies, Minds, Making, Sessions. Nothing below Ground changed relative to
anything else.

**Phase 2 is Things** and nothing in it has started. Its done-when is a wooden
box that is *one node* with a recipe describing six walls, which responds as one
box, loses a single wall to a hard enough strike, and returns to ~100 bytes when
nobody is watching — plus a rock that no `Program` made bouncing off a boulder.

What landed, in the order it was built — each of these has its own commit with
the measurement in the message:

1. **D3, adjacency.** `src/neighbourhood.rs` is the one primitive: a spatial
   hash over a node's occupants (bodies *and* promoted children, in one index),
   `pairs(within)`, and two laws over a pair — `exchange` for a conserved
   quantity crossing a boundary, `contact` for an overlap resolving as an
   impulse. Both are exact two-body solutions rather than `rate × dt`, and both
   are symmetric to the bit from either side. `engine.rs` calls them as
   `exchange_within` (Planetary and finer) and `contact_within`.
2. **D4, the promoted child.** It feels the force its parent's solver computed,
   through the mailbox, instead of being ballistic from the moment it was
   promoted.
3. **Tier follows size** — `retier`, and the spec that travels with it. See the
   trap above.
4. **The spread measurement.** `Spread::of(parts)` — centre, rms, furthest,
   count — and `occupancy(radius)`, which is what a node splitting will need and
   what `worst_occupancy` already reports.
5. **One second per second.** `pace_realtime`, `PaceMode::Fixed` as what a world
   *is*, and the throttle that no longer applies in it. See the trap above.
6. **`G_EARTH` deleted.** Gravity is `Tree::gravity_at` — shell theorem over
   what a node is inside, with its own mass subtracted — cached on the node,
   persisted, and carried in the recipe blob so a client regenerates the same
   structure. Three consumers, all structural; see the leaf entry in the backlog
   for the fourth that does not exist yet.
7. **§3.3, dispatch reads state.** `Node::structural_mask` partitions a node's
   contents from its topology's joint radii, and the tier solver is handed the
   disordered remainder only. A building, a wolf and a boulder are all
   `Continuum` and none of them is a fluid. Hydro also substeps to its Courant
   limit now, which it never had.
8. **§3.7, the resolution floor is reported.** `resolution_floor(signal_speed)`
   and `resolution_floor_of(node)`, and a node crossed by its ensemble says so
   in `Stats::ensembled` instead of doing it quietly.

Phase 1's done-when list, all four, are tests:
`two_promoted_things_collide_and_rebound`,
`a_hot_node_beside_a_cold_one_equilibrates`, `a_branch_lands_on_the_next_tree`,
and `a_node_holds_ordered_and_disordered_contents_at_once`. Suite at the end of
Phase 1: **349 passed, 1 ignored** (`no_node_flings_its_bodies_out_of_itself`),
plus 6 Postgres, and five demos run.

### What Phase 2 will meet first

Left deliberately undone, each with a measurement and a trigger in
`docs/BACKLOG.md`. Read those entries before touching any of it:

- ~~**The ball-in-box test still runs at `Tier::Galactic`.**~~ **Done**, and it
  was not a tier swap. A `Node` is already what this engine means by a rigid
  body — one velocity, one spin, and `apply_contact` has always written to both
  — so what was missing was a *shape* for it to present. `src/shape.rs` and
  `Node::collision_shape` supply one, partitioned by §3.3's own
  `structural_mask`. The box is now a node with capsule walls and takes a ball's
  momentum as 48 tonnes rather than as one 500 kg panel. It found three things
  on the way; see its backlog entry, which is kept for them.
- **Derived gravity is in the parent's axes**, because nothing composes
  orientation anywhere in the tree. Terrain on a sphere is the scenario that
  makes it bite.
- **Exchange has a radiative coefficient and no conductive one.** D3 names the
  law and does not specify it; heat conduction through ground or water needs it,
  and picking a thermal conductivity is a `PHYSICS.md`-weight decision.
- **Only a built thing has a surface**, so only a built thing collides. This is
  a *provenance test standing in for a state measurement*, and §3.3 already made
  the same call correctly one layer down. The derivation the backlog recorded
  does not work: `sound_speed()` is the gas formula, not an elastic wave speed,
  and `density()` is bulk, so `E = rho c^2` comes out 58-82x low. Strength has
  no law at all. Deferred to D11.
- **Growth accumulates internal energy nothing sheds** — 231× thermal after
  forty years, reading back as 67,000 K while `temperature` says 291.
- **The sampler inflates anything bound by chemistry by 4.3×10⁵.**
- **Collision geometry is a sphere** — *mostly closed*. Both sides of a contact
  now carry a hull baked from spheres, so a member is the capsule its `base`,
  `tip` and joint radius always described. What is left is **flatness**: a row
  of capsules is not a plane, and making one needs either a rule for which
  members share a convex piece or a planar primitive. Neither is decided.
- **A node cannot split.** Phase 1 built the measurement it needs; the splitting
  itself is untouched.

`PLAY.md` D11 also finds the largest standing axiom violation in the codebase:
**`morph::Program` is a species table.** Six variants, fourteen dispatch sites,
and seven per-variant columns including a tabulated per-species decay rate.
Five of those columns are properties of the *material* or the *measured
environment* rather than of a species, and belong there. Do not add a seventh
variant — that is what D11 exists to prevent. Phase 2 moves the first five
columns off it, because a derived erosion rate cannot coexist with a tabulated
one.

One decision in `PLAY.md` still changes text written down elsewhere, so do not
treat the older text as current where they disagree: **`PathKey` stops being an
identity.** The other — one second per second at every tier — is no longer a
plan; it is what `World::new` does.
