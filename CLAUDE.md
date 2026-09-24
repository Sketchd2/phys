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

**Nothing is deferred without the owner's say-so.** Finding an issue part-way
through an item is normal and expected; deciding what happens to it is not
yours. Say what you found, say what finishing it would take, and ask.

The default answer is that it goes **on the plan** — `docs/PLAY.md` §7,
scheduled into a phase — rather than into `docs/BACKLOG.md`, which is where a
*measured* gap waits for a stated trigger and not where an item goes to be
forgotten. Which of the two it is, is the owner's call, and so is doing it now
instead.

This covers the quiet forms, which are the ones that actually happen: a `TODO`,
a test marked `#[ignore]`, a case the code handles by doing nothing, a shortcut
with a comment explaining why it is temporary, a done-when reported met with a
caveat attached, and a commit message saying something was "left for later".
Each of those is a decision about what the project owes, taken without asking —
and by the time it is written down, the deferral has already happened.

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
one that grows; `emplace` states one that is already there; **`assemble` states
one that was made**, out of parts. Programs are `Tree`, `Coral`, `Tower`,
`Wall`, `Terrain`, `Settlement` — and a program's output bodies are the next
level's nodes, which is how a moon becomes a town becomes a building.

An **assembly** (`src/assembly.rs`) is the third case and does not have a
program of its own: it is a parts list the engine *generated* by assessing what
something is made of, and it decides the shape whenever it is present.
`Program` then records only provenance — what the material is. A composite is
one node with a recipe; a part becomes a node only when it detaches, and
`World::join` puts it back.

**A part is a `Body`.** There is one type for a thing inside a node, not one
for things that were sampled and another for things that were made. A `Body`
carries an `orientation`, `half`-extents (`Vec3::ZERO` meaning "a sphere of
`radius`") and the `substance` it is one of, which is everything a recipe needs
to state; an `Assembly` is `Vec<Body>` plus a parallel `Vec<Join>`, because a
join is a relationship between two things rather than a property of either.

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
move, a coarsen or a reload. `environments` is keyed by name, so `reparent` does
not touch it. Use `identify` on a write path (it issues), `identity` to read (it
does not).

**`mixtures` is gone.** D17 put the `Mixture` on `Matter`, so what a node is
made of travels with the node and needs no protection from a move at all — which
is what keying on a name was contriving to imitate. Speciation is now read with
`World::mixture_of` or straight off `node.matter.mixture`, and it is carried by
`promote` (down) and `coarsen` (up, blended by mass) like every other conserved
quantity.

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

**A node owns more than its radius, and a crossing measures three lengths.**
`Tree::domain` is the matter radius or the Hill radius about the parent,
whichever reaches further — without the second, anything standing on a planet is
ejected into interplanetary space, because a surface is exactly where a sphere's
boundary is. `Tree::contents_reach` is the third: `matter.radius` is the
equivalent uniform sphere, so a concentrated draw legitimately reaches past it,
and a node has only left when it is beyond *everything else its parent holds*.
The first implementation used the radius alone and flattened every ladder in
the engine in four frames.

**A crossing must not name a node.** Which nodes a frame advances depends on a
wall-clock allowance, so anything a crossing issues depends on how fast the
machine is. `World::cross` promotes a body to receive an arrival and
deliberately does not `identify` it.

**`settle` is not a general "refresh this node".** It re-summarises the whole
matter from the bodies, which is right on the way out and wrong in the middle:
a temperature somebody authored is replaced by whatever the sampled bodies
happen to carry. Measured, a 6000 K node and its 50 K neighbour both came back
at 3911 K. Where only the radius is meant, `Tree::remeasure_radius` is the one
to call — and only for a node whose detail is the authority (`pinned` or
`contains_edit`), because for a regenerable node the matter is the authority and
the bodies are a drawing of it. Letting the radius follow that drawing turned a
star into something else in four hundred frames.

**Every path that discards a body list must shed its children first.**
`Tree::shed_children`, and that includes `refine`: a reloaded node comes back
unmaterialised *with* promoted children, because the wire format writes
`children` and no bodies. Dropping them leaves a live node nothing can reach and
nothing will ever free.

**Wire format encodes enum *positions*.** Renaming a variant is safe; reordering
or inserting silently reinterprets old saves. Append only. `FORMAT_VERSION` in
`wire.rs`; `tests/persistence.rs` catches size changes but not reorderings.

**Temperature means two things.** Thermodynamic at Planetary and finer; a
velocity dispersion at Stellar and Galactic. Blackbody laws apply only to the
first — `evolve_matter` gates on tier for exactly this reason.

**A hand-set temperature no longer holds.** `evolve_matter` radiates it away
within a few frames. Drive thermal scenarios with light.

**`matter.radius` is the equivalent uniform sphere, not the bounding radius.**
`summarise` reports it from an rms and `sample` scales what it draws until it
comes back, with the factor `1/sqrt(3/5)` = 1.291 written down in both
`sampler::radius_scale` and `assembly::Assembly::extent`. Hand a recipe its
*bounding* radius and the sampler will helpfully scale the whole thing up until
the two numbers agree: a six-panel box's walls came out **1.87x too far out**,
silently, because every individual step was doing what it was told. If you are
writing anything that states a size, state that one, and use
`Assembly::bound` where you actually mean what fits inside a sphere.

**Collapsing is not forgetting.** `World::forgettable` asks whether the detail
may be *redrawn from the equilibrium ensemble*, and a broken thing must refuse
forever because the draw would mend it. `World::collapsible` asks whether the
detail may be *released*, and a broken thing accepts, because its recipe
reproduces the break exactly. One flag answering both is why no grown or built
thing ever coarsened — a thirty-year oak held 2,232 bodies for the life of the
world against 288 bytes that regenerate them. `Node::pinned` is this node's own
unregenerable detail; `Node::contains_edit` is a change below it that a fresh
draw would undo. `Tree::pin` no longer walks the ancestry.

**An assembled thing is the size somebody made it.** `sample_structured` skips
the fully-stressed optimiser when `morph.is_assembled()`. That sizing is how a
*grown* structure proportions itself, and running it on a box re-sizes every
seam until it can carry its load — which is the same thing as the box never
coming apart. If a joint is coming out a different size from the one the recipe
states, this is why.

**A sphere is not a cube.** `Body::half` is `Vec3::ZERO` for anything round,
and that zero is the discriminator — `(r,r,r)` would be a cube whose bounding
radius is `r*sqrt(3)`. `Body::solid` is the only thing that writes `half` and
`radius` together, so they cannot drift apart; do not set either by hand.

**A body's orientation is identity unless something stated it.** A recipe's
part, a promoted child coming back through `sync_from_child` — those know which
way a thing is facing. A sampled gas parcel does not, and giving it a random
draw would change nothing except every bit-exactness test in the suite.

**A part is promoted, not removed.** A wall that comes off a box keeps its slot
in the recipe: six parts, five joins and a break. The promoted child *is* that
slot (`children` runs parallel to `bodies`), so clearing the body list to
regenerate a five-part recipe orphans the node the break just made. Whether a
part is still attached is a join; whether it is its own object is the tree.

**Packing is measured against the volume a recipe states, not against a
bounding sphere.** A slab much wider than it is deep is nothing like its own
sphere, and `Formation::of_matter` divides by a volume: a patch of ground
presented rock at **60 kg/m³** until `material_of` and `derive_bonds` were made
to ask `Recipe::solid_volume`. That one fix moved three of Phase 2's box
measurements, which is recorded in its commit rather than quietly absorbed.

**A recipe is a pure function of itself and `(built, progress)` — so discard its
detail only when those move.** `grow_node` used to discard every structure's
bodies every frame on the reasoning that what it would generate had changed. For
a tree, which grows every frame anyway, that is free; for ground, which grows
never, it folded a planet's whole surface tree back on every frame it was looked
at — 42 promotions over six frames with two live nodes at the end.

**`f64::min` returns the *other* operand when one is NaN.** Terrain's embodied
energy is zero, so the growth law's `gross = usable / energy_density` was
`0.0/0.0`, and `net.min(reservoir)` then handed back the whole reservoir. A
patch of ground quietly ate itself every step. `Morphology::advance` now
dispatches on the habit, because what time does to a thing — grow, build, or
neither — is a property of how it is laid out.

**A piece of a sphere is not a sphere.** `gravity_at` modelled a face of a
planet as a ball at its own centre of mass, 1.9 million metres underground,
which put a ninth of the planet's mass pulling sideways: 8.96 m/s² tilted 2.6°.
The shell theorem is what a shell segment's field actually is. The same fix
found the second half of it — **a field derived from offsets is in root-aligned
axes**, so the rotation into a node's frame is the whole composition to the root
and not the part relative to one ancestor.

**A structure's parts are macroscopic objects, not molecules.**
`sample_structured` gave them the node's thermal energy as bulk motion —
`3/2 N k T` over the molecules rather than `3/2 n k T` over the parts. A 4733 kg
tree at 290 K jittered at 210 m/s and a new planet's six surface patches left at
60 km/s. The remainder is not lost: it belongs in the bodies' own internal
accounts, which is where `close_books` already puts it.

**History is only recorded for a node at `Residency::Causal` or better.** So
anything derived from a node's recorded past silently does nothing for a node
nobody is watching, and the result becomes a fact about the observer. Measured:
a melt cooling with no observer never recorded one moment in 600 frames. Where a
derivation needs the past, hang it off the *event* — `react_all` runs for every
node with a mixture whatever its residency — rather than off `histories`.

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

**Phases 1 to 4 of `docs/PLAY.md` §7 are done, and §7 has been reordered.**
The engine used to model what happens *inside* a node very well and what
happens *between* nodes barely at all; Phase 1 closed that. Phase 2 gave a
thing a shape and something to be made of, neither of which depends any more on
how it was made.

**Read `PLAY.md` §2A before anything else.** It records Issue 1 — *the engine has
no representation for the shape of a solid, at any scale* — and D13 to D18
follow from it, with D19 on how a change is stored. Two rules from those are
worth knowing before reading any of it:
**a surface primitive is always a filled solid, never hollow** (a box is six
slabs, and a hull over the whole box would enclose its own cavity), and **a
node's mixture is part of its matter**, which is what lets any node say what it
is made of and therefore whether it is solid at all. Two phases were inserted ahead of Ground as a result, so the
numbering below the insertion has moved: **Things**, **Crossings**, then Ground,
Water, Bodies, Minds, Making, Sessions. Nothing below Ground changed relative to
anything else.

**Phase 2 is done.** Its done-when was a wooden box that is *one node* with a
recipe describing six walls, which responds as one box, loses a wall to a hard
enough strike, and returns to a recipe when nobody is watching — plus a rock
that no `Program` made bouncing off a boulder. Measured:

```text
a box is one node        6 parts, 664 B of recipe persisted, 151.2 kg, one node
it responds as one box   six filled slabs, 8 spheres each, cavity empty
a part put in at 45 deg  its piece reaches 0.848 m along x, the square floor 0.600
struck at 20 m/s         utilisation 0.50, nothing comes off
struck at 700 m/s        utilisation 607, two panels become nodes
the seam decides         mortar 49.5 against cellulose 1.15 on the same box, 43x
nobody watching          1512 B of detail -> 664 B, regenerated at 0.0 m
a rock no Program made   rebounds at 0.0502 where a frame gives 0.0215
```

**The one decision in it that was the owner's**, because the plan did not
answer it: no `Program` variant describes a box and D11 forbids a seventh, so a
composite's recipe is **generated rather than selected**. A wooden box is a
tree, cut into pieces, carved and attached together; the engine assesses what
that produced and writes down a recipe for resampling it. `src/assembly.rs` is
that recipe, `Program` survives as *provenance* — what the material is, so a
box of oak planks weighs and burns like oak — and stops being the answer to
what shape a thing is. Do not read `Program` for geometry on a node that has an
assembly, and do not add a variant to describe an arrangement.

What landed in Phase 2, each with its measurement in its own commit:

1. **The two rigid-body defects** of §2A: `spin_rate` re-derived from
   `matter.spin` after a solve, and the collision shape composed with
   `motion.orientation`.
2. **`binding_energy` split** into `gravitational_binding` and
   `cohesive_binding`, before D17 rather than after.
3. **`Mixture` onto `Matter`** (D17). The `mixtures` side table is gone; read
   speciation with `World::mixture_of` or off `node.matter.mixture`.
4. **Material measured** (D14): `src/material.rs`, strength by Griffith with
   the flaw scale solved from the node's own cooling rate, and `rupture`
   deleted rather than kept as a reference.
5. **The recipe emits a surface** (D18): `shape::Surface` and `shape::Piece`,
   a union of filled convex solids with a material each, baked once and
   invalidated on `epoch`, and reconciled against the node's solid pools.
6. **Joining and breaking as one transform** (D15): `src/assembly.rs`,
   `Tree::assemble`, `World::join` and `World::detach`.
7. **A change is an edit** (D19): `Node::contains_edit`, and collapsing
   separated from forgetting — see the trap below.
8. **Collision runs against the surface**, at the level of detail the distance
   deserves.
9. **An equation of state outside its validity says so**, in `Stats`.
10. **A part and a body are one type** — `Body` gained an orientation,
    half-extents and a substance, `assembly::Part` was deleted, and `promote`
    hands a child the way its body was facing instead of `Quat::IDENTITY`.
    Added to the phase by the owner rather than deferred.
11. **A join is made of something, derived** — `Joint::bond` is a substance and
    `Topology::bonds` holds a material derived once per distinct substance.
    The same box with a mortar seam is loaded **43x** as hard as one with a
    cellulose seam.

Suite at the end of Phase 2: **400 passed, 1 ignored**
(`no_node_flings_its_bodies_out_of_itself`), plus 6 Postgres, and five demos
run, at `FORMAT_VERSION` 12.

**Phase 3 is done.** A boundary crossing is the event (D16). Its done-when was
the rocket — it leaves the forest, re-homes to the planet and then to the star,
and keeps correct gravity and correct neighbours throughout. Measured:

```text
leaves the forest      t = 0.150 s, 1.1 km up, re-homed to the planet
gravity there          9.815337 m/s2 against a closed form of 9.815335
reaches the star       t = 154002 s, 1.5088e9 m out — the planet's Hill radius
a limb shed by a tree  alive false, 3 live nodes not 4, energy moved 1e-16
a save taken mid-solve 1.52e-8 -> 4.3e-12 of the world's energy
a cloud in two minds   48.58 / 48.42, two nodes, energy moved 1.5e-16
a planet, 64 bodies    9 components, largest 86.5%, stays one node
the crossing pass      268 us for 8193 nodes, half a per cent of a frame
```

**The two decisions in it that were the owner's**, because the plan did not
answer them, are in `PLAY.md` §7 with their measurements: **what a node owns**
(its matter radius or its Hill radius about its parent, whichever reaches
further — the atmosphere is inverted for this and gives the large multiple to
the wrong body) and **the split's gate** (one component holding the bulk means
one region, applied twice so that dispersal is told apart from division).

What landed, each with its measurement in its own commit:

1. **A save summarises before it writes.** `Tree::settle` is `coarsen` without
   the destruction; `World::view` takes `&mut self` so the fix cannot be routed
   around.
2. **Two older defects the save's 2.35e-8 was hiding**, neither of which was
   the save: `sum_conserved` dropped `external_potential`, `cohesive_binding`
   and `chemical_energy` for every materialised node, and a stand-in `Body`
   could not carry what its child was holding (`Tree::unrepresented`).
3. **The crossing itself** — `World::cross_boundaries`, once a frame, after
   everything has moved.
4. **The orphan** — `Tree::shed_children`, on every path that discards a body
   list, `refine` included.
5. **The split and its inverse** — `Tree::components`, `resolve_extent`,
   `split_off`, `merge`, run over the nodes the frame advanced.

Suite at the end of Phase 3: **416 passed, 1 ignored**, plus 6 Postgres, and
five demos run. `FORMAT_VERSION` is 13 and `SCHEMA_VERSION` is 8.

**Phase 4 is closed, with three of its four done-when clauses met and the
fourth unmet by the owner's decision.** Ground, and one of D11's columns. A
planetary surface is a cubed-sphere parameterisation whose cells *are* the
node's children, so the surface tree is the scale tree with no adapter between
them, and a program is now a *generated rule* rather than a variant.

The clause that fails is **"`morph::Program` carries no column that is a
property of a material or a measured environment"**. `substrate` went; the
other four each want a law the engine does not have, and `design_flow` — a
table of a fluid density and a gust speed, which is precisely a measured
environment — is the one that is merely unwired. `PLAY.md` §7 records the
verdict clause by clause and what was left open with it. Do not read the
measurements below as the whole done-when met.

Measured:

```text
a vector between frames    9.8200 m/s2 along the node's own -z, oriented
a test that is watched     256 bodies cover 0.02-0.90 of the frame; 0.000 before
flatness, generated        a 2.172 m wall's rest height varies 0.0000 m against
                           0.2991 m as capsules; a patch's columns cover 1.0000
                           of their cell against 0.2848
Griffith on the worst flaw bedrock 0.91 of the retired table, not 77x low
a program is generated     five habits, stated size 1.000x drawn; before, out by
                           1.44 / 1.53 / 0.89 / 0.76 / 0.95
ground                     a planet's surface is 296 B of rule and 0 B of detail;
                           orbit to 0.55 m in 9 levels; 9.7055 m/s2 at 0.092 deg
                           off its own down; 10 km of walking crosses 41 patches
                           and stays within 2.1 m of 6.4098e6 m; look away and
                           back, worst displacement 0.0e0 m
a planet from origin       6.000e24 kg and 2.0000e7 m fell in for 1.000e5 s and
                           became a planet at frame 86: 6.7810e6 m, 1743 kg/m^3,
                           4.564 m/s^2 — and the run is filmable
a layout derived           molten silicate froze through 1845 K and laid down
                           1.662e-2 m grains; 0.25 m -> 4.53e-3, 4 m -> 6.09e-2;
                           a rock never molten gets no layout
forgetting is a decay      a 0.10 m mark, half-life 0.0116 s, gone; a channel
                           that took 34.6 kg still 0.3000 m deep a year later
one column retired         `substrate`, measured off the node's own feedstock
```

**Six decisions in it were the owner's**, because the plan did not answer them,
and they are in `PLAY.md` §7 with the rest: the terrain morphology is a
*generated program* rather than a variant; the flaw population is intragranular
and bounded below by the critical nucleus; the erosion done-when replaces the
tidal one, which moved to Phase 5; a planet generated from origin is a done-when
and has to be *watchable*; a floor plate goes in the general solver with a
generated program rather than a `Tower` type; and **there are no hardcoded
plans** — an actor designs and places outputs, the engine assesses what that
produced and writes a program for it, and a program is a function of the current
and historical data derived from the axioms rather than engineered.

What landed, each with its measurement in its own commit:

1. **A vector carried between frames is turned into them** — `Tree::axes_from`,
   `Tree::into_axes_of`, and `gravity_at` rotating by the composition **to the
   root** rather than relative to one ancestor.
2. **A test can be watched** — `src/render.rs`'s `Paint`/`Shot`/`Film` and
   `src/film.rs`, which gathers a scene without disturbing it. `PHYS_FILM=<dir>`
   or nothing happens at all.
3. **A generator that lays down a flat thing states it** — `Skeleton::push_box`
   and `push_plate`, so `Wall`, `Tower` and `Terrain` emit slabs the way an
   assembly does.
4. **Griffith on the worst flaw** — `Material::worst_flaw`, a weakest-link draw
   between the critical nucleus and the grain.
5. **A program is generated, not selected** — `src/recipe.rs`. `Habit` and
   `Recipe::generate`; `Morphology::recipe` replaced `assembly`.
6. **Ground** — `Recipe::Tiled`, the cubed sphere, `World::approach`.
7. **A planet from origin** — `World::assess_surface`, gated on random loose
   packing and then on central pressure against strength.
8. **A layout derived from what the node went through** — `Recipe::Granular`,
   written by the *freezing event* rather than by a state test.
9. **Forgetting is the deviation reaching zero** — `src/erode.rs`,
   `Morphology::field`, `Interaction::Mark`, `World::weather`.
10. **One of D11's five columns retired, four measured** — `substrate` is gone;
    `density`, `energy_density`, `material` and `maintenance` each hit a law the
    engine does not have, with numbers.

Suite at the end of Phase 4: **442 passed, 1 ignored**, plus 6 Postgres, and
five demos run. `FORMAT_VERSION` is 18.

**Three things were left open at the close rather than fixed**, and they are in
`PLAY.md` §7 and `BACKLOG.md` with their measurements: `tests/accretion.rs` runs
at a starved frame budget so the planet it derives differs run to run; the same
test prints a size and a density taken from two different radii; and a body
leaves the sampler with a temperature and an internal energy that disagree by
whatever the binding term is — 1739.50 K while holding 2430.79 K of energy on a
planet, which the first save turns into a 40% step.

### What Phase 5 will meet first

Left deliberately undone, each with a measurement and a trigger in
`docs/BACKLOG.md`. Read those entries before touching any of it:

- **An uncemented aggregate has no derived strength**, so nothing can erode
  sand. Every solid `Material::measured` returns is a bonded one: a poured pile
  of dry silica comes out at 3.26e7 Pa in tension against the 4.5e3 Pa water at
  3 m/s can press with. The decay machinery, the conservation criterion and the
  promotion rule are built and tested; this is the number that goes into them,
  and it is `PHYSICS.md`-weight.
- **`morph::Program`'s last four columns** each want a law the engine does not
  have, and each was measured rather than assumed — see the entry, and see
  `Program::maintenance`'s own doc comment for the forty-four orders.
- **Exchange has a radiative coefficient and no conductive one.** D3 names the
  law and does not specify it; Phase 5's free surface is its trigger and D17
  unblocks the choice.
- **Growth accumulates internal energy nothing sheds** — 231× thermal after
  forty years, reading back as 67,000 K while `temperature` says 291.
- **A detached fragment is neither a body nor a node.** Both of its blockers are
  gone, so what is left is the line of policy its entry describes.
- **Sideways handoff across a patch edge** has never been exercised in anger.
  The mechanism is Phase 3's and the tiling is Phase 4's; what is missing is a
  scenario, and the first one is Phase 6's quadruped.

**Scheduled rather than left**, and `PLAY.md` §7 is where they live — none of
them is a backlog entry to be picked up on a whim:

- **A 52/45 separation stays one node.** Phase 3's split gate refuses anything
  where one component holds the majority, and the boundary sits inside the
  working range. The alternative measured cleanly and was declined because it is
  a ratio the plan does not state; revisit when something built on the split
  needs the uneven case.
- **A planet stops at random loose packing**, 1743 kg/m^3 against Earth's 5514.
  What compacts a planet past the packing its own grains jam at is an equation
  of state for a solid, which is Phase 5's and is scheduled there.

`PLAY.md` D11 finds the largest standing axiom violation in the codebase:
**`morph::Program` is a species table.** Do not add a seventh variant — that is
what D11 exists to prevent. **Four of the five columns are still on it**, and
they are still there because Phase 4 measured each one against the laws the
engine actually has and found them out by four to forty-four orders. Read the
backlog entry and `Program::maintenance`'s doc comment before proposing a
derivation; the obvious ones have been tried and have numbers against them.

What Phase 2 took off `Program` is **geometry**, for anything assembled: a
composite's shape is a generated recipe and `Program` is provenance. Phase 4
took **the shape of everything else** — a habit is generated by
`recipe::generate` and a species names none of it — and **`substrate`**, which
is now the feedstock the node held. Those do not come back.

One decision in `PLAY.md` still changes text written down elsewhere, so do not
treat the older text as current where they disagree: **`PathKey` stops being an
identity.** The other — one second per second at every tier — is no longer a
plan; it is what `World::new` does.
