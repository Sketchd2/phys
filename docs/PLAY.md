# The play space

How the engine becomes a world people live in.

`DESIGN.md` says what the engine is. `PHYSICS.md` says what it models. This says
what has to be true before a person can stand on a planet, pick something up,
build something, break it, and share the result with somebody else — and it
records the decisions that were taken to get there, with the reasoning, so that
a later reader can tell a choice from an accident.

Nothing here is built yet. Everything here is decided.

---

## 1. What was actually missing

Three of the four things this looked like it needed already exist in skeleton,
and saying so matters, because the plan is much smaller if you know it.

- **Actor input exists.** `control.rs` is the inbound half of the seam:
  `Actor` with force, power, flow, reach and duration limits; `Act::{Push, Heat,
  Release, Sense}`; capability clamping with `Outcome::Clamped` reporting the
  fraction that survived; `Refusal`; wire encode and decode; a `Roster`.
  `tests/control.rs` holds the invariants — an act cannot state a fact, an act
  is bounded by the body, and what the world knows about an actor it measured.
  What is missing is not the vocabulary. It is that nothing consumes commands
  from many senders over time.

- **World generation exists.** `Program::Terrain` and `Program::Settlement`
  landed, and biomes fall out of insolation and water rather than a table. What
  is missing is not terrain as a type. It is that terrain is a *patch* with no
  parameterisation onto a planet, and nothing can stand on it.

- **Moving between areas exists.** `World::reparent` landed and `phys-rehome`
  demonstrates it. What is missing is anything that *triggers* it, and an
  identity that survives it.

**What is genuinely missing is narrower and harder.** The engine models what
happens *inside* a node extremely well and what happens *between* nodes barely
at all. `BACKLOG.md`'s coupling audit already reaches this conclusion: contact,
conduction, diffusion, radiative exchange and friction are absent, and they are
one missing primitive rather than five features. On top of that sit two facts
that make real-time interaction impossible today:

1. **A promoted node never feels a force.** `motion.velocity` is written at
   promotion and nowhere else. Measured: a promoted child and the parent body it
   stands for diverge by 0.79 of the child's own radius in forty frames. Two
   vehicles in one frame cannot collide, attract or perturb each other.
2. **Nothing rests on anything.** "The ground" is `ground_of`, one z-plane
   belonging to one structure. There is no surface in the world to stand on.

Those two, plus adjacency, plus an identity that survives a move, are the whole
foundation. Everything else in this document is built on them.

---

## 2. The decisions

### D1 — One clock, at every scale, always

**The world runs at one second per second.** Not within a band of tiers —
everywhere. `PaceMode::Fixed(1.0)` is what a shared world is, and `pace_to`
survives only as a tool for single-player exploration and for offline study.

This retires the engine's most quotable property. `README.md` says that zooming
into a nucleus does not slow the frame rate, it slows *time*. That is a
beautiful answer to a single-observer question and the wrong answer to a shared
one: a player inspecting a rifle bolt must not slow down the war.

**Why one clock at every tier is coherent rather than a wish.** A nucleus has a
stable timestep near 10⁻²³ s and cannot integrate one second of trajectory —
that would be 10²³ steps. It does not have to. `MAX_SUBSTEPS` already routes a
node that cannot be followed to its *equilibrium ensemble* instead of its
trajectory, and being carried to the instant statistically is both cheaper and
more nearly correct. So "every scale is at the same instant" is already what the
engine does. The only thing that changes is that the clock stops being dragged
slower by whoever is looking most closely.

**What gives way instead.** Overload no longer shows up as a slower world. It
shows up as a *staler* one: `stats.worst_lateness` rises and detail debt
accumulates. This makes the frame budget's job harder and more honest, and it
means the knapsack needs a rule it does not have today — a lateness ceiling for
anything an actor is interacting with, with resolution as the thing that is
surrendered to hold it.

**Time bubbles stay.** `Admin::Dilate` remains exactly as it is: localised,
unphysical, clamped, audited. That is the deliberate exception, and keeping it
in `Admin` rather than `Act` is what makes it one.

**Slow motion is a recording, not a pace.** A player who wants to watch water
freeze at a thousandth speed is asking to *review*, not to change the world's
clock — the freeze happened at one second per second and they want to look at it
again, closely. That is a replay, and the engine can already almost do it:

- regeneration is a pure function of `(world_seed, path_key, epoch, purpose, index)`;
- `causal::History` already keeps per-node rings of past state;
- the ledger already commits observed values permanently and never re-samples them.

What is missing is a **checkpoint plus an ordered input log**, from which a
bounded subtree can be re-run at an arbitrarily fine cadence into a buffer the
client scrubs, while the world clock never moves. The same two artefacts are
what multiplayer needs for rollback and what any audit of "what actually
happened" needs. One mechanism, three customers.

**The caveat, which must not be quietly dropped.** Replaying a region at finer
resolution than it originally ran gives *a* valid microhistory, not *the* one,
because detail that was never materialised was never decided. This is not a
defect to fix — it is axiom two, and observation-as-commitment, doing exactly
what they say. It also lands the right way round: the case a player cares about
is the one they were *watching*, and what they watched is in the ledger, pinned,
and replays exactly. Unobserved surroundings replay as a fresh draw from the
same distribution. The API must say which is which rather than presenting both
as the same kind of truth.

### D2 — Issued identity; `PathKey` demoted to an address

A node gains an `EntityId`: issued by the world, monotonic, persisted, never
reused. `PathKey` keeps the two jobs it is good at — deriving children and
seeding the sampler — and loses the one it is bad at, being a name.

**The prize is not tidiness.** Side tables (`mixtures`, `environments`,
`clocks`, `histories`, and everything gameplay adds) rekey to `EntityId`, and
the moment they do, `reparent` stops having to move them at all — they are keyed
by something a move does not change. `CLAUDE.md`'s standing trap, *add a table,
add a line to `World::reparent` or a moved object arrives without its
chemistry*, is dissolved rather than documented.

**It stays affordable** because most nodes never need a name. A regenerable node
is fully described by its address; only pinned, structured or touched things
need identity — and that is precisely the set the design already stores instead
of regenerating. The id table is bounded by the same argument as pinned detail.

Wire format is append-only, `FORMAT_VERSION` bumps, and `tests/persistence.rs`
catches size changes but not reorderings — so nothing existing may be reordered.

### D3 — Adjacency is one primitive, designed once

Contact, heat conduction, mass diffusion, sibling-to-sibling radiation, friction,
fire spread, flooding and debris landing on something other than its parent are
one sentence in different clothes: *two adjacent things exchange a conserved
quantity across the boundary they share*. The engine has no notion of adjacency
at all. `BACKLOG.md` warns that building the four consumers first would produce
four incompatible answers, and it is right.

**`Neighbourhood`: a per-node spatial index over that node's contents**, where
contents means its promoted children and its materialised bodies, treated
uniformly. Built lazily, cached, invalidated when the node's epoch moves.

- **Scale-free by construction.** Cell size is the node's radius over the cube
  root of its count — the same spacing SPH's smoothing length already uses — and
  everything is in node-relative units. The same code indexes a galaxy's arms
  and a nucleus's nucleons, because the numbers it sees are of order one either
  way. No per-tier variant, which is the test of whether this obeyed the axioms.
- **Across the hierarchy, not only between siblings.** Two adjacent things need
  not be siblings; a person on a hillside is a promoted node beside a terrain
  patch's contents. The query resolves at the lowest common ancestor — which
  `Tree::separation` already walks to — indexing there and descending into
  overlapping children. Precision then comes out right for free by the same
  argument `coords::Located` makes, and the work is bounded by the causal gate,
  because things too far apart to interact within a frame need not be
  enumerated.
- **Exchange is one function.** Given two adjacent things and a shared area, move
  a conserved quantity at a rate set by a transport coefficient derived from
  their `Mixture`s. Conduction, diffusion and radiative exchange are three calls
  to it with different coefficients. Delivery goes through `causal::Influence`,
  which is already a causally-ordered cross-node event.
- **Contact is the impulsive one.** Overlap resolves as an impulse pair through
  the existing path, with restitution and friction derived from the materials
  `topology.rs` already carries as data. This is where friction — currently
  absent above the Langevin thermostat — enters the engine.

### D4 — The promoted child is the authority

Of the two directions `BACKLOG.md` lays out, take the second: the child node is
real and the parent's body is its stand-in, so `sync_from_child` runs every
frame rather than only at `coarsen`. Cost is one write per promoted node per
frame. The parent's solver then sees the child's evolved position, and — the
half that does not exist at all today — the force the parent's solver computes
on that body is applied back to the child's `motion.velocity`.

This is what sibling interaction *is*. Without it there is no collision, no
mutual gravity between two ships, no debris that falls rather than floats, and
no reason to build any of the rest.

### D5 — Large rotation lives in the joint, not in the element

`dynamics.rs` is small-displacement linear about the reference geometry, exact
below about 0.1 rad of chord rotation and qualitatively wrong past 0.3. A
walking leg rotates through roughly a radian, so the obvious reading is that a
corotational element formulation is on the critical path for creatures.

**That reading is probably wrong, and it must be measured rather than assumed.**
`displacement_ratio` measures the transverse part of the relative displacement
*across each member*, over its length — the member bending, not the member being
carried. Its own doc makes the distinction: a twig at the end of a swaying
branch travels metres while rotating by almost nothing, because it is carried
rather than bent. A limb swinging about a hip is the same case: the large
rotation belongs to a *substructure moving rigidly*, and the members inside that
substructure barely deform.

So the decomposition is: **articulated segments carrying the large rotations,
the existing linear frame solver handling deformation within each segment.**
That is a corotational formulation in spirit — every segment carries its own
rotating frame — reached through substructuring rather than by rewriting the
element. And substructuring is already wanted twice over: `BACKLOG.md` needs it
for a vehicle whose load path spans promoted children, and D7 needs it so that a
scratch on a paw does not re-analyse a wolf.

`displacement_ratio` then stops being a limitation and becomes the *instrument*
that says whether the decomposition is adequate. Phase 0 probes it on a swinging
limb. If segments themselves turn out to exceed 0.1 rad, corotational elements
go back on the critical path — but that is a measurement, not a guess, and
guessing it is exactly the failure `CLAUDE.md` warns costs the most.

### D6 — Every planet walkable, by parameterisation not by type

A planetary surface is not a new kind of thing. It is a Planetary-tier node
which, on refinement, divides its sphere into **patches** — and a patch is a
node carrying `Program::Terrain`, which exists.

- **Cubed-sphere parameterisation.** Six faces, each patch a square in a local
  frame, each patch's children four squares. Chosen over HEALPix and geodesic
  subdivision for one reason that outranks their better area properties: the
  child relation is a clean four-way split, which is what `PathKey`'s
  child-index derivation consumes, so the surface tree *is* the scale tree with
  no adapter. Cubed-sphere's angular distortion is bounded and correctable at
  the patch level; a pentagon defect or a nested-ring indexing scheme is not.
- **Nothing is generated until approached.** A patch refines into sub-patches as
  an observer descends and coarsens behind them, and the terrain regenerates
  bit-identically — the invariant `tests/consistency.rs` already guards. A
  planet nobody has visited costs its `Matter` and nothing else.
- **Height derives, biomes already do.** Patch relief comes from address-derived
  noise conditioned on properties the node already carries; insolation and water
  already produce biomes without being told what a desert is.
- **Local gravity derives.** `drop_fragments` currently loads every falling piece
  with `G_EARTH`, a constant, while the engine computes real gravitational
  fields at every other tier. That constant goes; g comes from the planet's own
  mass and radius. A backlog entry closes as a side-effect.
- **Walking off the edge of a patch is a `reparent`**, triggered by leaving the
  patch's volume. That trigger does not exist, and it is the same spread
  measurement node-splitting needs — shared machinery, as `BACKLOG.md` predicted.

### D7 — A creature is a body plan, mechanisms, and derived shortcuts

**Not humanoid, and not a special case.** `Program::Creature` generates a
skeleton — segments, articulated joints, and attachment points — from a genome,
exactly as `Program::Tree` generates a tree. A biped, a quadruped and something
with six legs are three genomes, not three code paths. Wolves, pets and people
are the same machinery.

**Motion is mechanisms.** `Mechanism` gains one variant — an actuation, a force
between two joints along their connecting line. Muscles, hydraulic rams and
thrusters are then one thing, and the vocabulary that already carries wind,
snow, lightning and fire carries locomotion too. Standing is contact against
D3's adjacency relation. Nothing is told what a leg is.

**Shortcuts are derived, then stored** — axiom three, applied to behaviour
rather than to matter. Solving a gait from muscle activations every frame is not
the goal; solving it *once* and storing it is:

- A **gait** is an optimisation over actuations for a given tuple of body plan,
  mass distribution, local gravity and surrounding medium. Solved once, cached
  under that tuple, thereafter a lookup and a phase advance.
- The **invalidation condition is the tuple itself**, which is what makes this a
  derived shortcut and not a table. Break a leg, shoulder a heavy load, walk
  onto a low-gravity moon, wade into water — the tuple changes, the gait is
  re-derived, and the creature moves differently with nobody having written a
  rule about moons.
- A **grasp** is the same shape: feasible hold configurations solved once per
  manipulator topology and object shape class, cached, and a held object becomes
  a temporary joint in the carrier's topology. Its mass is then genuinely in the
  load path, so a creature can be overloaded and that falls out rather than
  being checked for.

**Interactable at any scale without re-analysing the whole animal.** This is
D5's substructuring doing its second job. A local interaction resolves inside
the affected segment against its condensed boundary, and only propagates to the
coarse solve when it is large compared to that condensation's own tolerance —
which is a measured criterion, not a policy about what counts as a big hit.

**The mind is not part of the creature.** See D8.

### D8 — A mind is an actor, and it lives outside the engine

Behaviour is authored, and none of it is a law. A mind — a wolf's, an NPC's, a
door's, a factory's — holds an `ActorId`, perceives only through `Sense`, and
acts only through capability-clamped `Act`s. It is the same door a player uses,
with the same clamps and the same refusals. It cannot state a fact, because
`Act` has no variant that states one; it cannot reach `Admin`, because it cannot
construct that variant. Those are type-level guarantees `control.rs` already
holds, inherited for free.

**The engine cannot tell a player from a wolf, and that is the property worth
having.** One mechanism covers NPCs, pets, scripted objects, automation, and the
long-standing `DESIGN.md` limitation that *nothing builds the buildings* —
`Program::Tower` already advances on a supplied labour rate, and a mind is what
supplies it.

**Where a mind runs is settled by a constraint, not a preference.** `Cargo.toml`
is explicit: the core crate has no dependencies and must keep building for
wasm32, which is why the Postgres store is an optional feature. A WebAssembly
script host is a large dependency and cannot go in the core. It does not need
to: a mind talks to the world through `Command`, which already encodes, decodes
and crosses a process boundary — so **a mind is a client**, exactly like a
player's, and the engine hosts nothing. Sandboxing, fuel and language choice
become properties of a separate host crate, and the engine's dependency-free
guarantee survives intact.

**Determinism is a constraint on that host, and it is not optional.** A mind's
outputs enter the input log, and replay (D1) and multiplayer (D10) both require
that the log replays identically. A guest may not read a clock, may not use
threads, and gets its randomness from the engine's addressed streams like
everything else.

### D9 — A built thing is a build log

A player-made object is a `topology` of members and joints, and the state that
persists is the **ordered build log**: what member, of what material, between
which points, placed by whom and when. That log is a genome in `morph`'s
existing sense — the geometry derives from it — so it stores in the hundreds of
bytes that a structure's developmental state already does.

Grid snapping is a *client* convention. The wire carries a placement; the engine
never learns what a cell is. This is what keeps a block game and a free-form
builder from being two engines.

Everything downstream comes free: the frame solver analyses it, damage and
fragmentation apply to it, `dynamics` shakes it, and it breaks the same way a
tree does.

**One thing is deliberately not free, and it is a divergence worth writing
down.** Generated structures run a fully-stressed design pass — the generator
places members and the pass sizes each one for the load it actually takes. A
player who places a member has *chosen* its section, and resizing it behind them
would be the engine overruling the builder. So player-built work is analysed but
never proportioned, which means players can build things that fall down. That is
the point; it is also the first place where two structures in the same world
were made by different rules, and it should stay the only one.

### D10 — One authoritative world, and a mandatory input log

Designed now, built after the first slice, so that nothing has to be unpicked.

- **Authority is the server's.** Clients send `Command`s and receive `Scene`s.
  Both already encode; `phys-headless` already proves the split is real rather
  than notional.
- **The input log is mandatory from the first line of code.** World state is
  `f(world_seed, ordered input log)`. Every accepted command is appended,
  timestamped and sequenced. This costs a write per command and buys replay,
  hindsight review, rollback, audit, and the ability to reproduce any bug
  exactly. It is the highest-value item in this document per line of code.
- **The client predicts only its own avatar, through the same mechanisms.** No
  second physics implementation, ever. Divergence is corrected by the next
  `Scene`.
- **Interest management already exists** as `Volume` queries with distance LOD,
  and bandwidth already tracks screen rather than world through recipes and
  per-client deltas.
- **The sharding unit is a subtree**, not built now, but `EntityId` and the input
  log are what make it possible later.
- **What must be true immediately** so this is not a rewrite: single-player runs
  a server and a client in one process, talking through the same encode and
  decode. The moment anything reaches across that seam, it will keep reaching.

---

## 3. What has not been measured

Per `CLAUDE.md`'s first trap, these are stated as unmeasured rather than
assumed, and each has a probe in Phase 0.

| Question | Why it matters | Probe |
|---|---|---|
| What does a neighbour query cost at 10³–10⁵ contents? | D3 is on every frame's critical path. If it is expensive, the frame budget arithmetic in `PERFORMANCE.md` changes. | Build the index over the existing scenarios and time it. |
| Do segments of a walking limb exceed 0.1 rad of chord rotation? | Decides whether D5's decomposition suffices or corotational elements go on the critical path. | Drive a two-segment limb through a stride and read `displacement_ratio`. |
| Does slaving the parent body every frame preserve `summarise(sample(m)) == m`? | D4 writes into the conserved set every frame. `IDEMPOTENT_TOLERANCE` is the contract. | Promote, run, coarsen, compare against the existing consistency harness. |
| How long does a gait optimisation take, and does it converge? | D7's shortcut is only a shortcut if deriving it is rare and bounded. | Solve one quadruped gait offline and time it. |
| How large is the checkpoint for an interactive subtree? | D1's replay and D10's rollback both pay for it. | Measure a populated patch's snapshot through the existing `persist` path. |

Two known defects will bite during this work and are scheduled rather than
discovered:

- **Tier is cached and never revisited.** `plant`, `emplace` and growth all change
  a node's size without updating its tier, so a node can be two tiers from what
  its radius says. An avatar walking across patches, and a patch refining under
  it, are both size changes. This must be fixed in Phase 1, not deferred.
- **A node cannot split when its contents spread out.** The spread measurement it
  needs is the same one D6 needs for patch handoff and `BACKLOG.md`'s fragment
  entry needs for promote-on-leaving. Build it once.

---

## 4. The order of work

Each phase ends with a test that fails today. A phase is not done because its
code exists.

**Phase 0 — Probes.** The five measurements above, each as a test that
demonstrates the thing it claims. Nothing is designed further until they are
numbers. *Done when:* `PERFORMANCE.md` carries five new measured rows and D5 is
either confirmed or replaced.

**Phase 1 — The primitives.** `EntityId` and the side-table rekey (D2);
`Neighbourhood`, boundary exchange and contact (D3); the promoted child made
authoritative and force-bearing (D4); tier revisited on size change; the spread
measurement; `PaceMode::Fixed(1.0)` as what a world is, and `G_EARTH` deleted in
favour of derived g. *Done when:* two promoted vehicles collide and rebound with
restitution derived from their materials; a hot node beside a cold one
equilibrates without either being told the other exists; a branch lands on the
next tree; and nothing in the existing suite regresses.

**Phase 2 — Ground.** Cubed-sphere parameterisation, patches as `Program::Terrain`
nodes, refinement and coarsening on approach, handoff by `reparent`, planetary
gravity. *Done when:* an observer descends from orbit to a square metre of any
planet in any scenario, travels ten kilometres across patch boundaries, and the
terrain behind them regenerates bit-identically.

**Phase 3 — Bodies.** Substructuring (D5); `Program::Creature` and its genome;
the actuation mechanism; derived-and-cached gait and grasp; local interaction
resolved inside a segment. *Done when:* a quadruped and a biped grown from two
genomes both walk, on two planets with different g, with no per-morphology code
anywhere; shouldering a load visibly changes the gait; and scratching a paw does
not re-analyse the animal.

**Phase 4 — Minds.** The actor-client host outside the core crate; the
determinism constraints; checkpoints and the input log. *Done when:* a wolf
pursues something, the engine cannot distinguish it from a player, and a
recorded session replays bit-identically from seed plus log.

**Phase 5 — Making things.** The build log as a genome; player-placed members;
analysis without proportioning; break and repair. *Done when:* a player builds a
bridge that holds and one that does not, and the engine was never told which was
which.

**Phase 6 — Sessions.** The server loop; two clients; prediction of one's own
avatar only; interest management under load; hindsight replay scrubbing on top
of Phase 4's checkpoints. *Done when:* two people share a world, one builds
while the other watches, and either can scrub back through a minute of it at a
thousandth speed without the world's clock moving.

---

## 5. The axioms, re-checked

Nothing above is worth building if it breaks the five things `CLAUDE.md` says
are not preferences.

- **Only physical law is axiom-side.** No decision here adds a named entity to
  the engine. A creature is a genome, a gait is an optimisation, a building is a
  log of placements, a biome is already a consequence. The one place a table was
  tempting — block types — was refused in D9, and the one place behaviour was
  tempting to make engine-side — authored devices — was refused in D8.
- **Measure, never be told.** D8 is this axiom applied to minds: a wolf learns
  what it can lift by trying and being told the fraction that got through.
  Ground contact is measured from geometry rather than from a grounded flag.
- **Derived, with shortcuts stored.** D7's cached gait is the clearest instance
  in the codebase: derived once, stored, and re-derived when the tuple it was
  derived under changes.
- **Detail exists where something is happening.** D6 generates a planet's surface
  only where somebody is standing, and D1's replay respects the same thing —
  which is exactly why an unobserved replay is a fresh draw and an observed one
  is exact.
- **No special cases.** D3 is the load-bearing one: contact, conduction,
  diffusion, radiation and debris landing become one primitive rather than five.
  If a later change needs a sixth mechanism to make one of them work, that is
  the signal the primitive was wrong.
