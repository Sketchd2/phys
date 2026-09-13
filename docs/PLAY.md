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
`Program::Tower` already advances on a supplied labour rate, and a mind actually
doing the work is what may supply it.

Note the boundary this draws, because §5.4 turns on it. A mind **building**
something is an actor spending real effort on a real act, and it is welcome. A
notional labour rate that **maintains** things nobody is working on is not a
mind, it is a number, and it is the thing §5.4 throws out. Construction by an
actor: yes. Ambient maintenance as a law: never.

**A mind is a client, and the reason is not the one it first appeared to be.**
`Cargo.toml` says the core crate has no dependencies and must keep building for
wasm32, and the first draft of this decision leaned on that: a script host is a
large dependency, therefore it cannot go in the core. That constraint is not
binding. The wasm build is real — it is what `viewer/` runs, and `src/wasm.rs`
is the C ABI it drives — but it is one target among several, not a rule that
forbids the core from ever taking a dependency. The reasoning has to stand on
its own, so here it is without the crutch:

A mind is an actor, an actor talks to the world through `Command`, and `Command`
already encodes, decodes and crosses a process boundary — `phys-headless` proves
it. Putting a mind on the far side of that boundary buys three things an
in-process host does not. The engine genuinely cannot tell a player from a wolf,
because both arrive the same way, so there is no second door to keep in step
with the first. Determinism and sandboxing are enforced by a serialised
interface rather than by trusting a guest sharing an address space. And the
sharding unit stays a subtree, because nothing about a mind is pinned to the
process holding the world.

**An in-process host is therefore an optimisation, not the architecture.** If
latency measurement later says a thousand wolves cannot afford a round trip, a
host crate behind a feature flag may run their guests in the engine's process —
but it must produce byte-identical `Command`s to the out-of-process path, and
the out-of-process path stays the definition of correct. That ordering is what
stops the fast path from quietly becoming a second, more privileged door.

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

### D11 — One deposition law; `Program` becomes a genome, not a species

**`Program` is a species table, and it is the largest standing violation of the
first axiom in the codebase.** Six variants, fourteen `match self.program` sites,
and an `impl` carrying seven per-variant columns:

```rust
is_planned()     Tower | Wall | Settlement
design_flow()    Tree (20.0, 1.225)  Coral (1.2, 1025.0)  Tower (42.0, 1.225) …
density()        WOOD_DENSITY | CORAL_DENSITY | BUILDING_DENSITY | ROCK_DENSITY
energy_density() BIOMASS_ENERGY | CONSTRUCTION_ENERGY | 0.0
substrate()      a Composition per variant
maintenance()    Tree 0.02/YEAR  Coral 0.05/YEAR  Tower 0.005/YEAR  Terrain 1e-5/YEAR
```

"Nothing in the engine knows what salt is, what a forest is" — but
`Program::Tree` knows what a tree is, knows it grows in air at 1.225 kg/m³,
knows it is made of wood, and knows it decays at 2% a year. That last row is a
tabulated decay rate, which is precisely what §5.2 says must be derived. And D7
was about to add a seventh row for creatures, with more columns.

**Five of the seven columns are not properties of a species at all.** They are
properties of two things the engine already measures:

| column | what it actually is | where it comes from |
|---|---|---|
| `density` | the material | `topology.rs` already has materials as data. A wooden tower and a wooden tree have one density between them. |
| `energy_density` | embodied energy per kg | the material again — and §5.8 needs this real anyway |
| `substrate` | what it is made of | the material's composition |
| `design_flow` | the fluid, and the gust to build against | **measured.** A coral is in water because its node's mixture is water; nobody tells it. The design gust is the gust the structure has met, which its own history already holds. |
| `maintenance` | the decay rate | §5.2: cohesion, local flux, geometry |

So the table is not needed, because the information is already in state. That is
"measure, never be told" applied to the thing that was telling.

The sixth, `is_planned`, is the interesting one. Physically, growth and
construction are the same operation — **material is deposited where a field says
to deposit it** — and they differ only in where the field comes from. A tree's
comes from light and load; a tower's is supplied by whoever is building. So
`is_planned` stops being a species flag and becomes a statement about the
deposition field's *source*, which is a property of the act rather than of the
thing.

**Decay is the same law with the flux reversed.** Deposition where it pays,
removal where it does not: a rotting log, a weathering wall, a dissolving bone
and an eroding dune are one expression over material resistance, local flux and
surface geometry. §5.2 already wrote it for terrain; there is no second version
for structures.

**And weather stops being a feature.** `Mechanism` is *already* the general load
vocabulary — `BodyAcceleration`, `FlowDrag`, `SurfaceAccretion`,
`ConductedEnergy`, `ThermalField` — and `PHYSICS.md` already says snow, wind,
lightning and fire are constructors over them rather than named weather. What is
missing is not generality, it is *reach*: mechanisms are handed to `shake` and
`damage` by the caller, so weather today is a test harness rather than a thing
the world has.

D3's adjacency is what closes it. Once a structure can measure what is next to
it, the adjacent air has a velocity, a temperature and a mixture with a water
fraction, and `FlowDrag` is constructed from the measurement rather than passed
in. Then **one environment measurement drives three consumers at once**: a tree
in a gale is loaded by the wind, grows thicker because it is loaded, and weathers
faster because of the driven rain. No weather system exists anywhere, and a storm
is a thing the air is doing.

**What this does not unify, and I will not claim it does.** The rate laws go to
zero species knowledge. The *geometry* does not. Branching, coursed masonry and a
subdivided street grid are genuinely different space-filling rules, and asserting
they collapse into one would be the third over-claim in this document rather than
the first.

What they do reduce to is a small set of **habits** — branching (trees, corals,
lungs, river deltas, lightning), coursed (walls, brickwork, strata), subdivided
plane (settlements, cracked mud, leaf venation) — selected and parameterised by
the genome instead of named by an enum. A creature needs a fourth, segmented and
bilateral, and that is an honest addition rather than a per-species renderer.
Three or four habits against six-and-climbing species is a real reduction.

**D12 moves this line**, for the derived-field half: branching in particular does
not have to be a habit, because it is what transport into an occluded field looks
like. Coursed masonry and the planned programs stay habits, because a wall is
placed rather than grown.

**The honesty test, stated before the work rather than after.** Can a genome
*nobody designed* produce something coherent? If every genome that works had to
be hand-tuned into working, then thirty coefficients have replaced six enum
variants, the special cases are still there, and they are now less legible than
when they were honest about being a table. Sample a hundred random genomes; if
the failures are ugly rather than impossible, it generalised.

**Where it lands.** The material-and-environment collapse goes in Phase 2, because
terrain decay needs `maintenance` derived before it can erode anything. The habit
refactor lands before Phase 4, because creatures would otherwise arrive as a
seventh species and make the table worse at exactly the moment it is hardest to
undo.

---

### D12 — Shape derives too, for anything grown

D11 stopped at habits: branching, coursed, subdivided plane, selected by genome.
**That line is in the wrong place, and this moves it** — for the half of D11's
split where the deposition field is *derived* rather than supplied.

**What decides a tree's shape today.**

```rust
Program::Tree  => self.render_branching(n, 0.62, 3),   // taper, splits
Program::Coral => self.render_branching(n, 0.72, 4),
```

`taper` and `splits` are constants at the call site, per species. `lean`, `twist`
and `spread` come from the genome, which is per-instance and right. And `render`
is pure in `(genome, age, built, progress, events)` — **`Environment` is not one
of its arguments.** A tree's shape cannot respond to where the light is. It
responds to its genome and its mass, and nothing else.

The tell is in the doc comment. It says the taper obeys da Vinci's rule, that
total cross-section is preserved across a branch point — a real law, which would
*supply* the number: `splits · taper³ = 1` gives 0.693 for three-way branching
and 0.630 for four-way. The constants are 0.62 and 0.72, giving 0.715 and 1.493.
One structure loses a quarter of its cross-section per level and the other gains
half. **A law was cited and two hand-picked numbers were used instead**, and
because both are constants nothing can notice they disagree with it.

**Branching does not have to be written down, because it is what transport into
an occluded field looks like.** Four mechanisms, all derived, all standard:

1. **A resource field with occlusion.** Tips compete; a branch shades what is
   under it; material goes where the unshaded resource is. This alone produces
   branching — it is what Laplacian growth, viscous fingering, Lichtenberg
   figures and river networks all are.
2. **Transport cost.** Murray's law — minimise pumping plus the cost of
   maintaining the conduit — *derives* the branching ratio. That is `taper`,
   computed rather than tabulated, and it is why the exponent is somewhere
   between the area-preserving 2 and the flow-optimal 3 rather than being
   whichever of 0.62 or 0.72 somebody typed.
3. **The mechanical constraint**, which the engine already has in full: the
   fully-stressed design pass already sizes every member for the load it
   actually takes, and already re-proportions against a load envelope.
4. **Lateral inhibition** for spacing, which is what `twist` and `spread`
   currently stand in for.

**Two ingredients are missing, and one of them is D3 again.** The engine has (3)
outright and has the resource field for growth *rate* but not for growth
*direction*. It needs **occlusion** — a part shading the parts behind it, which
is the adjacency relation applied inside a structure — and **transport cost**
along the structure, which nothing models. That is the whole gap.

**What the genome becomes, and this is the part that answers the question as
asked.** Not shape. DNA does not encode branch angles; it encodes proteins and
thresholds, and shape falls out of those meeting a particular patch of ground. So
the genome becomes physiological constants: transport efficiency, the material,
which resource is limiting, shade tolerance, the allocation between growing and
defending. Shape is then an *outcome* of those constants meeting a light field,
in the same sense that it is in a real tree.

**The consequence that makes it worth doing, and it is falsifiable.** The same
genome grown in a forest comes out tall and thin; grown in the open it comes out
short and spreading. That is what real trees do, and it is **impossible today** —
`taper` and `splits` are constants and `render` never sees the environment, so
the only way to get two shapes is two genomes or two species. If the derived
version cannot produce that difference without changing the genome, it did not
derive anything.

And the same law is the branching in river deltas, lightning, veins, lungs,
cracks and street networks. This does not delete `render_branching`; it deletes
the *category* that `render_branching` was the first member of, which is the
"no special cases" axiom collecting a large debt.

**Derive per species, instance per individual.** This is the third axiom applied
to morphogenesis, and it dissolves most of what follows.

The growth *response* is derived once per genome — offline, over a sampled space
of conditions — and stored. Every individual then evaluates the stored response
at its own conditions, which is a lookup rather than a simulation. An oak and a
birch are two genomes and therefore two response surfaces: same law, different
physiology, different tree. The expensive derivation happens once per species,
not once per tree and never per frame.

Note what makes this axiom-compliant rather than a table with better manners.
Axiom three's operative clause is "a tabulated constant *that was never derived*
is the thing to avoid" — so the objection is to provenance, not to storage.
`taper = 0.62` is forbidden because nothing produced it. A surface computed from
Murray's law and light competition is exactly what the axiom asks for, even
though both end up as stored numbers.

It is also the third use of one pattern: D7 caches a gait per
`(body plan, mass, gravity, medium)` and a grasp per `(manipulator, object
shape)`, and this caches a morphology response per genome. Three independent
arrivals at the same shape is a good sign it is the right one.

**Three costs, of which the first two are now answered.**

**Answered — it is iterative.** Growing into an occluded field is not one
`render` call, and the design's single biggest win is that "growth runs on the
aggregate, so a forest grows without any of its trees existing". A million trees
cannot each run a light-competition simulation. Under derive-once-per-species
they do not have to: a million trees is a million evaluations of a stored
surface. The simulation runs once, for the species, and never again.

**Answered — regenerability.** The worry was that shape would depend on the
environment a thing *grew in*, which is a history rather than a state, so a node
could no longer throw its detail away and rebuild it bit-identically. The
resolution is better than storing a summary of that history: **the environment
does not need storing, because it re-derives.** In a deterministic world the
conditions at a place are a pure function of that place and that time — the
terrain, the latitude, the climate all descend from the world seed. The same node
always had the same weather. So shape stays pure in `(genome, place, age)` and
the founding invariant is untouched.

The qualifier is the one that was already stated: *without significant triggering
factors*. A bushfire is a perturbation, and a perturbation is a deviation, which
§5 already handles as an event that either decays or is promoted. So a structure
is `derived_baseline(genome, place, age) + events` — which is exactly the
existing `structure = program(genome, age, events)` with `place` added and
`program` replaced by something that was derived rather than typed.

**The remaining risk, and it is now the only real one: is growth Markovian in
its conditions?** A response surface tabulates shape against *conditions*. If
shape instead depends on their *ordering* — a drought at age ten making a
different tree from a drought at age fifty — then the input is a history rather
than a point, and it does not tabulate at all. Real trees certainly record their
history, but mostly *internally*, in rings and reaction wood and scars, rather
than in gross form; height, spread and taper are plausibly driven by integrals.
Plausibly is not measured, and this one decides whether the cheap version exists.

Dimensionality is the second-order version of the same worry. Six conditions —
light, water, temperature, crowding, wind, soil — at five samples an axis is
about 15,600 growth runs per species, which is fine once and offline, and grows
badly if the list gets longer. Which conditions actually move the shape is worth
measuring before committing to a surface over all of them.

**Bounded: it does not apply to built things.** A wall is not grown into a field;
it is placed by an intention. D11's supplied-versus-derived split is exactly the
right boundary and it is not a species boundary — it is whether something meant
it. Coursed masonry stays a habit. A street grid is interestingly *both*, since
real street networks do follow transport optimisation, and the split is per
structure rather than per kind.

**The honesty test.** Grow one genome in three light fields and get three shapes
a botanist would call the same species in three situations. Grow two genomes in
one field and get two species. If the forest-versus-open difference needs a
genome change, nothing derived and the constants merely moved.

---

## 3. The tiers under these decisions

The ladder is the part of the architecture most likely to be assumed changed by
all of the above, so this says explicitly what happens to it. Short version:
**the ladder does not move, the play space occupies barely any of it, and one
clock exposes a boundary inside `Continuum` that observer-following pace was
hiding.**

### 3.1 The ladder does not move

Seven tiers, same boundaries, same meaning: a tier is a physics regime, not a
tree level, and `Tier::containing(metres)` still derives it from size. Nothing
in D1–D10 argues for a new tier, a moved boundary, or a different solver
assignment. A regime is a fact about physics and the decisions above are facts
about a game.

### 3.2 The play space is one tier and a bit

`Continuum` runs from 10⁻⁸ m to 10⁴ m — twelve orders of magnitude — and it
contains a grain of sand, a person, a wolf, a tree, a building, a settlement and
a ten-kilometre terrain patch. Above it, the bottom of `Planetary` (10⁴ to 10⁹ m)
holds the planet itself and a very large ship.

So **essentially the whole game is one tier**, and almost every level-of-detail
step in gameplay is a refinement *within* `Continuum` rather than a tier change.
That is exactly what `units.rs` says tiers are for — "many refinements happen
within one tier" — and it is good news: the game does not lean on the part of
the architecture that spans thirty-eight orders of magnitude. It leans on the
part that refines within one.

### 3.3 Dispatch must read state, not only size

`solvers::for_tier(Continuum)` is `Hydro`. A building, a wolf and a boulder are
all `Continuum`, and none of them is a fluid.

They escape SPH today, but implicitly: ordered matter enters through `shake` and
`damage`, which require `morphology` and `topology` and never consult
`for_tier`, while disordered matter enters through `advance_node`, which
consults nothing else. Dispatch is already state-dependent — it is just split
across two entry points that no single node can straddle.

The play space is made of nodes that must straddle them. A room with furniture
and air in it. A ship with a hull and an atmosphere. A creature standing in
water. So the two paths become one: **the tier says which regime the disordered
contents are in, and the node's own state says which contents are ordered**, in
one node, in one pass. This is not a new tier and not a new solver; it is
"measure, never be told" applied to solver selection, made explicit instead of
emergent.

### 3.4 One clock puts a resolution floor inside `Continuum`

This is the consequence of D1 that the first draft of this document missed, and
it is the most important thing in this section.

`node_dt` is `min(tier.dt(), dynamical_time/50, 0.25 · h/c_signal)` with
`h = radius / parts^(1/3)`. At one second per second and twenty updates per
second a frame covers 50 ms of world time, and `MAX_SUBSTEPS` is 256. Measured
by drilling the galaxy scenario from root to a half-millimetre node and reading
`node_dt` at every step:

```text
tier             radius        parts    node_dt (s)   substeps
galactic       4.629e20        20000       3.156e12        0.0
planetary       3.305e4         4000       1.429e-2        3.5
continuum       1.041e3         4000       4.785e-4      104.5
continuum       3.279e1         8000       1.104e-5     4529.4   <- ensemble
continuum      8.198e-1         8000       2.675e-7   186917.8   <- ensemble
continuum      5.124e-4         8000      1.727e-10 289551456.3  <- ensemble
```

**The trajectory path runs out inside `Continuum`.** Everything coarser is free —
the sky costs nothing, because a megayear-step regime crossing 50 ms is a coast.
`Planetary` is comfortable at a few substeps. And then within one tier the
requirement climbs through six orders of magnitude and falls off the end.

Where exactly it falls off depends on the material, not on the tier. Rearranging
the same expression, the finest a fluid can be resolved and still be *followed*
is `h ≥ 4 · frame_span · c_signal / MAX_SUBSTEPS`, which at a 50 ms frame is
`c_signal / 1280`:

| medium | signal speed | finest followed |
|---|---|---|
| air | 340 m/s | 0.27 m |
| water | 1500 m/s | 1.2 m |
| rock | 5000 m/s | 3.9 m |

The measured table above crosses over at a much coarser radius than the air row
suggests, because the galaxy scenario's `Continuum` gas is hot and its signal
speed is tens of km/s. Both say the same thing; the crossing point moves with
the material.

**Under observer-following pace none of this was visible**, because zooming in
slowed the clock until the substeps fit. Fixing the clock is what surfaces it.
That is not an argument against D1 — it is D1 doing what it was chosen to do,
which is to make the cost land somewhere honest instead of in a silently slower
world.

### 3.5 Solids escape the floor, and that is most of the game

`dynamics.rs` integrates structures with trapezoidal Newmark-beta, whose own
module doc says it "removes the stability limit entirely". `shake` substeps for
*accuracy* — the structure's own period — not for stability, capped at
`MAX_SHAKE_STEPS = 240`.

So the floor applies to **free fluid only**. Everything a player is, touches,
builds or breaks is ordered matter on the unconditionally-stable path: a
creature resolved to the centimetre, a building to the member, a ship to its
frame, all stepped across a 50 ms frame without a CFL condition anywhere. The
constraint bites on air in a room, water in a lake, smoke, and the blast from an
explosion — bulk flow, where a quarter-metre cell is coarse but not absurd, and
splashes and flames, where it is.

### 3.6 Below `Continuum`, replay is the only access

At 1 s/s a `Molecular` node needs ~5×10¹³ substeps to cross a frame, `Atomic`
~5×10¹⁶, `Nuclear` ~5×10¹⁹. They are always crossed by their ensemble. That was
already true and is not a regression — but with the clock fixed it becomes
permanent, and it has a consequence for what a player can ever see:

**A player cannot watch chemistry happen in real time.** They can watch its
consequences — the mixture changes, a phase changes, something freezes — because
those are `Continuum` facts. If they want to watch the mechanism, they use D1's
hindsight replay, which re-runs a bounded subtree at whatever cadence it likes
precisely because the world clock is not moving.

So the replay facility is not a luxury feature bolted onto the time model. Once
the clock is fixed, **it is the only way the fine tiers are ever directly
observable at all**, and that is the argument for building it rather than the
convenience of scrubbing.

### 3.7 What is decided, and the one thing that is not

Decided: the ladder stands; dispatch reads state as well as size; the resolution
floor is real, is derivable, and the engine should *report* it rather than
silently drop a node to its ensemble — the same discipline `displacement_ratio`
already applies to the small-displacement regime.

Open, and deliberately not decided here: **whether `Continuum` eventually gets
an unconditionally-stable fluid option**, so that free fluid stops being
CFL-bound the way solids already are not. Three ways to respond to the floor,
in the order they should be tried:

1. **Accept it.** This is D1 working as intended — detail gives way, not the
   clock — and metre-scale bulk air may simply be adequate. Costs nothing.
2. **Reduce the signal speed rather than the timestep.** The beach test (§4.2)
   found this after the list was first written, and it belongs second because it
   is far cheaper than (3) and buys most of what (3) buys. Weakly-compressible
   SPH replaces the physical sound speed with an artificial one about ten times
   the flow speed, chosen so density varies under a percent: for a 2 m/s wave,
   20 m/s instead of 1500, and a scene needing 30,000 substeps needs 80. It is
   an approximation with a known error bound rather than a fudge, and it is the
   standard treatment for exactly the free-surface case the play space is full
   of. It does nothing for a genuinely stiff medium.
3. **Raise `MAX_SUBSTEPS` for the play space.** Available immediately and
   directly buys resolution, but 256 is already a number nothing derives, and a
   larger one spends frame budget on exactly the nodes with the most of it.
   Phase 3 expects to need 512.
4. **Give `Continuum` an implicit integrator.** The principled answer and much
   the largest piece of work.

Try them in that order, and let a measurement decide when to move on. Choosing
(4) now would be the same mistake as assuming corotational elements were needed
in D5 — and note that the beach test already moved the play space past (1),
which is what a good acceptance scenario is for.

---

## 4. The beach test

One scene, checked link by link, because a plan that cannot be falsified by a
concrete question is not a plan:

> An actor is standing on a beach. They carve a 5 cm channel in the sand. Do the
> waves flow through it, in real time?

**As written, no.** Not because of hardware — the particle counts are
comfortable — but because two subsystems the scene needs are absent from the
engine and, worse, absent from every phase below. Recording it here because
finding that out is what the question was worth.

### 4.1 The chain, link by link

| # | What the scene needs | Status |
|---|---|---|
| 1 | A planet with a beach on it | Phase 2 |
| 2 | An actor standing on it, not falling through | Phases 1 and 3 |
| 3 | Sand that can be **carved** and stays carved | **not designed** |
| 4 | An ocean that has a **surface** | **absent** — see the coupling audit |
| 5 | Waves, i.e. free-surface gravity waves with a driver | **absent**, needs 4 |
| 6 | Water followed at ~1 cm at 1 s/s | **blocked on a scheme choice** |
| 7 | Water meeting arbitrary 5 cm sand geometry | **not designed** |
| 8 | Coarse ocean feeding the fine channel | **not designed** |

### 4.2 The timestep, measured

The interesting link is 6, and it is not the one that looked hardest.
`node_dt`'s CFL term is `0.25 · h / c_signal`, so for a 50 ms frame:

```text
water: rho = 1000.0 kg/m^3
  pressure()            = 4.0366e8 Pa
  velocity_dispersion() = 635.3 m/s
  sound_speed()         = 820.2 m/s   <- what node_dt uses
  real water            = 1500 m/s

signal speed                  h (m)      dt (s)   substeps
engine sound_speed()         0.0500   1.524e-5     3280.9  over
engine sound_speed()         0.0100   3.048e-6    16404.7  over
real water 1500              0.0100   1.667e-6    30000.0  over
WCSPH, 10x a 2 m/s wave      0.0500   6.250e-4       80.0  ok
WCSPH, 10x a 2 m/s wave      0.0156   1.950e-4      256.4  over
WCSPH, 10x a 2 m/s wave      0.0100   1.250e-4      400.0  over

particles, 2 m x 2 m x 0.1 m sheet of swash:
  h = 0.0500 m -> 3.200e3
  h = 0.0100 m -> 4.000e5
```

Three things fall out.

**The wall is a scheme choice, not hardware.** 4×10⁵ particles at 1 cm is inside
the 10⁵–10⁷ band `PERFORMANCE.md` projects. What fails is the timestep, and the
timestep is set by the signal speed the solver assumes. Weakly-compressible SPH
— the standard treatment for free-surface flow — replaces the physical sound
speed with an artificial one about ten times the flow speed, chosen so density
varies by under a percent. At a 2 m/s wave that is 20 m/s instead of 1500, and
the same scene goes from 30,000 substeps to 80.

**So the answer is a near miss rather than a fantasy.** With WCSPH, 5 cm cells
fit in 80 substeps; 1.56 cm sits exactly on the 256 cap; 1 cm needs 400. To see
water flow *through* a 5 cm channel you want three to five cells across it, so
~1 cm — which needs `MAX_SUBSTEPS` at 512 and about 4×10⁵ particles in the local
patch. That is a defensible configuration, not a wish.

**There is no condensed-matter equation of state.** `pressure()` returns
4×10⁸ Pa for water at room conditions, because the equation of state is ideal
gas everywhere; `sound_speed()` then comes out at 820 m/s rather than 1500,
capped by `1.3 · velocity_dispersion`. `chem` tracks phase fractions, but
nothing at the `Matter` level knows a liquid is not a gas. For a beach this
matters, and it is a gap the audit did not name.

### 4.3 What this exposes about the plan

The phases as first written delivered solids, creatures, terrain, contact,
construction and sessions, and **fluid at play resolution appeared in none of
them.** That was an omission rather than a deferral — "water behaves like water" is not a garnish on
the stated play space, it is most of a beach, a river, a flooded compartment, a
wake and a rainstorm.

The subsystem, stated once so it can be sequenced rather than rediscovered:

1. **A free surface.** The audit's own row: a node holds one `Mixture` with
   phase fractions, not an interface. No water level, no buoyancy, no sloshing.
   Everything else here waits on it.
2. **A liquid equation of state**, so pressure and sound speed stop being ideal
   gas for condensed matter.
3. **Weakly-compressible or projection-based SPH**, which is §3.7's open
   question with this scene as its trigger.
4. **Boundary conditions against arbitrary geometry**, so water meets a carved
   channel rather than a sphere.
5. **Multi-resolution transport** — particles crossing between a coarse ocean
   and a fine channel. `sample` and `summarise` are exactly the right
   abstraction, conserving across a scale change by construction; they have
   simply never been asked to do it for a flowing fluid crossing a live
   interface each frame.

And one thing that is not fluid at all: **terrain has to be editable.** A carve
is a persisted deviation over a derived base — pinned, because it was touched.
That is the same object as D9's build log seen from the other side, additive
there and subtractive here, so it should be one mechanism and not two.

**Both are now sequenced.** The five fluid pieces are Phase 3, immediately after
ground and before creatures; editable terrain moves into Phase 2 alongside the
surface it edits. The beach test is Phase 3's completion criterion, which is the
point of having written it down.

---

## 5. What the world keeps

The beach test asks whether the channel *fills*. This asks what happens to it
afterwards, and it turns out to be the same question as one of the engine's
founding axioms — asked at a resolution the axiom was not written for.

### 5.1 The binary that should be a rate

"Detail exists where something is happening" is implemented as: most of the tree
is deleted every frame, anything regenerable regenerates bit-identically, and
anything *touched* is pinned and persisted instead. `mixing_time` says it
literally —

```rust
if n.pinned || n.morphology.is_some() || n.topology.is_some() {
    return f64::INFINITY;
}
```

— so a user's fingerprint is exempt from forgetting **outright and forever**. At
the scales the axiom was written for that is right: a tree that lost a branch has
lost it, and a broken thing does not un-break. At play resolution it is wrong,
because it makes a footprint as permanent as a felled trunk. `DESIGN.md` already
flags the mixing time as "a discriminator, not a derivation"; this is where the
discrimination stops being adequate.

**The fix is to make the exemption a decay.** A deviation carries an amplitude,
and the amplitude falls at a rate the physics sets. When it drops below the
resolution anything could observe it at, the deviation is dropped and the node
goes back to being purely regenerable. Forgetting stops being the deletion of a
memory and becomes **the deviation reaching zero**, which is a different and much
better-behaved thing.

### 5.2 Erosion is not a mechanic, it is the deviation's own physics

A footprint fills in because grains move under wind, rain, gravity and traffic.
That rate derives from the material's cohesion, the local flux, and the feature's
own geometry — a sharp narrow notch has a steeper gradient than a broad shallow
dish, drives more flux, and goes faster. Nothing needs a table: wet sand at the
tide line goes in hours, dry sand above it in days, a rut in clay in months, a
scar in granite in millennia, and all four are the same expression with different
material and flux.

It is D7's pattern applied to terrain rather than to gaits: derive the rate once
for a `(material, flux, geometry)` tuple, store it, re-derive when the tuple
changes. A channel that the tide starts reaching erodes at a different rate that
afternoon, and nobody wrote a rule about tides.

**The honesty test for this, when it is built:** can one expression produce the
granite case and the wet-sand case with only material and flux differing? If it
needs a per-material correction, it is a table wearing a derivation's clothes and
it has failed.

### 5.3 Consequence is promotion to baseline, and the criterion is `summarise`

The sharp half of the beach example is that a squiggle should vanish while a
channel that breaks a river through should permanently change the region — and
that even then, the *fine* geometry can still be dropped.

That is not two behaviours. It is one, and the engine's central operation already
expresses it:

> **A deviation is absorbed into the baseline exactly when `summarise` of the
> edited node differs from `summarise` of the unedited node by more than the
> coarse level can represent. Otherwise it decays.**

The squiggle does not change the patch's summarised state, so it decays. The
breakthrough changes where the water goes — drainage, mass distribution, the
patch's own terrain parameters — so it is promoted, and the fine channel geometry
is then **dropped**, because it is no longer a deviation at all. The river running
here is the new derived baseline. Permanent large-scale change, and zero
fine-detail storage.

**Say what "differs" means, because §5.8 turns on it:** the comparison is over
the *conserved set* — energy, momentum, angular momentum, charge, baryon and
lepton number — and then over whatever the coarse level represents beyond it.
Stating it as "drainage and mass distribution" was an example standing in for the
rule, and the example is what lets a reader think a deviation that walked off
with a node's mass might be droppable. It cannot be.

Two properties worth naming. The criterion is a **measurement, not a judgement** —
nothing decides which edits are important, and no edit is tagged as significant
when it is made. And the threshold is not free to choose: it must be the coarse
level's own representational resolution, in the shape `IDEMPOTENT_TOLERANCE`
already has. A tuned number here is where "no special cases" would quietly die.

### 5.4 Nothing repairs anything; forgetting does the work

An earlier draft of this section had damage lower a planned program's progress
and a labour rate raise it back, so a city healed because its inhabitants were
notionally working on it. **That was wrong, and it was wrong in the way the
axioms are meant to catch.**

`env.labour` is a number nobody derives. A "target" state encodes intent. Real
masonry does not grow back, and the only thing that ever repairs a wall is
somebody deciding to and doing work — which is *intention*, and intention is not
physical law. "Only physical law is axiom-side" therefore forbids it outright,
and building a labour economy into the engine would have been exactly the
`if is_forest` the axiom warns about, wearing an economics costume.

The alternative is not to write repair agents either. It is to notice that
**nothing has to be repaired for the damage to stop being there.**

A bullet hole is a deviation over a derived base (§5.1). While it is remembered
it is remembered exactly. When it is dropped, the node **regenerates from its
program** — and the program describes a maintained wall, because a maintained
wall is what that program builds. No process ran. Nothing was mended. The
specific damage simply stopped being remembered, and what came back is what the
program says is there.

That is the same mechanism as the beach squiggle, applied one level up, and it
costs nothing new.

**Why the program may legitimately describe a maintained wall.** This is not
teleology sneaking back in. `README.md`'s fourth idea already establishes the
category: matter that is *ordered* cannot be regenerated by max-entropy sampling
because its configuration is contingent, so for those the generator becomes a
developmental program. A planned program encodes intent because **intent is what
made the thing** — `Program::is_planned` means "target known, progress-driven",
and that was accepted long before any of this. Using it for the building's
maintained condition is consistent with what it already is, not a new exception.

**Three outcomes, from one rule, with no agent anywhere:**

| what happened | does `summarise` notice? | outcome |
|---|---|---|
| bullet holes in a wall | no | dropped when convenient; regenerates as the program describes |
| the roof has come down | yes | promoted to baseline; permanent, and the program now generates a wreck |
| nobody lives here | — | the settlement's own state says derelict, so regeneration gives a ruin |

The last row matters: **an abandoned place does not heal**, and not because a
rule says so. Its program generates what it generates, and a settlement with no
population generates a ruin. That falls out rather than being enforced.

**And the hand-wave was the axiom all along.** "It got fixed off screen" is not a
fudge under this scheme — it is precisely correct, and it is the engine's ninth
idea. An unobserved quantity has no committed value. A wall nobody watched has no
fact of the matter about its bullet holes, so it comes back as its program
describes it. Observation-as-commitment is what makes that honest rather than
convenient, and `README.md` is blunt that this is not an approximation of quantum
mechanics but the thing itself.

**The rule that keeps it from being a visible pop:** never drop a deviation that
is currently resolved for an observer. That is already what `pinned` means; it
just has to be honoured here too. The consequence is real and is the next
section.

### 5.5 The busy square is a budget problem, not a semantics problem

The worry is that a place under constant observation never gets an unobserved
window in which to tidy up, so damage accumulates until something forces a
repair.

**Decay and repair run on the world clock, not on the absence of an audience.**
Wind, rain, feet and stallholders do not care whether anyone is watching. What
observation changes is not the rate — it is whether the intermediate detail has
to be *stored* rather than regenerated.

**So a watched place does not recover while it is watched.** Damage in a busy
market accumulates: it cannot be dropped, because dropping it in front of
somebody is the pop that §5.4 forbids. It goes from pristine, to pristine with
five holes, to pocked — and never back to pristine until it is out of sight. That
is the honest behaviour, and it is also the right one: a place under constant use
*should* look used.

So it is a budget question, and it has a bound. The stored deviation count in a
busy square is arrival rate times mean lifetime: footprints arriving ten a second
and lasting five minutes is three thousand stored, which is nothing; the same
footprints lasting a week is six million, which is not affordable. **The decay
rate is what bounds the memory**, it is derived rather than chosen, and whether a
given square is affordable is a measurement rather than a hope. It converges
instead of growing without bound, which is the property the design actually
needed.

### 5.6 When it does not converge, summarise rather than repair

If arrival outruns decay, something has to give, and D1 already says what: detail
gives way. But the right way is not deletion, and not a repair effect either —
it is §5.3 again, with many small edits instead of one consequential one.

**A thousand footprints summarise into trampled ground**: compacted, grass gone,
different roughness, different cohesion, different erosion rate. Which is exactly
what a desire path is. The individual prints are dropped, the coarse consequence
is promoted to baseline, and the world keeps the information that mattered
instead of the information that was merely recorded. It is physically right, it
is visually right, and it is the same operation as the river.

That is strictly better than popping an object back to a simplified state behind
an effect, because nothing is lost that a person would notice was lost.

The structural case is the same sentence: **a hundred bullet holes summarise
into a pocked wall.** The individual holes go; the surface's condition — its
roughness, its weathering, a little less section where the section was lost —
becomes part of the member's derived description and regenerates from it. And
note what that means at the coarse level: "a wall patched a hundred times" and "a
wall weathered by a hundred impacts" are the same description. At the resolution
where the individual holes stopped being tracked, the distinction between having
been repaired and having been worn stops existing. Nobody has to choose which
one it was.

**Where a visible repair is legitimate** is when a mind actually performs one. A
stallholder mending their stall is an actor issuing acts (D8) — real, observable,
and worth rendering because it is happening, not as cover for bookkeeping. That
is the only repair in the design, it is never required for the world to look
lived-in, and nothing has to be written for a city to avoid rotting.

### 5.7 Two representations, which turn out to be one scale transform

**Decided: a field delta for terrain, an edit list for structures.**

Terrain is always many small edits over a continuous surface, and a field
superposes for free — a thousand footprints cost what one costs, and storage
scales with area rather than with how eventful the ground has been. A structure
is a discrete thing built of discrete acts, the first bullet hole is
individually memorable, and an edit list is both far cheaper when edits are
sparse and the only representation in which "this member was placed by that
person" survives at all.

§4.3 said these "should be one mechanism and not two", and §5.7's first draft
retracted that as an over-claim. Both were wrong, and the truth is better than
either: **they are the fine and coarse ends of one scale transform.**

An edit list *summarises into* a field delta. A field delta *samples into*
plausible individual edits. That is `summarise` and `sample`, which is the
engine's central pair, applied to deviations instead of to bodies — and it is
what §5.6 is actually describing when a hundred bullet holes become a pocked
wall. The list is the resolved representation, the field is the coarse one, and
the merge nobody wanted to write is the scale transform that already exists.

So the two are not two mechanisms and not one representation. They are one
mechanism with two representations and a transform between them, which is
exactly the shape everything else in this engine has. Terrain mostly lives at
the field end because ground is never eventful one grain at a time; structures
start at the list end and slide toward the field as damage accumulates.

The conservation requirement comes along for free and must be honoured: the
transform has to conserve what a deviation *means* — lost section stays lost
section when a hundred holes become pocking — in the same way `summarise` already
conserves the conserved tuple.

---

### 5.8 Forgetting must be conservative

**The exploit.** Take a window out of a building. Walk away. Nobody is watching,
so the deviation is dropped, the wall regenerates from its program, and the
window is back tomorrow — while the glass is still in your pack. Repeat. Free
material, for ever.

**This is not a balance problem, and treating it as one is the mistake.** It is a
conservation bug: mass and energy appear from nowhere. The engine has the
strongest possible machinery for that already — `summarise(sample(m)) == m` on
the conserved set, held to 5.8 × 10⁻¹⁶ over 126 configurations — and the fix is
to hold forgetting to the same standard that materialisation is held to.

> **A deviation may be dropped only if dropping it leaves the conserved tuple
> unchanged.** Forgetting is not permitted to be a source or a sink.

§5.3 already said this and I stated it badly, by example — "drainage, mass
distribution" — rather than as what it is. The criterion is the conserved set:
energy, momentum, angular momentum, charge, baryon number, lepton number. Mass
walking out of a node in somebody's pack is the most obvious possible change to
it, which is why the exploit is caught by the rule that was already written
rather than needing a new one.

**The cases, walked through.**

| what the actor did | conserved tuple | outcome |
|---|---|---|
| takes the window away | baryon and rest-mass energy both drop | **promoted to baseline** — the window stays gone |
| breaks it, cullet on the floor | unchanged; the glass never left | droppable — re-glazed from its own cullet, and mass-neutral |
| dismantles it brick by brick | each removal changes it | every one promoted; permanently dismantled |
| swaps in an equal mass of sand | see below | **not caught today** — needs embodied energy derived rather than carried |

**The clever version deserves its own paragraph**, because it is the one that
nearly works. Match the mass. Glass is SiO₂ and sand is SiO₂, so against eight
lumped `CoarseElement` buckets their composition is very nearly identical, and
`DESIGN.md` is already honest that the elemental account is too coarse for
molecular diversity. Mass matches, composition matches, so the deviation looks
droppable and the swap is free.

It is *not* caught by the mechanism the first draft of this section claimed, and
the claim was made in exactly the way this project's first trap describes:
asserted from a plausible reading of a doc comment, and disproved by three
minutes of grep.

`Matter::chemical_energy` is indeed embodied energy — its own doc says "a steel
frame holds its embodied energy. Destroying the structure releases it" — and
`non_rest_energy` does fold it into the total that `Conserved::energy` reports.
All true, and none of it sufficient, because of how it crosses a scale
transition:

```rust
// state.rs, summarise()
// Not knowable from the children alone; the caller reinstates these.
external_potential: 0.0,
chemical_energy: 0.0,
```
```rust
// sampler.rs
back.chemical_energy = matter.chemical_energy;
```

**Embodied energy is carried, not derived.** `summarise` zeroes it and the caller
copies the stored value back unchanged, so it is opaque to the scale transform.
The §5.3 criterion compares `summarise` of the edited node against `summarise` of
the unedited one — and embodied energy is zero on both sides of that comparison.
It cannot catch anything. The swap works, and the exploit is open.

**The fix is to stop carrying it and start deriving it, for the case where it
can be.** The comment is right that a *body list* cannot report a structure's
embodied energy: it lives in the arrangement, not in the parts. But a structure
does not only have a body list. It has a `topology` — members, joints, materials
as data — and from those the embodied energy **is** derivable: it is the sum over
members of what that material cost to form into that member.

So for ordered matter, `chemical_energy` becomes a derived quantity with a stored
shortcut, which is the third axiom exactly: derive it once from the topology,
store it, re-derive when the topology changes. Then removing a window removes a
member, the derived total drops by that member's share, `summarise` differs, and
the deviation is promoted. **Nothing has to remember to debit anything**, which
is the property worth having — the alternative is a rule that every removal path
must observe, which is the same shape as the `reparent` side-table trap and would
fail the same way.

With that in place the substitution argument works as originally described, and
its closure is a satisfying one: to defeat the check you would need a substitute
matching the window in mass, in eight-bucket composition, *and* in embodied
energy — at which point you have not forged glass, you have made some. The
economy conserves because the physics does.

**A second hole in the same place, and it is the bigger one.** All of this checks
the node the window left. It says nothing about where the glass went. If taking
the window appends a row to an item list, then the world lost mass, the check
fires correctly, and the glass in the pack is matter that no node accounts for —
conservation broken at the boundary, with every node-side test still green.

So: **there is no inventory. There is matter that happens to be carried.** A
stowed object is matter inside the actor's node; a held one is jointed to them
(D7's grasp), so its mass is genuinely in the load path and a creature can be
overloaded. Taking the window is then a transfer between two nodes, checked on
both sides by machinery that already exists, and an actor's carrying capacity
stops being a number in a design document and becomes what their body can bear.

**Observation is irrelevant to this, which is what makes it un-gameable.** A
removal is never droppable, watched or not. There is no looking away, no waiting
for the node to coarsen, no server-restart trick — the test is on the conserved
tuple, not on who was present.

**Harvest is not theft, and the same account says why.** Fruit taken from a tree
regrows, and that is legitimate income rather than duplication, because growth is
a process with a real input: the sun's flux enters the energy balance and pays
for the fruit. Regeneration-from-forgetting has no input, so it may not be a
source. One account separates renewable resources from money-printing, and
nobody has to write a list of which is which.

**Two things that were already blocking neighbouring exploits**, worth naming so
they are not re-solved. Re-rolling a sample by coarsening and refining until the
contents are favourable does not work: regeneration is deterministic in
`(matter, spec, world_seed, path_key, epoch)`, and `epoch` moves only on a
recorded interaction. And a measured quantity is committed to the ledger and
never re-sampled, so measuring-then-rerolling is not available either.

**The cost is not a problem, because promotion is absorption rather than
retention.** A statue that souvenir-hunters chip at for a year does not
accumulate a year of stored chips; it becomes a smaller, more worn statue, folded
into the baseline by §5.6's summarising. Storage tracks how much genuine
modification happened, not how many times somebody poked it.

**Why not the obvious fixes.** Ownership flags on windows, loot cooldowns, "this
object cannot be taken twice" — each is a special case the engine has to be told,
each is a rule an ingenious player will find the edge of, and collectively they
are the table of melting points the axioms exist to prevent. Conservation is a
law, it is already enforced to machine epsilon, and there is no edge to find.

**One implementation trap, straight out of the code's own warning.** The
droppability check must difference `non_rest_energy()`, **never**
`total_energy()`. Rest mass exceeds every other term by roughly 10¹⁶, so a
difference of totals leaves about seven significant digits and none of them
reliable — `state.rs` says exactly this, about exactly this arithmetic. A
droppability test written against `total_energy()` would find every deviation
conservative, pass its tests, and leave the exploit wide open.

**And the weakness to hold onto.** Eight lumped elements cannot tell silica from
silica, so the sand-for-glass case rests on embodied energy and on nothing else.
That means two things must both be true, and neither is true today: embodied
energy must be **derived from the topology** rather than carried opaquely through
the scale transform, and every program must **book what building cost** so there
is something to derive from. A `Program::Tower` that raises a wall without
recording what raising it cost leaves that wall forgeable no matter how good the
check is. Both are testable invariants rather than hopes, and both belong in the
suite before any player can reach a building.

There is also a precision question, which is the same one that bites everywhere
else in this engine: one window's embodied energy against a whole building's is a
small number differenced from a large one. `non_rest_energy` removes the 10¹⁶ of
rest mass, but a building's own embodied energy may still dominate a single
member's by five or six orders. Whether the check has the digits it needs is a
measurement, not an assumption.

---

### 5.9 Which window, and what happens at the sixty-fifth edit

Two questions that turn out to have one good answer and one bad one, and the bad
one is a defect in code that ships today rather than a gap in this plan.

**Does the taken window always sample the same? Yes, and the design already
guarantees it.** A structure's geometry is not sampled statistically. `morph.rs`
states it twice — `structure = program(genome, age, events)` in the module doc,
and "Pure in `(genome, age, built, progress, events)`" on the generator — and
**`epoch` is not in that tuple.** That matters more than it looks: taking a
window is a recorded interaction, so it bumps the node's `epoch`, and everything
sampled statistically in that node redraws. If structures were epoch-seeded, a
player removing one window would watch the other nine jump. They do not, because
a structure's variation comes from a genome derived once from the path key and
then stored.

The taken window itself is no longer generated by anything. It is matter with an
`EntityId` (D2), carried as matter rather than as an item row (§5.8), and its
properties were committed at the moment somebody took it — which is the ledger
doing exactly its job. Before that interaction there was no fact of the matter
about that window's specifics; afterwards there is, permanently.

**Which of ten windows is missing is `Event { site }`, and that mechanism already
exists.** `morph::Event` is documented as "the generalisation of `Node::epoch`,
which can only say 'the procedural detail is stale' and not what changed… so the
deviations are logged and replayed rather than discarded", `EventKind::Severed`
is "a limb removed — pruned, snapped, demolished", and `Skeleton::site` is
described as "a program-stable name for this part, so an event can refer to it
and mean the same thing after the structure is regenerated" — which is precisely
the property the question asks about. The edit list §5.7 chose is not a thing to
invent. Whether it is a thing already built is the next paragraph, and the answer
is only partly.

**Then it was measured, and it is worse than reading the code suggested — in a
different place.** Severing sites one at a time on each program and re-rendering
after every one:

```text
program      intact   after 64 severed   after 65   severed sites present again
tree            512          0 parts      512 parts        33 of 65
coral           512          0 parts      512 parts        33 of 65
tower           296        296 parts      296 parts        64 of 64  (from the first)
wall            504        504 parts      504 parts        64 of 64  (from the first)
terrain         484        484 parts      484 parts        64 of 64  (from the first)
settlement      256        256 parts      256 parts        64 of 64  (from the first)
```

**Two separate defects, and the second is the one that matters here.**

**One: the log cap resurrects a destroyed structure.** For the branching programs
the event log works exactly as documented up to the cap — a tree severed at the
trunk goes to zero parts and stays there. Then the sixty-fifth event drops the
oldest thirty-two, and the tree **comes back whole**, all 512 parts, with only
the 33 most recent severances still honoured. A felled tree returns. `built` was
decremented each time, so the mass is right and the geometry is not, which is a
structure contradicting its own conserved state. This is live today and has
nothing to do with the play space.

**Two: four of the six programs never honour a severance at all.** Tower, wall,
terrain and settlement return exactly the same part count from the first
severance onward, with every severed site still present. The skip —

```rust
// A severed limb and everything above it is simply absent.
if self.events.iter().any(|e| e.kind == Severed && e.site == s.id) { continue; }
```

— lives in `render_branching`, which only `Tree` and `Coral` call. The planned
programs, which are exactly what buildings are made of, do not consult the event
log when generating geometry. So a demolished wall section does not survive
sixty-four edits; **it does not survive one.** The window is back on the next
regeneration, with no cap involved.

That makes the earlier sections' confidence misplaced in a specific way worth
recording: `morph.rs` documents that "deviations are logged and replayed rather
than discarded", and for two-thirds of the programs the logging happens and the
replay does not. The doc describes an intent that one renderer implements. §5.7
chose the edit list for structures on the strength of that doc; the choice is
still right, but it is a thing to finish rather than a thing to use.

**The fix is the rule from §5.8, applied one level down.** Ask of each event the
same question asked of each deviation: *did it change the conserved tuple?*

- **`Severed` changed it** — mass left the structure. It may never be compacted
  away. Either the event stays, or the removal is **promoted into the program's
  own description** so the program stops generating that member at all, which is
  §5.3's promotion done properly and keeps the log bounded without anything
  coming back.
- **`Damaged` and `Suppressed` did not** — nothing left. These are exactly the
  events that may merge into a surface field, which is §5.6 and §5.7 already.

So the compaction policy stops being a count and becomes the conservation test
that governs everything else here, and the split falls exactly where §5.7 put it:
**damage merges into a field, removal does not.** Raising `MAX_EVENTS` would only
move the wall; nothing derives 64, and nothing would derive a larger number
either.

Worth noting in passing that `genome[7]` is currently doing double duty as a sink
for accumulated history, which makes a per-instance variation slot mean two
things at once. Whatever replaces compaction should leave the genome alone.

---

### 5.10 How a region remembers

D12 says a grown shape depends on the conditions it grew in, and asserts those
conditions re-derive rather than needing storage. That is the load-bearing claim
and it deserves spelling out, because "the environment it grew in" sounds like a
time series and a time series per region would be ruinous.

**A region is a node. There is no other kind of thing.** A terrain patch, a
settlement, a continent and a planet are all nodes at different tiers, and
`Tier::containing(metres)` already maps a size to a regime. Adding a `Region`
type would be precisely the special case the axioms forbid, and it would
immediately need its own persistence, its own side table and its own line in
`reparent`.

**But a node holds an instant, not a history.** `Matter` is state at a moment.
`causal::History` does exist and must not be repurposed for this: it is a short
fixed-capacity ring of `Moment { position, velocity, mass, luminosity,
temperature }` whose job is retarded-light lookups into the past light cone. It
forgets by design and it carries none of the quantities a climate needs.
Conflating "where was this three seconds ago, for light travel" with "was there a
drought in 1340" would wreck both.

**Most of a climate needs no storage, because it re-derives.** Insolation from
latitude, axial tilt and orbital phase. Prevailing wind from rotation and the
equator-to-pole thermal gradient. Rainfall from orography and moisture
advection. Each is a pure function of place and time, and place and time descend
from the world seed — so the baseline is *computed* rather than remembered, and
the same node always had the same weather. This is the same argument that keeps
a regenerated terrain patch bit-identical, applied to the weather over it.

**So what is stored is only what deviated from that** — which is §5's model one
level up, and the shape proposed in the question is the right one. A drought is a
bounded episode with a start, an end, a magnitude and an extent: an *event*. A
slow warming is a *field*. Those are §5.7's two representations, so this is the
same mechanism rather than a parallel one, including the transform between them:
individual storms summarise into "that decade ran wetter than derived" exactly as
a hundred bullet holes summarise into a pocked wall.

**The event picks its own node, by extent.** Nobody declares a region. A
disturbance has a size, and it is stored on the coarsest node whose radius covers
it — the same move `Tier::containing` makes when it turns a size into a tier. A
continental drought lands on a continental node; a late frost lands on one patch.
Finer nodes read up the ancestor chain, which `environment_at` already does for a
parent's luminosity.

That is where the economy comes from: **one drought is stored once and read by a
million trees.** Storage is O(events), not O(events × individuals).

**And a tree does not need the history — it needs the integrals.** D12's response
surface takes conditions as its input, so what a tree consumes over its life is
integrated light, mean water, crowding, peak wind: a handful of accumulating
numbers. Because the baseline re-derives, the baseline's integrals re-derive too,
analytically or by cheap quadrature over a known function. Only the deviations'
contribution has to be accumulated, and that is a sum over a short event list.
**There is no time series anywhere in this.**

**One detail that will break a transplant if it is got wrong:** the integrals
belong to the *individual* and the climate belongs to the *region*. They must be
two different homes. A tree moved to another continent should meet its new
region's weather from now on while keeping the growth it has already done — and
if its accumulated integrals lived on the region node, `reparent` would silently
rewrite its history. So integrals ride with the entity (D2), climate deviations
ride with the region node, and `reparent` moves one and not the other.

**Whether old events can be compacted is decided by D12's probe, not separately.**
A drought can eventually be folded into the baseline — a decade that ran drier is
a shift in the derived baseline for that decade — which preserves the integral
and discards the ordering. That is safe *if and only if* growth is Markovian in
the integrals, which is exactly the question §D12 already sends to a probe. One
measurement decides both whether the response surface exists and whether regional
history can be compacted at all.

And §5.9's lesson governs the compaction when it happens: **it must be
conservative in the quantity that matters.** There, folding severances into a
mean magnitude lost the fact of the severance and the structure regrew. Here the
quantity is the integral, so the invariant is that compaction preserves it. A
compaction that changes what a tree would have grown into is the same bug wearing
a climate's clothes.

---

## 6. What has not been measured

Per `CLAUDE.md`'s first trap, these are stated as unmeasured rather than
assumed, and each has a probe in Phase 0. The substeps-per-tier question that
§3.4 answers was one of these and has been measured; the numbers there are from
a scratch probe that Phase 0 should commit properly rather than from arithmetic.

| Question | Why it matters | Probe |
|---|---|---|
| What does a neighbour query cost at 10³–10⁵ contents? | D3 is on every frame's critical path. If it is expensive, the frame budget arithmetic in `PERFORMANCE.md` changes. | Build the index over the existing scenarios and time it. |
| ~~Do segments of a walking limb exceed 0.1 rad of chord rotation?~~ | **Answered: no, and D5 is confirmed more strongly than it claimed.** The member ruptures at 700 N having travelled 1.3% of its length, with chord rotation only 0.018. A limb has no elastic path to a stride at all, so the large rotation must live in a joint between substructures and corotational elements are *off* the critical path — they would permit a bend the material does not. | done — `tests/probes.rs` |
| Does slaving the parent body every frame preserve `summarise(sample(m)) == m`? | D4 writes into the conserved set every frame. `IDEMPOTENT_TOLERANCE` is the contract. | Promote, run, coarsen, compare against the existing consistency harness. |
| How long does a gait optimisation take, and does it converge? | D7's shortcut is only a shortcut if deriving it is rare and bounded. | Solve one quadruped gait offline and time it. |
| ~~How large is the checkpoint for an interactive subtree?~~ | **Answered: 181 bytes per pinned body.** Only pinned nodes write bodies, and an interactive subtree is pinned by definition, so 10⁵ bodies is ~18 MB a checkpoint. Regenerable detail costs 45× less. | done — `tests/probes.rs` |
| ~~Does metre-scale bulk fluid actually hurt?~~ | Answered by §4: yes, for anything at play scale — a 5 cm channel is two orders below the floor. §3.7's option (1) is not the end of the matter, and option (2) is what Phase 3 adopts. | — |
| Does a busy square's stored-deviation count converge, and to what? | §5.5 argues arrival rate times mean lifetime is a bound rather than a hope. If the number is millions, the decay rate is wrong or §5.6's summarising is load-bearing much earlier than expected. | Simulate arrivals at a plausible footfall against a derived sand/paving erosion rate and count what is held. |
| ~~Can one genome give a forest tree and an open-grown tree?~~ | **Answered: no, exactly as D12 predicted.** Three light fields grow masses differing 500-fold; at equal age and mass the skeletons are bit-identical. Conditions reach form only through how much mass they grew. | done — `tests/probes.rs` |
| ~~Is growth Markovian in its conditions?~~ | **Answered: no. D12's stated mechanism does not work.** Same light integral, opposite ordering: rising +29.7%, falling −26.3% in final mass against flat. A per-species surface over *integrated* conditions cannot reproduce a life. **This needs review before Phase 4** — see the note below. | done — `tests/probes.rs` |
| Which conditions actually move the shape? | Decides the response surface's dimensionality, and it grows badly. Six axes at five samples is ~15,600 runs per species. | Vary each condition alone and rank by how much the gross form moves. |
| Can a grown shape stay regenerable when it depends on a history? | D12's fatal risk. If the environment history will not summarise into ~200 bytes, shape stops regenerating bit-identically and the founding invariant breaks. | Grow a tree, summarise its environment history, regrow from the summary, and diff the skeletons. |
| ~~Can an undesigned genome make something coherent?~~ | **Answered: yes, 600 of 600.** A hundred randomised genomes per program, six programs, all finite and non-degenerate. Caveat: this tests the robustness of the *existing* genome slots, not of D11's proposed unified law. | done — `tests/probes.rs` |
| ~~At what edit count does a structure start regrowing removed parts?~~ | **Answered: 65 for the branching programs, 1 for the other four.** Tree and coral suppress all 64 then resurrect 32 on the sixty-fifth; tower, wall, terrain and settlement never suppress a severance at all. | done — `tests/probes.rs` |
| Can embodied energy be derived from a topology, and does removing one member show up? | §5.8's fix depends on it entirely. Carried opaquely, as today, the check is blind. | Derive it for a walled structure, remove one member, and difference `non_rest_energy` against the unedited regeneration. |
| ~~Does the embodied-energy check have the digits?~~ | **Answered: yes, comfortably — 12.0 to 13.5 digits** after differencing a whole structure to find its smallest member, against the ~6 the check needs. The precision worry is retired. | done — `tests/probes.rs` |
| ~~Does every program book the embodied energy of what it builds?~~ | **Answered: yes.** Tree and coral 1.7e7 J/kg, tower, wall and settlement 2.5e6 J/kg. Terrain books zero, correctly — rock was not manufactured, and that contrast is what §5.8 relies on. | done — `tests/probes.rs` |
| Does an edit list summarise into a field without losing what the damage meant? | §5.7 makes this the same guarantee `summarise` already carries for the conserved tuple. If lost section does not survive the merge, a wall repairs itself by being forgotten in the wrong sense. | Shoot a member a hundred times, summarise, and compare the member's strength against the un-summarised case. |
| Can one erosion expression give both granite and wet sand? | §5.2's honesty test. If it needs a per-material correction it is a table, and the axiom is broken. | Derive the rate for four materials from cohesion and flux alone and compare against observed rates. |
| Where does the §3.4 crossover fall on *terrestrial* material? | The measured table used the galaxy scenario, whose `Continuum` gas is hot and fast. A room is not that. | Re-run the same drill on a planetary-surface scenario once Phase 2 exists. |

### The one Phase 0 result that changes a decision

**D12's mechanism, as written, does not work, and this needs review before
Phase 4 rather than a unilateral rewrite.**

D12 says the growth response is derived once per species over a sampled space of
conditions and stored, then evaluated per individual. The probe says the space it
would be sampled over is the wrong one: growth is *not* Markovian in the
integrated conditions. Same light integral, opposite ordering, over a hundred
simulated years —

```text
history          final mass    against flat
flat 400          3.594e4 kg        —
rising 80->720    4.663e4 kg     +29.7%
falling 720->80   2.648e4 kg     -26.3%
```

Ordering is not a correction here, it is a third of the answer. A surface indexed
by "how much light did it get" cannot tell a tree that was shaded young from one
shaded old, and those are different trees.

What this does *not* kill is the axiom D12 rests on. Axiom three's second clause
already says a stored rule "stays a function of its conditions and never a frozen
outcome" — and an outcome surface over integrated conditions is exactly a frozen
outcome. So the axiom forbade what the probe disproved, which is a good sign for
the axiom and a bad one for the paragraph.

The obvious repair is to store a **rate law over the current state** rather than
an outcome over an integral: a surface mapping `(current size, current
conditions)` to an increment, applied step by step. That is Markovian in the
*state*, which is a different and much weaker claim than being Markovian in the
integral, and it is what `advance` already does. It would keep D12's economy —
derive once per species, evaluate per individual — while dropping the part the
probe just falsified.

**That repair is not made here.** It changes a decision, and the plan says
decisions get made before code rather than during it.

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

## 7. The order of work

Each phase ends with a test that fails today. A phase is not done because its
code exists.

**Phase 0 — Probes.** The measurements above, each as a test that demonstrates
the thing it claims — including committing the substeps-per-tier probe that
§3.4 reports, since it is currently a scratch run. Nothing is designed further
until they are numbers. *Done when:* `PERFORMANCE.md` carries the new measured
rows and D5 is either confirmed or replaced.

**Phase 1 — The primitives.** `EntityId` and the side-table rekey (D2);
`Neighbourhood`, boundary exchange and contact (D3); the promoted child made
authoritative and force-bearing (D4); tier revisited on size change; the spread
measurement; `PaceMode::Fixed(1.0)` as what a world is, and `G_EARTH` deleted in
favour of derived g. Plus the two tier corrections from §3: solver dispatch that
reads ordered-versus-disordered state as well as size (§3.3), and a reported
resolution floor so a node dropping to its ensemble says so instead of doing it
quietly (§3.7). *Done when:* two promoted vehicles collide and rebound with
restitution derived from their materials; a hot node beside a cold one
equilibrates without either being told the other exists; a branch lands on the
next tree; one node holds both a structure and loose contents and steps both
correctly in one pass; and nothing in the existing suite regresses.

**Phase 2 — Ground.** Cubed-sphere parameterisation, patches as `Program::Terrain`
nodes, refinement and coarsening on approach, handoff by `reparent`, planetary
gravity, and terrain as an **editable deviation over a derived base** (§4.3) —
carving shares D9's build log's *concept* but not its representation (§5.7),
and both are settled here — along with the decay machinery of §5: a deviation
with an amplitude, a derived erosion rate, and promotion to baseline when
`summarise` says the coarse level noticed. Plus D11's first half — `density`,
`energy_density`, `substrate`, `design_flow` and `maintenance` moved off
`Program` onto the material and the measured environment — because a derived
erosion rate cannot coexist with a tabulated one. The surface representation is designed with
Phase 3 as a named consumer, because a surface that cannot hold a puddle is one
that gets rebuilt. *Done when:* an observer descends from orbit to a square metre of any
planet in any scenario, travels ten kilometres across patch boundaries, and the
terrain behind them regenerates bit-identically. **And:** a squiggle drawn in
sand is gone by the next tide while a channel that redirects drainage is still
there a year later, with neither having been tagged as important when it was
made.

**Phase 3 — Water.** The five pieces §4.3 names, in dependency order: a **free
surface**, so a node holds an interface and not only phase fractions; a **liquid
equation of state**, so `pressure()` stops returning 4×10⁸ Pa for a bucket of
water; **weakly-compressible SPH**, which is what makes centimetre flow fit a
frame at all; **boundary conditions against arbitrary geometry**, so water meets
a carved channel; and **multi-resolution transport**, particles crossing between
a coarse ocean and a fine channel through `sample` and `summarise`, which
conserve across a scale change by construction and have simply never been asked
to do it for a flowing fluid.

It sits here, immediately after ground and before creatures, for one reason: a
beach, a river and rain are most of what makes a planet feel like a place, and
water is the first thing anybody tests. Designing it against a surface that
exists — rather than retrofitting it to one built without a consumer — is the
whole argument for the position. It is also the largest single subsystem in this
document, and putting it third is a deliberate acceptance of that cost.

*Done when:* **the beach test passes.** An actor carves a 5 cm channel in wet
sand and the swash runs through it, at one second per second, at roughly 1 cm
cells — which §4.2 measures as `MAX_SUBSTEPS` at 512 and about 4×10⁵ particles
in the local patch. The ocean beyond the patch stays coarse, and nothing between
the two loses mass.

**Phase 4 — Bodies.** D11's habit refactor first, so creatures arrive as a genome
rather than a seventh species; substructuring (D5); the creature genome;
the actuation mechanism; derived-and-cached gait and grasp; local interaction
resolved inside a segment. *Done when:* a quadruped and a biped grown from two
genomes both walk, on two planets with different g, with no per-morphology code
anywhere and no enum variant naming either of them; shouldering a load visibly changes the gait; and scratching a paw does
not re-analyse the animal.

**Phase 5 — Minds.** The actor-client host outside the core crate; the
determinism constraints; checkpoints and the input log. *Done when:* a wolf
pursues something, the engine cannot distinguish it from a player, and a
recorded session replays bit-identically from seed plus log.

**Phase 6 — Making things.** The build log as a genome; player-placed members;
analysis without proportioning; break and repair. Plus §5's deviation machinery at the structure end
(§5.7): the edit list, and its summarising into a member's surface condition.
*Done when:* a player builds a bridge that holds and one that does not, and the
engine was never told which was which; a wall shot a hundred times becomes a
pocked wall rather than a hundred stored holes; and a market stall nobody
watched comes back as its program describes it while an identical stall in an
abandoned village comes back a ruin — with no labour rate, no repair process and
no agent anywhere in the path. **And the window cannot be farmed:** an actor who
removes one and leaves finds it still missing, an actor who breaks one and
leaves finds it re-glazed, and an actor who backfills the hole with an equal mass
of sand is caught by the embodied energy that was never in the sand (§5.8).

**Phase 7 — Sessions.** The server loop; two clients; prediction of one's own
avatar only; interest management under load; hindsight replay scrubbing on top
of Phase 5's checkpoints. *Done when:* two people share a world, one builds
while the other watches, and either can scrub back through a minute of it at a
thousandth speed without the world's clock moving.

---

## 8. The axioms, re-checked

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
- **Derived, with shortcuts stored.** This axiom was *extended* by D12 rather
  than merely obeyed by it: a shortcut may now hold a rule and not only a
  constant, derived once for a kind and run for each individual of it. D7's
  cached gait, D7's cached grasp and D12's morphology response are three
  instances of the one pattern. The second clause is the one that bites — a
  stored rule stays a function of its conditions, so a regional drought reaches
  every tree it touches; a response surface that a drought cannot move is a
  snapshot, and D12's whole argument fails if it becomes one.
- **Detail exists where something is happening.** D6 generates a planet's surface
  only where somebody is standing, and D1's replay respects the same thing —
  which is exactly why an unobserved replay is a fresh draw and an observed one
  is exact.
- **No special cases.** D3 is the load-bearing one: contact, conduction,
  diffusion, radiation and debris landing become one primitive rather than five.
  If a later change needs a sixth mechanism to make one of them work, that is
  the signal the primitive was wrong.
