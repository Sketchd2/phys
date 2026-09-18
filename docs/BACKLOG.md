# Backlog

Deferred work, with the reason it was deferred. Not a roadmap — the roadmap is
the architecture review. This is the smaller list of things noticed while
building something else, parked deliberately rather than forgotten.

An entry earns its place by being **specific about the trigger**: what would
make it worth doing, or what will go wrong if it is not. "Would be nice" is not
a trigger.

Several entries below now have their trigger pulled and their direction chosen
by `docs/PLAY.md`, which is the plan for making the world inhabitable — the
adjacency relation, the promoted child that never feels a force, node splitting,
the `PathKey` that does two jobs, terrain and constant-gravity debris, and
structures that cannot span promoted children. The entries stay as written
because what they *measured* is still the evidence; `PLAY.md` says what is being
done about it and in what order.

---

## ~~Pace control has no honest "manual" mode~~ — done

**Noticed:** building `phys-persist` (Phase 1). **Closed:** building
`phys-bubble`, which was the third caller the trigger predicted.

`refresh_pace` runs at the top of every frame, so an assignment to `world.pace`
was silently discarded on the next one, and the only way to drive the clock by
hand was to set `paced_to = NodeIdx::NONE` first — which worked by accident of
an early return, read as a bug at every call site, and depended on a return
nothing stopped a refactor from removing.

`PaceMode::{Follow, Fixed}` makes the choice explicit, `World::pace_fixed`
is the honest counterpart to `World::pace_to`, and the two are mutually
exclusive by construction. It persists, because a fixed pace is a property of
the world rather than of the process that set it.

---

## Structures cannot span promoted children

**Noticed:** scoping the play space to vehicle scale (§09 of the review).
**Where:** `solvers/structure.rs` — `analyse`, `analyse_with` and friends all
take `bodies: &[Body]`, one node's body list. `docs/DESIGN.md` records the
matching limitation.

A 160 km vehicle at human resolution is ~5×10^14 bodies. It cannot be one body
list, so it must be a node with promoted children — decks, sections,
compartments — and then its load path spans the hierarchy where the solver
cannot see.

**The fix** is substructuring: condense each section's stiffness onto its
boundary degrees of freedom, solve the coarse frame with sections as
super-members, then solve interiors with the coarse solution as boundary
conditions. Static condensation is ordinary finite-element practice, and here it
is the same operation `summarise` already performs on mass and energy, applied to
a stiffness matrix instead.

**Trigger:** the first structure that does not fit in one node. Nothing built so
far comes close, so this is not urgent — but it is the one item in the review
that is a genuine architectural gap rather than unwritten code.

---

## `Material.name` cannot round-trip

**Noticed:** writing `persist.rs` (Phase 1).
**Where:** `topology.rs` — `Material.name` is `&'static str`.

A saved material recovers its label by matching against `Material::PRESETS` and
falls back to `"custom"`. Every one of the thirteen numbers round-trips exactly;
only the name of a material somebody *built* is lost.

**Trigger:** goes away by itself when materials become data (Track C2 of the
review), at which point the name is data too. Not worth fixing before then —
nothing in the physics reads `name`.

---

## Pacing to an unmaterialised node runs the clock 10^10x too fast

**Noticed:** measuring the volume query's neighbourhood (Phase 2 bandwidth work).
**Where:** `engine.rs` — `node_cadence`, `pace_to`. `state.rs` —
`characteristic_speed`.

`pace_to(idx)` sets the world clock from the target's characteristic time, and
`node_cadence` derives that from body speeds when the node is materialised and
from the matter when it is not. A promoted child's matter has
`momentum = ZERO` by construction — `promote` says so explicitly, and is right
to: inside a node's own frame the net momentum *is* zero, which is what a rest
frame means. But that leaves nothing for `characteristic_speed` to measure, so
an unmaterialised node reports a characteristic time far longer than it has.

Measured on a galaxy drilled to a planetary node with 128 promoted neighbours,
pacing to one of the neighbours:

```text
    unmaterialised target:   pace 2.597e11 s/frame,  worst lateness 3.6e8
    same target, refined:    pace 6.187e0  s/frame,  worst lateness 1.2
```

A factor of 4x10^10 in the world clock, decided by whether the thing being
watched happens to have been materialised yet.

This is on the client path, which is what makes it urgent rather than curious:
`pace_to` follows what is being watched, a `ViewRequest` may name any node, and
recipes now mean a client can be looking at scenery the engine never
materialised. `tests/view.rs::an_unpinned_reload_comes_back_coarse` records the
same mechanism as an oddity of reloading; it is not an oddity, it is this.

**The fix** is that a node with no bodies has no measurable internal motion, so
its cadence must come from something else — the tier's own timestep, or the
speed the parent's body list says it has — and never from a rest-frame
matter that is zero by construction. `pace_to` should also refuse a target
whose cadence it cannot measure rather than silently accepting a nonsense one.

**Trigger:** now. Anything measured against a badly paced world is measuring
the pace.

---

## ~~Crowd lateness is tier-dependent, and unexplained~~ — explained

**Closed:** measured. It is not tier-dependent and there is nothing to fix; it
is the single-clock constraint, restated.

Lateness came out at 6.8 for a crowd of promoted planetary nodes and 1.000 for
the same crowd at continuum tier, which looked like a fortyfold difference
between tiers. It is neither fortyfold nor about tiers:

```text
    planetary:  pace 1.163e3 s, worst node's cadence 1.709e2 s  ->  ratio 6.8
    continuum:  pace 1.008e-2 s, worst node's cadence 1.008e-2 s ->  ratio 1.0
```

Lateness is `pace / cadence of the fastest live node`, exactly. The pace is
taken from the node being *watched*; its promoted children have their own
cadences, and at planetary tier those happened to be 6.8x shorter while at
continuum they happened to match. Change which node is watched and the number
moves with it.

**What it does say** is that the world clock is set by what someone is looking
at, and anything live with a shorter cadence falls behind by exactly that
ratio. `refresh_pace` could instead bound the pace by the fastest *live* node's
cadence, which would guarantee lateness at or under one everywhere at the cost
of a clock that runs slower than the viewer asked for. That is a policy
decision, not a bug fix, and it wants deciding rather than defaulting — a time
bubble is the escape hatch either way.

---

## Recipes assume the client's sampler matches the server's

**Noticed:** building `view::Recipe` (Phase 2 bandwidth work).
**Where:** `view.rs` — `Recipe::build`, and `sampler.rs` underneath it.

A recipe is ~300 bytes that regenerate a node's detail on the client, and it is
only ever sent for detail the server does not itself hold. That is what makes it
safe: there is no server-side truth for a client's version to disagree with, so
two clients sampling untouched scenery a few ulps apart is a difference nobody
can observe. `build` also checks the generated mass against the matter, which
catches a blob from a different build of the engine.

What it does *not* catch is a sampler that differs subtly — a different
`libm`, a different FMA contraction, x87 excess precision — producing bodies
that conserve mass and sit in slightly different places. Today that is
cosmetic.

**Trigger:** the moment a client's regenerated detail is allowed to *matter* —
client-side hit detection against recipe scenery, say, or a client authoring
into a node it generated. At that point the sampler needs a determinism
guarantee (soft-float for the sampling path, or a server-supplied hash of the
positions), not just a mass check.

---

## A node whose bodies are 10^20 radii outside it — cause found, partly fixed

**Noticed:** the first complete debug run of the suite (the naming pass).
**Where:** surfaces at `solvers/hydro.rs::key_of`; the cause is upstream.

`tests/simultaneity.rs` has two tests that fail in a **debug** build with
`attempt to add with overflow`, and pass in release. They are not new — the same
two fail identically at `1d0d6e5`, before the pace fix, the relativity work, the
chemistry layer and the naming pass. They had never been seen because
`frames_stay_within_budget` fails first in debug on a slow container and
fail-fast stopped the run before reaching them. That masking failure has since
been explained and fixed — it was the unoptimised test profile, not the
container: `frames_stay_within_budget` measures wall-clock frame time against a
50 ms target, and unoptimised the engine took 1228 ms. It passes at
`opt-level = 2`, which is now what `[profile.test]` sets.

**The immediate cause.** `key_of` buckets a body into the neighbour grid with
`(p.x / s).floor() as i64`. Measured, at the moment it breaks:

```text
    p = (6.564e10, -2.411e10, -6.473e10) m      grid spacing s = 1e-9 m
    p.x / s = 6.6e19                            i64 tops out at 9.2e18
```

so the cast saturates to `i64::MAX` and the neighbour walk's `cx + 1` overflows.

**Why it matters in release.** There is no crash: `i64::MAX + 1` wraps to
`i64::MIN`, the lookup consults an arbitrary far-away cell, and that body's SPH
forces are computed against whatever neighbours happen to be there. Silently
wrong physics is worse than the panic.

**The real cause is upstream, and it is not what this entry first said.** The
first guess here was that something puts parent-frame coordinates into a node's
body list. Probing disproved that. The bodies start where they belong — a few
times 10^-11 m from the node origin — and are *flung* there by an unstable
integration inside a single `advance_node` call:

```text
    PROBE MD Atomic dt=1.000e-18 stable=1.000e-24 wanted=1.000e6
             substeps capped to 64 -> h=1.563e-20
             far=9.208e-12 fastest=2.811e4
```

`engine.rs` hands the Atomic node a scheduler step of 1e-18 s.
`configuration_dt` answers that the Lennard-Jones force field is stable at
1e-24 s, so a million substeps are wanted. `.clamp(1, 64)` then caps that
silently, and the node integrates at 1.6x10^4 times the stable step. The
comment directly above the clamp says what happens next — "a molecular system
handed a step longer than its own vibrational period does not integrate
inaccurately, it detonates" — and then the code does it anyway. Bodies leave at
2.8x10^4 m/s and keep going; twenty tiers later the grid index overflows.

**Three fixes, and they are separate.**

1. **The clamp must stop lying.** A ceiling on substeps is reasonable — 10^6
   substeps in one frame is not affordable — but capping and continuing is not.
   The node should integrate the span it *can* integrate stably and fall
   behind, which is the engine's own philosophy everywhere else: nodes fall
   behind honestly, and their lateness says so. Silently running an unstable
   step is the one thing the scheduler is built to avoid.

2. **The grid should saturate, not wrap.** A bucket index that cannot be
   represented should clamp, whatever put the body there. Release currently
   wraps `i64::MAX + 1` to `i64::MIN`, consults an arbitrary far cell, and
   computes that body's SPH forces against whatever is in it. Silently wrong
   physics is worse than the panic debug gives.

3. **A spread check would have caught it at the source.** A `debug_assert`
   that a node's bodies lie within some multiple of its own radius fails on the
   first frame after the detonation, not twenty tiers later inside a hash. This
   is the same measurement node splitting needs, so the two share a detector.

Note also that `MdParams::default().cutoff` is a fixed 1e-9 m and does not
scale with node radius, unlike gravity's softening and hydro's `h` (both
`radius / count^(1/3) * k`). That is a separate latent inconsistency.

**Closed** by all three. The clamp now takes `min(dt / substeps, stable)`, so a
node that cannot afford the whole span covers the part it can integrate stably,
reports the shortfall in `dt_used`, and has its clock advanced by what was
actually integrated rather than by what was asked for. `advance_to` re-asks the
unreachability question on the way out using the step the loop *demonstrated*
rather than the one `node_dt` predicted, and thermalises or counts
`Stats::unreachable` accordingly. `key_of` clamps to `KEY_LIMIT` and sends
non-finite coordinates off-grid instead of to cell zero.

`no_node_flings_its_bodies_out_of_itself` is the spread check, as a test, and
it is `#[ignore]`d because **it still fails**. That is the useful part. The
clamp fix cut the overshoot from 4.5x10^8 radii to 4.0x10^7 and stopped there,
which is how we know the clamp was never the cause either.

**Re-measured after the binding split**, and the first thing it catches has
moved — which is worth recording, because the two failures are different faults
and the old one is gone:

```text
  before:  frame 0, node 18 (Atomic),    r=5.917e-12 m, body 47 at 2.371e-4 m — 4.0e7 radii
  after:   frame 1, node 14 (Molecular), r=6.720e-9  m, body  3 at 6.300e2  m — 9.4e10 radii
```

Frame *zero* is the tell: the old failure was the sampler inflating the node at
birth, before anything was integrated, and it is fixed. What is left fires on
frame one, during a solve, in a node holding 8,000 molecules at a spacing equal
to their own Lennard-Jones sigma — a liquid-density parcel with no container,
handed a span a galactic pace supplies. It expands, and nothing is wrong with
the expansion except its rate. Still ignored, still the reminder.

**The cause. Measured this time, and it is not what the two paragraphs above
first replaced the earlier guess with either — the packing is not unphysical,
the solver assignment is wrong.**

```text
  node  14  tier=Molecular  MolecularDynamics  n=8000  Molecule  r=6.720e-9   sep=1.540e-11
  node  15  tier=Atomic     MolecularDynamics  n=   8  Atom      r=1.680e-10  sep=1.331e-10
  node  16  tier=Atomic     MolecularDynamics  n=  56  Nucleon   r=5.292e-11  sep=5.704e-12
  node  17  tier=Atomic     MolecularDynamics  n=  56  Nucleon   r=2.264e-11  sep=1.342e-12
  ...
  node  22  tier=Atomic     MolecularDynamics  n=  56  Nucleon   r=2.762e-14  sep=2.090e-15
```

Node 15 is a real atom and integrates perfectly. Nodes 16 onward hold **56
nucleons in a Woods-Saxon profile** — which is exactly right, a nucleus *is*
fifty-six nucleons — and are labelled `Tier::Atomic`, so `solvers::for_tier`
hands them **`MolecularDynamics`**. Lennard-Jones, with a sigma of about
1e-10 m, applied to nucleons a few femtometres apart. The contents are correct
physics; the solver is for a scale five orders of magnitude coarser.

**Why the tier is wrong.** `promote` derives the child's tier from
`Tier::containing(body.radius)`, and the body being promoted out of node 15 is
an *atom*, so its radius sits squarely in the Atomic band. The spec it is
filled with, `default_spec(Atomic.finer())`, is the Nuclear one. Tier from the
radius, contents from the spec, and nothing reconciles them.

**Why it then runs away.** `drill` stops when `tier >= to_tier`. The tier never
reaches Nuclear, so it promotes the most massive *nucleon* into a smaller node,
fills it with 56 more nucleons, and repeats — twenty-four levels of nucleons
splitting into nucleons, each one narrower than the last, until the radius
finally falls below `Tier::Atomic.floor()`.

**The guard for this is already written, and is half of one.** `tree.rs`
lines 404-427 replace a spec meant for a *coarser* scale with the tier's own
policy, and deliberately keep one meant for a finer scale: "asking to split an
atom into nucleons is a deliberate step down and not a mistake, and overriding
it would leave the ladder unable to reach its own bottom." That is right. But
having kept the finer spec, it leaves the node's tier — and therefore its
solver, its timestep and its cadence — describing the coarser scale it took
the radius from.

**The fix**, then, is narrower than the two options this entry previously
offered, both of which were aimed at the wrong target: when the spec is
deliberately finer than the radius-derived tier, the **tier should follow the
spec**, not the radius. `tier = tier_of(spec.kind).max(...)` rather than
`Tier::containing(body.radius).max(parent_tier)` in that branch. Node 16 then
becomes Nuclear, takes the `Statistical` solver, no force field touches a
nucleon, and `drill` terminates where it was asked to.

Two things to check before believing that is all of it: a node of radius
5.3e-11 m holding a nucleus is carrying its *atom's* radius, which is
defensible for a node that represents the atom but means the Woods-Saxon
profile is being scaled to the atom rather than to the nucleus (hence
`sep=5.7e-12` where a real nucleus is ~1e-15 m across); and the ladder's shape
changes, since the path shortens from 24 to about 17.

**Trigger:** now, in the sense that the ignored test is the reminder. Nothing in
the play space reaches Atomic tier yet, so it blocks nothing — but the debug
suite cannot be a gate until it is green, and every MD result below 10^-10 m is
meaningless until then. The rest of what stopped it being a gate is gone: the
suite runs in 123 s and is otherwise green, so this ignored test is the only
thing left between here and using it as one.

The `MdParams::default().cutoff` inconsistency noted above is untouched and
still latent.

---

## Growth accumulates internal energy that nothing sheds

**Noticed:** measuring why the loose contents of a grown tree leave it at
14 km/s, while implementing `PLAY.md` §3.3.
**Where:** `engine.rs` — the growth transaction, `node.matter.internal_energy
+= txn.heat_released`.

Measured, on a 900 kg tree planted at 291 K and grown a year at a time:

```text
    after plant     U = 4.1782e8 J      thermal at 291 K = 4.1782e8   (exact)
    after  1 year   U = 6.0499e8 J      mass 900.0 kg, built   2.2 kg
    after 10 years  U = 7.4845e9 J      mass 900.0 kg, built  47.7 kg
    after 40 years  U = 9.6555e10 J     mass 900.0 kg, built 528.8 kg
```

231 times its own thermal energy, at unchanged mass. Read back as a
temperature that is about 67,000 K — while `matter.temperature` still says 291,
because growth adds to `internal_energy` and never touches the field the rest of
the engine reads.

**The intent is already written down and is not what happens.** The line
carries the comment "Only the thermalised share stays. What was re-radiated has
left the node, and adding it here would cook a forest in a season." It is
cooking the forest; it is taking forty years rather than one.

**Why it had not been seen.** `evolve_matter` is what sheds heat and it works
from `matter.temperature`, which growth leaves alone — so the one mechanism that
would have caught it is looking at the wrong field. And nothing asked the node
for its energy: a grown tree is normally *damaged* or *shaken*, paths that read
the topology and the morphology, never the thermal state.

**What it breaks now.** `sample_structured` converts `internal_energy` into the
parts' kinetic energy, faithfully, so a forty-year-old tree materialises with
its branches and its litter moving at 14 km/s. §3.3's dispatch keeps the members
out of the tier solver and the Courant substepping stops the parcels being
accelerated further, but neither addresses the speed they are *born* at. A tree
should materialise still.

**Candidate causes, not yet separated.** Either the growth transaction is
thermalising a share it should be re-radiating, or the share is right and
nothing is radiating it afterwards because `evolve_matter` cannot see it. The
second is the more likely and the cheaper test: make growth move `temperature`
with `internal_energy` and see whether `evolve_matter` then carries it away.
That is one probe and it has not been run.

**Trigger:** before anything materialises a grown structure and expects it to be
still — which is the first time an observer looks at a tree. It is also a
prerequisite for the ball-in-box test moving to Continuum, since that test is
about things that are and are not moving.

---

## Derived gravity is in the parent's axes, because nothing composes orientation

**Noticed:** building `Tree::gravity_at` for `PLAY.md` Phase 1's "`G_EARTH`
deleted in favour of derived g".
**Where:** `tree.rs` — `gravity_at`, and `offset_from` underneath it.

`gravity_at` walks the ancestor chain and adds what each one pulls with, which
gives 9.82 m/s^2 on an Earth-mass, Earth-radius body and 10^-13 in a galaxy. The
*direction* is right in the parent's axes and is not necessarily right in the
node's own: `offset_from` composes offsets and not rotations, and `Motion`
carries an `orientation` that nothing consults.

So a node on the `+x` side of a planet is told gravity points along `-x`. A
structure is generated with `+z` up — `ground_of` returns a `z` — so such a node
would carry its weight sideways through its own geometry.

**Why it does not bite yet.** Nothing orients a node against the body it is on.
There is no terrain, no surface, and no reason for a node's frame to be anything
but its parent's. Both tests that needed a planet were placed on its `+z` axis,
where the question does not arise, and they say so.

**Why it will.** A cubed-sphere terrain patch is oriented by construction: that
is what a patch *is*, a square of surface with a local up. `PLAY.md` Phase 2
builds them, and the first structure emplaced on one at any latitude but the
pole reads its own weight in the wrong direction. Actors standing on a planet
are the same problem one level down.

**What it needs.** `offset_from` to compose orientation as well as offset, or a
`gravity_at` that rotates the accumulated field into the node's frame on the way
out. The second is smaller and is probably wrong: the same gap affects any
vector carried across frames — a velocity, a wind direction, an impulse — so the
fix belongs where frames are composed rather than at one consumer. `coords.rs`
already carries the error-bound argument for offsets and would be the place.

**Trigger:** the first oriented node, which is Phase 2's first terrain patch.
Before that there is nothing to get wrong.

---

## ~~The ball-in-box test runs at the wrong tier~~ — done, and it found three things

**Was:** both tests built metre-scale solids at `Tier::Galactic`, the
collisionless-gravity regime, because `for_tier(Continuum)` is SPH and SPH reads
a solid through a gas equation of state. Measured then: this box's `pressure()`
was 1.0x10^8 Pa against green wood's 4.5x10^7 tensile strength, so it burst from
its own equation of state before anything touched it.

**Confirmed still true before the rebuild**, since the entry was a year of
commits old: 1.04x10^8 Pa, and a box at `Continuum` with the old arrangement
reached 1159..2857 m from a 3 m box in 100 frames at 2.86x10^3 m/s.

**What the rebuild actually was.** Not a tier swap. The box is now a *node*,
because a node is already what this engine means by a rigid body — one velocity,
one spin, and `apply_contact` has always written to both. What a node lacked was
a **shape**, so `Node::collision_shape` gives it one: a capsule per structural
member, partitioned by §3.3's own `structural_mask` rather than by a second
opinion about what a node is. The root stays `Galactic` because the root is deep
space; the *solids* are what moved to `Continuum`, and `tier_for` puts them there
from their own radii without being asked.

```text
                       before (96 panels, Galactic)      now
    momentum drift           8e-7                      1.0055e-14
    angular momentum         5e-7                      7.3675e-10
    box response        one 500 kg panel at 0.128    48 t at 1.1468e-3 m/s
    box afterwards      grew by ~6 m                 still a box
    parts                   96 spheres                48 capsules
```

The box velocity is the rigidity check and it is arithmetic: the ball carries
53.5724 kg m/s, so a rigid 48-tonne box takes 1.1161e-3 m/s. Measured 1.1468e-3,
2.8% over because the ball ends travelling slowly backwards.

**Three things it found on the way, which is why the entry is worth keeping.**

1. **A node made entirely of ordered matter never moved.** The `count == 0`
   early return in `advance_node` sat *before* the node's own clock, motion,
   contact and exchange, so a solid object — every body a member, nothing loose
   — was frozen in space. Measured: a 48-tonne box given 1 m/s for one second of
   frames moved 0.000000 m with its clock reading 0.0 s after a hundred steps;
   the same box with a single loose body added moved 1.000000 m. **Fixed** —
   skipping the solver is not the same statement as skipping the node.

2. **A `Continuum` node whose bodies are all stand-ins detonates.** See the new
   entry below. Not fixed.

3. **`E = rho c^2` does not give a surface.** See the surface entry. Not fixed.

---

## ~~A struck node banks angular momentum and never turns~~ — done

**Noticed:** answering what a node's velocity and spin actually are, while
scoping `PLAY.md` §2A.
**Where:** `coords::Motion::spin_rate` against `state::Matter::spin`.

A node carries two spin quantities and nothing keeps them in step:

```text
    matter.spin        angular momentum      kg m^2/s
    motion.spin_rate   angular velocity      rad/s
```

`spin_rate` is derived from `matter.spin` at **node construction and nowhere
else** — two sites, both creation paths. A collision adds to `matter.spin`
through `apply_contact` and `spin_rate` never hears about it, so the node's
`orientation` never moves. The ball-in-box test measures a box accumulating
`matter.spin` of 7.4865 over five off-centre strikes while remaining, as far as
`motion` is concerned, perfectly still.

This is the concrete mechanism behind "nothing composes orientation anywhere in
the tree", which until now was recorded only as a consequence for derived
gravity's axes.

**Second defect, same area.** `Node::collision_shape` translates a child's
pieces to where the child is and never rotates them by `motion.orientation`. So
even with the above fixed, a spinning box's walls would stay where they were
built. Both are `PLAY.md` Phase 2 item 1, and neither is more than a few lines.

**Trigger:** immediately, as the first item of Phase 2 — the rigid-body claim
that `a_ball_loose_in_a_box...` rests on is not true until both are done.

**Done**, as Phase 2's first item, and both halves are held by a test that was
checked against the defect it claims to catch.

`Node::sync_spin_rate` re-derives the angular velocity from the angular
momentum, and is called from the three places where the momentum, the mass or
the radius can have moved — after a solve, after a contact, and after a
coarsen folds a node's own detail back into it. Deliberately *not* from
coasting, which asserts that nothing changed. `Hull::placed` turns a hull by an
orientation and then moves it, and `contact_within` uses it for a promoted
child, which is the same composition `Motion::body_to_parent` uses.

Measured on `a_plank_struck_off_centre_turns_and_its_surface_turns_with_it`, a
1-tonne plank struck near one end by a 50 kg ball:

```text
                        before              after
  matter.spin           835.2953 kg m^2/s   835.2953 kg m^2/s
  motion.spin_rate      0 rad/s             0.23203 rad/s about z
  orientation.angle()   0 rad               0.44317 rad
  far end of the capsule moved                1.3187 m
```

The angular momentum was always right; nothing read it. The chord a point at
3 m swings through 0.44317 rad is `2 * 3 * sin(0.44317/2) = 1.31866`, which is
the measured figure to seven digits, so the shape really is turning with the
node and not merely translating somewhere new.

`a_turned_plank_is_struck_where_its_wall_now_is` is the second half on its own:
a plank turned a quarter turn before anything happens, struck 2.4 m along where
its unturned self does not reach. Against a `contact_within` that translates
without rotating, the ball passes straight through and the run records zero
contacts.

---

## A Continuum node whose bodies are all stand-ins detonates — visible now, fixed in Water

**Noticed:** building the ball-in-box scene, where a root holds two promoted
children and nothing else.
**Where:** `engine.rs::advance_node`, `state::Matter::pressure`.

A node whose every body is the stand-in for a promoted child has no ordered
contents of its own, so `structural_mask` returns `None` and the whole lot goes
to the tier solver. At `Continuum` that is SPH, and SPH prices the node's
`Matter` through the gas-plus-radiation equation of state — which, for a node
that is mostly vacuum with two solid objects in it, is nonsense. The force lands
on the stand-ins and D4 hands it straight to the children.

Measured, on a 12 m root at `Continuum` holding a 48-tonne box and a 10 kg ball,
with **zero collisions** in the run:

```text
    ball speed, m/s      5.36 -> 7.96 -> 20.9 -> 63.3 -> 194 -> 566 -> ...
                         (roughly x3 per frame, geometric)
    with the root not advanced:  5.3572 every frame, exactly
```

The control is the point: the blow-up is the parent's own solver, not contact.

**Why the obvious fix is wrong.** Excluding stand-ins from the solver would
break D4 outright — `stand_in_velocities` before and `apply_body_forces` after
are *how* a promoted child feels its parent's forces, so a stand-in that is not
integrated is a child that feels nothing. The stand-in has to be solved; what is
wrong is the equation of state it is solved under.

**So this is the liquid-EOS gap seen from the parent's side**, and `PLAY.md`
Phase 3's second piece names it: "a liquid equation of state, so `pressure()`
stops returning 4x10^8 Pa for a bucket of water". There is no EOS anywhere for a
node that is mostly empty space.

**Worked around, not fixed**, in `a_ball_loose_in_a_box...`: the root is
`Galactic`, which is what deep space is, and the solids are `Continuum` on their
own radii. That is honest for a scene in space and is no help at all to a room
with furniture in it.

**Trigger:** the first `Continuum` node that contains promoted children and is
not empty space — a room, a vehicle interior, a crate of objects. Phase 3 at the
latest.

**Made visible**, as Phase 2's ninth item, and deliberately not fixed:
correcting it needs a liquid and a solid equation of state, which is `PLAY.md`
Phase 5's second piece and is named there in as many words.

`Matter::gas_law_applies` is the reading, and D17 is what made it possible — a
node carrying a mixture knows its own phase, so "is this matter a gas" stopped
being a guess. `Stats::eos_outside_validity` counts the node-steps where the
fluid solver priced matter the gas law does not describe, and
`eos_outside_validity_at` says which node, on §3.7's precedent.

Two ways to be outside it, and both are reported:

- **The matter is not a gas.** A bucket of water prices at 4x10^8 Pa.
- **The node is mostly vacuum with solids in it**, which is this entry's own
  case: every body a stand-in for a promoted child, so there are no contents of
  its own for a fluid solver to be about.

**Matter nobody has described is deliberately not reported.** "No information"
is not the same answer as "measured and wrong", and almost every node in a
galaxy genuinely is a gas — `a_gas_law_asked_about_a_solid_says_so` in
`tests/dispatch.rs` holds all three cases, including that one.

---

## ~~The narrow phase runs against every piece a recipe emits~~ — done

**Noticed:** taking Phase 2's entry baseline, and again wiring D18's stored
surface into contact.
**Where:** `shape::closest_of`, `engine::contact_within`.

`closest_of` is N x M over the pieces two sides present, and the cost is flat
per piece, so it is linear in the piece count. Measured, in `PERFORMANCE.md`:

```text
  pieces on one side   closest_of   per piece   whole contact
                   1      0.31 us     0.31 us         0.40 us
                   6      1.94 us     0.32 us         2.04 us
                  64     20.1  us     0.31 us        20.2  us
                 512    161    us     0.31 us       162    us
```

A wooden box of six walls is 2 us a contact. **A generated tree is 3,400
members**, which at the same rate is 1.1 ms for one pair — a fiftieth of a
frame, for one tree touching one thing.

D18 is explicit that the generator states the pieces and nothing infers a
decomposition, so the count is the recipe's to choose. What is missing is the
other half: `PLAY.md` §7's item 8, "collision runs against the surface **at the
level of detail the distance deserves**". A tree a hundred metres away does not
need its twigs, and the broad phase already knows the distance.

**Not done, and deliberately separated from the bake.** The bake is what D18
asks for and it is correct; picking a level of detail is a second decision with
its own measurement, and doing both at once would have made it impossible to
tell which one moved the numbers.

**Measured cost of the bake itself**, which is the part that is done: four
surfaces over a twenty-five-frame run of `frames_stay_within_budget`, because
`epoch` is what invalidates one and an undisturbed node never moves it.

**Trigger:** the first scene where a tree or a settlement is close enough to be
collided with. `PLAY.md` §7 item 8.

**Done**, and the saving is **exact rather than approximate**, which is what
makes it a level of detail rather than a fudge: a hull's `bound` contains it, so
a piece further from the other side's bounding sphere than the two bounds
together cannot be the nearest pair, and skipping it changes no answer.

Re-measured on a slower moment of the same container — the sphere-against-sphere
row moved from 0.106 µs to 0.174, so the comparison is between columns:

```text
  pieces on one side   touching   out of reach   reaching a few
                   1    0.56 us        0.03 us          0.56 us
                   6    1.20 us        0.09 us          1.89 us
                  64    1.63 us        0.40 us          2.21 us
                 512    5.15 us        2.99 us          6.32 us
```

**512 pieces went from 161 µs to 5.15 µs**, about fifty times once the machine
is accounted for, and per piece is no longer flat — it falls from 0.56 µs at one
piece to 0.010 at five hundred, because the count actually *tested* stops
growing.

`pieces_out_of_reach_are_not_tested_and_the_gap_is_a_lower_bound` in
`tests/shape.rs` pins the one semantic change: when nothing can touch, the query
answers from the bounds and returns a gap that is a **lower bound** on the true
separation rather than the separation itself. Measured, a capsule against a ball
40 m off its end: 34.4 m from the bounds against an exact 39.4. Every caller in
the engine uses it to decide whether a pair is in contact, and a lower bound is
sound for that; one that wanted the true distance would have to ask for the
descent.

**What is left is O(n) rather than O(1):** the bounding sphere over a side is
recomputed per query, which is the 2.99 µs in the middle column, and it could be
cached on the surface alongside the pieces. Nothing measured needs it yet.

---

## ~~Collision geometry is a sphere~~ — mostly done; the residual is flatness

**Was:** the general contact path knew a thing by a centre and a radius, so a
26 m settlement beam presented a 13.5 cm bead at its midpoint and a ball thrown
at a timber frame passed through it unless it happened to hit one.

**What landed.** `src/shape.rs`: a `Hull` over a *set of spheres*, chosen
because a sphere is what the engine already holds everywhere — a sampled body, a
member's two ends, a fluid parcel. So the hull of two spheres **is** a capsule
rather than a faceted stand-in for one, and nothing is stored that was not
already there. The narrow phase is GJK on the hulls' cores with the radii
restored at the witness points, and there is deliberately no EPA: `contact`
resolves from a normal, a contact point, two masses and a closing velocity, and
never reads a penetration depth.

Of the three options this entry listed, it is the first — **a shape per
occupant** — arrived at from the third. Both sides of a contact now carry
geometry:

- A **promoted child** presents `Node::collision_shape`: a capsule per
  structural member, or the sphere of its own radius when its contents are
  disordered. Partitioned by §3.3's `structural_mask`, so it is the same
  ordered-versus-disordered question the solver dispatch already asks.
- A **body of the node** presents its member capsule too, which is the side a
  limb landing on a tree strikes.

Measured on materialised programs: every structural member now presents a
capsule — Tower 248/248, Tree 3400/3400, Wall 704/704, Settlement 256/256.

**What is left, and it is real.**

- **A wall is scalloped, not flat.** A row of capsules is not a plane. Measured
  on the ball-in-box walls, on axis: the inner face sits at 2.500 m where a real
  wall would be at 3.000 m, against 1.875 m for the ninety-six-sphere version —
  the error halves, 1.125 m to 0.500 m — and the scallop between neighbouring
  capsules falls from 0.287 m to 0.169 m. Flat needs either a hull *grouping*
  rule (which members share one convex piece) or a genuinely planar primitive,
  and **neither is decided**. Grouping is the harder half: a support subtree is
  not convex for a tree, and `site` names a failure rather than a convex piece.
- **Terrain is untouched.** The signed-distance option this entry lists is still
  what a cubed-sphere field wants, and Phase 2 will meet it first.
- **One contact per pair per frame.** `closest_of` takes the minimum-gap pair of
  pieces, so a ball wedged into a corner resolves one wall this frame and the
  other next. Resolving both at once over-corrects, which is the standard reason
  a solver iterates; that is a scheduling question of the kind §3.3 left open.
- **Cost is still unmeasured.** `PERFORMANCE.md` has no row for contact. The
  narrow phase went from subtracting two radii to a GJK over small point sets,
  and nobody has priced it.

**Trigger for the remainder:** flatness, when something has to rest or stand on
a built surface rather than bounce off it — Phase 4, or Phase 2 for terrain.

---

## A struck structure is not rigid — resolved for a node, open for a member

**Noticed:** writing `a_ball_loose_in_a_box_conserves_momentum_and_angular_momentum`.
**Where:** `engine.rs::apply_contact`.

**The half that is closed.** A ball no longer rebounds from one 500 kg panel of
a two-tonne box. The answer was not a rigid-body group bolted onto contact: a
`Node` is *already* the rigid body — one velocity, one spin, and `apply_contact`
has always written to both — so what was missing was a shape for it to present,
which `Node::collision_shape` now supplies. Measured: a 48-tonne box takes a
10 kg ball's 53.5724 kg m/s and moves at 1.1468e-3 m/s against an arithmetic
1.1161e-3, with momentum conserved to 1.0e-14 and angular momentum to 7.4e-10.

This also means §3.3 needed no revisiting. "The members are left where they are"
is right, because a struck box's motion belongs to the *node* and not to its
panels — which is exactly why the frozen-member problem never appears once the
box is a node rather than a heap of bodies in the root.

**The half that is open.** Where the struck thing is a **member of the node
doing the striking** — a limb landing on the tree it fell off — the impulse
still lands on `bodies[k].vel` and that member alone. The member now presents
its capsule, so the *geometry* is right and the *response* is not.

`solvers::structure` still has the machinery and contact still does not call it:
`Mechanism::PointImpulse` and `World::damage` resolve an impulse against the
whole frame, and `drop_fragments` uses exactly that. The reasons it was not
wired up are unchanged and are design rather than plumbing — `damage`
regenerates the structure and renumbers its members, so contact would have to
batch as `drop_fragments` does, and **which of three outcomes the contact path
should produce** (absorbed rigidly, a joint broken, a loose member knocked off)
is still not decided anywhere.

**Trigger:** the first time a *member* of a node has to respond as part of its
structure rather than on its own — a branch that should sway the tree instead of
flying off. Unavoidable by Phase 6.

---

## Exchange has a radiative coefficient and no conductive one

**Noticed:** building D3's transport half.
**Where:** `neighbourhood::radiative_conductance` is the only coefficient there is.

`exchange` takes a conductance and moves a conserved quantity across a boundary
at that rate. `PLAY.md` D3 says "conduction, diffusion and radiative exchange
are three calls to it with different coefficients", and one of the three exists.

**Radiation was taken first because it is the one that is derivable with what is
already here.** `sigma A (T_a^4 - T_b^4)` factors exactly into a conductance —
`(T_a + T_b)(T_a^2 + T_b^2)` — so it needs no linearisation and no new
property, and the exchange area comes from the two radii and their separation.

**Conduction has no coefficient anywhere in the codebase.** Checked:

- `topology::Material` carries `specific_heat`, `thermal_onset`, `thermal_gone`,
  `destruction_enthalpy` and an **electrical** `resistivity` — "ohm-metres,
  sets how a conducted discharge distributes its energy between members". No
  thermal conductivity.
- `chem::analyse::Properties` carries `unit_mass`, `molar_mass`,
  `cohesive_energy`, `ionicity`, `polarity`, `dipole`, `hydrogen_bonds`,
  `lattice_binding_ev`, `density`, `melting_point`, `boiling_point`. No
  transport quantity of any kind.
- `solvers::structure`'s `Mechanism::ThermalField` takes a **convective**
  coefficient `h` as authored data — a number the caller supplies for a fire or
  an immersion, not a property of the material.

**What the axioms will and will not allow.** A table of conductivities per
material is what axiom one forbids. The derivable candidate is the
Einstein-Cahill-Pohl minimum conductivity, `k_min ~ (pi/6)^(1/3) k_B n^(2/3)
v_s`, which needs a number density and a sound speed — and `Matter` already has
both, `number_density()` and `sound_speed()`, so it needs no new stored
property at all. It is a *lower bound* on a real crystal's conductivity, right
to within a factor of a few for a disordered solid or a liquid and an order or
more low for a good crystal or a metal, where electrons carry most of the heat.

Whether that bound is the right answer, whether an electronic term should be
derived alongside it from `resistivity` through Wiedemann-Franz, or whether
something else entirely, is a `PHYSICS.md`-weight decision and is not made.

**What this blocks.** Two things touching do not conduct, only radiate, which
is wrong by orders of magnitude for anything in contact: a hand on cold metal,
a pan on a hob, heat spreading through a wall. Radiation is the correct and
dominant mechanism across a gap, so nothing is wrong for things that are merely
near each other. Diffusion is in the same position and is wanted by Phase 3's
water rather than by anything now.

**D17 unblocks the choice, and `PLAY.md` §7 schedules it in Water.** The
Einstein-Cahill-Pohl bound needs a number density and a sound speed, both of
which `Matter` has; the objection was that picking it is `PHYSICS.md`-weight.
With a mixture on every node there is also a real *phase*, so the coefficient can
differ for solid, liquid and gas rather than being one bound stretched across all
three — which is most of what made the single number hard to defend.

**Trigger:** the first scenario where two touching things have to reach the same
temperature — which is contact heating, cooking, or anything a hand rests on.
**Water** at the latest, since a free surface exchanging with what it sits on is
the same function.

---

## ~~Only a built thing has a surface, so only a built thing collides~~ — done

**Noticed:** building D3's contact half.
**Re-opened:** the owner, asking why a wooden ball is not a "built" thing, and
observing that "grown" and "built" are themselves two special cases where the
axioms allow only one kind of thing.
**Where:** `engine.rs::World::surface_of`.

Contact needs three numbers — density, Young's modulus and a strength — and they
come from `topology::Material`. `surface_of` reads `topology.material` for a
materialised structure or `morphology.material()` for one not yet materialised,
and returns `None` for everything else, so `neighbourhood::contact` refuses the
pair.

**That is a provenance test standing in for a state measurement**, which is how
this entry should have read and did not. The engine has already made the same
decision correctly once: `PLAY.md` §3.3's `Node::structural_mask` partitions a
node's contents by the *measured* joint radii it carries, not by whether it has
a `Program`. A building, a wolf and a boulder are all `Continuum` and the
dispatcher decides from state. One layer over, contact decides from a birth
certificate. The question a surface answers is not "was this built" but "does
this hold together", and that is measurable.

**What the ball-in-box test does about it is the tell.** `tests/adjacency.rs`
hand-authors `Topology { material: GREEN_WOOD, ..Default::default() }` onto both
the ball and the box, with no joints, for no reason except to make `surface_of`
return something. A hand-set value on a path that exists to measure is what
`Interaction::Author` is audited for, and here it is load-bearing for the test to
run at all.

**The derivation this entry used to record does not work, measured.** It said
Young's modulus falls out exactly, because "the longitudinal wave speed of a
solid *is* `sqrt(E/rho)`" and `Matter::sound_speed()` already exists. Both halves
fail:

- `Matter::sound_speed()` is **the gas formula** — `sqrt(gamma n k T / rho)`,
  capped by `velocity_dispersion` — and not an elastic wave speed at all. It is
  the same equation-of-state mismatch that makes SPH burst a solid.
- `Matter::density()` is the **bulk** density of a node, and a surface is a
  property of the material at the point of contact. A hollow box is mostly air.

Measured on the ball-in-box assembly, against what green wood actually is
(600 kg/m^3, 1.0x10^10 Pa):

```text
                    density()              E = rho c^2
  box node      5.31e1   11.3x low      1.73e8   57.7x low
  ball node     3.73e1   16.1x low      1.22e8   82.1x low
```

So **none of the three numbers is currently derivable from a node's own
`Matter`**, which is a stronger statement than this entry made before. Density
and stiffness are not missing laws; they are being asked of the wrong object.
They are properties of *what the matter is made of*, and the place that knows is
`chem::analyse::Properties` — `density`, `cohesive_energy`, `lattice_binding_ev`
— which is reachable for a node with a registered `Mixture` and not from a bare
`Composition`.

**Strength is the one with no law at all.** Theoretical strength is about `E/10`
for a perfect crystal and real materials are one to three orders below that
because of dislocations, which the engine does not represent. The derivable
candidate is a cohesive energy density — `cohesive_energy` over a formula unit's
volume, in Pa — an honest *upper* bound of the same shape as the
Einstein-Cahill-Pohl minimum conductivity this file proposes for the thermal gap.
**Deliberately not chosen**, and deferred to D11 with the owner asked again
there.

**The owner's wider point, recorded rather than resolved.** "There should only be
'things' that everything is" — and grown versus built is the same provenance
split one level up. `PLAY.md` D11 already goes most of the way: growth and
construction "are the same operation — material is deposited where a field says
to deposit it — and they differ only in where the field comes from", and five of
`Program`'s seven columns move onto the material and the measured environment.
But D11 explicitly stops short of unifying geometry: "Branching, coursed masonry
and a subdivided street grid are genuinely different space-filling rules, and
asserting they collapse into one would be the third over-claim in this document."
Whether the *habits* are themselves a species table in disguise is the owner's
question and is not settled here.

**Trigger:** when anything that is not a built structure has to be collided with
— the first loose object an actor can pick up, kick or trip over. Not before D11,
unless something needs it sooner.

**Done**, as Phase 2's fourth item, and the entry's own diagnosis was right:
none of the three numbers is derivable from a node's own `Matter`, because they
are properties of *what the matter is made of*. D17 put a `Mixture` there, and
`Material::measured` reads the solid pools of it.

`World::surface_of` now measures first and inherits second: a node that carries
a mixture answers from it whatever made it, and `Morphology::material` is the
fallback for a node whose chemistry nobody has stated. A node with no solid pool
gets `None`, which is D13's own line — a liquid's surface belongs to its
container and a gas has none.

The hand-authored `Topology { material: GREEN_WOOD }` this entry complained
about is gone from `a_ball_loose_in_a_box...`: the box and the ball are given a
cellulose mixture and the material is measured.

**And the strength gap this entry left open is closed**, by `PLAY.md` D14's
Griffith law rather than by the cohesive-energy-density upper bound it was
deferring to. See the new entry below for what that reproduces and what it does
not.

---

## Griffith on a grain is not Griffith on the worst flaw — bedrock is 77x low

**Noticed:** building `PLAY.md` D14, measuring the derived strengths against the
`rupture` column they replace.
**Where:** `material.rs` — `Material::strength`, `grain_scale`.

Seven of the eight presets land within a factor of 2.3 of the retired table, and
the ordering is the table's own except for two positions. The measurements, none
of which the derivation was shown:

```text
  material           a            derived     retired    ratio
  green wood       3.0e-5 m       5.16e7      4.5e7      1.15
  dry timber       1.2e-5         6.52e7      7.0e7      0.93
  aragonite        3.0e-3         1.35e7      1.2e7      1.13
  masonry          7.0e-2         2.07e6      2.0e6      1.04
  ice              3.8e-2         3.85e6      1.7e6      2.27
  steel            8.2e-3         1.94e8      4.0e8      0.49
  reinforced       1.3e-3         7.90e7      1.8e8      0.44
  bedrock          2.1e-1         1.66e6      1.3e8      0.013
```

**Bedrock is the one real miss**, and the reason is measurable rather than
mysterious: a granite at 130 MPa with the stiffness and surface energy this
derives implies a **34 µm** crack, against the **21 cm** grain the nucleation
solve gives it. The cracks that matter in rock are *inside* the grains, not
around them — D14 says as much in a line ("a grain is also not the same as the
worst flaw") and does not say what the intragranular scale is.

**The frozen branch's grain sizes are plausible and the nucleation solve is not
the problem.** Measured: steel 8.2 mm (a slowly cooled ingot), ice 3.8 cm (lake
ice is centimetre-grained), bedrock 21 cm (coarse, and plutonic grains do reach
centimetres). What is missing is the step from a grain to the worst crack in it.

**Two further approximations are recorded here rather than hidden.** A metal's
cohesive energy comes out about 2.3x low, because iron really has eight nearest
neighbours and the valence model allows three — so steel's stiffness derives at
7.5e10 against 2.0e11, and its yield at 1.94e8 against 4.0e8, both low by the
same root cause. And ice's density derives at 1653 kg/m^3 against 917, because
`analyse`'s van der Waals packing estimate is a correlation rather than a
structure.

**Trigger:** when something brittle and crystalline has to break *correctly*
rather than merely break — a rock face that spalls, a stone wall that a siege
engine has to beat down. Until then the direction of the error is the safe one:
bedrock is weaker than it should be, so anything standing on it is
conservatively judged.

---

## ~~The sampler inflates anything bound by chemistry by 4.3x10^5~~ — done, and it found two more

**Noticed:** building D3's exchange pass. `Neighbourhood::pairs` returned zero
for a refined granite block, which should be the easiest case in the engine.
**Where:** `sampler.rs` lines 246-259, the relaxation loop.

A second, unrelated cause of the same symptom as *"A node whose bodies are
10^20 radii outside it"* above — that one is an MD detonation during
integration, this one happens at sample time, before anything is integrated.

**Measured**, `sampler::sample` on each scenario at `count = 512`:

```text
  scenario          tier        R          U_int      E_bind      relax  max|pos|/R
  Spiral galaxy     Galactic    4.629e20   +1.086e48  -2.896e47   0      3.109e0
  Molecular cloud   Stellar     6.171e17   +2.138e42  -2.566e42   0      3.569e0
  The Sun           Planetary   6.957e8    +1.138e41  -2.276e41   0      3.512e0
  Rocky planet      Planetary   6.371e6    +2.242e31  -2.242e32   0      3.512e0
  Granite block     Continuum   8.660e-1   +4.154e8   -6.392e10   33     4.396e5
  Water vapour      Molecular   3.000e-9   +1.279e-17 -7.891e-16  33     4.592e5
  Carbon atom       Atomic      7.000e-11  +8.251e-17 -1.650e-16  33     5.158e5
  Iron nucleus      Nuclear     4.591e-15  +1.776e-10 -7.887e-11  0      1.224e0
```

`1.5^32 = 4.3x10^5`, which is the whole of the discrepancy.

**The cause.** The loop exists for a real case: a configuration too tightly
bound to hold the energy it claims must be bigger, so it is scaled up by 1.5
until `internal_energy + binding_energy - phi` turns positive, up to 32 times.
That is right for anything **self-gravitating**, because spreading it out is
exactly what releases the binding — `phi` is the gravitational potential, and it
is the only term the scaling moves.

It is wrong for anything bound by **chemistry**. A granite block's
`binding_energy` is a silicate cohesive energy, `-5 eV` per atom; its `phi` is
about `1e-4 J` and utterly negligible at that size. Scaling the geometry cannot
make the budget positive because it does not touch the term that is negative.
So the loop runs all 32 iterations, fails, falls through to the
`random_ke_target = internal_energy.abs()` branch — **and leaves the 4.3x10^5
inflation in place.** Nothing undoes it.

The nucleus is the control: its binding is nuclear, equally non-gravitational,
and it never triggers because its Fermi energy exceeds it.

**Why nothing caught it.** The conserved-set tests check that the books close,
and they close either way — the energy budget absorbs whatever `phi` comes out
to, which is the sampler's design and is correct. Nothing had asked a
*geometric* question of a Continuum node until adjacency existed. `sample`
reports it honestly in `report.radius_overridden` and `report.relaxations`;
no caller reads either.

**Blast radius, today and tomorrow.** Today it is the three scenarios that set
a real chemical binding — granite, vapour, carbon. A node reached by descending
the galaxy ladder has `binding_energy == 0` from the second level down, so the
demo path never triggers it. Tomorrow it is everything: the trigger is
`internal_energy + binding_energy < 0`, which is every solid and every liquid
with a real cohesive energy, and §5A of `PLAY.md` is about making matter carry
exactly that. The Continuum tier is the whole play space.

Nothing geometric works on an affected node: adjacency, contact, the exchange
pass, and hydro's own neighbour finding all see contents scattered 4.3x10^5
radii from a node they are supposed to be inside.

**Not yet fixed, because the fix is a physics decision rather than a patch.**
The loop needs to know which part of `binding_energy` is gravitational — only
that part is released by expansion. The engine does not currently separate
them; `binding_energy` is one number and `scenario.rs` puts gravitational,
cohesive, covalent, electronic and nuclear binding into it. Options, none of
them chosen:

- **Split the field.** `binding_energy` becomes gravitational-only and a second
  carries the rest. Most honest, touches the conserved set and the wire format.
- **Compare against `phi` before relaxing.** If `phi` is negligible against the
  deficit, expansion cannot fix it, so do not try — recognise the budget as
  non-gravitational and take the fallback branch immediately, *without* the
  inflation. Smallest change, and it makes the loop's own precondition explicit
  rather than assumed.
- **Undo the scaling when the loop fails.** Narrowest of all, and it leaves the
  loop still wrong about what it is doing for 32 iterations.

**Now chosen, and promoted to a Phase 2 blocker.** `PLAY.md` D17 puts a
`Mixture` on every `Matter`, which is exactly what gives every solid a real
cohesive energy — so the "tomorrow it is everything" above is what D17 causes.
**Option one is taken: split the field.** It was the most honest and was held
back for touching the conserved set and the wire format, and D17 moves
`FORMAT_VERSION` regardless, so the two ride one migration instead of two.

**Trigger:** *now* — before D17, and therefore before anything else substantive
in Phase 2. Phase 2's own done-when is a rock bouncing off a boulder, and the
rock is this granite block. Previously: before anything at Continuum tier is
asked a question about where its contents are, which is D3's contact half, the
beach test, and every structure that has to sit on a surface.

**Done.** `binding_energy` is now `gravitational_binding` — both halves renamed,
so no call site could keep compiling while meaning something else — and
`cohesive_binding` carries the rest. Only the gravitational half is in the
relaxation loop's budget, which is the loop's own precondition made explicit
rather than assumed. `FORMAT_VERSION` moves to 9.

Measured, `sample` at `count = 512` on the whole shelf, before and after:

```text
  scenario          relax before / after    max|pos|/R before / after
  Spiral galaxy          0        0              3.354      3.354
  Molecular cloud        0        0              3.574      3.574
  The Sun                0        0              3.515      3.515
  Rocky planet           0        0              3.515      3.515
  Granite block         33        0          4.502e5        1.043
  Water vapour          33        0          4.387e5        1.017
  Carbon atom           33        0          4.738e5        1.098
  Iron nucleus           0        0              1.080      1.080
```

The four gravitationally bound scenarios are unchanged to every digit, which is
the control: the loop still does the job it was written for. `sample`'s cost is
unchanged — measured at 1.115 and 1.070 µs/body at 10^5 and 5x10^5 with the new
pass compiled out, against 1.116 and 1.041 with it in.

`every_scenario_samples_its_contents_inside_itself` in `tests/scenarios.rs` is
the guard, and it is the first *geometric* assertion any test has made about a
Continuum node — which is why nothing caught this: the conserved-set tests
close the books either way, because the energy budget absorbs whatever `phi`
comes out to and that is the sampler's design. `it_finds_the_two_defects_already_on_the_list`
in `tests/spread.rs` has its first half inverted rather than deleted.

### The two things it exposed, both of which it had been hiding

**1. An independent draw is an ideal-gas draw.** With the vapour node's
molecules correctly inside it, the closest pair sat at 2.55x10^-11 m against a
Lennard-Jones sigma of 3.12x10^-10 — eight per cent of contact, where
`(sigma/r)^12` is 8.6x10^12. The engine throttled to 2.83x10^-20 s and the pair
still left at 4.5x10^9 m/s. `sample` chooses every position without reference
to the others, and a real interacting system's pair correlation vanishes below
contact.

Fixed in the sampler rather than defended in the solver, which was the owner's
call: bodies whose kind means *one object* are pushed apart until their surfaces
no longer overlap. See `PHYSICS.md` §3 for the method and for the two bounds
that make it terminate. `SampleReport::worst_overlap` reports what could not be
resolved, on §3.7's precedent.

**Open, and recorded here rather than fixed:** an iron nucleus reports
`worst_overlap = 0.45` every time it is sampled. `child_radius` gives a nucleon
the 1.2 fm that is the radius *per nucleon* in `R = r0 A^(1/3)`, so fifty-six of
them fill their own nucleus exactly and no arrangement of spheres reaches a
packing fraction of one. A real nucleon's charge radius is 0.84 fm, which would
put the nucleus at 0.30 and inside the bound. Not changed here because it moves
the geometry of every nuclear sample and nothing currently depends on it.
**Trigger:** when anything reads `worst_overlap` to make a decision, or when
nuclear geometry is next worked on.

**2. Lennard-Jones was being applied to nucleons.** The carbon atom scenario is
a node the size of an atom holding its twelve nucleons — the right contents, and
what the scenario's own doc comment says it is for. `solvers::for_tier(Atomic)`
hands it molecular dynamics, and `dominant()` reads carbon off the composition:

```text
    body kind          Nucleon, radius 1.200e-15 m
    min separation     1.9016e-11 m
    sigma applied      3.431e-10 m
    (sigma/r)^12       1.190e15
```

The contents left at 150 c. Fixed by `md::has_electron_cloud`: the van der Waals
term applies between bodies that have electron clouds, and a nucleon, an
electron and a photon do not. Same shape as the bonded-pair exclusion beside it,
and for the same reason — a term that does not describe a pair is removed rather
than tuned.

**This is deliberately the *body* half and not the tier half.** The entry "A
node whose bodies are 10^20 radii outside it" above records the other one: a
node taking its solver from a radius while its contents came from a spec five
orders finer. That remains open and this does not close it.

---

## The idle floor is fine to ~10^4 live nodes and dominates past ~3x10^4

**Noticed:** measuring whether a town, then a city, then a forest fits.
**Where:** `engine.rs` — `survey`, `coast_to`, `evolve_matter`,
`record_histories`, each walking every live node every frame.

**This entry has been wrong twice and is now measured properly.** It first
claimed 5–25 us a node growing as n^1.5, from a probe that averaged frame time
across frames containing a Barnes-Hut solve over the root's 16,384 bodies and
divided by node count. It was then retracted outright, which over-corrected: the
retraction generalised one clean number at 8,193 nodes into "the floor is fine".

Measured across a range, counting **only frames where the plan accepted no
task**, so the solve is excluded rather than averaged in:

```text
live nodes  idle frame   us / node   share of 50 ms
       513    0.181 ms      0.3528            0.4%
      2049    0.882 ms      0.4304            1.8%
      8193    3.694 ms      0.4509            7.4%
     32769   27.149 ms      0.8285           54.3%
    131073  144.148 ms      1.0998          288.3%
```

**Flat at about 0.45 us a node to 8k, then the per-node cost itself climbs** —
0.83 us at 32k, 1.10 us at 131k. Overall n^1.32 from 8k to 131k, not the n^1.5
first claimed and not the linear the retraction implied.

**So the answer depends on the scene, and the wall sits between 10^4 and 10^5:**

* A town, or a city street bounded by an interest volume — a few thousand live
  nodes — costs under 2% and is free.
* **A forest of 10^4 individually promoted trees costs 9–18% of the frame doing
  nothing.** Affordable, and worth noting a forest does not normally need it:
  growth runs on the aggregate, so a forest is one node until the trees matter
  separately.
* **10^5 live nodes is 288% of the frame before anything happens.** A city that
  keeps that many nodes live is not viable, and the fourth axiom is what is
  supposed to stop it.

**What is not measured, and must not be guessed:** why the per-node cost climbs.
`survey` dominates at 8k (1.5–1.8 ms against `coast_to`'s 1.0 and
`evolve_matter`'s 0.4), and one candidate is memory rather than algorithm —
131k nodes at 576 bytes is 75 MB, far past any L3 — but the parent-chain walks
in `time_rate_of` are another and the two are distinguishable only by measuring.

**The fix is the lazy-catch-up rule** `docs/PLAY.md` §5A.5a states: nothing that
can be advanced in closed form is advanced by ticking, and the span settled is
the node's own elapsed proper time. It takes the floor from per live node to per
node being looked at.

**Trigger:** a scene that holds more than about 10^4 live nodes at once. A
forest with every tree promoted reaches it; a crowded city square plausibly
does. Below that this is not urgent, and the reason to do it anyway is
correctness rather than speed — ticking cannot cover a three-year absence
however cheap each tick is.

---

## The frame's cost model under-estimates a large gravity step, increasingly

**Noticed:** diagnosing a retracted claim about per-node cost.
**Where:** `budget.rs` — `SolverKind::cost`, `n * n.max(2.0).log2() * 1.4` for
Gravity. `solvers/gravity.rs` for what it actually does.

Measured, the one task the budget accepts on those frames:

```text
bodies    estimate      actual    ratio
  2048   23,552 us   48,076 us     2.0x
 16384  188,416 us  738,559 us     3.9x
```

The model is optimistic and **the error grows with n**, so it is not a wrong
constant — the shape is wrong, or the constant was calibrated at a body count
far below where it is being used. `docs/PERFORMANCE.md` says the constants come
from its own measurements, so either those were taken at small n or something
has changed since.

**Why it matters.** The budget is the engine's central promise: frame rate is
the invariant and detail gives way. A cost model that under-reads by 4x means
the knapsack fits work that does not fit, and the frame overruns instead of the
detail giving way. It is the one number the whole scheduling argument rests on.

**Trigger:** pulled. Re-derive the Gravity coefficient against measured
Barnes-Hut cost across three decades of body count, and make
`budget.observe_frame` — which already sees planned against actual — report the
ratio so a drift like this cannot go unnoticed again.

---

## One task can be fifteen times the frame budget, by design

**Noticed:** the same diagnosis.
**Where:** `budget.rs` — the plan takes the best task even when it exceeds the
whole frame. `docs/DESIGN.md` §3.7 states the rule outright: "A plan that
accepts nothing is worse than one that runs late, so the best task is taken even
when it costs more than the whole frame."

Measured: at 8,193 live nodes the plan accepts *one* task about every fifth
frame — a `Step` on the galactic root — and that frame takes 738 ms against a
50 ms target. Every other frame is 3.7 ms.

The rule is right for a single-observer explorer, where a long frame is a pause.
It is wrong for the play space: at one second per second a 738 ms frame is a
fifteen-frame hitch and the world falls 0.69 s behind real time in a single
step, which is exactly the "staler world" `docs/PLAY.md` D1 predicts and does
not want to arrive in lumps.

**The fix is task splitting, not a smaller budget.** A Barnes-Hut step over
16,384 bodies is divisible — half the bodies this frame, half the next — and a
task that can be cut into frame-sized pieces never forces the choice the rule
was written to resolve. What cannot be split (an indivisible solve) should still
be taken, so the rule survives for the case it was written for.

**Trigger:** the first play-space scene under a fixed clock, since D1 removes
the option of absorbing the overrun by slowing time. Not before task granularity
is decided, because splitting changes what a `Task` is.

---

## ~~A node can hold eight substances, and the play space needs more~~ — sized, and the loss is no longer silent

**Noticed:** asked whether "a galaxy is not made of anything you could put in a
beaker" — the comment justifying chemistry as a sparse side table — survives
contact with the play space.
**Where:** `chem/registry.rs` — `MIXTURE_SLOTS`, `Mixture::add`. `engine.rs` —
`environment_at`, `set_mixture`.

**It does not, and the justification is currently self-fulfilling.** Grepped:
`set_mixture` is called from **nowhere in `src/`** — only from tests. Chemistry
is entirely author-supplied, so it is rare because nothing makes it, not because
matter is rarely made of anything. In a galaxy-shaped node population the
sparsity claim is true; in a Continuum-shaped one almost every node of note is a
specific substance, and the claim inverts.

**Three consequences, one of them measured.**

**1. The slot count is a hard ceiling on what a node may contain.**
`MIXTURE_SLOTS` is 8, and `Mixture::add` displaces the smallest pool when full,
or returns `false` if the newcomer is smaller. Nothing at the call sites checks
that return. Measured, twelve equal substances into one mixture:

```text
   added   accepted   total fraction           lost
       8          8         0.666667       0.000000
       9          8         0.666667       0.083333
      12          8         0.666667       0.333333
```

A third of the *described* fraction, gone. **Corrected from this entry's first
version, which implied conservation was broken: it is not.** Mass, elemental
composition and chemical energy are scalars on `Matter`, and a `Mixture` is a
descriptive overlay whose fractions sum to `explained` — so dropping a pool
lowers `explained` honestly and the node is merely less described.

What is wrong is that *which* description survives depends on insertion order,
and that `add`'s `false` return — documented as telling the caller its trace
species did not make the cut — is read by nobody. And underneath both:
**mixtures do not aggregate at all.** `Composition::blend` exists for the
elemental account; there is no `Mixture` equivalent, and neither `summarise` nor
`coarsen` touches chemistry, so a wood cannot know its trees contain sugar.
`docs/PLAY.md` §5A is the plan for that.

A room as one node is already over: air is five substances (N₂, O₂, Ar, CO₂,
H₂O), a wooden table three (cellulose, lignin, water), a beaker of brine two.

**2. At Continuum, the `Mixture` is what the matter is, and `Composition` is its
summary — which is the wrong way round today.** Eight lumped elements cannot
tell brine from sodium metal plus chlorine plus water. Speciation is the fine
account and the elemental one is what survives coarsening, which is exactly
`sample` and `summarise` — a star has no use for salt. The engine already gates
a `Matter` field's meaning on tier for the same kind of reason (`evolve_matter`
does it for temperature), so this is an existing pattern applied to an account
that does not yet use it, not a new mechanism.

**3. `water = 1.0` is honest now and a lie later.** `environment_at` falls back
to unlimited water for a node with no mixture, and says so plainly: "there is
nothing to measure and the fallback is unlimited". Correct while chemistry is
author-only. The moment terrain generates patches that ought to be dry, every
undescribed node silently reads as infinitely wet, and "measure, never be told"
fails quietly. `Program::Terrain::substrate()` already returns a `Composition`;
producing a `Mixture` instead would let a patch measure its own water. That is
the direction D11 points anyway.

**The decision this forces, and it is not made here.** Either `MIXTURE_SLOTS`
grows to whatever a play-space node needs, or node granularity is constrained by
what a node can describe. The second is more interesting and more in keeping:
the engine already refines a node whose *gravity* has stopped being
representable (`needs_refinement`, on the Jeans length), and "refine when the
representation can no longer describe the contents" is the same rule in a
different account. It would make chemistry a driver of subdivision rather than a
passenger.

At minimum, `Mixture::add` must stop losing mass silently — a full mixture is
the signal to split the node, not to quietly forget the smallest thing in it.

**Trigger:** pulled for the silent loss, which is wrong now. The granularity rule
is Phase 2, when terrain starts generating nodes that are made of something.

**Done, as far as Phase 2 goes**, and the parts that are not are named rather
than left implied.

`MIXTURE_SLOTS` is **twelve**, sized against `PLAY.md` §5A.4's own worked case —
a room is air (five), a wooden table (three) and a beaker of brine (two), so
ten — with margin so the granularity rule fires on a real judgement rather than
an off-by-two. Twelve rather than sixteen because the cost is linear and now
universal: a `Pool` is 16 bytes, so a `Mixture` is 200 rather than 136, and it
sits on every `Matter` instead of in a side table. Measured: `Matter` 248 -> 456,
`Node` 1,136 -> 1,344.

**The silent loss is gone**, which was this entry's "at minimum". `Mixture::blend`
is what `coarsen` uses to combine a promoted child's description with its
parent's, and it returns the mass fraction it could not fit. Which pools survive
is decided by **size** rather than by the order the caller happened to add things
in, which was the other half of the complaint. `TreeStats::over_described` and
`worst_description_lost` carry it, on §3.7's precedent that a node crossed by its
own ensemble says so rather than doing it quietly.

**Point 2 is answered by D17**: the `Mixture` is on `Matter` now and travels by
`promote` and `coarsen` like every other conserved quantity, so the elemental
account really is the summary and the speciation really is the fine account.

**Point 3 is not.** `environment_at` still falls back to unlimited water for a
node with no mixture, and that is still honest only while most nodes have none.
D17 makes every node *able* to carry one; it does not make terrain generate one,
which is Phase 4's `Program::Terrain` work. `Matter::is_described` exists so a
caller can tell "measured and dry" from "never measured", which is the
distinction the fallback needs and did not have.

**Still not built: §5A.3's derived merge criterion**, which would blend two pools
into an interned substance rather than drop the smaller, and §5A.4's
subdivide-when-it-will-not-fit. `over_described` is the measurement standing
where that rule goes. **Trigger:** when a real scene drives it non-zero.

---

## ~~Identity allocation depends on the frame budget, and reaches the save file~~ — done

**Noticed:** asked whether an id reallocated across a sample/summarise cycle
would break the bit-identical regeneration rule.
**Where:** `engine.rs` — `advance_node` calls `issue_identity` to key its clock,
and which nodes `advance_node` runs on is chosen by the frame budget from a
wall-clock allowance.

**It does not break regeneration.** `sample` is deterministic in
`(matter, spec, world_seed, path_key, epoch)` and takes no identity, so
regenerated bodies are bit-identical whatever the ids are. That axiom is intact.

**It breaks save determinism.** Measured, same scenario, forty frames, varying
only the wall-clock budget:

```text
 budget us  next_entity     bytes             checksum
     50000            2    182196     784948e4801d0fc7
      5000            5    182196     ce4db2dd4b519220
       200            5    182196     ce4db2dd4b519220
         1            5    182196     ce4db2dd4b519220
```

Same length, different contents. A slower machine advances a different set of
nodes, issues a different number of identities, and saves a different world.
`next_entity` is persisted, so the divergence is durable.

**Why it matters more than it looks.** `docs/PLAY.md` D10 makes the input log
mandatory and states world state as `f(world_seed, ordered input log)`, and D1's
hindsight replay rests on the same thing. The moment an input names an
`EntityId` — "actor picks up entity 4471" — a replay on a different machine
binds that name to a different node. The failure would be silent and would look
like a physics bug.

**The cause is one line in the wrong place.** Identity is issued on a path the
*scheduler* drives rather than a path something *happened* on. A clock is
bookkeeping the budget created, not a fact about the world.

**Three ways out**, and they are not equivalent:

1. **Issue only on paths where something happened** — chemistry set, environment
   authored, node pinned, actor interacted. Clocks and histories stop naming
   nodes. This is what `docs/PLAY.md` D2's economy argument already claimed was
   true ("only pinned, structured or touched things need identity") and what the
   code does not do. Leaves the question of what keys `clocks` and `histories`.
2. **Derive the id from the address** — deterministic, no counter, and it
   re-breaks the thing D2 exists to fix, because a move would rename the node.
3. **Make the counter deterministic some other way** — for instance advancing it
   only on `epoch` bumps. Keeps clocks named, at the cost of a second rule about
   when ids may be handed out.

**Fixed** in Phase 1, by (1). Identity is issued only where something happened —
chemistry set, an environment authored, a node pinned. `clocks` and `histories`
went back to being keyed by address, which is what they should always have been:
both are bookkeeping the frame budget creates, and `persist.rs` already listed
them as transient. `World::reparent` migrates those two along with the identity
index, because a clock that changed rooms did not un-tick.

Measured after: `next_entity` is 1 across a 50,000x range of budget, where it
was 2 to 5 before. `tests/reparent.rs` asserts it, with a control that the
budget really did change which nodes were advanced — the clock count still moves
from 1 to 4, which is the point of keying clocks by address.

**One correction to what this entry first claimed.** The save file still differs
across budgets, and that part is *not* a defect. The residual difference is the
world instant: 7.11e13 s against 7.63e13 s after forty frames, because a
generous budget resolves more nodes, resolving shortens the pace, and less
simulated time therefore passes. That is the documented pace mechanism working,
and `time_throttle` exists to report it.

It also makes an argument for `docs/PLAY.md` D1 that the plan did not make:
**a fixed one-second-per-second clock is what makes a save reproducible across
machines.** While pace follows what is being watched, two machines running the
same scenario legitimately reach different instants. Once the clock is fixed,
they do not, and the only remaining source of divergence would be a real bug.
That is worth having as a test the moment `PaceMode::Fixed(1.0)` is the default.

---

## Nothing prunes the identity index, or the clocks and histories it names

**Noticed:** reviewing the `EntityId` change, asked to justify it.
**Where:** `engine.rs` — `identities`, `clocks`, `histories`. Grepped: the only
`remove` on any of the three is the one `reparent` does, and it re-inserts.

`clocks` gains an entry for every node that is ever advanced and `histories` for
every node that ever reaches Causal residency. Neither is ever removed, so both
grow with the number of nodes a world has *ever had*, not the number it has.
A node coarsened away leaves its clock behind.

That is pre-existing. What the identity change did is add a third map with the
same shape: `identities` gains an entry whenever a clock or a history names a
node, and is pruned only when *persisting*, never at runtime. So the leak is now
about twice the size it was.

Measured on the drilled galaxy scenario after twenty frames: 8 live nodes,
4 clocks, 0 histories — small, because that scenario coarsens nothing. A world
that repeatedly refines and coarsens the same region is the case that grows, and
nothing here has measured one yet.

**The fix is one rule, not three.** An entry in any of the three is worth keeping
exactly as long as the thing it names can still come back — which is what
`forgettable` and `mixing_time` already decide for detail. Prune on the same
criterion, in one place, rather than three ad-hoc sweeps. Note that `identities`
must be pruned *last*: dropping an address-to-name entry while a table still
holds that name orphans the row, which is the exact failure the change was made
to prevent.

**Trigger:** the first long-running session, or the first world that cycles a
region in and out repeatedly — the play space does both constantly. Measure the
three map sizes over a few thousand frames of a region being entered and left
before deciding how aggressive the rule needs to be.

---

## A regenerable snapshot still grows with body count

**Noticed:** measuring checkpoint size for the Phase 0 probes.
**Where:** `persist.rs` — `put_node`, and whatever else `encode` walks.

`persist.rs` is explicit that materialised bodies of unpinned nodes are
transient, "regenerated from address and epoch, bit-for-bit", and the code
agrees: `put_node` writes bodies only `if n.pinned`. So a snapshot of a
regenerable world should not care how many bodies it stands for. Measured, it
does:

```text
unrefined world, 1 live node                            736 B
8 live nodes, spec.count 512,   28,512 bodies       118,211 B
8 live nodes, spec.count 4096,  32,096 bodies       132,547 B
8 live nodes, spec.count 16384, 44,384 bodies       181,699 B
```

One node costs 736 bytes; eight cost 118 KB and up, and the figure tracks the
body count of detail that is not being written. Something on the path scales
with what was materialised — a side table, a spec, or a per-node array — and
it is not the bodies.

**Why it may not matter yet.** 182 KB is small in absolute terms and nothing is
wrong with the world that comes back; `tests/persistence.rs` passes. What it
undermines is the *claim* that regenerable detail is free to store, which the
design leans on in several places and `docs/PLAY.md` D10 leans on for sharding.

**And now `PLAY.md` D15 leans on it too.** D15's argument for a composite being
one node with a recipe is a storage one — a house of ~50 parts at under 1 KB
rather than ~29 KB, a town at 500 nodes rather than 25,000 — and that argument
assumes regenerable detail is nearly free to store. These numbers say it is not.
Re-measuring is part of validating D15 rather than a separate errand.

**Trigger:** before anything sizes a shard or a replay buffer against "coarse
nodes are nearly free", and before the first world large enough for 14 KB of
overhead per drilled node to matter. Sooner if D15's cost case is to be trusted.

---

## A structure regrows what was removed from it

**Noticed:** asked, while planning the play space, whether a window taken from a
building stays taken when the structure is summarised and re-sampled.
**Where:** `morph.rs` — `record`, `compact_events`, `MAX_EVENTS`,
`render_branching`, and the four renderers that are not it.

`morph::Event` is the mechanism for this and its doc is explicit: "A branch
breaking must survive coarsening — the whole point of a structure is that its
history is visible — so the deviations are logged and replayed rather than
discarded." `Skeleton::site` is "a program-stable name for this part, so an event
can refer to it and mean the same thing after the structure is regenerated."

**The logging happens. The replay mostly does not.** Measured by severing sites
one at a time and re-rendering after each:

```text
program      intact   after 64 severed   after 65   severed sites present again
tree            512          0 parts      512 parts        33 of 65
coral           512          0 parts      512 parts        33 of 65
tower           296        296 parts      296 parts        64 of 64  (from the first)
wall            504        504 parts      504 parts        64 of 64  (from the first)
terrain         484        484 parts      484 parts        64 of 64  (from the first)
settlement      256        256 parts      256 parts        64 of 64  (from the first)
```

**Two defects.**

**1. `MAX_EVENTS` resurrects.** The cap is 64; on overflow `compact_events` drops
the oldest half, folding only their mean magnitude into `genome[7]`. A segment is
suppressed only while its own `Severed` event survives, so the sixty-fifth
severance brings back everything the first thirty-two named. A tree severed at
the trunk sits at zero parts through 64 events and returns to all 512 on the
sixty-fifth. `built` was decremented each time, so mass is correct and geometry
is not: the structure contradicts its own conserved state.

**2. Only `render_branching` honours a severance.** The skip is written once, in
the renderer `Tree` and `Coral` share. `render_tower`, `render_wall`,
`render_terrain` and `render_settlement` never consult `events`, so for them a
severance has no geometric effect at all — not after 64, but immediately. A
demolished wall section is back on the next regeneration.

**Why it has not been noticed.** Every structural test uses a tree, and every
tree test breaks far fewer than sixty-four joints. `tests/fragments.rs` and
`tests/topology.rs` exercise breakage against the branching renderer, where the
skip exists and the cap is never reached; nothing anywhere severs a wall and
re-renders it.

**The fix, and it is a rule rather than a bigger number.** Ask of each event what
`docs/PLAY.md` §5.8 asks of each deviation: did it change the conserved tuple? A
`Severed` did — mass left the structure — so it may never be compacted away, and
must either persist or be promoted into the program's own description so the
program stops generating that member at all. `Damaged` and `Suppressed` did not,
and those are the events that may legitimately merge into an aggregate. Raising
`MAX_EVENTS` only moves the wall; nothing derives 64 and nothing would derive a
larger number either. The second defect is simpler: the four planned renderers
have to consult the log, and the skip wants to live somewhere all six share.

Also worth cleaning up while there: `genome[7]` currently doubles as a sink for
compacted history, so a per-instance variation slot means two things at once.

**Trigger:** pulled. Any actor who can remove part of a building reaches both of
these immediately, and the second one on the first edit. `docs/PLAY.md` §5.9
carries the same measurement and the reasoning behind the fix.

---

## A node cannot split when its contents spread out

**Noticed:** chasing the overflow above, which is *not* an instance of it.
**Where:** `engine.rs` — `refine`, `coarsen`, `promote` change a node's
resolution; nothing changes a node's **extent** or its **count**.

The tree can make a node's contents finer or coarser in place, and it can push
a child up a tier. It cannot say "these bodies are no longer one neighbourhood"
and hand them to two nodes. A node's radius is fixed when it is created, so
contents that legitimately expand — a gas cloud, dispersing debris, an
explosion, anything with a positive velocity divergence — either stay inside a
radius that no longer describes them, or leave it and are tracked by a node
that claims a volume they are not in.

Every downstream consumer of node radius is then wrong by the same factor: the
SPH smoothing length `h = radius / count^(1/3) * 1.2`, the gravity softening,
the LOD's angular size, the volume query's node selection, and the neighbour
grid spacing. None of them fail loudly; they all quietly describe a
neighbourhood that has stopped existing.

**What it needs.** A spread measurement per node — RMS distance of bodies from
the centre of mass, or the principal axes of their second moment — evaluated on
the same cadence as the node itself. Three outcomes: within the radius, do
nothing; larger than the radius but still one clump, grow the radius and
re-derive everything that depends on it; genuinely bimodal, split into two
nodes, each with its own centre, radius and body list, and `summarise` the pair
back to the parent so the conserved quantities still add up. The inverse merge
belongs with it, or two clumps that fall back together stay two nodes forever.

**The measurement now exists.** `state::Spread` and `Tree::spread`, built for
Phase 1 on its own, since three separate things want it and it is to be written
once. Centre of mass, mass-weighted RMS, and the furthest occupant *surface*,
over bodies and promoted children together; `Spread::occupancy` is the ratio to
what the node claims. `World::advance_node` measures it on the node it has just
touched — the same cadence, and a cache line that is already warm — and reports
the worst in `Stats::worst_occupancy`, with the `PathKey` it was seen at.

The **principal axes are deliberately not built**. They are what *splitting*
needs, to say which way to cut, and none of the three Phase 1 consumers can use
them; a symmetric eigensolver written for a caller that does not exist is what
this list is for avoiding. They go with the split.

**A baseline, for whoever later wants a threshold.** A healthy node does *not*
come out at or below one. A Plummer sphere's tail legitimately reaches three to
four radii and the scenario shelf measures 1.5 to 4.0, so "outgrew its radius"
is not the fault signal — orders of magnitude are. The two faults already on
this list read 4.5x10^5 (the sampler's chemical-binding inflation) and
4.2x10^11 (the ladder flinging a node's bodies out), which is the separation to
design against.

**Still to build:** the three outcomes above, and the merge.

**It does not fix the overflow entry above.** That was checked. The bodies there
are flung across twenty orders of magnitude *inside one* `advance_node` call, so
a split evaluated afterwards would faithfully split corrupt state into two nodes
of corrupt state. The spread measurement is still worth having as the *detector*
for that class of fault — it fires on the first frame, in the node that caused
it, rather than in a hash function twenty tiers away.

**Trigger:** the first simulation whose contents are meant to expand and are
meant to be measured afterwards — a detonation, a vented compartment, an
ablating surface. Nothing built so far expands; every test either holds a bound
configuration or collapses one. This is the reason it has not bitten yet, and
the reason it will.

---

## A coarse node does not age

**Noticed:** asked how a planet revisited after a hundred years can be right
when `sample` regenerates it from a static seed.
**Where:** `engine.rs` — `coast_to`, `survey`, `react_all`. `tree.rs` — `refine`.

The seed is not the problem, and the answer to the question as asked is that
`sample` is deterministic in `(matter, spec, world_seed, key, epoch)`, of which
`matter` and `epoch` both move. Regenerating gives *a* sample of the node's
**current** matter rather than the *same* sample — which is exactly the
rule `tree.rs` already states, that past a mixing time a stored sample is "no
longer *that* state, only *a* state" — and detail somebody touched is exempt
anyway, because it is pinned and comes back verbatim from the store.

**The problem is what moves `matter` while the node is coarse.** Measured:

* `coast_to` advances **`motion` only** — position, velocity, orientation.
  Not temperature, not composition, not internal energy.
* `TaskKind::Grow` runs on any node with a `morphology`, materialised or not.
  The comment there is right and is the model for everything below: "growth
  advances whether or not anything is materialised — in fact especially when
  nothing is. This is the payoff of the matter representation."
* `react_all` runs on any node with a `mixture`, materialised or not, on the
  matter's temperature.
* `TaskKind::Step` needs bodies, so it does nothing for a coarse node.

So **a coarse node evolves if and only if it has a morphology or a mixture.**
Everything else is frozen but moving. A planet with neither coasts a century
and comes back at the same temperature, with the same composition and the same
internal energy, having only changed position.

The sharpest instance: `matter.luminosity` is computed from Stefan-Boltzmann
and *read* — for illumination in `environment_at`, for flux in `observe.rs` —
but nothing anywhere subtracts `luminosity * dt` from `internal_energy`. Every
star in the world radiates into every scene and never spends anything.

**What it needs** is a *matter evolution* law: `advance_matter` standing to a
node's `Matter` as `advance_node` stands to its bodies, run from `survey` on
the same "materialised or not" basis growth already uses. The candidates are the processes that are slow, monotone and
depend only on the matter itself — radiative cooling, radioactive decay of the
composition, tidal and orbital evolution, accretion and mass loss. Structurally
this is the trick growth already proves works, applied to the quantities growth
does not own; the cost argument is the same one, that 10^4 coarse nodes cost 10^4
ODE steps whatever they stand for.

**How this differs from growth, which is the obvious thing to mistake it for.**
They look alike — both advance a coarse node, both write `internal_energy`,
`radius` and `luminosity` — and they are opposites in the way that matters.

| | `grow` | matter evolution |
|---|---|---|
| applies to | a node with a `morphology` somebody planted | every node, because it is physics |
| driven by | a `Program` — a developmental rule or a construction plan | a law with no parameters to choose |
| state | its own, in `Morphology`: segments, extent, stored energy. Not derivable from the matter | none beyond the matter itself |
| under coarsen/refine | persists; it *is* the state | must be idempotent, or looking changes the rate |
| fine-solver counterpart | none — the morphology is the model at every resolution | must agree with it on the conserved quantities |

`grow` is a program a node **runs**. Matter evolution is a law a node **cannot
escape**. That is why `grow` may own state that only it can produce, and why
matter evolution must own none: the moment it has private state, a node that was
materialised and re-coarsened evolves differently from one that was not, and
the observer has changed the physics.

**They will collide, and the collision is specific.** `grow_node` already
writes `internal_energy` (the thermalised share only), `entropy`,
`entropy_exported`, `radius` and `luminosity`, and its `GrowthStep` books
`energy_radiated` explicitly — "a leaf absorbs the whole solar flux and stores
about 0.3% of it; the other 99.7% leaves again". A cooling law that subtracts
`luminosity * dt` on every node would double-count that outflow on any node
with a morphology. So the two need one energy account between them, not two,
and `grow`'s existing `validate` is the right place to keep them honest.

Two more things to get right rather than assume. It has to be **consistent with
the fine solver**: a node cooled as matter for a century and then materialised must
land where materialising it and integrating for a century would have — at
least in the conserved quantities, which is the same guarantee
`summarise(sample(m)) = m` already carries. And it has to be **idempotent
across coarsen/refine cycles**, or a node that is looked at repeatedly evolves
at a different rate from one that is not, which would make observation change
the physics.

**Trigger:** the first thing whose *state* rather than position is expected to
differ after being left alone — a star that should have dimmed, a reactor slug
that should have decayed, a body that should have cooled. It is invisible until
someone looks twice and compares, and it is the mechanism by which "detail
exists where something is happening" stays honest: a node that is not happening
still has to *become* the thing you would find when you look again.

---

## ~~The test suite took twenty minutes and appeared to hang~~ — done

**Noticed:** asked why the tests are slow and why they sometimes never finish.

Two unrelated causes, and the second was not the tests at all.

**Slow: the test profile was unoptimised.** `Cargo.toml` set `[profile.release]`
and nothing for tests, so the whole numeric engine ran at `opt-level = 0`.
Measured per test, one accounted for most of it:

```text
    stepping_does_not_heat_a_node   405 s      11 scenarios x 200 solver passes
    every_scenario_steps             71 s
    the other three                 ~2 s
```

Full release was 60 s for that test against 405 s. But release also drops the
overflow checks that found the neighbour-grid bug, so the useful measurement
was the middle one: `opt-level = 2` **with** `debug-assertions` and
`overflow-checks` forced back on gives 66.7 s — within 10% of full release
while keeping every check. The suite went from about twenty minutes to **123 s**.

`overflow-checks` defaults *off* once `opt-level` rises, which would have
silently turned a panicking overflow into a wrapping one — the exact failure
that hid the grid bug in release. `tests/profile_guard.rs` asserts both flags
are still in force, so the saving cannot quietly cost the checks later.

**It also fixed `frames_stay_within_budget`**, which had been failing for the
whole session and which I repeatedly reported as pre-existing and unrelated. It
was pre-existing, and it was not unrelated: it measures wall-clock frame time
against a 50 ms target and the unoptimised engine took 1228 ms. Nothing about
the frame budget was wrong. The suite is now green.

**Appeared to hang: it did not — the watcher did.** No test loops on a
condition (the only `while` in `tests/` is bounded by `steps < 200`), every
tree walk terminates on a parent chain that ends at `NodeIdx::NONE`, `execute`
caps substeps at `MAX_SUBSTEPS` deliberately so the schedule does not depend on
machine speed, and the Postgres mutex recovers from poisoning rather than
deadlocking.

What did not terminate were shell loops of the form
`until ! pgrep -f "cargo test"; do sleep; done`, waiting for a command whose
name appears in the waiting shell's *own* command line. The condition is
permanently false and the loop never exits. Five were still running hours
later. Match on the test binary (`deps/<name>-<hash>`) or a PID, never on a
pattern the watcher itself contains.

The remaining honest cause is that `cargo test` has no per-test timeout, so a
405-second test is indistinguishable from a hang. At 123 s for the suite that
matters much less, and `timeout` around the command turns the remaining case
into an error rather than a wait.

---

## ~~A promoted child never feels a force~~ — done

**Noticed:** auditing what the engine does and does not couple, after the
neutron-bombardment question.
**Where:** `tree.rs` — `promote`, `sync_from_child`. `coords.rs` —
`Motion::advance`. `engine.rs` — `coast_to`, `advance_node`.

`promote` sets a child's `Motion` from the body it stands for — `offset =
body.pos`, `velocity = body.vel` — and sets `matter.momentum = ZERO`, because
the child's *frame* now carries the bulk motion. That is right. What is missing
is the other half: nothing ever changes `motion.velocity` again.

`Motion::advance` is `offset += velocity * dt` and a constant spin. Grepped
across the crate, `motion.velocity` is written in exactly two places —
`promote`, and decoding a saved world — and read everywhere else. **A promoted
node moves ballistically in its parent's frame for the rest of its life.**

Meanwhile the parent still holds the body it was promoted from, in the same
slot, and the parent's own solver goes on integrating it. So the same object has
two representations moving under different laws. Measured, promoting the most
massive body of a galaxy and running forty frames:

```text
    child node velocity  1.523667e4 m/s -> 1.523667e4 m/s   (unchanged)
    parent body velocity 1.523667e4 m/s -> 1.520674e4 m/s   (gravity acting)
    divergence after 40 frames: 6.7e18 m — 0.79 of the child's own radius
```

They reconcile only on `coarsen`, which calls `sync_from_child` — the single
call site — and overwrites the parent's body from the child, discarding
whatever the parent's solver did to it.

**Why it has not been noticed.** Every test either promotes and then looks
(where a fraction of a radius is invisible), or promotes and then coarsens
(where the sync hides it). Nothing yet promotes two siblings and expects them to
interact, which is the case where it becomes obvious: two vehicles in one frame
do not attract, collide, or perturb each other at all, because neither one's
node can be moved by anything.

**The fix is a direction to choose, not a line to write.** Either the parent's
body is the authority and the child's `Motion` is slaved to it each frame — cheap,
but then a promoted node cannot have its own dynamics; or the child is the
authority and the parent's body is slaved to the child, which is `sync_from_child`
run every frame instead of only at coarsen, and costs one write per promoted node
per frame. The second is more consistent with what promotion means everywhere
else — the child is the real thing, the body is the stand-in — and it is what
makes the parent's solver see the child's evolved position, which is what
sibling interaction needs.

**Done** in Phase 1, taking the second direction: the child is the authority and
the parent's body is slaved to it. `Tree::sync_children` runs at the top of
`advance_node` so the solver sees where the child actually is, and
`Tree::apply_body_forces` hands each child the velocity change its stand-in just
received — the half that did not exist at all.

Position is deliberately not copied back, because the child owns where it is.
So the two differ between syncs by the solver's own second-order term: leapfrog
integrates `x + v·dt + ½a·dt²` and the child coasts linearly. Measured over four
hundred frames that oscillates between 0.13 and 0.29 of the child's radius and
does not climb, where the old behaviour passed 0.79 in forty frames and kept
going. Bounded and reset rather than accumulated is the whole of the fix.

`tests/promoted.rs` holds it, each assertion verified against a mutation:
removing the force transfer, or removing the sync, fails the tests that name
them.

---

## Sparse discrete transport has no solver shape

**Noticed:** asked how "a region of water bombarded by neutrons" would be
handled.
**Where:** nowhere, which is the point. `solvers::for_tier` offers Gravity,
GravityHydro, Hydro, MolecularDynamics and Statistical.

**Every solver in the engine is a dense-interaction solver.** Gravity: every
body pulls every body. Hydro: every neighbour within `h`. MD: every neighbour
within `cutoff`. All three assume a particle interacts with everything nearby,
continuously, every step.

A neutron is the opposite. It is a Nuclear-tier object (10^-15 m) whose *mean
free path in water is centimetres* — Continuum tier, five tiers coarser than
itself. It crosses ~10^23 molecules' worth of matter without touching any of
them, then interacts once, discretely. Put it in a Continuum node as a `Body`
today and `for_tier` hands it SPH, which would give it a smoothing length of
centimetres and have it interact with everything continuously: modelling as
dense the one thing whose entire physics is that it is rare.

The same shape covers photons through a medium, cosmic rays, and any beam.

**What exists.** The water is fine: `chem` builds H2O from first principles, a
node carries it as a `Mixture`, `react_all` runs each frame. `Isotope::Neutron`
decays correctly with an 878.4 s half-life in `advance_statistical`. And the
delivery mechanism is already the right shape — `causal::Influence` is a
discrete cross-node event with a causally-ordered arrival time, and
`InfluenceKind::Radiation` already exists.

**What does not.** Grepped for cross-section, mean free path, capture,
moderation, attenuation, radiolysis: nothing.

* `solvers::nuclear` is entirely *stellar* — Gamow peak, pp-chain, CNO,
  triple-alpha, SEMF, opacity. Those are thermonuclear rates for a hot plasma
  in equilibrium, not (n,elastic) or (n,gamma) for a beam.
* `advance_statistical` moves bodies ballistically and samples decay. It never
  asks what a body is passing *through*.
* `apply_influence` can only `add_heat` and add momentum. Every `InfluenceKind`
  collapses to those two, so `1H(n,gamma)2H` has no route in.
* `Interaction::Inject` takes a `Composition`, which is the eight
  `CoarseElement` buckets. There is no free neutron in that account.
* `chem::react` moves mass between *phases only*, and its own module doc says
  that is what makes it safe to run everywhere: "the elemental account cannot
  move, because dissolving a salt does not transmute anything." Radiolysis
  produces H., OH., e-aq, H2, H2O2 — new species. It is precisely what that
  pass is built not to do.

**What the shape would be**, and it fits the engine's principles rather than
fighting them:

1. **A travelling particle is an influence in flight, not a body.** The mailbox
   already models a scheduled discrete event with a causally correct arrival
   time. What is missing is that `apply_influence` cannot change composition.
2. **Schedule by the mean free path, not the particle's size** — the same
   lesson as `drill_to`. Lambda is a *length*, and `Tier::containing(1e-2)`
   puts the flight at Continuum where it belongs.
3. **Most cross sections derive.** Sigma = n·sigma, and the number density
   comes from `Matter`'s mass, radius and composition, which every node already
   carries. Elastic scattering off hydrogen is two-body kinematics — a neutron
   loses half its energy per collision with a proton, which is *why* water
   moderates, and ~18 collisions take 2 MeV to thermal. Only radiative capture
   needs real data (1/v plus resonances): a legitimate derived-and-stored
   shortcut rather than a table of everything.
4. **Almost no molecule ever needs to exist.** 10 cm^3 of water is 3x10^23
   molecules; a 10^9 n/s beam over a microsecond is ~10^3 interactions. Each is
   one place to `drill_to`, do the physics, and let go. This scenario is close
   to the ideal advertisement for the design.

**Trigger:** any radiation, any beam weapon, any dosimetry, any reactor, and the
"radiation damage to tissue" item that has been on the informal list since the
morphology work. Also the honest answer to "what happens if I shine this at
that", which is a question the play space will ask constantly.

---

## What else is not coupled — an audit

**Noticed:** asked to check for other interactions that are not handled, having
found the two above.
**Method:** grepped for each mechanism by name across `src/`, then checked
whether what was found is *read* by anything rather than merely stored.

Recorded together because they share a cause: the engine models what happens
*inside* a node very well and what happens *between* nodes barely at all. The
only inter-node couplings that exist are the causal influence mailbox, the
parent's luminosity reaching a child through `environment_at`, and
`sync_from_child` at coarsen.

| interaction | state | note |
|---|---|---|
| **Contact / collision** | only `Fragment` | `Fragment::contacts` handles a detached piece against the structure it fell from and the ground. Two vehicles cannot collide; nothing can rest on anything. |
| **Electromagnetism above molecular tier** | absent | `md.rs` has a shifted-force Coulomb term, so charge is real at Molecular/Atomic tier. `Body::charge` is read by **no other solver** — an ion at Continuum tier has a charge that does nothing. |
| **Magnetism** | dead state | `Matter::magnetic_energy` is persisted and initialised and **read by nothing**. No field, no Lorentz force, no MHD. |
| **Radiative transfer** | one hop only | `environment_at` gives a child the flux from its *parent's* luminosity. No sibling-to-sibling light, no shadowing, no opacity attenuation. The opacity functions exist but serve stellar burning. |
| **Heat conduction between nodes** | absent | A hot node beside a cold one never equilibrates. Conduction exists *within* a structure (`Mechanism::ConductedEnergy`) and within SPH neighbours, never across a node boundary. |
| **Mass diffusion between nodes** | absent | Chemistry runs within one node's `Mixture`. Nothing crosses a boundary, so a solute cannot spread from one node into the next. |
| **Free surfaces and interfaces** | absent | A node holds one `Mixture` with phase fractions, not a *surface* between them. So no buoyancy, no capillarity, no surface tension, no sloshing, no "water level". |
| **Reactions beyond phase change** | absent | Combustion, corrosion, acid/base, redox. `react` may not move the elemental account by construction. |
| **Transmutation by reaction** | stellar only | `burn` changes composition through fusion channels; nothing else can. |
| **Friction between bodies** | absent | `MdParams` has Langevin friction (a thermostat, not contact) and `Mechanism::FlowDrag` is a fluid load on a structure. Sliding, rolling and static friction between two objects do not exist. |

**The pattern worth naming.** Several of these are one mechanism wearing
different clothes: conduction, diffusion, radiative exchange and contact are all
"two adjacent nodes exchange a conserved quantity across their shared boundary".
The engine has no notion of adjacency at all — nodes know their parent and their
children, and `separation` can measure any two, but nothing enumerates *what is
next to what*. A neighbour relation between sibling nodes is the single missing
primitive underneath four of the rows above, and it is worth designing once
rather than four times.

**Trigger:** each row has its own, but the neighbour relation should be designed
before the first of them is built, or it will be built four incompatible ways.

---

## A detached fragment is neither a body nor a node

**Noticed:** asked whether a broken tree branch should become its own node.
**Where:** `engine.rs` — `World::falling`, `drop_fragments`, `MAX_FALLING`,
`MAX_FALL_SECONDS`. `solvers/structure.rs` — `Fragment`, `detach`.

A branch that snaps lives in a third place: `World::falling`, a
`Vec<(NodeIdx, Fragment)>` keyed by the node it fell from. Its mass stays in
that node — `GrowthStep::mass_detached` states the policy outright, "mass that
left the structure but stayed in the node; a fallen limb is litter, not an
absence" — while the `Fragment` carries its own topology and tumbling dynamics
for at most twelve seconds and is then dropped.

**The policy, when this is built: promote on leaving the volume, not on
detaching.** Detachment is a topological fact; deserving a node is a spatial
one. A branch that settles two metres down is still in the tree's
neighbourhood — same frame, same volume, litter rather than structure — and
promoting it buys nothing. A branch blown off a cliff, or knocked off a vehicle
in flight, has left the parent's neighbourhood, and that is what a node is for.
The test is the spread measurement that node splitting already needs, so the
two share machinery.

**Five things that are wrong today regardless.**

1. `MAX_FALL_SECONDS = 12` is a terrestrial assumption. A piece falling from a
   200 m tower, or in low gravity, is written off in mid-air and its mass
   silently reverts to being inside the parent.
2. `MAX_FALLING = 24`, and past it the loop simply `break`s. A vehicle at the
   top of the stated play-space scale shedding debris gets twenty-four pieces.
3. A fragment can only strike the node it fell from — `contacts` is called with
   that node's bodies. A branch cannot land on the next tree, on a person, or
   on a vehicle. That is the missing adjacency relation again.
4. While it falls, its mass is in the wrong place: the parent's `matter.mass`,
   `com` and `radius` are unchanged, so gravitationally and thermodynamically
   the branch is still inside the tree while visually ten metres away.
5. Fragments do not survive a save — `persist.rs` says in-flight solver state
   is not stored. Mass is safe because it stayed in the node, but a save taken
   mid-collapse loses the debris.

**Two blockers, both already in this file, and the order matters.** A promoted
fragment would *float, not fall*, because a promoted node never feels a force.
And `promote(node, slot, spec)` takes a single slot, while a fragment is a *set*
of members — expressing it needs the sibling-from-a-subset operation that node
splitting is about. Fix those two and this stops being an architectural question
and becomes one line of policy in the spread check.

**Trigger:** the first debris that has to outlive twelve seconds, land on
something other than its own parent, or survive a save.

---

## `World` accumulates what could not decide whose property it was

**Noticed:** asked whether the world is not just a node with children, and
whether the design had crept.

**It has, in one specific way, and the core has not.** `Tree { nodes, root }`
*is* the world, the root node is the universe, and `World::new` paces to the
root immediately so a world is watchable the moment it exists. That much is
intact. What `World` adds falls into three groups and only the third is creep.

**Legitimate — the apparatus for simulating.** `budget`, `gate`, `observers`,
`time`, `pace`, `pace_mode`, `stats`, `gpu`, `history_depth`. A universe does
not have a frame budget; the process watching it does. These belong off the
node, and the useful test is: *would this still be true if nobody were
simulating?*

**Defensible but drifting — sparse per-node data in side tables.**
`histories`, `clocks`, `environments`, `mixtures`, all `HashMap<PathKey, _>`.
Each has a real justification, and `mixtures` states it well: `Matter` is `Copy`
and about two hundred bytes, a `Mixture` is another hundred and forty, and most
nodes have no chemistry at all — "a galaxy is not made of anything you could put
in a beaker". Keying by `PathKey` so they survive coarsening is the same trick
pinned detail uses, and is right.

What has drifted is that a node's identity now lives in five places and
**nothing enumerates them together**. No table is ever pruned: grep finds no
`remove` or `retain` on `histories`, `clocks` or `environments`, and `mixtures`
has exactly one, in `set_mixture`, for the unrelated case of a mixture emptying.

**This entry first called that a leak and prescribed sweeping the tables when a
node dies. Both halves were wrong**, and checking before implementing is what
caught it.

*There is no death path.* `release_subtree` has exactly two callers — `coarsen`,
and its own recursion. Every release is a coarsening, and `PathKey::child` is a
pure deterministic hash of the parent key and slot, so a re-promoted path
recovers the same key. Surviving release is the entire point of keying by path,
and sweeping there would destroy a coarsened node's chemistry and history the
moment nobody was looking at it.

*And the four tables are not one kind of thing.* `persist.rs` sorts them:

| table | persisted | what it is |
|---|---|---|
| `mixtures` | yes | what a node is *made of*. Authored, durable, must survive |
| `environments` | yes | authored per-node overrides. Same |
| `clocks` | no | "rebuilt as the simulation runs" |
| `histories` | no | retained observations for viewing at a distance |

The two that persist are durable state, and their growing with the number of
places the world has been given chemistry is *correct* rather than a leak. The
two that do not persist are process state, and those are the ones with no
eviction at all — `history_depth` bounds each `History`'s length, nothing bounds
how many exist.

So the real defect is narrower than first written: **nothing enumerates the
tables**, so a sixth can be added and silently missed by anything that does
need to sweep; and **`clocks` and `histories` have no eviction policy**, which
their not being persisted is the tell for.

**Creep — `falling` and `shaking`.** These are world state: a branch really is
falling whether or not anyone is simulating it, which is exactly the test above,
and it fails. They sit on the simulator as `Vec<(NodeIdx, T)>`, capped at
twenty-four, unpersisted, and keyed by **`NodeIdx` rather than `PathKey`** while
`release_subtree` recycles indices through `Tree::free`. A fragment that
outlives its node would then deliver its impulses to whatever node next takes
that slot. *Not demonstrated* — it needs a promoted node damaged, released while
its pieces are still in the air, and its index reused inside twelve seconds —
but the shape is wrong regardless, and it is why a fragment can only strike the
node it fell from: its identity is a *pair* rather than a thing in the world.

**What to do, corrected.**

1. **Enumerate the tables in one place** — a struct or a macro listing them — so
   a sixth cannot be added and silently missed. This is the part that holds
   whatever policy is chosen later, and it is cheap now.
2. **Present them as node state.** Accessors (`World::mixture_of`,
   `environment_of`) rather than public `HashMap` fields. The side table is a
   memory optimisation and should read like one rather than like a separate
   concept a caller has to know about.
3. **Give `clocks` and `histories` an eviction policy.** Not the persisted two.
   The engine already has the right idea for it — a path nobody has observed
   for a mixing time does not need its history kept — but it is a policy
   decision rather than an obvious default, so it wants choosing rather than
   assuming.
4. **Move `falling` and `shaking` into the tree**, which is the fragment entry
   above and wants doing with it.

**Trigger:** (1) and (2) whenever the next side table is added, which is when
the cost is lowest. (3) when a session runs long enough for unobserved paths to
outnumber observed ones.

---

## There is no terrain — and gravity for debris was a constant, which is now fixed

**Since Phase 1**, two of the three things below have moved and the entry is
kept because the first has not. `G_EARTH` is deleted: `drop_fragments` derives
the field from `Tree::gravity_at` on the node the piece is falling inside, so
debris weighs what the place says. And a falling piece is no longer restricted
to the structure it came off — D3's adjacency index gives the parent's siblings
as candidates, each with its offset, which is what makes
`a_branch_lands_on_the_next_tree` pass. **What is unchanged is the first
paragraph: there is still no ground.** Everything a piece can land on is a
structure with a topology, and a planet's surface is not one yet. That is
Phase 2.

**Noticed:** asked whether a falling branch can hit the ground, and whether
terrain exists.
**Where:** `engine.rs` — `ground_of`, `drop_fragments`. `solvers/structure.rs`
— `G_EARTH`.

**There is no terrain, no heightmap, no surface of any kind.** "The ground" is
`ground_of`, which returns a single `z` — the lowest unsupported joint of the
structure itself, in that structure's own frame. It exists so debris does not
fall through the floor of a structure recentred on its own centre of mass, and
for that it is exactly right. It is not a world feature; it is a property of one
structure, and a falling piece can strike only that plane and members of the
structure it came from.

`DESIGN.md` states the debris half honestly — "debris collides with the
structure it fell from and with the ground, not with other debris" — and
`PHYSICS.md` explains why the ground is not at zero. So this is a documented
frontier rather than drift. What is not written down is the consequence: a
branch cannot land on the earth, on another tree, on a person, or on a vehicle,
because none of those are things it can be tested against.

**Second, smaller, and undocumented — and since fixed, see the note at the top
of this entry:** `drop_fragments` loaded every falling
piece with `st::G_EARTH.scale(m)`, and `G_EARTH` was the constant
`(0, 0, -9.80665)`. Debris therefore falls at Earth gravity along its own
structure's negative z wherever the node actually is — on a ship under thrust,
in orbit, on a body of any other mass. The engine computes real gravitational
fields at every other tier and then ignores them here.

**What terrain should be, given the design.** Not a special case: a planetary
surface is a `Continuum` node with a topology, which is a thing the engine can
already build. What is missing is not terrain as a type but the ability for a
falling piece to be tested against a node that is not its own parent — the
adjacency relation from the audit. Terrain is the first customer for it rather
than a separate feature.

**Trigger:** the first scene where something falls onto anything other than the
structure it came from. That is essentially the whole play space.

---

## ~~Nothing can change parent~~ — done, with one thing it exposed

**Noticed:** asked whether a branch thrown into space could leave the world and
collide with things outside it.
**Where:** `tree.rs` — `Node::parent` is written in exactly two places, `NONE`
for the root and `parent: i` inside `promote`. Nothing else assigns it.

**A node's place in the hierarchy is fixed at creation for life.** There is no
re-parenting, no hand-off between frames, no transfer. There is also no
"outside the world" to leave to — everything descends from one root — so the
question's second half is really "can it reach a different part of the
hierarchy", and the answer is also no, for the same reason plus the missing
adjacency relation.

So "throw a branch into space" is not expressible. As it rises it would have to
leave the tree's frame, enter the planet's, and then the system's; each of those
is a change of parent, and none of them can happen.

**Why this was not a mistake.** The engine was built galaxy-downward, and
galactic containment does not change on any timescale that matters: a star does
not leave its cluster during play. Building the hierarchy static was right for
what had been built. It stops being right at the stated play scale, where
containment changes constantly — pick up a rock and it enters your frame, throw
it and it enters the ground's, walk into a ship and you enter the ship's.

**The four gaps found in this session are one gap.** The tree is a *static
containment hierarchy*, and a game-engine platform needs it to be a *dynamic
spatial index*:

| gap | what it prevents |
|---|---|
| a promoted child never feels a force | anything moving under physics once promoted |
| a node cannot split | making a node out of part of another |
| no adjacency relation | finding what is next to what |
| **nothing can change parent** | **anything moving between frames** |

Re-parenting is the foundational one. It is what turns the tree from a record of
what owns what into an index of what is where, and the other three are much
easier to reason about once a node's parent is allowed to change.

**What it needs.** Moving a node between parents is a frame change, so its
`Motion` must be re-expressed in the new parent's frame — the engine already has
that arithmetic in `separation`, `offset_from` and `velocity_from`, which walk
to a common ancestor. Conservation must hold across the move: the old parent
loses the mass and the new one gains it, which is `summarise` in both directions
and is the guarantee `IDEMPOTENT_TOLERANCE` already covers. And the slot the
node occupied in its old parent's body list has to be vacated rather than left
holding a stale body.

**Trigger:** the first object that moves between containers, which is the first
object a player picks up.

---

## ~~A `PathKey` is doing two jobs, and re-parenting made it show~~ — done

**Closed** by `docs/PLAY.md` D2, in Phase 1. `ids::EntityId` is issued once and
never reused; `mixtures`, `environments`, `clocks` and `histories` are keyed by
it, so `World::reparent` no longer enumerates them. One index — `World::identities`,
address to name — is migrated instead, because a node discarded and rebuilt
recovers its name from its address and nothing else. `tests/reparent.rs` holds
the three properties: a move changes the address and not the name, a name is
issued once and asking does not issue, and the tables keyed by name are not
touched by a move.


**Noticed:** building `Tree::reparent`. **Where:** `ids.rs`, and every table
keyed by `PathKey`.

`PathKey`'s own doc calls it "stable identity for nodes across materialisation
cycles", and it is — as long as nothing moves. It is a rolling hash of the
child-index path from the root, so it answers two different questions with one
number:

* **Where is this?** It derives children, seeds the sampler, and is the address
  a node's procedural contents are generated from.
* **Which thing is this?** The ledger, pinned detail, chemistry and environments
  are all filed under it.

Moving a node changes the answer to the first and must not change the answer to
the second. `reparent` handles that by rekeying the subtree and migrating every
table, which is correct and is tested — but the identity is *reconstructed* at
the new address rather than carried, and two consequences follow that the
`phys-rehome` demo shows plainly:

**A round trip does not return the original key.** Put an object down where it
was picked up and it gets a new key, because it takes a new slot. Nothing
outside the migrated tables can recognise it as the same object — an external
reference held across a move is dangling, and there is no way to ask "is this
the thing I was holding".

**Vacated slots are permanent.** A slot cannot be reused or removed: `children`
is parallel to `bodies`, and a sibling's key is derived from its slot index, so
removing an element would renumber every sibling after it and change what each
one *is*. So the vacated body is zeroed in place, and a parent that has had a
hundred objects pass through it carries a hundred dead slots forever. Every
solver iterates them.

**The fix is to stop conflating the two jobs.** A node that has been
individuated — pinned, authored, moved — needs an identity that is *issued*
rather than derived, and stays with it wherever it goes. Procedural scenery
does not: for a node nobody has touched, the address genuinely is the identity,
which is what makes recipes work at all. That split is the same one the engine
already draws between regenerable and pinned detail, so it has a precedent to
follow rather than a new concept to invent.

With an issued identity, vacated slots also stop being permanent, because slot
index would no longer determine what a sibling *is*.

**Trigger:** the first thing that holds a reference to an object across a move —
an actor's inventory, a target lock, a command naming what to act on. Which is
to say, the first thing built on top of re-parenting.

---

## A node's tier is decided once and never revisited

**Noticed:** the moon/town test. The terrain patch prints as `Planetary` while
being 1.4 km across, which is `Continuum` by the table.
**Where:** `tree.rs` — `promote` sets `tier` from the promoted body's radius and
nothing sets it again.

`Tier::containing(body.radius).max(parent_tier)` runs once, at promotion. Every
later thing that changes a node's size leaves the tier behind: `emplace` and
`plant` both set `matter.radius` from the program's `extent`, growth changes it
every step, and `Morphology::extent` is explicitly the authority on how big a
structure is. So a node can be two tiers away from what its own radius says.

Measured: a surface patch promoted out of a moon inherits a body radius in the
planetary band, is then emplaced as terrain 1.4 km across, and stays
`Planetary`. `solvers::for_tier(Planetary)` is `GravityHydro`.

This is the same shape as the nucleons-marked-`Atomic` entry above, and the
same root: **a tier is derived from a radius at one instant and then cached**,
so the two drift apart the moment anything changes size. The difference is that
this one is benign so far, because a structural node's behaviour comes from its
morphology and topology rather than from `for_tier`.

**The fix is probably not to recompute it everywhere.** A tier that changed
under a node would change its solver, its timestep and its cadence mid-flight.
More likely: recompute it at the points that already change a node's size —
`plant`, `emplace`, and the growth step — and make `promote` share that path so
there is one place that decides.

**Trigger:** the first time a structural node's tier matters to something. It
does not today, which is exactly why it should be recorded rather than fixed on
the spot.

---

## An authored environment is all-or-nothing

**Noticed:** writing the biome tests, where a scenario wants to control the
light on a node and still have its water measured from what the node is made
of.
**Where:** `engine.rs` — `environment_at` returns `self.environments[key]`
wholesale when one is present.

An override replaces the entire `Environment`, so authoring any one field
silently discards the derivation of all the others. Setting `light_flux` to
stage a lit surface also pins `water` at whatever the author happened to leave
in the struct — usually `Default`'s `1.0` — and a patch that should have frozen
into a desert goes on measuring as well watered forever.

This was worse before `plant` and `emplace` took `Option<Environment>`: they
inserted an override unconditionally, so *every* structure in the world had one
and none of them ever felt the weather. That much is fixed. What remains is that
an author cannot say "this much light, and derive the rest".

**The fix** is a partial override — `Option` per field, merged over the derived
values rather than replacing them. It touches the persisted `environments`
table, so it is a wire-format change, and it wants doing before scenarios start
depending on the all-or-nothing behaviour.

**Trigger:** the first scenario that authors one field and is surprised by
another. The biome tests dodged it by building a star and dimming that instead,
which is more honest anyway — but it is a workaround, not a preference.

---

## A leaf cannot fall: nothing sheds one, nothing slows one, and nothing is under it

**Noticed:** asked directly for a test of a leaf falling and hitting the ground.
**Where:** `engine.rs` — `drop_fragments` and `ground_of`; `morph.rs`, where a
leaf is an area rather than a part; `solvers/structure.rs` — `Mechanism::FlowDrag`.

The scenario is one sentence — a leaf comes off a tree, flutters down, lands,
stays — and four separate things are missing. Recording which, because three of
them are wanted by far more than this and only one is Phase 2's by name.

**1. A leaf is not a thing.** Leaves enter the engine as `capture_area` and as a
line in the light budget; they are never parts. What a tree's body list holds is
members with joint radii, and then `sample_structured`'s "unstructured
remainder: litter, air, rubble" — a mass fraction of loose bodies with no shape
and nothing attaching them. So nothing can shed a leaf, because nothing has one.
What *does* come off a tree is a fragment: a detached run of members, produced
by `damage` or `shake` and tracked in `World::falling`.

**2. Loose contents feel no gravity.** The derived field —
`Tree::gravity_at`, cached as `Node.gravity` since D6 — has exactly three
consumers: `damage` (`engine.rs:2453`), `shake` (`2613`) and `drop_fragments`
(`2754`). All three are structural. A node's *loose* bodies are handed to the
tier solver instead, and `hydro` and `gravity` give them only their mutual
attraction, which across a six-metre node is 10^-8 m/s^2. A leaf modelled as a
loose body would hang in the air. Making it fall is a fourth consumer of a field
that already exists, and it is the same gap a dropped tool, a thrown stone and
settling dust all have.

**3. A falling thing feels no drag.** `drop_fragments` builds its load as `g·m`
per joint and nothing else. `Mechanism::FlowDrag` exists, takes the fluid's
density and a drag coefficient so nothing in it is specific to air — and is
applied only by the standing-load paths, never to a piece in flight. For a limb
that is nearly right. For a leaf it is the entire physics. The arithmetic, for a
0.3 g leaf of 30 cm^2 at `Cd` 1.2 in air at 1.2 kg/m^3 — terminal velocity
1.17 m/s:

```text
    fall      arrives (drag)   arrives (as coded)   too fast by   energy
    10 m      1.17 m/s, 8.7 s   14.0 m/s, 1.4 s         12x         144x
    30 m      1.17 m/s, 25.8 s  24.3 m/s, 2.5 s         21x         432x
```

`MAX_FALL_SECONDS` is 12, so the 30 m case is written off as litter in mid-air,
about 14 m up. A limb from the same branch arrives in 2.5 seconds and is fine.

**4. There is nothing under it.** `ground_of` returns the lowest unsupported
joint of *the structure being struck*, in that structure's own frame, and says
why: "what anchors a structure is what it is standing on". There is no world
surface anywhere in the engine. A leaf that misses the tree it came off has
nothing beneath it at all, and `MAX_FALL_SECONDS` is what catches it. This is
the debris end of "There is no terrain, and gravity for debris is a constant" —
whose second half is now done and whose first half is Phase 2.

**What the test should assert, once it can run.** Written down now, while the
reasons are fresh:

- It arrives at roughly terminal velocity rather than at `sqrt(2gh)`. That ratio
  is the whole point; everything else about the scenario is scenery.
- It lands on the terrain, not on the foundation plane of whatever structure
  happens to be nearest, and it stays there.
- The mass arrives: what left the tree is what reached the ground.
- A second leaf released from a different height reaches the *same* speed. That
  is the signature of a terminal velocity and it cannot be faked by a tuned
  constant, which is what makes it the assertion worth having.
- One control, or the test proves nothing: the same leaf with drag suppressed
  must arrive about twelve times faster. A landing test that still passes
  without `FlowDrag` is a test of gravity.

**Trigger:** Phase 2's terrain, for the ground to exist. Items 2 and 3 come due
earlier and for more callers — the first dropped or thrown object needs the
derived field on loose contents, and the first object whose shape matters more
than its mass needs drag in flight.
