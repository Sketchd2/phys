# Backlog

Deferred work, with the reason it was deferred. Not a roadmap — the roadmap is
the architecture review. This is the smaller list of things noticed while
building something else, parked deliberately rather than forgotten.

An entry earns its place by being **specific about the trigger**: what would
make it worth doing, or what will go wrong if it is not. "Would be nice" is not
a trigger.

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

## A promoted child never feels a force

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

**Trigger:** the first time two promoted things are meant to affect each other.
That is most of the stated play space, so this is nearer than its position in
this file suggests.

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

## There is no terrain, and gravity for debris is a constant

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

**Second, smaller, and undocumented:** `drop_fragments` loads every falling
piece with `st::G_EARTH.scale(m)`, and `G_EARTH` is the constant
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

## A `PathKey` is doing two jobs, and re-parenting made it show

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
