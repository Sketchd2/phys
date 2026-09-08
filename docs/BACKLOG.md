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
fail-fast stopped the run before reaching them.

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
meaningless until then.

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
**current** bulk state rather than the *same* sample — which is exactly the
rule `tree.rs` already states, that past a mixing time a stored sample is "no
longer *that* state, only *a* state" — and detail somebody touched is exempt
anyway, because it is pinned and comes back verbatim from the store.

**The problem is what moves `matter` while the node is coarse.** Measured:

* `coast_to` advances **`motion` only** — position, velocity, orientation.
  Not temperature, not composition, not internal energy.
* `TaskKind::Grow` runs on any node with a `morphology`, materialised or not.
  The comment there is right and is the model for everything below: "growth
  advances whether or not anything is materialised — in fact especially when
  nothing is. This is the payoff of the bulk representation."
* `react_all` runs on any node with a `mixture`, materialised or not, on the
  bulk temperature.
* `TaskKind::Step` needs bodies, so it does nothing for a coarse node.

So **a coarse node evolves if and only if it has a morphology or a mixture.**
Everything else is frozen but moving. A planet with neither coasts a century
and comes back at the same temperature, with the same composition and the same
internal energy, having only changed position.

The sharpest instance: `matter.luminosity` is computed from Stefan-Boltzmann
and *read* — for illumination in `environment_at`, for flux in `observe.rs` —
but nothing anywhere subtracts `luminosity * dt` from `internal_energy`. Every
star in the world radiates into every scene and never spends anything.

**What it needs** is a bulk evolution law: the coarse-state counterpart to the
solvers, run from `survey` on the same "materialised or not" basis growth
already uses. The candidates are the processes that are slow, monotone and
depend only on the bulk tuple — radiative cooling, radioactive decay of the
composition, tidal and orbital evolution, accretion and mass loss. Structurally
this is the trick growth already proves works, applied to the quantities growth
does not own; the cost argument is the same one, that 10^4 bulk nodes cost 10^4
ODE steps whatever they stand for.

**How this differs from growth, which is the obvious thing to mistake it for.**
They look alike — both advance a coarse node, both write `internal_energy`,
`radius` and `luminosity` — and they are opposites in the way that matters.

| | `grow` | bulk evolution |
|---|---|---|
| applies to | a node with a `morphology` somebody planted | every node, because it is physics |
| driven by | a `Program` — a developmental rule or a construction plan | a law with no parameters to choose |
| state | its own, in `Morphology`: segments, extent, stored energy. Not derivable from the bulk tuple | none beyond the bulk tuple itself |
| under coarsen/refine | persists; it *is* the state | must be idempotent, or looking changes the rate |
| fine-solver counterpart | none — the morphology is the model at every resolution | must agree with it on the conserved quantities |

`grow` is a program a node **runs**. Bulk evolution is a law a node **cannot
escape**. That is why `grow` may own state that only it can produce, and why
bulk evolution must own none: the moment it has private state, a node that was
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
the fine solver**: a node cooled as bulk for a century and then materialised must
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
