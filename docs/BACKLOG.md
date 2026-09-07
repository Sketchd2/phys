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

## A node whose bodies are 10^20 radii outside it

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

**The real cause is upstream.** The spacing is `2h` and `h` comes from
`radius / count^(1/3) * 1.2`, so `s = 1e-9` implies a node about seven
nanometres across — and its bodies are sitting 6.6x10^10 m away, twenty orders
of magnitude outside it. Bodies are node-relative by construction and should be
within a few radii of the origin. Something is putting parent-frame coordinates
into a node's body list, or a node's radius is being set without its contents.
The failing tests both drill the full ladder ("24 tiers, galaxy to nucleus"), so
a scale transition is the place to look.

**Two fixes, and they are separate.** The grid should use saturating arithmetic
regardless — a bucket index that cannot be represented should clamp, not wrap,
whatever put the body there. And the thing that put the body there needs
finding, which is the actual bug; a `debug_assert` that a node's bodies lie
within some multiple of its own radius would have caught it at the source
rather than twenty tiers later.

**Trigger:** before trusting any SPH result at a tier boundary, and before the
debug suite can be used as a gate. It is not blocking the release suite, which
is green.
