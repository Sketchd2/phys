# Backlog

Deferred work, with the reason it was deferred. Not a roadmap — the roadmap is
the architecture review. This is the smaller list of things noticed while
building something else, parked deliberately rather than forgotten.

An entry earns its place by being **specific about the trigger**: what would
make it worth doing, or what will go wrong if it is not. "Would be nice" is not
a trigger.

---

## Pace control has no honest "manual" mode

**Noticed:** building `phys-persist` (Phase 1).
**Where:** `engine.rs` — `World::refresh_pace`, `World::pace`, `World::paced_to`.

`refresh_pace()` runs at the top of every frame and recomputes `pace` from
`paced_to`, which is correct and is what makes "zooming in slows time" fall out
as arithmetic rather than policy. But it means an assignment to `world.pace` is
silently discarded on the next frame, and the only way to drive the pace by hand
is:

```rust
w.paced_to = NodeIdx::NONE;   // obscure: "nothing is being watched"
w.pace = 90.0 * 24.0 * 3600.0;
```

That works because `refresh_pace` returns early on a none index, but it reads as
a bug rather than an intent, and nothing stops a later refactor from removing
the early return and quietly breaking every caller doing this.

**The fix** is a small enum making the choice explicit:

```rust
pub enum PaceMode {
    /// Follow a node's characteristic time. The default, and the reason
    /// resolution and time rate stay coupled.
    Follow(NodeIdx),
    /// A fixed span per frame, set by the caller. Growth demos, tests,
    /// deterministic replay, anything scripted.
    Fixed,
}
```

**Trigger:** do it before anything else needs to script the clock — the headless
worker (Phase 2) and replay both will. Doing it after there are several callers
setting `paced_to = NONE` means finding all of them.

**Cost:** small. One enum, one branch in `refresh_pace`, and updating the two
current callers. Deliberately not done inside a persistence commit, because
changing scheduler semantics there would have muddled what that change was for.

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
is the same operation `restrict` already performs on mass and energy, applied to
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

## The pacer cannot serve a crowd

**Noticed:** building the volume query (Phase 2 bandwidth work).
**Where:** `engine.rs` — `refresh_pace`, `pace_to`, `lateness`.

`phys-headless serve` promotes two hundred neighbours around the node it is
watching and paces the clock to that node. Worst lateness settles around 300:
two hundred live nodes cannot all be re-solved at the cadence of the finest one
inside a 50 ms frame, so most of them fall behind by a couple of orders of
magnitude and stay there. The scheduler is behaving correctly — it ranks by
lateness and spends what it has — but "correctly" here means every node is
equally, badly late.

It shows up in `tests/view.rs::a_neighbourhood`, which has to call `pace_to`
explicitly or the neighbours travel 10^8 node radii between two frames.

**The fix** is probably not more throughput. It is that a crowd of nodes at the
same tier doing nothing in particular should not each be a scheduling unit:
they want to be *one* unit — a coarse node with them as bodies — until
something happens to one of them, which is the promotion rule running in
reverse. `coarsen` already exists and `mixing_time` already says when detail may
be released; what is missing is the same judgement applied to a *sibling group*
rather than to one node's interior.

**Trigger:** the first scene that needs more than a few dozen live nodes at
once, which is any city. Not urgent while the play space is one node deep, and
deliberately not attempted inside a bandwidth change — the two would have been
impossible to tell apart in the numbers.

---

## Recipes assume the client's sampler matches the server's

**Noticed:** building `view::Recipe` (Phase 2 bandwidth work).
**Where:** `view.rs` — `Recipe::build`, and `prolong.rs` underneath it.

A recipe is ~300 bytes that regenerate a node's detail on the client, and it is
only ever sent for detail the server does not itself hold. That is what makes it
safe: there is no server-side truth for a client's version to disagree with, so
two clients sampling untouched scenery a few ulps apart is a difference nobody
can observe. `build` also checks the generated mass against the aggregate, which
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
