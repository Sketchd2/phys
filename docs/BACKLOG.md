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
