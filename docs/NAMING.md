# Naming pass

A survey of names that are wrong, colliding, or needlessly opaque, with what
they should become and what each will cost. Nothing here has been changed yet:
this is the thing to argue with before any of it is done.

The standard being applied is the one the rest of the project already holds
itself to. A name should say what the thing *is*, in words someone who has not
read the file would recognise, and no two things should share one.

## What a rename cannot break

Worth stating first, because it decides how nervous to be.

- **Saved worlds are safe.** Every enum crosses the wire as a *position*, not a
  name: `put_tier` writes `t.index()`, `put_property` writes
  `PROPERTIES.iter().position(..)`. Renaming a variant changes no byte. Only
  *reordering* one would, and nothing here reorders anything.
- **The database is safe.** No SQL column is named after any type below.
- **What people see is separate.** `view::TIER_NAMES`, `Species::name`,
  `Program::name` and friends are string tables the physics never reads. They
  are not affected by a type rename and can be changed independently — or not
  at all.

So this is a compile-time-checked, mechanical change throughout. The risk is
churn, not breakage.

---

## Tier 1 — collisions

Two or three different things wearing one name. These are the ones worth doing
whatever else is decided, and most of them are recent: they arrived with the
chemistry work and have not had time to set.

| Name | Where | What it actually is | Becomes |
|---|---|---|---|
| `Bond` | `chem::arrange` | A chemical bond: two atoms, an order | **`Bond`** — keep, this is the true one |
| `Bond` | `solvers::md` | A Morse potential between two particles | **`MorseBond`** |
| `Bond` | `topology` | *"One joint between two parts"* — its own doc | **`Joint`** |
| `Element` | `chem::elements` | A chemical element | **`Element`** — keep |
| `Element` | `solvers::frame` | *"One member of a frame"* — its own doc | **`Member`** |
| `Frame` | `coords` | Where a node is, how fast, which way | **`Motion`** |
| `Frame` | `solvers::frame` | A truss ready to analyse | **`Framework`** |
| `Snapshot` | `persist` | A saved world | **`Snapshot`** — keep, idiomatic |
| `Snapshot` | `causal` | One node's state at one past instant | **`Moment`** |

Two notes on the harder ones.

**`Frame` has three meanings, not two.** Beyond the two types there is the
render frame — `step_frame`, `frame_dt`, `FrameBudget`, `stats.frames` — which
is the most colloquial use and should keep the word. So `coords::Frame` gives
it up. `Motion` is what that struct holds: offset, velocity, orientation, spin
rate, proper time. `n.motion.advance(dt)` reads better than `n.frame.advance(dt)`
and stops the reader wondering which kind of frame is being advanced.

**`solvers::frame::Frame.nodes` is a fourth collision** — it is a `Vec<Vec3>` of
joint positions, nothing to do with tree nodes. Its own doc calls them joints,
so: `.nodes` → `.joints`, and the module `solvers/frame.rs` → `solvers/truss.rs`.

Also colliding, as functions rather than types:

| Name | Where | What it does | Becomes |
|---|---|---|---|
| `settle` | `chem::react` | Run to chemical equilibrium | **`equilibrate`** |
| `settle` | `persist::Snapshot` | Carry every node up to the world instant | **`catch_up`** |
| `settle` | `engine::World` | Clear shaking and falling state | **`stop_dynamics`** |

All three mean different things and none of them means "settle". The third is
the worst: it sounds like dust settling and it means *stop integrating*.

---

## Tier 2 — the name says something the thing is not

### `Species` — the headline, and the one that needs a decision

204 references. It means "one of eight lumped element buckets": H, He, C, N, O,
Si, Fe, and everything-else. It is not a species in any sense a biologist or a
chemist would recognise, and since the chemistry work it sits beside
`chem::Element`, which *is* an element — so the codebase now has two things that
both mean "what matter is made of" and neither name says which is which.

Three ways to go, and this is the one I would like a view on rather than a
default:

1. **`CoarseElement`** *(recommended)*. Honest and unmistakable next to
   `chem::Element`; the relationship between the two accounts is legible from
   the names alone. `NSPECIES` → `COARSE_ELEMENTS`. Slightly clunky to type,
   which matters less than it reads.
2. **`ElementGroup`**. Smoother, and true — each bucket *is* a group of
   elements. Very slightly less clear that it is the coarse one.
3. **`Nuclide`**. The most technically apt: the buckets exist for nuclear
   purposes (`binding_per_nucleon_mev`, burning rates) and each carries a
   representative mass number — H=1, He=4, C=12, Fe=56 — so each really does
   behave as one nuclide standing in for many. Against it: `Other` at A=65 is
   not a nuclide, and it is the least friendly of the three, which cuts against
   the point of the pass.

`Composition` stays whatever is chosen — it is a composition.

### The rest of Tier 2

| Name | Where | Why it is wrong | Becomes |
|---|---|---|---|
| `prolong` | `prolong` | Multigrid jargon; and the tree already calls this axis `refine`/`coarsen`, so there are two vocabularies for one idea | **`sample`**, module → `sampler.rs` |
| `restrict` | `state` | Same, other direction. Nothing is being restricted | **`summarise`** |
| `ProlongSpec` / `ProlongReport` | | follow the above | **`SampleSpec`** / **`SampleReport`** |
| `Transaction` | `morph` | Reads as a database transaction; it is one growth step's energy and mass accounting | **`GrowthStep`** |
| `Located` | `coords` | A `Vec3` carrying an error bound. "Located" says nothing about the error, which is the entire point of the type | **`Bounded`** |

On `prolong`/`restrict`: the docs already describe prolongation as *"a
maximum-entropy sample of the same conserved tuple"*, so `sample` is not a
softening — it is what the function is. `Tree::refine` and `Tree::coarsen` keep
their names and stay the tree-level verbs; `sample` and `summarise` become the
state-level ones they call. Two layers, two verbs each, no synonyms.

---

## Tier 3 — accurate but opaque

Worth doing if the pass is going that far; harmless to leave.

| Name | Where | Note |
|---|---|---|
| `Aggregate` | `state` | Its own doc opens *"Bulk state. Always present — this is what a node **is**."* → **`Bulk`** pairs with `summarise` and is shorter than what it replaces. 100 references. |
| `Residency` | `tree` | Why a node is holding detail. Not wrong, just abstract. → `HoldReason`, or leave. |
| `Dof` | `solvers::structure` | Degrees of freedom. Standard in the field; an abbreviation everywhere else. Leave. |

---

## Tier 4 — deliberately left alone

Listed so the pass does not churn them on a second look.

- **`Tier`** (185 refs) — short, unambiguous, and `Scale` would collide with the
  existing `Scales`.
- **`Body`, `Node`, `Tree`, `World`, `Mixture`, `Substance`, `Phase`, `Pool`** —
  say what they are.
- **`Speck`** — a drawable body in a scene. Evocative and unique; nothing else
  in the codebase could be mistaken for it.
- **`Ledger`, `Fact`, `Mailbox`, `Influence`, `Epoch`** — established, and each
  is the ordinary word for what it holds.
- **`pace`, `cadence`, `lateness`, `mixing_time`, `characteristic_time`** — the
  scheduling vocabulary is the clearest part of the codebase and should not be
  touched.
- **`encode`/`decode`, `step`, `analyse`** repeated across modules — module-scoped
  and idiomatic; not collisions in practice.

---

## What each step costs

References counted across `src/` and `tests/`, so an upper bound on the lines a
rename touches.

| Step | References | Notes |
|---|---|---|
| `Bond` ×3 → `Bond` / `MorseBond` / `Joint` | 62 | Split three ways; the compiler finds every one |
| `Element` ×2 → `Element` / `Member` | 66 | |
| `Frame` ×2 → `Motion` / `Framework` | 176 | Largest of the collisions; `Motion` is currently unused as a name |
| `Snapshot` ×2 → `Snapshot` / `Moment` | 39 | |
| `settle` ×3 | 3 definitions | Trivial |
| `Species` → chosen name | 204 | Largest single rename |
| `prolong` / `restrict` | 80 | Plus a module file rename |
| `Transaction` → `GrowthStep` | 27 | |
| `Located` → `Bounded` | 18 | |
| `Aggregate` → `Bulk` | 100 | Tier 3, optional |

## How to run it

One rename per commit, each mechanical and each with the suite green, so a
bisect lands on a single name if anything goes wrong. Order matters only in
that the collisions should go first — they are the ones that make the code
actively misleading today.

1. Tier 1 collisions (nine types, one field, one module file, three functions).
2. `Species` → whichever of the three is chosen. Largest single change at 204
   references; entirely mechanical.
3. `prolong`/`restrict` → `sample`/`summarise`, with the module rename.
4. Tier 2 remainder, then Tier 3 if wanted.

Each step also touches `docs/`, which describes several of these by name, and
the doc comments — which is most of the value: a comment saying "restrict the
bodies" reads differently once the function is called `summarise`.
