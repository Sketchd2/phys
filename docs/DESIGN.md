# Design

## 1. The arithmetic that forces everything else

| | |
|---|---|
| Baryons in a 10⁹-star galaxy | 1.2 × 10⁶⁶ |
| Bodies an RTX 2060 can step per 50 ms frame | 10⁵ – 10⁷ |
| Ratio | ~10⁶⁰ |

Every design decision below follows from that ratio. It is worth being precise
about how hopeless it is: if you had one atom of memory per particle you would
need more atoms than exist in 10³⁷ galaxies. This is not an engineering problem
with a clever solution. It is a category error to attempt at all.

The way out is to change the goal. We do not need a galaxy; we need something
no observer inside it can distinguish from a galaxy. And "no observer can
distinguish" is a *finite* requirement, because observers are finite: they have
apertures, integration times, photon budgets, and — at the bottom — the
uncertainty principle, which puts a hard ceiling on how much detail any region
can be asked to have (`solvers::quantum::max_distinguishable_states`).

So the engine's real specification is:

> For any observation an observer can physically perform, return the answer a
> real galaxy would have returned, within the observation's own error bars, in
> bounded time.

That is achievable, and everything below is about achieving it.

## 2. Existing approaches, and why each is not enough alone

The engine is built out of established techniques. What is new is the way they
are combined and the guarantee that combination is made to satisfy.

| Technique | Field | What it gives | Why it is not enough |
|---|---|---|---|
| Adaptive mesh refinement | CFD, cosmology | Resolution where gradients are steep | Refines on *physics*, unbounded when everything is interesting; no notion of an observer |
| Barnes-Hut / FMM | N-body | O(n log n) gravity | Still linear in *stored* particles; 10⁶⁶ of them |
| Sub-grid models | Galaxy simulation | Unresolved physics as effective terms | One-way: you can never zoom in and see the real thing |
| Procedural generation | Games | Infinite worlds from a seed | No conservation, no dynamics; regenerating a region loses whatever happened there |
| L-systems, growth models | Graphics, biology | Plausible organic structure | Not conservative, not thermodynamic; growth is free |
| Level of detail / impostors | Rendering | Cheap distant geometry | Visual only; the physics is not consistent with what you see |
| Conservative PDES | Distributed simulation | Correct parallel event ordering | Needs a lookahead, which is hard to obtain and usually tiny |
| Multiscale / QM-MM coupling | Chemistry | Different physics in different regions | Fixed, hand-placed regions; does not move with an observer |
| Renormalisation group | Physics | Principled coarse-graining, running couplings | A theory, not a runtime |
| Density-matrix / ensemble methods | Quantum | Statistical description of unobserved DOF | Applies only at the bottom of the ladder |

The gaps are consistent: the LOD techniques are not conservative, the
conservative techniques are not lazy, and none of them is driven by what an
observer can actually see.

## 3. What this engine adds

### 3.1 Constraint-projected materialisation

*The problem.* Procedural generation and physical simulation are normally
incompatible. Generation invents; simulation conserves. If you throw away a
cloud's million parcels and regenerate them later, the regenerated set has
different total energy and momentum, and the difference accumulates every time
the camera moves.

*The solution.* Sample from the maximum-entropy distribution consistent with the
coarse state — that makes it look right — then project exactly onto the
constraint surface. The projection is built so that each correction lives in the
null space of the constraints already satisfied:

```
centre positions        =>  Σ m r = 0
subtract mean velocity  =>  Σ m v = 0
subtract rigid rotation =>  Σ m r × v = 0     (leaves Σ m v = 0)
scale the residual by s =>  both still 0      (scaling is linear)
add rigid rotation ω    =>  L = L_target      (adds no net momentum)
add bulk drift v_b      =>  P = P_target      (adds no L about the com)
```

Nothing is ever undone, so this is a single pass rather than a fixed point. And
because the residual field has zero angular momentum by construction, it is
*energetically orthogonal* to the rigid rotation — `Σ m δ·(ω×r) = ω·Σ m r×δ = 0`
— which is what lets the energy scale be solved in closed form.

Two details make it exact rather than merely good:

- **Energy is closed algebraically, not by iteration.** `restrict` computes the
  parent's internal energy as `E_kin + Σ U_i − K_bulk`, so `Σ U_i` is a free
  parameter that appears linearly and disturbs neither momentum nor angular
  momentum. Solving for it makes total energy exact for *any* configuration,
  including the awkward ones (extreme mass ratios, degenerate geometry) where
  the velocity projection converges poorly.
- **Angular momentum that the geometry cannot carry becomes intrinsic spin.**
  Two point masses have no moment of inertia about their own axis; the inertia
  solve correctly reports this as singular. The angular momentum still has to go
  somewhere, and the honest place is the children's intrinsic spin — which is
  exactly what the parent's spin *was*, one level down.

*The guarantee.* `restrict(prolong(s)) = s` on the conserved set, to machine
precision, at every tier. Measured worst case: 5.8 × 10⁻¹⁶
(`tests/consistency.rs`).

### 3.2 Causality as the source of scheduling lookahead

Conservative parallel discrete-event simulation needs a lookahead: a guarantee
that no message will arrive with a timestamp earlier than `now + L`. Obtaining a
useful one is the classical difficulty of the field, and without it conservative
schemes deadlock or crawl.

Here it is free and it is physical. Nothing influences anything else faster than
light, so a region whose nearest active neighbour is at distance *d* has a
lookahead of *d/c*, always, with no analysis and no annotation. Better, the
lookahead is largest exactly where it is most useful: cosmological distances
give millennia of independence, and the lockstep requirement only bites at
femtometre separations where the systems are tiny anyway.

The corollary that keeps it correct: **the constraint applies between disjoint
regions, not between a node and its own ancestor.** A nucleus is not a separate
system sitting zero metres from the galaxy containing it — it is *part of* it,
and the galaxy's aggregate already accounts for it. Applying the light-speed
rule across that relationship would force the galaxy to advance at the nucleus's
zeptosecond timestep, which is precisely the catastrophe the scheme exists to
avoid (`Tree::sibling_separations`).

The same principle gates materialisation: a region only needs fine detail if
fine detail there could reach an observer within the horizon. Everything else
stays coarse no matter how interesting it is, because its influence has not
arrived yet.

### 3.3 Developmental state for matter that is ordered

Everything in §3.1 rests on the generated detail being *ergodic*: one
max-entropy sample of a gas cloud is as good as another, because no observation
can tell them apart. That interchangeability is what licenses throwing detail
away.

A tree is not like that. Its branch structure is low-entropy and historically
contingent — not a typical sample from anything, but the specific record of
which branch got shaded in year three — and the difference is *observable*,
because someone who saw it yesterday will notice if handed a different one. The
conserved tuple pins mass and momentum exactly and pins none of what matters.

So for structured matter the generator changes, and only the generator:

```text
    structure = program(genome, age, events)
```

a pure function, addressed exactly as everything else in `rng.rs` is, so it
regenerates bit-for-bit. Measured: 104 bytes of developmental state standing in
for 6,000 rendered parts — **10,615× smaller** — and the compression grows with
the structure, because the state does not.

Three consequences make this fit rather than bolt on:

**Growth runs on the aggregate.** A forest does not grow by integrating 10⁹
trees; it grows by advancing one ODE on a forest node, at O(1) per node. That is
cheap enough to run on the entire world every frame while the fine structure
stays unbuilt — so the laziness the rest of the engine works for is simply free
here. The demo grows a tree for 160 simulated years across 1,920 growth steps
without materialising it once.

**Growth is a transaction, not an exemption.** Building order out of disorder
costs free energy and exports entropy, so every step returns a `Transaction`
that is validated before it is applied: energy in equals energy stored plus heat
released plus radiation, and local entropy plus exported entropy is
non-negative. A program that tried to grow too efficiently is *refused*, not
smoothed over.

**The same projection applies.** A generated oak goes through the identical
`close_books` routine a sampled gas cloud does — the momentum, angular-momentum
and energy projection, the algebraic energy close, the intrinsic-spin residual.
Adding morphology needed no second conservation story, only a second way of
choosing where the parts go. Measured worst case across four programs, three
masses and three budgets: 3.6 × 10⁻¹⁶.

Two things the implementation had to get right that are easy to get wrong.
`restrict` is not entitled to an opinion about a structure's entropy: it sees an
unstructured heap of parts and reports the entropy of the same mass as a gas,
erasing precisely the order that makes the thing a structure. `Body` carries no
topology, so the information is not there to recover — the developmental state
is the authority. And a node holding a structure is not made *only* of that
structure; assuming so works until a limb is severed, at which point the
structure's mass drops while the node's does not, and the missing mass is
silently redistributed into the surviving branches — a tree that grows heavier
every time you prune it.

### 3.4 Observation as commitment

An unobserved quantity has no value. A measured one is recorded in a ledger and
returned identically forever after (`observe::Ledger::get_or_sample`). This is
what makes the deception undetectable *over time* rather than merely at one
instant: measure a decay time twice and you get the same answer, because the
first measurement made it a fact.

At the subatomic tier this stops being a trick. An unobserved system genuinely
does not have a definite position; a measurement genuinely does sample from a
distribution and leave the system in the sampled state. The engine's
"generate on demand, commit on observation" machinery is, at the bottom of the
ladder, a literal implementation of the physics rather than an approximation of
it. The approximation runs the other way — classical trajectories are the
approximation, and `solvers::quantum::regime` says when they are safe.

The ledger is the only structure in the engine whose size grows with what users
*do* rather than with the size of the universe. In the demo it holds one fact in
60 bytes against 20 MB of regenerable detail.

### 3.5 The frame budget as the invariant

A conventional simulation decides what to compute and takes however long it
takes. This one is given 50 ms and decides what fits. Every candidate piece of
work carries an estimated cost and an estimated value:

```
value = salience × urgency × error × (1 + novelty)
```

where salience is the solid angle the work subtends for some observer, urgency
is how close it is to the causal horizon, and error is the estimated
inaccuracy of *not* doing it (an unresolved Jeans length, a dynamical time
shorter than the frame step). A greedy knapsack by value density fills the
frame; the rest is reported as detail debt.

Greedy rather than exact on purpose: the 0/1 knapsack is NP-hard, greedy is
within a factor of two, and that is far inside the error of the cost estimates
themselves. Spending frame time to plan the frame better than the estimates
justify would be self-defeating.

The cost model calibrates itself from measured frames, asymmetrically — fast to
back off, slow to push — because being late is visible to the user and being
early is not.

## 4. Structure

### 4.1 The scale tree

A node is a region at a tier holding a bulk `Aggregate`. It has two independent
finer representations:

- **materialised bodies**, produced by `prolong` — cheap to make, cheap to
  destroy, regenerable bit-for-bit;
- **promoted children**, full nodes standing in for individual bodies, created
  only for the few bodies someone is actually looking at.

Materialising takes you from "a molecular cloud" to "a million gas parcels";
promoting takes you from "one of those parcels" to "a protostar with its own
interior". A galaxy-to-nucleus path promotes about twenty times and materialises
about twenty times, so the live node count along any single zoom is in the
dozens and the body count in the hundreds of thousands — not 10⁶⁶.

Crucially, **the tier is derived from physical size, not from tree depth**
(`Tier::containing`). A tier is a physics regime; many refinements happen within
one. Conflating the two is a tempting simplification that gives seven levels of
refinement between a galaxy and a nucleus when the mass ratio alone demands
more than twenty.

### 4.2 No global coordinates

The engine spans 10²¹ m to 10⁻¹⁵ m. Expressing a femtometre offset at galactic
radius needs ~36 significant digits; `f64` has 15.95, and fixed-point would need
120 bits per axis — 45 bytes of position per particle, three times the rest of
the state and unusable on a GPU.

So no global coordinate is ever formed. Every position is an offset from its
parent, and two positions are compared by walking to their lowest common
ancestor (`Tree::separation`). The precision you get is then exactly the
precision the question deserves: two nucleons in one nucleus are located
relative to each other to ~10⁻³¹ m; a nucleon and a star across the galaxy to
~10⁵ m. That is not a limitation, it is the same statement as the causal gate —
things that are far apart cannot interact sharply, so they need not be located
sharply relative to one another. `coords::Located` carries the accumulated
round-off so the engine can refuse to evaluate a process at a separation it
cannot resolve, rather than reporting a number finer than its own arithmetic
supports.

### 4.3 Determinism contract

Regeneration is a pure function of an address:

```
value = f(world_seed, path_key, epoch, purpose, index)
```

- **path_key** is a 128-bit rolling hash of the child-index path from the root,
  so it survives the node being destroyed and rebuilt (arena indices do not).
- **epoch** increments only when a recorded interaction changes a node's
  contents, so an undisturbed region regenerates identically forever.
- **purpose** decorrelates streams, so adding a new physical process never
  disturbs the numbers an existing one draws — old scenarios keep replaying.
- **index** makes draws order-independent, so a GPU kernel computing draw *i* in
  any lane order matches the CPU exactly.

Two consequences the implementation had to be careful about. First, anything
drawn repeatedly *over time* — a thermostat, a scattering kernel — must include
the step in its address, or "deterministic" quietly becomes "identical every
step", which for a Langevin thermostat means a constant force and unbounded
heating (`rng::Stream::split`). Second, reductions are pairwise with a fixed
tree shape (`math::det_sum`) so that CPU and GPU agree bit for bit.

Coarsening is also made *idempotent*: if the restricted state agrees with the
stored aggregate to within 10⁻¹², the coarse state is left exactly as it was.
Without that, every visit perturbs the aggregate in its last bits, the next
materialisation samples from a marginally different distribution, and a region a
user visits a thousand times slowly drifts away from itself.

## 5. What the user can do

Every interaction is expressed as a change to a conserved quantity delivered at
a place and a time, so none of them can violate the invariants the rest of the
engine depends on.

| Action | Behaviour |
|---|---|
| **Observe** | A demand for resolution over a solid angle. Decides what gets materialised. Returns retarded, Doppler-shifted, shot-noise-limited readings. |
| **Measure** | An interaction, because measurement disturbs. Locating a proton to 1 fm deposits 5.2 MeV, and the engine applies it. |
| **Impulse / Deposit / Extract** | Momentum and energy, delivered after `d/c`. Pins the target. |
| **Inject** | Adds matter with a composition; rebalances baryon and lepton number. |
| **Pin** | Forces detail to persist even when nobody is looking. |
| **Author** | Sets a bulk property directly. The one path that can break conservation — so it records exactly how much it broke it by, in an audit log. |
| **Time control** | `time_rate` scales simulated seconds per wall second. Because the timestep is fixed by accuracy, zooming into a nucleus does not slow the frame rate — it slows *time*. |

## 6. Honest limitations

- **The GPU path is designed, not written.** The cost model carries an explicit
  60× assumption for it, stated as an assumption. Everything measured in
  `docs/PERFORMANCE.md` is single-core CPU.
- **The CPU tree code is 3–10× off production quality.** It uses single-body
  leaves rather than buckets and per-particle rather than grouped traversal;
  both are known, standard optimisations, quantified in `docs/PERFORMANCE.md`.
- **Nuclear and electronic structure are rate- and table-based**, not first
  principles. This is the same choice every stellar evolution code makes, and it
  is not a compromise — the tabulated rates *are* the experimental facts.
- **General relativity is post-Newtonian only.** Adequate to a few hundred
  Schwarzschild radii; a genuine metric solver would be needed closer.
- **Chemistry is eight lumped species.** Enough for burning, cooling, opacity
  and gross chemistry; not enough for real molecular diversity.
- **Turbulence is a prescribed solenoidal field**, not a solved cascade. It
  produces the right correlations at one scale, not the right spectrum across
  scales.
- **Structures have no topology.** `Body` carries no bonds, so materialising a
  tree and running molecular dynamics on it would let the trunk fall apart.
  Emitting constraints alongside bodies is the next substantial piece of work,
  and until it exists a structure is geometry that happens to hold still.
- **Structures cannot straddle nodes.** A building spanning several nodes, or
  roots reaching into the soil node, would need cross-links, which the strictly
  hierarchical tree deliberately forbids — that hierarchy is what makes the
  precision and causal-gating arguments work. For now a structure must fit
  inside one node.
- **Nothing builds the buildings.** Planned construction advances at a supplied
  labour rate; there are no agents with goals, plans or logistics behind it.
- **The conservation guarantee is about the conserved tuple, not about
  trajectories.** Regenerated detail is a valid sample from the right
  distribution, not the specific configuration that would have evolved had the
  detail been tracked continuously. For unobserved matter this is exactly the
  right standard — it is the standard statistical mechanics uses — but it is a
  weaker claim than "the same simulation, computed lazily", and it should not be
  advertised as the stronger one.
