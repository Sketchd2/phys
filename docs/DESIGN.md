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
| Finite element analysis | Engineering | Exact internal forces | A sparse solve per structure per load case; far too slow to run on everything |
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

### 3.4 Cohesion, and why the load path being a tree matters

A generated structure was, until topology existed, geometry that happened to
hold still: nothing in a `Body` says this segment is attached to that one, so
materialising a tree and handing the parts to the molecular dynamics solver
would let the trunk fall apart.

The joints come from the same program as the geometry — a branching generator
knows perfectly well which segment grew out of which, and was simply discarding
it — so cohesion is regenerated rather than stored, exactly like shape.

The load path of every structure the engine generates is a **tree**: a branch
hangs off a branch, a floor stands on the floor below, a course of bricks rests
on the course beneath. That is worth a great deal. A real framed building is a
redundant lattice whose true internal forces need a stiffness matrix and a
sparse solve — thousands of unknowns, iterative, awkward to budget for at 20
frames a second. On a tree it is exact in one pass: accumulate force and moment
from the leaves inward, and every joint's load is known in O(n) with no solve at
all. Since parts are emitted parents-first, a single reverse iteration over the
array does it.

A redundant structure gets a real solve instead, and the fast path is not an
approximation of it — it is the same code with an empty tie list.

#### Beam elements, and why springs were not enough

The redundant path was a network of anisotropic springs: stiff along each
member's axis, soft across it, with the transverse stiffness set to `3EI/L³`.
That is the tip stiffness of a cantilever, so the obvious validation — a
cantilever's `PL³/3EI` deflection — passes and proves nothing.

It is wrong for everything else, because a spring pair carries no *moment*
between its ends. A member's rotation is invisible to its neighbours, so the
model cannot tell a fixed end from a pinned one. The discriminating case is a
beam built in at both ends, whose midspan deflection is `PL³/192EI` — sixty-four
times stiffer than the cantilever, entirely because the fixed ends resist
rotation. A translation-only model gets that wrong by that factor, and braced,
portalised and continuous structures are exactly the ones that need the
redundant solver in the first place.

`solvers::frame` is a standard three-dimensional Euler-Bernoulli formulation:
six degrees of freedom per node, twelve per element, matrix-free and
Jacobi-preconditioned. Preconditioning is not a nicety here — translational
stiffnesses of order `EA/L` and rotational ones of order `EI/L` differ by
`L²/r²`, four orders of magnitude for a slender member — and the claim is
measurable rather than asserted: the same problem takes 242 iterations
preconditioned and does not converge at all without.

Two failure modes ride on top of the elastic solve. **Euler buckling**, because
a slender member in compression fails far below its material strength and no
stress check will ever notice — it is a stability failure, not a strength one,
and omitting it means a hundred-metre column reporting 57% utilised while it
folds up. And **elastic-perfectly-plastic redistribution** by secant stiffness,
because a ductile material sheds load off an overstressed member onto its
neighbours and a brittle one does not. That single number, `Material::ductility`,
is the difference between a steel frame that sags, redistributes and warns you
and a masonry wall that is standing one moment and rubble the next.

#### The same tree, twice

There is one more consequence of the load path being a tree, and it is what
makes the dynamics affordable at all.

A tree has a **perfect elimination ordering**. Eliminate the leaves first and
each one's only surviving coupling is to its parent, so factorising the
stiffness matrix produces no fill-in whatever and costs one 6×6 inverse per
joint. For a determinate structure that factorisation *is* the inverse, and
preconditioned conjugate gradient converges in a single iteration — at any
timestep and any size. For a braced one the braces are the only edges outside
the forest, so it stays a good approximation and the iteration count reflects
how redundant the structure is rather than how big it is.

This matters because Jacobi preconditioning is hopeless on a tree. A slender
chain of `n` beam elements has a condition number growing like `n⁴`, so a
2000-member tree took 3645 iterations and 692 ms for one dynamic substep. The
same substep with the factorisation takes one iteration and 5 ms.

Two details decide whether it works. The forest has to be a *maximum* spanning
forest by stiffness rather than a breadth-first one, or a soft stay skipping
along a chain gets kept and a load-carrying member dropped. And the
factorisation is declined when the structure is redundant enough that it stops
paying for itself, because applying it costs about an order of magnitude more
than a diagonal division. Both are measured rather than assumed, in
`docs/PERFORMANCE.md`.

#### Motion

Knowing what a bond can take is not knowing what it makes its neighbours *do*.
Until `solvers::dynamics`, a tree in a gale either stood exactly still or
snapped, with nothing in between.

The dynamic step runs on the *same operator*. `Frame` grew a lumped mass vector
and two scalar coefficients, turning its operator from `K` into `s_m M + s_k K`;
Newmark-beta with Rayleigh damping is then a static solve against that shifted
operator, so the conjugate gradient, the preconditioner, the element forces and
the failure criteria are all shared. A structure cannot move according to one
stiffness and break according to another, because there is only one.

Not backward Euler, which is the `γ = 1, β = 1` corner of the same family: it is
unconditionally stable but removes 93% of a structure's energy in four cycles
for reasons that have nothing to do with the material. A tree under it deflects
into the wind and stops dead. Trapezoidal Newmark conserves energy exactly for a
linear system, so what damping there is, is the damping that was asked for.

The payoff is a number a quasi-static analysis cannot produce. A load that
arrives suddenly deflects a structure about *twice* as far as the same load
standing still — the dynamic load factor for a step load is exactly 2, and the
solver measures 1.985. A gust that a static check passes at 60% utilisation
breaks the same member outright.

What this buys is that damage stops being scripted. Nothing in the engine says
"lightning destroys a tree" or "wet snow breaks branches". Snow settles on
upward-facing area, wind adds drag to projected area, lightning deposits
enthalpy along the conduction path, fire raises temperature and consumes
material — and then the *same* stress calculation decides what survives. A limb
comes down because the moment at its base exceeded what its cross-section could
carry, which is also why real limbs come down.

Three calibrations had to be right for that to be true rather than merely
plausible, and each was wrong at first in a way only a physical check caught:

* **Member radii follow from density, not from the position scale.** Scaling
  radii geometrically alongside positions produced a 13 m tree with a half-metre
  trunk radius — three times the whole tree's volume in the trunk alone. Section
  modulus goes as `r^3`, so a factor of two in radius is a factor of eight in
  apparent strength, and every structural conclusion drawn from it was wrong.
* **Snow load is bounded by the crown's silhouette, not the sum of member
  areas.** Branches shade one another and snow falling between them reaches the
  ground; summing per-member areas over-counts by the crown's area index and
  made 100 mm of snow destroy a tree that would not have noticed it.
* **Interception capacity depends on wetness.** Dry powder barely adheres and
  blows off; wet snow near freezing bonds to the bark, and that is the snow that
  brings limbs down. There is no single value that is right for both — treating
  all snow alike either makes powder lethal or makes wet snow harmless.

#### Generalising it

The first version of this had the weather *in* the solver: an `Insult` enum with
arms for snow, wind, lightning and fire. That is the wrong shape. Every new load
case was a new arm inside the physics, and the four cases that happened to be
implemented were the four that existed.

What the solver knows now are *mechanisms* — a body-force field, drag in a
moving fluid, mass accreting on upward-facing surfaces, energy conducted along a
path, a thermal field, a point impulse. Snow, wind, lightning and bushfire are
constructors in `solvers::structure::weather` that produce mechanisms, and a
user can write their own without touching the solver. The same accretion law
covers snow, rime and volcanic ash; the same drag law covers air and water at
eight hundred times the density. Materials went the same way: a closed enum of
four became a struct of numbers anyone can construct.

The determinacy question got the same treatment. A forest load path is exactly
solvable by statics; anything with a brace in it is not, and how the load
divides depends on relative stiffness. So the structure decides which solver
runs — and the fast path is not an approximation of the general one, it is the
general one with an empty tie list. `tests/topology.rs` checks that they agree
to 10^-6 on a determinate structure, and validates the redundant path against
the analytic three-bar truss: 585.8 N in the middle bar against an analytic
585.8 N, in three conjugate-gradient iterations.

Two failures worth recording, because in both the solver was right and the
expectation was wrong. Bracing a tower does *not* lower the peak utilisation
over all members — load taken off the columns goes into the braces, and the peak
may legitimately move there. It does not necessarily relieve any *particular*
column either, because under a lateral load bracing transfers force between the
windward and leeward sides. What it does is reduce the total stress carried,
which is what the test asserts now.

### 3.5 Observation as commitment

An unmeasured quantity has no value. A measured one is recorded in a ledger and
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

### 3.6 One instant, and what it costs to stay on it

The world has a single instant and everything in it is at that instant. This is
not obviously affordable — a nucleus needs 10⁻²³ s steps and a galaxy arm is
happy with a million years — and the naive reading of "simulate every scale at
once" is that the fast thing sets the pace for everything. That is exactly what
an earlier version of this engine did, by taking a global minimum over every
live node's stability limit, and it meant that resolving anything small stopped
everything large.

The way out is that **being at an instant and being solved at it are different
things**. A node's offset under constant velocity and its orientation under
constant angular velocity are closed-form solutions, so between two solves a
node can be asked where it is at any moment and answer without approximation.
Carrying a node forward therefore costs one add and introduces no error at all.
The expensive question is not where a node *is*, it is what it is *doing* — and
that only needs re-deriving when it has changed.

So each frame:

1. every node's **lateness** is computed against the new instant — how long
   since it was last re-derived, divided by its own characteristic time;
2. the budget takes the most overdue work that fits;
3. everything else is **coasted** to the same instant in closed form.

Step 3 is nearly every node, nearly every frame.

**The characteristic time.** τ = ℓ / v, where ℓ is the resolution the node is
currently drawn at and v the fastest thing it is doing: moving through its
parent, turning, or rearranging inside itself. One expression, and it spans the
ladder without a table of special cases — an isolated Earth at one radian of
turn per update, the same Earth in orbit at 3.6 minutes, a swimming bacterium
at half a second, a nucleus at 10⁻²³ s.

Two things are deliberately in or out of v, and both were found by getting them
wrong:

- **Rotation is in.** Jupiter's equator runs at 12.3 km/s against about 1 km/s
  of internal signal speed, so a cadence chosen from the signal speed alone
  samples the planet once every 1.45 revolutions and aliases it completely.
- **The equilibrium sound speed is out.** It is the fastest speed anything
  inside a node is travelling at, and it is the wrong quantity: a body in
  thermal equilibrium is not changing, however fast its molecules are going.
  Scheduling on it asks the engine to re-solve a bacterium every five
  nanoseconds and produce the same answer every time. The internal term is the
  *stirring* speed — the internal energy the stored temperature does not
  account for. The sound speed keeps its real job, bounding the timestep of a
  node that is already being solved.

**When a span cannot be integrated.** A resolved nucleus asked to cross a
millisecond would need 10¹⁹ steps. It also does not need them: over that span it
has sampled its accessible states 10¹⁹ times, and where it ends up is a draw
from its equilibrium ensemble, not the endpoint of a trajectory. So past a
sub-step ceiling the node is *thermalised* — restricted to its bulk state,
carried across in closed form, and drawn again at the far end. Both halves are
things the engine already guarantees: restriction is conservative to within
`IDEMPOTENT_TOLERANCE`, and prolongation is a maximum-entropy sample of the same
conserved tuple, which is exactly what a fresh draw from the ensemble means.
Detail somebody has touched, and detail something finer has been built on, is
exempt and falls behind honestly instead.

### 3.7 The frame budget as the invariant

A conventional simulation decides what to compute and takes however long it
takes. This one is given 50 ms and decides what fits. Every candidate piece of
work carries an estimated cost and an estimated value:

```
value = lateness × urgency × error
```

where lateness is how many of its own characteristic times the node has gone
unsolved, urgency is how close the work is to the causal horizon, and error is
the estimated inaccuracy of *not* doing it (an unresolved Jeans length, a
dynamical time shorter than the frame step). A greedy knapsack by value density
fills the frame; the rest is reported as detail debt.

Value used to be led by observer salience — the solid angle the work subtended
— with a small novelty bonus to stop a region starving. That put the camera in
charge of the physics, which is the wrong way round: a tree falls whether or not
anyone is pointing at it. Lateness needs no anti-starvation term of its own,
because a node passed over grows more overdue every frame until it outranks
whatever kept beating it. What the camera legitimately controls is *resolution*,
and resolution feeds back into lateness on its own: a node drawn more finely has
a shorter τ and comes due more often.

Greedy rather than exact on purpose: the 0/1 knapsack is NP-hard, greedy is
within a factor of two, and that is far inside the error of the cost estimates
themselves. Spending frame time to plan the frame better than the estimates
justify would be self-defeating.

Two rules keep the invariant from degenerating into a world that is perfectly on
time and completely still:

- **A plan never accepts nothing.** Work is indivisible at this level — a solver
  pass over a node is one pass or none — so a world whose cheapest node costs
  more than a whole frame would defer every task forever. The best one is taken
  and the overrun is recorded.
- **When the work does not fit, simulated time gives way.** The signal is
  shortfall on work that actually ran: a node the frame accepted, integrated,
  and which still did not reach the instant. That says the span was too long,
  and shortening it fixes it. Work deferred outright says something else —
  more of the world is resolved than can be simulated — and slowing the clock
  does not help that, it just stops it while the debt stays where it was.

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
  only for the few bodies something is actually happening to.

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
| **Pin** | Marks detail as non-regenerable, so it is stored rather than re-drawn. |
| **Author** | Sets a bulk property directly. The one path that can break conservation — so it records exactly how much it broke it by, in an audit log. |
| **Time control** | `time_rate` scales simulated seconds per wall second, on top of a pace taken from whatever is being watched (`pace_to`). Zooming into a nucleus does not slow the frame rate — it slows *time*, and that is now arithmetic rather than policy, because materialising a node shortens its characteristic time and the pace is re-read every frame. |

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
- **The refinement table is a choice, not a derivation.** How many children a
  node has and how they are arranged comes from a table in `prolong.rs`, and the
  table is only *consistent* — a node's radius, its contents' count and their
  interaction radii have to agree, and where they do not the materialised
  configuration is unphysical. Two cases where they did not have been fixed
  (eight thousand molecules in an atom-sized node, sixty-four atoms in a
  molecule-sized one); there is no proof there are no others, only tests that
  step every scenario and check nothing heats up.
- **Structural dynamics is small-displacement and linear.** The restoring force
  is exact for members whose chord rotates by well under a tenth of a radian,
  which covers a building in a storm and a trunk in a gale but not a sapling
  bent double. `StepReport::displacement_ratio` reports the worst chord rotation
  every step, so the regime is measured rather than assumed; a corotational
  element formulation would remove the limit and has not been written.
- **Chemistry is valence-limited, not bond-order.** Bonds form and break, and
  what stops a hydrogen acquiring five neighbours is a valence count rather than
  a many-body bond order that weakens with over-coordination. A Tersoff or
  Brenner potential is the better answer where the energetics of intermediate
  coordination matter; the question asked here is only whether a molecule holds
  together and can react.
- **Falling debris is bounded and one-way.** At most twenty-four pieces fall at
  once and each is written off after twelve seconds, because a crown fire breaks
  thousands of joints and simulating every twig's descent would spend the whole
  budget on litter. Debris collides with the structure it fell from and with the
  ground, not with other debris.
- **The design pass optimises against an envelope, and an envelope is a
  choice.** A vertical overload at 2.5 g and the program's design flow from a
  few directions is not every load a structure will meet. What it buys is that
  the structure is not *brittle* off-design, which is checked; what it cannot
  buy is optimality against a load nobody listed.
- **Structures cannot straddle nodes.** A building spanning several nodes, or
  roots reaching into the soil node, would need cross-links, which the strictly
  hierarchical tree deliberately forbids — that hierarchy is what makes the
  precision and causal-gating arguments work. For now a structure must fit
  inside one node.
- **The mixing time is a discriminator, not a derivation.** Detail is released
  once a node has had time to forget what put it there, and "time to forget" is
  taken as the crossing time of the *random* part of the internal motion, with
  ordered rotation removed and with structured matter — anything carrying a
  morphology, a topology, promoted children or a user's fingerprints — exempted
  outright. That gets the cases that matter right (a gas parcel forgets in
  milliseconds, a rotating disc never, a broken tree never) and it is not a
  first-principles relaxation time. A collisionless stellar system's real
  relaxation time is a two-body calculation this does not do.
- **The sub-step ceiling is a number.** 256 sub-steps is where following a
  trajectory stops being the cheaper way to answer the question and the node is
  crossed by its ensemble instead. Nothing derives 256; it is set high enough
  that anything watchable is integrated properly and low enough that a nucleus
  resolved inside a galaxy cannot stall a frame.
- **Sub-step allowance is shared equally, not by need.** Once the plan has
  chosen which nodes run, each gets the same number of passes, because the
  alternative — spending until a wall clock runs out — would make the schedule
  depend on how fast the machine happened to be that frame, and replay would
  diverge. A node that needs more passes than its share falls behind and says so
  through its lateness.
- **Nothing builds the buildings.** Planned construction advances at a supplied
  labour rate; there are no agents with goals, plans or logistics behind it.
- **The conservation guarantee is about the conserved tuple, not about
  trajectories.** Regenerated detail is a valid sample from the right
  distribution, not the specific configuration that would have evolved had the
  detail been tracked continuously. For unobserved matter this is exactly the
  right standard — it is the standard statistical mechanics uses — but it is a
  weaker claim than "the same simulation, computed lazily", and it should not be
  advertised as the stronger one.
