# Physics

## 1. The scale ladder

A **tier** is a physics regime, not a tree level: many refinements happen within
one tier, and the tier changes only when the characteristic size crosses a
boundary where a different description becomes appropriate
(`units::Tier::containing`).

| Tier | Size range | Timestep | Light crossing | Representation | Solver |
|---|---|---|---|---|---|
| Galactic | > 10¹⁸ m | 100 kyr | 10.6 kyr | Collisionless super-particles + static halo | Barnes-Hut |
| Stellar | 10¹² – 10¹⁸ m | 100 yr | 1.06 yr | Clouds, clusters, individual stars | Barnes-Hut + SPH |
| Planetary | 10⁴ – 10¹² m | 100 s | 3.3 s | Stars, planets, orbits, interiors | Barnes-Hut + 1PN + SPH |
| Continuum | 10⁻⁸ – 10⁴ m | 1 ms | 3.3 ns | Fluid parcels, grains, continuum matter | SPH |
| Molecular | 3×10⁻¹⁰ – 10⁻⁸ m | 1 fs | 3.3 as | Molecules | Molecular dynamics |
| Atomic | 10⁻¹⁴ – 3×10⁻¹⁰ m | 1 as | 0.3 zs | Atoms, electronic structure | MD + level ensembles |
| Nuclear | < 10⁻¹⁴ m | 1 zs | — | Occupation numbers, not trajectories | Statistical |

The demo's descent from a disc galaxy to a nucleus crosses all seven in 23
refinement levels, spanning 10³⁰ in length and 10⁶⁰ in mass, with conservation
error never exceeding 4 × 10⁻¹⁶ at any level.

## 2. The conserved set

Six quantities are preserved exactly across every scale transition, in both
directions:

**energy · momentum · angular momentum · charge · baryon number · lepton number**

Baryon and lepton number are what make the subatomic tier consistent with the
galactic one. They are additive integers in disguise, which is why the
bookkeeping can span 60 orders of magnitude: you cannot fuse hydrogen in a star,
coarsen the star, refine it again and find the protons back.

### Energy decomposition

```
E = M c²  +  K_bulk(M, P)  +  U_internal  +  Φ_gravitational  +  Φ_cohesive  +  Φ_external
```

- `K_bulk` is the **exact** relativistic form, `√((Mc²)² + (Pc)²) − Mc²`, not
  `p²/2M`. The decomposition has to be invertible to the last bit — `summarise`
  recovers `U` from `(M, P, E, Φ)` — and a Newtonian bulk term makes the round
  trip lossy at the 10⁻⁵ level for anything moving at galactic rotation speeds,
  which shows up as visible energy drift when a user pans across a disc.
- **`Φ_gravitational` and `Φ_cohesive` are separate because expansion releases
  the first and not the second.** Gravitational binding is held by the node's
  own long-range field, so spreading the contents out is exactly what pays it
  back; cohesive binding is held by bonds at a scale far below the node, and
  moving the bodies apart releases none of it. `sampler::sample` has a
  relaxation loop whose whole premise is the first — a configuration too
  tightly bound to hold the energy it claims must be bigger — and run against
  the sum it inflated a cubic metre of granite by `1.5^32 = 4.3 × 10⁵`,
  scaling for thirty-two iterations to move a term that was not the negative
  one, then giving up and leaving the inflation in place. Measured: the block's
  contents sampled `4.502 × 10⁵` radii outside the node they are inside,
  against 3.35 for a spiral galaxy, which is what a correctly relaxed
  configuration looks like. Split, they sample at 1.04.

  `Φ_cohesive` is carried through both directions of a scale transition
  untouched, like `Φ_external`: `summarise` measures where the bodies ended up,
  and a bond is far below the scale of a body, so it has nothing to report.
- `Φ_external` is separate from both because it is *not recoverable from
  the node's own contents*. Refine a galaxy's baryons and their mutual potential
  is nine times smaller than the dark halo's grip on them. Folding the two
  together makes every refinement demand a thermal budget that does not exist,
  and forces the sampler either to invent energy or to violate the virial
  theorem. This was a real bug: it produced a systematic 1.3 × 10⁻⁸ energy error
  at the root node, identical for every profile and every particle count, which
  is exactly the signature of a modelling error rather than a numerical one.

Nuclear binding energy is deliberately **not** in this list. It is already
inside the rest mass. It enters the invariant set only through the composition,
which is conserved exactly; fusion releases energy by *changing* the composition,
and the burning solver takes the difference. Treating it as an available energy
pool — an easy mistake, since it has units of energy — injects ~10¹³ J/kg of
spurious heat.

### Measuring conservation error

Errors are measured against the **natural scale** of each quantity (the sum of
the magnitudes of its contributions), not against its own net value.

This is forced by the physics. A hot cloud with almost no net rotation has
particles each carrying angular momentum of order `m r v`, but the sum nearly
cancels — the net can be 20 orders of magnitude below any individual term.
Double precision holds 16 digits, so that net is *not representable*, and no
finite arithmetic would do better. Dividing by it reports 100% and tells you
nothing; dividing by the total angular momentum content tells you what is true
and what matters, namely that the bookkeeping is good to one part in 10¹³ of
everything in the system. No observer inside the simulation can make a
measurement finer than that.

## 3. Solvers and their validation

Every claim below is a passing assertion in `tests/solvers.rs`.

### Gravity — Barnes-Hut with retardation

Sources are evaluated at their retarded position. The naive version of this idea
is famously wrong: aberrating the position alone produces a tangential force
that does net work on a bound orbit and unbinds it — Laplace's objection to
finite-speed gravity. The physical field of a uniformly moving source points at
its *instantaneous* position; the retardation is cancelled by the
velocity-dependent part of the field, and what survives is the gravitomagnetic
term proportional to relative radial velocity.

| Test | Result |
|---|---|
| Earth–Sun, 50 orbits, retarded | dE/E = 3.1 × 10⁻⁸, dL/L = 2.1 × 10⁻¹⁴, radius drift 5.3 × 10⁻⁶ |
| Leapfrog convergence order | error ratio 4.0 per halving (second order, confirmed) |

Leapfrog rather than a higher-order non-symplectic scheme on purpose: a
higher-order method has smaller error per step but *secular* drift, which over
the 10⁵ steps between a user's visits shows up as a galaxy that slowly
evaporates. Leapfrog's error oscillates instead of accumulating.

### Hydrodynamics — SPH

Cubic-spline kernel, Monaghan artificial viscosity, optically-thin cooling.
Meshless on purpose: it composes with the tiers on either side without a
remeshing step, and interpolation between grids is exactly where conservation
dies.

| Test | Result |
|---|---|
| Momentum conservation, 20 steps | drift 4 × 10⁻¹⁶ of total momentum content |
| Compressive flow | temperature rises, never falls |
| Cooling curve | positive everywhere; metals cool faster below 10⁷ K; line peak at 10⁵ K exceeds bremsstrahlung at 10⁷ K |

Momentum is exact because the pressure force is written in the symmetric form,
so the force on *i* from *j* is precisely minus the force on *j* from *i*.

### Molecular dynamics

Velocity Verlet, Lennard-Jones with Lorentz-Berthelot mixing, screened and
shifted Coulomb, cell lists, deterministic Langevin thermostat.

| Test | Result |
|---|---|
| LJ minimum position and depth | at 2^(1/6)σ, depth −ε, to 10⁻⁹ |
| Thermostat, 2000 steps at target 300 K | 253 K, stable, no divergence |
| Thermostat noise decorrelation | consecutive steps differ |

That last test exists because of a real bug: the thermostat seeded its noise
stream from the node address alone, so every step drew the *same* numbers. A
"random" force that is identical every step is a constant force, and the system
heated to 81,000 K. The fix threads the step index into the address
(`rng::Stream::split`); the test would catch a regression.

#### What Lennard-Jones is *for*

The van der Waals term applies between bodies that have electron clouds. The
repulsive half is two clouds refusing to overlap and the attractive half is
their induced dipoles, and a nucleon has neither — so applying it to one is a
category error rather than an approximation. Measured, on the carbon atom
scenario, a node the size of an atom holding its twelve nucleons:

```text
    body kind          Nucleon, radius 1.200e-15 m
    min separation     1.9016e-11 m
    σ applied          3.431e-10 m   (carbon, taken from the composition)
    (σ/r)^12           1.190e15
```

The contents left at 150 c. This was invisible for as long as the sampler
inflated the node by 4.7 × 10⁵ — everything sat beyond the 1 nm cutoff and felt
nothing at all — and splitting the binding energy is what exposed it. The rule
is the same shape as the bonded-pair exclusion below and is there for the same
reason: a term that does not describe a pair is removed rather than tuned.

#### Excluded volume, and why an independent draw is not enough

`sample` draws each position without reference to the others, which is an
*ideal gas* draw. That is right for stars in a galaxy and for parcels of gas,
whose SPH kernels are supposed to overlap, and it is wrong for anything solid:
a real interacting system's pair correlation vanishes below contact, because
the repulsion that makes two things two things has already turned them around.

Measured, once the granite and vapour nodes stopped being inflated: the water
vapour node's closest pair sat at `2.55 × 10⁻¹¹ m` against a σ of
`3.12 × 10⁻¹⁰ m` — eight per cent of contact, where `(σ/r)^12` is
`8.6 × 10¹²`. No timestep rescues that. The engine throttled to
`2.83 × 10⁻²⁰ s` and the pair still left at `4.5 × 10⁹ m/s`.

So the configuration is corrected rather than the solver defended. Bodies whose
kind means *one object* — a grain, a molecule, an atom, a nucleon, a star, a
planet — are pushed apart until their surfaces no longer overlap, Jacobi-style
so the answer does not depend on the order the neighbour grid visits pairs, and
re-centred and re-scaled to the node's radius after each pass. A fixed point
always exists for an unstructured body: `child_radius` gives it half the mean
spacing, which fills exactly one eighth of the node's volume whatever the
count, well under the 0.64 of a random close packing.

Where it does *not* exist the sampler says so rather than pretending, in
`SampleReport::worst_overlap`. An iron nucleus is the standing case: a nucleon
is given the 1.2 fm that is the radius *per nucleon* in `R = r₀A^(1/3)`, so
fifty-six of them fill their own nucleus exactly, and a packing fraction of one
is above what any arrangement of spheres reaches. It reports 0.45 and is left
alone.

The pass is skipped where the expected number of overlapping pairs —
`n(n−1)/2 · (2r̄/R)³` — is below `10⁻⁶`, which is most samples. For a galaxy it
is `3.6 × 10⁻²⁴`: a star is 10⁹ m across and its neighbours are 10¹⁶ m away.

#### Covalent bonds

Lennard-Jones plus screened Coulomb describes a gas well and a *molecule* not at
all: nothing in it distinguishes two hydrogens that are bonded from two that
happen to be near each other. A water molecule handed to this tier came apart
the moment anything warmed it, because the bond holding it was a van der Waals
well two hundred times too shallow.

Bonds are Morse, `V(r) = D_e [1 − e^{−a(r−r₀)}]²` with `a = sqrt(k / 2D_e)`.
A harmonic bond has the right stiffness near equilibrium and is infinitely
strong — it can be stretched across the box and still pull back — so
dissociation has to be bolted on as a threshold nobody can defend. Morse has the
same curvature at the bottom of the well and flattens out at `D_e`, so a
molecule given more than its dissociation energy comes apart because the
potential ran out, not because a branch fired. Angles are harmonic, since a
molecule loses its shape by breaking a bond rather than by opening an angle to
infinity.

The constants are spectroscopic, and that is what makes them a test rather than
a fit: bond length, dissociation energy and vibrational frequency are not
independent, so fixing any two fixes the third.

| Test | Result |
|---|---|
| H₂ vibrational fundamental | 4403 /cm against an observed 4401 |
| Period against `2π√(μ/k)` | within 0.06% |
| 75% of the well depth | turns around within 2% of the Morse turning point |
| 130% of the well depth | dissociates |
| Net bonded force and torque | 0.000 of both |
| Energy over 50 000 steps | drift 3.5 × 10⁻¹³ |
| H–O–H bend released 30° off rest | settles to 104.50° |

#### Bonds that form

The hard part of reactive chemistry is not deciding when a bond exists. It is
that adding a term to the potential ordinarily *changes* the potential, and a
simulation that gains energy every time two atoms meet is worthless however
plausible its chemistry.

Every bond is multiplied by a smooth switch — one at short range, zero beyond,
half a cosine between. Bonds are created and destroyed only where that switch is
zero, so **the moment of creation and the moment of destruction are both
energetically invisible**. Conservation across a topology change is not patched
up afterwards; there is nothing to patch.

| Test | Result |
|---|---|
| Potential moved by any formation or breaking | 0.000 J |
| Total drift over a run that forms a bond | 3.6 × 10⁻⁶ of a well depth |
| Well converted to motion on approach | 4.22 eV of a 4.75 eV well |
| Valence saturation (O with four H nearby) | O takes 2, each H takes 1 |
| Three loose atoms with a heat sink | O–H at 0.958 Å, H–O–H at 104.5° |

What the bond releases is real: an exothermic reaction warms the gas, and the
conserved tuple shows exactly where the heat came from. A two-atom collision
then comes *apart* again, which is not a failure — two atoms approaching from
infinity have positive total energy and nothing in a two-body collision can take
any of it away. Real recombination needs a third body, which is why the water
test has one.

Three things had to be right and none was at first. The switch enters the force
as well as the energy, and its term is *added*: where the switch is closing it
lifts a negative potential towards zero, which is extra restoring force, and
with the sign wrong a molecule with three quarters of its dissociation energy
escaped. Dispersion has to be faded out over the same range rather than switched
off, because the same wall that would tear a molecule apart also stops two atoms
ever reaching the distance at which they could bond. And the handover has to
happen where neither term is doing anything: at 2.4 bond lengths it is deep
inside the van der Waals wall and leaves a 0.026 eV barrier that is nobody's
chemistry.

Bonded pairs are excluded from the nonbonded sum, 1-2 and 1-3, as any force
field does. This is not a small double-count: two hydrogens sit 0.74 Å apart
with a Lennard-Jones σ of 2.57 Å, so the repulsive term between them is enormous
and entirely spurious. It tore every molecule apart within a few hundred
femtoseconds, and it is what these tests found first. The second thing they
found was a sign error in the angle force — the bend was anti-restoring. The
force-balance test could not see it, since an anti-restoring force balances just
as well as a restoring one; the shape test could.

### Nuclear

Rate-based, from measured cross sections and decay constants. No attempt to
integrate the strong interaction — that is lattice QCD, at supercomputer-months
per femtometre — and this is not a compromise: the tabulated rates *are* the
experimental facts.

| Quantity | Engine | Reference |
|---|---:|---:|
| pp-chain rate at solar centre (1.57×10⁷ K, 1.5×10⁵ kg/m³) | 9.6 × 10⁻⁴ W/kg | ~10⁻³ |
| pp temperature sensitivity, d ln ε / d ln T | 3.76 | ~4 |
| CNO temperature sensitivity | 16.4 | 16–20 |
| H → He energy release | 0.759% of rest mass | 0.71% |
| Mean neutron lifetime | 1270.6 s | 1267.3 s |

The rate coefficients were another real bug: the standard formulae are quoted in
cgs (erg g⁻¹ s⁻¹, g cm⁻³) and the engine is SI. Using them raw made the Sun
100,000 times too bright — at least an obvious failure.

Iron is the floor, and the engine gets that from the binding-energy curve rather
than from a special case: fusing anything lighter *to* iron releases energy, and
fusing past it costs.

### Quantum

Below the decoherence scale the engine stops storing trajectories and stores
occupation numbers, because that is what exists. `regime()` decides, from the
thermal de Broglie wavelength, which description is valid.

| Test | Engine | Reference |
|---|---:|---:|
| Lyman-α (n=2→1, hydrogen) | 121.50 nm | 121.57 nm |
| Mean blackbody photon energy | 2.688 kT | 2.701 kT |
| Uncertainty relation, enforced | ΔxΔp ≥ ħ/2 | — |
| Locating a proton to 1 fm | deposits 5.2 MeV | ~5–20 MeV |

That last row is the engine's cleanest illustration of measurement-as-
interaction: you cannot watch a nucleus without changing it, and the engine
charges you for the privilege.

The uncertainty principle also supplies the engine's most useful bound. A
phase-space volume contains at most `V p³/h³` distinguishable states
(`max_distinguishable_states`), so "simulate down to the subatomic level" is a
*finite* demand per unit volume. An observer cannot request infinite detail,
because there is no infinite detail to request.

## 3.5 Growth and construction

Structured matter is advanced by a developmental program rather than a force
law, and every step is a transaction that must balance before it is applied
(`morph::Transaction::validate`).

| Check | Result |
|---|---|
| First law, per step: in + released = stored + warmed + radiated | exact to 10⁻⁹ relative |
| Second law: local + exported entropy | never negative, asserted every step |
| Node energy change vs net boundary flux, over a decade | mismatch 2.1 × 10⁻¹⁶ |
| Mass, composition, baryon number under growth | unchanged exactly |
| Structural conservation, 4 programs × 3 masses × 3 budgets | worst 3.6 × 10⁻¹⁶ |

Two calibrations are worth recording because the obvious values are wrong by
large factors.

**Photosynthetic yield is an ecosystem number, not a leaf number.** The
laboratory quantum efficiency is around 3%; a temperate forest actually fixes
~1.2 kg of dry matter per m² per year under ~200 W/m² of mean insolation, which
is 0.32% of incident energy. Using the leaf figure makes trees grow about a
hundred times too fast — plausible for a few frames, absurd after a simulated
decade.

**The whole incident flux has to be on the books, not just the usable part.** A
leaf absorbs the full solar flux and stores 0.3% of it; the rest leaves as heat
and infrared. Booking only the fraction that gets used describes a perfectly
efficient converter, and the second-law check correctly rejects it — a device
that turns all its input into stored free energy while lowering its own entropy
is not allowed. This was caught by the validator rather than by inspection.

Carrying capacity is not a constant anyone typed in. Light capture scales with
crown *area* and maintenance with *mass*, so the two balance at a finite size
and the ceiling emerges from the allometry — which is the actual biological
reason trees stop growing. Height follows McMahon's elastic-similarity result,
`H ∝ V^(1/4)`, calibrated so a one-tonne tree stands about 15 m. The engine
grows a tree from 2 kg to 8.5 t and 30 m over 160 simulated years, then
saturates near 31 t and 40 m.

Growth also responds to conditions without anything being arranged: a structure
held below freezing loses mass to maintenance, and one in drought or shade grows
slowly. That is what makes tree rings.

### A third case: something that was made

`PLAY.md` D15. A composite — a box, a crate, a shelf — neither grows nor is
constructed by a program. It is a **parts list the engine generated** by
assessing what somebody put together, stored so the thing can be sampled again.
It does not grow, it has no light budget and no allometric ceiling, and it
ages only because weathering still applies to it.

Two things follow from "somebody made it" rather than "it grew", and both are
physics rather than bookkeeping:

**Its size does not follow from its mass.** A grown thing derives its extent
from what it has accumulated; an assembled one is the size its parts are, and
its mass follows from them. The arrow between geometry and mass points the
other way.

**Its joins are the size they were made, not the size they need to be.**
Fully-stressed sizing — putting material where the load is, in proportion to
what each member carries — is what thirty years of growth does to a tree, and
regenerating one without it would regenerate a tree that could not stand up. A
box of 25 mm planks glued along their edges has had no such process applied. Re-
sizing its seams until they could carry their load would make every join as
strong as it needed to be, which is indistinguishable from the box never coming
apart. So the seam is stated: a plank glued along one edge is held by a joint of
its thickness by its width, and it comes off long before it snaps.

That is also why joining and breaking are one transform rather than two
features. A weld, a glue line and a grown-together seam differ in what the join
is made of and how much of it there is, and the strength follows from the same
Griffith law as everything else. Break the join and the part becomes a thing of
its own, carrying its own recipe — still a wall, not a lump of mass.

## 3.5A What a thing is made of

`PLAY.md` D13 and D14. A material used to be thirteen numbers somebody typed,
reached by asking what *generated* a node; §2A measured what that cost, which
was that a boulder had no surface and could not collide. Every number is now
derived from `chem::Properties`, which is derived from atoms and bonds.

| property | law |
|---|---|
| density | the substance's own, times how much of the volume the process filled |
| stiffness `E` | cohesive energy density, `3U/Ω`, times the packing squared (Gibson–Ashby) |
| surface energy `γ` | a quarter of an atom's cohesive energy over the area it presents |
| strength | Griffith for a brittle solid, Hall–Petch for a metal, weighted by how metallic it is |
| flaw scale `a` | **the formation history** — a grain if it froze, a layer if it was laid down |
| specific heat | Dulong and Petit, `3R` per mole of atoms |
| thermal limits | `0.4 T_m` to `T_m`, the homologous-temperature rule |
| ductility | metallic from electronegativity, polymeric from molecule size |
| tensile ratio | whether the solid has cleavage planes to open and close |
| resistivity, combustibility | **no law**; stated, and said to be stated |

### Strength is a property of how well a thing was made

Theoretical strength is about `E/10` and real materials are one to three orders
weaker, because they fail at defects rather than everywhere at once. Chemistry
cannot answer that and history can, which is why the one stored number is a
*length with physical meaning* rather than a stress somebody chose.

A solid that **froze** has grains, and the grain size comes out of the
competition between nucleation and growth. The undercooling is solved for rather
than stated — a melt cools past its freezing point, nuclei accumulate, and the
transformation runs away when the latent heat the growing grains release
overtakes the heat being extracted. Turnbull's relation supplies the
solid–liquid interface energy from Richard's rule, and reproduces iron's
0.20 J/m² without being shown it.

A solid that was **laid down** has no grains. Its flaw is the increment — a
growth ring's cells, a course of blocks, a deposited layer — and the generator
knows it, because the generator is what laid it. A tree lays down thirty-micron
cells, and Griffith on thirty microns gives green wood 5.16 × 10⁷ Pa against the
4.5 × 10⁷ the retired table held.

**Two substances, four materials.** Masonry and bedrock are both silicate;
green wood and dry timber are both cellulose. A per-material stress could only
ever record those as four materials, and a formation history tells them apart:
a brick is fired and laid in courses, a pluton crystallises over a geological
age, a fast summer ring has wide cells and a slow one has fine ones.

### What it reproduces, and the one it does not

Seven of eight land within a factor of 2.3, and the ordering is the retired
table's own except that bedrock comes out last where the table put it third.
That is the intragranular-crack gap — a granite at 130 MPa implies a 34 µm
crack against its 21 cm grain — and it is recorded in `BACKLOG.md` with the
measurement rather than papered over. `tests/material.rs` holds the ordering,
and holds bedrock's position in the direction it actually comes out, so that
closing the gap *fails* the test and gets read.

### What a solid presents

`PLAY.md` D18. A solid's boundary is a **union of solid convex primitives**,
each carrying its own material, emitted by whatever generated it.

**A primitive is always a filled solid. Never a shell, never hollow.** That rule
is what makes the vocabulary unambiguous, and it is the general statement of a
measured failure: one convex hull over a box encloses its own cavity, so
anything inside reads as deeply interpenetrating on every frame. A hollow box is
six solid slabs generated together, and a void is not represented at all — it is
simply where no primitive is.

**The generator emits the pieces; nothing infers a decomposition.** That is what
keeps convex decomposition — normally the hard, unsolved half of this problem —
from arising. A wall with a doorway emits four boxes around the opening, because
the recipe is what put the opening there. There are three generators and no
fourth: an **assembly** states its parts, each a filled box with its own
material; a structure states its members, each a capsule; and unstructured
solid matter states one piece, because a rock is one filled solid with no
cavity in it.

The two primitives are both sphere sets, so there is one narrow phase over
both. A capsule is the hull of its two end spheres exactly. A **slab** is the
hull of eight spheres at the inset corners of a box — which is that box with
its edges rounded by the inset, a convex solid in its own right rather than a
discretisation of a sharp one. Flatness is what a capsule cannot do and a panel
needs: a probe standing off a 1.2 m panel reads the same gap at its centre, two
thirds out and across the diagonal, where a row of capsules scallops by 0.169 m
between members.

The surface is **derived once and stored**, and regenerated when the node's
`epoch` moves — which is precisely the definition of when its arrangement
changed. An undisturbed tree computes its boundary once and reuses it for a
thousand frames.

**Material attaches per primitive**, which only an assembly has ever needed: a
tree is wood everywhere and a sampled rock is one substance throughout, so a
single measured material was enough until a box of oak panels on a steel frame
existed. A contact reads the material of the piece it actually struck.

**Its materials are a partition of the node's solid pools**, checked at bake
time rather than documented. A node whose mixture is all water while its
primitives present steel would be telling a contact something its own bulk
contradicts, which is the second axiom failing by way of having two answers to
one question. `Stats::surface_mismatches` reports it.

## 3.6 Structural failure

Joints carry a cross-section, a material and a remaining integrity. Loads
accumulate from the leaves inward in one O(n) pass; peak fibre stress is the
bending moment over the section modulus, `4M / (pi r^3)`, plus the axial term.

| Load case | Result |
|---|---|
| A 13 m, 900 kg tree under its own weight | peak utilisation 0.285 — safety factor 3.5 |
| Wind, 15 and 25 m/s | nothing breaks |
| Wind, 40 m/s | limbs come down |
| Wind, 60 m/s | crown destroyed |
| 600 mm dry powder snow (100 kg/m³) | utilisation 0.39, nothing breaks |
| 300 mm settled snow (200 kg/m³) | utilisation 0.89, nothing breaks |
| 100 mm wet snow (400 kg/m³) | utilisation 1.90, limbs come down |
| Lightning, 10⁷ / 10⁸ / 10⁹ J | 3 / 4 / 6 members destroyed along the channel |
| Brief low ground fire | nothing consumed; trunk reaches 305 K |
| Sustained crown fire | fine fuel consumed; trunk reaches 526 K, not 1100 K |

Every row is a passing assertion in `tests/topology.rs`.

Redundant structures do not use that O(n) pass. They go to a three-dimensional
Euler-Bernoulli frame solver — six degrees of freedom per node, matrix-free,
Jacobi-preconditioned — with Euler buckling and elastic-perfectly-plastic
redistribution on top. Every case in `tests/frame.rs` has a closed-form answer:

| Case | Solver | Closed form |
|---|---|---|
| Axial extension | `PL/EA` | exact to six figures |
| Cantilever tip | 0.108650 m | `PL³/3EI` = 0.108650 |
| Simply supported midspan | 4.420971 × 10⁻² m | `PL³/48EI` = 4.420971 × 10⁻² |
| **Fixed-fixed midspan** | 1.105243 × 10⁻² m | `PL³/192EI` = 1.105243 × 10⁻² |
| Redundant three-bar truss | 585.79 / 292.89 N | `P/(1+2cos³θ)` = 585.79 / 292.89 |
| Buckling at ½, 0.95, 1.2, 3× `P_cr` | utilisation 0.500, 0.950, 1.200, 3.000 | |
| Ductile vs brittle load spread | 1.00× vs 2.40× | |
| Preconditioned vs not | 242 iterations vs no convergence | |

The fixed-fixed row is the one that matters. The previous solver was a spring
network, which carries no moment between a member's ends and so cannot tell a
fixed end from a pinned one; it would report that beam sixty-four times too
flexible. The cantilever row, by contrast, proves nothing — the spring model was
*built* from the cantilever's tip stiffness and passes it by construction. A
stress check alone would also read only 0.0197 at the Euler critical load, a
factor of 51 of false margin, which is why buckling is a separate criterion.

### 3.7 Design, and coming apart

A generated structure is analysed and proportioned the moment it is
materialised. The generator decides where members go; it cannot know what any of
them will carry, so its radii are a shape scaled as a group to match the
structural mass. Fully-stressed design fixes that: analyse, give every member
the section that brings it to the utilisation every other member is at, rescale
the set back to the volume it started with, repeat.

| Structure | Peak utilisation | Spread | Structural volume |
|---|---|---|---|
| Tree, 900 kg | 0.628 → 0.343 | 0.103 → 0.057 | error 1.0 × 10⁻¹⁵ |
| Tree, 6 t | 0.797 → 0.368 | 0.131 → 0.068 | error 1.4 × 10⁻¹⁴ |
| Tower, 3000 t, braced | 0.080 → 0.014 | 0.013 → 0.005 | error 2.1 × 10⁻¹⁵ |

The material does not change; where it is does. Members are sized against an
*envelope* — a vertical overload standing in for anything that settles on a
structure, plus the program's design flow from three directions — because sizing
against one case gives a structure that is optimal for that case and brittle in
every other. A wind from six directions none of the cases used reaches 0.417
utilisation on a tree designed for 20 m/s.

Breaking apart is then a re-analysis. A break is a member ceasing to be
*supported*, and everything hanging off it comes away as its own object with its
own roots — which matters because the static analysis walks the support forest
towards its roots, and a branch still carrying its old support index is still,
as far as any analysis knows, being held up by the trunk it fell off. The piece
falls, tests its members as capsules against what is below, and hands over the
reaction through the ordinary mechanism vocabulary. On a 4 t tree in a 38 m/s
gale: 35 joints break into 24 pieces, 704 kg falls from up to 13 m, 178 members
are struck, and four fail under impact at a peak utilisation of 2.26.

Five things had to be right there, and the failures were all silent. A
fragment's roots are *free*, not fixed — the distinction is not in the topology,
which has no opinion about the planet — and built wrong every piece hangs in the
air exactly where it broke. The ground is not at zero, because a structure is
recentred on its own centre of mass. An impulse is not a force. The whole piece
is behind the contact, not the joint that touched, so using the local lumped
mass under-reads the force by the number of joints in the piece. And impulses
have to accumulate into one analysis per structure per step, because each
analysis regenerates the structure and renumbers the members the next impulse
was aimed at.

### 3.8 Structural dynamics

Newmark-beta on the same operator, with lumped mass and rotational inertia from
the members' own geometry. `tests/dynamics.rs`:

| Case | Solver | Closed form |
|---|---|---|
| Cantilever period, integrated | 0.472675 s | `2π/((1.8751/L)²√(EI/ρA))` = 0.472050 |
| Same, by Rayleigh quotient | 0.466388 s | within 1.2% |
| Dynamic load factor, step load | 1.985 | exactly 2 |
| Rigid-body translation | 5.2 × 10⁻¹⁶ J of strain | 0 |
| Undamped energy, four cycles | 99.3% retained | 100% |
| Damping ledger | closes to 1 part in 10⁶ | |
| Tower under held wind | settles to the static answer to 0.000% | |
| Tree, 18 m/s gust released | swings to +1.47 m, back to −1.24 m, four crossings | |

The dynamic load factor is the row worth dwelling on. A load that arrives
suddenly deflects a structure twice as far as the same load standing still, so a
quasi-static analysis of a gust under-reads the stress in it by a factor of two
— the difference between a member at 60% utilisation and one that has already
failed. No amount of care in a static solver recovers that number.

The tower row is the consistency check: hold a load steady long enough and the
dynamics must settle to exactly what the static analysis predicted. Running the
dynamic step on the static operator is what makes that a theorem rather than a
coincidence.

The same consistency is checked the other way in `tests/topology.rs`: on a
1500-member determinate tree the frame solver and the O(n) statics pass agree to
**7 × 10⁻¹¹**, member by member. A determinate structure's internal forces are
fixed by equilibrium alone, so that is not "close enough" — it is the only
answer a correct beam model can give.

That test found two real things. It had been passing because the redundant
solver could not converge and was silently falling back to the statics it was
being compared against. And once it did converge, the two disagreed by up to
10% on the least-loaded members: lumping a member's own load half to each end
makes the *end moments* right but halves the axial force, because the base
section carries all of that load and the tip carries none. Adding the missing
half back is what closes the gap. A safety factor of 3.5
against self-weight matches measured values for real trees, which is the check
that makes the rest of the table meaningful — a model that could not stand up
would fail everything else for the wrong reason.

The wind and snow rows are the ones worth dwelling on. A tree that survives a
25 m/s gale and loses limbs at 40 m/s is behaving as trees do, and the
distinction between dry and wet snow is not a tuning knob but the actual reason
snow damage happens: interception capacity scales steeply with how well the snow
adheres. Lightning follows the support chain to ground and deposits energy in
proportion to each member's resistance, so thin members — highest resistance per
kilogram — are destroyed while the trunk survives. Fire heats members on a time
constant set by their thermal mass, which is why a ground fire takes the
understory and scorches but does not fell mature trees.

## 4. Generated detail is statistically correct

Conservation is necessary but not sufficient: a cloud whose parcels conserve
energy but follow the wrong velocity distribution is detectable by anyone with a
spectrograph. `tests/statistics.rs`:

| Property | Engine | Expected |
|---|---:|---:|
| Velocity kurtosis ⟨v⁴⟩/⟨v²⟩² | 1.667 | 1.6667 (Maxwell-Boltzmann) |
| Velocity anisotropy between axes | 0.17% | 0 |
| Kroupa IMF high-mass slope (per octave) | m^−1.29 | m^−1.3 |
| Plummer profile r₅₀/r₂₅ | 1.6049 | 1.6086 |
| Poisson(25) counting variance | 25.24 (mean 24.98) | 25 |
| Uniform draws, 64 bins | χ² = 64.0 (63 dof) | ≈ 63 |
| Composition scatter, mass-weighted mean | exact to 10⁻⁹ | exact |

The kurtosis test is deliberately sharp: it fails immediately for a uniform or
top-hat distribution, both of which would pass a naive "the mean speed looks
right" check.

## 5. Modelling choices worth knowing about

- **Dark matter is a static halo potential**, not particles mixed into the
  baryonic composition. Both are standard; only the first keeps the books
  straight, since dark matter carries no baryon number and folding it into the
  composition would have the engine report 10⁶⁷ baryons for a galaxy containing
  10⁶⁶.
- **Temperature is derived from random kinetic energy**, never averaged from
  children. Averaging temperatures across unequal masses is wrong in general and
  wrong by orders of magnitude across this dynamic range.
- **Entropy may only increase** under coarse-graining — but the monotonic
  quantity is the *total*, local plus exported, not the local part. A growing
  structure legitimately lowers its local entropy, and the original one-line
  `max` in `Tree::coarsen` would have clamped it back up and silently
  unbalanced the books the moment anything in the world became more ordered.
- **Cooling is a three-regime fit** — molecular/fine-structure below 10⁴ K, the
  Lyman-α peak to 10⁷ K, bremsstrahlung above — with metallicity scaling the
  line cooling. This is why the first generation of stars forms differently from
  later ones, which the engine reproduces for free.
