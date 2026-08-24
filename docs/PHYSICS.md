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
| Continuum | 10⁻⁸ – 10⁴ m | 1 ms | 3.3 ns | Fluid parcels, grains, bulk matter | SPH |
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
E = M c²  +  K_bulk(M, P)  +  U_internal  +  Φ_binding  +  Φ_external
```

- `K_bulk` is the **exact** relativistic form, `√((Mc²)² + (Pc)²) − Mc²`, not
  `p²/2M`. The decomposition has to be invertible to the last bit — `restrict`
  recovers `U` from `(M, P, E, Φ)` — and a Newtonian bulk term makes the round
  trip lossy at the 10⁻⁵ level for anything moving at galactic rotation speeds,
  which shows up as visible energy drift when a user pans across a disc.
- `Φ_external` is separate from `Φ_binding` because it is *not recoverable from
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
