# Liquid solutions: the water-law checklist and the plan for a general solution

*Written 2026-10-18 at the owner's request, as part of Phase 6 (`PLAY.md` E8e).
It is a plan, not a result: nothing in section 3 onward has been run. Reference
values for experiments are the author's recollection of the literature and are
marked as such; each is to be looked up and cited when the test that uses it is
written, because a test that compares against a remembered number is a test
against nothing.*

Contents: 1. why this document exists · 2. the water-law checklist · 3. what the
plan (`PLAY.md`) already says that this builds on · 4. the general solution ·
5. deriving the phase diagram · 6. the test set and what each test can resolve ·
7. stages, in order, with done-whens · 8. questions the plan does not answer
(the owner's) · 9. cost.

---

## 1. Why

The water law (`runs/e8d/core/law-core-B5.txt` and its refits) is the only
intermolecular law the engine has derived to a liquid. It is close in density
(equilibrium 1.003–1.007 g/cm³ against 0.997, kappa_T 4.3e-5 per bar against
4.5e-5) and short in structure (g(OO) at 3.4 Å 0.98 against 0.8; at 4.4 Å 0.98
against 1.1; no tetrahedral gap). A solution is that law plus a second kind of
molecule, so every fault of the water law is inherited and every term has to be
defined for unlike pairs. The checklist comes first because it is the list of
places where a solution can go wrong that water already shows, and the plan
after it is built so that each place is *measured* rather than assumed.

## 2. The water-law checklist

Compiled from a literature list of what sets inter- and intramolecular forces in
a liquid (general physical-chemistry and simulation knowledge: Stone's
*Theory of Intermolecular Forces*, the SAPT and HIPPO/MEDFF/AVDO/DensIP papers
already cited in `PLAY.md` E8b, Gillan–Alfè–Michaelides 2016 for water, and the
nuclear-quantum-effects reading in `runs/many/E8d-notes.md` §15).

Status key: **D** derived from the engine's own electrons and tested ·
**F** fitted to the engine's own MP2 energies (a stored shortcut under the E8
conditions) · **P** present in part · **A** absent · **?** not measured. "Test"
is what would settle the row; a row is closed by a number, not by a judgement.

### 2.1 Intermolecular

| # | Factor | In the law now | Status | Evidence in the project | What settles it |
|---|---|---|---|---|---|
| 1 | Exchange repulsion: overlap of closed-shell densities | Gaussian-cloud overlap `K exp(-r²/(sa²+sb²))` between the clouds, plus a Born–Mayer core on O–H (exponent 7) that closes the short-range holes. | F | E8b (m), E8d §10.6; holes gone, pair error 0.17 held out. | Out-of-sample slopes at 3.1–3.8 Å are still 0.07–0.22 kcal/mol/bohr too attractive for every law (E8d §9). Needs the anisotropy test below and a repulsion derived from the density, not fitted per site pair (AVDO form). |
| 2 | Repulsion anisotropy (lone-pair and bond-pair directions) | Isotropic in each site's cloud; zero-charge lone-pair repulsion sites were driven to zero. | A | E8d §10.4, and the antenna test: energy error correlates 0.00 with the O–H···O angle. | Not supported by the pair data. Re-test on bent H···O contacts at 95–125° using MP2 slopes, where the shared error sits. |
| 3 | Permanent electrostatics: multipoles | Gaussian charge clouds (O 5.44 e, bisector 1.65 e, H 0.45 e), derived from the monomer's HF potential and exact electrostatic energies, dipole scaled to MP2's 1.865 D. | D (fit to the engine's own exact first-order energies) | E8b (k)–(l): 0.26 fitted / 0.41 held out against 1.9 for point charges. | Closed to 0.009 kcal/mol beyond 5.5 Å. Remaining: the quadrupole against the engine's MP2 value (not reported). |
| 4 | Charge penetration | Included: overlapping Gaussian clouds, nuclei as points. | D | Same. Point charges were 1.1 / 2.0 / 4.0 kcal/mol too weak at 2.8–3.1 / 2.5–2.8 / <2.5 Å. | Closed. |
| 5 | Long-range treatment of the sum (Ewald or truncation) | Ewald with Gaussian smearing; tin-foil boundary. | D, validated | E8d §7: energy, forces, torques and virial to 1e-8; truncation cost 1.2 kbar and 0.75 kcal/mol/molecule and cost the structure nothing. | Tin-foil against vacuum boundary for the slab geometry: **?** (needed before a coexistence run; the slab is not isotropic). |
| 6 | Induction: polarisability | Per-site isotropic polarisabilities (O 9.234, H 0.854 bohr³) from the finite-field response of the molecule's HF and MP2 energy; molecular 10.05 against MP2 9.908. | D | E8a. | Anisotropy of the molecular tensor against MP2's. **?** |
| 7 | Induction: non-additivity | Explicit induced dipoles solved self-consistently. Hartree–Fock three-body terms: the law's induction is about 1.5× too strong per triple (slope 0.64, correlation 0.75). | P | E8d §1, §6.1. | The pressure budget says each 5% of polarisability is ~1,100 bar. Scan done in the liquid (E8d §16): structure improves with induction and pressure and cohesion worsen by the same step, so one scale does not give both. |
| 8 | Higher polarisabilities, field-dependent response (hyperpolarisability) | None. | A | — | Compare dipole-in-field of the molecule at the liquid's fields (~0.1 au) with the linear response. **?** |
| 9 | Induction damping at short range (Thole-type) | Self-consistent: Gaussian widths set by the polarisability, `sigma = (sqrt(2/pi) alpha/3)^(1/3)`. | D | E8a. Widening to 1.2× made the liquid looser, not better (E8d §17). | Closed in that direction. |
| 10 | Exchange–induction quenching (polarisation weakens at contact) | None. | A | Triples say induction is weaker than modelled at contacts (E8d). | The physics a point-dipole model lacks; test by fitting an induction damping that depends on overlap and checking whether triples, pressure and g(r) move together. |
| 11 | Charge transfer | None. | A | Not modelled; bundled into the pair fit. | Energy decomposition of the MP2 dimer (SAPT-like, or constrained-DFT charge transfer). Not run. |
| 12 | Dispersion C6 | Derived from the vdW-DF kernel split by the Hirshfeld partition (O–O 17.48, O–H 6.33, H–H 2.38; molecule 52.3). Free fit moved it to 35 (O–H only) and 48–51 once the repulsion had a second shape. | D, with the fit disagreeing | E6/E7 and E8b. | Which to trust is open; the molecule–molecule C6 measured ≈ 45–50 is the check (E8b Casimir–Polder from the MP2 frequency-dependent polarisability would be an independent derivation). **?** |
| 13 | C8, C10 | None. | A | A free C8 earlier absorbed the repulsion's error (E8b). | Derive from the same kernel or from higher multipoles; do not fit. |
| 14 | Dispersion damping at short range | Tang–Toennies, exponent fixed. | P | — | Sensitivity of the pressure budget to damping. **?** |
| 15 | Many-body dispersion (Axilrod–Teller, screening) | None. | A | HF triples cannot see it (no correlation). | MP2 three-body correlation part is a few hundredths of a kcal/mol on compact triples (E8a stage 3); larger clusters unmeasured. |
| 16 | Dispersion tail beyond the cut-off | Added analytically (−117 to −214 bar, −0.05 kcal/mol/molecule) | D | E8d §6.8. | Closed for bulk; slab needs the inhomogeneous form. |
| 17 | Hydrogen bonding: directionality and cooperativity | Emergent from items 3, 4, 6, 7; no separate term. | P | Straight bonds fit to 0.23 kcal/mol, bent contacts 0.7 (E8b j); 3.5–3.6 hydrogen bonds per molecule against ~3.7 for TIP4P-like control; q 0.57 against 0.72. | Cooperativity: MP2 energies of the 4–6 molecule ring and chain clusters against the law (hexamer not yet run). |
| 18 | Three-body (and higher) non-additivity of the total | Only through induction. | P | HF triples −0.12 ± 0.08 (n = 3) to −0.35 ± 0.18 (n = 6); MP2 triples −0.10 ± 0.08; compact clusters cannot sample it (non-additive induction −2.43 kcal/mol/molecule in the liquid). | Cluster cut from the liquid, 10+ molecules, law against reference. The one centre taken to 10 went +1.11 at n = 10; **?** whether that is signal. |

### 2.2 Intramolecular

| # | Factor | In the law now | Status | Evidence | What settles it |
|---|---|---|---|---|---|
| 19 | Bond stretch, angle bend | Rigid at the ground-state average (0.97894 Å, 104.078°) in the working law; options 2 and 3 flexible on the MP2 surface (35-term quartic, rms 3.1 cm⁻¹; harmonics 1626 / 3847 / 3973 against 1649 / 3832 / 3943). | D | E8c. | The flexible liquids' own table; the working law is rigid. |
| 20 | Anharmonicity | Quartic surface. | D | E8c. | Fundamentals (3657, 3756, 1595 measured) not yet checked against the surface. |
| 21 | Cross terms (stretch–bend) | In the quartic. | D | E8c. | — |
| 22 | Environment changes the geometry (O–H lengthens, angle opens in the liquid) | Option 2 and 3 only: the cloud follows the angle (`GeoCloud`); bonds held in option 2. | P | E8c: the angle-following cloud gives 105.4° where a fixed cloud collapses it to 95°. | Option 3 run. |
| 23 | Dipole enhancement 1.85 → 2.4–3 D | Through induction. | ? | `phys-bulk` does not report the mean dipole. | Report it. One line of output; compare with the 2.4–3.0 D of ab initio liquid water (a literature range, to be cited). |
| 24 | Intramolecular non-bonded and conjugation | Not applicable to water (no 1–4). | n/a | — | Becomes real for methanol and every larger molecule: a torsion (methyl rotor) is the first test of treatment (2). |
| 25 | Electronic-state effects (open shell, spin) | Closed shell. | n/a | — | Not for water; ions here are closed shell. |

### 2.3 Condensed-phase and statistical

| # | Factor | In the law now | Status | Evidence | What settles it |
|---|---|---|---|---|---|
| 26 | State point: density, temperature | Fixed-density runs at 298 K, a density scan 0.95 / 1.0 / 1.05. | D | E8d §10.7: equilibrium density 1.003–1.007. | NPT runs and the temperature of maximum density (277 K measured, absent in the plan's fitted law). |
| 27 | Thermal sampling and length | 30–60 ps, 216 molecules. Pressure scatter ±600 bar instantaneous; ±89 bar after 60 ps. | P | E8d §10.7. | Sampling error stated per observable before a comparison counts. |
| 28 | Nuclear quantum effects | Classical nuclei. Literature (Reddy; E8d §15): density shifted ~0.3 %, structure softened: lower first peak, broader correlations. | A, estimated | E8d §15 reads as small on density, wrong direction to explain the missing gap. | Path-integral run on the flexible law (planned, `PLAY.md` E8c). Heat of vaporisation needs the intermolecular zero-point (~2.8 kcal/mol in ice) added before it can be compared with 10.5. |
| 29 | Environmental screening (effective dielectric, effective polarisability) | Self-consistent induction supplies it. | P | — | Dielectric constant from the dipole fluctuation, against the measured ~78 (to be cited): derived, not tuned, and a strong test of items 3, 6, 7. **?** |
| 30 | Finite size, periodic images | 216 molecules, tin-foil. | P | E8d §6.8: pressure oscillates with the cut-off, converged with Ewald. | Repeat at 512 and 1000 molecules; extrapolate as the owner's solution rule requires. |
| 31 | Dynamics and transport (collective response, diffusion, viscosity) | Not measured. | ? | — | Self-diffusion (2.3e-5 cm²/s) and viscosity (0.89 mPa·s) are *dynamical* checks on the law; they are not in any done-when, so are listed only as a source of evidence. |

### 2.4 Method (what the reference and the fit can get wrong)

| # | Factor | Status | Evidence / action |
|---|---|---|---|
| 32 | Reference: MP2 basis 5 % short in the dimer | Known | Completing the basis (g on O, f on H) takes the dimer from −4.778 to −4.826; not the reason the second shell is missing (E8b h). |
| 33 | Reference: MP2 itself (overbinds dispersion for large polarisable systems; CCSD(T) is the benchmark) | Known | S22 within 4–10 % of CCSD(T)/CBS, under-bound. A CCSD(T) check on the 29 trimer pairs would bound it. **Not run.** |
| 34 | Fit weighting (energies fit; the pressure is a sum of slopes) | Measured | A law good to 0.1–0.3 kcal/mol can be 25 % wrong in the mid-range force (E8b). Slope data in the fit was the cure; out-of-sample slopes are the audit. |
| 35 | Fit conditioning (degenerate starts, collapsed laws) | Measured | `cloudfit.log` "after" line must read under ~1 kcal/mol, dipole ≈ 1.87 D, no charge of tens of electrons (E8c). |
| 36 | Repulsive-wall holes | Closed | O–H core, exponent 7: +33.5 … +2.6 kcal/mol at 0.6–1.4 Å. |
| 37 | Simulator correctness | Closed | A TIP4P/2005-like law gives 3.12 at 2.75 Å and coordination 4.25 (E8b). |
| 38 | Heat of vaporisation compared like with like | Rule | Classical rigid value above 10.5 by ~2.8 kcal/mol is right; agreement with 10.5 is not evidence (E8z list). |

### 2.5 What the checklist says to do first

Rows where the law is **absent or unmeasured and the evidence says it could
matter**, in the order the project's own evidence ranks them: **7 / 10**
(induction strength and its quenching at contact; the 650 bar per 10 % trade
cannot be one scale), **18** (a cluster of 10+ molecules against the reference),
**23** and **29** (dipole and dielectric constant: cheap, derived, no fit), **12**
(which C6), **28** (nuclear quantum effects, expected small), **11** (charge
transfer, which becomes large for an anion). Rows 5 (slab boundary) and 30
(size) are prerequisites for the phase diagram and are listed again below.

The checklist is a **plan for diagnosis**. Whether a row becomes a new term in
the law is a `PHYSICS.md`-weight decision and is the owner's; this document does
not add one.

---

## 3. What the plan already says that this builds on

All from `PLAY.md`, restated so the design below can be checked against it.

- **E8**: the four substances water, methanol, ammonia, methane (373 / 338 / 240
  / 112 K); rigid at E5 shapes; coexistence by a liquid slab against its vapour;
  the boiling point where the vapour reaches one atmosphere. **None of the four
  has had its coexistence run with the working law.** `phys-coexist` exists. The
  only liquids simulated with a derived law are water (216 molecules, many laws)
  and methane (E8d §12: heat of vaporisation 2.01 against 1.96, pressure −222
  bar). Ammonia and methanol have growth states (`grow-ammonia.state`,
  `grow-methanol.state`) and no pair data.
- **E8 route A conditions**: only engine-computed energies go in; every term a
  physical mechanism; residual measured on pairs the fit never saw; what the form
  misses measured against E6's site-site potential and a molecule in a shell (B).
- **E5b, dissolved compounds (the owner)**: clusters in the solution's own ratio,
  all components solved together, ions and nearest shells in full with the rest
  through the pair-and-three law, extrapolated in size; the ratio **bucketed**,
  the bucket as wide as the ratio can move without the derived behaviour
  changing by more than the tolerance; done for salt in water over its range and
  for mixtures of the 45 substances; buckets measured first on representative
  mixtures of three or more components. The thorough derivation runs behind a
  quick one (`BACKLOG.md`).
- **E11 (the owner, 2026-10-11)**: a compound's stored rule is its whole phase
  diagram — melting curve, vaporisation curve to the critical point, sublimation
  curve below the triple point — built from measured heats of fusion and
  vaporisation and volume changes and integrated by Clapeyron; the triple point is
  where the derived melting and vaporisation curves intersect; a mixture's
  depression and elevation follow from the pure compounds' rules; the simulation
  runs once per compound and pressure range and never in the frame loop; every
  shortcut carries the pressure and temperature range it was derived over and
  says when it is outside it. **Checks, not inputs**: carbon dioxide's triple point
  and its sublimation at 194.7 K, nitrogen on Pluto, ammonia and water ice.
- **D17**: a node's mixture is part of its matter. **Backlog items 5 and 7**:
  dissolved solute carries no volume (seawater ~3 % light); water's expansion is
  half the ordinary rate and the density maximum at 277 K is absent.
- **The axioms** (`CLAUDE.md`): nothing the engine knows is named; nothing is told
  that can be measured; a stored rule stays a function of its conditions.

## 4. The general solution

### 4.1 The principle

A solution is not a new kind of calculation. It is the same law, in which every
term is a function of *two sites' derived properties*, applied to a pair that is
no longer two copies of the same thing. The solution's law is therefore
**the pure components' laws plus nothing that names a substance**:

1. **Every species supplies a monomer record**, derived once at intern and
   stored (E9): geometry (E5b), charge cloud and dipole (E8b), site
   polarisabilities (E8a), the dispersion blocks from the vdW-DF kernel and
   partition (E6/E7), the energy surface and its geometry-following terms (E8c).
   A neutral molecule and an ion differ in the record's *values* (net charge,
   closed-shell density, polarisabilities), never in its form.
2. **A pair term is a function of two records**, with no mixing rule typed in:
   - electrostatics and penetration: the overlap of the two clouds (exact for
     the form);
   - induction: each site's polarisability in the field of the other's charges and
     induced dipoles, solved self-consistently (explicit many-body);
   - dispersion: the cross dispersion block from the kernel evaluated on the two
     monomers' densities; Lorentz–Berthelot-type combining rules are **not used**,
     and the kernel's cross value is *checked against* the geometric and
     London combinations as a measure of how wrong a combining rule would have
     been;
   - repulsion: the overlap of the two clouds with the strengths the fit assigns.
     **This is the one place a fit lives**, as in E8: strengths and the core are
     fitted to the engine's own MP2 energies of unlike pairs. Transferability is
     measured, not assumed (§4.4).
3. **Ewald, induction, the integrator and the thermostat are unchanged.** `liquid`
   already carries several `Kind`s (`Kind::from_atoms`), so a mixture of
   species in one box needs a species index and a pair table, not a new
   simulator. This is a statement about the design to be checked in the code
   before it is relied on.
4. **A solution's behaviour is measured, never entered**: the ideal solution is a
   *baseline* to subtract, not an assumption to build from. The excess properties
   (G^E, H^E, V^E and the activity coefficients) are what the simulation
   measures, and the phase diagram is read from them (§5).

### 4.2 What an ion needs that a molecule does not

Sodium and chloride are single closed-shell atoms. What changes:

- A net charge, so the monomer record has a charge and the Ewald sum a neutralising
  partner. A solution is built electroneutral by construction (equal numbers of
  Na⁺ and Cl⁻).
- **Diffuse basis functions for the anion** (chloride holds its extra electron
  weakly). E5b's per-molecule growth finds these for itself; ion records are
  grown, not taken from a table, and checked against PBE's converged ion
  energy as the atoms were at E4.
- **Reference method**: MP2 on ion–water clusters, with the same counterpoise
  rule, because the pair data is what fixes the repulsion strengths. PBE's
  self-interaction error is largest for the anion, which is one more reason
  the reference is wavefunction theory (E8a). The size of the many-body part is
  the finding to expect: ion–water induction is non-additive and charge transfer
  (row 11) is largest here.
- **The data set**: Na⁺–H₂O and Cl⁻–H₂O pairs on the biased draw (`phys-pairs
  --bias`) and wall draw; Na⁺–(H₂O)₂ and Cl⁻–(H₂O)₂ trimers and the ion with its
  first shell, up to the six-coordinate Na⁺(H₂O)₆ and the chloride shell, cut
  from the simulated solution as for water; and the contact and solvent-separated
  Na⁺···Cl⁻ pair in a water shell, which is where the owner's rule (all components
  solved together, in the solution's own ratio) applies.
- **The long range**: the ion's field falls as 1/r, so water beyond the solved
  shell still screens it. Ewald handles the sum in the simulation; in the
  reference, the shell is solved in full and the rest enters through the law, as
  the owner decided for dilute solutions.

### 4.3 Which pair data every solution needs

For each unlike pair of species A, B in a solution:

| Data | Purpose |
|---|---|
| A···B pairs, biased draw from a liquid of the mixture | the repulsion strengths, the cross dispersion check |
| A···B wall pairs (0.6–2.3 Å nearest atoms) | the core, so there are no holes (the E8c/E8d finding) |
| A···B slopes (finite difference along the line of centres) | the force, which fixes the pressure (E8b) |
| A···A···B and A···B···B triples cut from the mixture | the three-body part, induction in particular |
| exact first-order electrostatics of A···B | verify the cloud model on the unlike pair |

Whether the strengths fitted on like pairs predict the unlike energies is the
**transferability test** and is run first (§7, S1).

### 4.4 The transferability test, the experiment that decides how general this is

For each of water–methane, water–ammonia and water–methanol: compute 40 MP2 pairs
with the unlike biased draw; evaluate them with the law formed from the two pure
laws **with no unlike data** (the cloud overlaps and the cross dispersion already
follow from the form) and report the error beside the pure-component held-out
errors. Three outcomes, each defined before the result is known:

- error within 1.5× the pure held-out error: the form transfers, unlike data is
  used only to audit;
- error 1.5–4×: fit the unlike repulsion strengths only, with the rest held, and
  report what moved;
- worse: the form lacks a term for the interaction, and the finding goes to the
  owner before anything is fitted.

### 4.5 A hierarchy of derivations (quick behind thorough)

Following the backlog's rule, the world uses the cheapest rule that is still a
function of its conditions:

1. **Quick**: pure-component curves (from the E8/E11 derivations) combined by
   Raoult's law with an ideal-solution baseline, and the *measured* excess terms
   stored per ratio bucket as activity-coefficient parameters;
2. **Thorough**: the full free-energy calculation of §5, running behind the quick
   one and replacing its excess terms when it completes.

The quick rule's accuracy is judged by comparison with the thorough one, which is
the same use the owner puts the 0.1 % tolerance to.

## 5. Deriving the phase diagram

The plan's E11 already says what is stored: a compound's rule is its **phase
diagram**, not two curves. For a solution the diagram is a function of
temperature, pressure and composition. The derivation in order of what is
available now.

### 5.1 Pure component: liquid–vapour

- `phys-coexist` (direct coexistence, a liquid slab against its vapour) gives
  `rho_L(T)`, `rho_V(T)`, the vapour's pressure from its density, and the heat of
  vaporisation from the energies; Clausius–Clapeyron gives the boiling point at
  any pressure. **To do first: the boiling points of the four substances**,
  which E8 listed as the first test and which has not been run with the working
  law.
- Slab corrections that bulk runs did not need: the dispersion tail for an
  inhomogeneous system, the Ewald sum with a slab-appropriate boundary (checklist
  row 5), and a cut-off long enough that the interfacial tension is not
  truncated. Each is measured, not assumed.
- Critical point from `rho_L − rho_V ∝ (T_c − T)^β` and the rectilinear diameter
  law. Finite size biases it near T_c, so it is quoted with the box size.
- **Pressure rule:** pressure comes from the virial the code already has; a vapour
  at these pressures is treated by its second virial coefficient, computed from the
  same law by direct integration of the pair potential, not the ideal-gas law.

### 5.2 Pure component: solid–liquid and solid–vapour (E11)

Needs a crystal and a free-energy method (Einstein-crystal / Frenkel–Ladd
thermodynamic integration, or coexistence), for ice Ih (proton disorder: a
decision), solid methane (phase I is a plastic crystal), ammonia (phase I) and
methanol (α). E11's order puts molecular ices first, so these are in scope and are
sequenced after the liquid–vapour work they depend on (§7). The triple point is
where the derived melting and vaporisation curves meet; sublimation is integrated
from it by Clapeyron with the Hess heat of sublimation.

### 5.3 Binary solution: liquid–vapour

The general method is by **chemical potentials**, because the composition of a
direct-coexistence slab drifts and a Gibbs-ensemble insertion fails in a dense
hydrogen-bonded liquid:

1. For each component, `mu_i^ex(T, P, x)` in the liquid by thermodynamic integration
   (coupling the molecule on, a staged soft-core path) or Bennett's acceptance
   ratio, at the composition grid of the buckets.
2. The vapour: second virial coefficients `B_ij(T)` from the law, giving the vapour
   fugacities directly.
3. Equilibrium: `mu_i^L = mu_i^V` at fixed T, P determines x and y; the bubble and
   dew curves and the activity coefficients `gamma_i(x)` follow.
4. **Self-consistency checks, which are the tests of the machinery**: Gibbs–Duhem
   (`sum x_i d ln gamma_i = 0`) satisfied by the measured `mu_i`; the pure-limit
   `mu_i^ex` equal to the pure-component calculation of §5.1; and the same point
   reached by direct coexistence of a mixed slab at one composition.

### 5.4 Binary solution: liquid–liquid (miscibility)

From the curvature of the Gibbs energy of mixing `G^mix(x)`: spinodal where
`d²G/dx²` changes sign and binodal by the common tangent. Water–methane is the
case where it matters (the second liquid is nearly immiscible); it is a
solubility calculation, `mu_CH4` in water against pure methane, and it is a test of
the hydrophobic hydration the law gives, which is a known weak point of
pair-additive laws.

### 5.5 Binary solution: solid–liquid and the eutectic

The freezing point of a solution whose solute does not enter the ice is where the
solvent's chemical potential in the liquid equals the ice's: `mu_w^L(T, x) =
mu_ice(T, P)`. So it needs the activity of water in the solution (§5.3) and the
ice's free energy (§5.2). The same chain gives the eutectic for salt–water. Hydrates and
clathrates (ammonia hydrates; methane clathrate, stable at several MPa and low
temperature) are phases the plan has not named; whether they are in scope is for the owner
(§8).

### 5.6 Dilute electrolyte (the saline case)

At 1 % by mass the solution is the dilute limit, where a different route is
cheaper than free energy of the solvent. The arithmetic the design rests on:

```text
1 g NaCl : 99 g water   -> 0.017111 mol : 5.4953 mol, one NaCl per 321.2 waters
smallest exact box      -> 321 H2O + 1 Na+ + 1 Cl-   (1.0006 % by mass)
molality 0.1728 mol/kg, concentration ~0.17 mol/L
ion-pair separation ~ 19 A at this concentration, against a 22 A box
```

One ion pair in a 321-molecule box is dominated by finite size, so the solution is
run at 321, 642 and 1284 waters with 1, 2 and 4 pairs and extrapolated, as the
owner's rule for solutions says. What the sizes can and cannot resolve:

```text
quantity (author's recollection, to be cited)            size            resolvable?
water activity: ln a_w ~ -phi nu m M_w, phi ~ 0.93        -0.0058         -> mu_w shifts 0.0034 kcal/mol: NO by direct measurement
boiling-point elevation  (K_b 0.512)                      ~0.16 K         via activity, not by heating
freezing-point depression (K_f 1.86)                      ~0.60 K         via activity and ice mu
osmotic pressure                                          ~8 bar          YES: a pressure, resolvable at +-100 bar? marginal
density increment (1.0053 against 0.9982 g/cm3)           +0.7 %          marginal: a ~150 bar signal in a pressure that scatters by ~100
apparent molar volume of NaCl (~16.6 cm3/mol)             -               YES with long runs
Na-O and Cl-O first peaks (~2.4 and ~3.1 A), coordination ~5.5-6   -    YES, structural
single-ion hydration shells, ion pair PMF                 -               YES
```

So the plan is: **derive the solute's excess chemical potential by thermodynamic
integration of the ion pair and obtain the solvent's activity by Gibbs–Duhem**,
because the solvent's own chemical potential change (0.003 kcal/mol) is well below
what direct sampling resolves; resolve **structure, apparent molar volume and the
ion pair's potential of mean force** directly; and **state the resolution** of
each comparison before reading it. The colligative values are checked
against the activity route's output, which is the honest test, and the direct
heating or freezing is not attempted.

## 6. The test set

The user's list is "water, methane, ammonia, alcohol, and a 1 % saline solution".
The tests have two layers: the **pure liquids** (every liquid is a trivial
one-component solution and their diagrams are the foundation of every mixture) and
**solutions** built from them. Reference values are from the author's memory and
must be verified from sources before a test uses them.

### 6.1 Pure components

| Substance | T_b (K) | T_m (K) | T_c (K) | P_c (MPa) | note |
|---|---|---|---|---|---|
| water | 373.15 | 273.15 | 647.1 | 22.06 | triple 273.16 K, 611.66 Pa; density maximum 277 K; ice floats |
| methane | 111.7 | 90.7 | 190.6 | 4.60 | triple 90.69 K, ~11.7 kPa |
| ammonia | 239.8 | 195.4 | 405.4 | 11.33 | triple 195.4 K, ~6.1 kPa |
| methanol | 337.8 | 175.5 | 512.6 | 8.10 | the substance in E8's set |
| ethanol | 351.4 | ~159 | 514 | 6.15 | the alcohol most people mean |

### 6.2 Solutions

| Solution | what the diagram must show | why it is a test |
|---|---|---|
| NaCl in water, 1 % by mass | a < 1 K elevation and < 1 K depression; density 1.0053; ion structure | dilute electrolyte, charge, finite size |
| water + methanol | miscible at all ratios; mild positive deviation from Raoult; no azeotrope at 1 atm | polar–polar, both hydrogen bonded |
| water + ammonia | miscible; strongly negative deviation | unlike donor/acceptor balance; N has one lone pair for three N–H |
| water + methane | near-immiscible; solubility x ~ 2.5e-5 at 1 atm and 298 K | hydrophobic hydration; the case a pair law is least likely to get |
| methanol + methane, ammonia + methane | partially miscible at low T | dispersion-only partner against a polar one |
| water + ethanol, at the azeotrope and at x = 0.10 and 0.50 | azeotrope near 95.6 wt % (78.2 °C); the boiling curve and activity coefficients at the two other compositions | an azeotrope derived from the law and nothing else is a sharp test; the other two ratios show where the curve departs from Raoult |

## 7. Stages

Each stage ends with a measurement; none is done because its code exists. Numbers
are the done-whens, to be agreed with the owner; where one is an assumption it
says so.

**S0. Close the water checklist's cheap rows.** Report the mean molecular dipole in
`phys-bulk`; compute the dielectric constant from the dipole fluctuations; the
molecular quadrupole and polarisability anisotropy against MP2; rows 23, 29, 6, 3.
*Done when* the rows read as numbers in the checklist. Rows 7, 10, 11, 12, 18, 28
are then run as the evidence directs. *Nothing here adds a term to the law.*

**S1. The four boiling points.** `phys-coexist` for water, methane, ammonia and
methanol with the working law. Ammonia and methanol need pair data and a law:
the E8 recipe (cloud from the monomer, MP2 pairs, biased draw, wall, slopes, the
fit conditions). Methane has a law; water has several. *Done when* a table gives
`T_b`, `dH_vap`, `rho_L(T_b)` and the coexistence curve for each, with the
sampling error stated; E8's own "done when" (the 45) stays at the end of Phase 6
as decided.

**S2. Cross terms and the transferability test (§4.4)** for water with each of the
other three. *Done when* the three outcomes are reported with their errors and the
owner has seen which branch each took.

**S3. Chemical potentials.** Thermodynamic integration of a molecule in a pure
liquid, against `phys-coexist`'s `P_sat` for the same law. *Done when* the two
routes to the pure vapour pressure agree within their errors. This is the check on
the machinery all mixing depends on.

**S4. Binary liquid–vapour and liquid–liquid.** Water + methanol first (the
mildest), then water + ammonia, water + methane. **S4b (after S4, the owner's
decision 1):** water + ethanol at x = 0.10, 0.50 and the azeotrope's composition. *Done when* each gives `gamma_i(x)`
satisfying Gibbs–Duhem within error, the bubble and dew curves, and miscibility or
the gap, and one composition reproduced by direct mixed-slab coexistence.

**S5. Ions.** Records for Na⁺ and Cl⁻ (grown bases), ion–water pairs, trimers and
shells at MP2, the ion–ion pair in a water shell, a law fitted under the E8
conditions, `phys-bulk` with a net-neutral mixture. *Done when* single-ion
hydration structure (Na–O, Cl–O peaks and coordination) and the ion pair's
potential of mean force are reported with the law's pair error on unseen ion–water
pairs.

**S6. Saline at 1 % by mass.** 321 / 642 / 1284 waters with 1 / 2 / 4 ion pairs;
density, apparent molar volume, structure; excess chemical potential of the ion pair
by TI and the water activity by Gibbs–Duhem; boiling elevation and freezing
depression from the activity. *Done when* each result is given with its resolution
and compared with its reference.

**S7. Solid phases and the full diagram (E11 for these substances).** Crystal free
energies for ice Ih, methane I, ammonia I and methanol; melting curves, triple
points and sublimation by Clapeyron; then the solution freezing curve (§5.5).
*Done when* the stored rule for each pure substance reproduces its own coexistence
points and the triple point is where the curves cross.

**S8. Buckets and the stored rule.** The ratio bucket width, measured as the owner
decided: the widest step in composition over which the derived excess terms change
by less than the tolerance, on the representative mixtures above, then a three-
component mixture. Store the quick rule (§4.5) per bucket with its derivation range
and a thorough derivation running behind it. *Done when* a world can ask a node for
its mixture's boiling point and gets the rule, flagged when outside its range.

## 8. Questions the plan does not answer — the owner's

These came up while designing this; none has been decided. They are listed rather
than chosen between, per `CLAUDE.md`.

1. **Which alcohol. DECIDED (the owner, 2026-10-18): methanol first, ethanol
   after** (option C). Methanol is E8's substance and goes through S1-S4 as
   planned. Ethanol is added once the mixing law works, as the azeotrope test.
   **And the ethanol-water solution is tested at two ratios besides the
   azeotrope's.** Read as: the azeotrope composition (about 95.6 wt %, mole fraction
   about 0.89; to be verified from sources) *plus* two other compositions. *The two
   are the author's choice inside the item and the owner's to overrule:* mole
   fraction of ethanol 0.10 (about 22 wt %, a water-rich solution where water's
   network is still mostly intact) and 0.50 (about 72 wt %, the middle of the range).
   Together with the azeotrope they sample both sides of the boiling curve's
   extremum; the pure ends are already tested (§6.1). In a 216-molecule box 0.10 is 22
   ethanol + 194 water and 0.50 is 108 + 108; the azeotrope's 0.89 is 192 + 24.
   Ethanol needs everything methanol does (growth, MP2 pairs, fit) plus a second
   rotor, so its torsions follow the treatment (2) result from methanol.
2. **Solutions beyond the one asked for.** The request names only the saline
   solution as a solution; the binaries in §6.2 are included because a general
   mixing law cannot be judged on one electrolyte. Confirm the binaries are in
   scope, and which.
3. **Solid phases in this item.** E11 puts the crystal free energies after the liquid
   and the owner placed it at the end; the phase diagram asked for here includes them. S7
   schedules them last, on E11's order. Say if they should come earlier.
4. **Hydrates and clathrates** (ammonia hydrates; methane clathrate at several MPa)
   are phases of these mixtures. Not in the plan. In or out?
5. **The ice structure.** Ice Ih is proton-disordered; the choice of how to sample or
   enumerate the disorder, or whether to use ice XI or another ordered proxy for the
   free energy, is a modelling decision the plan does not make.
6. **What the law owes a mixture.** If the transferability test (§4.4) falls in the
   third branch, a term is missing from the form — a `PHYSICS.md`-weight addition.
   This document commits to reporting that and not to adding it.
7. **The tolerance for a derived phase diagram.** The plan gives 0.1 % on a
   geometry and the 25 % boiling-point bound in a failing test; it states no
   tolerance for an excess property or a miscibility gap. Boiling points are the
   measure the owner chose for the basis tolerance; the same use is proposed here
   (a looser tolerance that leaves the diagram unchanged is enough), to be agreed.
8. **Pressure-dependence.** E11 says the range is determined by the world, not fixed
   in advance. Whether a solution's diagram is first derived at 1 atm only, or over
   a fixed set of pressures, is the owner's.
9. **Quantum nuclei in the diagram.** A classical rigid liquid's heat of vaporisation
   and its melting point are off by zero-point terms (about 2.8 kcal/mol in ice). Whether
   the stored rule includes that correction (a path-integral run or a harmonic
   correction) is not decided; E8c plans the path-integral comparison only for the
   water liquid.
10. **Simulation speed.** Free-energy runs need nanoseconds (§9). Porting the Ewald
    reciprocal sum to the GPU in double precision, in the manner of the induction and
    pair loops, is the large enabling item; it is engineering for shorter compute, which
    the owner has said is acceptable, and must keep results to the tolerance the CPU
    path holds. Confirm.

## 9. Cost, from what has been measured

Measured: a 216-molecule Ewald liquid takes 72 ms per force evaluation on the CPU
(E8d §11); 30 ps took 70–95 minutes in the induction scan with another job
running; the GPU induction solve is 8.8 ms and the pair loop 8.7 ms, 22.3 ms a
step for the flexible liquid. Estimates, flagged as such:

- one nanosecond of the rigid Ewald liquid at 216 molecules is on the order of 40 hours
  on the present code, so a thermodynamic integration of twelve windows at 1 ns each
  is weeks, and the 1284-molecule saline run is several times that;
- the reciprocal-space sum is the next thing to move to the GPU after the induction
  and pair loops, and 1,284 molecules scale roughly with the number of waves;
- ammonia, methanol and the ions need pair data first at about 80–110 s per MP2 pair on
  the GPU: 400 pairs a species-pair is 9–12 hours, and the ion–water clusters cost
  more.

None of this is a reason to drop a stage; it is the reason S0 and the speed item in §8
come first.

## 10. What this document does not do

It does not add a term to the law, choose a functional, change a default, or
schedule anything for later without saying so. The checklist rows marked **A** and **?**
are the places the owner may decide a term is needed, and the §8 list is where the plan
and this design disagree about what is in scope.
