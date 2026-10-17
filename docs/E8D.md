# E8d: where the liquid's pressure comes from

The detailed record behind the E8d entry in `PLAY.md` section 7: the many-body
probe, the pressure cross-check, the literature, the designs tried and the runs,
as they were made (2026-10-14 to 16). Sections 1-6 are the probe, the pressure
budget, the refits and the cut-off study; section 7 is the Ewald sum. Numbers are
as measured at the time they were written; where a run was still going the text
says so, and the table at the end of section 7 is the place the finished numbers
go. Everything is reproducible from the tools named here (`phys-cluster-gpu
--manybody`, `phys-manybody`, `phys-virial`, `phys-ewald`, `phys-bulk` with
`PHYS_EWALD`, `phys-fit --bm`) and the scripts under `runs/many`, `runs/virial`
and `runs/e8d`; the data are not committed (they are 60 MB and regenerable).

# E8d — many-body probe, pressure cross-check, literature, designs (running log)

Written 2026-10-15 (AEST, UTC+10). Working notes, **not** a PLAY.md entry: nothing here
replaces a decided method. It is evidence for the owner's decision on E8c's open
question (the −1.9 / −2.9 kbar pressure). Everything below is reproducible from
`runs/many/` and `runs/virial/`.

## 1. What was measured (the many-body probe)

Clusters of the N nearest molecules cut from `bulk-water-298.snap`; non-additivity =
cluster interaction less the sum of its counterpoise-corrected pairs. Quantum: Hartree-
Fock only (`phys-cluster-gpu --manybody --method hf`), pairs + triples (three-body
expansion), because a direct 4-mer's first field alone was 2348 s and spilled to disk.
Classical: the fitted law with induced dipoles, every dipole solved together
(`phys-manybody`, law `runs/e8c/law-bias.txt`, 216 centres).

| n | centres | HF three-body, kcal/mol (mean ± s.e.) | law many-body, same centres | law many-body, 216 centres |
|---|---|---|---|---|
| 3 | 10 | −0.123 ± 0.079 | −0.075 | −0.016 ± 0.022 |
| 4 | 10 | −0.127 ± 0.060 | +0.020 | +0.009 ± 0.028 |
| 5 | 10 | −0.087 ± 0.133 | +0.151 | −0.003 ± 0.031 |
| 6 | 10 | −0.352 ± 0.181 | −0.018 | −0.056 ± 0.043 |
| 7 | 1 (c0) | −0.156 | +0.252 | −0.013 ± 0.052 |
| 8 | 1 (c0) | +0.185 | +0.475 | +0.027 ± 0.060 |
| 9 | 1 (c0) | +0.330 | +1.016 | +0.003 ± 0.069 |
| 10 | 1 (c0) | +1.113 | +1.306 | −0.055 ± 0.079 |

Centre 0 at n = 10: 120 triples, sum +1.11 = +3.95 repulsive − 2.84 attractive, each
triple sd 0.096. The big contributors are (0,5,9) +0.48, (0,4,6) −0.37, (0,4,5) −0.29.
Validation: the HF trimer reproduces the stored MP2 file's HF numbers (centre 0: pair sum
1.4621, cluster 1.4556, 3-body −0.0065 kcal/mol).

Reading: small (a few % of the pair sum) at n ≤ 6, attractive on average, then the one
cluster that was extended turns repulsive as the shell fills. NOT a trend from one cluster.

## 2. Pressure cross-check (what the three-body energy implies for P)

Conversions: 1 kcal/mol/Å^3 = 69,477 bar; ρ = 216 / 35.20^3 bohr^3 = 0.03343 Å^-3.

(a) **What a three-body term would have to be.** For a short-ranged three-body energy per
molecule u3 ∝ ρ², P3 = ρ² ∂u3/∂ρ = 2 ρ u3. To cancel −1872 bar: u3 = +0.403 kcal/mol per
molecule (+0.61 for −2868). That is 3.6 % (5.5 %) of the −11.2 binding energy per molecule.

(b) **What the probe measured** (cluster energy per molecule, crude: assumes the cluster's
triples per molecule stand for the liquid's, the ρ² scaling, and no surface effects):
- n = 6, ten centres: −0.352/6 = −0.059 ± 0.030 kcal/mol/molecule → P3 = **−270 ± 140 bar**.
  Relative to the law's own induction (which is already in the liquid): −0.334/6 → −260 bar.
- centre 0, n = 10: +1.113/10 = +0.111 → **+510 bar** (one cluster, n=8: +0.023→+110;
  n=9: +0.037 → +170 bar).
So the measured sign is wrong or unresolved for n ≤ 6 and right but 4x too small at n = 10.
The growth of centre 0 with n is the interesting thing: the repulsion lives in the
non-hydrogen-bonded close contacts, which a 6-molecule shell has few of and a liquid with 9
neighbours inside 4.15 Å (see §3) has many. This is a hypothesis the probe suggests, not a
measurement of the liquid.

(c) **Direct virial check** (running, `runs/virial/`): the HF three-body energy of the 10
triples of centre 0 / n = 5 at molecule-centre separations scaled by 0.99 and 1.01
(`--scale`), giving dE3/dlnλ without assuming the ρ² scaling. Compare with the same
finite difference for the law's induction (cheap, classical). → results in §6 when done.

(d) **Things the simulator leaves out, by arithmetic** (law `runs/e8c/avg/refit1/law-opt1.txt`,
molecule-molecule C6 = 47.66 hartree·bohr^6, nearly all on O–H; C8 ≈ 0):
- **No dispersion tail beyond the cutoff** (`bulk.rs`: r_cut = min(0.49 box, 16) bohr,
  switch from 14; no correction anywhere in `liquid.rs`).
  **Corrected (22:3x AEST):** the sharp-cutoff formula −(4π/3)ρ²C6/r_c³ (−350…−525 bar)
  overstates it, because the force is switched (quintic, 14→16 bohr) and the switch carries
  its own part of the integral. With the actual switch: P_tail = −2π ρ² C6 ∫(1−S) r⁻⁴ dr,
  ∫ = 9.91e-5 bohr⁻³: **−214 bar for the wall-fit law** (C6_mol-mol 47.66, U −0.092
  kcal/mol/molecule), **−117 bar for the no-wall law** (C6 25.95). So the wall-fit law's
  pressure with its own tail is ≈ −1872 − 214 = **−2086 bar** (no-wall: −1007 − 117 =
  −1124). Derived from the law's own C6, not a free parameter; the sign is the wrong
  way to help.
- Electrostatics: molecular-centre switch at 14–16 bohr, no Ewald, neutral molecules. Effect
  not yet measured (design A below).
- Nuclear quantum effects: MB-pol classical vs quantum density at 298 K is 1.004 vs
  1.001 g/cm³ (Reddy et al. 2016), i.e. ≈ 0.3 % ≈ 70 bar by 4.5e-5/bar. Cannot be the
  explanation of a −2 kbar error.
- Basis shortfall (pair energies ~5 % short, E8a): the wrong direction. Scaling pair
  energies by 1.05 made the pressure −41 → −613 bar (PLAY.md E8b (h)).

## 3. The liquid these numbers come from is a collapsed one

E8c, 1 g/cm³, wall in: first minimum of g(OO) at 4.15 Å with **8.99–9.02 neighbours**
(water: ~4.4–4.5 inside 3.3–3.5 Å), q = 0.51–0.56 (TIP4P-like control 0.72), 3.4 H-bonds a
molecule. Bukowski et al. (Science 2007) found the same thing for the ab initio pair-only
CC-pol: coordination number 11.9 until three-body terms were added. Our pair-only law is
a pair-only liquid of the same disease, milder.

## 4. Literature (what was found; sources at the end)

1. **Pair-only ab initio water collapses.** CC-pol pair potential alone: coordination 11.9;
   the three-body terms restore tetrahedral coordination (Bukowski et al., Science 315,
   1249, 2007). "The neglect of accurate short-range 3-body interactions invariably leads to
   an incorrect description of the liquid structure" (Babin, Medders, Paesani, HBB2-pol,
   2012), which adds an explicit short-range three-body exchange-repulsion/charge-transfer
   term to induction.
2. **Size of the three-body term.** ~20 % of the interaction energy in water clusters
   (MB-pol; Reddy et al. 2016); 15–20 % in cyclic clusters. MB-pol's 3B polynomials were fit
   to 12,347 CCSD(T)/CBS trimers; 2B to 42,508 dimers. "Local structure is determined by a
   delicate balance between 2-body and 3-body interactions."
3. **Sign structure** (Gillan, Alfè, Bartók, Csányi 2013): beyond-2-body errors in PBE come
   from "the close approach of monomers that are not H-bonded to each other" (ice VIII,
   compact non-amers), while B2B overbinding in BLYP comes from H-bond cooperativity.
   Both effects are present; they have opposite sign, which is why compact nearest-N
   clusters average to something small and a liquid does not.
4. **Quantum nuclei** shift the density by ~0.3 % for MB-pol (1.004 → 1.001, Reddy et al.); they
   broaden g(r) and weaken H-bonds (Wang, Ceriotti, Markland 2014). By the liquid's compressibility
   0.3 % is ~70 bar: my inference from the density shift, not a statement of those papers.
5. **MP2 for the liquid**: MP2/TZ-class baselines over-structure water (first peak ~11 %
   high); results depend strongly on basis convergence (2026 reassessment: "MP2 with a
   triple-ζ basis demonstrates poor performance indicative of insufficient convergence"; an
   MP2/RPA isobaric-isothermal Monte Carlo study by Del Ben, VandeVondele and Hutter reports RPA
   0.994 g/cm³ and an over-structured MP2). The engine's basis is derived and ~5 % short in the dimer.
6. **Long-range dispersion/cutoff** artefacts bias pressure and density (the study of Ewald-summed
   dispersion, arXiv:0704.0474, finds them secondary to Coulomb in water). The dispersion tail from the
   law's own C6 is −117 to −214 bar (section 2d corrects the first, sharp-cut-off estimate).

## 5. Hypotheses, ranked by what the evidence supports

H1. A missing **short-range three-body exchange repulsion** among non-H-bonded close
    neighbours (a pair-only law over-compacts; literature 1–3; our probe's centre-0 growth).
    Needed size: +0.4 kcal/mol/molecule.
H2. **Mid-range pair force 20–25 % too strong at 3.4–4.7 Å** (E8b: 16 neighbours × 0.045
    kcal/mol/bohr ≈ 2000 bar); the slopes-fit cured it for the rigid law (−448 bar) and the
    wall data re-opened it (−1872). Pair-form limit.
H3. **Truncation errors** in the simulator: tail (−350…−525 bar, wrong direction to help),
    electrostatic cutoff (unmeasured).
H4. Quantum nuclei (≲ 100 bar), basis shortfall (wrong sign): ruled out as the main cause.
H5. Pair-law holes inside 1.0 Å, rigid monomer, no flexibility cost: second order.

## 6. Run log

(see the end of this file; appended as runs finish)

### 6.1 Per-triple analysis (2026-10-15 22:3x AEST) — `runs/many/triples.py`, `law-triples.txt`

300 triples (10 centres × 20 at n≤6, + the 100 extra of centre 0) with both the HF three-body
energy and the law's induction three-body energy for the *same triple*:
- HF3B mean −0.008, sd 0.147; law3B mean +0.001, sd 0.171; residual mean −0.009, sd 0.115.
  corr(HF, law) = 0.75, slope 0.64: the law's induction reproduces ~56 % of the variance.
- By the closest O–O in the triple: <2.90 Å (n=75) HF −0.015, resid sd 0.190; 2.90–3.06 Å
  HF −0.020; **3.06–3.39 Å HF −0.002 (resid sd 0.059); >3.39 Å HF +0.003, sd 0.006.**
  All the three-body energy is in triples that contain an H-bonded pair. Triples of
  non-H-bonded close neighbours carry ≈ 0 at the Hartree-Fock level.
- Product-of-exponentials (three-body exchange form S_ab S_bc S_ca) fitted to the HF
  three-body: R² 0.01–0.02, to the residual R² 0.03–0.05, **negative (attractive) sign**, for
  every decay length 0.3–1.0 Å. No repulsive core is visible in the data.
- Needed for P = 0: ≈ +0.018 kcal/mol per active triple (assuming ~22 active triples per
  molecule); measured −0.008 ± 0.0085 → H1 (a missing three-body exchange repulsion) is
  disfavoured at ~3σ *in this structure*. Caveat: the structure is the over-compact one;
  the non-H-bonded contacts that the literature blames sit at 2.5–3.2 Å and are rare here.

### 6.2 The first-shell force: where the wall data moved the pressure (22:3x AEST)

`phys-diag --slopes` (law's slope of pair energy along the line of centres, minus MP2's,
kcal/mol per bohr; 66 MP2 pairs; in-sample), run in `runs/e8c/avg/`:

| O–O (Å) | n | MP2 mean | no-wall law − MP2 | wall-fit law (law-opt1) − MP2 |
|---|---|---|---|---|
| 0–2.8 | 7 | −8.65 | +0.02 (rms 0.59) | −0.54 (rms 1.24) |
| 2.8–3.1 | 17 | −1.22 | −0.00 (0.24) | **−0.27 (0.54)** |
| 3.1–3.4 | 9 | −0.22 | +0.05 | +0.04 |
| 3.4–3.8 | 22 | +0.21 | −0.01 | +0.00 |
| 3.8–4.7 | 11 | +0.28 | +0.02 | +0.04 |

The wall data cost the law 0.27 kcal/mol/bohr in the **first shell** (2.8–3.1 Å, where a
water molecule has ~4 neighbours). E8b's virial arithmetic (neighbours × ½ r × slope error
/ 3V) turns that into ≈ 4 × ½ × 5.6 bohr × 0.27 / (3 × 202 bohr³) = 0.0050 kcal/mol/bohr³
≈ **2,300 bar**, the size of the −1007 → −1872 move. The same data moved the law's
dispersion from C6 25.9 to 47.7 hartree·bohr⁶ (physical water–water ≈ 46, E6 derived 52.3).
Cause: the repulsion in these laws is **entirely** the Gaussian cloud overlap
K exp(−r²/(σa²+σb²)) with widths fixed by the electrostatic clouds (every Born–Mayer
amplitude A is 1e-300 in the fitted file): one fixed shape has to serve 0.6–2.3 Å (the
wall) and 2.8–3.5 Å (contact). It cannot, so the fit trades them.

**Variants on the existing data** (`runs/e8d/fit/variant.sh`; 20 s a fit; PBE-exchange
form; judged on reference data only — never on the liquid's pressure):

| variant | held-out RMS (≤3 kcal/mol of lowest) | slope err 2.8–3.1 Å | 3.1–3.4 | 3.4–3.8 | 3.8–4.7 | inner 0–2.8 | C6 total |
|---|---|---|---|---|---|---|---|
| no wall (control) | 0.178 | −0.00 | +0.05 | −0.01 | +0.02 | +0.02 | 6.5* |
| wall, slopeW 3, wallW 0.05 | 0.304 | −0.32 | +0.11 | +0.02 | +0.04 | −0.99 | 9.8* |
| wall, wallW 0.02 | 0.235 | −0.23 | +0.08 | −0.01 | +0.01 | −0.70 | 7.5* |
| wall, wallW 0.01 | 0.207 | −0.17 | +0.07 | −0.02 | −0.00 | −0.50 | 6.5* |
| wall, slopeW 10, wallW 0.05 | 0.273 | −0.24 | +0.09 | +0.00 | +0.03 | −0.71 | 8.9* |
| wall, slopeW 30, wallW 0.02 | 0.201 | −0.12 | +0.04 | −0.01 | +0.02 | −0.19 | 7.2* |

(* sum of the pair-table C6 entries only, not molecule–molecule; these fits start from the
default, not from the E8c start, so they are NOT the E8c laws — `law-opt1.txt` was not
reproducible from the recorded command and I have not found what produced it.)
Reading: there is a smooth trade-off; the best of these (slopeW 30 / wallW 0.02) halves the
first-shell error without losing the wall. It does not remove it: a fixed-shape repulsion
cannot, which points to a more flexible repulsion (a second overlap shape) as the
form-level remedy, and to longer-range force data for the dispersion.

### 6.3 Run queue (set 22:2x AEST)
- GPU: `runs/virial/run.sh` (finite-difference virial of HF three-body, c0/n5, λ = 0.99,
  1.01) → `runs/e8d/g1/run.sh` (100 MP2 pairs drawn from the wall-law liquid, population-
  weighted; then their MP2 radial slopes to 5.2 Å).
- CPU: `runs/e8d/scan/run.sh` (density scan of law-opt1 at 0.95, 1.05, 1.10 g/cm³, 60 ps)
  → `runs/e8d/liquids.sh` (liquids at 1.0 g/cm³ of the refit variants w001, s30w002, w002).

### 6.4 Freeing Born–Mayer beside the overlap (22:3x AEST) — `phys-fit --bm`, no code change

The fit already has a `--bm` switch that adds Born–Mayer exponentials A exp(−B r) per atom
pair beside the Gaussian cloud overlaps (the E8c fits ran with it off: `no_born_mayer`).
That gives the repulsion two shapes, which is what 6.2 said was missing. Fits on the same
existing data (20 s each), PBE-exchange form, scored on reference data only:

| variant (slopeW, wallW, --bm) | held-out RMS | slope err 0–2.8 / 2.8–3.1 / 3.8–4.7 Å | C6 mol–mol (hartree·bohr⁶) |
|---|---|---|---|
| E8c-style, no bm (30, 0.02) | 0.201 | −0.19 / −0.12 / +0.02 | (P_tail basis 47.7 for law-opt1) |
| bm (3, 0.02) | 0.200 | +0.24 / −0.22 / +0.03 | — |
| bm (30, 0.05) | 0.204 | +0.19 / −0.17 / +0.04 | 57.3 |
| bm (30, 0.02) | 0.180 | +0.14 / −0.13 / +0.05 | 47.7 |
| **bm (30, 0.01)** | **0.166** | **+0.10 / −0.11 / +0.04** | **50.4** |
| **bm (100, 0.02)** | **0.170** | **+0.09 / −0.10 / +0.05** | **51.1** |
| bm (100, 0.01) | 0.331 | degenerate (RMS 26, C6 98) — discarded | 98.3 |
(no-wall control: 0.178; inner/first-shell slope errors ≈ 0.)

The dispersion lands at 48–51 hartree·bohr⁶ molecule–molecule with no C6 pinned — E6's
derived value is 52.3 and the physical water–water C6 is ≈ 46. The no-bm wall fits had
26–48 depending on the weights.

**Hole scan** (400,000 random dimers, `tests/zz_scan.rs`, lowest energy by nearest-atom
distance 0.6–0.8 / 0.8–1.0 / 1.0–1.2 / 1.2–1.4 / 1.4–1.6 / 1.6–1.8 / 1.8–2.0 / 2.0+ Å,
kcal/mol):
- E8c wall-fit law-opt1: −297.8, −88.9, +11.6, +3.7, −1.1, −3.6, −4.7, −4.8 (holes inside 1 Å)
- bm (30, 0.01): +34.9, +18.9, +8.7, +1.8, −1.7, −3.8, −4.7, −4.6
- bm (100, 0.02): +32.8, +18.0, +8.3, +1.6, −1.9, −3.9, −4.7, −4.6
- deepest dimer −4.65 / −4.72 kcal/mol at O–O 2.92 Å (MP2 −5.0 in the engine's basis).
**E8c's open limit ("holes inside 1.0 Å … a short-range core term the law does not have")
is met by two shapes in the repulsion that the fit already supports.** Their liquids are the
test of whether the contact-shell force is also right (queue: `runs/e8d/liquids.sh`).

### 6.5 Direct three-body virial (23:2x AEST) — `runs/virial/`, `analyse.py`

Finite difference of the three-body energy under isotropic scaling of the molecule centres
(each molecule keeping its shape and orientation), λ = 0.99 against 1.00, centre 0, n = 5
(10 triples, HF; a one-sided difference — the λ = 1.01 run was dropped to free the GPU).
P3 = −ρ/(3n) · ΣdE3/dlnλ, ρ = 0.03343 Å⁻³, 1 kcal/mol/Å³ = 69,477 bar:

- **HF three-body: ΣdE3/dlnλ = +0.93 kcal/mol → P3 = −144 bar** for this cluster.
- law induction (same triples): −5.29 kcal/mol → **+819 bar**. Dominated by triple (0,3,4):
  E3 HF +0.196, law +0.383; dE/dlnλ HF −1.40, law −3.79 — the law's induction is ~2×
  the reference in energy and ~2.7× in its derivative on the one triple that matters.
- the ρ² scaling estimate (§2b) for the same cluster: +109 bar. Same order, wrong sign:
  the scaling assumption is only good for order of magnitude.
All three are O(100) bar against the −1,900 … −2,900 bar to be explained. The three-body
energy is not where the pressure is. (One cluster, one sign of strain: a cross-check, not
a measurement of the liquid's three-body pressure.)

### 6.6 The pressure budget, and where it is sensitive (23:3x AEST) — `phys-virial`

`phys-virial law snapshot` takes every pair within the cut-off, moves one molecule along the
line of centres, and sums −d·dU/ds over pairs (/3V); then it adds the many-body induction by
finite-differencing the induction solve (all dipoles together, same links and switches as
`Liquid::add_induction`) under uniform scaling of the molecule centres. **Validation:** on
the final snapshot of the 80 ps wall-law liquid it gives −1766 bar against the run's
−1872 ± 73, and U/N = −8.660 (pairs) − 2.43 (non-additive induction) = −11.09 against
−11.18.

| term (wall-fit law, 1 snapshot) | pressure, bar |
|---|---|
| translational kinetic N kT/V | +1,375 |
| repulsion (pair) | +118,233 |
| permanent electrostatics (pair) | −64,994 |
| dispersion (pair) | −27,595 |
| induction, pair sum | −14,847 |
| induction, many-body beyond the pairs | **−13,939** |
| total | **−1,766** (with the tail beyond 16 bohr: −1,980) |

By shell: the first shell (2.65–3.2 Å, ~400 pairs) carries +18,000 of the virial net; shells
beyond 3.7 Å carry −9,000. Nothing here is small: **the pressure is a ~1-in-15 residual of
terms of order 10⁴ bar**, so a 7 % error in any one of them is the whole error.

**Non-additive induction in the liquid is −2.43 kcal/mol/molecule** (many-body solve −4.89 vs
pair sum −2.46), 22 % of U/N. The 10-molecule clusters of §1 see ≈ 0 of it — non-additive
induction is long-ranged, so a compact cluster is not a sample of the liquid for it. Which
means the cluster probe could not have found where the law's non-additivity matters, and
the probe's per-triple comparison (below) is the only handle on it.

**Sensitivity** (same snapshot, α scaled, nothing refit): α × 1.00 → −1,766; × 0.95 → −706;
× 0.90 → +431; × 0.80 → +2,937. **Each 5 % of polarisability is ≈ 1,100 bar; the whole
error is an ~8 % polarisability error.** The law's molecular α from its sites is 10.05
bohr³ (MP2 finite field 9.908; +1.4 %), so the *gas-phase* polarisability is right; what
the reference data can say is how much of it acts at liquid contacts.

**What the reference says about the law's induction** (300 HF triples, `ind_scan.py`):
law three-body vs HF three-body per triple — RMS 0.115 / slope 0.64 at α × 1.00,
0.102 / 0.75 at × 0.90, 0.093 / 0.87 at × 0.80, **0.089 / 1.04 at × 0.70**. The triples want
the law's induction ~30 % weaker (α × 0.70), i.e. the
law's induction is stronger than the reference at the contacts the liquid has. Exchange-
induction quenching and charge penetration are the physics point-dipole induction lacks.
**Pair data cannot see it**: refit with α × 1.0 / 0.9 / 0.8 / 0.7 (Born–Mayer + overlap,
slopeW 30, wallW 0.01), held-out pair error 0.166 / 0.164 / 0.164 / 0.164 — the repulsion
absorbs the change. So the triples are the only discriminating reference, and they point
the right way (weaker induction → pressure up), but by more than the pressure needs (0.92
for P = 0, 0.70 from triples), without a refit of the law around it.

**Instantaneous trap:** the same snapshot gives −8,899 bar with the Born–Mayer refit
(α × 1.0) against −1,766 for the E8c law, though the two have the same held-out pair error
(0.166 vs 0.304) — pair accuracy to 0.2–0.3 kcal/mol does not pin the virial, a repeat of
E8b's conclusion. The Born–Mayer refit's own liquid (80 ps, running to 25 %) reads
−2,500 … −3,100 bar, U/N −11.7: a **more** collapsed liquid than E8c's. Not a pressure fix
on its own; it is a hole fix.

### 6.7 Screening liquids queued (23:3x AEST)
`runs/e8d/liquids2.sh`: 40 ps at 1.0 g/cm³, Born–Mayer refits with α × 0.8, 0.7, 0.9
(`runs/e8d/fit/alpha_*`). G1 (`runs/e8d/g1/`): 100 MP2 pairs from the wall-law liquid
(78 s each) then their radial slopes — the long-range shells' reference.

### 6.8 Truncation: the pressure is not converged in the cut-off (23:4x AEST) — the main finding so far

`phys-bulk` cuts all pair terms off at 16 bohr with a quintic switch on the molecules' centres
(box 35.2 bohr, so 17.6 is the most a 216-molecule box can offer). `bulk.rs`'s own header
calls this "a known error of this prototype". It had not been measured. With `phys-virial`
(pair virial + many-body induction virial, dispersion tail included), the **same** snapshot
(final frame of the 80 ps wall-law liquid), the law fixed, only the cut-off changed:

| r_cut (bohr), switch from r_cut−2 | box | pressure, bar (incl. tail) | non-additive induction, kcal/mol/molecule |
|---|---|---|---|
| 10 | 1×216 | +2,685 | −1.31 |
| 12 | 1×216 | −3,015 | −2.01 |
| 14 | 1×216 | −4,433 | −2.29 |
| **16 (the run)** | 1×216 | **−1,980** | −2.43 |
| 17.5 (half box) | 1×216 | +499 | −2.28 |
| 20 | 2×2×2 replicated (1,728) | +778 | −2.11 |
| 24 | 2×2×2 | +460 | −1.98 |
| 28 | 2×2×2 | −1,167 | −2.06 |
| 32 | 2×2×2 | −424 | −2.12 |

(The replication — the same configuration repeated periodically, `wall-liquid-x8.snap` —
reproduces the 16-bohr value exactly, −1,766 without the tail.) The pressure **oscillates by
±1–2 kbar with the cut-off and is still doing so at 32 bohr**; the mean over 20–32 is ≈ −90 bar.
So: (i) the −1,872 bar of the E8c table is, for this configuration, an artefact of the
truncation to within the scatter of the truncation itself; (ii) a switched cut-off cannot
resolve this pressure to better than ~±1 kbar at any cut-off a box this size allows;
(iii) every pressure in E8b/E8c (and all my screening liquids) carries that error, and the
structure they generated (9 neighbours, q 0.5) may be shaped by it as well. The remedy is
the standard one, long-range sums (Ewald for the Gaussian charges and the induced dipoles),
and it is general: it applies to every polar molecule of the class, which no choice of
pair-law numbers does.

Caveats: one configuration; its periodic images repeat the 216-molecule structure, so the
sum is the periodic crystal's, not an infinite liquid's; the dipolar sums are conditionally
convergent (surface term, the cell's net dipole), which is part of why a spherical
truncation oscillates. This is evidence for the cause, not yet the fix.

### 6.9 Queue (23:4x AEST)
CPU (BelowNormal): `alpha_0.80` (running) → `rc/17.5` and `rc/12` (the E8c wall-fit law at
other cut-offs, 40 ps: does the *structure* depend on the truncation?) → `alpha_0.70`,
`alpha_0.90`. GPU: G1 pairs (78 s each) then slopes.

**Second configuration, second law (no-wall law, `runs/e8c/avg/law-final.txt`, final frame of
its own 100 ps liquid, MD average −1,007 ± 116 bar), same scan, with tail:** 16 → −2,583;
20 → +1,419; 24 → +634; 28 → −722; 32 → −57 bar (mean of 20–32: +319). The same pattern: the
16-bohr value is ~2.9 kbar below the mean of the long cut-offs, and the long cut-offs
oscillate by ±1 kbar. (Single frames: the 16-bohr values, −2,583 and −1,980, are within a
frame's scatter of the runs' averages.) `runs/e8d/trunc.log`, `trunc2.log`.

## 7. Ewald (built 16 Oct, the owner's instruction: "Try the Ewald method, and review the results")

### 7.1 What was built (`src/ewald.rs`, `induction::solve_ewald`, `liquid::{EwaldLaw, add_ewald}`, `phys-ewald`)

Every object in the liquid (permanent charge, induced dipole) is a Gaussian of width w and
two of them interact through derivatives of erf(r/s)/r with s² = w_a² + w_b², so one
**common** width s′ ≥ s for all pairs splits every interaction as [kernel(s) − kernel(s′)]
(real space, nearest images, the engine's own `dd_terms`/`charge_terms` differenced) +
kernel(s′) (reciprocal space, point multipoles, 4π/k² exp(−k²s′²/4)).
`E_rec = (4π/V) Σ_half f(k)|S(k)|²`, `S = Σ q e^(ik·r) + Σ i(k·μ) e^(ik·r)`. Subtracted: each
charge's and each dipole's self term, the long-range part of the permanent charges' pairs
inside a molecule (a constant), and the long-range field of a molecule's own charges on its
own dipoles (the model excludes it; done as real-space "exclusion" items with an ordinary
pair's force). The dipoles are solved with real-space remainder and reciprocal coupling
together (CG). Forces on every charge and dipole site (→ COM force and torque), reciprocal
virial by central difference under rigid scaling. Tin-foil. Default behaviour of every
existing tool is unchanged; `PHYS_EWALD=<s′ bohr>` turns it on in `phys-bulk`.

### 7.2 Validation (`phys-ewald law state snapshot s′,… check`), E8c wall-fit law on its own
final frame, 216 molecules, r_cut 16

- Reconstructing the liquid from the snapshot reproduces the truncated run: U/N −11.092
  and P −1,782 bar (the run: −11.18, −1,872 ± 73; `phys-virial`: −1,766).
- **Independence of the splitting width** — U/N and P with Ewald: s′ = 3.4: −10.5170,
  170.0 bar; 3.8: −10.5170, 169.5; 4.2: −10.5175, 166.9 (4.8 and 5.4 drift, −10.521 / 148
  and −10.532 / 89, because erfc(16/s′) is no longer negligible at the cut-off).
- **Forces and torques vs central differences of the re-minimised energy**, 3 molecules ×
  3 axes each: agree to ~1e-8 (e.g. force +2.509788e-3 vs +2.509776e-3; torque −2.983456e-3
  vs −2.983431e-3), i.e. 1e-5 … 1e-6 relative.
- **Virial vs −dU/dln λ** under rigid scaling of box and centres: agrees to 1.8e-5 hartree
  (0.0 bar) for both the truncated and the Ewald versions.
- cost: 0.2 s per force evaluation (0.1 s truncated), threads over disjoint output slices.

### 7.3 Static Ewald pressure of the truncated runs' final frames
| law (frame) | truncated, this frame | Ewald (s′ 3.8 / 4.2) | U/N truncated → Ewald |
|---|---|---|---|
| E8c wall-fit (80 ps frame) | −1,782 | +170 / +167 bar | −11.09 → −10.52 |
| no-wall (100 ps frame) | −2,481 | +510 / +507 bar | −11.13 → −10.48 |
The truncated sum over-binds by ~0.6 kcal/mol a molecule and its pressure is 2–3 kbar too
low on these frames. But these frames were *made* by truncated dynamics.

### 7.4 Truncated dynamics is self-consistent at each cut-off (so the static scan of 6.8 is
not the whole story): E8c wall-fit law, 40 ps, 216 molecules, 1.0 g/cm³ —
r_cut 17.5: P −1,667 ± 67 bar, U/N −11.12, 9.01 neighbours, first minimum 4.15 Å
(r_cut 16, 80 ps: −1,872 ± 73, −11.18, 8.99, 4.15 Å). The structure and the pressure do not
move when the cut-off does; what moves is the pressure of a *fixed* configuration. Only
Ewald *dynamics* can say whether the structure itself is a truncation artefact. Born–Mayer
refit with α × 0.8, truncated, 40 ps: P −2,146 ± 139, U/N −10.58, 7.66 neighbours, first
peak 2.52 (so: weaker induction decollapses the liquid somewhat, and does not fix the
pressure at this truncation).

### 7.5 Running (02:1x AEST): Ewald dynamics `runs/e8d/ewald/{opt1,nowall}` (s′ = 4.2, 40 ps
each, ~2–3 h each).

### 7.6 Result: Ewald dynamics of the E8c wall-fit law (finished 05:4x AEST, 16 Oct)

216 molecules, 1.0 g/cm3, 298 K, s' = 4.2 bohr, k accuracy 1e-5 (1535 vectors; 1e-5 changes
the static pressure by 5 bar against 1e-7), r_cut 16, 30 ps with 24 ps sampled
(`runs/e8d/ewald/opt1`), against the truncated runs of the same law:

| | Ewald, 30 ps | truncated 16 bohr, 80 ps (E8c) | truncated 17.5 bohr, 40 ps |
|---|---|---|---|
| pressure, bar | **-701** (-916 with the dispersion tail of -214) | -1,872 +- 73 | -1,667 +- 67 |
| U/N, kcal/mol | -10.427 +- 0.013 (-10.519 with the tail) | -11.181 | -11.119 |
| heat of vaporisation | 11.02 (11.11 with the tail) | 11.77 | 11.71 |
| g(OO) first peak | 2.68 at 2.85 A | 2.83 at 2.85 | 2.79 at 2.85 |
| first minimum | 0.86 at 4.15 A | 0.88 at 4.15 | 0.88 at 4.15 |
| neighbours inside the minimum | **9.07** | 8.99 | 9.01 |
| tetrahedral q (one frame) | **0.541** | 0.558 | 0.547 |
| hydrogen bonds a molecule (one frame) | 3.45 | 3.43 | 3.50 |

The block standard error printed for the Ewald pressure (20 bar) is not believable: the
samples run -1,934 ... +211 with a typical successive scatter of ~1,000 bar; +- 200 is the
honest figure for 24 ps. **The truncation took 1.2 kbar off the pressure and 0.75 kcal/mol a
molecule off the binding, and took nothing off the structure.** The collapsed 9-neighbour
liquid is the law's.

## 8. The owner's idea: ball-and-stick molecules, atom-atom forces by brute force (16 Oct)

*"Treat all molecules as stick-and-ball models, a spring between atoms with the appropriate
stiffness; put many in the scene; brute-force the forces between each atom; step to steady
state; see where it falls short."* Built as `phys-sticks` (`src/bin/sticks.rs`) with the
fit in `runs/e8d/sticks/fit_sticks.py`. Everything is derived from the engine's own
reference: nothing is taken from a published force field.

### 8.1 The springs (derived, and checked)
The Hessian of the engine's monomer surface (`surface-water.txt`, HF + RI-MP2 over r1, r2,
θ) at its minimum (r_OH 0.96154 Å, θ 104.145°): H = [[0.5495, −0.0089, 0.0319], [−0.0089,
0.5495, 0.0319], [0.0319, 0.0319, 0.1582]] (hartree per bohr², bohr·rad, rad²). It gives
harmonic frequencies **1626, 3847, 3973 cm⁻¹** — the numbers `phys-vib` reports for the
same surface — so "the appropriate stiffness" is derived and agrees with the engine. The
springs sit at the pair data's averaged shape (0.97894 Å, 104.078°), so that the
intermolecular law and the monomer agree on the geometry.

### 8.2 The atom–atom law, fitted to the 410 MP2 pair energies the engine has
Point charges on the atoms (q_O = −2 q_H), an exponential repulsion and a Tang–Toennies
damped C6 for each of O–O, O–H, H–H, nothing else (no clouds, no extra site, no induction),
weighted to the low-energy pairs as the project does; every fifth pair held out.

| variant | q_H | C6 (OO, OH, HH) | held-out RMS within 3 kcal/mol of the lowest | all pairs < 20 kcal/mol |
|---|---|---|---|---|
| gas-phase charge (from the monomer's MP2 dipole 1.873 D) | 0.3238 | 133.9, 11.4, **−5.0** | 0.72 | 1.21 |
| charge fitted to the pairs | 0.2973 | 103.5, 23.0, **−8.5** | 0.72 | 1.21 |
| dispersion held at the engine's derived C6 | 0.3238 | 17.48, 6.33, 2.38 | 0.74 | 1.36 |
| (the project's cloud + induction law, for scale) | — | — | 0.24–0.30 | — |
When the C6s are free the fit takes them to 3× the physical molecule–molecule value (160
against ~50 hartree·bohr⁶) and one of them negative: dispersion is being used as a fudge
for what point charges and pair additivity cannot do. With them pinned to the derived
values the pair error is no worse. In every case the pair error is 2.5–3× the project's law.

### 8.3 A short test (64 molecules, 3 ps, gas-charge variant)
Stable, ~16 ms a step (648 → 192 atoms brute force), 2000 steps a picosecond at dt = 0.48 fs.
The liquid **collapses**: 12.8 neighbours inside the first minimum (water: 4.4), first peak
2.15 at 2.95 Å, P −2,837 bar (−4,940 with the dispersion tail), U/N −8.6 kcal/mol. Spring
energy 1.08 kcal/mol a molecule (equipartition 0.89: the classical hydrogens are hot,
T 315 K on a 298 K bath in a short run).

### 8.4 216 molecules, 20 ps: `runs/e8d/sticks/run-{e6disp,gascharge,freecharge}` (running)

**e6disp finished (216 molecules, 1.0 g/cm3, 298 K, 20 ps with 16 ps sampled, 47 min on a loaded machine; dispersion held at the engine's derived C6, gas-phase charge q_H 0.3238):**

| | ball-and-stick, e6disp | rigid + induction, Ewald (E8c law) | measured water |
|---|---|---|---|
| pressure at 1 g/cm3 | **+379 +- 147 bar** (+143 with the dispersion tail -235) | -701 (-916 with tail) | ~0 |
| intermolecular U/N | -6.73 kcal/mol (-6.83 with tail) | -10.43 (-10.52) | |
| heat of vaporisation | 7.3 (7.4 with tail) | 11.0 (11.1) | 10.5 |
| g(OO) first peak | 2.20 at 2.85 A | 2.68 at 2.85 | 2.6-2.75 at 2.80 |
| first minimum, neighbours inside | 0.84 at 4.45 A, **12.01** | 0.86 at 4.15, 9.07 | ~0.8 at ~3.4, ~4.4 |
| tetrahedral q | **0.479** | 0.541 | ~0.7 (TIP4P-like control) |
| hydrogen bonds a molecule | 2.71 | 3.45 | ~3.5-3.8 |
| spring energy a molecule | 1.13 (equipartition 0.89) | n/a | |

The pressure is right (within ~1 sigma of zero) and everything else is wrong: a dense simple liquid
(12 neighbours is what a close-packed fluid has, and what Bukowski et al. found for the
pair-only CC-pol potential, 11.9), 3.7 kcal/mol under-bound, with fewer hydrogen bonds and
a less tetrahedral arrangement than even the pair-only rigid law. **A near-zero pressure is
therefore not evidence that a liquid is water.** The hydrogens' spring energy is 27 % above
equipartition although the kinetic temperature is 297.7 K. I have not found the cause (the Langevin
thermostat's weak coupling in production, the short sampling, or an energy leak from the time step
are the candidates); it is flagged, not investigated, and does not touch the pressure's sign or the structure.

**gascharge finished (same box and length; gas-phase charge q_H 0.3238, C6 free: OO 133.9, OH 11.4, HH -5.0, molecule-molecule ~160):**
pressure **-4,777 +- 132 bar** (-5,495 with the dispersion tail of -718, which is large because the
fitted C6 is 3x physical), intermolecular U/N -9.27 kcal/mol (-9.58), heat of vaporisation 9.86 (10.17
with the tail), g(OO) first peak 2.25 at 2.95 A, first minimum 0.81 at 4.45 A with **12.17** neighbours
inside it. It lands near the measured heat of vaporisation (10.5) because the inflated dispersion
supplies the binding that point charges and pair additivity do not, and it pays for it with a
pressure of -5 kbar and the same close-packed structure. The two variants bracket the answer: the
honest dispersion (e6disp) gives the right pressure and 3.7 kcal/mol too little binding; the
fudged dispersion gives roughly the right binding and a pressure 5 kbar too low. Neither is
water's structure (12 neighbours).

### 7.7 Result: Ewald dynamics of the no-wall law (finished 08:2x AEST, 16 Oct; 30 ps, 24 sampled)

| | Ewald | truncated 16 bohr (E8c, 100 ps) |
|---|---|---|
| pressure, bar | **-21 +- 98** (-138 with the dispersion tail of -117) | -1,007 +- 116 |
| U/N, kcal/mol | -10.30 (-10.32 with tail) | -11.13 |
| g(OO) first peak | 2.56 at 2.85 A | 2.73 at 2.75 |
| first minimum, neighbours inside | 0.90 at 4.05 A, **8.12** | 0.88 at 3.95 A, 7.49 |
So for the law without the wall data, Ewald takes the pressure to zero within error, and (as for the
wall-fit law) the structure does not become more water-like: slightly more neighbours, not fewer.

## 9. Out-of-sample slopes from the liquid (G1, 83 MP2 radial slopes on 100 pairs drawn from the E8c liquid)

law - MP2 mean slope (kcal/mol/bohr), by O-O shell; the pairs were in no fit.
| O-O (A) | n | MP2 mean | E8c wall-fit | no-wall | Born-Mayer 30/0.01 | Born-Mayer 100/0.02 |
|---|---|---|---|---|---|---|
| 0-2.8 | 4 | -4.16 | -1.55 | +0.07 | +0.01 | +0.02 |
| 2.8-3.1 | 11 | -0.22 | **-0.335** | +0.03 | +0.05 | +0.08 |
| 3.1-3.4 | 8 | +0.28 | +0.22 | +0.22 | +0.19 | +0.20 |
| 3.4-3.8 | 11 | -0.47 | +0.07 | +0.11 | +0.16 | +0.14 |
| 3.8-4.7 | 29 | +0.13 | -0.005 | -0.02 | +0.005 | +0.007 |
The long-range shell (3.8-4.7, 29 pairs) is right for every law (the +0.04 suspected from 11 in-sample
pairs was noise). The wall-fit law's contact shell is off by -0.335 (confirmed out of sample). All four
laws are too attractive in 3.1-3.8 A by 0.07-0.22: by the virial arithmetic (156 pairs at 3.18-3.44 A
x 0.19, 179 at 3.44-3.70 A x 0.15) about -0.7 and -0.7 kbar. A refit that adds the even-numbered
liquid slopes (`runs/e8d/fit2/`) leaves the withheld odd ones at +0.08..+0.11 and +0.12..+0.16 in
those two shells (before: +0.08, +0.18): the isotropic pair form does not remove it, with 4-5
withheld pairs a bin.

**freecharge finished (charge fitted to the pairs, q_H 0.2973; C6 free, 103.5 / 23.0 / -8.5):** pressure
-5,388 +- 71 bar (-6,115 with the tail of -726), intermolecular U/N -9.30 (-9.61), heat of vaporisation
9.9 (10.2), g(OO) first peak 2.26 at 2.95 A, first minimum 0.80 at 4.55 A, **12.85** neighbours. Same
picture as the gas-charge fit: free dispersion buys the binding and costs the pressure.

### 7.8 Ewald dynamics of the Born-Mayer refit (bm 30/0.01; at 67 % when this was written)
Production samples read -2,050 ... -2,960 bar (mean about -2,400), U/N about -10.8. With Ewald on,
the law that closes the dimer holes and fixes the contact shell's slope has a pressure *three times
as negative* as the E8c wall-fit law's -700 and binds 0.4 kcal/mol a molecule more.

## 10. The owner's instruction: "come up with the best solution, even a novel one" (16 Oct, afternoon)

What was tried, in the order tried, with what it showed. All checked against the engine's own reference
(pair energies, the 83 out-of-sample liquid slopes, the 300 HF triples) before any liquid was run.

### 10.1 Widening the dipoles (short-range damping of induction), `PHYS_IND_WIDTH` in `induction::sigma_of`
Scaling the smeared dipoles' width by f and refitting the site polarisabilities to the same MP2 tensor
(`phys-polar --target 10.2675,10.1771,9.7044,-0.18`): f = 1.0 / 1.2 / 1.4 / 1.6 gives alpha_O / alpha_H =
9.23/0.85, 8.25/1.47, 6.34/2.47, 3.22/3.94 (polarisability moves from oxygen to hydrogen, the tensor stays
matched to f = 1.4). **Against the 300 HF triples it gets worse**: rms 0.115 -> 0.134 -> 0.161 -> 0.281 kcal/mol,
correlation 0.75 -> 0.64 -> 0.36 -> -0.21. The triples do not support short-range damping by width. (The
distribution of polarisability that Thole-damped literature models use is the f = 1.4-1.6 one; for *this*
induction model and *this* reference it fits the three-body energies worse.)

### 10.2 A fixed Born-Mayer exponent from the density's tail, `phys-fit --bm-exponent`
Exchange repulsion falls as exp(-2 kappa R), kappa = sqrt(2 IP), about 2.0 per bohr from the engine's own HF HOMO,
while free fits chose 2.7-3.4 for O-H and H-H. Holding B at 2.0 / 2.2 / 1.8 / 1.6 and fitting amplitudes:
held-out pair error 0.231 / 0.187 / 0.393 / 1.05 (free: 0.169), held-out liquid slope at 3.4-3.8 A +0.14 / +0.16 /
+0.25 / +0.49 (free: +0.157). **It does not remove the shell error** and costs pair accuracy. Not the cause.

### 10.3 Where the shared 3.1-3.8 A error comes from (`phys-slopeparts`)
On the 83 liquid pairs the law's slope splits into electrostatics, pair induction, repulsion, dispersion.
In 3.4-3.8 A: electrostatics -0.08, induction +0.10, repulsion -0.96, dispersion +0.58 (MP2 total -0.47, law
-0.31): the error (+0.16) is a small fraction of the repulsion/dispersion balance; induction is not it.
Per pair, the large errors (+0.2 ... +0.45) are on *bent contacts*: shortest H..O 2.6-3.2 A at O-H..O angles
of 95-125 deg (pairs 10, 25, 26, 44, 52, 75), where the law is not repulsive enough. The law has
essentially no O...O repulsion at 3.4 A (the Gaussian overlaps fall as exp(-r^2/s^2); the O-O Born-Mayer
term is ~1e-3 kcal/mol there).

### 10.4 Repulsion-only lone-pair sites (zero charge, Boys positions 0.564 bohr from O, +-0.499 out of plane)
Added as a site type with a start of non-zero overlap strength (the fit works in logs, a start of zero cannot
move; my first attempt had that artefact and is discarded): the fit sends every overlap involving the new
sites to zero at widths 0.8, 1.1 and 1.5 bohr. **The data do not want an out-of-plane lone-pair repulsion.**

### 10.5 Population-weighted pair fits (weighting temperature 1000 -> 100,000 K)
Held-out slope at 3.4-3.8 A goes +0.157 -> +0.108; the energy bias at 3.1-3.4 A stays -0.19 ... -0.22 kcal/mol.
The weighting is not the cause; the isotropic form cannot represent that shell.

### 10.6 What does work: a hole-free law with the best liquid-range pressure (`runs/e8d/core/make_core.py`)
The no-wall law has the best Ewald liquid (P -21 +- 98 bar, q 0.60, 3.57 hydrogen bonds) and attractive holes
inside 1 A (-577 kcal/mol at 0.6-0.8 A). The plan named the remedy ("a short-range core term the law does not
have"). Held: everything in the no-wall law; added: a steep Born-Mayer core (exponent 7.0 per bohr, the
dispersion damping exponent written out as 2.0, its implicit value) whose amplitudes are fitted by non-
negative least squares to what the law leaves on the 510 pair energies (`phys-lawenergy`). Only O-H is
needed (A_OO = A_HH = 0). It changes pairs below 10 kcal/mol by at most 0.06 kcal/mol and the wall pairs
to within 17 % relative. **Hole scan, 400,000 random dimers, lowest energy by nearest-atom distance
(0.6-0.8 ... 2.0+ A): no-wall -577.4, -260.3, +4.3, +2.6, -1.4, -3.8, -4.6, -4.6; with the core +33.5, +20.6,
+10.2, +2.6, -1.4, -3.8, -4.6, -4.6**; deepest dimer -4.63 kcal/mol at O-O 2.92 A. The liquid range is
identical by construction.

### 10.7 Ewald dynamics of the hole-free law (no-wall law + O-H core, exponent 7): finished 18:2x AEST
216 molecules, 1 g/cm3, 298 K, s' 4.2, 30 ps (24 sampled), `runs/e8d/ewald/core_B7`:
pressure **-88 +- 97 bar (-204 with the dispersion tail of -117)**, U/N -10.214 (-10.264), heat of vaporisation 10.8
(10.86; measured 10.5), g(OO) first peak **2.57 at 2.85 A** (measured 2.57-2.75 at 2.80), first minimum 0.91 at 3.95 A,
**7.56** neighbours inside it, tetrahedral q 0.572, hydrogen bonds 3.61 a molecule. Against the same law without
the core (Ewald, `ewald/nowall`): -21 +- 98 bar, 8.12 neighbours, q 0.603, 3.57. The core changes nothing in the
liquid beyond noise, as it was built not to, and removes every hole. This is the best law so far on the combination
of pressure, energy, hydrogen bonds and hole-free dimers; it still has no first-shell gap in g(OO) (minimum 0.91).

## 11. Making the Ewald liquid cheaper (16 Oct, the owner asked)

Profile of one force evaluation, 216 molecules, s' 4.2, k accuracy 1e-5 (1535 vectors), on a quiet machine
(`phys-ewald ... profile`, the other jobs suspended for the measurement):

| | before (184a556) | after |
|---|---|---|
| real-space pair loop | 10 ms | 10 ms |
| + charge reciprocal sum | 62 ms | 34 ms |
| truncated with induction (phys-bulk's old path) | 38 ms | 38 ms |
| **Ewald with induction (the production run)** | **117 ms** | **72 ms** |

Three exact changes, none altering a result beyond rounding (energy -10.4757, P 507.2 bar, force/torque
finite differences and the virial check give the same digits, `tests/ewald.rs` and 21+5+7 other tests pass):
1. **Analytic reciprocal virial** (`Ewald::virial_analytic`) in place of the central difference, which recomputed
   every site's phase factors twice more. dE/dlnL = -3E + (4pi/V) sum f (k^2 s^2/2 + 2)|S|^2 +
   (8pi/V) sum f Re[S* dS/dlnL]. It agrees with the numerical one to the printed digits
   (`PHYS_EWALD_NUMERIC_VIRIAL=1` keeps the old one as an oracle).
2. **Coincident charges merged in reciprocal space** (a core and its cloud on one atom): 7 sites a molecule -> 4;
   the structure factor, forces on the molecule and virial are unchanged; self and same-molecule terms keep the
   original charges.
3. **`erf_fast`** (7th-order Taylor about a table of nodes every 1/64): used only in the Ewald paths; agrees with
   the series to 2e-15 below x = 4 (the series itself carries ~1e-14 beyond). It did not matter for speed
   (erf was not the bottleneck) and is kept because it is exact and free.
Measured and found not worth it: relaxing the dipole solve's tolerance from 1e-10 to 1e-7 (`PHYS_IND_TOL`;
91 -> 89 ms: warm starts already keep it at 0-3 iterations).
Not done, with what each would give: a persistent thread pool (the ~15 parallel regions a step each spawn 12
threads; ~10 % at most); a time step of 80 au instead of 40 (2x, but needs an energy-drift check first and
changes the dynamics, so it is the owner's call); k accuracy 1e-4 (28 % fewer vectors, shifts the pressure ~40 bar).
**Contention is the larger factor**: with two other jobs running the same step took 85-300 ms.

## 12. A second molecule: methane (16 Oct, evening)

The recipe — derive the dispersion from the engine's own electrons, fit a site-site law to the engine's MP2 pair
energies, run the liquid — applied to methane. 80 MP2 pairs (random, `phys-pairs-gpu methane 80 --mp2`, 225 s each
on a loaded machine, ~100 s quiet); `phys-fit methane`; `phys-bulk methane` at 111 K, 0.4226 g/cm3, 216 molecules,
40 ps (32 sampled), no charges so no Ewald question. Measured liquid methane at its normal boiling point (111.7 K):
density 0.4226 g/cm3, pressure 1 atm, heat of vaporisation 8.2 kJ/mol = 1.96 kcal/mol, g(C-C) first peak ~2.5 at
~4.1 A, first minimum at ~5.6-5.8 A, ~12 neighbours.

**The engine's derived dispersion** (`phys-c6`, new: the non-local functional's C6 split among atoms by the density
partition): C-C 20.17, C-H 7.85, H-H 3.12 hartree bohr^6, molecule-molecule **132.9** (measured ~130). The same tool
on water gives 17.84 / 6.56 / 2.51, total 54.1 (E6's 17.48 / 6.33 / 2.38 at the averaged shape; physical ~46).

| | as fitted (nothing held) | dispersion held at the derived C6 | measured |
|---|---|---|---|
| molecule-molecule C6 | 33.5 (the fit put the rest in C8 = 200) | 132.9 | ~130 |
| held-out pair error, kcal/mol | 0.040 | 0.057 | |
| pressure, bar | +29 +- 8 (-5 with the tail) | -222 +- 24 (-357 with the tail of -135) | ~1 |
| U/N, kcal/mol | -1.44 (-1.47) | -1.785 (-1.908) | |
| heat of vaporisation, kcal/mol | 1.66 (1.69) | **2.01 (2.13)** | 1.96 |
| first g(C-C) peak | 2.61 at 4.05 A | **2.59 at 4.05 A** | ~2.5 at ~4.1 |
| first minimum, neighbours inside | 5.55 A, 11.5 | 5.75 A, **12.4** | ~5.7 A, ~12 |

Reading: with nothing held, the as-fitted law gets the pressure right for the wrong reason (a C6 a quarter of the
physical value, the near field carried by C8) and is 14 % under-bound. With the derived dispersion the cohesion is
right to 2-9 %, the structure is right, and the pressure at the measured density is -357 bar, i.e. the law's liquid
wants to be denser by roughly 9 % (methane's isothermal compressibility ~2.5e-4 per bar). **The recipe holds on a
second molecule where the many-body physics is small; what water adds (hydrogen bonding, strong induction) is
where it does not yet.** Not done for methane: a smaller pair-fit error with more or liquid-sampled pairs; the
7 % of the pressure that a measured-density MP2 pair-additive methane would still miss.

## 13. Density scan of the hole-free law under Ewald (finished 00:00 AEST, 17 Oct)

Law `runs/e8d/core/law-core-B7.txt` (no-wall law + O-H core), 216 molecules, 298 K, Ewald s' 4.2, 30 ps (24 sampled).

| density, g/cm3 | pressure, bar | with the dispersion tail | U/N | heat of vaporisation | g(OO) first peak | first minimum | neighbours inside it | q | H-bonds |
|---|---|---|---|---|---|---|---|---|---|
| 0.95 | -1,045 +- 41 | -1,151 | -10.146 | 10.74 | 2.67 at 2.85 A | 0.88 at **3.45 A** | **4.65** | **0.634** | 3.50 |
| 1.00 | -88 +- 97 | -204 | -10.214 | 10.81 | 2.57 at 2.85 | 0.91 at 3.95 | 7.56 | 0.572 | 3.61 |
| 1.05 | +1,301 +- 97 | +1,172 | -10.350 | 10.94 | 2.55 at 2.75 | 0.90 at 4.05 | 8.70 | 0.566 | 3.58 |
(the +- are 5-block errors; the 1.0 value is the noisier, ~+-200 realistic)

**Equilibrium density and compressibility from the pressure** (linear in density): P = 0 at 1.003 g/cm3 (no tail) to
1.007 (with the tail) against 0.997 measured; dP/drho over 0.95-1.05 is 23,200 bar per g/cm3 (with the tail), so
kappa_T = 1/(rho dP/drho) = **4.3e-5 per bar** against 4.5e-5 measured (0.95-1.00: 5.4e-5; 1.00-1.05: 3.6e-5).
**The hole-free law's bulk thermodynamics are right**: density within 1 %, compressibility within 5 %, heat of
vaporisation 10.8 against 10.5.

**The structure is not hysteretic; the real deficiency is the shell gap.** Counting neighbours inside the first
g(OO) minimum made 0.95 and 1.0 g/cm3 look like two structures (4.65 against 7.6), but the minimum is so shallow
that its location jumps. Compared at fixed radius, runs started from the open and from the dense configuration agree
to about 0.03 in g(r) (g(3.4) 0.91/0.93 at 0.95; 0.98/0.97 at 1.0), neighbours inside 3.3 A are 4.47 from either
start, and the energies differ by 0.06-0.09 kcal/mol. The law's first shell is about right; its minimum is too
shallow (g(3.4) ~0.98 against ~0.8) and its second shell too weak (g(4.4) ~0.98 against ~1.1), the signature of the
bent-contact over-attraction at 3.1-3.8 A (section 10.3).


## 14. A direction-dependent repulsion: tested, not supported

Idea (the owner's): give each atom's repulsion an antenna-style strength map, a size that varies with direction
(an anisotropic atom-atom potential), as the cure for the shared 3.1-3.8 A over-attraction on bent H..O contacts.
Test before any build, on the hole-free law (`runs/e8d/aniso_test.py`, `phys-slopeparts` with `VERBOSE=1`):

- Energies, 47 liquid-sampled pairs with O-O 2.7-4.3 A: rms error 0.18 kcal/mol against 2.12 of pair energy;
  correlation with the O-H..O angle 0.00, with the approach side to the acceptor's bisector 0.02, with H..O
  distance -0.20. At O-O 3.2-3.9 A the mean error is +0.085 +- 0.144 (90-115 deg, n 7), +0.018 +- 0.120
  (115-135, n 9), -0.052 +- 0.056 (135-180, n 4).
- Slopes along the line of centres, 32 pairs at 3.0-3.9 A: error +0.171 +- 0.131 (90-115 deg, n 7),
  +0.105 +- 0.154 (115-135, n 10), +0.119 +- 0.169 (135-180, n 15) kcal/mol per bohr; correlation with angle -0.20.

The error is a nearly constant radial offset, not an angular pattern, so a direction-dependent term fitted here would
mostly act as a radial adjustment, which the Born-Mayer and refit searches already exhausted. Not built.
Candidates left: induction strength, many-body effects beyond three bodies, nuclear quantum effects.


## Sources

Pair-only ab initio water collapses (coordination 11.9 for the CC-pol pair
potential alone) and three-body terms restore tetrahedral coordination:
Bukowski, Szalewicz, Groenenboom, van der Avoird, Science 315, 1249 (2007),
https://pubmed.ncbi.nlm.nih.gov/17332406/ . The neglect of short-range
three-body exchange-repulsion and charge-transfer terms leads to an incorrect
liquid structure: Babin, Medders, Paesani (HBB2-pol), arXiv:1210.7022 .
MB-pol: three-body about 20% of the interaction energy; classical against
quantum density 1.004 against 1.001 g/cm3 at 298 K: Reddy et al., arXiv:1609.02884 .
Beyond-two-body errors of PBE come from close approach of monomers that are not
hydrogen-bonded to each other, BLYP's from hydrogen-bond cooperativity: Gillan,
Alfe, Bartok, Csanyi, arXiv:1303.0751 . Nuclear quantum effects weaken the
hydrogen bond and broaden g(r): Wang,
Ceriotti, Markland, arXiv:1405.7366 (the pressure size quoted in the notes is my inference from the density shift). MP2 with a triple-zeta basis over-structures
liquid water and depends on basis convergence: arXiv:2603.20451; Del Ben,
VandeVondele, Hutter (MP2 and RPA isobaric-isothermal Monte Carlo; RPA 0.994 g/cm3).
Ewald summation of long-range dispersion, and cut-off artefacts in pressure and density: arXiv:0704.0474.
