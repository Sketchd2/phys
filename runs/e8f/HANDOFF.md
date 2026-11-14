# E8f handoff: boiling points from the engine's own laws

Written so that this work can be continued without the session that did it.
Everything here is on branch `claude/laughing-keller-sagqde`. Nothing in it has
been written into `docs/PLAY.md` or `docs/PHYSICS.md`: the owner approves final
results first. Read `CLAUDE.md` before changing anything; its rule on asking the
owner when the plan is silent applies to every "owner decision" marked below.

## 1. What is established

| Result | Number | Where |
|---|---|---|
| Free-energy route validated (methane, two independent checks) | -0.685 +- 0.012 against -0.661 +- 0.038 kcal/mol | `ti-methane/`, `check-slab/` |
| Stored water law boils at | 333.3 K (real 373.1) | `ti-water/path-result.txt` |
| One-constant water law boils at | 327.5 K | `newlaw/water/path-result.txt` |
| One-constant methane law (unweighted fit) boils at | 125.9 K (real 111.7) | `newlaw/methane/path-result.txt` |
| Methane law refitted with liquid-weighted pairs and fitted damping boils at | 98.2 K | `newlaw/methane2/path-result.txt` |
| System size, water, 64 against 216 molecules | within about 4 K | `ti-water/n216/` |
| A law scaled by x at T is the law at T / x | 5% more binding: 348.7 K; 12.7%: 373 K | commit `71e9722` |
| Water dimer in the engine's basis | counterpoise -4.778, uncorrected -5.161, accurate -5.020 | `shortfall/water.txt` |
| Methane dimer | counterpoise -0.475, uncorrected -0.660, accurate -0.530 | `shortfall/methane.txt` |
| Repulsion constant against orbital energy | K = 5.435 / abs(HOMO energy), three species | commit message of "basis bracket" |

What these say together: the route from a law to a boiling point is sound and
cheap; the laws reproduce the pair energies they were fitted to; the pair
energies themselves are short of binding because the basis is incomplete (about
5% for water, far more in relative terms for methane); and about 3 to 5% of
water's 12.7% is still unexplained.

## 2. The tools, and exactly how to run them

All commands are run from the repository root unless a directory is given.
Build once: `cargo build --release`.

**A boiling point from a law (the free-energy path).** A state is one liquid
simulation at one temperature and density.

```sh
# water: 64 molecules, Ewald, started from stored snapshots
cd runs/e8f/newlaw && ./run-states.sh water water/part1.txt
# any other molecule: lattice start, plain cut-off
cd runs/e8f/newlaw && NAME=methane ./run-states3.sh methane2 methane2/part1.txt 108
# then, from runs/e8f, the analysis ("T:B2" pairs come from phys-b2, cm3/mol)
python3 path.py newlaw/water 18.015 900 1.0 "300:-722.6,320:-543.7,..."
```

A list file has one state per line: `T density picoseconds [snapshot]`. A
state is finished when its directory holds a file `done`; the scripts skip
finished states, so after any stop the same command is simply run again.

**The vapour side:** `phys-b2 NAME law.txt T1 T2 ... --orientations 20000`
(needs `grow-NAME.state` in the working directory; copy it from the root).

**A law from the monomer record and one constant:** `newlaw/build2.py` writes
the law; `newlaw/fitk4.py` fits the constant with liquid-weighted pairs;
`newlaw/scan.sh` scans the dispersion damping; `newlaw/core.py` fits a
short-range core from close-contact pairs. Each script's first lines say how
to call it.

**The basis bracket of a dimer:** `phys-shortfall water 0.9 1.0 1.1`
(also `ammonia`, `methane`; `--more Z:l:exponent,...` adds functions).

**Pair energies split into parts:** `phys-overlap NAME PAIRS --out FILE`
(columns: index, separation, total, Hartree-Fock, exact electrostatics, three
density overlaps, Heitler-London energy).

## 3. Checks to make every time (each of these went wrong once)

1. **A law must have no hole.** Run `phys-b2` at one temperature before any
   liquid run. If it prints a number with hundreds of digits, the law has an
   attractive hole at contact (induced dipoles against point charges) and needs
   a core. Water's first one-constant law had this.
2. **A liquid under tension may have cavitated.** If the pressure does not rise
   with density at one temperature (methane: -400 to -670 bar from 0.36 to 0.44
   g/cm3), the box has opened a cavity and those states are not a liquid. Use
   higher densities; do not integrate through them.
3. **An average pair error hides the first shell.** Always print the error by
   separation (`fitk4.py` does). Methane's 0.09 kcal/mol overall hid a 30% bias
   at nearest-neighbour distances.
4. **`phys-lawenergy` reads `pairs-water<tag>.txt` whatever the molecule.** Copy
   or name the pair file that way in the working directory.
5. **A law's sites follow its own monomer shape.** Evaluating a law at a dimer
   whose monomers have other bond lengths gave -3.76 where the law's own
   minimum is -4.6.
6. **Cloud machines are reclaimed when a session is idle.** A long run must be
   watched in the foreground (short waits, one after another) or it dies; the
   run scripts are restartable for this reason.
7. **`PHYS_START_SNAP` is for three-atom molecules only.**
8. **Do not quote a boiling point** from a path whose states used different
   law files. Compare each state's `law.txt` with the intended one (`cmp`).

## 4. Designs to carry on with

Ordered by the evidence behind them. Each has a done-when that can fail.

### D1. Complete the basis until the bracket closes (attacks the shortfall)

*Idea.* The counterpoise and uncorrected energies bracket the complete-basis
answer. Add functions to the basis, by the engine's own growth rule, until the
two agree within a tolerance; then every pair energy is the reference the laws
deserve. One rule, any species, nothing empirical.

*Steps.*
1. `phys-shortfall water 1.0 --more ...` with candidate functions; record the
   bracket width after each addition (section 5 has the first measurements).
2. Find which kind of function closes it fastest (higher angular momentum,
   diffuse, or both) and make that the growth criterion.
3. Recompute a small set of pairs (20 is enough for the constant, see D2) in
   the closed basis; refit; rerun the path.

*Done when:* the water dimer's bracket is under 2% of its binding, its middle
is within 2% of -5.020, and the same recipe closes methane's to within 5% of
-0.530. *Fails if* the bracket stops closing or closes away from the accurate
value (then the method, not the basis, is short).

*Cost:* a water dimer is about 2.5 minutes on 4 cores in the present basis and
grows steeply with basis size. Methanol and larger need the owner's machine.

*Owner decision:* the tolerance, and whether pairs are recomputed for all
species or only where the bracket is wide.

### D2. The repulsion constant from the monomer's own orbital energy

*Idea.* Exchange repulsion is a constant times the overlap of valence
densities, and the constant differs by species (10.70 water, 12.49 ammonia,
10.40 methane). It follows the highest occupied orbital energy: K = 5.435 /
abs(e) gives 10.66, 12.62, 10.00. If that holds, no pair energy is needed to
set the repulsion of a new species.

*Steps.* For each new species: `phys-density NAME` (orbital energies are the
`orbital` lines), `phys-overlap` on about 30 close pairs, fit K, compare with
5.435 / abs(e).

*Done when:* five or more species beyond the three agree to 5%. *Fails if* one
is out by more than 10%; then look at whether the valence definition (which
orbitals are left out) is what breaks.

*Caveat:* three species and one fitted number is weak evidence. This constant
is for the Hartree-Fock exchange alone; the law's net constant also absorbs
short-range attraction and is smaller (5.8 for water).

### D3. One path for any species: fewer states, no cavitation

*Idea.* Because a law scaled by x at T is the law at T / x, the isochore up to
high temperature is already a coupling-strength integration. Choose the
isochore at a density the liquid holds at every temperature on it (above its
zero-pressure density), so no state is under tension.

*Steps.* (1) Three short runs at the lowest temperature to find the density
where the pressure is zero. (2) Isochore at 5% above that density from the
lowest temperature to twice the expected critical temperature, 10 states.
(3) Isotherm at the top from that density to the dilute gas, 10 states.
(4) One density leg per temperature of interest, 3 states each.

*Done when:* it reproduces 327.5 K for the one-constant water law within 3 K
using at most 30 states. *Fails if* the isotherm at high density does not
converge (pressure errors above 5%).

### D4. A general density model for the law's sites

*Idea.* The one-constant law needs the molecule's density as a few Gaussians.
Water's came from a fit to exact electrostatic pair energies; methane's from
the monomer alone; ammonia's from a least-squares fit of the density held to
its own charge, dipole and second moments (`gaussfit.py --moments`,
`newlaw/build3.py`). One procedure for all is missing.

*What was measured:* ammonia with 16 Gaussians fits bound pairs to 0.24
kcal/mol unbiased but is 40% soft at contact and its liquid collapses (-7.0
kcal/mol a molecule, -3200 bar). With 32 the electrostatics improve (0.29 rms)
but the simulation is about 25 times slower than water's.

*Next step:* fit the Gaussians to the density *and* to the exact overlap of a
handful of pairs at once, so the tail that sets the repulsion is right;
`gaussfit.py --check` already prints model against exact overlap.

*Done when:* ammonia's liquid at 240 K and 0.68 g/cm3 holds between -5.5 and
-4.7 kcal/mol a molecule with a pressure within 500 bar of zero.

### D5. A placeholder with no liquid simulation

The second virial coefficient from the law (seconds to compute) gives an energy
and size scale; corresponding states turns those into a boiling point. Expect
10% for simple liquids and worse for water. Useful only so that a substance has
*some* derived boiling point while the real one is computed. *Owner decision:*
whether a placeholder of that quality is wanted at all.

## 5. Measurements made for D1

Water dimer at its minimum, kcal/mol (accurate value -5.020):

| Basis | Functions | Counterpoise | Uncorrected | Gap |
|---|---|---|---|---|
| the engine's | 428 | -4.778 | -5.161 | 0.383 |
| + diffuse s, p, d on O and s, p on H | 464 | -4.781 | -5.119 | 0.338 |
| + one g on O (2.0) and one f on H (1.6) | 498 | -4.826 | -5.045 | 0.219 |
| + both of the above | 534 | -4.828 | -5.023 | 0.195 |
| + a second g and f, an f on O and a d on H | 648 | -4.870 | -5.021 | 0.151 |

Higher angular momentum closes the bracket (43% narrower for one shell, 9
minutes instead of 2.5); diffuse functions barely move it. The two ends are
closing on about -4.95, which is where MP2 at the complete basis is expected,
a little above the accurate -5.020. With two further levels the uncorrected
value settles at -5.02 (-5.045, -5.023, -5.021) and the counterpoise one keeps
climbing toward it: for water the complete basis gives the accurate answer,
and the uncorrected energy gets there first.

Ammonia dimer at its minimum: counterpoise -3.049, uncorrected -3.236,
accurate -3.171 (65% across the gap; water 63%, methane 30%).

Along the water dimer's approach the gap is 5.5 to 8% of the binding from 2.9
to 4.7 A; along methane's it is 39% at the minimum and 58% at 4.6 A. So the
order of work for D1 is methane first: it is where the reference is worst.

What D1 should deliver, by the scaling identity (an estimate, not a run): the
re-weighted methane law boils at 98.2 K; scaled by 0.530 / 0.475, the dimer's
accurate over counterpoise binding, it boils at 108.3 K against methane's
111.7. The same exercise for water (5% at the dimer) gives 348.7 K against
373.1, so water has a further shortfall that the dimer minimum does not show.

## 6. Open, for the owner

- Three tests fail on the base commit and on this branch alike
  (`boiling_points_are_right_to_about_a_quarter`,
  `a_branch_lands_on_the_next_tree`, `rain_on_a_basin_stands_as_a_lake`); the
  owner said to ignore them for this work. `s22_subset` did not finish in 70
  minutes on a cloud machine.
- The fitted boiling law in `src/chem/analyse.rs` is still in place. Replacing
  it with any of the above is the owner's decision.
- Dispersion damping is a second fitted constant per species (methane 1.3,
  ammonia 1.4, water flat between 1.5 and 2.0). Whether it can be derived
  (it should follow the repulsion's own decay) has not been tried.
- The second virial coefficient of real water quoted in the commits (about
  -450 cm3/mol at 373 K) is from memory; check it against Harvey and Lemmon
  (2004) before relying on it.
