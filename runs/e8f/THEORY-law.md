# Why the fleet's laws fail, and how to get good ones without breaking an axiom

Written by the fleet orchestrator on 2026-11-27 from the fleet's own data (17 laws, 4 boiling
points) and four probes run on those laws. Nothing in the engine, `law.py`, `PLAY.md` or
`PHYSICS.md` was changed: what form a law takes is the owner's decision. Every number here is
reproducible from `runs/e8f/species/*/` and `target/release/phys-b2`.

## What the data say

**1. The one fitted repulsion constant decides everything.** `law.py` fits K (and a damping
exponent) to the pairs' total energies, weighting a repulsive pair by exp(-E/e), so the contact
geometries where a hole lives carry almost no weight. Every law with a hole has K <= 5.89; every
law without one has K >= 6.59. Among the four with a liquid, the boiling-point error falls in
order of K over the exchange the electrons imply:

| species | K fit / K orbital | boils at | real |
|---|---|---|---|
| ammonia | 0.55 | 577 K | 240 K |
| ethane | 0.59 | 316 K | 185 K |
| acetonitrile | 0.66 | 415 K | 355 K |
| ethanol | 0.76 | 148 K | 351 K |

**2. The exchange the fit should have found is already measured.** `full-NAME.txt` carries, for
every pair, the exact electrostatics, the density overlap S and the Heitler-London energy, so
exchange = E_hl - E_es is measured directly. Exchange / S is 6.8 to 10.6 (median per species) for
all 17 species, and exchange / S_valence is about 10.4 (median of medians, 8.2 to 13.4): close to
one constant, and close to HANDOFF D2's 5.435 / |HOMO| (10.0 to 14.4). The fits took 23 to 76% of it.

**3. Three short-range faults, separated by probes** (B2 at 300 K, 1500 orientations, each law
otherwise unchanged):

| species | as fitted | measured K | fitted K, no induction | measured K, no induction |
|---|---|---|---|---|
| methane | hole | hole | -278 | +87 |
| ethylene | hole | hole | hole | -87 |
| dichloromethane | hole | hole | hole | +103 |
| methanol | hole | hole | -1.1e5 | -1,008 |
| dimethyl ether | hole | hole | hole | -1,572 |
| water | hole | -1.75e4 | -1.1e4 | -7,780 |
| methylamine | hole | -3,521 | -5,175 | -2,350 |
| acetaldehyde | hole | -1.4e4 | hole | -9,847 |
| acetic acid | hole | -4.2e5 | hole | -9.3e4 |
| formic acid | hole | hole | -1.2e15 | -4.7e5 |
| hydrazine | hole | hole | hole | -9.9e5 |
| hydrogen peroxide | hole | hole | -7.9e10 | -1.6e8 |
| chloroform | hole | hole | hole | hole |

(cm3/mol; "hole" = infinite or hundreds of digits. In methane, switching dispersion off instead of
induction changes nothing.)

- **Repulsion too weak** (fault A). The measured K alone closes 4 of 13 holes.
- **Point polarisabilities with no short-range damping** (fault B). Induction is the whole hole
  in methane, and with the measured K and no induction 8 of the 13 become ordinary.
- **Bare point nuclei against compact clouds** (fault C). The hydrogen-bonding species (formic acid,
  hydrazine, hydrogen peroxide, acetic acid) stay 100 to 100,000 times over-bound with neither
  fault A nor B, and chloroform keeps a true hole: a nucleus with no cloud of its own reaches into
  the partner's tight negative cloud with nothing to stop it. The hand-built water law that boiled
  at 325.6 K had an O-H core term (`pair 0 1 1238.76 ...`); `law.py` writes none.

The same faults explain the boiling points, not only the holes: ammonia's law has B2(300 K) about
-1e5 cm3/mol against roughly -265 for real ammonia (from memory), i.e. wells far deeper than any
of its 12 reference pairs, in orientations the pairs never sampled. A hole and an over-bound liquid
are one defect at two strengths.

## The theory

The long-range parts of every law (charges, C6, polarisabilities) are derived from the molecule's
own electrons and are fine. **Every failure is at short range, where three things the electrons
determine are instead fitted to a total or left out**: the strength of exchange repulsion, the
damping of induction, and the shape of the charge where molecules interpenetrate. A fit to total
energies over a dozen near-equilibrium pairs lets one constant absorb all three; errors cancel
where the pairs are and add up where the liquid goes.

So, in the spirit of *measure, never be told* and *derived, with shortcuts stored*:
**derive each term of the law from its own measured component, fit nothing to the total, and use
the total only as the test.** The pair files already hold each component.

| law term | measured against | replaces |
|---|---|---|
| electrostatics (clouds) | E_es, the exact electrostatics of the two densities, per pair | clouds fitted to the monomer only |
| exchange repulsion | E_hl - E_es against the model overlap (one K, from the pairs or 5.435/abs(HOMO)) | K fitted to totals |
| induction | E_hf - E_hl, with the induced dipoles responding to the partner's clouds, not to point charges (damping width from the clouds' own sigma) | undamped point polarisabilities |
| dispersion | E_total - E_hf (MP2 correlation), damped with the exchange's own decay | damping exponent fitted to totals |

Zero numbers fitted to totals. Then, to stop validating only where the pairs already are: **ask
the law where it is least sure.** The B2 integrator already finds the law's deepest
configurations; compute those as the next reference pairs (cheap for 2-3 heavy atoms), and add
them until B2 stops moving. That is measurement, not instruction, and it is one rule for every
species.

## Tests, each of which can fail

1. **Component fits hold up.** For water, methane and ammonia, each term against its own component
   is within 10% at every reference pair. *Fails if* one term cannot match its component with the
   form it has (that term's form is then wrong, not its constant).
2. **No holes.** With all four terms derived, B2 is finite for all 17 species. *Fails if* any
   stays holed (then a fourth short-range fault exists; chloroform is the likeliest).
3. **B2 is right.** Derived B2(300 K) within 30% of the measured value for water, ammonia,
   methane and ethane (take the measured values from a reference such as Harvey and Lemmon
   2004, not from memory). *Fails if* the laws are now all under-bound, meaning the
   fitted K was also standing in for real short-range attraction that no term carries.
4. **Boiling points move the right way.** Ammonia and ethane fall, ethanol rises. *Fails if*
   ethanol does not rise: its law has K at 0.76 and its reference density is 0.435 against
   0.789, so its error may be a different fault (it uses the MP2 mean shared equally among
   heavy atoms, not phys-polar's split).

## What this does not settle

- The model overlap of Gaussian clouds is assumed to track the exact S; HANDOFF D4 measured
  ammonia's 40% soft at contact. Test 1 checks this.
- B2 from 1500 random orientations is a rare-event estimate: ammonia's 600 K value moved the
  wrong way when K rose. Numbers in the 1e4 to 1e5 range are indicative, not precise.
- The basis shortfall (HANDOFF D1: 5% for water, 30-60% for methane) is separate and still applies
  to the dispersion component.
- The fleet's large species are blocked on memory, not on this; nothing here helps them.

## Owner decisions

Changing the law's form (damped induction, a derived core, component fits) is `PHYSICS.md`
weight. So is whether the reference pair set may grow by asking the law where it is unsure.
