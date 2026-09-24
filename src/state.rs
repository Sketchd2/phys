//! Matter state: what the engine stores when it is *not* storing particles.
//!
//! A node holds its contents as `Matter`. Everything below its own
//! resolution is absent — not approximated, absent — and is regenerated on
//! demand by `sampler.rs`. For that to be legitimate, that matter must
//! carry every quantity that the missing detail is *not allowed to change*:
//! the conserved set. If refinement and re-coarsening return exactly the same
//! conserved tuple, no experiment performed at the coarse scale can tell
//! whether the detail was ever there.
//!
//! That is the engine's central correctness claim, and `Conserved` is the
//! object it is stated in terms of.

use crate::math::{det_sum, Mat3, Quat, Vec3};
use crate::units::*;

/// The invariant set. Every scale transition preserves this exactly (to
/// round-off), at every tier, in both directions.
///
/// Baryon and lepton number are here because they are what make the *subatomic*
/// tier consistent with the galactic one: you cannot fuse hydrogen in a star,
/// coarsen the star, refine it again and find the protons back. The bookkeeping
/// spans 60 orders of magnitude in scale precisely because these are additive.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Conserved {
    /// Total energy including rest mass, J.
    pub energy: f64,
    pub momentum: Vec3,
    /// About the node's centre of mass.
    pub angular_momentum: Vec3,
    /// Coulombs.
    pub charge: f64,
    /// Net baryon number (protons + neutrons, minus antibaryons).
    pub baryon: f64,
    /// Net lepton number.
    pub lepton: f64,
}

impl Conserved {
    pub fn zero() -> Self {
        Self::default()
    }

    pub fn add(self, o: Self) -> Self {
        Conserved {
            energy: self.energy + o.energy,
            momentum: self.momentum + o.momentum,
            angular_momentum: self.angular_momentum + o.angular_momentum,
            charge: self.charge + o.charge,
            baryon: self.baryon + o.baryon,
            lepton: self.lepton + o.lepton,
        }
    }

    pub fn sub(self, o: Self) -> Self {
        Conserved {
            energy: self.energy - o.energy,
            momentum: self.momentum - o.momentum,
            angular_momentum: self.angular_momentum - o.angular_momentum,
            charge: self.charge - o.charge,
            baryon: self.baryon - o.baryon,
            lepton: self.lepton - o.lepton,
        }
    }

    /// Worst discrepancy against a reference, measured against the *natural
    /// scale* of each quantity rather than against its own net value.
    ///
    /// This distinction is not pedantry, it is forced by the physics. Consider
    /// a hot molecular cloud with almost no net rotation: each of its particles
    /// carries angular momentum of order `m r v`, but the sum very nearly
    /// cancels, and the net can easily be 20 orders of magnitude smaller than
    /// any individual term. Double precision holds 16 digits. The net angular
    /// momentum of such a cloud is therefore *not representable* — not because
    /// the algorithm is careless, but because the information is below the
    /// noise floor of the arithmetic, and would be below the noise floor of any
    /// finite arithmetic.
    ///
    /// Dividing the error by the net value in that situation reports 100% and
    /// tells you nothing. Dividing by the total angular momentum *content* —
    /// `sum |r_i x p_i|` — tells you what is actually true and what actually
    /// matters: that the engine's bookkeeping is good to one part in 10^13 of
    /// everything in the system. An observer cannot detect an error smaller
    /// than that, because no measurement they can make inside the simulation
    /// has access to a finer distinction.
    ///
    /// `Scales::of` computes those denominators from a materialised set.
    pub fn error_against(&self, reference: &Conserved, scales: &Scales) -> f64 {
        fn rel(a: f64, b: f64, scale: f64) -> f64 {
            let d = (a - b).abs();
            if scale > 0.0 {
                d / scale
            } else if d > 0.0 {
                1.0
            } else {
                0.0
            }
        }
        let e = rel(self.energy, reference.energy, scales.energy);
        let p = {
            let d = (self.momentum - reference.momentum).norm();
            if scales.momentum > 0.0 { d / scales.momentum } else { 0.0 }
        };
        let l = {
            let d = (self.angular_momentum - reference.angular_momentum).norm();
            if scales.angular_momentum > 0.0 { d / scales.angular_momentum } else { 0.0 }
        };
        let q = rel(self.charge, reference.charge, scales.charge);
        let b = rel(self.baryon, reference.baryon, scales.baryon);
        let le = rel(self.lepton, reference.lepton, scales.lepton);
        e.max(p).max(l).max(q).max(b).max(le)
    }

    /// Worst *relative* discrepancy against a reference, per component.
    ///
    /// Relative rather than absolute because the same tuple has to be checked
    /// at 10^41 J (a galaxy) and 10^-13 J (a nucleus); an absolute tolerance is
    /// meaningless across that range.
    pub fn max_relative_error(&self, reference: &Conserved) -> f64 {
        fn rel(a: f64, b: f64, scale: f64) -> f64 {
            let d = (a - b).abs();
            let s = scale.max(b.abs()).max(a.abs());
            if s > 0.0 {
                d / s
            } else {
                0.0
            }
        }
        let e = rel(self.energy, reference.energy, 0.0);
        let p = {
            let d = (self.momentum - reference.momentum).norm();
            let s = reference.momentum.norm().max(self.momentum.norm());
            // Momentum can legitimately be ~0 in a node's own frame; compare
            // against the energy scale (E/c is the natural momentum unit).
            let floor = reference.energy.abs() / C;
            if s.max(floor) > 0.0 {
                d / s.max(floor)
            } else {
                0.0
            }
        };
        let l = {
            let d = (self.angular_momentum - reference.angular_momentum).norm();
            let s = reference
                .angular_momentum
                .norm()
                .max(self.angular_momentum.norm());
            if s > 0.0 {
                d / s
            } else {
                0.0
            }
        };
        let q = rel(self.charge, reference.charge, E_CHARGE);
        let b = rel(self.baryon, reference.baryon, 1.0);
        let le = rel(self.lepton, reference.lepton, 1.0);
        e.max(p).max(l).max(q).max(b).max(le)
    }

    pub fn is_finite(&self) -> bool {
        self.energy.is_finite()
            && self.momentum.is_finite()
            && self.angular_momentum.is_finite()
            && self.charge.is_finite()
            && self.baryon.is_finite()
            && self.lepton.is_finite()
    }
}

/// Natural magnitudes of each conserved quantity in a system: the denominators
/// that make a conservation error meaningful. See `Conserved::error_against`.
#[derive(Debug, Clone, Copy)]
pub struct Scales {
    pub energy: f64,
    pub momentum: f64,
    pub angular_momentum: f64,
    pub charge: f64,
    pub baryon: f64,
    pub lepton: f64,
}

impl Scales {
    /// Sum of the *magnitudes* of every contribution — the total amount of each
    /// quantity present, as opposed to the net that survives cancellation.
    pub fn of(bodies: &[Body]) -> Scales {
        let n = bodies.len();
        if n == 0 {
            return Scales::unit();
        }
        let mass = det_sum_by(n, &|i| bodies[i].mass);
        let com = if mass > 0.0 {
            det_sum_v3_by(n, &|i| bodies[i].pos.scale(bodies[i].mass)).scale(1.0 / mass)
        } else {
            Vec3::ZERO
        };
        let momentum = det_sum_by(n, &|i| bodies[i].momentum().norm());
        let angular = det_sum_by(n, &|i| {
            let b = &bodies[i];
            (b.pos - com).cross(b.momentum()).norm() + b.spin.norm()
        });
        let charge = det_sum_by(n, &|i| bodies[i].charge.abs());
        let baryon = det_sum_by(n, &|i| bodies[i].mass * bodies[i].composition.nucleons_per_kg());
        Scales {
            energy: (mass * C2 + det_sum_by(n, &|i| bodies[i].internal_energy.abs())).max(1e-300),
            momentum: momentum.max(mass * 1e-30),
            angular_momentum: angular.max(1e-300),
            charge: charge.max(E_CHARGE),
            baryon: baryon.max(1.0),
            lepton: baryon.max(1.0),
        }
    }

    pub fn unit() -> Scales {
        Scales {
            energy: 1.0,
            momentum: 1.0,
            angular_momentum: 1.0,
            charge: 1.0,
            baryon: 1.0,
            lepton: 1.0,
        }
    }
}

/// Mass fractions by coarse element. Always sums to 1 for a non-empty node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Composition(pub [f64; COARSE_ELEMENTS]);

impl Default for Composition {
    fn default() -> Self {
        Composition::primordial()
    }
}

impl Composition {
    /// Big Bang nucleosynthesis output: the initial condition for gas that has
    /// never been through a star.
    pub fn primordial() -> Composition {
        let mut c = [0.0; COARSE_ELEMENTS];
        c[CoarseElement::Hydrogen as usize] = 0.75;
        c[CoarseElement::Helium as usize] = 0.25;
        Composition(c)
    }

    /// Roughly solar (Asplund 2009 mass fractions, lumped into our buckets).
    pub fn solar() -> Composition {
        let mut c = [0.0; COARSE_ELEMENTS];
        c[CoarseElement::Hydrogen as usize] = 0.7381;
        c[CoarseElement::Helium as usize] = 0.2485;
        c[CoarseElement::Carbon as usize] = 0.0024;
        c[CoarseElement::Nitrogen as usize] = 0.0007;
        c[CoarseElement::Oxygen as usize] = 0.0057;
        c[CoarseElement::Silicon as usize] = 0.0007;
        c[CoarseElement::Iron as usize] = 0.0013;
        c[CoarseElement::Other as usize] = 0.0026;
        Composition(c).normalised()
    }

    /// Stated, element by element, and normalised. What a scenario uses when
    /// it is saying what something is made of — the composition equivalent of
    /// stating a mixture, and the honest replacement for a table of
    /// per-species substrates.
    pub fn of(parts: &[(CoarseElement, f64)]) -> Composition {
        let mut c = [0.0; COARSE_ELEMENTS];
        for (e, f) in parts {
            c[*e as usize] += f.max(0.0);
        }
        Composition(c).normalised()
    }

    /// Carbon, hydrogen and oxygen in the proportions of a carbohydrate —
    /// `CH2O`, 44/6/50 by mass.
    ///
    /// **A scenario's statement, not a program's property.** It sits beside
    /// [`Composition::primordial`] and [`Composition::solar`] for the same
    /// reason they do: a scenario has to say what its matter is made of, and
    /// saying so is not a table the engine reads. `docs/PLAY.md` D11 retired
    /// `Program::substrate`, which was the version of this that the *engine*
    /// consulted; what a thing grows into is now the feedstock it grew from.
    pub fn organic() -> Composition {
        Composition::of(&[
            (CoarseElement::Carbon, 0.44),
            (CoarseElement::Hydrogen, 0.06),
            (CoarseElement::Oxygen, 0.50),
        ])
    }

    /// Crustal silicate rock: oxygen and silicon with iron and the rest lumped.
    ///
    /// The same kind of statement as [`Composition::organic`].
    pub fn crustal() -> Composition {
        Composition::of(&[
            (CoarseElement::Oxygen, 0.46),
            (CoarseElement::Silicon, 0.28),
            (CoarseElement::Iron, 0.09),
            (CoarseElement::Other, 0.17),
        ])
    }

    /// Nothing stated: every bucket zero.
    ///
    /// Distinct from any real composition, and it means "this has not been
    /// built out of anything yet" rather than "this is made of nothing".
    pub fn none() -> Composition {
        Composition([0.0; COARSE_ELEMENTS])
    }

    pub fn is_none(&self) -> bool {
        self.0.iter().all(|v| *v <= 0.0)
    }

    /// Pure one element — used when the user drills into a specific atom.
    pub fn pure(s: CoarseElement) -> Composition {
        let mut c = [0.0; COARSE_ELEMENTS];
        c[s as usize] = 1.0;
        Composition(c)
    }

    /// Cold dark matter: gravitationally active, chemically inert. Lives in the
    /// `Other` bucket but is flagged separately by the tier's solver.
    pub fn dark() -> Composition {
        Composition::pure(CoarseElement::Other)
    }

    pub fn normalised(mut self) -> Composition {
        let s = det_sum(&self.0);
        if s > 0.0 {
            for v in self.0.iter_mut() {
                *v /= s;
            }
        }
        self
    }

    pub fn get(&self, s: CoarseElement) -> f64 {
        self.0[s as usize]
    }

    /// Metallicity: everything heavier than helium.
    pub fn metallicity(&self) -> f64 {
        det_sum(&self.0[2..])
    }

    /// Mass-weighted blend of two compositions.
    pub fn blend(a: Composition, ma: f64, b: Composition, mb: f64) -> Composition {
        let t = ma + mb;
        if t <= 0.0 {
            return a;
        }
        let mut c = [0.0; COARSE_ELEMENTS];
        for i in 0..COARSE_ELEMENTS {
            c[i] = (a.0[i] * ma + b.0[i] * mb) / t;
        }
        Composition(c)
    }

    /// Mean molecular mass in kg, assuming full ionisation above 10^4 K and
    /// neutral below. This one number sets the pressure, the sound speed, the
    /// Jeans mass and the thermal velocity, so it is worth getting right.
    pub fn mean_molecular_mass(&self, temperature: f64) -> f64 {
        let ionised = temperature > 1.0e4;
        let mut inv = 0.0;
        for s in CoarseElement::ALL {
            let x = self.get(s);
            if x <= 0.0 {
                continue;
            }
            let particles_per_nucleus = if ionised { 1.0 + s.z() } else { 1.0 };
            inv += x * particles_per_nucleus / s.a();
        }
        // mu = AMU / sum(x_i * particles_per_nucleus_i / A_i).
        //
        // Written directly in terms of AMU rather than routed through
        // Avogadro's number and a gram-to-kilogram factor, because that route
        // is where a stray 10^3 hides — and a factor of 1000 here silently
        // scales every temperature, pressure, sound speed and Jeans length in
        // the engine.
        if inv > 0.0 {
            AMU / inv
        } else {
            M_PROTON
        }
        .max(M_ELECTRON)
    }

    /// Electrons per nucleon — needed for opacity and for charge bookkeeping.
    pub fn electrons_per_nucleon(&self) -> f64 {
        let mut n = 0.0;
        for s in CoarseElement::ALL {
            n += self.get(s) * s.z() / s.a();
        }
        n
    }

    /// Nucleons per kilogram.
    /// Mean mass of one atom of this mixture, kg.
    ///
    /// Not the same as the mean molecular mass, which is per free particle and
    /// depends on how ionised the gas is. This is per *atom*, which is what
    /// bounds how many atoms a given mass can be split into.
    pub fn mean_atomic_mass(&self) -> f64 {
        let mut mass = 0.0;
        let mut number = 0.0;
        for s in CoarseElement::ALL {
            let f = self.get(s);
            if f <= 0.0 {
                continue;
            }
            let a = s.a();
            mass += f;
            number += f / a;
        }
        if number > 0.0 {
            mass / number * AMU
        } else {
            AMU
        }
    }

    pub fn nucleons_per_kg(&self) -> f64 {
        let mut n = 0.0;
        for s in CoarseElement::ALL {
            n += self.get(s) / s.mass_kg() * s.a();
        }
        n
    }

    /// Nuclear binding energy per kg relative to free nucleons (negative).
    ///
    /// This is *not* part of the energy budget: it is already inside the rest
    /// mass. It appears in the invariant set only through the composition,
    /// which is conserved exactly. Fusion releases energy by *changing* this
    /// number — the burning solver takes the difference and adds it as heat.
    /// Treating it as an available energy pool (an easy mistake, since it has
    /// units of energy) injects ~10^13 J/kg of spurious heat.
    pub fn nuclear_energy_per_kg(&self) -> f64 {
        let mut e = 0.0;
        for s in CoarseElement::ALL {
            let x = self.get(s);
            if x <= 0.0 {
                continue;
            }
            // per kg of this element: (A nucleons / mass) * B/A
            e -= x * (s.a() / s.mass_kg()) * s.binding_per_nucleon_mev() * MEV;
        }
        e
    }
}

/// What a node is made of: everything it knows about itself without
/// materialising its children.
///
/// 248 bytes. A galaxy's worth (a few million live nodes) is well under a
/// gigabyte, which is what makes the whole approach fit in a 6 GB card
/// alongside the materialised working set.
///
/// # Why this is not just fields on `Node`
///
/// It is `Copy` and `Node` is not — `Node` holds body lists, children, a
/// morphology — and that difference is load-bearing. The engine copies the
/// matter out of a node to release the borrow on the tree before doing work
/// that needs `&mut self`, which it could not do with something holding a
/// `Vec`.
///
/// It is also the unit the scale transform is defined over: `summarise` takes
/// bodies and returns this, `sample` takes this and returns bodies, and neither
/// has any business seeing a node's address or its scheduling timestamps. The
/// engine's central guarantee, `summarise(sample(m)) = m` on the conserved set,
/// is a statement about a *value*, so it has to be one. Fifty-odd places build
/// one with no node behind it at all — the scenario builders, and `promote`,
/// where a single body is reinterpreted as the matter of a new node.
///
/// # Why it was called `Aggregate`
///
/// Because it can be produced by aggregating children. But that is the rarer
/// direction: usually there are no children, and they are produced from *this*.
/// The name claimed a direction the type does not have, and `Bulk` — the other
/// candidate — would have claimed a resolution it does not have either, since a
/// materialised node carries one of these kept in step with its bodies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matter {
    pub mass: f64,
    /// Centre of mass offset from the node's own origin (metres, node frame).
    /// Kept near zero by construction; drift here is a diagnostic.
    pub com: Vec3,
    /// Bulk momentum in the node's own frame.
    pub momentum: Vec3,
    /// Intrinsic (spin) angular momentum about the centre of mass.
    pub spin: Vec3,
    /// Random/thermal kinetic energy plus chemical and nuclear binding, J.
    /// Excludes rest mass and excludes bulk motion.
    pub internal_energy: f64,
    /// Binding held by the node's *own long-range field* — gravity, and at fine
    /// tiers the electrostatic self-energy of a charged ball. Negative for a
    /// bound object. Tracked explicitly so that coarsening does not silently
    /// destroy it.
    ///
    /// **This is the half of the binding that spreading the contents out
    /// releases**, and that is the whole reason it is a field of its own rather
    /// than one number covering every kind of binding. `sampler::sample` has a
    /// relaxation loop whose premise is exactly that: a configuration too
    /// tightly bound to hold the energy it claims must be bigger, so scale it
    /// up until the budget turns positive. Scaling moves `phi`, and `phi` is
    /// this term.
    ///
    /// See [`Matter::cohesive_binding`] for the half it does not move, and why
    /// running the loop against the sum inflated every solid by 4.3x10^5.
    pub gravitational_binding: f64,
    /// Binding held by *bonds* — chemical, electronic, nuclear. Negative for a
    /// bound object, and zero for anything held together only by its own
    /// gravity.
    ///
    /// Split out of what used to be one `binding_energy` because expansion
    /// does not release it. A granite block's binding is a silicate cohesive
    /// energy of about 5 eV per atom; its self-gravity at that size is around
    /// 10^-4 J and utterly negligible against it. The relaxation loop in
    /// `sampler::sample` scaled the block up by 1.5 thirty-two times trying to
    /// make a budget positive by moving a term that was not the negative one,
    /// then gave up and left the inflation in place — measured, a cubic metre
    /// of granite sampled its contents 4.396x10^5 radii outside the node they
    /// are supposed to be inside, against 3.1 for a spiral galaxy.
    ///
    /// The blast radius was three scenarios while only three scenarios set a
    /// real cohesive energy. `PLAY.md` D17 puts a `Mixture` on every `Matter`,
    /// which is what gives every solid in the world one — so this had to be
    /// split first, and `PLAY.md` §7 puts it before D17 for that reason.
    ///
    /// Carried through both directions of a scale transition untouched, like
    /// [`Matter::external_potential`] and [`Matter::chemical_energy`]: the
    /// bonds are at a scale far below the node, so nothing the sampler does to
    /// the arrangement of its bodies can be read back out of them.
    pub cohesive_binding: f64,
    /// Energy this node holds by virtue of sitting in someone *else's*
    /// potential — a dark matter halo, a parent star, an applied field.
    ///
    /// Kept separate from `binding_energy` because it is not recoverable from
    /// the node's own contents: refine a galaxy's baryons and their mutual
    /// potential is nine times smaller than the halo's grip on them. Folding
    /// the two together makes every refinement demand a thermal budget that
    /// does not exist, and the sampler is then forced to either invent energy
    /// or violate the virial theorem. This field is carried through scale
    /// transitions untouched, by both directions.
    pub external_potential: f64,
    /// Characteristic radius, m.
    pub radius: f64,
    pub temperature: f64,
    pub composition: Composition,
    pub charge: f64,
    pub baryon_number: f64,
    pub lepton_number: f64,
    /// Thermodynamic entropy of the node's *contents*, J/K.
    ///
    /// This may legitimately fall: a growing organism or a building under
    /// construction becomes more ordered. What may never fall is the total —
    /// this plus everything dumped into the surroundings, tracked below.
    pub entropy: f64,
    /// Cumulative entropy exported to the environment as waste heat, J/K.
    ///
    /// Without this account the second law cannot be checked at all, because
    /// the interesting processes are precisely the ones that lower local
    /// entropy while raising the total. Growth is not a violation; it is a
    /// transaction, and this is the other side of it.
    pub entropy_exported: f64,
    /// Free energy stored in chemical bonds and ordered structure, J.
    ///
    /// Distinct from `internal_energy`, which is thermal. Biomass holds about
    /// 17 MJ/kg here; a steel frame holds its embodied energy. Destroying the
    /// structure releases it — which is what makes a forest fire an energy
    /// source rather than a rendering effect.
    ///
    /// Like `external_potential`, this is not recoverable from the children
    /// alone, so both directions of a scale transition carry it through
    /// unchanged and the caller reinstates it.
    pub chemical_energy: f64,
    /// Magnetic energy density integrated over the node, J. Drives the ISM.
    pub magnetic_energy: f64,
    /// Bolometric luminosity, W — what the node emits, for observation.
    pub luminosity: f64,
    /// What this matter is made of, by substance. `docs/PLAY.md` D17.
    ///
    /// # Why this is state and not a lookup
    ///
    /// It lived in a side table keyed by `EntityId`, for the handful of nodes
    /// somebody had described. That made "what is this made of" unanswerable
    /// for almost every node in the world, and D13's measured solidity and
    /// D14's cohesive energy both need an answer.
    ///
    /// There is no derivation available. `Registry::intern` takes an
    /// `Arrangement` — actual atoms and bonds — and `chem::analyse` works only
    /// from one; a `Composition` is eight mass fractions of coarse element
    /// buckets, and **there is no path between them**: wood and coal have
    /// nearly the same composition, and no amount of carbon-hydrogen-oxygen
    /// arithmetic yields "oak". Mineralogy from bulk composition is non-unique
    /// in the same way — the same composition is basalt or granite depending on
    /// cooling history.
    ///
    /// So `Composition` is the wrong resolution to carry substance identity,
    /// and that is what makes putting `Mixture` here axiom-consistent rather
    /// than a lookup. *Which substance something is* **is state**. Nothing is
    /// told anything: `summarise` blends mixtures when detail collapses,
    /// `promote` carries one down when a body becomes a node, and the answer
    /// travels from wherever the matter came from by the same transform that
    /// carries mass and energy. `SubstanceId::UNSPECIATED` stays a real answer
    /// rather than a failure — matter at ten million kelvin genuinely has no
    /// molecular identity.
    ///
    /// # What it settles
    ///
    /// Solidity, which D13 named and did not specify. Phase is already a law
    /// the engine has, derived from temperature against a melting point
    /// `chem::analyse` computes, so `mixture.in_phase(Phase::Solid)` is the
    /// solid mass fraction and D13's "measured solidity" is a reading.
    ///
    /// # What it costs
    ///
    /// 200 bytes against `Matter`'s own 248 and a `Node`'s 1,136. `chem::react`
    /// runs at 0.069 µs per node for one substance and 0.934 for eight, which
    /// is 1.4% to 18.7% of a frame at 10^4 nodes — so it is gated on a mixture
    /// that is actually non-trivial rather than run over `UNSPECIATED`.
    pub mixture: crate::chem::Mixture,
}

impl Default for Matter {
    fn default() -> Self {
        Matter {
            mass: 0.0,
            com: Vec3::ZERO,
            momentum: Vec3::ZERO,
            spin: Vec3::ZERO,
            internal_energy: 0.0,
            gravitational_binding: 0.0,
            cohesive_binding: 0.0,
            external_potential: 0.0,
            radius: 1.0,
            temperature: 2.725, // CMB floor: nothing in the engine is colder
            composition: Composition::primordial(),
            charge: 0.0,
            baryon_number: 0.0,
            lepton_number: 0.0,
            entropy: 0.0,
            entropy_exported: 0.0,
            chemical_energy: 0.0,
            magnetic_energy: 0.0,
            luminosity: 0.0,
            mixture: crate::chem::Mixture::new(),
        }
    }
}

impl Matter {
    /// Neutral matter of the given composition: charge zero, and baryon/lepton
    /// numbers implied by the mass. This is the normal way to create matter.
    pub fn neutral(mass: f64, radius: f64, temperature: f64, composition: Composition) -> Matter {
        let nucleons = mass * composition.nucleons_per_kg();
        let electrons = nucleons * composition.electrons_per_nucleon();
        let mut a = Matter {
            mass,
            radius,
            temperature,
            composition,
            baryon_number: nucleons,
            lepton_number: electrons,
            charge: 0.0,
            ..Default::default()
        };
        a.internal_energy = a.thermal_energy();
        a.entropy = a.estimate_entropy();
        a
    }

    /// Set a net charge by removing (or adding) electrons, keeping lepton
    /// number consistent.
    ///
    /// Charge and lepton number are not independent: an iron nucleus stripped
    /// to +26e has no electrons left, and saying otherwise creates 26 leptons
    /// out of nothing. Setting `charge` directly is therefore a trap, so the
    /// supported path adjusts both together.
    pub fn with_charge(mut self, charge: f64) -> Matter {
        let electrons_removed = charge / E_CHARGE;
        self.charge = charge;
        self.lepton_number -= electrons_removed;
        self
    }

    /// Check the internal consistency the invariants assume. Returns the worst
    /// relative violation; zero for a well-formed state.
    ///
    /// Used on the authoring path, where a user can set fields directly, and in
    /// the tests. The engine never produces a state that fails this, so a
    /// non-zero result always points at an external cause.
    pub fn validate(&self) -> f64 {
        let expected_lepton =
            self.baryon_number * self.composition.electrons_per_nucleon() - self.charge / E_CHARGE;
        let scale = self.lepton_number.abs().max(expected_lepton.abs()).max(1.0);
        let lepton_err = (self.lepton_number - expected_lepton).abs() / scale;
        let expected_baryon = self.mass * self.composition.nucleons_per_kg();
        let baryon_err = (self.baryon_number - expected_baryon).abs()
            / self.baryon_number.abs().max(expected_baryon.abs()).max(1.0);
        lepton_err.max(baryon_err)
    }

    /// Number of constituent particles at the current temperature.
    pub fn particle_count(&self) -> f64 {
        let mu = self.composition.mean_molecular_mass(self.temperature);
        if mu > 0.0 {
            self.mass / mu
        } else {
            0.0
        }
    }

    /// (3/2) N k T for an ideal gas.
    pub fn thermal_energy(&self) -> f64 {
        1.5 * self.particle_count() * K_B * self.temperature
    }

    /// Joules per kelvin. The derivative of [`Self::thermal_energy`] with
    /// respect to temperature, and derived rather than tabulated for the
    /// obvious reason: a specific heat looked up per material is a table of
    /// answers, and equipartition is the law underneath it.
    ///
    /// It is the quantity [`crate::neighbourhood::exchange`] needs to know how
    /// far two things can move each other before they meet in the middle. Note
    /// that it depends on temperature through the ionisation in
    /// `mean_molecular_mass`, which is correct and is why it is a method.
    pub fn heat_capacity(&self) -> f64 {
        1.5 * self.particle_count() * K_B
    }

    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * std::f64::consts::PI * self.radius.powi(3)
    }

    pub fn density(&self) -> f64 {
        let v = self.volume();
        if v > 0.0 {
            self.mass / v
        } else {
            0.0
        }
    }

    pub fn number_density(&self) -> f64 {
        let v = self.volume();
        if v > 0.0 {
            self.particle_count() / v
        } else {
            0.0
        }
    }

    /// Ideal gas + radiation pressure.
    pub fn pressure(&self) -> f64 {
        let gas = self.number_density() * K_B * self.temperature;
        let rad = A_RAD * self.temperature.powi(4) / 3.0;
        gas + rad
    }

    /// Whether [`Matter::pressure`] describes this matter at all.
    ///
    /// `docs/PLAY.md` §7's ninth Phase 2 item, and §3.7's own precedent: a node
    /// crossed by its ensemble reports it rather than doing it quietly.
    ///
    /// `pressure` is an **ideal gas plus radiation**, and there is no other
    /// equation of state in the engine. Handed something condensed it is not
    /// merely inaccurate — a bucket of water prices at 4x10^8 Pa and a solid
    /// handed to SPH bursts from its own pressure before anything else happens.
    /// Measured, on a forty-year-old tree advanced for a twentieth of a second
    /// before §3.3's dispatch kept its members away from the fluid solver: they
    /// reached 1.9x10^8 m/s, 64% of light, across a tree 6.3 m wide.
    ///
    /// D17 is what makes this answerable rather than a guess: a node carrying a
    /// mixture knows its own phase. **Matter nobody has described is not
    /// reported**, because "no information" is not the same answer as "measured
    /// and wrong", and almost every node in a galaxy genuinely is a gas.
    ///
    /// Phase 2 makes it visible. **Water** fixes it, with a liquid and a solid
    /// equation of state — which is `PLAY.md` Phase 5's second piece, named
    /// there in as many words.
    pub fn gas_law_applies(&self) -> bool {
        !self.is_described() || self.mixture.in_phase(crate::chem::Phase::Gas) >= 0.5
    }

    /// Speed at which a disturbance crosses this node's contents, m/s.
    ///
    /// The gas formula, capped by what the internal energy can actually
    /// support. The cap is not a safety clamp — it is the correction for an
    /// assumption the formula makes and does not state: that the internal
    /// energy is heat, held in molecules at `temperature`, colliding often
    /// enough to carry a pressure wave.
    ///
    /// Where that is false the formula answers with the sound speed of a gas
    /// that is not there. A galaxy's internal energy is the motion of its
    /// stars, and `n k T` over that many notional hydrogen atoms came back with
    /// 0.577c — which then set the galaxy's timestep, demanded three thousand
    /// sub-steps a frame, and priced every stellar solve out of the budget.
    ///
    /// [`Matter::velocity_dispersion`] is the same quantity computed from
    /// the energy that is present rather than from a temperature. For a real
    /// gas the two agree within 30% — `c_s = 1.29 sigma` exactly, for an ideal
    /// monatomic one — so the cap costs nothing where the formula is valid and
    /// rescues it where it is not.
    pub fn sound_speed(&self) -> f64 {
        let rho = self.density();
        if rho <= 0.0 {
            return 0.0;
        }
        let gas = (1.6667 * self.pressure() / rho).sqrt();
        gas.min(1.3 * self.velocity_dispersion()).min(C * 0.577)
    }

    /// Free-fall / dynamical time, `1/sqrt(G rho)`. Sets the natural timestep
    /// for a self-gravitating node and hence its scheduling priority.
    pub fn dynamical_time(&self) -> f64 {
        let rho = self.density();
        if rho <= 0.0 {
            return f64::INFINITY;
        }
        1.0 / (G * rho).sqrt()
    }

    /// Moment of inertia about any axis through the centre, kg m^2.
    ///
    /// A uniform sphere, `2/5 M R^2`. Centrally condensed bodies are stiffer
    /// than this — the Earth's is 0.33 rather than 0.4 — but the matter does
    /// not carry a density profile, and inventing one to get a 20% correction
    /// on a quantity used for scheduling would be false precision.
    pub fn moment_of_inertia(&self) -> f64 {
        0.4 * self.mass * self.radius * self.radius
    }

    /// Angular velocity implied by the stored angular momentum, rad/s.
    pub fn angular_velocity(&self) -> Vec3 {
        let i = self.moment_of_inertia();
        if i > 0.0 {
            self.spin.scale(1.0 / i)
        } else {
            Vec3::ZERO
        }
    }

    /// Speed of coherent internal motion, m/s.
    ///
    /// The internal energy of a node covers two different things. Most of it is
    /// thermal: molecules going nowhere in particular, at a temperature the
    /// node already records. The rest is *stirring* — convection, a shock still
    /// crossing, gas that has been hit and has not settled — and only that part
    /// moves the node's contents from one place to another.
    ///
    /// The difference decides how often the node is worth re-solving, so it has
    /// to be drawn somewhere. `thermal_energy` is what the stored temperature
    /// accounts for; anything above it is motion the temperature does not
    /// explain, and is taken to be coherent.
    pub fn stirring_speed(&self) -> f64 {
        if self.mass <= 0.0 {
            return 0.0;
        }
        let excess = self.internal_energy - self.thermal_energy();
        if !(excess > 0.0) {
            return 0.0;
        }
        (2.0 * excess / (3.0 * self.mass)).sqrt().min(C * 0.999)
    }

    /// Speed at which this node's state is actually changing, m/s.
    ///
    /// Three things change it: the node moves, the node turns, and the node's
    /// contents rearrange. Rotation is the term that was missing and it is not
    /// a small correction — Jupiter's equator runs at 12.3 km/s against about
    /// 1 km/s of internal signal, so a cadence blind to it would sample the
    /// planet once every two revolutions and alias it completely.
    ///
    /// What is deliberately *not* here is the sound speed. It is tempting,
    /// because it is the fastest speed anything inside the node is travelling
    /// at, and it is the wrong quantity: a body in thermal equilibrium is not
    /// changing, however fast its molecules are going. Scheduling on the sound
    /// speed asks the engine to re-solve a granite block every 174 microseconds
    /// and a bacterium every 5 nanoseconds, and to produce the same answer
    /// every time. The sound speed still bounds the *timestep* once a node is
    /// being solved — that is [`Matter::signal_crossing`], and it is a
    /// stability limit, not a cadence.
    pub fn characteristic_speed(&self) -> f64 {
        let bulk = if self.mass > 0.0 {
            self.momentum.norm() / self.mass
        } else {
            0.0
        };
        let surface = self.angular_velocity().norm() * self.radius;
        bulk.max(surface).max(self.stirring_speed())
    }

    /// How long before this node's state has changed by `resolution` metres.
    ///
    /// The engine's update cadence. `resolution` is how finely the node is
    /// currently being represented — the size of its materialised children, or
    /// of whatever is looking at it — because that is what decides how far
    /// something may move before the difference is real. A planet treated as a
    /// point may be left alone for half an hour; the same planet resolved into
    /// parcels, or watched at a kilometre, may not.
    pub fn characteristic_time(&self, resolution: f64) -> f64 {
        let v = self.characteristic_speed();
        if !(v > 0.0) || !v.is_finite() {
            return f64::INFINITY;
        }
        resolution.max(0.0) / v
    }

    /// Time for a sound wave to cross one resolution element, in seconds.
    ///
    /// The other limit on a timestep, and the one that has nothing to do with
    /// gravity. A hot, diffuse node has a long dynamical time and a short sound
    /// crossing time, and a solver stepped on the first while ignoring the
    /// second is unstable in exactly the way that looks like physics: parcels
    /// gain energy every step and the node heats up.
    pub fn sound_crossing(&self, parts: usize) -> f64 {
        self.signal_crossing(parts, 0.0)
    }

    /// As [`Matter::sound_crossing`], for a flow that is already moving.
    ///
    /// Information crosses a resolution element at the greater of the sound
    /// speed and the flow speed, and a supersonic flow is the case where the
    /// difference matters. Using the sound speed alone is what lets a node
    /// whose gas is already moving faster than its own sound speed keep taking
    /// steps that the gas has long since outrun.
    pub fn signal_crossing(&self, parts: usize, flow: f64) -> f64 {
        let signal = self.sound_speed().max(flow.abs());
        if !(signal > 0.0) || !signal.is_finite() {
            return f64::INFINITY;
        }
        let h = self.radius / (parts.max(1) as f64).cbrt();
        h / signal
    }

    /// Jeans length: below this, pressure wins and the node will not fragment;
    /// above it, the node *must* be refined or the engine misses collapse.
    pub fn jeans_length(&self) -> f64 {
        let rho = self.density();
        let cs = self.sound_speed();
        if rho <= 0.0 || cs <= 0.0 {
            return f64::INFINITY;
        }
        cs * (std::f64::consts::PI / (G * rho)).sqrt()
    }

    /// Root-mean-square speed of this node's constituents relative to its own
    /// frame, m/s.
    ///
    /// Taken from the greater of the booked internal energy and the thermal
    /// energy the stored temperature implies. The floor matters: a node may be
    /// handed an internal energy of zero and a temperature of 300 K, and it is
    /// still 300 K — its molecules are moving whether or not anybody wrote the
    /// energy down.
    pub fn velocity_dispersion(&self) -> f64 {
        if self.mass <= 0.0 {
            return 0.0;
        }
        let ke = self.internal_energy.max(self.thermal_energy()).max(0.0);
        (2.0 * ke / (3.0 * self.mass)).sqrt().min(C * 0.999)
    }

    /// Nuclear binding energy of the current composition. Bookkeeping only —
    /// see `Composition::nuclear_energy_per_kg`. Differences of this quantity
    /// across a composition change are real energy; the quantity itself is not.
    pub fn nuclear_energy(&self) -> f64 {
        self.mass * self.composition.nuclear_energy_per_kg()
    }

    /// Sackur-Tetrode-ish ideal gas entropy. Absolute value is not meaningful
    /// at this level of modelling; *differences* are, and the engine only ever
    /// uses differences (to assert coarse-graining never decreases entropy).
    pub fn estimate_entropy(&self) -> f64 {
        let n = self.particle_count();
        if n <= 0.0 || self.temperature <= 0.0 {
            return 0.0;
        }
        let v = self.volume().max(1e-300);
        let mu = self.composition.mean_molecular_mass(self.temperature);
        let lambda = H_PLANCK / (2.0 * std::f64::consts::PI * mu * K_B * self.temperature).sqrt();
        let arg = (v / (n * lambda.powi(3))).max(1e-300);
        n * K_B * (arg.ln() + 2.5)
    }

    /// Local entropy plus everything exported. This is the quantity the second
    /// law constrains, and the only one worth asserting monotonicity on.
    pub fn total_entropy(&self) -> f64 {
        self.entropy + self.entropy_exported
    }

    /// Total energy, including rest mass. The `Conserved.energy` slot.
    ///
    /// ```text
    ///   E = M c^2  +  K_bulk(M, P)  +  U_internal  +  Phi_binding
    /// ```
    ///
    /// The bulk term is the *exact* relativistic one, not `p^2/2M`, so that the
    /// decomposition is invertible to the last bit: given `(M, P, E, Phi)` you
    /// recover `U` exactly, which is what `summarise` does. A Newtonian bulk
    /// term makes the round trip lossy at the 10^-5 level for anything moving
    /// at galactic-rotation speeds, which is enough to be visible as energy
    /// drift when a user pans across a disk.
    pub fn total_energy(&self) -> f64 {
        self.mass * C2 + self.non_rest_energy()
    }

    /// Total energy *excluding* rest mass.
    ///
    /// Rest mass exceeds every other term by roughly 10^16 for ordinary matter,
    /// so any process that moves energy around — growth, heating, radiation,
    /// construction — is completely invisible in a difference of
    /// `total_energy()` at double precision. A decade of a tree's growth is
    /// 10^10 J against a rest mass of 10^19 J: differencing the totals leaves
    /// about seven significant digits of the answer and none of them reliable.
    ///
    /// So anything auditing an energy *flow* must difference this instead. It
    /// is the same lesson as measuring conservation against natural scales
    /// rather than net values (see `Conserved::error_against`): when two terms
    /// differ by more orders of magnitude than the arithmetic carries, the
    /// small one has to be tracked on its own.
    pub fn non_rest_energy(&self) -> f64 {
        bulk_kinetic(self.mass, self.momentum)
            + self.internal_energy
            + self.gravitational_binding
            + self.cohesive_binding
            + self.external_potential
            + self.chemical_energy
    }

    /// Orbital angular momentum contributed by this node sitting at `offset`
    /// and moving with `momentum` in its parent's frame, plus its own spin.
    pub fn angular_momentum_about(&self, offset: Vec3) -> Vec3 {
        offset.cross(self.momentum) + self.spin
    }

    pub fn conserved(&self) -> Conserved {
        Conserved {
            energy: self.total_energy(),
            momentum: self.momentum,
            angular_momentum: self.spin,
            charge: self.charge,
            baryon: self.baryon_number,
            lepton: self.lepton_number,
        }
    }

    /// Fraction of this matter's mass that is in a solid phase.
    ///
    /// `docs/PLAY.md` D13 says a node's boundary is derived from *measured
    /// solidity* rather than from what generated the node, and left the
    /// measurement unspecified because there was nothing to measure it from.
    /// D17 supplies it: a mixture on the matter, and phase against a melting
    /// point `chem::analyse` derives, so this is a **reading** rather than a
    /// law still to be written.
    ///
    /// Zero for matter nobody has described, which is the honest answer to no
    /// information and not a claim that it is a gas. A caller deciding whether
    /// something has a surface has to distinguish "measured as not solid" from
    /// "not described", and [`Matter::is_described`] is how.
    pub fn solid_fraction(&self) -> f64 {
        self.mixture.in_phase(crate::chem::Phase::Solid)
    }

    /// Whether anything has said what this matter is made of.
    ///
    /// The difference between "measured and found not solid" and "never
    /// measured", which every reading off the mixture has to keep apart.
    pub fn is_described(&self) -> bool {
        !self.mixture.is_empty()
    }

    /// Set the temperature and rebalance internal energy to match. Used when a
    /// solver decides a node has heated up.
    pub fn set_temperature(&mut self, t: f64) {
        self.temperature = t.max(2.725);
        self.internal_energy = self.thermal_energy();
        self.entropy = self.estimate_entropy();
    }

    /// Add heat, letting the temperature follow. Returns the new temperature.
    pub fn add_heat(&mut self, joules: f64) -> f64 {
        let n = self.particle_count();
        if n <= 0.0 {
            self.internal_energy += joules;
            return self.temperature;
        }
        let dt = joules / (1.5 * n * K_B);
        self.internal_energy += joules;
        self.temperature = (self.temperature + dt).max(2.725);
        self.entropy = self.estimate_entropy();
        self.temperature
    }

    pub fn is_finite(&self) -> bool {
        self.mass.is_finite()
            && self.com.is_finite()
            && self.momentum.is_finite()
            && self.spin.is_finite()
            && self.internal_energy.is_finite()
            && self.radius.is_finite()
            && self.temperature.is_finite()
    }
}

/// One materialised body: the fine-grained representation.
///
/// Deliberately compact. On the GPU this is 4 x `vec4<f32>` for the hot fields
/// (position, velocity, mass+radius+charge+flags) with the cold fields in a
/// parallel array, which is the layout the bandwidth budget in
/// `docs/PERFORMANCE.md` is computed against.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Body {
    /// Position in the *parent node's* frame, metres.
    pub pos: Vec3,
    pub vel: Vec3,
    pub mass: f64,
    pub radius: f64,
    pub charge: f64,
    /// Internal (thermal + binding) energy carried by this body.
    pub internal_energy: f64,
    pub temperature: f64,
    pub composition: Composition,
    pub spin: Vec3,
    /// Which way this thing is facing, in its parent's frame.
    ///
    /// `spin` says how fast it is turning and about what axis; this says where
    /// it has turned *to*, and the two are as different as velocity and
    /// position. Nothing in the engine had it: a body was a point with a size,
    /// so a plank and a boulder of the same mass were the same object seen from
    /// every direction, and `promote` had to hand its child `Quat::IDENTITY`
    /// because there was nothing else to hand it.
    ///
    /// Identity for anything sampled. A gas parcel has no orientation worth
    /// carrying and giving one a random draw would change nothing except every
    /// bit-exactness test in the suite. It is stated by whatever *knows* the
    /// answer — a recipe's part, a promoted child coming back — which is the
    /// same rule the rest of `Body` follows.
    pub orientation: Quat,
    /// Half-extents of this thing's own box, along its own axes, metres.
    ///
    /// **`Vec3::ZERO` means "a sphere of `radius`"**, which is what everything
    /// sampled is and what `radius` alone has always described. Non-zero means
    /// a box, and then `radius` is that box's bounding radius — derived, stored
    /// alongside, and written by the same constructor so the two cannot
    /// disagree. See [`Body::solid`].
    ///
    /// A sphere cannot be written as half-extents (`(r,r,r)` is a cube whose
    /// bounding radius is `r*sqrt(3)`), which is why the zero case carries the
    /// discriminator rather than a third field or an enum.
    pub half: Vec3,
    /// Which substance this body is one of, or `UNSPECIATED`.
    ///
    /// D17 decided a body carries no *mixture* of its own, and that stands: a
    /// `Mixture` is 200 bytes on a 184-byte struct and every body of a sampled
    /// node is drawn from one matter anyway. A single id is four bytes and is a
    /// different question — *which* of the node's pools this particular thing
    /// is made of — which only something with named parts can answer and which
    /// a box of oak panels on a steel frame has to.
    pub substance: crate::chem::SubstanceId,
    /// Index of this body within its parent's materialised set. Also the index
    /// into the parent's random streams, which is what makes regeneration
    /// order-independent.
    pub slot: u32,
    pub kind: BodyKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodyKind {
    /// A statistical stand-in for many real objects (dark matter, ISM parcel,
    /// star cluster). Its mass is not the mass of one thing.
    Super,
    Star,
    CompactObject,
    Planet,
    GasParcel,
    Grain,
    Molecule,
    Atom,
    Nucleus,
    Nucleon,
    Electron,
    Photon,
}

impl Default for Body {
    fn default() -> Self {
        Body {
            pos: Vec3::ZERO,
            vel: Vec3::ZERO,
            mass: 0.0,
            radius: 0.0,
            charge: 0.0,
            internal_energy: 0.0,
            temperature: 2.725,
            composition: Composition::primordial(),
            spin: Vec3::ZERO,
            orientation: Quat::IDENTITY,
            half: Vec3::ZERO,
            substance: crate::chem::SubstanceId::UNSPECIATED,
            slot: 0,
            kind: BodyKind::Super,
        }
    }
}

/// Rounding on a solid body's edges, as a fraction of its smallest half-extent.
///
/// Not needed for correctness — a sharp box is convex — but [`crate::shape::Hull::slab`]
/// builds the box out of corner spheres, and a zero radius collapses all eight
/// onto the corners: the support function is still the box's, but a contact
/// near an edge gets a corner's normal rather than a face's. A tenth of the
/// thinnest dimension is a plank's arris.
const EDGE_ROUNDING: f64 = 0.1;

impl Body {
    /// A solid piece of stated size, substance and orientation — what a recipe
    /// names and what `sample` then fills in the rest of.
    ///
    /// `docs/PLAY.md` D15's part *is* one of these, which is the point: there is
    /// one type for a thing inside a node and not one for things that were
    /// sampled and another for things that were made. Two types would be the
    /// engine knowing how a thing came to be there, which is the same failure
    /// as "only a built thing has a surface" one layer up.
    ///
    /// `radius` is set from `half` here and nowhere else, so the bounding
    /// radius and the extents cannot drift apart.
    pub fn solid(
        pos: Vec3,
        half: Vec3,
        mass: f64,
        substance: crate::chem::SubstanceId,
        slot: u32,
    ) -> Body {
        let half = Vec3 { x: half.x.abs(), y: half.y.abs(), z: half.z.abs() };
        Body {
            pos,
            mass: mass.max(0.0),
            radius: half.norm(),
            half,
            substance,
            slot,
            kind: BodyKind::Grain,
            ..Body::default()
        }
    }

    /// Is this body a box rather than a sphere?
    #[inline]
    pub fn is_boxed(&self) -> bool {
        self.half != Vec3::ZERO
    }

    /// Half the longest diagonal — how far this body reaches from its own
    /// centre, whichever shape it is.
    #[inline]
    pub fn reach(&self) -> f64 {
        if self.is_boxed() {
            self.half.norm()
        } else {
            self.radius
        }
    }

    /// The volume the body occupies, m³.
    pub fn volume(&self) -> f64 {
        if self.is_boxed() {
            8.0 * self.half.x * self.half.y * self.half.z
        } else {
            4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3)
        }
    }

    /// The filled convex solid this body presents, in its parent's frame.
    ///
    /// One place that turns a body into geometry, so a caller never has to ask
    /// what kind of body it is holding. A sphere is a sphere; a box is a
    /// [`crate::shape::Hull::slab`] turned by the body's own orientation.
    pub fn hull(&self) -> crate::shape::Hull {
        if !self.is_boxed() {
            return crate::shape::Hull::sphere(self.pos, self.radius);
        }
        let thinnest = self.half.x.min(self.half.y).min(self.half.z);
        crate::shape::Hull::slab(Vec3::ZERO, self.half, thinnest * EDGE_ROUNDING)
            .placed(self.orientation, self.pos)
    }

    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.vel.norm2()
    }

    pub fn momentum(&self) -> Vec3 {
        // Relativistically correct: p = gamma m v.
        self.vel.scale(self.mass * crate::coords::gamma(self.vel))
    }

    pub fn angular_momentum(&self) -> Vec3 {
        self.pos.cross(self.momentum()) + self.spin
    }

    /// Constituent particle count, by the same route [`Matter::particle_count`]
    /// takes. One derivation, used at both resolutions: a body and the matter
    /// it summarises into must agree about how many things they are made of, or
    /// `summarise(sample(m)) == m` would not hold on energy.
    pub fn particle_count(&self) -> f64 {
        let mu = self.composition.mean_molecular_mass(self.temperature);
        if mu > 0.0 {
            self.mass / mu
        } else {
            0.0
        }
    }

    /// Joules per kelvin. See [`Matter::heat_capacity`].
    pub fn heat_capacity(&self) -> f64 {
        1.5 * self.particle_count() * crate::units::K_B
    }

    /// Put heat into this body, moving its temperature with it.
    ///
    /// The same law as [`Matter::add_heat`], written once for each resolution
    /// because the two hold their state in different structs and not because
    /// the physics differs. Negative joules take heat out, which is half of
    /// what any exchange needs.
    ///
    /// The 2.725 K floor is the microwave background: nothing in the universe
    /// is colder, and an exchange that ran a body below it would be taking heat
    /// from somewhere that has none to give. `internal_energy` is *not* floored
    /// with it — the joules are the conserved quantity and they are recorded
    /// exactly, so the pair stays auditable at the point where the floor bites.
    pub fn add_heat(&mut self, joules: f64) -> f64 {
        let n = self.particle_count();
        self.internal_energy += joules;
        if n > 0.0 {
            let dt = joules / (1.5 * n * crate::units::K_B);
            self.temperature = (self.temperature + dt).max(2.725);
        }
        self.temperature
    }
}

/// Reduce a materialised set back to `Matter`.
///
/// This is the *summarising* operator R. Together with sampling P it must
/// satisfy `R(P(s)) = s` on the conserved set — the property that lets the
/// engine throw detail away safely. See `tests/consistency.rs`.
pub fn summarise(bodies: &[Body], mutual_potential: f64) -> Matter {
    if bodies.is_empty() {
        return Matter::default();
    }
    let n = bodies.len();
    let mass = det_sum_by(n, &|i| bodies[i].mass);
    if mass <= 0.0 {
        return Matter::default();
    }
    let com = det_sum_v3_by(n, &|i| bodies[i].pos.scale(bodies[i].mass)).scale(1.0 / mass);
    let momentum = total_momentum(bodies);

    // Angular momentum about the centre of mass. The mass-weighted share of the
    // bulk momentum contributes nothing (`sum m_i (r_i - com) = 0` by
    // definition of `com`), so this single expression is both the spin and the
    // total angular momentum about the com.
    let spin = total_spin(bodies, com);

    // Energy, exactly. `kinetic_energy_of` is relativistic, so this is valid for
    // a 20 K gas parcel and for a 10 GeV cosmic ray in the same expression.
    let e_kin = kinetic_energy_of(bodies);
    let child_internal = det_sum_by(n, &|i| bodies[i].internal_energy);
    let internal = e_kin + child_internal - bulk_kinetic(mass, momentum);

    let radius = {
        let r2 = det_sum_by(n, &|i| {
            let b = &bodies[i];
            b.mass * (b.pos - com).norm2()
        }) / mass;
        (r2.max(0.0).sqrt() * RMS_TO_RADIUS).max(1e-30)
    };

    let mut comp = [0.0; COARSE_ELEMENTS];
    for s in 0..COARSE_ELEMENTS {
        comp[s] = det_sum_by(n, &|i| bodies[i].mass * bodies[i].composition.0[s]) / mass;
    }
    let composition = Composition(comp).normalised();
    let charge = det_sum_by(n, &|i| bodies[i].charge);

    let mut matter = Matter {
        mass,
        com,
        momentum,
        spin,
        internal_energy: internal,
        gravitational_binding: mutual_potential,
        // Not knowable from the children alone; the caller reinstates these.
        // `cohesive_binding` joins the list for the reason its own doc gives:
        // the bonds are far below the scale of the bodies, so measuring where
        // the bodies ended up says nothing about them.
        cohesive_binding: 0.0,
        external_potential: 0.0,
        chemical_energy: 0.0,
        entropy_exported: 0.0,
        radius,
        temperature: 2.725,
        composition,
        charge,
        baryon_number: mass * composition.nucleons_per_kg(),
        lepton_number: 0.0,
        entropy: 0.0,
        magnetic_energy: 0.0,
        luminosity: det_sum_by(n, &|i| match bodies[i].kind {
            BodyKind::Star | BodyKind::CompactObject => {
                stefan_boltzmann(bodies[i].radius, bodies[i].temperature)
            }
            _ => 0.0,
        }),
        // A `Body` carries no speciation, deliberately: a `Mixture` is 200
        // bytes against a body's 184, and `docs/PLAY.md` §5A.5's own costing
        // puts one on a `Matter` and not on a body. So `summarise` has nothing
        // to blend *from the bodies*, and reports none — the caller reinstates
        // it, exactly as it does for `chemical_energy` and `external_potential`.
        //
        // Nothing is lost by that, because `sample` distributes one mixture to
        // every body it makes: a set of bodies drawn from one matter is made of
        // what that matter was made of, and blending identical mixtures is the
        // identity. The case that *does* need a blend is a node whose slots
        // hold promoted children with chemistries of their own, and those are
        // nodes rather than bodies — `Tree::coarsen` blends them with
        // `Mixture::blend` before this is called.
        mixture: crate::chem::Mixture::new(),
    };
    matter.lepton_number =
        matter.baryon_number * composition.electrons_per_nucleon() - charge / E_CHARGE;

    // Temperature is *derived* from the random kinetic energy, never averaged
    // from the children: averaging temperatures across unequal masses is wrong,
    // and across this dynamic range it is wrong by orders of magnitude.
    matter.temperature = temperature_of(mass, composition, internal);
    matter.entropy = matter.estimate_entropy();
    matter
}

/// The temperature a mass of this composition has when it holds this much
/// internal energy — the inverse of `Matter::thermal_energy`.
///
/// **One convention, in one place.** Energy is the conserved quantity and
/// temperature is derived from it, so the derivation has to give the same
/// answer at every resolution or the two disagree across a scale transform.
/// It did: a body left `sampler::close_books` carrying an internal energy
/// written to close the books exactly and a temperature copied from its parent
/// before that energy was known, and the two were out by a factor of 1.3974 on
/// a self-gravitating planet — the body said 1739.50 K while holding the energy
/// of 2430.79 K, and the next `summarise` believed the energy and stepped the
/// node's temperature by 40%.
///
/// `mean_molecular_mass` is evaluated at `1e4` and not at the answer, which is
/// deliberate and is the convention `summarise` has always used: the function
/// is a step at 10^4 K (`ionised = temperature > 1.0e4`), so evaluating it at
/// the temperature being solved for makes this implicit for no gain below the
/// step. **Above the step it is wrong**, and equally wrong at both resolutions,
/// which is the property worth having — one answer, and its error reported
/// rather than differing by resolution. See `docs/BACKLOG.md`.
pub fn temperature_of(mass: f64, composition: Composition, internal: f64) -> f64 {
    let mu = composition.mean_molecular_mass(1e4);
    let particles = if mu > 0.0 { mass / mu } else { 0.0 };
    if particles > 0.0 {
        (2.0 * internal.max(0.0) / (3.0 * particles * K_B)).max(2.725)
    } else {
        2.725
    }
}

/// RMS radius to equivalent-uniform-sphere radius: a uniform ball of radius R
/// has `<r^2> = 3R^2/5`, so `R = sqrt(5/3) * rms`.
pub const RMS_TO_RADIUS: f64 = 1.290_994_448_735_805_6;

/// Exact relativistic momentum of a materialised set.
pub fn total_momentum(bodies: &[Body]) -> Vec3 {
    det_sum_v3_by(bodies.len(), &|i| bodies[i].momentum())
}

/// Exact angular momentum about `centre`, including intrinsic spins.
pub fn total_spin(bodies: &[Body], centre: Vec3) -> Vec3 {
    det_sum_v3_by(bodies.len(), &|i| {
        let b = &bodies[i];
        (b.pos - centre).cross(b.momentum()) + b.spin
    })
}

/// Exact relativistic kinetic energy `sum (gamma - 1) m c^2`.
pub fn kinetic_energy_of(bodies: &[Body]) -> f64 {
    det_sum_by(bodies.len(), &|i| {
        let b = &bodies[i];
        (crate::coords::gamma(b.vel) - 1.0) * b.mass * C2
    })
}

/// Kinetic energy of a body of rest mass `m` carrying total momentum `p`,
/// exactly: `sqrt((mc^2)^2 + (pc)^2) - mc^2`. Reduces to `p^2/2m` for small p
/// but stays correct — and, crucially, stays *invertible* — for large p.
pub fn bulk_kinetic(mass: f64, momentum: Vec3) -> f64 {
    if mass <= 0.0 {
        return momentum.norm() * C;
    }
    let mc2 = mass * C2;
    let pc = momentum.norm() * C;
    // Numerically stable form: for pc << mc2 the naive difference of two nearly
    // equal huge numbers loses every significant digit.
    let x = pc / mc2;
    if x < 1e-4 {
        mc2 * x * x * (0.5 - x * x / 8.0)
    } else {
        (mc2 * mc2 + pc * pc).sqrt() - mc2
    }
}

pub fn stefan_boltzmann(radius: f64, temperature: f64) -> f64 {
    4.0 * std::f64::consts::PI * radius * radius * SIGMA_SB * temperature.powi(4)
}

/// Inertia tensor of a materialised set about `centre`.
/// How far a set of contents actually extends from its own centre.
///
/// A node's radius is fixed when the node is created and nothing has ever
/// checked that it still describes what the node holds. Contents with a
/// positive velocity divergence — a dispersing cloud, debris, an explosion —
/// either outgrow it or leave it, and every consumer of node radius is then
/// wrong by the same factor: SPH's smoothing length `radius / count^(1/3)`,
/// the gravity softening, the LOD's angular size, the volume query's node
/// selection, the neighbour grid's spacing. None of them fail loudly. They all
/// quietly describe a neighbourhood that has stopped existing.
///
/// This is the measurement that notices. `docs/PLAY.md` Phase 1 asks for it on
/// its own — three separate things need it and it is to be built once: node
/// splitting, D6's patch handoff, and promoting a fragment that has left its
/// parent. It is also the detector `BACKLOG.md` wanted for the class of fault
/// where a solver flings a node's bodies across twenty orders of magnitude: it
/// fires on the first frame, in the node that caused it, rather than twenty
/// tiers away inside a hash function.
///
/// # What is measured, and what is not
///
/// Mass-weighted RMS distance and the furthest surface. **Not** the principal
/// axes of the second moment, which `BACKLOG.md` offers as the alternative.
/// The axes are what *splitting* needs — they say which way to cut — and
/// splitting is not in Phase 1; none of the three consumers above can use them.
/// Writing a symmetric eigensolver for a caller that does not exist is the
/// thing this project's backlog discipline exists to prevent, so the axes go
/// with the split that wants them.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Spread {
    /// Centre of mass of the contents, in the node's own frame. Not assumed to
    /// be the origin: a node whose contents have drifted off-centre is exactly
    /// one of the states worth noticing.
    pub centre: Vec3,
    /// Mass-weighted RMS distance from `centre`. The characteristic size of the
    /// contents, and robust to one outlier in a way `furthest` is not.
    pub rms: f64,
    /// Distance from `centre` to the furthest occupant's *surface*, not its
    /// centre. A body is not inside a volume its own bulk sticks out of.
    pub furthest: f64,
    pub count: usize,
}

impl Spread {
    /// How much of what the node claims its contents actually occupy.
    ///
    /// One means they exactly fill it. Below one they rattle around inside it.
    /// Above one they have outgrown it, and every length derived from the
    /// radius is wrong by this factor.
    pub fn occupancy(&self, radius: f64) -> f64 {
        if radius > 0.0 && self.furthest.is_finite() {
            self.furthest / radius
        } else {
            0.0
        }
    }

    /// Measure a set of positions carrying a mass and a size each.
    ///
    /// Taken as an iterator rather than a `&[Body]` because a node's contents
    /// are its bodies *and* its promoted children, and those are the same thing
    /// at two resolutions. `Tree::spread` supplies both through one call.
    pub fn of(parts: impl IntoIterator<Item = (Vec3, f64, f64)>) -> Spread {
        let parts: Vec<(Vec3, f64, f64)> = parts.into_iter().collect();
        let count = parts.len();
        if count == 0 {
            return Spread::default();
        }
        let mass = crate::math::det_sum_by(count, &|i| parts[i].1);
        let centre = if mass > 0.0 {
            crate::math::det_sum_v3_by(count, &|i| parts[i].0.scale(parts[i].1)).scale(1.0 / mass)
        } else {
            // Massless contents still have a geometry, and answering with the
            // origin would report a spread about a point nothing is near.
            crate::math::det_sum_v3_by(count, &|i| parts[i].0).scale(1.0 / count as f64)
        };
        let rms = if mass > 0.0 {
            (crate::math::det_sum_by(count, &|i| parts[i].1 * (parts[i].0 - centre).norm2()) / mass).max(0.0).sqrt()
        } else {
            (crate::math::det_sum_by(count, &|i| (parts[i].0 - centre).norm2()) / count as f64).max(0.0).sqrt()
        };
        let mut furthest = 0.0f64;
        for (p, _, r) in &parts {
            furthest = furthest.max((*p - centre).norm() + r.max(0.0));
        }
        Spread { centre, rms, furthest, count }
    }
}

pub fn inertia_tensor(bodies: &[Body], centre: Vec3) -> Mat3 {
    let mut m = Mat3::zero();
    for b in bodies {
        let r = b.pos - centre;
        let r2 = r.norm2();
        let outer = r.outer(r);
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                m.0[i][j] += b.mass * (r2 * delta - outer.0[i][j]);
            }
        }
    }
    m
}

/// Mutual gravitational potential energy of a materialised set — direct O(n^2),
/// used for exactness checks and for small sets. Large sets get the tree
/// estimate from `solvers::gravity`.
pub fn mutual_gravitational_energy(bodies: &[Body], softening: f64) -> f64 {
    let n = bodies.len();
    let mut terms = Vec::with_capacity(n * (n.saturating_sub(1)) / 2 + 1);
    for i in 0..n {
        for j in (i + 1)..n {
            let d = (bodies[i].pos - bodies[j].pos).norm2() + softening * softening;
            terms.push(-G * bodies[i].mass * bodies[j].mass / d.sqrt());
        }
    }
    det_sum(&terms)
}

use crate::math::{det_sum_by, det_sum_v3_by};
