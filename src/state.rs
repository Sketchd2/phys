//! Aggregate state: what the engine stores when it is *not* storing particles.
//!
//! A node holds a bulk description of its contents. Everything below its own
//! resolution is absent — not approximated, absent — and is regenerated on
//! demand by `prolong.rs`. For that to be legitimate, the bulk description must
//! carry every quantity that the missing detail is *not allowed to change*:
//! the conserved set. If refinement and re-coarsening return exactly the same
//! conserved tuple, no experiment performed at the coarse scale can tell
//! whether the detail was ever there.
//!
//! That is the engine's central correctness claim, and `Conserved` is the
//! object it is stated in terms of.

use crate::math::{det_sum, Mat3, Vec3};
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

/// Mass fractions by species. Always sums to 1 for a non-empty node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Composition(pub [f64; NSPECIES]);

impl Default for Composition {
    fn default() -> Self {
        Composition::primordial()
    }
}

impl Composition {
    /// Big Bang nucleosynthesis output: the initial condition for gas that has
    /// never been through a star.
    pub fn primordial() -> Composition {
        let mut c = [0.0; NSPECIES];
        c[Species::Hydrogen as usize] = 0.75;
        c[Species::Helium as usize] = 0.25;
        Composition(c)
    }

    /// Roughly solar (Asplund 2009 mass fractions, lumped into our buckets).
    pub fn solar() -> Composition {
        let mut c = [0.0; NSPECIES];
        c[Species::Hydrogen as usize] = 0.7381;
        c[Species::Helium as usize] = 0.2485;
        c[Species::Carbon as usize] = 0.0024;
        c[Species::Nitrogen as usize] = 0.0007;
        c[Species::Oxygen as usize] = 0.0057;
        c[Species::Silicon as usize] = 0.0007;
        c[Species::Iron as usize] = 0.0013;
        c[Species::Other as usize] = 0.0026;
        Composition(c).normalised()
    }

    /// Pure one species — used when the user drills into a specific atom.
    pub fn pure(s: Species) -> Composition {
        let mut c = [0.0; NSPECIES];
        c[s as usize] = 1.0;
        Composition(c)
    }

    /// Cold dark matter: gravitationally active, chemically inert. Lives in the
    /// `Other` bucket but is flagged separately by the tier's solver.
    pub fn dark() -> Composition {
        Composition::pure(Species::Other)
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

    pub fn get(&self, s: Species) -> f64 {
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
        let mut c = [0.0; NSPECIES];
        for i in 0..NSPECIES {
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
        for s in Species::ALL {
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
        for s in Species::ALL {
            n += self.get(s) * s.z() / s.a();
        }
        n
    }

    /// Nucleons per kilogram.
    pub fn nucleons_per_kg(&self) -> f64 {
        let mut n = 0.0;
        for s in Species::ALL {
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
        for s in Species::ALL {
            let x = self.get(s);
            if x <= 0.0 {
                continue;
            }
            // per kg of this species: (A nucleons / mass) * B/A
            e -= x * (s.a() / s.mass_kg()) * s.binding_per_nucleon_mev() * MEV;
        }
        e
    }
}

/// Everything a node knows about itself without materialising its children.
///
/// Roughly 200 bytes. A galaxy's worth of these (a few million live nodes) is
/// well under a gigabyte, which is what makes the whole approach fit in a
/// 6 GB card alongside the materialised working set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aggregate {
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
    /// Self-gravitational (and, at fine tiers, electrostatic) binding energy.
    /// Negative for a bound object. Tracked explicitly so that coarsening does
    /// not silently destroy it.
    pub binding_energy: f64,
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
}

impl Default for Aggregate {
    fn default() -> Self {
        Aggregate {
            mass: 0.0,
            com: Vec3::ZERO,
            momentum: Vec3::ZERO,
            spin: Vec3::ZERO,
            internal_energy: 0.0,
            binding_energy: 0.0,
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
        }
    }
}

impl Aggregate {
    /// Neutral matter of the given composition: charge zero, and baryon/lepton
    /// numbers implied by the mass. This is the normal way to create matter.
    pub fn neutral(mass: f64, radius: f64, temperature: f64, composition: Composition) -> Aggregate {
        let nucleons = mass * composition.nucleons_per_kg();
        let electrons = nucleons * composition.electrons_per_nucleon();
        let mut a = Aggregate {
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
    pub fn with_charge(mut self, charge: f64) -> Aggregate {
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

    pub fn sound_speed(&self) -> f64 {
        let rho = self.density();
        if rho <= 0.0 {
            return 0.0;
        }
        (1.6667 * self.pressure() / rho).sqrt().min(C * 0.577)
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

    /// Velocity dispersion implied by the internal energy.
    pub fn velocity_dispersion(&self) -> f64 {
        if self.mass <= 0.0 {
            return 0.0;
        }
        let ke = self.internal_energy.max(0.0);
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
    /// recover `U` exactly, which is what `restrict` does. A Newtonian bulk
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
            + self.binding_energy
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
            slot: 0,
            kind: BodyKind::Super,
        }
    }
}

impl Body {
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
}

/// Reduce a materialised set back to a bulk description.
///
/// This is the *restriction* operator R. Together with prolongation P it must
/// satisfy `R(P(s)) = s` on the conserved set — the property that lets the
/// engine throw detail away safely. See `tests/consistency.rs`.
pub fn restrict(bodies: &[Body], mutual_potential: f64) -> Aggregate {
    if bodies.is_empty() {
        return Aggregate::default();
    }
    let n = bodies.len();
    let mass = det_sum_by(n, &|i| bodies[i].mass);
    if mass <= 0.0 {
        return Aggregate::default();
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

    let mut comp = [0.0; NSPECIES];
    for s in 0..NSPECIES {
        comp[s] = det_sum_by(n, &|i| bodies[i].mass * bodies[i].composition.0[s]) / mass;
    }
    let composition = Composition(comp).normalised();
    let charge = det_sum_by(n, &|i| bodies[i].charge);

    let mut agg = Aggregate {
        mass,
        com,
        momentum,
        spin,
        internal_energy: internal,
        binding_energy: mutual_potential,
        // Not knowable from the children alone; the caller reinstates these.
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
    };
    agg.lepton_number =
        agg.baryon_number * composition.electrons_per_nucleon() - charge / E_CHARGE;

    // Temperature is *derived* from the random kinetic energy, never averaged
    // from the children: averaging temperatures across unequal masses is wrong,
    // and across this dynamic range it is wrong by orders of magnitude.
    let mu = composition.mean_molecular_mass(1e4);
    let particles = if mu > 0.0 { mass / mu } else { 0.0 };
    agg.temperature = if particles > 0.0 {
        (2.0 * internal.max(0.0) / (3.0 * particles * K_B)).max(2.725)
    } else {
        2.725
    };
    agg.entropy = agg.estimate_entropy();
    agg
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
