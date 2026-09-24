//! Morphology: matter that has a history rather than a temperature.
//!
//! # Why the rest of the engine cannot represent a tree
//!
//! Everything `sample` regenerates is *ergodic*. One max-entropy sample of a
//! gas cloud is as good as another, because no observation can tell them apart,
//! and that interchangeability is what licenses throwing the detail away.
//!
//! A tree is not like that. Its branch structure is low-entropy and
//! historically contingent — not a typical sample from anything, but the
//! specific record of which branch got shaded in year three. And the difference
//! is *observable*: someone who saw the tree yesterday will notice if they are
//! handed a different one. The conserved tuple, which pins mass and momentum
//! exactly, pins none of what actually matters here.
//!
//! # What replaces it
//!
//! For structured matter the generator changes from "sample the max-entropy
//! distribution" to "run a developmental program to a given age". The stored
//! state is then not a million vertices but a few dozen bytes:
//!
//! ```text
//!     structure = program(genome, age, events)
//! ```
//!
//! — a pure function, addressed the same way as everything else in `rng.rs`, so
//! it regenerates bit-for-bit. A forest of 10^9 trees costs nothing until
//! someone walks into it.
//!
//! # Why this fits the architecture better than the physics does
//!
//! Growth runs on the node's *matter*, never on the fine structure. A forest does
//! not grow by integrating 10^9 trees; it grows by advancing one ordinary
//! differential equation on a forest node, at O(1) per node. That is cheap
//! enough to run on the entire world every frame, coarse or not — so the
//! laziness that the rest of the engine works hard for is simply free here.
//!
//! # Growth is a transaction, not an exemption
//!
//! Building order out of disorder costs free energy and exports entropy. Every
//! step returns a [`GrowthStep`] that has to balance before it is applied:
//! energy in equals energy stored plus heat released, and local entropy plus
//! exported entropy is non-negative. A program that tried to grow too
//! efficiently would be rejected by `GrowthStep::validate` rather than
//! quietly minting free energy.

use crate::math::{v3, Vec3};
use crate::rng::{Purpose, Stream};
use crate::state::{BodyKind, Composition};
use crate::units::*;

/// Free energy density of dry biomass, J/kg. Wood is ~17-19 MJ/kg.
pub const BIOMASS_ENERGY: f64 = 17.0e6;
/// Embodied energy of reinforced concrete and steel construction, J/kg.
pub const CONSTRUCTION_ENERGY: f64 = 2.5e6;
/// Bulk density of a framed building including voids, kg/m^3.
///
/// **The last of the tabulated densities, and it is not a material's.** A
/// framed building is mostly air: its walls, floors and frame are dense and the
/// rooms between them are not, so this is a statement about an arrangement
/// rather than about a substance — which is exactly the distinction D11 draws,
/// and exactly why it cannot be derived from a substance. What would derive it
/// is the recipe's own parts against the volume they enclose, and that is
/// circular while the recipe is what turns the mass into the volume. See
/// `World::rewrite_recipe` for the measurement.
pub const BUILDING_DENSITY: f64 = 250.0;

/// Effective conversion of incident radiation into stored biomass.
///
/// Not the photosynthetic quantum efficiency (~3%), which is a laboratory
/// number for a leaf. This is the ecosystem figure: a temperate forest fixes
/// roughly 1.2 kg of dry matter per square metre per year under ~200 W/m^2 of
/// mean insolation, which is 0.32% of the incident energy. Using the leaf
/// number instead makes trees grow about a hundred times too fast, which looks
/// plausible for a few frames and absurd after a simulated decade.
pub const PHOTOSYNTHETIC_YIELD: f64 = 0.0032;

/// Fraction of stored free energy that must be paid as ordering entropy.
///
/// A crude stand-in for the entropy of polymerisation: assembling monomers into
/// an ordered polymer lowers the configurational entropy roughly in proportion
/// to the free energy stored. The exact coefficient does not matter much — what
/// matters is that it is non-zero, so the second-law check has something to
/// bite on.
pub const ORDERING_FRACTION: f64 = 0.3;

/// What a structure is made of and how it is put together.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Program {
    /// Recursive branching under elastic self-similarity. Light-limited.
    Tree,
    /// Radial branching, carbonate skeleton. Limited by dissolved mineral.
    Coral,
    /// Planned: a framed tower built floor by floor.
    Tower,
    /// Planned: a wall laid course by course.
    Wall,
    /// A patch of ground. Columns of bedrock standing on nothing, with a
    /// surface height that varies across the patch.
    ///
    /// Terrain is a program for the same reason a tree is: what a node holds is
    /// a rule and a seed, and the geometry is derived on demand and thrown away
    /// again. It differs from the living programs in that it does not grow —
    /// its `advance` weathers rather than builds — and from the planned ones in
    /// that nobody is constructing it.
    Terrain,
    /// A settlement. Plots on a street grid, each one a footprint that a
    /// building can be promoted out of.
    ///
    /// The level *between* terrain and a building: it does not describe walls
    /// and floors, it describes where the buildings are. Promoting one of its
    /// plots gives a node that carries [`Program::Tower`] and builds itself.
    Settlement,
}

impl Program {
    pub fn name(self) -> &'static str {
        match self {
            Program::Tree => "tree",
            Program::Coral => "coral",
            Program::Tower => "tower",
            Program::Wall => "wall",
            Program::Terrain => "terrain",
            Program::Settlement => "settlement",
        }
    }

    /// Which space-filling rule a thing made this way follows.
    ///
    /// **One line per provenance, and not a plan.** D11's reduction is that
    /// geometry is not a species property: the *rule* is one of a handful, and
    /// its numbers are measured. This selects among the rules the way
    /// `solvers::for_tier` selects a solver, and everything that decides what
    /// the thing actually looks like is in `crate::recipe::generate`.
    ///
    /// Where an actor places a design, the habit comes from the design and this
    /// is not consulted — which is what `Recipe::Placed` is, and why a wooden
    /// box needed no variant here.
    pub fn habit(self) -> crate::recipe::Habit {
        use crate::recipe::Habit;
        match self {
            // Growing into a field it cannot see all of.
            Program::Tree | Program::Coral => Habit::Branching,
            // Laid down course on course.
            Program::Wall => Habit::Coursed,
            // Framed, storey on storey.
            Program::Tower => Habit::Framed,
            // A piece of a surface.
            Program::Terrain => Habit::Tiled,
            // A plan divided into plots.
            Program::Settlement => Habit::Subdivided,
        }
    }

    /// Emergent growth (target unknown, rule-driven) versus planned
    /// construction (target known, progress-driven). The two need different
    /// state and different advance laws, and conflating them is how you end up
    /// with buildings that grow organically towards the light.
    pub fn is_planned(self) -> bool {
        matches!(self, Program::Tower | Program::Wall | Program::Settlement)
    }

    /// The load a structure is proportioned against when it is created.
    ///
    /// A design load, not a survival load. A tree grows against the wind it
    /// meets most days and comes down in the storm it does not; an engineered
    /// frame is proportioned against a code gust with margin on top. Returning
    /// the *fluid* the structure lives in as well as the speed is what lets a
    /// coral be designed against a current at a thousand times the density of
    /// air without any of this having to know what a coral is.
    pub fn design_flow(self) -> (f64, f64) {
        match self {
            // Metres per second, and the fluid's density in kg/m^3.
            Program::Tree => (20.0, 1.225),
            Program::Coral => (1.2, 1025.0),
            Program::Tower => (42.0, 1.225),
            Program::Wall => (34.0, 1.225),
            // Ground is not proportioned against wind; it is proportioned
            // against what stands on it. The flow is what erodes it.
            Program::Terrain => (25.0, 1.225),
            // A settlement is designed to the same gust its buildings are.
            Program::Settlement => (42.0, 1.225),
        }
    }

    /// How heavy a cubic metre of the structure is, kg/m^3.
    ///
    /// **Read off the material, not off the program.** `docs/PLAY.md` D11: a
    /// wooden tower and a wooden tree have one density between them, and the
    /// place that knows is the material. This stays a method on `Program` only
    /// because a program is still what picks the material — which is D11's
    /// remaining work, not this — and it is now one line rather than a column.
    pub fn density(self) -> f64 {
        self.material().density
    }

    /// Free energy stored per kilogram of structure, J/kg.
    ///
    /// **Deliberately still a column**, and `docs/PLAY.md` D11 lists it as a
    /// material property. It is not one, and trying to retire it here is what
    /// showed why: run off `destruction_enthalpy`, terrain came out holding
    /// 8.2 MJ/kg of free energy, which is 820 J in a small hill that nothing
    /// put there.
    ///
    /// What this measures is not how tightly the material is bound but **how
    /// far uphill the making pushed it**, and that is a difference between the
    /// substance and the feedstock it was made from. A tree fixes cellulose out
    /// of carbon dioxide and water, which is a long way uphill; a hill's
    /// silicate came from silicate and is exactly where it started. The engine
    /// has cohesive energies and no formation enthalpies, so the comparison
    /// cannot be made yet.
    ///
    /// D11 puts the rest of its columns in **Ground**, where a derived erosion
    /// rate cannot coexist with a tabulated one; this belongs with them.
    pub fn energy_density(self) -> f64 {
        match self {
            Program::Tree | Program::Coral => BIOMASS_ENERGY,
            Program::Tower | Program::Wall | Program::Settlement => CONSTRUCTION_ENERGY,
            // Bedrock is already at the bottom of its own energy landscape.
            // There is no free energy stored in a hill.
            Program::Terrain => 0.0,
        }
    }

    /// What this program builds out of, as a material.
    ///
    /// D11's remaining species column, and the one this does not close: a
    /// program still names a material. What has changed is that the material it
    /// names is *derived* rather than tabulated — see `material.rs` — so the
    /// table is one of substances and histories rather than of thirteen numbers
    /// each. The column goes away when a generator seeds its node's mixture and
    /// the material is measured from that, which is Ground's work.
    pub fn material(self) -> crate::material::Material {
        match self {
            Program::Tree => crate::material::Material::green_wood(),
            Program::Coral => crate::material::Material::aragonite(),
            Program::Tower => crate::material::Material::reinforced_frame(),
            Program::Wall | Program::Settlement => crate::material::Material::masonry(),
            Program::Terrain => crate::material::Material::bedrock(),
        }
    }

    /// Maintenance cost per kilogram per second — respiration for the living
    /// programs, weathering and depreciation for the built ones.
    ///
    /// **Deliberately still a column, and this is the phase that measured
    /// why.** `docs/PLAY.md` D11 calls it a tabulated per-species decay rate
    /// that must be derived, and Phase 4 put both of the engine's candidate
    /// laws against it:
    ///
    /// * **Erosion** — §5.2's expression, now built in [`crate::erode`]. Green
    ///   wood's grain-scale cohesion is 5.16e7 Pa and air at 20 m/s presses
    ///   with 245 Pa, four orders short of the threshold, so the derived rate
    ///   for a tree standing in wind is **exactly zero**. A tree that pays no
    ///   maintenance grows without bound, and the emergent carrying capacity —
    ///   capture scaling with area against upkeep scaling with mass — goes with
    ///   it.
    /// * **Thermal degradation** — the other channel the engine has, an
    ///   attempt frequency from the lattice against the cohesive energy per
    ///   atom. Measured, per atom and at 291 K:
    ///
    /// ```text
    ///   cellulose   3.859 eV   E/kT 153.9   nu 3.48e13   1.6e-46 /yr
    ///   silica      5.678 eV   E/kT 226.4   nu 2.25e13   3.3e-78 /yr
    ///   aragonite   5.485 eV   E/kT 218.7   nu 2.29e13   7.5e-75 /yr
    ///   iron        1.835 eV   E/kT  73.2   nu 7.82e12   4.2e-12 /yr
    /// ```
    ///
    ///   against the 0.02/yr this column holds for a tree. **Forty-four orders
    ///   out**, and no choice of attempt frequency closes a gap that size.
    ///
    /// The reason is not a missing coefficient: a tree's maintenance is
    /// *metabolic*, the cost of running enzymatic turnover, and it is a
    /// property of being alive rather than of cellulose. The engine has no law
    /// for that, and inventing one here would be the `if is_forest` the axioms
    /// forbid wearing a rate constant's clothes. So the column stays, with its
    /// numbers measured rather than assumed, until the plan says what governs
    /// the upkeep of a living thing.
    pub fn maintenance(self) -> f64 {
        match self {
            Program::Tree => 0.02 / YEAR,
            Program::Coral => 0.05 / YEAR,
            Program::Tower | Program::Wall => 0.005 / YEAR,
            Program::Settlement => 0.008 / YEAR,
            // Erosion. Slow enough that a hill outlasts everything on it, and
            // not zero, because it is the same account as everything else.
            Program::Terrain => 1.0e-5 / YEAR,
        }
    }
}

/// Something that happened to this structure and is not in the nominal program.
///
/// The generalisation of `Node::epoch`, which can only say "the procedural
/// detail is stale" and not what changed. A branch breaking must survive
/// coarsening — the whole point of a structure is that its history is visible —
/// so the deviations are logged and replayed rather than discarded.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Event {
    pub at: f64,
    pub kind: EventKind,
    /// Which part of the structure. Interpreted by the program.
    pub site: u32,
    pub magnitude: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventKind {
    /// A limb removed — pruned, snapped, demolished.
    Severed,
    /// Local damage that heals or is repaired over time.
    Damaged,
    /// Growth suppressed here (shade, obstruction).
    Suppressed,
    /// Construction milestone reached.
    Completed,
}

/// The developmental state of one structure.
///
/// About 200 bytes against the 10^6 vertices it stands for. This is the whole
/// trick: the *state* is small even though the *structure* is not, because the
/// structure is recomputed from it on demand.
#[derive(Debug, Clone)]
pub struct Morphology {
    pub program: Program,
    /// Per-instance variation. Derived from the node's path key, so two trees
    /// in the same forest differ but the same tree is always itself.
    pub genome: [f32; 8],
    /// Developmental clock, seconds. Not the same as the node's coordinate
    /// time: a structure that spent a decade in shade is younger than its age.
    pub age: f64,
    /// Structural mass accumulated, kg.
    pub built: f64,
    /// Completion fraction for planned programs, 0..1. Unused for growth.
    pub progress: f64,
    /// Deviations from the nominal program, replayed at render time.
    pub events: Vec<Event>,
    /// Baked state, so the event log does not grow without bound.
    pub checkpoint_age: f64,
    /// Design mass for planned programs, kg. Growth programs leave it zero and
    /// discover their own ceiling from the allometry.
    pub design_mass: f64,
    /// **The generated program**: the rule the engine wrote down for placing
    /// this thing's bodies. See [`crate::recipe::Recipe`].
    ///
    /// `program` records what the thing is made *of* — provenance — and this
    /// records how it is put together. Nothing reads `program` for geometry:
    /// D11's whole argument is that a species table cannot be the answer to
    /// what shape a thing is, and D15 already proved it for a composite, whose
    /// recipe no variant describes.
    ///
    /// `None` only between a morphology being made and its recipe being
    /// generated, which is one call later. Everything that draws goes through
    /// it.
    pub recipe: Option<crate::recipe::Recipe>,
    /// **The field delta**: what has happened over this thing's surface that
    /// its own rule does not describe. `docs/PLAY.md` §5.
    ///
    /// The coarse half of §5.7's pair — `events` is the fine half. A field
    /// superposes, so a thousand footprints are a thousand of these summed
    /// rather than a heightmap, and each one decays at the rate its own
    /// material, flux and geometry set. See [`crate::erode`].
    pub field: Vec<crate::erode::Deviation>,
    /// **What this thing is actually made of**, by coarse element.
    ///
    /// `docs/PLAY.md` D11's third column, retired. A structure used to be built
    /// out of whatever its species was declared to be built out of; it is now
    /// built out of **what was there** — `Environment::feedstock`, measured off
    /// the node's own matter — and this is the mass-weighted blend of
    /// everything that has gone into it.
    ///
    /// Stored rather than re-derived, which is the third axiom: it is a fact
    /// about what happened, and what happened does not come back from the
    /// node's present composition once the structure is part of it.
    /// [`Composition::none`] until the first kilogram goes in.
    pub substrate: Composition,
}

impl Morphology {
    /// Seed a new structure. `genome` is derived from the address, so this is
    /// as reproducible as everything else.
    pub fn new(program: Program, world_seed: u64, path_key: u128, epoch: u32) -> Morphology {
        let mut st = Stream::at(world_seed, path_key, epoch, Purpose::Structure);
        let mut genome = [0.0f32; 8];
        for g in genome.iter_mut() {
            *g = st.uniform() as f32;
        }
        let mut m = Morphology {
            program,
            genome,
            age: 0.0,
            built: 0.0,
            progress: 0.0,
            events: Vec::new(),
            checkpoint_age: 0.0,
            design_mass: 0.0,
            recipe: None,
            field: Vec::new(),
            substrate: Composition::none(),
        };
        m.regenerate(&Environment::default(), &program.material());
        m
    }

    /// What this thing is made of, as far as anything knows.
    ///
    /// The blend of everything that has gone into it. A structure nobody has
    /// built anything into yet answers with the feedstock's default rather than
    /// with nothing, because a mass leaving it has to have a composition.
    pub fn made_of(&self, node: Composition) -> Composition {
        if self.substrate.is_none() {
            node
        } else {
            self.substrate
        }
    }

    /// Take `mass` of `feedstock` into the structure, and say what went in.
    ///
    /// Blended by mass, the same way `coarsen` blends anything else: a tree
    /// that spent its first decade in one soil and its second in another is
    /// made of both in the proportions it took them.
    fn absorb(&mut self, mass: f64, feedstock: Composition, standing: f64) -> Composition {
        if !(mass > 0.0) {
            return self.made_of(Composition::primordial());
        }
        let taken = if feedstock.is_none() { self.made_of(Composition::primordial()) } else { feedstock };
        self.substrate = if self.substrate.is_none() || !(standing > 0.0) {
            taken
        } else {
            Composition::blend(self.substrate, standing, taken, mass)
        };
        taken
    }

    /// What the recipe is a function of, right now.
    pub fn growth(&self) -> crate::recipe::Growth {
        crate::recipe::Growth::new(self.built, self.progress)
    }

    /// Write down the rule for placing this thing's bodies.
    ///
    /// **The analysis step.** `docs/PLAY.md` D11 and the owner's call on Phase
    /// 4: a program is *generated* by assessing what is actually there, not
    /// selected from a table of species. What is assessed is the mass, the
    /// material it is made of, and the fluid and light it is in — all of which
    /// the node already measures — and what comes out is one of a handful of
    /// space-filling rules with its numbers filled in.
    ///
    /// Re-run whenever the conditions change, because a recipe is a rule and
    /// not a frozen outcome: a wall somebody decides to build twice as long is
    /// a different rule, and a tree that has doubled its mass is the same one.
    pub fn regenerate(&mut self, env: &Environment, material: &crate::material::Material) {
        self.recipe = Some(crate::recipe::generate(
            self.program.habit(),
            &self.genome,
            self.design_mass.max(self.built),
            env,
            material,
        ));
    }

    /// A recipe generated from parts somebody put together.
    ///
    /// `program` is provenance, not species: it says what the parts are made
    /// of, so a box of oak planks burns and weighs like oak. The shape comes
    /// from the parts themselves, and no `Program` variant describes it — which
    /// is the whole of D11's argument applied to D15's composite.
    ///
    /// `built` is the parts' own mass rather than an accumulated one: an
    /// assembled thing did not grow into its size, it was made at it.
    pub fn assembled(
        program: Program,
        parts: crate::assembly::Assembly,
        world_seed: u64,
        path_key: u128,
    ) -> Morphology {
        let mut m = Morphology::new(program, world_seed, path_key, 0);
        m.built = parts.mass();
        m.design_mass = m.built;
        m.progress = 1.0;
        m.recipe = Some(crate::recipe::Recipe::Placed(parts));
        m
    }

    /// Is this thing's shape a list of parts somebody placed?
    pub fn is_assembled(&self) -> bool {
        self.assembly().is_some_and(|a| !a.is_empty())
    }

    /// The parts, where the recipe is a parts list.
    pub fn assembly(&self) -> Option<&crate::assembly::Assembly> {
        self.recipe.as_ref().and_then(|r| r.placed())
    }

    /// The parts, to be changed.
    pub fn assembly_mut(&mut self) -> Option<&mut crate::assembly::Assembly> {
        self.recipe.as_mut().and_then(|r| r.placed_mut())
    }

    /// Start, or replace, the parts list this thing is made of.
    pub fn set_assembly(&mut self, parts: crate::assembly::Assembly) {
        self.recipe = Some(crate::recipe::Recipe::Placed(parts));
    }

    /// Bulk density, kg/m³ — measured from the parts when there are parts, and
    /// read off the program when there are not.
    ///
    /// The first of D11's tabulated columns to become a measurement rather than
    /// a constant for anything that has the state to measure it. A box of oak
    /// planks with air inside is not oak's density and never was; its parts
    /// say what volume they occupy and its matter says what they weigh.
    pub fn density(&self) -> f64 {
        match self.recipe.as_ref().map(|r| r.density()) {
            Some(d) if d > 0.0 => d,
            // A recipe that has not been written yet, which is one call's
            // worth of window. See `Recipe::density` for why nothing else
            // reads the program for this.
            _ => self.program.density(),
        }
    }

    pub fn planned(program: Program, design_mass: f64, world_seed: u64, path_key: u128) -> Morphology {
        let mut m = Morphology::new(program, world_seed, path_key, 0);
        m.design_mass = design_mass;
        // The rule, against the conditions this constructor can see. A caller
        // with a world measures them and writes it again — see
        // `World::rewrite_recipe`.
        m.regenerate(&Environment::default(), &m.program.material());
        m
    }



    /// Characteristic size of the structure, metres.
    ///
    /// The node's radius is kept equal to this, so that geometry and matter
    /// agree by construction rather than by correction. **It is the bounding
    /// radius of what the recipe will actually draw**, which is one meaning
    /// rather than two: a structure used to state one size and be drawn at
    /// another, by up to 53%, because the sampler rescaled a unit skeleton
    /// until `summarise` reported the stated radius back. See
    /// `crate::recipe`'s module documentation for the measurement.
    pub fn extent(&self) -> f64 {
        match self.recipe.as_ref() {
            Some(r) => r.extent(self.growth()).max(1e-9),
            None => 1e-9,
        }
    }

    /// Overall height of the structure, whatever kind it is.
    pub fn height(&self) -> f64 {
        match self.recipe.as_ref() {
            Some(r) => r.height(self.growth()),
            None => 0.0,
        }
    }

    /// Light-intercepting area, m^2. Crown projection, not total leaf area —
    /// what limits a thing growing into a field is the ground it shades, not
    /// the foliage it carries. Zero for anything that is not growing into one.
    pub fn capture_area(&self) -> f64 {
        match self.recipe.as_ref() {
            Some(crate::recipe::Recipe::Branching(b)) => b.capture_area(self.growth()),
            _ => 0.0,
        }
    }

    /// Advance the developmental state, returning the transaction that has to
    /// balance for the step to be legitimate.
    ///
    /// Runs on the node's matter. The fine structure is never touched, never
    /// materialised, and does not need to exist.
    pub fn advance(&mut self, dt: f64, env: &Environment) -> GrowthStep {
        if dt <= 0.0 {
            return GrowthStep::none();
        }
        // **The habit decides what time does to this thing**, and the three
        // answers are different. Something growing into a field accumulates
        // what it can catch; something being built advances towards a design;
        // and something that was *made* or is simply *there* does neither.
        //
        // A box of oak planks has the oak's chemistry and none of a tree's
        // development: no light budget, no allometric ceiling, no progress to
        // advance. A patch of ground has less than that — nobody grew it and
        // nobody built it. Both still age, because weathering applies to them,
        // and that is the whole of what time does to them here.
        //
        // Running the growth law on ground was not merely pointless: terrain's
        // embodied energy is zero, so the gross-growth division was `0.0/0.0`,
        // and `f64::min` returns the *other* operand when one is NaN — so a
        // patch of ground quietly incorporated the whole of its own reservoir
        // every step.
        let txn = match self.recipe.as_ref() {
            Some(crate::recipe::Recipe::Branching(_)) => self.advance_growth(dt, env),
            Some(
                crate::recipe::Recipe::Coursed(_)
                | crate::recipe::Recipe::Framed(_)
                | crate::recipe::Recipe::Subdivided(_),
            ) => self.advance_construction(dt, env),
            _ => GrowthStep::none(),
        };
        self.age += dt * env.suppression();
        txn
    }

    fn advance_growth(&mut self, dt: f64, env: &Environment) -> GrowthStep {
        let program = self.program;
        // Seed mass: a structure has to start somewhere, and a zero-mass tree
        // has zero capture area and can never grow.
        if self.built <= 0.0 {
            self.built = SEED_MASS;
        }

        // The whole incident flux crosses the boundary. Only a sliver of it is
        // convertible; the rest is waste heat and re-radiation, and it has to
        // be on the books for the second law to have anything to check.
        let energy_absorbed = env.light_flux * self.capture_area() * dt;
        let usable = energy_absorbed
            * PHOTOSYNTHETIC_YIELD
            * env.thermal_factor()
            * env.water
            * (1.0 - env.crowding).max(0.0);

        // Gross new structure, before maintenance.
        let gross = usable / program.energy_density().max(1e-30);
        // Maintenance is paid out of standing structure: respiration for a
        // tree, dissolution for a coral. It scales with mass while capture
        // scales with area, so the two balance at a finite size and the
        // carrying capacity is *emergent* rather than an imposed constant.
        let upkeep = program.maintenance() * self.built * dt;
        let net = gross - upkeep;

        let limited = net.min(env.reservoir_mass.max(0.0));
        let before = self.built;
        self.built = (self.built + limited).max(0.0);
        let actual = self.built - before;

        // Energy locked into new structure, and energy freed by structure that
        // was respired away. Only one of the two is ever non-zero.
        let energy_stored = actual.max(0.0) * program.energy_density();
        let energy_released = (-actual).max(0.0) * program.energy_density();
        // What it is made of is what it was made *from*, blended in as it goes.
        let took = self.absorb(actual.max(0.0), env.feedstock, before);
        GrowthStep::build(
            actual.max(0.0),
            took,
            energy_absorbed,
            energy_stored,
            energy_released,
            THERMALISED_FRACTION,
            env.temperature,
        )
    }

    fn advance_construction(&mut self, dt: f64, env: &Environment) -> GrowthStep {
        let program = self.program;
        if self.design_mass <= 0.0 || self.progress >= 1.0 {
            return GrowthStep::none();
        }
        // Progress is limited by whichever of labour and materials runs out
        // first — the honest bottleneck on any real site.
        let by_labour = env.labour * dt;
        let remaining_mass = self.design_mass * (1.0 - self.progress);
        let by_material = if self.design_mass > 0.0 {
            env.reservoir_mass.max(0.0) / self.design_mass
        } else {
            0.0
        };
        let step = by_labour.min(by_material).min(1.0 - self.progress).max(0.0);

        let mass = self.design_mass * step;
        self.progress = (self.progress + step).min(1.0);
        self.built += mass;
        if self.progress >= 1.0 {
            self.events.push(Event {
                at: self.age,
                kind: EventKind::Completed,
                site: 0,
                magnitude: 1.0,
            });
        }
        let _ = remaining_mass;

        // Construction is inefficient: most of the energy poured in is lost as
        // process heat, and only the embodied energy ends up in the structure.
        let energy_stored = mass * program.energy_density();
        let energy_absorbed = energy_stored / CONSTRUCTION_EFFICIENCY;
        let took = self.absorb(mass, env.feedstock, self.built - mass);
        GrowthStep::build(
            mass,
            took,
            energy_absorbed,
            energy_stored,
            0.0,
            // A building site dumps most of its waste heat locally rather than
            // radiating it: kilns, curing concrete, machinery.
            0.6,
            env.temperature,
        )
    }

    /// Record a deviation from the nominal program.
    ///
    /// Returns the transaction it implies, so that severing a limb is booked
    /// like every other change: the structure loses mass, the free energy in
    /// that mass is released, and the matter itself stays in the node as
    /// litter. Nothing is created or destroyed — it just stops being part of
    /// the structure.
    pub fn record(&mut self, event: Event, temperature: f64, node: Composition) -> GrowthStep {
        let mut txn = GrowthStep::none();
        if event.kind == EventKind::Severed {
            let lost = (self.built * event.magnitude.clamp(0.0, 1.0)).max(0.0);
            self.built -= lost;
            txn = GrowthStep::build(
                -lost,
                self.made_of(node),
                0.0,
                0.0,
                lost * self.program.energy_density(),
                1.0,
                temperature,
            );
            txn.mass_incorporated = 0.0;
        }
        self.events.push(event);
        // Bound the log. Beyond this the oldest deviations are folded into the
        // genome as a permanent bias — a checkpoint, so replay stays finite for
        // a structure someone interacts with for a very long time.
        if self.events.len() > MAX_EVENTS {
            let drop = self.events.len() - MAX_EVENTS / 2;
            let bias: f64 = self.events[..drop].iter().map(|e| e.magnitude).sum::<f64>()
                / drop as f64;
            self.genome[7] = (self.genome[7] as f64 * (1.0 - 0.1 * bias)).clamp(0.0, 1.0) as f32;
            self.compact_events();
        }
        txn
    }

    /// Sever several sites at once, losing `fraction` of the structural mass
    /// between them.
    ///
    /// Separate from `record` because recording n breaks individually would
    /// compound the mass loss n times — each call taking a fraction of what the
    /// previous one left. A storm that breaks two hundred joints does not
    /// remove two hundred successive fractions of the tree.
    pub fn sever_many(&mut self, sites: &[u32], fraction: f64, node: Composition) -> GrowthStep {
        let lost = (self.built * fraction.clamp(0.0, 1.0)).max(0.0);
        self.built -= lost;
        for &site in sites {
            self.events.push(Event {
                at: self.age,
                kind: EventKind::Severed,
                site,
                magnitude: 0.0,
            });
        }
        self.compact_events();
        // The limb is on the ground, not gone: its free energy is still locked
        // in the wood. Nothing is released until something decomposes or burns
        // it, which is a separate process.
        GrowthStep {
            mass_incorporated: 0.0,
            composition: self.made_of(node),
            ..GrowthStep::none()
        }
        .with_detached(lost)
    }

    /// Consume structural mass outright — burned, vaporised — releasing the
    /// free energy that was holding it together.
    pub fn consume(&mut self, mass: f64, temperature: f64, node: Composition) -> GrowthStep {
        let lost = mass.clamp(0.0, self.built);
        self.built -= lost;
        GrowthStep::build(
            0.0,
            self.made_of(node),
            0.0,
            0.0,
            lost * self.program.energy_density(),
            0.5,
            temperature,
        )
    }

    fn compact_events(&mut self) {
        if self.events.len() > MAX_EVENTS {
            self.checkpoint_age = self.age;
            let drop = self.events.len() - MAX_EVENTS / 2;
            self.events.drain(..drop);
        }
    }



    // -- geometry ---------------------------------------------------------

    /// Place the structure's bodies, in metres, in the node's own frame.
    ///
    /// Pure in `(recipe, built, progress)`, and **`epoch` is not in that
    /// tuple**: taking a window out of a building bumps the node's epoch and
    /// everything sampled statistically in it redraws, and if a structure were
    /// epoch-seeded a player removing one window would watch the other nine
    /// jump. They do not, because a structure's variation comes from a recipe
    /// generated once from the path key and then stored.
    pub fn render(&self, budget: usize) -> Skeleton {
        let mut sk = match self.recipe.as_ref() {
            Some(r) => r.render(budget.max(1), self.growth(), &self.field),
            None => Skeleton::default(),
        };
        // A severed limb and everything above it is simply absent. Applied
        // here rather than inside a habit, because it is true of every habit
        // and was only ever implemented for one: `docs/PLAY.md` §5.9 measured
        // that four of the six programs never honoured a severance at all, and
        // the skip lived inside `render_branching`.
        if self.events.iter().any(|e| e.kind == EventKind::Severed) {
            sk = self.without_severed(sk);
        }
        sk
    }

    /// Drop every part the event log says is gone, and everything held on by
    /// one.
    ///
    /// `docs/PLAY.md` §5.9's second defect, which was measured and is closed
    /// here: "four of the six programs never honour a severance at all. Tower,
    /// wall, terrain and settlement return exactly the same part count from the
    /// first severance onward, with every severed site still present. The skip
    /// lives in `render_branching`, which only `Tree` and `Coral` call." A
    /// demolished wall section did not survive sixty-four edits; it did not
    /// survive one.
    ///
    /// Done over the emitted skeleton rather than inside each habit, so there
    /// is one implementation and a new habit cannot forget it.
    fn without_severed(&self, sk: Skeleton) -> Skeleton {
        let gone: std::collections::HashSet<u32> = self
            .events
            .iter()
            .filter(|e| e.kind == EventKind::Severed)
            .map(|e| e.site)
            .collect();
        let n = sk.len();
        let mut keep = vec![true; n];
        for i in 0..n {
            if gone.contains(&sk.site[i]) {
                keep[i] = false;
            }
        }
        // Anything held on by a part that is gone is gone with it. Emitted
        // order is support-before-supported for every habit here, so one
        // forward pass closes the transitive case.
        for i in 0..n {
            let sup = sk.support[i];
            if sup != NO_SUPPORT && (sup as usize) < n && !keep[sup as usize] {
                keep[i] = false;
            }
        }
        if keep.iter().all(|k| *k) {
            return sk;
        }
        let mut out = Skeleton::with_capacity(n);
        let mut remap = vec![NO_SUPPORT; n];
        for i in 0..n {
            if !keep[i] {
                continue;
            }
            remap[i] = out.len() as u32;
            let sup = sk.support[i];
            let support = if sup != NO_SUPPORT && (sup as usize) < n {
                remap[sup as usize]
            } else {
                NO_SUPPORT
            };
            out.push_joined(
                sk.base[i],
                sk.tip[i],
                sk.mass[i],
                sk.radius[i],
                support,
                sk.site[i],
                sk.joint_radius[i],
                sk.joint_bond[i],
            );
            *out.half.last_mut().unwrap() = sk.half[i];
            *out.free.last_mut().unwrap() = sk.free[i];
            *out.orientation.last_mut().unwrap() = sk.orientation[i];
        }
        for (a, b, f) in &sk.ties {
            let (a, b) = (*a as usize, *b as usize);
            if a < n && b < n && keep[a] && keep[b] {
                out.tie(remap[a], remap[b], *f);
            }
        }
        out
    }

    /// What the structure is physically made of, for the failure analysis.
    ///
    /// One dispatch site fewer than it looks: it is the program's, and the
    /// program's is one line per variant naming a *substance and a history*
    /// rather than thirteen numbers. See `Program::material`.
    pub fn material(&self) -> crate::material::Material {
        self.program.material()
    }

    /// The body kind the structure's parts should be tagged with.
    pub fn body_kind(&self) -> BodyKind {
        match self.program {
            Program::Tree | Program::Coral => BodyKind::Grain,
            Program::Tower | Program::Wall => BodyKind::Grain,
            Program::Terrain | Program::Settlement => BodyKind::Grain,
        }
    }

    /// Free energy currently locked in the structure, J.
    pub fn stored_energy(&self) -> f64 {
        self.built * self.program.energy_density()
    }

    /// Approximate byte cost of the developmental state — the number that has
    /// to be compared against the millions of vertices it replaces.
    pub fn state_bytes(&self) -> usize {
        std::mem::size_of::<Morphology>()
            + self.events.len() * std::mem::size_of::<Event>()
            + self.recipe.as_ref().map(|r| r.state_bytes()).unwrap_or(0)
    }
}

const SEED_MASS: f64 = 1e-4;
/// Fraction of waste energy that warms the structure rather than leaving as
/// radiation. A leaf runs only a few kelvin above ambient, so most of it goes.
const THERMALISED_FRACTION: f64 = 0.02;
const MAX_EVENTS: usize = 64;
/// Fraction of construction energy that ends up embodied rather than wasted.
const CONSTRUCTION_EFFICIENCY: f64 = 0.35;

/// Raw geometry in units of the structure's extent, together with the
/// connectivity that holds it together.
///
/// The connectivity was always implicit in the generators — a branching program
/// knows perfectly well which segment grew out of which — and was simply being
/// discarded. Keeping it costs three arrays and is what makes the difference
/// between a cloud of parts that happens to be tree-shaped and a structure that
/// can be loaded, stressed and broken.
#[derive(Debug, Clone, Default)]
pub struct Skeleton {
    /// Midpoint of each part.
    pub pos: Vec<Vec3>,
    pub mass: Vec<f64>,
    pub radius: Vec<f64>,
    /// Index of the part that supports this one; `NO_SUPPORT` for a part
    /// anchored to the ground.
    pub support: Vec<u32>,
    /// Program-stable name for this part, so an event can refer to it and mean
    /// the same thing after the structure is regenerated.
    pub site: Vec<u32>,
    /// Endpoints. The base is where the joint to the supporting part is, and
    /// therefore where the bending stress is highest and where things break.
    pub base: Vec<Vec3>,
    pub tip: Vec<Vec3>,
    /// The box each part presents, as half-extents along the part's own axes,
    /// or `Vec3::ZERO` for a part that is a tube of `radius`.
    ///
    /// **The generator emits the slab; nothing infers one.** `docs/PLAY.md`
    /// Phase 4: a grown or coursed structure emitted one capsule per member, so
    /// a masonry wall was a row of beads with a 0.169 m scallop between them
    /// and a patch of ground was a bed of cylinders with gaps at every corner
    /// of its grid. D18 rules out a grouping pass over the member list — a
    /// generator never infers a decomposition — so the generators that lay down
    /// flat things state them, which is what an assembly has always done.
    ///
    /// Zero is the discriminator, exactly as it is on [`crate::state::Body`]:
    /// `(r,r,r)` would be a cube whose bounding radius is `r*sqrt(3)`, so a
    /// sphere cannot be written as half-extents and the zero case has to carry
    /// it.
    pub half: Vec<Vec3>,
    /// Which axes of a boxed part the cross-section correction may move. See
    /// [`Skeleton::free_axes`].
    pub free: Vec<u8>,
    /// Which way each part is facing, in the node's own frame.
    ///
    /// Identity for everything laid out on an axis, which is every habit but
    /// one: a patch of ground on a *sphere* has cells whose own ups are not its
    /// up, because that is what being on a sphere means. A body carries an
    /// orientation for exactly this reason and it had nowhere to come from for
    /// a generated part.
    pub orientation: Vec<crate::math::Quat>,
    /// Radius of the *joint* at each part's base, metres, where that is not
    /// simply the member's own.
    ///
    /// Zero means "the member's", which is right for everything generated: a
    /// branch meets its parent across its own cross-section, and so does a
    /// course of masonry. An **assembled** thing is the case where it is not —
    /// `docs/PLAY.md` D15 is explicit that weld, glue and grown-together differ
    /// in what the join is made of and therefore in its strength, and the first
    /// half of that is how much of the part is actually stuck down. A 25 mm
    /// plank glued along one edge is held by a seam of its thickness, not by
    /// its whole face, and it comes off long before it snaps.
    pub joint_radius: Vec<f64>,
    /// What each part's *joint* is made of, where that is not what the part is
    /// made of. `UNSPECIATED` for everything generated: a branch meets its
    /// parent in wood. See [`crate::topology::Joint::bond`].
    pub joint_bond: Vec<crate::chem::SubstanceId>,
    /// Redundant connections `(a, b, fraction)` beyond the support forest —
    /// bracing, ties, anything giving load a second route to ground.
    ///
    /// The third value is the tie's cross-section as a *fraction* of the
    /// smaller member it joins, not an absolute area. Absolute areas are
    /// meaningless in a skeleton that gets scaled to whatever mass the
    /// structure has grown to, and a brace whose stiffness is out of proportion
    /// to its members makes the linear system impossible to condition — the
    /// solve then runs to its iteration cap and returns noise.
    pub ties: Vec<(u32, u32, f64)>,
}

/// A part anchored to the ground rather than to another part.
pub const NO_SUPPORT: u32 = u32::MAX;

/// Axis masks for [`Skeleton::free_axes`].
pub const FREE_X: u8 = 1;
pub const FREE_Y: u8 = 2;
pub const FREE_Z: u8 = 4;
pub const FREE_ALL: u8 = FREE_X | FREE_Y | FREE_Z;

impl Skeleton {
    pub fn with_capacity(n: usize) -> Skeleton {
        Skeleton {
            pos: Vec::with_capacity(n),
            mass: Vec::with_capacity(n),
            radius: Vec::with_capacity(n),
            support: Vec::with_capacity(n),
            site: Vec::with_capacity(n),
            base: Vec::with_capacity(n),
            tip: Vec::with_capacity(n),
            joint_radius: Vec::with_capacity(n),
            joint_bond: Vec::with_capacity(n),
            half: Vec::with_capacity(n),
            free: Vec::with_capacity(n),
            orientation: Vec::with_capacity(n),
            ties: Vec::new(),
        }
    }

    /// Add a redundant connection between two existing parts, sized as a
    /// fraction of the smaller member's cross-section.
    pub fn tie(&mut self, a: u32, b: u32, fraction: f64) {
        self.ties.push((a, b, fraction.clamp(0.0, 4.0)));
    }

    /// Add a part that is a segment between two points.
    pub fn push_segment(&mut self, base: Vec3, tip: Vec3, m: f64, r: f64, support: u32, site: u32) {
        self.pos.push((base + tip).scale(0.5));
        self.mass.push(m.max(1e-12));
        self.radius.push(r.max(1e-12));
        self.support.push(support);
        self.site.push(site);
        self.base.push(base);
        self.tip.push(tip);
        self.joint_radius.push(0.0);
        self.joint_bond.push(crate::chem::SubstanceId::UNSPECIATED);
        self.half.push(Vec3::ZERO);
        self.free.push(FREE_ALL);
        self.orientation.push(crate::math::Quat::IDENTITY);
    }

    /// Add a part that is a filled box rather than a tube.
    ///
    /// The member is still a member — `base`, `tip` and a cross-section radius
    /// are what the topology and the structural solve read, and they are
    /// unchanged — and `half` is the solid it *presents*. A block of masonry
    /// spans from one end to the other and is a box while it does it; the two
    /// descriptions are of the same thing and neither is inferred from the
    /// other.
    #[allow(clippy::too_many_arguments)]
    pub fn push_box(
        &mut self,
        base: Vec3,
        tip: Vec3,
        half: Vec3,
        free: u8,
        m: f64,
        r: f64,
        support: u32,
        site: u32,
    ) {
        self.push_segment(base, tip, m, r, support, site);
        *self.half.last_mut().unwrap() = v3(half.x.abs(), half.y.abs(), half.z.abs());
        *self.free.last_mut().unwrap() = free & FREE_ALL;
    }

    /// Which of a boxed part's own axes the cross-section correction may move.
    ///
    /// **Stated, never inferred.** The sampler scales member sections as a
    /// group so that the volume they enclose agrees with the structural mass at
    /// the material's density — a tube has one free dimension and the question
    /// never arises, and a box has three and the answer is different for every
    /// flat thing there is. A course of masonry may get thicker and may not get
    /// longer, or the courses stop meeting; a column of ground may get deeper
    /// and may not get narrower, or the patch stops being a surface. D18 says a
    /// generator never infers a decomposition, and this is the same sentence
    /// about the same information: the generator knows, so it says.
    ///
    /// Bit 0 is the part's local x, bit 1 y, bit 2 z. Zero means the box is
    /// fully stated and takes no correction at all.
    pub fn free_axes(&self, i: usize) -> u8 {
        self.free.get(i).copied().unwrap_or(FREE_ALL)
    }

    /// The axis a part runs along: its own direction where it has one, and its
    /// thinnest where it does not.
    ///
    /// What a cross-section is *across*. A beam's section is the two axes that
    /// are not its span; a plate's is its footprint, and the axis it is thin
    /// along is the one that plays the part of a span.
    pub fn axial_axis(&self, i: usize) -> usize {
        let d = self.tip[i] - self.base[i];
        if d.norm2() > 0.0 {
            let (dx, dy, dz) = (d.x.abs(), d.y.abs(), d.z.abs());
            return if dx >= dy && dx >= dz {
                0
            } else if dy >= dz {
                1
            } else {
                2
            };
        }
        let h = self.half.get(i).copied().unwrap_or(Vec3::ZERO);
        if h.x <= h.y && h.x <= h.z {
            0
        } else if h.y <= h.z {
            1
        } else {
            2
        }
    }

    /// Add a part held on by a joint narrower than the part itself — the
    /// assembled case. See [`Skeleton::joint_radius`].
    #[allow(clippy::too_many_arguments)]
    pub fn push_joined(
        &mut self,
        base: Vec3,
        tip: Vec3,
        m: f64,
        r: f64,
        support: u32,
        site: u32,
        joint_radius: f64,
        joint_bond: crate::chem::SubstanceId,
    ) {
        self.push_segment(base, tip, m, r, support, site);
        *self.joint_radius.last_mut().unwrap() = joint_radius.max(0.0);
        *self.joint_bond.last_mut().unwrap() = joint_bond;
    }

    /// Add a flat part that rests on several supports at once.
    ///
    /// **A plate is not a beam, and the difference is measurable.** A floor
    /// spanning its bay as a single member came out 77.6 m long with a 0.33 m
    /// section — slenderness 827 against the frame's next worst of 30 — and the
    /// static solve stopped converging. What a plate actually is, to a frame
    /// solver, is a *grillage*: strips of the plate running to each of the
    /// supports it rests on, each carrying the section a strip of it has. That
    /// is expressed here in the vocabulary the solver already has — one member
    /// with the plate's own section, tied to every other support — so any
    /// generator that lays down a flat thing gets it, and nothing in the solver
    /// has to know what a floor is.
    ///
    /// The member's section is the *equal-area* circle of the strip, not the
    /// plate's half-diagonal: handing the solver a bounding radius gave one
    /// element more stiffness than the rest of the structure put together.
    pub fn push_plate(
        &mut self,
        centre: Vec3,
        half: Vec3,
        m: f64,
        supports: &[u32],
        site: u32,
    ) {
        let h = v3(half.x.abs(), half.y.abs(), half.z.abs());
        // Span along the widest plan axis; the section is the other two.
        let (span, a, b) = if h.x >= h.y {
            (v3(h.x, 0.0, 0.0), h.y, h.z)
        } else {
            (v3(0.0, h.y, 0.0), h.x, h.z)
        };
        let radius = (4.0 * a * b / std::f64::consts::PI).sqrt().max(1e-12);
        let first = supports.first().copied().unwrap_or(NO_SUPPORT);
        let me = self.len() as u32;
        self.push_box(centre - span, centre + span, h, FREE_Z, m, radius, first, site);
        // Every other support carries its own strip. A tie is exactly what the
        // solver calls a second route to ground, and a plate on four columns is
        // four routes.
        for &sup in supports.iter().skip(1) {
            if sup != NO_SUPPORT {
                self.tie(me, sup, 1.0);
            }
        }
    }

    /// Add a part with no extent of its own — a block, a slab, a parcel.
    pub fn push(&mut self, p: Vec3, m: f64, r: f64) {
        let half = v3(0.0, 0.0, r);
        let site = self.pos.len() as u32;
        self.push_segment(p - half, p + half, m, r, NO_SUPPORT, site);
    }

    /// Add a part supported by another.
    pub fn push_supported(&mut self, p: Vec3, m: f64, r: f64, support: u32) {
        let half = v3(0.0, 0.0, r);
        let site = self.pos.len() as u32;
        self.push_segment(p - half, p + half, m, r, support, site);
    }

    pub fn len(&self) -> usize {
        self.pos.len()
    }
    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }
    /// Length of each part along its own axis.
    pub fn length(&self, i: usize) -> f64 {
        (self.tip[i] - self.base[i]).norm()
    }
    /// Unit direction of each part.
    pub fn direction(&self, i: usize) -> Vec3 {
        (self.tip[i] - self.base[i]).unit()
    }
}

/// Conditions the structure grows in. Derived from the node's own matter
/// plus whatever is arriving from outside.
#[derive(Debug, Clone, Copy)]
pub struct Environment {
    /// Incident radiation, W/m^2.
    pub light_flux: f64,
    pub temperature: f64,
    /// Water availability, 0..1.
    pub water: f64,
    /// Competition for the same resource, 0..1.
    pub crowding: f64,
    /// Feedstock available in the node this step, kg.
    pub reservoir_mass: f64,
    /// Construction rate, fraction of the design per second.
    pub labour: f64,
    /// Density of the fluid this thing is standing in, kg/m^3.
    ///
    /// **Measured, not named.** D11: "A coral is in water because its node's
    /// mixture is water; nobody tells it." This is the mass of the node's own
    /// non-solid fraction over the volume it occupies, so a thing growing in
    /// water is proportioned against water and one in air against air, and the
    /// two species the table used to hold become one habit with a measurement
    /// in it. Air at sea level where nothing has been described.
    pub fluid_density: f64,
    /// The flow this thing has actually met, m/s.
    ///
    /// D11 again: "The design gust is the gust the structure has met, which its
    /// own history already holds." A structure is proportioned against what it
    /// lives through, so a tree on a headland is stouter than one in a valley
    /// without either of them being told which it is.
    pub flow_speed: f64,
    /// **What is here to build out of**, by coarse element.
    ///
    /// `docs/PLAY.md` D11's `substrate` column, made a measurement: a structure
    /// is built out of what its node holds, so a tree in a silicate world is
    /// made of silicate and nobody tells it which. The node's own
    /// `matter.composition`, and [`Composition::none`] where nothing has been
    /// said.
    pub feedstock: Composition,
}

impl Default for Environment {
    fn default() -> Self {
        Environment {
            light_flux: 200.0,
            temperature: 288.0,
            water: 1.0,
            crowding: 0.0,
            reservoir_mass: f64::INFINITY,
            labour: 0.0,
            fluid_density: 1.225,
            flow_speed: 20.0,
            feedstock: Composition::none(),
        }
    }
}

impl Environment {
    /// Temperature response: a broad optimum around 298 K, falling to zero at
    /// freezing and at protein denaturation. Growth stops in winter, which is
    /// what makes tree rings.
    pub fn thermal_factor(&self) -> f64 {
        let t = self.temperature;
        if !(273.0..=323.0).contains(&t) {
            return 0.0;
        }
        let x = (t - 298.0) / 20.0;
        (1.0 - x * x).max(0.0)
    }

    /// How much of elapsed time counts as developmental time. A structure in
    /// the dark does not age towards maturity.
    pub fn suppression(&self) -> f64 {
        (self.thermal_factor() * self.water * (1.0 - self.crowding)).clamp(0.0, 1.0)
    }
}

/// One growth or construction step, as a set of books that must balance.
///
/// The point of routing every change through this type is that the second law
/// becomes a precondition rather than an aspiration. `validate` is called
/// before the transaction is applied, so a program cannot silently mint free
/// energy or order.
#[derive(Debug, Clone, Copy, Default)]
pub struct GrowthStep {
    /// Mass moved from the surrounding reservoir into the structure, kg.
    pub mass_incorporated: f64,
    /// Mass that left the structure but stayed in the node, kg. A fallen limb
    /// is litter, not an absence.
    pub mass_detached: f64,
    /// What that mass is made of.
    pub composition: Composition,
    /// Energy crossing the node boundary inwards, J.
    pub energy_absorbed: f64,
    /// Energy now held as chemical or structural free energy, J.
    pub energy_stored: f64,
    /// Energy thermalised locally, J. Stays in the node as internal energy.
    pub heat_released: f64,
    /// Energy re-radiated back out across the node boundary, J.
    ///
    /// A leaf absorbs the whole solar flux and stores about 0.3% of it. The
    /// other 99.7% leaves again, mostly as thermal infrared. Booking only the
    /// fraction that gets used — the obvious simplification — describes a
    /// perfectly efficient converter, and `validate` rightly refuses it,
    /// because a device that turns all of its input into stored free energy
    /// while lowering its own entropy is a second-law violation.
    pub energy_radiated: f64,
    /// Free energy liberated by structure that was lost this step, J.
    /// Respiration, dissolution, demolition.
    pub energy_released: f64,
    /// Entropy change of the structure itself, J/K. Negative when it orders.
    pub entropy_local: f64,
    /// Entropy delivered to the surroundings, J/K. Never negative.
    pub entropy_exported: f64,
}

impl GrowthStep {
    pub fn none() -> GrowthStep {
        GrowthStep {
            composition: Composition::primordial(),
            ..Default::default()
        }
    }

    fn with_detached(mut self, mass: f64) -> GrowthStep {
        self.mass_detached = mass;
        self
    }

    #[allow(clippy::too_many_arguments)]
    fn build(
        mass: f64,
        composition: Composition,
        energy_absorbed: f64,
        energy_stored: f64,
        energy_released: f64,
        thermalised_fraction: f64,
        temperature: f64,
    ) -> GrowthStep {
        let t = temperature.max(2.725);
        // Whatever is not stored is waste. A little of it warms the structure;
        // the rest leaves as radiation.
        let waste = (energy_absorbed + energy_released - energy_stored).max(0.0);
        let heat_released = waste * thermalised_fraction;
        let energy_radiated = waste - heat_released;
        GrowthStep {
            mass_incorporated: mass,
            mass_detached: 0.0,
            composition,
            energy_absorbed,
            energy_stored,
            heat_released,
            energy_radiated,
            energy_released,
            // Net ordering: building lowers local entropy, decomposing raises it.
            entropy_local: ORDERING_FRACTION * (energy_released - energy_stored) / t,
            entropy_exported: (heat_released + energy_radiated) / t,
        }
    }

    /// Does this step obey the first and second laws?
    ///
    /// Returns the reason it does not, so a failure is diagnosable rather than
    /// merely a rejection.
    pub fn validate(&self) -> Result<(), &'static str> {
        let scale = self
            .energy_absorbed
            .abs()
            .max(self.energy_stored.abs())
            .max(1e-30);
        // First law across the structure: what came in, plus what was freed by
        // decomposition, equals what was stored, warmed and radiated away.
        let inflow = self.energy_absorbed + self.energy_released;
        let outflow = self.energy_stored + self.heat_released + self.energy_radiated;
        if (inflow - outflow).abs() > 1e-9 * scale {
            return Err("energy in does not equal energy stored, warmed and radiated");
        }
        if self.energy_radiated < -1e-9 * scale {
            return Err("negative radiation");
        }
        if self.heat_released < -1e-9 * scale {
            return Err("negative heat release: the step is refrigerating for free");
        }
        if self.entropy_exported < 0.0 {
            return Err("negative entropy export");
        }
        if self.entropy_local + self.entropy_exported < -1e-12 * self.entropy_exported.abs().max(1e-30) {
            return Err("total entropy decreased: second law violated");
        }
        if self.mass_incorporated < 0.0 {
            return Err("negative mass incorporated");
        }
        Ok(())
    }

    pub fn total_entropy_change(&self) -> f64 {
        self.entropy_local + self.entropy_exported
    }

    /// Net energy the node gains: what crossed the boundary inwards minus what
    /// left. This is the quantity the engine's books have to match.
    pub fn net_boundary_flux(&self) -> f64 {
        self.energy_absorbed - self.energy_radiated
    }
}
