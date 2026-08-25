//! Morphology: matter that has a history rather than a temperature.
//!
//! # Why the rest of the engine cannot represent a tree
//!
//! Everything `prolong` regenerates is *ergodic*. One max-entropy sample of a
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
//! Growth runs on the *aggregate*, never on the fine structure. A forest does
//! not grow by integrating 10^9 trees; it grows by advancing one ordinary
//! differential equation on a forest node, at O(1) per node. That is cheap
//! enough to run on the entire world every frame, coarse or not — so the
//! laziness that the rest of the engine works hard for is simply free here.
//!
//! # Growth is a transaction, not an exemption
//!
//! Building order out of disorder costs free energy and exports entropy. Every
//! step returns a [`Transaction`] that has to balance before it is applied:
//! energy in equals energy stored plus heat released, and local entropy plus
//! exported entropy is non-negative. A program that tried to grow too
//! efficiently would be rejected by `Transaction::validate` rather than
//! quietly minting free energy.

use crate::math::{v3, Vec3};
use crate::rng::{Purpose, Stream};
use crate::state::{BodyKind, Composition};
use crate::units::*;

/// Free energy density of dry biomass, J/kg. Wood is ~17-19 MJ/kg.
pub const BIOMASS_ENERGY: f64 = 17.0e6;
/// Bulk density of wood, kg/m^3.
pub const WOOD_DENSITY: f64 = 600.0;
/// Bulk density of coral skeleton (aragonite), kg/m^3.
pub const CORAL_DENSITY: f64 = 2700.0;
/// Embodied energy of reinforced concrete and steel construction, J/kg.
pub const CONSTRUCTION_ENERGY: f64 = 2.5e6;
/// Bulk density of a framed building including voids, kg/m^3.
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
}

impl Program {
    pub fn name(self) -> &'static str {
        match self {
            Program::Tree => "tree",
            Program::Coral => "coral",
            Program::Tower => "tower",
            Program::Wall => "wall",
        }
    }

    /// Emergent growth (target unknown, rule-driven) versus planned
    /// construction (target known, progress-driven). The two need different
    /// state and different advance laws, and conflating them is how you end up
    /// with buildings that grow organically towards the light.
    pub fn is_planned(self) -> bool {
        matches!(self, Program::Tower | Program::Wall)
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
        }
    }

    pub fn density(self) -> f64 {
        match self {
            Program::Tree => WOOD_DENSITY,
            Program::Coral => CORAL_DENSITY,
            Program::Tower | Program::Wall => BUILDING_DENSITY,
        }
    }

    /// Free energy stored per kilogram of structure.
    pub fn energy_density(self) -> f64 {
        match self {
            Program::Tree | Program::Coral => BIOMASS_ENERGY,
            Program::Tower | Program::Wall => CONSTRUCTION_ENERGY,
        }
    }

    /// What the structure is built out of.
    pub fn substrate(self) -> Composition {
        let mut c = [0.0; NSPECIES];
        match self {
            // Cellulose, CH2O: roughly 44% C, 6% H, 50% O by mass.
            Program::Tree => {
                c[Species::Carbon as usize] = 0.44;
                c[Species::Hydrogen as usize] = 0.06;
                c[Species::Oxygen as usize] = 0.50;
            }
            // Calcium carbonate, lumped: carbon, oxygen, and heavier cations.
            Program::Coral => {
                c[Species::Carbon as usize] = 0.12;
                c[Species::Oxygen as usize] = 0.48;
                c[Species::Other as usize] = 0.40;
            }
            // Concrete and steel: silicates, oxygen, iron.
            Program::Tower | Program::Wall => {
                c[Species::Oxygen as usize] = 0.46;
                c[Species::Silicon as usize] = 0.27;
                c[Species::Iron as usize] = 0.12;
                c[Species::Carbon as usize] = 0.03;
                c[Species::Other as usize] = 0.12;
            }
        }
        Composition(c).normalised()
    }

    /// Maintenance cost per kilogram per second — respiration for the living
    /// programs, weathering and depreciation for the built ones.
    pub fn maintenance(self) -> f64 {
        match self {
            Program::Tree => 0.02 / YEAR,
            Program::Coral => 0.05 / YEAR,
            Program::Tower | Program::Wall => 0.005 / YEAR,
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
        Morphology {
            program,
            genome,
            age: 0.0,
            built: 0.0,
            progress: 0.0,
            events: Vec::new(),
            checkpoint_age: 0.0,
            design_mass: 0.0,
        }
    }

    pub fn planned(program: Program, design_mass: f64, world_seed: u64, path_key: u128) -> Morphology {
        let mut m = Morphology::new(program, world_seed, path_key, 0);
        m.design_mass = design_mass;
        m
    }

    /// Genome value in `[lo, hi]`.
    fn gene(&self, i: usize, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.genome[i % 8] as f64
    }

    /// Characteristic size of the structure, metres.
    ///
    /// The aggregate's radius is kept equal to this, so that geometry and bulk
    /// state agree by construction rather than by correction. The developmental
    /// state is the authority on how big the thing is; the aggregate follows.
    pub fn extent(&self) -> f64 {
        match self.program {
            Program::Tree => self.tree_height().max(1e-3) * 0.5,
            Program::Coral => {
                let v = self.built / CORAL_DENSITY;
                (v / 0.3).cbrt().max(1e-3)
            }
            Program::Tower => {
                let (floors, side) = self.tower_design();
                let h = floors as f64 * FLOOR_HEIGHT * self.progress.max(1e-3);
                0.5 * (h * h + 2.0 * side * side).sqrt()
            }
            Program::Wall => {
                let (len, h, _) = self.wall_design();
                0.5 * ((len * self.progress.max(1e-3)).powi(2) + h * h).sqrt()
            }
        }
    }

    /// Overall height of the structure, whatever kind it is.
    ///
    /// `tree_height` applies quarter-power allometry, which is meaningful for
    /// something that grew under its own self-loading and meaningless for a
    /// wall.
    pub fn height(&self) -> f64 {
        match self.program {
            Program::Tree => self.tree_height(),
            _ => self.extent() * 2.0,
        }
    }

    /// Tree height from biomass, by elastic self-similarity.
    ///
    /// McMahon's buckling criterion gives trunk radius ∝ height^1.5, so volume
    /// ∝ height^4 and height ∝ volume^(1/4) — the same quarter-power that runs
    /// through all of allometry. Calibrated so a one-tonne tree stands ~15 m.
    pub fn tree_height(&self) -> f64 {
        let v = self.built / WOOD_DENSITY;
        let shape = self.gene(0, 0.8, 1.25);
        (v / 3.3e-5).max(0.0).powf(0.25) * shape
    }

    /// Light-intercepting area, m^2. Crown projection, not total leaf area —
    /// what limits a tree is the ground it shades, not the foliage it carries.
    pub fn capture_area(&self) -> f64 {
        match self.program {
            Program::Tree => {
                let h = self.tree_height();
                let crown = 0.3 * h * self.gene(1, 0.75, 1.3);
                std::f64::consts::PI * crown * crown
            }
            Program::Coral => {
                let r = self.extent();
                std::f64::consts::PI * r * r * 1.5
            }
            _ => 0.0,
        }
    }

    /// Advance the developmental state, returning the transaction that has to
    /// balance for the step to be legitimate.
    ///
    /// Runs on the aggregate. The fine structure is never touched, never
    /// materialised, and does not need to exist.
    pub fn advance(&mut self, dt: f64, env: &Environment) -> Transaction {
        if dt <= 0.0 {
            return Transaction::none();
        }
        let txn = if self.program.is_planned() {
            self.advance_construction(dt, env)
        } else {
            self.advance_growth(dt, env)
        };
        self.age += dt * env.suppression();
        txn
    }

    fn advance_growth(&mut self, dt: f64, env: &Environment) -> Transaction {
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
        let gross = usable / program.energy_density();
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
        Transaction::build(
            actual.max(0.0),
            program.substrate(),
            energy_absorbed,
            energy_stored,
            energy_released,
            THERMALISED_FRACTION,
            env.temperature,
        )
    }

    fn advance_construction(&mut self, dt: f64, env: &Environment) -> Transaction {
        let program = self.program;
        if self.design_mass <= 0.0 || self.progress >= 1.0 {
            return Transaction::none();
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
        Transaction::build(
            mass,
            program.substrate(),
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
    pub fn record(&mut self, event: Event, temperature: f64) -> Transaction {
        let mut txn = Transaction::none();
        if event.kind == EventKind::Severed {
            let lost = (self.built * event.magnitude.clamp(0.0, 1.0)).max(0.0);
            self.built -= lost;
            txn = Transaction::build(
                -lost,
                self.program.substrate(),
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
    pub fn sever_many(&mut self, sites: &[u32], fraction: f64) -> Transaction {
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
        Transaction {
            mass_incorporated: 0.0,
            composition: self.program.substrate(),
            ..Transaction::none()
        }
        .with_detached(lost)
    }

    /// Consume structural mass outright — burned, vaporised — releasing the
    /// free energy that was holding it together.
    pub fn consume(&mut self, mass: f64, temperature: f64) -> Transaction {
        let lost = mass.clamp(0.0, self.built);
        self.built -= lost;
        Transaction::build(
            0.0,
            self.program.substrate(),
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

    /// Total mass severed from the structure by logged events.
    fn severed_fraction(&self) -> f64 {
        let lost: f64 = self
            .events
            .iter()
            .filter(|e| e.kind == EventKind::Severed)
            .map(|e| e.magnitude)
            .sum();
        lost.clamp(0.0, 0.95)
    }

    // -- geometry ---------------------------------------------------------

    /// Generate the structure's geometry, in units of `extent()`.
    ///
    /// Pure in `(genome, age, built, progress, events)`. Returns positions,
    /// relative masses and per-part radii; `prolong` scales them so the totals
    /// match the aggregate exactly, exactly as it does for a sampled cloud.
    pub fn render(&self, budget: usize) -> Skeleton {
        let n = budget.max(1);
        match self.program {
            Program::Tree => self.render_branching(n, 0.62, 3),
            Program::Coral => self.render_branching(n, 0.72, 4),
            Program::Tower => self.render_tower(n),
            Program::Wall => self.render_wall(n),
        }
    }

    /// Recursive branching, shared by tree and coral.
    ///
    /// Segments are laid down depth-first with a mass budget: each level takes
    /// `taper^3` of its parent's cross-section, so the structure obeys
    /// da Vinci's rule (total cross-section is preserved across a branch point)
    /// and therefore looks like a tree rather than like a fractal.
    fn render_branching(&self, budget: usize, taper: f64, splits: usize) -> Skeleton {
        let mut sk = Skeleton::with_capacity(budget);
        let severed = self.severed_fraction();
        let lean = self.gene(2, -0.12, 0.12);
        let twist = self.gene(3, 0.0, std::f64::consts::TAU);
        let spread = self.gene(4, 0.45, 0.85);

        // Depth grows with the structure: a seedling is a stick, a mature tree
        // has six or seven orders of branching. This is what makes the same
        // program produce a plausible object at every age.
        let max_depth = (2.0 + (self.built.max(SEED_MASS) / SEED_MASS).log10() * 1.4)
            .clamp(1.0, 8.0) as usize;

        struct Seg {
            base: Vec3,
            dir: Vec3,
            len: f64,
            rad: f64,
            depth: usize,
            id: u32,
            /// Segment id of the part this one grew out of.
            parent: u32,
        }
        // Breadth-first, not depth-first. With a stack the budget is spent
        // rendering one branch down to its finest twigs while the rest of the
        // tree is simply absent — which looks broken at any budget below
        // saturation, and makes the level-of-detail meaningless. A queue spends
        // the budget on the structurally significant parts first, so a coarse
        // render is a coarse *tree* rather than a fragment of one.
        let mut queue = std::collections::VecDeque::from(vec![Seg {
            base: v3(0.0, 0.0, -1.0),
            dir: v3(lean, lean * 0.5, 1.0).unit(),
            len: 0.55,
            rad: 0.055,
            depth: 0,
            id: 0,
            parent: NO_SUPPORT,
        }]);
        let mut next_id = 1u32;
        // Segment id to emitted index. Breadth-first guarantees a parent is
        // emitted before its children, so every lookup here succeeds.
        let mut emitted: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();

        while let Some(s) = queue.pop_front() {
            if sk.len() >= budget {
                break;
            }
            // A severed limb and everything above it is simply absent.
            if self
                .events
                .iter()
                .any(|e| e.kind == EventKind::Severed && e.site == s.id)
            {
                continue;
            }
            let tip = s.base + s.dir.scale(s.len);
            let support = if s.parent == NO_SUPPORT {
                NO_SUPPORT
            } else {
                match emitted.get(&s.parent) {
                    Some(&idx) => idx,
                    // The parent was pruned or fell outside the budget, so this
                    // segment has nothing to hang from and is not emitted.
                    None => continue,
                }
            };
            emitted.insert(s.id, sk.len() as u32);
            sk.push_segment(s.base, tip, s.rad * s.rad * s.len, s.rad, support, s.id);
            if s.depth >= max_depth {
                continue;
            }
            let child_rad = s.rad * taper;
            let child_len = s.len * (0.62 + 0.12 * self.gene(5, 0.0, 1.0));
            for k in 0..splits {
                let phi = twist + k as f64 * std::f64::consts::TAU / splits as f64
                    + s.depth as f64 * 0.7;
                let axis = v3(phi.cos(), phi.sin(), 0.0);
                let dir = (s.dir + axis.scale(spread)).unit();
                queue.push_back(Seg {
                    base: tip,
                    dir,
                    len: child_len,
                    rad: child_rad,
                    depth: s.depth + 1,
                    id: next_id,
                    parent: s.id,
                });
                next_id += 1;
            }
        }
        let _ = severed;
        sk
    }

    fn tower_design(&self) -> (usize, f64) {
        let floors = (6.0 + self.gene(0, 0.0, 34.0)).round() as usize;
        let side = 0.6 + self.gene(1, 0.0, 0.5);
        (floors.max(1), side)
    }

    /// A framed tower, built floor by floor.
    ///
    /// Planned construction is the easier case: the target is known in advance,
    /// so materialising a half-built tower is the finished design masked by the
    /// completion fraction, with the topmost floor partial.
    ///
    /// The members are real segments — columns spanning a storey, beams
    /// spanning between columns. An earlier version emitted each element as a
    /// point with a nominal radius, which was adequate while the parts were
    /// only drawn and became nonsense the moment they had to carry load: a
    /// zero-length member has no bending stiffness, and the density correction
    /// inflated its radius until the tower rendered as a smear of vertical
    /// streaks.
    fn render_tower(&self, budget: usize) -> Skeleton {
        let (floors, side) = self.tower_design();
        let mut sk = Skeleton::with_capacity(budget);
        let built_floors = self.progress * floors as f64;
        let storey = 2.0 / floors as f64;
        const COLS: usize = 4;
        let corner = |c: usize, z: f64| {
            let a = std::f64::consts::TAU * c as f64 / COLS as f64 + std::f64::consts::FRAC_PI_4;
            v3(side * a.cos(), side * a.sin(), z)
        };

        let mut below = [NO_SUPPORT; COLS];
        let mut site = 0u32;
        for f in 0..floors {
            let complete = (built_floors - f as f64).clamp(0.0, 1.0);
            if complete <= 0.0 || sk.len() + 2 * COLS > budget {
                break;
            }
            let z0 = -1.0 + storey * f as f64;
            let z1 = z0 + storey * complete;

            // Columns: one storey each, standing on the column below.
            let mut here = [NO_SUPPORT; COLS];
            for c in 0..COLS {
                here[c] = sk.len() as u32;
                sk.push_segment(corner(c, z0), corner(c, z1), 1.0, 0.030, below[c], site);
                site += 1;
            }
            // Beams: the floor plate, spanning between column heads.
            if complete >= 0.999 {
                for c in 0..COLS {
                    let n = (c + 1) % COLS;
                    sk.push_segment(corner(c, z1), corner(n, z1), 0.7, 0.020, here[c], site);
                    site += 1;
                }
                // Cross-bracing between adjacent columns, and diagonally to the
                // storey below. This is what makes a frame a frame rather than
                // a stack of posts — and it makes the structure statically
                // indeterminate, so the redundant solver has a generated case
                // to work on and not only a test rig.
                for c in 0..COLS {
                    let n = (c + 1) % COLS;
                    sk.tie(here[c], here[n], 0.33);
                    if below[c] != NO_SUPPORT {
                        sk.tie(here[c], below[n], 0.25);
                    }
                }
            }
            below = here;
        }
        if sk.is_empty() {
            sk.push_segment(v3(0.0, 0.0, -1.0), v3(0.0, 0.0, -0.9), 1.0, 0.03, NO_SUPPORT, 0);
        }
        sk
    }

    fn wall_design(&self) -> (f64, f64, usize) {
        let len = 1.6;
        let h = 0.3 + self.gene(0, 0.0, 0.5);
        let courses = (4.0 + self.gene(1, 0.0, 12.0)).round().max(1.0) as usize;
        (len, h, courses)
    }

    /// A wall, laid course by course.
    ///
    /// Each block is a horizontal member resting on the course beneath, offset
    /// by half a block on alternate courses — which is what stops a wall being
    /// a set of independent vertical columns of brick, and is why the load path
    /// is a tree that runs diagonally down to the footing.
    fn render_wall(&self, budget: usize) -> Skeleton {
        let (len, h, courses) = self.wall_design();
        let mut sk = Skeleton::with_capacity(budget);
        let per_course = (budget / courses.max(1)).clamp(2, 64);
        let laid = self.progress * courses as f64;
        let block = 2.0 * len / per_course as f64;
        let mut prev_start = 0u32;
        let mut prev_count = 0usize;
        let mut site = 0u32;
        for c in 0..courses {
            let complete = (laid - c as f64).clamp(0.0, 1.0);
            if complete <= 0.0 || sk.len() >= budget {
                break;
            }
            let z = -h + 2.0 * h * (c as f64 + 0.5) / courses as f64;
            let blocks = ((per_course as f64 * complete).round() as usize).max(1);
            let offset = if c % 2 == 0 { 0.0 } else { 0.5 };
            let start = sk.len() as u32;
            for b in 0..blocks {
                let x0 = -len + block * (b as f64 + offset);
                let support = if c == 0 || prev_count == 0 {
                    NO_SUPPORT
                } else {
                    prev_start + (b.min(prev_count - 1)) as u32
                };
                sk.push_segment(
                    v3(x0, 0.0, z),
                    v3(x0 + block, 0.0, z),
                    1.0,
                    block * 0.35,
                    support,
                    site,
                );
                site += 1;
                if sk.len() >= budget {
                    break;
                }
            }
            prev_start = start;
            prev_count = sk.len() as usize - start as usize;
        }
        if sk.is_empty() {
            sk.push_segment(v3(-0.1, 0.0, -1.0), v3(0.1, 0.0, -1.0), 1.0, 0.05, NO_SUPPORT, 0);
        }
        sk
    }

    /// What the structure is physically made of, for the failure analysis.
    pub fn material(&self) -> crate::topology::Material {
        match self.program {
            Program::Tree => crate::topology::Material::GREEN_WOOD,
            Program::Coral => crate::topology::Material::ARAGONITE,
            Program::Tower => crate::topology::Material::REINFORCED_FRAME,
            Program::Wall => crate::topology::Material::MASONRY,
        }
    }

    /// The body kind the structure's parts should be tagged with.
    pub fn body_kind(&self) -> BodyKind {
        match self.program {
            Program::Tree | Program::Coral => BodyKind::Grain,
            Program::Tower | Program::Wall => BodyKind::Grain,
        }
    }

    /// Free energy currently locked in the structure, J.
    pub fn stored_energy(&self) -> f64 {
        self.built * self.program.energy_density()
    }

    /// Approximate byte cost of the developmental state — the number that has
    /// to be compared against the millions of vertices it replaces.
    pub fn state_bytes(&self) -> usize {
        std::mem::size_of::<Morphology>() + self.events.len() * std::mem::size_of::<Event>()
    }
}

const SEED_MASS: f64 = 1e-4;
/// Fraction of waste energy that warms the structure rather than leaving as
/// radiation. A leaf runs only a few kelvin above ambient, so most of it goes.
const THERMALISED_FRACTION: f64 = 0.02;
const MAX_EVENTS: usize = 64;
const FLOOR_HEIGHT: f64 = 3.2;
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

/// Conditions the structure grows in. Derived from the node's own aggregate
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
pub struct Transaction {
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

impl Transaction {
    pub fn none() -> Transaction {
        Transaction {
            composition: Composition::primordial(),
            ..Default::default()
        }
    }

    fn with_detached(mut self, mass: f64) -> Transaction {
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
    ) -> Transaction {
        let t = temperature.max(2.725);
        // Whatever is not stored is waste. A little of it warms the structure;
        // the rest leaves as radiation.
        let waste = (energy_absorbed + energy_released - energy_stored).max(0.0);
        let heat_released = waste * thermalised_fraction;
        let energy_radiated = waste - heat_released;
        Transaction {
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
