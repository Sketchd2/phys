//! The generated program: the rule the engine wrote down for placing a thing's
//! bodies.
//!
//! # Why a program is generated rather than selected
//!
//! `docs/PLAY.md` D11 finds the largest standing axiom violation in the
//! codebase: [`crate::morph::Program`] was a species table with a renderer per
//! variant, so the engine knew what a tree was, what a tower was, and what a
//! town was, and each of them had its geometry written into a `match`. D15
//! already answered that for a composite — a wooden box is **generated**, by
//! assessing what the parts are, and no variant describes it — and Phase 4
//! carries the same answer to everything else.
//!
//! So this module holds the *rules*, and nothing in it is a species. A recipe
//! is produced by an analysis of what is actually there:
//!
//! * somebody **placed** parts, and the engine wrote down where they are;
//! * something is **growing** into an occluded field and branches, whether it
//!   is a tree in air or a coral in water — the fluid is measured, not named;
//! * something is being **laid** course on course, or framed storey on storey,
//!   or **subdivided** across a plan, which are the space-filling rules D11
//!   keeps as *habits* rather than as species;
//! * something is a piece of a surface and **tiles** it.
//!
//! Six habits serve every kind of thing the engine has, and a seventh kind of
//! thing needs no seventh variant — it needs a measurement that picks one of
//! these and parameterises it. That is the reduction D11 asks for.
//!
//! # A recipe is a rule, not a frozen outcome
//!
//! `CLAUDE.md`: "A stored rule stays a function of its conditions and never a
//! frozen outcome, so a drought still reaches every tree in the region it
//! touches." A branching recipe stores a taper, a split count and a lean — the
//! shape of the rule — and derives its *size* from the mass the thing has
//! accumulated, every time it is asked. A planned recipe states metres, because
//! a planned thing's size is a decision somebody made rather than a consequence
//! of how much it has grown, and it is regenerated when that decision changes.
//!
//! # It emits metres
//!
//! **This is the fix for a measured defect.** A generated skeleton used to be
//! emitted in units of the structure's own extent and then rescaled by
//! `sampler::radius_scale` until `summarise` reported the node's radius back —
//! which meant the size a program *stated* and the size it was *drawn* at were
//! two different numbers. Measured, as the factor the drawing was inflated by:
//! tree 1.44, coral 1.53, tower 0.89, wall 0.76, settlement 0.95. A tree was
//! drawn 44% taller than `tree_height` said it was.
//!
//! A recipe states where its parts are, in metres, and the sampler places them
//! there. `matter.radius` for a structure is then the bounding radius of what
//! was actually drawn, which is one meaning rather than two.

use crate::assembly::Assembly;
use crate::math::{v3, Quat, Vec3};
use crate::morph::{Skeleton, FREE_ALL, FREE_Y, FREE_Z, NO_SUPPORT};

/// How much of a thing there is, at the moment it is being drawn.
///
/// The conditions a recipe is a function of. Separated from the recipe so the
/// recipe cannot accidentally store one of them and become a frozen outcome.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Growth {
    /// Structural mass accumulated, kg.
    pub built: f64,
    /// Completion fraction for a planned recipe, 0..1.
    pub progress: f64,
}

impl Growth {
    pub fn new(built: f64, progress: f64) -> Growth {
        Growth { built: built.max(0.0), progress: progress.clamp(0.0, 1.0) }
    }
}

/// The generated program. See the module documentation.
#[derive(Debug, Clone, PartialEq)]
pub enum Recipe {
    /// Parts somebody placed, and the joins that hold them. `docs/PLAY.md` D15.
    Placed(Assembly),
    /// Recursive branching under a transport field: a tree, a coral, a delta.
    Branching(Branching),
    /// Laid course on course: a wall, brickwork, strata.
    Coursed(Coursed),
    /// A frame of columns, beams and floors: anything with storeys.
    Framed(Framed),
    /// A tiled surface: ground.
    Tiled(Tiled),
    /// A subdivided plane: a settlement, cracked mud, leaf venation.
    Subdivided(Subdivided),
    /// A packing of grains: a melt that froze.
    Granular(Granular),
}

impl Recipe {
    /// Place the bodies, in metres, in the node's own frame.
    ///
    /// `field` is what has happened to this thing that its own rule does not
    /// describe — `docs/PLAY.md` §5's deviations, superposed over the derived
    /// baseline. Only [`Tiled`] reads it, because §5.7 decides that terrain
    /// carries a field and a structure carries an edit list, and the edit list
    /// is applied one level up in `Morphology::render`.
    pub fn render(&self, budget: usize, g: Growth, field: &[crate::erode::Deviation]) -> Skeleton {
        match self {
            Recipe::Placed(a) => a.render(),
            Recipe::Branching(b) => b.render(budget, g),
            Recipe::Coursed(c) => c.render(budget, g),
            Recipe::Framed(f) => f.render(budget, g),
            Recipe::Tiled(t) => t.render(budget, field),
            Recipe::Subdivided(s) => s.render(budget, g),
            Recipe::Granular(x) => x.render(budget),
        }
    }

    /// The radius the node claims, metres.
    ///
    /// Closed-form per habit rather than measured off a render, because growth
    /// asks for it every frame on every structure in the world.
    ///
    /// **This is the size the recipe states, and the sampler no longer argues
    /// with it.** What it means differs by habit and that is deliberate: for a
    /// parts list it is the equivalent uniform sphere of the parts, which is
    /// what `matter.radius` means everywhere else in the engine and what every
    /// composite has been measured as since D15; for a generated habit it is
    /// the half-diagonal of what the habit lays out, which is what a structure
    /// has always claimed. Reconciling the two is a separate question from the
    /// one this module fixes — see the module documentation — and moving a
    /// structure's radius moves what `Formation::of_matter` measures its
    /// packing as, which moves its strength.
    pub fn extent(&self, g: Growth) -> f64 {
        match self {
            Recipe::Placed(a) => a.extent(),
            Recipe::Branching(b) => b.extent(g),
            Recipe::Coursed(c) => c.extent(g),
            Recipe::Framed(f) => f.extent(g),
            Recipe::Tiled(t) => t.extent(),
            Recipe::Subdivided(s) => s.extent(g),
            Recipe::Granular(x) => x.extent(),
        }
    }

    /// Overall height, metres — the dimension a thing is measured by when
    /// somebody asks how tall it is.
    pub fn height(&self, g: Growth) -> f64 {
        match self {
            Recipe::Branching(b) => b.height(g),
            Recipe::Coursed(c) => c.height(g),
            Recipe::Framed(f) => f.height(g),
            Recipe::Tiled(t) => t.depth,
            Recipe::Subdivided(s) => s.height,
            Recipe::Granular(x) => x.height(),
            Recipe::Placed(a) => 2.0 * a.bound(),
        }
    }

    /// Bulk density, kg/m^3 — what the analysis measured and turned this
    /// thing's mass into its size with.
    ///
    /// **D11's first column, gone.** `Program::density` was one value per
    /// species; this is the number the material actually has, carried on the
    /// rule that used it, so the sampler's volume check and the recipe cannot
    /// disagree. They did: a patch generated against a measured cellulose and
    /// checked against a tabulated rock came out three times too shallow,
    /// because the correction was making up the difference between two
    /// densities that were describing the same thing.
    pub fn density(&self) -> f64 {
        match self {
            Recipe::Placed(a) => {
                let v: f64 = a.parts.iter().map(|p| p.volume()).sum();
                let m = a.mass();
                if v > 0.0 && m > 0.0 {
                    m / v
                } else {
                    0.0
                }
            }
            Recipe::Branching(b) => b.density,
            Recipe::Coursed(c) => c.density,
            Recipe::Framed(f) => f.density,
            Recipe::Tiled(t) => t.density,
            Recipe::Subdivided(s) => s.density,
            Recipe::Granular(x) => x.density,
        }
    }

    /// The volume of solid this recipe actually lays down, m^3.
    ///
    /// **What a node's packing should be measured against.** A node's bulk
    /// density is its mass over the sphere it claims, which is right for a
    /// node that *is* the material and wrong for one that is an arrangement of
    /// it: a crate of solid planks reads as 4.8% packed against its bounding
    /// sphere, and a patch of ground — a slab much wider than it is deep —
    /// reads at 60 kg/m^3 where the rock it is made of is 2600.
    ///
    /// A recipe knows better, because the volume it lays down is the volume it
    /// turned the mass into a size with. A parts list states it outright; a
    /// habit's is its own mass over its own density, which is the same number
    /// the analysis used.
    pub fn solid_volume(&self, g: Growth) -> f64 {
        match self {
            Recipe::Placed(a) => a.parts.iter().map(|p| p.volume()).sum(),
            _ => {
                let d = self.density();
                if d > 0.0 {
                    g.built.max(0.0) / d
                } else {
                    0.0
                }
            }
        }
    }

    /// The parts, where this recipe is a parts list.
    pub fn placed(&self) -> Option<&Assembly> {
        match self {
            Recipe::Placed(a) => Some(a),
            _ => None,
        }
    }

    pub fn placed_mut(&mut self) -> Option<&mut Assembly> {
        match self {
            Recipe::Placed(a) => Some(a),
            _ => None,
        }
    }

    /// Byte cost of the rule, for the detail budget.
    pub fn state_bytes(&self) -> usize {
        match self {
            Recipe::Placed(a) => a.state_bytes(),
            _ => std::mem::size_of::<Recipe>(),
        }
    }

    /// A short name for logs. Not read by any physics.
    pub fn habit(&self) -> &'static str {
        match self {
            Recipe::Placed(_) => "placed",
            Recipe::Branching(_) => "branching",
            Recipe::Coursed(_) => "coursed",
            Recipe::Framed(_) => "framed",
            Recipe::Tiled(_) => "tiled",
            Recipe::Subdivided(_) => "subdivided",
            Recipe::Granular(_) => "granular",
        }
    }
}

// ---------------------------------------------------------------------------
// branching
// ---------------------------------------------------------------------------

/// Recursive branching under elastic self-similarity.
///
/// D12: branching is what transport into an occluded field looks like, so this
/// is one habit and not two. A tree and a coral differ in the numbers below and
/// in nothing else — and those numbers are measured, because the fluid a thing
/// grows in is its node's own mixture and the load it is proportioned against
/// is the flow it has met.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Branching {
    /// Cross-section carried into each child, as a fraction of the parent's.
    pub taper: f64,
    /// How many children a branch point has.
    pub splits: u8,
    pub lean: f64,
    pub twist: f64,
    pub spread: f64,
    /// Slenderness: how tall this individual stands for its mass.
    pub slenderness: f64,
    /// Bulk density of the structure, kg/m^3 — measured, not a column.
    pub density: f64,
    /// Volume of the smallest thing this habit produces, m^3. The constant in
    /// `height ~ (V/k)^(1/4)`, which is McMahon's buckling criterion and the
    /// same quarter power that runs through all of allometry.
    pub allometry: f64,
}

/// Mass a structure starts from, kg. Below this there is nothing to draw.
pub const SEED_MASS: f64 = 1e-4;

impl Branching {
    /// Height from mass, by elastic self-similarity.
    ///
    /// McMahon's buckling criterion gives trunk radius proportional to
    /// `height^1.5`, so volume goes as `height^4` and height as `volume^(1/4)`.
    pub fn height(&self, g: Growth) -> f64 {
        let v = g.built.max(0.0) / self.density.max(1e-9);
        (v / self.allometry.max(1e-30)).max(0.0).powf(0.25) * self.slenderness
    }

    /// How many orders of branching this much mass supports.
    fn depth(&self, g: Growth) -> usize {
        (2.0 + (g.built.max(SEED_MASS) / SEED_MASS).log10() * 1.4).clamp(1.0, 8.0) as usize
    }

    /// Half the height it stands at — what a thing that grew into a field
    /// claims as its radius.
    pub fn extent(&self, g: Growth) -> f64 {
        0.5 * self.height(g)
    }

    /// Light-intercepting area, m^2. Crown projection, not total leaf area —
    /// what limits a thing that grows into a field is the ground it shades.
    pub fn capture_area(&self, g: Growth) -> f64 {
        let h = self.height(g);
        let crown = 0.3 * h * (0.75 + 0.55 * self.spread.clamp(0.0, 1.0));
        std::f64::consts::PI * crown * crown
    }

    /// Segments laid down breadth-first with a mass budget.
    ///
    /// Each level takes `taper^3` of its parent's cross-section, so the
    /// structure obeys da Vinci's rule — total cross-section is preserved
    /// across a branch point — and therefore looks like a tree rather than like
    /// a fractal.
    ///
    /// Breadth-first, not depth-first. With a stack the budget is spent
    /// rendering one branch down to its finest twigs while the rest is simply
    /// absent, which looks broken at any budget below saturation and makes the
    /// level of detail meaningless.
    pub fn render(&self, budget: usize, g: Growth) -> Skeleton {
        let unit = 0.5 * self.height(g);
        let mut sk = Skeleton::with_capacity(budget);
        let max_depth = self.depth(g);
        let splits = self.splits.max(1) as usize;

        struct Seg {
            base: Vec3,
            dir: Vec3,
            len: f64,
            rad: f64,
            depth: usize,
            id: u32,
            parent: u32,
        }
        let mut queue = std::collections::VecDeque::from(vec![Seg {
            base: v3(0.0, 0.0, -1.0),
            dir: v3(self.lean, self.lean * 0.5, 1.0).unit(),
            len: 0.55,
            rad: 0.055,
            depth: 0,
            id: 0,
            parent: NO_SUPPORT,
        }]);
        let mut next_id = 1u32;
        let mut emitted: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();

        while let Some(s) = queue.pop_front() {
            if sk.len() >= budget {
                break;
            }
            let tip = s.base + s.dir.scale(s.len);
            let support = if s.parent == NO_SUPPORT {
                NO_SUPPORT
            } else {
                match emitted.get(&s.parent) {
                    Some(&idx) => idx,
                    // The parent fell outside the budget, so this segment has
                    // nothing to hang from and is not emitted.
                    None => continue,
                }
            };
            emitted.insert(s.id, sk.len() as u32);
            sk.push_segment(
                s.base.scale(unit),
                tip.scale(unit),
                s.rad * s.rad * s.len,
                s.rad * unit,
                support,
                s.id,
            );
            if s.depth >= max_depth {
                continue;
            }
            let child_rad = s.rad * self.taper;
            let child_len = s.len * (0.62 + 0.12 * self.spread);
            for k in 0..splits {
                let phi = self.twist
                    + k as f64 * std::f64::consts::TAU / splits as f64
                    + s.depth as f64 * 0.7;
                let axis = v3(phi.cos(), phi.sin(), 0.0);
                let dir = (s.dir + axis.scale(self.spread)).unit();
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
        sk
    }
}

// ---------------------------------------------------------------------------
// coursed
// ---------------------------------------------------------------------------

/// Laid course on course: a wall, brickwork, a stratum.
///
/// Planned, so its size is stated rather than grown into: somebody decided how
/// long the wall is. `progress` is how much of it has been laid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coursed {
    /// Metres, the finished thing.
    pub length: f64,
    pub height: f64,
    pub thickness: f64,
    /// Height of one course, metres — the thickness the process lays down at a
    /// time. **Not a count**: how many courses are drawn depends on how finely
    /// anyone is looking, and a rule that stored the count would be storing a
    /// level of detail.
    pub course: f64,
    /// Bulk density, kg/m^3 — what the analysis turned the mass into a size
    /// with, so the sampler's volume check and the recipe agree by
    /// construction.
    pub density: f64,
}

impl Coursed {
    pub fn height(&self, _g: Growth) -> f64 {
        self.height
    }

    pub fn extent(&self, g: Growth) -> f64 {
        let laid = (self.length * g.progress.max(1e-3)).max(1e-9);
        0.5 * (laid * laid + self.height * self.height + self.thickness * self.thickness).sqrt()
    }

    /// Each block rests on the course beneath, offset by half a block on
    /// alternate courses — which is what stops a wall being a set of
    /// independent vertical columns of brick, and is why the load path runs
    /// diagonally down to the footing.
    ///
    /// **A block is a box.** A capsule over the same two endpoints is a bead,
    /// and a row of beads has a scallop between them that something walking
    /// along the top falls into. The wall's thickness is the one axis the
    /// sampler's cross-section correction may move: it may not get longer or
    /// the courses stop meeting, and it may not get taller or they stop
    /// stacking.
    pub fn render(&self, budget: usize, g: Growth) -> Skeleton {
        // **How many courses is a level of detail, not a shape.** A block is
        // about twice as long as its course is tall, which is what a brick is;
        // that fixes the aspect, and the budget fixes how many of them there
        // are. Coarsened, a wall is fewer bigger blocks laid the same way.
        let want_courses = (self.height / self.course.max(1e-9)).max(1.0);
        let want_per = (self.length / (2.0 * self.course.max(1e-9))).max(2.0);
        let fit = (budget as f64 / (want_courses * want_per)).min(1.0).sqrt();
        let courses = (want_courses * fit).round().max(1.0) as usize;
        let mut sk = Skeleton::with_capacity(budget);
        let per_course = ((want_per * fit).round() as usize).clamp(2, 256);
        let laid = g.progress.clamp(0.0, 1.0) * courses as f64;
        let block = self.length / per_course as f64;
        let course_half = 0.5 * self.height / courses as f64;
        let mut prev_start = 0u32;
        let mut prev_count = 0usize;
        let mut site = 0u32;
        for c in 0..courses {
            let complete = (laid - c as f64).clamp(0.0, 1.0);
            if complete <= 0.0 || sk.len() >= budget {
                break;
            }
            let z = -0.5 * self.height + self.height * (c as f64 + 0.5) / courses as f64;
            let blocks = ((per_course as f64 * complete).round() as usize).max(1);
            let offset = if c % 2 == 0 { 0.0 } else { 0.5 };
            let start = sk.len() as u32;
            for b in 0..blocks {
                let x0 = -0.5 * self.length + block * (b as f64 + offset);
                let support = if c == 0 || prev_count == 0 {
                    NO_SUPPORT
                } else {
                    prev_start + (b.min(prev_count - 1)) as u32
                };
                sk.push_box(
                    v3(x0, 0.0, z),
                    v3(x0 + block, 0.0, z),
                    v3(0.5 * block, 0.5 * self.thickness, course_half),
                    FREE_Y,
                    1.0,
                    0.5 * self.thickness,
                    support,
                    site,
                );
                site += 1;
                if sk.len() >= budget {
                    break;
                }
            }
            prev_start = start;
            prev_count = sk.len() - start as usize;
        }
        if sk.is_empty() {
            let t = 0.5 * self.thickness.max(1e-6);
            sk.push_segment(v3(-t, 0.0, 0.0), v3(t, 0.0, 0.0), 1.0, t, NO_SUPPORT, 0);
        }
        sk
    }
}

// ---------------------------------------------------------------------------
// framed
// ---------------------------------------------------------------------------

/// Columns, beams, bracing and a floor at every storey.
///
/// Planned construction: the target is known in advance, so a half-built frame
/// is the finished design masked by the completion fraction with the topmost
/// storey partial.
///
/// The members are real segments — columns spanning a storey, beams spanning
/// between column heads. An earlier version emitted each element as a point
/// with a nominal radius, which was adequate while the parts were only drawn
/// and became nonsense the moment they had to carry load: a zero-length member
/// has no bending stiffness, and the density correction inflated its radius
/// until the frame rendered as a smear of vertical streaks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Framed {
    /// Bulk density, kg/m^3. See [`Coursed::density`].
    pub density: f64,
    pub floors: u16,
    /// Height of one storey, metres.
    pub storey: f64,
    /// Half-width of the plan, metres.
    pub side: f64,
    /// Thickness of a floor, metres.
    pub plate: f64,
}

/// Columns per storey. Four corners is the smallest frame that is a frame.
pub const COLUMNS: usize = 4;

impl Framed {
    pub fn height(&self, g: Growth) -> f64 {
        self.storey * self.floors.max(1) as f64 * g.progress.clamp(0.0, 1.0).max(1e-3)
    }

    pub fn extent(&self, g: Growth) -> f64 {
        let h = self.height(g);
        // The plan's corners sit at `side` from the axis, so the plan's
        // diagonal is `2 side`.
        0.5 * (h * h + 4.0 * self.side * self.side).sqrt()
    }

    pub fn render(&self, budget: usize, g: Growth) -> Skeleton {
        let floors = self.floors.max(1) as usize;
        let mut sk = Skeleton::with_capacity(budget);
        let total = self.storey * floors as f64;
        let built_floors = g.progress.clamp(0.0, 1.0) * floors as f64;
        let side = self.side;
        let corner = |c: usize, z: f64| {
            let a = std::f64::consts::TAU * c as f64 / COLUMNS as f64 + std::f64::consts::FRAC_PI_4;
            v3(side * a.cos(), side * a.sin(), z)
        };

        let mut below = [NO_SUPPORT; COLUMNS];
        let mut site = 0u32;
        for f in 0..floors {
            let complete = (built_floors - f as f64).clamp(0.0, 1.0);
            if complete <= 0.0 || sk.len() + 2 * COLUMNS + 1 > budget {
                break;
            }
            let z0 = -0.5 * total + self.storey * f as f64;
            let z1 = z0 + self.storey * complete;

            let mut here = [NO_SUPPORT; COLUMNS];
            for (c, slot) in here.iter_mut().enumerate() {
                *slot = sk.len() as u32;
                sk.push_segment(corner(c, z0), corner(c, z1), 1.0, 0.03 * side, below[c], site);
                site += 1;
            }
            if complete >= 0.999 {
                for c in 0..COLUMNS {
                    let n = (c + 1) % COLUMNS;
                    sk.push_segment(corner(c, z1), corner(n, z1), 0.7, 0.02 * side, here[c], site);
                    site += 1;
                }
                // **The floor.** A frame of columns and beams is a wireframe:
                // there is nothing between the beams, so anything set down on a
                // storey falls through it to the footing. `docs/PLAY.md` Phase
                // 4 asks the generators that lay down flat things to emit them,
                // and a floor plate is the flat thing a frame lays down.
                //
                // It is a *plate* and not a beam, and `solvers::frame` knows
                // the difference: a plate spanning its bay as a beam came out
                // with a slenderness of 827 and the static solve stopped
                // converging. See `Skeleton::push_plate`.
                let plan = side * std::f64::consts::FRAC_1_SQRT_2;
                sk.push_plate(
                    v3(0.0, 0.0, z1 - 0.5 * self.plate),
                    v3(plan, plan, 0.5 * self.plate),
                    2.0,
                    &here,
                    site,
                );
                site += 1;
                // Cross-bracing between adjacent columns, and diagonally to the
                // storey below. This is what makes a frame a frame rather than
                // a stack of posts — and it makes the structure statically
                // indeterminate, so the redundant solver has a generated case
                // to work on and not only a test rig.
                for c in 0..COLUMNS {
                    let n = (c + 1) % COLUMNS;
                    sk.tie(here[c], here[n], 0.33);
                    if below[c] != NO_SUPPORT {
                        sk.tie(here[c], below[n], 0.25);
                    }
                }
            }
            below = here;
        }
        if sk.is_empty() {
            let h = 0.1 * self.storey.max(1e-6);
            sk.push_segment(v3(0.0, 0.0, -h), v3(0.0, 0.0, h), 1.0, 0.03 * side, NO_SUPPORT, 0);
        }
        sk
    }
}

// ---------------------------------------------------------------------------
// tiled
// ---------------------------------------------------------------------------

/// A square of surface: ground.
///
/// Stated in metres, because a patch's side is the tiling's rather than a
/// consequence of how much it weighs. A standalone patch has its side derived
/// from its own mass when the recipe is generated; a patch of a planet has it
/// from the parameterisation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tiled {
    /// Bulk density, kg/m^3. See [`Coursed::density`].
    pub density: f64,
    /// Radius of the sphere this patch is a piece of, metres, or zero for a
    /// flat patch that is a thing on its own.
    pub sphere: f64,
    /// Which face of the cube, or [`WHOLE_BALL`] for the ball itself.
    pub face: u8,
    /// How many times the face has been divided to reach this patch.
    pub level: u8,
    /// Where on the face, at that level. Sixty-four bits because the address
    /// is exact integer arithmetic and a descent to a millimetre on an Earth
    /// is twelve levels of eight, which is `8^12` — past what thirty-two bits
    /// hold, and the overflow is a panic rather than a wrong answer because
    /// `[profile.test]` keeps overflow checks on.
    pub u: u64,
    pub v: u64,
    /// How many cells across this patch divides into. Fixed by the recipe
    /// rather than by whoever is looking, because the address has to mean the
    /// same thing at every level of detail.
    pub cells: u8,
    /// Side of the square, metres.
    pub side: f64,
    /// How deep the slab goes below its own mean surface, metres.
    pub depth: f64,
    /// The height field, as the amplitudes and phases of a short sum of
    /// sinusoids. See [`Tiled::surface_height`].
    pub relief: [f32; 8],
}

impl Tiled {
    pub fn extent(&self) -> f64 {
        if self.is_ball() {
            // The ball is the size it is. Its cells are pieces of its surface
            // and its radius is not theirs.
            return self.sphere;
        }
        if self.on_sphere() {
            // **A patch on a sphere is curved, and its bound is not a flat
            // slab's.** The half-diagonal of a face laid out flat is 6.59e6 m
            // on an Earth — larger than the planet the face is part of — so a
            // promoted face came out bigger than its own parent. The honest
            // answer is the furthest of its own corners from its own centre of
            // mass, which curvature makes smaller rather than larger.
            let com = self.centre_of_mass_from_planet();
            let (a, b, half) = self.face_span();
            let mut worst: f64 = 0.0;
            for (da, db) in [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)] {
                let dir = cube_to_sphere(self.face, a + da * half, b + db * half);
                for depth in [0.0, self.depth] {
                    worst = worst.max((dir.scale(self.sphere - depth) - com).norm());
                }
            }
            return worst.max(1e-30) + self.relief_amplitude();
        }
        let r = self.relief_amplitude();
        0.5 * (2.0 * self.side * self.side + (self.depth + 2.0 * r).powi(2)).sqrt()
    }

    fn gene(&self, i: usize, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.relief[i % 8] as f64
    }

    /// How far the surface departs from its own mean, metres.
    ///
    /// A fraction of the slab's depth rather than of its width: ground varies
    /// by a part of how deep it is, not by a part of how wide it is.
    pub fn relief_amplitude(&self) -> f64 {
        self.depth * self.gene(0, 0.15, 0.55)
    }

    /// A deterministic height at a point on the patch, in units of the patch's
    /// own half-side.
    ///
    /// Not noise from a library and not a table: a short sum of sinusoids whose
    /// frequencies and phases come from the recipe, which came from the node's
    /// address. So the same patch of ground is the same shape every time it is
    /// regenerated, two patches differ, and nothing has to be stored. The
    /// octaves halve in amplitude and roughly double in frequency, which is
    /// what makes a landscape read as a landscape rather than as a sine wave:
    /// most of the relief is in the largest feature and the rest is detail on
    /// it.
    pub fn surface_height(&self, x: f64, y: f64) -> f64 {
        let tilt = self.gene(1, -0.15, 0.15);
        let mut h = tilt * x;
        let mut amp = 1.0;
        let mut freq = self.gene(2, 1.1, 2.2);
        for o in 0..4 {
            let px = self.gene(3 + o % 4, 0.0, std::f64::consts::TAU);
            let py = self.gene((5 + o) % 8, 0.0, std::f64::consts::TAU);
            h += amp * ((freq * x + px).sin() * (freq * y + py).cos());
            amp *= 0.5;
            freq *= 2.07;
        }
        h
    }

    /// A grid of columns, each standing on nothing.
    ///
    /// `NO_SUPPORT` is exactly right here and is not a shortcut. It means "this
    /// part is anchored, load stops here", which is what bedrock is — the
    /// terrain is what everything else's load path terminates in.
    ///
    /// **A square prism, not a cylinder.** Touching at their midlines is not
    /// the same as tiling: round columns on a square grid leave a gap at every
    /// corner of it, and something walking across the patch drops into each
    /// one. The cell is square, so the column that fills it is.
    pub fn render(&self, budget: usize, field: &[crate::erode::Deviation]) -> Skeleton {
        if self.on_sphere() {
            return self.render_on_sphere(field);
        }
        let mut sk = Skeleton::with_capacity(budget);
        let n = ((budget as f64).sqrt().floor() as usize).clamp(2, 64);
        let step = self.side / n as f64;
        let half = 0.5 * self.side;
        let amp = self.relief_amplitude();
        // Relief about the patch's *own mean*, measured over the same grid it
        // is drawn on. Without that the mean height is an arbitrary offset and
        // the drawn volume is not the stated one.
        let mut mean = 0.0;
        for i in 0..n {
            for j in 0..n {
                let x = -1.0 + 2.0 * (i as f64 + 0.5) / n as f64;
                let y = -1.0 + 2.0 * (j as f64 + 0.5) / n as f64;
                mean += self.surface_height(x, y);
            }
        }
        mean /= (n * n) as f64;

        let mut site = 0u32;
        for i in 0..n {
            for j in 0..n {
                if sk.len() >= budget {
                    break;
                }
                let u = -1.0 + 2.0 * (i as f64 + 0.5) / n as f64;
                let v = -1.0 + 2.0 * (j as f64 + 0.5) / n as f64;
                // The derived baseline, plus whatever has happened here that
                // the baseline does not describe. A field superposes, so a
                // thousand footprints cost what one costs — §5.7.
                let mut top = amp * (self.surface_height(u, v) - mean);
                for d in field {
                    top += d.height_at(u, v, half);
                }
                let base = v3(u * half, v * half, -self.depth);
                let tip = v3(u * half, v * half, top);
                let len = (top + self.depth).max(1e-9);
                // The cell is the cell: a column may get deeper and may not get
                // narrower, or the patch stops tiling and the ground has holes
                // in it.
                sk.push_box(
                    base,
                    tip,
                    v3(0.5 * step, 0.5 * step, 0.5 * len),
                    FREE_Z,
                    step * step * len,
                    0.5 * step,
                    NO_SUPPORT,
                    site,
                );
                site += 1;
            }
        }
        sk
    }
}

// ---------------------------------------------------------------------------
// subdivided
// ---------------------------------------------------------------------------

/// Plots on a street grid: a settlement, cracked mud, leaf venation.
///
/// This habit deliberately does not describe walls, floors or roofs. It
/// describes *where the buildings are*, and each plot it emits is a body that
/// can be promoted into a node of its own — which is the whole shape of world
/// generation here: a recipe's output bodies are the next level's nodes, and
/// the ladder from a moon to a room is recipes all the way down rather than one
/// generator that knows about everything.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subdivided {
    /// Bulk density, kg/m^3. See [`Coursed::density`].
    pub density: f64,
    /// Side of the plan, metres.
    pub side: f64,
    /// How tall the things standing on it are, metres.
    pub height: f64,
    /// Fraction of a block taken by the street rather than the plot.
    pub street: f64,
    /// Per-instance variation, so a town is not uniform and is not different
    /// every time it is drawn.
    pub variation: f32,
}

impl Subdivided {
    pub fn extent(&self, _g: Growth) -> f64 {
        0.5 * (2.0 * self.side * self.side + self.height * self.height).sqrt()
    }

    /// Ordered from the middle outward, because towns fill in from their
    /// centre: a partly built one is a core with edges missing rather than a
    /// scatter of lone houses.
    pub fn render(&self, budget: usize, g: Growth) -> Skeleton {
        let mut sk = Skeleton::with_capacity(budget);
        let blocks = ((budget as f64).sqrt().floor() as usize).clamp(2, 16);
        let step = self.side / blocks as f64;
        let plot = 0.5 * step * (1.0 - self.street.clamp(0.0, 0.9));
        let built = g.progress.clamp(0.0, 1.0);
        let half = 0.5 * self.side;

        let mut plots: Vec<(f64, usize, usize)> = Vec::with_capacity(blocks * blocks);
        for i in 0..blocks {
            for j in 0..blocks {
                let x = -half + step * (i as f64 + 0.5);
                let y = -half + step * (j as f64 + 0.5);
                plots.push((x * x + y * y, i, j));
            }
        }
        plots.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let wanted = ((plots.len() as f64) * built).round() as usize;

        let mut site = 0u32;
        for (_, i, j) in plots.into_iter().take(wanted.min(budget)) {
            if sk.len() >= budget {
                break;
            }
            let x = -half + step * (i as f64 + 0.5);
            let y = -half + step * (j as f64 + 0.5);
            let vary = frac(
                (i as f64 * 12.9898 + j as f64 * 78.233 + self.variation as f64 * 43.5).sin()
                    * 43758.5453,
            );
            let h = self.height * (0.6 + 1.8 * vary);
            let base = v3(x, y, -self.height);
            let tip = v3(x, y, -self.height + h);
            sk.push_box(
                base,
                tip,
                v3(plot, plot, 0.5 * h),
                FREE_ALL,
                plot * plot * h,
                plot,
                NO_SUPPORT,
                site,
            );
            site += 1;
        }
        sk
    }
}

/// Fractional part, for the small deterministic hashes the flat habits use.
fn frac(x: f64) -> f64 {
    x - x.floor()
}

// ---------------------------------------------------------------------------
// the analysis
// ---------------------------------------------------------------------------

/// A space-filling rule, before its numbers are filled in.
///
/// **Not a species.** D11: "Branching, coursed masonry and a subdivided street
/// grid are genuinely different space-filling rules, and asserting they
/// collapse into one would be an over-claim. What they do reduce to is a small
/// set of **habits** — selected and parameterised by the genome instead of
/// named by an enum." This is that set. A tree and a coral are one habit with
/// different numbers; a wall and a stratum are one habit; a settlement and
/// cracked mud are one habit.
///
/// What selects a habit is what the thing is *doing*, and where an actor has
/// placed a design the habit comes from the design rather than from here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Habit {
    /// Growing into an occluded transport field.
    Branching,
    /// Being laid course on course.
    Coursed,
    /// Being framed, storey on storey.
    Framed,
    /// A piece of a surface.
    Tiled,
    /// A plan divided into plots.
    Subdivided,
}

/// Air at sea level, kg/m^3 — the reference the measured fluid is read against.
const AIR_DENSITY: f64 = 1.225;

/// Write down the rule for placing a thing's bodies.
///
/// **This is the analysis, and it is the only place a habit's numbers come
/// from.** Everything it uses is measured or derived: the mass the thing has,
/// the material it is made of — which is read off the node's own mixture — and
/// the fluid it is standing in, which is the node's own mixture again.
///
/// D11's own example is the one that proves the point: "A coral is in water
/// because its node's mixture is water; nobody tells it." A branching thing
/// growing in a dense fluid closes up and splits more, because the fluid both
/// carries more load and delivers more of what it is building with; one in air
/// opens out. That is the whole of the difference between the two branching
/// species the table used to hold, and it is a measurement now.
pub fn generate(
    habit: Habit,
    genome: &[f32; 8],
    mass: f64,
    env: &crate::morph::Environment,
    material: &crate::material::Material,
) -> Recipe {
    let gene = |i: usize, lo: f64, hi: f64| lo + (hi - lo) * genome[i % 8] as f64;
    let density = material.density.max(1e-6);
    let mass = mass.max(0.0);
    match habit {
        Habit::Branching => {
            // How dense the fluid is, against air. Zero for air, one for water.
            let wet = (env.fluid_density.max(AIR_DENSITY) / AIR_DENSITY).log10()
                / (1025.0f64 / AIR_DENSITY).log10();
            let wet = wet.clamp(0.0, 1.0);
            Recipe::Branching(Branching {
                // A thing carrying its own weight in air tapers hard; one held
                // up by the fluid it is in barely has to.
                taper: 0.62 + 0.10 * wet,
                // More of what it is building with arrives from more
                // directions, so a branch point has more children.
                splits: (3.0 + wet.round()) as u8,
                lean: gene(2, -0.12, 0.12),
                twist: gene(3, 0.0, std::f64::consts::TAU),
                spread: gene(4, 0.45, 0.85),
                slenderness: gene(0, 0.8, 1.25),
                density,
                // McMahon's buckling criterion, calibrated once so that a
                // one-tonne thing of wood's density stands about fifteen
                // metres. A universal of the allometry rather than a column:
                // the same number gives a coral its proportions, because what
                // differs between them is the density and the fluid.
                allometry: 3.3e-5,
            })
        }
        Habit::Coursed => {
            // The course is the material's own deposition increment: a mason
            // lays a course of the thickness the process puts down at a time,
            // which is the very same number D14's flaw scale is. Nothing here
            // states a block size.
            // The course is the material's own deposition increment: a mason
            // lays a course of the thickness the process puts down at a time,
            // which is the very same number D14's flaw scale is. Nothing here
            // states a block size. Bounded below by a millimetre, because a
            // material whose increment is a cell wall is being used as
            // something it was not laid down as, and a wall of 30 um courses
            // is a level of detail nothing can draw.
            let course = material.flaw_size.clamp(1e-3, 1.0);
            let volume = mass / density;
            // A wall is much longer than it is tall and much taller than it is
            // thick. The proportions are the genome's; the size is the mass's.
            let slender = gene(0, 2.0, 5.0);
            let stoutness = gene(1, 0.12, 0.30);
            // volume = length * height * thickness, with
            // length = slender * height and thickness = stoutness * height.
            let height = (volume / (slender * stoutness)).max(0.0).cbrt().max(1e-6);
            Recipe::Coursed(Coursed {
                length: slender * height,
                height,
                thickness: stoutness * height,
                course: course.min(height),
                density,
            })
        }
        Habit::Framed => {
            let storey = gene(2, 2.6, 4.0);
            let floors = (6.0 + gene(0, 0.0, 34.0)).round().max(1.0);
            // The plan follows from the mass: a frame of this much material,
            // this many storeys tall, covers this much ground. The plate is a
            // fixed fraction of a storey, which is what a floor is.
            let volume = mass / density;
            let plate = 0.06 * storey;
            // Each storey costs its own floor plus its columns and beams; the
            // floor dominates, so the plan is what the volume buys.
            let area = (volume / (floors * plate)).max(1e-6);
            let side = (area.sqrt() * std::f64::consts::FRAC_1_SQRT_2).max(1e-3);
            Recipe::Framed(Framed { density, floors: floors as u16, storey, side, plate })
        }
        Habit::Tiled => {
            // A patch is much wider than it is deep, and this is how much.
            let volume = mass / density;
            let side = (volume / SLAB_ASPECT).max(0.0).cbrt().max(1e-6);
            let mut relief = [0.0f32; 8];
            relief.copy_from_slice(genome);
            Recipe::Tiled(Tiled {
                density,
                sphere: 0.0,
                face: 0,
                level: 0,
                u: 0,
                v: 0,
                cells: 2,
                side,
                depth: side * SLAB_ASPECT,
                relief,
            })
        }
        Habit::Subdivided => {
            // A town is an *area*, not a volume, and that is the difference
            // that matters. Deriving its size from its mass the way a solid
            // body's is derived — the cube root of a volume — describes a town
            // cast as one lump of concrete, and gives a six-thousand-tonne
            // settlement fifty metres across with fifty-metre buildings in it.
            let height = gene(2, 8.0, 16.0);
            // What a town actually is: things a few storeys tall covering some
            // of the ground, and streets and yards for the rest.
            let street = gene(0, 0.22, 0.38);
            let cover = (1.0 - street) * (1.0 - street);
            let areal = density * height * cover;
            let side = (mass / areal.max(1e-9)).max(0.0).sqrt().max(1e-6);
            Recipe::Subdivided(Subdivided { density, side, height, street, variation: genome[1] })
        }
    }
}

/// Depth of a patch of ground as a fraction of its side. Terrain is much wider
/// than it is deep, and this is how much.
pub const SLAB_ASPECT: f64 = 0.125;

// ---------------------------------------------------------------------------
// the cubed sphere
// ---------------------------------------------------------------------------

/// A [`Tiled`] recipe that is the whole ball rather than a square on it.
pub const WHOLE_BALL: u8 = 6;

/// The three axes of one face of the cube, as `(right, up, out)`.
///
/// Six faces, in the order `+x -x +y -y +z -z`. The handedness is consistent
/// across all six — `right × up = out` — so a patch's local frame is a rotation
/// of its parent's and never a reflection, which is what lets `Tree::axes_from`
/// compose them.
pub fn face_axes(face: u8) -> (Vec3, Vec3, Vec3) {
    match face % 6 {
        0 => (v3(0.0, 1.0, 0.0), v3(0.0, 0.0, 1.0), v3(1.0, 0.0, 0.0)),
        1 => (v3(0.0, 0.0, 1.0), v3(0.0, 1.0, 0.0), v3(-1.0, 0.0, 0.0)),
        2 => (v3(0.0, 0.0, 1.0), v3(1.0, 0.0, 0.0), v3(0.0, 1.0, 0.0)),
        3 => (v3(1.0, 0.0, 0.0), v3(0.0, 0.0, 1.0), v3(0.0, -1.0, 0.0)),
        4 => (v3(1.0, 0.0, 0.0), v3(0.0, 1.0, 0.0), v3(0.0, 0.0, 1.0)),
        _ => (v3(0.0, 1.0, 0.0), v3(1.0, 0.0, 0.0), v3(0.0, 0.0, -1.0)),
    }
}

/// A point on the unit sphere from a point on one face of the unit cube.
///
/// `docs/PLAY.md` D6 chose the cubed sphere over HEALPix and geodesic
/// subdivision for one reason that outranks their better area properties: **the
/// child relation is a clean square split**, which is what `PathKey`'s
/// child-index derivation consumes, so the surface tree *is* the scale tree
/// with no adapter. A pentagon defect or a nested-ring indexing scheme is not.
///
/// `(a, b)` run over `[-1, 1]` on the face. The plain normalisation is used
/// rather than one of the equal-area or tangent-warped variants: the angular
/// distortion is bounded at 1.3 and is corrected at the patch level, which is
/// exactly what D6 says to do with it, and a warp would put a transcendental
/// function inside the address arithmetic that has to regenerate bit-for-bit.
pub fn cube_to_sphere(face: u8, a: f64, b: f64) -> Vec3 {
    let (right, up, out) = face_axes(face);
    (out + right.scale(a) + up.scale(b)).unit()
}

/// How deep ground goes, as a fraction of how wide the piece of it is.
///
/// One rule at every level, which is what makes the ladder a ladder: a patch
/// ten thousand kilometres across is the outer eighth of a planet, and a patch
/// a metre across is the top twelve centimetres of soil. The number is
/// [`SLAB_ASPECT`], shared with a flat patch, because it is the same statement.
impl Tiled {
    /// Is this the whole ball, tiled into its six faces?
    pub fn is_ball(&self) -> bool {
        self.face == WHOLE_BALL
    }

    /// Is this a square on a sphere rather than a flat slab?
    pub fn on_sphere(&self) -> bool {
        self.sphere > 0.0
    }

    /// How many cells across this patch divides into.
    pub fn cells(&self) -> usize {
        (self.cells as usize).max(2)
    }

    /// The half-width of this patch on its face, in cube coordinates.
    ///
    /// Level zero is the whole face; every level halves it. `u` and `v` are the
    /// cell's index on the face at that level, so the address is exact integer
    /// arithmetic all the way down and a patch twenty-four levels deep
    /// regenerates bit-for-bit.
    pub fn face_span(&self) -> (f64, f64, f64) {
        let n = (self.cells() as f64).powi(self.level as i32);
        let half = 1.0 / n;
        let a = -1.0 + (2.0 * self.u as f64 + 1.0) * half;
        let b = -1.0 + (2.0 * self.v as f64 + 1.0) * half;
        (a, b, half)
    }

    /// The direction of this patch's centre from the planet's own centre.
    pub fn centre_direction(&self) -> Vec3 {
        if self.is_ball() {
            return v3(0.0, 0.0, 1.0);
        }
        let (a, b, _) = self.face_span();
        cube_to_sphere(self.face, a, b)
    }

    /// Which way is up here, as a rotation of the *planet's* axes.
    ///
    /// A patch is oriented by construction — that is what a patch is, a square
    /// of surface with a local up — and this is the rotation that says so. It
    /// is what `Tree::axes_from` composes and what makes gravity arrive along a
    /// standing thing's own `-z` at any latitude rather than along the planet's
    /// `-x`.
    pub fn frame(&self) -> crate::math::Quat {
        let up = self.centre_direction();
        let z = v3(0.0, 0.0, 1.0);
        let axis = z.cross(up);
        let s = axis.norm();
        if s < 1e-12 {
            return if up.z >= 0.0 {
                crate::math::Quat::IDENTITY
            } else {
                crate::math::Quat::from_axis_angle(v3(1.0, 0.0, 0.0), std::f64::consts::PI)
            };
        }
        crate::math::Quat::from_axis_angle(axis.scale(1.0 / s), s.atan2(z.dot(up)))
    }

    /// The recipe for one of this patch's cells.
    ///
    /// **The child relation, and the whole of what makes a surface a tree.**
    /// A cell of the ball is a face; a cell of a face is a quarter of it, or an
    /// `n`-th, and so on down. Nothing is stored per patch except the address,
    /// so a planet nobody has visited costs its `Matter` and nothing else.
    ///
    /// The last cell of a patch on a sphere is the *substrate* — what is under
    /// the surface rather than part of it — and it has no recipe, because it is
    /// not a piece of surface. See [`Tiled::render`].
    pub fn child(&self, cell: usize) -> Option<Tiled> {
        if !self.on_sphere() {
            return None;
        }
        if self.is_ball() {
            if cell >= 6 {
                return None;
            }
            let side = self.sphere * FACE_SIDE;
            return Some(Tiled {
                density: self.density,
                sphere: self.sphere,
                face: cell as u8,
                level: 0,
                u: 0,
                v: 0,
                cells: self.cells,
                side,
                depth: side * SLAB_ASPECT,
                relief: self.relief,
            });
        }
        let n = self.cells();
        if cell >= n * n {
            return None;
        }
        let (i, j) = (cell % n, cell / n);
        let side = self.side / n as f64;
        // The address has a bottom. Below it the arithmetic would wrap and a
        // patch would silently be a different patch, so the division stops
        // instead — which is honest: nothing can be addressed finer than the
        // address goes.
        let (u, v) = (
            self.u.checked_mul(n as u64)?.checked_add(i as u64)?,
            self.v.checked_mul(n as u64)?.checked_add(j as u64)?,
        );
        Some(Tiled {
            density: self.density,
            sphere: self.sphere,
            face: self.face,
            level: self.level.checked_add(1)?,
            u,
            v,
            cells: self.cells,
            side,
            depth: side * SLAB_ASPECT,
            relief: self.relief,
        })
    }

    /// Where a cell's centre of mass sits, relative to this patch's own.
    ///
    /// **In the planet's axes, not the patch's**, and that is deliberate: a
    /// node's `offset` is a position in the frame its parent's offset is
    /// expressed in, and nothing in the tree rotates an offset on the way up.
    /// Rotating these into the patch's own axes put a cell at 1.0e7 m from the
    /// centre of a planet 6.4e6 m across, because the offset and the one it was
    /// added to were in two different frames.
    ///
    /// What *is* relative is the cell's own facing, which is what
    /// `Tree::axes_from` composes. The two conventions are different and each
    /// is consistent: a position is stated once in one frame, and a facing is
    /// stated against the thing it is inside.
    fn cell_offset(&self, cell: usize) -> Option<(Vec3, f64)> {
        let child = self.child(cell)?;
        let here = self.centre_of_mass_from_planet();
        let there = child.centre_of_mass_from_planet();
        Some((there - here, child.side))
    }

    /// Where this patch's own centre of mass sits, from the planet's centre.
    ///
    /// **Integrated over the patch, not taken as its middle.** A patch is a
    /// curved wedge, and the centre of a curved thing is not on its surface's
    /// midpoint: a whole face of the cube has a centroid well inside the
    /// direction of its centre, and even a small patch is a little inward of
    /// it.
    ///
    /// Getting this wrong is not cosmetic. The sampler recentres a node's
    /// bodies on their own centre of mass, so if a recipe states its cells
    /// about a point that is *not* that centre, every level of a descent
    /// shifts by the difference — and they accumulate. Measured before this: a
    /// descent seven levels deep left the patch's centre 93 km from the
    /// observer it was descending towards, and stopped there.
    ///
    /// The integral separates: the radial factor is `int r^3 dr / int r^2 dr`
    /// over the patch's own depth, and the angular factor is the mean of the
    /// cube map's direction over the cell weighted by its Jacobian,
    /// `(1 + a^2 + b^2)^(-3/2)`. Eight points a side is far finer than the
    /// answer needs — the integrand has no structure — and it is deterministic,
    /// which is what regenerating bit-for-bit requires.
    pub fn centre_of_mass_from_planet(&self) -> Vec3 {
        if self.is_ball() {
            return Vec3::ZERO;
        }
        self.region_com(self.sphere - self.depth, self.sphere)
    }

    /// The centre of mass of the shell of this patch between two radii.
    /// The solid angle this patch covers, sr — the cube map's Jacobian,
    /// `(1 + a^2 + b^2)^(-3/2)`, over its span, by the same quadrature as
    /// `region_com`, so that the two agree about where its mass is.
    pub fn solid_angle(&self) -> f64 {
        let (ca, cb, half) = self.face_span();
        const Q: usize = 8;
        let cell = (2.0 * half / Q as f64).powi(2);
        let mut total = 0.0;
        for i in 0..Q {
            for j in 0..Q {
                let a = ca - half + 2.0 * half * (i as f64 + 0.5) / Q as f64;
                let b = cb - half + 2.0 * half * (j as f64 + 0.5) / Q as f64;
                total += (1.0 + a * a + b * b).powf(-1.5) * cell;
            }
        }
        total
    }

    pub fn region_com(&self, r0: f64, r1: f64) -> Vec3 {
        let (r0, r1) = (r0.min(r1).max(0.0), r0.max(r1).max(0.0));
        let radial = if r1 > r0 {
            0.75 * (r1.powi(4) - r0.powi(4)) / (r1.powi(3) - r0.powi(3)).max(1e-300)
        } else {
            r1
        };
        let (ca, cb, half) = self.face_span();
        const Q: usize = 8;
        let (right, up, out) = face_axes(self.face);
        let mut acc = Vec3::ZERO;
        let mut weight = 0.0;
        for i in 0..Q {
            for j in 0..Q {
                let a = ca - half + 2.0 * half * (i as f64 + 0.5) / Q as f64;
                let b = cb - half + 2.0 * half * (j as f64 + 0.5) / Q as f64;
                let w = (1.0 + a * a + b * b).powf(-1.5);
                acc = acc + (out + right.scale(a) + up.scale(b)).unit().scale(w);
                weight += w;
            }
        }
        if weight <= 0.0 {
            return self.centre_direction().scale(radial);
        }
        acc.scale(radial / weight)
    }
}

/// Side of one face of a cubed sphere, as a fraction of the sphere's radius.
///
/// A cube face covers a sixth of the sphere's area, `4 pi R^2 / 6`, so the
/// square with that area has a side of `R sqrt(2 pi / 3)`. Used to turn a
/// planet's radius into the side of the six patches that cover it, so a face
/// and the surface it stands for have the same area rather than the same
/// angular width.
pub const FACE_SIDE: f64 = 1.447_202_699_454_372_5;

impl Tiled {
    /// A piece of a planet's surface: its cells, and what is under them.
    ///
    /// **The cells are the children.** `docs/PLAY.md` D6: "a patch refines into
    /// sub-patches as an observer descends and coarsens behind them", and the
    /// bodies a patch holds *are* those sub-patches, so the surface tree is the
    /// scale tree with no adapter between them. The ball's cells are its six
    /// faces; a face's cells are its `n²` squares; and so on down to a square
    /// metre, which on Earth is twenty-four levels.
    ///
    /// # Why there is one more body than there are cells
    ///
    /// Each cell is only as deep as it is wide — [`SLAB_ASPECT`] of its own
    /// side, the same rule at every level — so the cells of a patch account for
    /// `1/n` of its volume and not all of it. The rest is *under* them, and it
    /// is one body: the substrate. Without it the mass would not conserve, and
    /// with it a descent sheds depth as it sheds width, which is what makes the
    /// thing you finally stand on a shallow patch of ground rather than a
    /// column reaching to the centre of the planet.
    fn render_on_sphere(&self, field: &[crate::erode::Deviation]) -> Skeleton {
        let n = self.cells();
        let count = if self.is_ball() { 6 } else { n * n };
        let mut sk = Skeleton::with_capacity(count + 1);
        // Volume of this patch, and of the cells that cover it.
        let mine = self.volume();
        let mut covered = 0.0;
        let here = self.centre_of_mass_from_planet();
        let mut cells: Vec<(Vec3, f64, f64, Quat)> = Vec::with_capacity(count);
        for c in 0..count {
            let Some(child) = self.child(c) else { continue };
            let Some((at, side)) = self.cell_offset(c) else { continue };
            let v = child.volume();
            covered += v;
            // The cell's own frame, expressed in this patch's: out of the
            // cell's axes by its frame, then into the patch's by the inverse of
            // the patch's. `a.then(b)` applies `b` first (`Quat::then`).
            let turn = self.frame().conjugate().then(child.frame());
            cells.push((at, side, v, turn));
        }
        let n_f = self.cells() as f64;
        for (i, (at, side, v, turn)) in cells.iter().enumerate() {
            let half = 0.5 * side;
            let depth = 0.5 * side * SLAB_ASPECT;
            // A cell is a slab, turned the way its own surface faces. Its
            // plan is the tiling's and may not move; its depth is what gives.
            let up = turn.rotate(v3(0.0, 0.0, 1.0));
            // What has happened here that the rule does not describe, at this
            // cell's own place on the patch. A patch on a sphere is divided the
            // same way a flat one is, so the coordinates mean the same thing.
            let mut lift = 0.0;
            if !field.is_empty() && !self.is_ball() {
                let u = -1.0 + 2.0 * ((i % self.cells()) as f64 + 0.5) / n_f;
                let w = -1.0 + 2.0 * ((i / self.cells()) as f64 + 0.5) / n_f;
                for d in field {
                    lift += d.height_at(u, w, 0.5 * self.side);
                }
            }
            let at = &(*at + up.scale(lift));
            sk.push_box(
                *at - up.scale(depth),
                *at + up.scale(depth),
                v3(half, half, depth),
                FREE_Z,
                *v,
                half,
                NO_SUPPORT,
                i as u32,
            );
            *sk.orientation.last_mut().unwrap() = *turn;
        }
        // What is under them. On a ball, one body — its whole interior, at its
        // centre, which every face rests on straight down. On a patch, **a
        // column under each cell**, reaching from the cell's underside to the
        // patch's floor at the column's own centre of mass, so that each cell
        // rests on what is directly beneath it. As one body at the patch's
        // centre instead, a cell at a face's corner rested on it along a line
        // mostly across the face rather than down, nothing held it up or down
        // but its compressed neighbours, and the face buckled like a shell at
        // its corners, 1.2e-6 s^-2.
        let under = (mine - covered).max(0.0);
        if under > 0.0 && !self.is_ball() && covered > 0.0 {
            let skin = self.side / self.cells() as f64 * SLAB_ASPECT;
            let depth = 0.5 * (self.depth - skin).max(0.0);
            // Each column's share by the solid angle it covers: a cell of a
            // cube-mapped sphere near a face's corner covers less of it than
            // one in the middle. Shared by the cells' flat areas instead, the
            // columns' masses disagreed with where their centres are, every
            // level of a descent was recentred by the difference, and an
            // observer nine levels down read 9.1954 m/s^2 against 9.7055.
            let angles: Vec<f64> = (0..cells.len()).map(|i| self.child(i).map(|c| c.solid_angle()).unwrap_or(0.0)).collect();
            let all: f64 = angles.iter().sum();
            // `under` is the patch's shell less its cells' shells; each
            // column takes its share by the solid angle over it, which is the
            // same shell of rock between the cells' floor and the patch's.
            for (i, (_, side, _, turn)) in cells.iter().enumerate() {
                let Some(child) = self.child(i) else { continue };
                if !(all > 0.0) {
                    continue;
                }
                let at = child.region_com(self.sphere - self.depth, self.sphere - skin) - here;
                let up = turn.rotate(v3(0.0, 0.0, 1.0));
                sk.push_box(
                    at - up.scale(depth),
                    at + up.scale(depth),
                    v3(0.5 * side, 0.5 * side, depth),
                    FREE_Z,
                    under * angles[i] / all,
                    0.5 * side,
                    NO_SUPPORT,
                    (count + i) as u32,
                );
                *sk.orientation.last_mut().unwrap() = *turn;
            }
            return sk;
        }
        if under > 0.0 {
            let skin = if self.is_ball() {
                self.sphere * FACE_SIDE * SLAB_ASPECT
            } else {
                self.side / self.cells() as f64 * SLAB_ASPECT
            };
            let at = if self.is_ball() {
                // A ball's substrate is its whole interior, centred on it.
                Vec3::ZERO
            } else {
                self.region_com(self.sphere - self.depth, self.sphere - skin) - here
            };
            let r = (0.75 * under / std::f64::consts::PI).cbrt();
            sk.push_segment(
                at - v3(0.0, 0.0, r),
                at + v3(0.0, 0.0, r),
                under,
                r,
                NO_SUPPORT,
                count as u32,
            );
        }
        sk
    }

    /// The volume this patch stands for, m^3: the shell it cuts from its
    /// sphere, `Omega (R^3 - (R - d)^3) / 3` over the solid angle it covers.
    ///
    /// **The volume of the shell, not of a flat slab of its side.** A cell
    /// of a cube-mapped sphere near a face's corner covers a fifth of the
    /// solid angle of one in its middle, and the same side. Weighed flat, the
    /// cells' masses disagreed with where their centres are, and — once the
    /// columns under them were weighed by their shells — a corner cell stood on
    /// a column a fifth of the weight it expected, and a face under air
    /// buckled at 9e-7 s^-2.
    pub fn volume(&self) -> f64 {
        if self.is_ball() {
            return 4.0 / 3.0 * std::f64::consts::PI * self.sphere.powi(3);
        }
        let inner = (self.sphere - self.depth).max(0.0);
        self.solid_angle() * (self.sphere.powi(3) - inner.powi(3)) / 3.0
    }
}

/// How two pieces of ground touch: side by side across a face they share, one
/// resting on the other, or corner to corner across a diagonal of the grid.
///
/// **The diagonal is what carries shear.** Springs along a grid's rows and
/// columns alone let a square of it lean into a rhombus for nothing, and a
/// ground carrying its own weight is in compression, which leans it further:
/// measured, a face of a turning Earth holding its weight as a stress between
/// neighbours slipped over its ground from 4 to 21 m/s in two hours. A
/// diagonal spring of `G t` — the shear modulus over the slab's depth — gives
/// the square lattice the rock's own resistance to shear.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Touch {
    Beside,
    On,
    Across,
}

impl Tiled {
    /// The first of this ground's pieces that rest on what is outside it — the
    /// columns at a patch's floor, which the rest of the planet holds up. A
    /// ball rests on nothing, and every piece of it is its floor.
    pub fn floor(&self) -> usize {
        if !self.on_sphere() || self.is_ball() {
            0
        } else {
            self.cells() * self.cells()
        }
    }

    /// Which pieces of this ground touch which, by the index of the body each
    /// is drawn as (`Tiled::render_on_sphere`), and how.
    ///
    /// A patch's cells touch the cells beside them in its grid and rest on
    /// the substrate under them, which is the body after the last cell. A
    /// ball's six faces each touch the four that are not opposite them and
    /// rest on its interior. Only ground on a sphere is laid out this way; a
    /// flat patch is drawn as columns and has no joints here.
    pub fn joints(&self) -> Vec<(usize, usize, Touch)> {
        let mut out = Vec::new();
        if !self.on_sphere() {
            return out;
        }
        if self.is_ball() {
            for a in 0..6usize {
                for b in (a + 1)..6 {
                    if a / 2 != b / 2 {
                        out.push((a, b, Touch::Beside));
                    }
                }
                out.push((a, 6, Touch::On));
            }
            return out;
        }
        // Two layers of the same lattice: the cells, and the columns under
        // them (`render_on_sphere`), each cell resting on its own column.
        let n = self.cells();
        let count = n * n;
        for layer in [0, count] {
            for c in 0..count {
                let (i, j) = (c % n, c / n);
                if i + 1 < n {
                    out.push((layer + c, layer + c + 1, Touch::Beside));
                }
                if j + 1 < n {
                    out.push((layer + c, layer + c + n, Touch::Beside));
                }
                if i + 1 < n && j + 1 < n {
                    out.push((layer + c, layer + c + n + 1, Touch::Across));
                }
                if i > 0 && j + 1 < n {
                    out.push((layer + c, layer + c + n - 1, Touch::Across));
                }
            }
        }
        // Each cell on its own column, and across to the columns under its
        // neighbours: the diagonals that carry shear between the two layers.
        // Without them a layer of cells resting on compressed columns is a
        // sheet of inverted pendulums, and slid over the columns under it at
        // 2e-5 s^-2.
        for c in 0..count {
            out.push((c, count + c, Touch::On));
            let (i, j) = (c % n, c / n);
            if i + 1 < n {
                out.push((c, count + c + 1, Touch::Across));
                out.push((c + 1, count + c, Touch::Across));
            }
            if j + 1 < n {
                out.push((c, count + c + n, Touch::Across));
                out.push((c + n, count + c, Touch::Across));
            }
        }
        out
    }
}

impl Tiled {
    /// Which cell of this patch a direction from the planet's centre falls in.
    ///
    /// **The inverse of the parameterisation, and the reason there is one.**
    /// Finding the cell an observer is over by searching for the nearest cell
    /// body is a plausible thing to write and it does not work: the cells of a
    /// face are spread over a curved square nine thousand kilometres across,
    /// and the nearest *centre* to a point above the middle of it is not the
    /// cell below that point. Measured, a descent that searched drifted to a
    /// corner cell within three levels and stopped 93 km from the observer.
    ///
    /// Inverting the cube map is exact and costs a divide: the face is the axis
    /// the direction leans on hardest, and `(a, b)` are the other two
    /// components over that one. `None` when the direction is not over this
    /// patch at all, which is how a walk knows it has left.
    pub fn cell_of_direction(&self, dir: Vec3) -> Option<usize> {
        let d = dir.unit();
        if !d.is_finite() {
            return None;
        }
        if self.is_ball() {
            // The ball's cells are the six faces, and a direction is over
            // exactly one of them.
            let (x, y, z) = (d.x, d.y, d.z);
            let face = if x.abs() >= y.abs() && x.abs() >= z.abs() {
                if x >= 0.0 { 0 } else { 1 }
            } else if y.abs() >= z.abs() {
                if y >= 0.0 { 2 } else { 3 }
            } else if z >= 0.0 {
                4
            } else {
                5
            };
            return Some(face);
        }
        let (right, up, out) = face_axes(self.face);
        let o = d.dot(out);
        if o <= 1e-12 {
            // On the far side of the cube from this face.
            return None;
        }
        let (a, b) = (d.dot(right) / o, d.dot(up) / o);
        let (ca, cb, half) = self.face_span();
        let n = self.cells() as f64;
        // Where in this patch, as a fraction of its own span.
        let fa = (a - (ca - half)) / (2.0 * half);
        let fb = (b - (cb - half)) / (2.0 * half);
        if !(0.0..1.0).contains(&fa) || !(0.0..1.0).contains(&fb) {
            return None;
        }
        let i = (fa * n).floor().clamp(0.0, n - 1.0) as usize;
        let j = (fb * n).floor().clamp(0.0, n - 1.0) as usize;
        Some(j * self.cells() + i)
    }

    /// Where a direction lands on this patch, in the patch's own normalised
    /// coordinates: -1 to 1 across it, in the same frame `surface_height` and a
    /// [`crate::erode::Deviation`] use.
    ///
    /// `None` when the direction is not over this patch at all, which is the
    /// honest answer and is what makes a deviation belong to exactly one node.
    pub fn local_of_direction(&self, dir: Vec3) -> Option<(f64, f64)> {
        let d = dir.unit();
        if !d.is_finite() || self.is_ball() {
            return None;
        }
        let (right, up, out) = face_axes(self.face);
        let o = d.dot(out);
        if o <= 1e-12 {
            return None;
        }
        let (a, b) = (d.dot(right) / o, d.dot(up) / o);
        let (ca, cb, half) = self.face_span();
        let x = (a - ca) / half;
        let y = (b - cb) / half;
        if !(-1.0..=1.0).contains(&x) || !(-1.0..=1.0).contains(&y) {
            return None;
        }
        Some((x, y))
    }

    /// The direction of the point on this patch's surface nearest a direction.
    ///
    /// Used to ask where a walk has got to when it has left this patch: the
    /// answer is the edge it went out through.
    pub fn clamped_direction(&self, dir: Vec3) -> Vec3 {
        if self.is_ball() {
            return dir.unit();
        }
        let (right, up, out) = face_axes(self.face);
        let d = dir.unit();
        let o = d.dot(out);
        if o <= 1e-12 {
            return self.centre_direction();
        }
        let (ca, cb, half) = self.face_span();
        let a = (d.dot(right) / o).clamp(ca - half, ca + half);
        let b = (d.dot(up) / o).clamp(cb - half, cb + half);
        cube_to_sphere(self.face, a, b)
    }
}

// ---------------------------------------------------------------------------
// granular
// ---------------------------------------------------------------------------

/// A packing of grains: what a melt that froze is laid out as.
///
/// **Nothing selects this habit; it is derived.** A node whose own recorded
/// past shows it crossing its melting point while cooling *froze*, and what a
/// freezing melt produces is grains — at a size the competition between how
/// fast nuclei appear and how fast they grow decides, which `grain_scale`
/// already solves from the cooling rate. So the layout, and the one number in
/// it, are both consequences of the node's history rather than of anything
/// anybody stated.
///
/// The grain is the *physics*; the cells this draws are the *resolution*. A
/// cubic metre of granite has 10^8 grains and nobody is going to draw them, so
/// a cell stands for as many of them as the budget requires — which is the
/// third axiom again, a rule derived once and run at whatever detail is asked
/// for.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Granular {
    /// The grain the freezing left, metres. Derived from the node's own
    /// cooling rate, and the thing a crack in this rock has to get around.
    pub grain: f64,
    /// Bulk density, kg/m^3. See [`Coursed::density`].
    pub density: f64,
    /// Side of the cube this fills, metres.
    pub side: f64,
}

impl Granular {
    pub fn extent(&self) -> f64 {
        0.5 * self.side * 3.0f64.sqrt()
    }

    pub fn height(&self) -> f64 {
        self.side
    }

    /// How many grains across this piece is.
    pub fn grains_across(&self) -> f64 {
        (self.side / self.grain.max(1e-30)).max(1.0)
    }

    /// A cubic packing of cells filling the piece.
    pub fn render(&self, budget: usize) -> Skeleton {
        // As many cells across as the budget and the grain both allow: never
        // finer than the grains themselves, because below that there is nothing
        // there to draw.
        let by_budget = (budget as f64).cbrt().floor().max(1.0);
        let n = by_budget.min(self.grains_across()).max(1.0) as usize;
        let step = self.side / n as f64;
        let half = 0.5 * step;
        let origin = -0.5 * self.side + half;
        let mut sk = Skeleton::with_capacity(n * n * n);
        let mut site = 0u32;
        let mass = 1.0;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let at = v3(
                        origin + step * i as f64,
                        origin + step * j as f64,
                        origin + step * k as f64,
                    );
                    // A packing fills what it fills; none of its axes is free,
                    // because a cell that grew would overlap the one next to it
                    // and one that shrank would leave a hole.
                    sk.push_box(at, at, v3(half, half, half), 0, mass, half, NO_SUPPORT, site);
                    site += 1;
                }
            }
        }
        sk
    }
}
