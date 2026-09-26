//! The scale tree: the object that is 10^68 particles without storing them.
//!
//! A `Node` is a region of space at a tier, holding a `Matter`. It has
//! two optional finer representations:
//!
//! * **materialised bodies** — a `Vec<Body>` produced by `sample`. Cheap to
//!   make, cheap to throw away, regenerable bit-for-bit.
//! * **promoted children** — full `Node`s standing in for individual bodies,
//!   created only for the handful of bodies someone is actually looking at.
//!
//! The pair gives two independent LOD axes. Materialising is how you get from
//! "a molecular cloud" to "a million gas parcels"; promoting is how you get
//! from "one of those parcels" to "a protostar with its own internal
//! structure". A path from the galaxy to a nucleus promotes about seven times
//! and materialises about seven times, so the live node count along any single
//! zoom is in the thousands, not the billions.
//!
//! # Discarding detail is the normal case
//!
//! Every frame, most of the tree is *deleted*. That is not a cache eviction
//! policy bolted on the side; it is the whole design. The invariant that makes
//! it safe is stated in `sampler.rs` and enforced in `tests/consistency.rs`:
//! summarising after sampling returns the same conserved tuple. Detail that
//! has been *touched* — measured, or hit by something — is different, and is
//! pinned (see `Node::pinned` and `observe::Ledger`).

use crate::coords::{Motion, Bounded};
use crate::ids::{NodeIdx, PathKey};
use crate::math::Vec3;
use crate::sampler::{SampleReport, SampleSpec};
use crate::state::{summarise, Matter, Body};
use crate::units::Tier;
use std::collections::HashMap;

/// How closely summarising must agree with the stored matter before the
/// engine treats the two as the same state. Set well above the round-off floor
/// (~10^-16) and far below anything physically detectable.
pub const IDEMPOTENT_TOLERANCE: f64 = 1e-12;

/// Why a node currently holds fine detail. Drives the eviction order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Residency {
    /// Nothing is looking at it; free to discard at any time.
    Speculative,
    /// Inside an observer's interest volume this frame.
    Observed,
    /// Inside the causal past of something observed — must be resolved even
    /// though nobody is looking at it directly, because its influence will
    /// arrive at an observer within the horizon.
    Causal,
    /// Touched by a user interaction or a recorded measurement. Its detail is
    /// no longer derivable and must be persisted, not regenerated.
    Pinned,
}

impl Residency {
    /// Higher survives eviction.
    pub fn rank(self) -> u8 {
        match self {
            Residency::Speculative => 0,
            Residency::Causal => 1,
            Residency::Observed => 2,
            Residency::Pinned => 3,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Node {
    pub key: PathKey,
    pub parent: NodeIdx,
    /// Index of this node in its parent's materialised body list.
    pub slot: u32,
    pub depth: u32,
    pub tier: Tier,

    /// The node's matter. Always present — this is what a node *is*.
    pub matter: Matter,
    /// Position and velocity relative to the parent node's frame.
    pub motion: Motion,

    /// Fine detail, if currently materialised.
    pub bodies: Vec<Body>,
    /// Self-potential the materialisation was built against. Summarising must
    /// use this same number (see `SampleReport::potential`).
    pub potential: f64,
    /// The gravitational field this node sits in, in its own frame.
    ///
    /// Derived by [`Tree::gravity_at`] from what the node is inside, and stored
    /// for the same reason `potential` is: a structure is *proportioned*
    /// against its own weight, so its member radii depend on this number, and a
    /// client regenerating the structure has to use the one the geometry was
    /// built with rather than whatever it would work out for itself. It travels
    /// in the node's wire payload with `matter` and `spec` because it is an
    /// input to regeneration exactly as they are.
    ///
    /// Refreshed where a node is materialised, which is where it matters. A
    /// node that has not been sampled has no geometry for it to be wrong about.
    pub gravity: Vec3,
    /// Children promoted from `bodies`; `NodeIdx::NONE` where not promoted.
    /// Parallel to `bodies`, and empty when nothing is promoted.
    pub children: Vec<NodeIdx>,

    /// How this node is to be split when refined.
    pub spec: SampleSpec,
    /// Bumped whenever a recorded interaction changes the node's contents.
    /// Detail regenerated at the same epoch is identical; a new epoch means the
    /// old detail is gone for good.
    pub epoch: u32,
    /// World instant this node's state is valid at.
    ///
    /// Every live node is brought to the world instant every frame, whether or
    /// not it was solved: position and orientation under constant velocity and
    /// spin are closed-form, so carrying a node forward costs one add and is
    /// exact. That is what makes a shared "now" affordable across thirty-eight
    /// orders of magnitude — the expensive thing is not knowing where a node
    /// is, it is working out what it is doing.
    pub time: f64,
    /// World instant at which something last *happened* here: the node was
    /// resolved, hit, measured, or otherwise put into a state the sampler did
    /// not produce.
    ///
    /// Detail is kept for a mixing time after this and then released, because
    /// past a mixing time the stored sample is no longer *that* state, only *a*
    /// state of the same matter — which the sampler can draw for free.
    pub last_disturbed: f64,
    /// World instant at which this node's dynamics were last re-derived, as
    /// opposed to merely carried forward.
    ///
    /// `time - last_solved`, divided by the node's characteristic time, is its
    /// lateness, and lateness is what the scheduler ranks on.
    pub last_solved: f64,
    /// World instant this node's morphology was last advanced to. Separate
    /// from `last_solved` because growth runs on the matter and costs O(1),
    /// so it keeps up on nodes whose dynamics cannot.
    pub last_grown: f64,
    pub residency: Residency,
    /// Set when the node's detail has been altered away from what `sample`
    /// would produce, so it must be stored rather than regenerated.
    ///
    /// **This node's own detail**, and not its descendants'. See
    /// [`Node::contains_edit`], which is what the ancestry walk sets and which
    /// `docs/PLAY.md` D19 separates from this.
    pub pinned: bool,
    /// Set when something *below* this node has been changed away from what the
    /// sampler would produce.
    ///
    /// `docs/PLAY.md` D19. `Tree::pin` used to set `pinned` on the whole
    /// ancestry, on the reasoning that "a changed child means the parent's
    /// materialisation no longer matches what `sample` would produce" — which
    /// is true, and does not imply what it was being used to imply. `pin` is
    /// one-way and walks to the root, and the scheduler refuses to coarsen a
    /// pinned node at all, so **felling one tree made a whole planet
    /// permanently un-coarsenable**, and then its star, and then its galaxy.
    ///
    /// An ancestor that merely contains an edit can still collapse, because its
    /// recipe plus its descendants' own edits is all it needs. What it may not
    /// do is *forget*: a fresh draw from the ensemble would quietly mend the
    /// change, which is the one thing the flag is for.
    pub contains_edit: bool,
    /// Density the node's condensed matter has at rest, kg/m^3, or zero for a
    /// node whose matter the gas law describes or nobody has described.
    ///
    /// Derived by `World` from the mixture and its registry (`eos.rs`) and
    /// stored, for the reason `gravity` is: it is an **input to regeneration**.
    /// A liquid or a packed solid is drawn at the spacing its own density
    /// implies rather than as a Poisson cloud — a random draw of water put its
    /// parcels at 0.86 to 2.3 times rest density, which Tait turns into 10^11
    /// Pa — and a client regenerating the node holds no registry to derive it
    /// from.
    pub rest_density: f64,
    /// The largest net acceleration the node's last solve left on its loose
    /// contents, m/s^2, or infinity for "not yet measured".
    ///
    /// Only ever not zero for contents that are held in a field by something
    /// in their own node — water in a bucket, litter on the ground — since that
    /// is the only case where contents at rest can be out of balance: see
    /// `World::node_cadence`. Marked unknown where such a node is drawn or has
    /// something promoted into it, and measured by its next solve. Not
    /// persisted; a reloaded node is drawn again and measured again.
    pub unrest: f64,
    /// The ocean this node carries at its own scale, if it has one: its
    /// dynamic tide and the sea the wind has raised on it. See `ocean.rs`.
    ///
    /// State, not a drawing: a tide depends on its history, so it travels
    /// with the node and is persisted. Boxed, because almost no node has one.
    pub ocean: Option<Box<crate::ocean::Ocean>>,
    /// The world instant this node's *motion* has been carried to, s.
    ///
    /// **Not the same clock as `time`**, which is how far the node's contents
    /// have been solved. Where a node is can be carried to the world instant
    /// every frame, in closed form, and is — that is what keeps one instant
    /// everywhere. What it holds cannot: the owner's rule for Phase 5 is that
    /// no node is ever left unsolved, only scheduled, so a materialised node the
    /// frame did not reach keeps its contents at the time they were solved to
    /// and covers the whole of what it is owed when it runs. Coasting `time`
    /// along with the motion used to lose that span from the contents for good:
    /// measured, a moon advanced 300 km per 3.6-day solve instead of 313,000.
    /// See `Tree::carry` and `Tree::position_at`.
    pub carried: f64,
    /// The angular velocity of the frame this node is carried in between its
    /// parent's solves, rad/s, root-aligned: its parent's spin where the parent
    /// *holds* it, and zero where it does not.
    ///
    /// **Held things turn with what holds them** — the owner's decision for
    /// Phase 5. A straight line is the exact carry for a thing no force acts
    /// on, and wrong for one whose forces supply the centripetal of a turning
    /// body: measured, a parcel of air neutrally buoyant over a turning Earth
    /// that nobody solved went along its tangent to 31,800 km in a day. The
    /// turning carry is the account of that centripetal, so a held thing's
    /// solves apply only what its forces do beyond it.
    ///
    /// **What it does not carry is Coriolis.** Between solves a held thing
    /// keeps its velocity relative to the ground, so wind and current over a
    /// turning planet are not deflected by its turning. A first version counted
    /// the centripetal at each solve as well — the fluid's own acceleration in
    /// its push, and the carry's velocity change handed back — so that over any
    /// span the forces alone decided; that is exact only where every solve
    /// applies the true forces, and a patch of ground's solve applies no
    /// gravity to what it holds. Measured there: a parcel of air at its own
    /// ambient density went inward at 2.0 m/s in the first minute and 13.5 m/s
    /// after the second solve. Set each frame by `World::hold`; see
    /// `Tree::carry`.
    pub turning: Vec3,
    /// What this node's ground is, where it is ground: derived on its first
    /// solve after it is drawn. See `solvers::ground`.
    pub ground: Option<Box<crate::solvers::ground::Cache>>,
    /// A deliberate, unphysical multiplier on how fast this node's *interior*
    /// runs, and its whole subtree's with it.
    ///
    /// One for everything the engine builds. An administrator sets it to watch
    /// a century of growth in an afternoon, or a reaction go to completion
    /// while the world barely moves — testing and balancing, not physics. It
    /// multiplies into `dilation::TimeRate` alongside the two relativistic
    /// factors, so there is exactly one place that decides how fast anything
    /// evolves; and it is reported separately from them, so a log can always
    /// say which part was the universe and which part was somebody's thumb.
    ///
    /// It scales the interior only. Position and orientation still advance on
    /// coordinate time, because an object that *travelled* a hundred times
    /// faster is an object with a hundred times the velocity — which the engine
    /// already models — and because a node whose position ran ahead of the
    /// world clock would outrun the influences it had itself emitted.
    pub bubble: f64,
    pub alive: bool,
    /// Developmental state, when this node is a structure rather than a
    /// statistical population. Its presence switches materialisation from
    /// max-entropy sampling to program-driven generation, and makes the node's
    /// geometry, entropy and stored free energy the morphology's business
    /// rather than the sampler's.
    pub morphology: Option<crate::morph::Morphology>,
    /// Joints holding the materialised parts together. Present only while the
    /// node is materialised, and regenerated with the geometry.
    pub topology: Option<crate::topology::Topology>,
    /// Number of solver steps this node has taken. Part of the address for any
    /// per-step randomness (see `rng::Stream::split`).
    pub steps_taken: u64,
    /// The boundary this node presents, baked once and kept. `docs/PLAY.md`
    /// D13 and D18.
    ///
    /// **Derived, with the shortcut stored** — the third axiom read literally,
    /// applied to shape, which is the one place the engine never applied it.
    /// `collision_shape` rebuilt a proxy from the member list on every frame
    /// and kept nothing, which is a derivation with no shortcut. An undisturbed
    /// tree now computes its boundary once and reuses it for a thousand frames;
    /// a felled one recomputes, because its arrangement changed.
    ///
    /// Not persisted and not part of the wire format: it is regenerable from
    /// the node's own contents, which is the same reason `last_report` is not.
    pub surface: Option<crate::shape::Surface>,
    /// The `epoch` the surface above was baked at. A different one means the
    /// arrangement moved and the surface is stale.
    ///
    /// `epoch` is exactly the right trigger and not an approximation of one:
    /// it is "bumped whenever a recorded interaction changes the node's
    /// contents", which is the definition of when a boundary stops being the
    /// boundary.
    pub surface_epoch: u32,
    pub last_report: SampleReport,
}

/// What a move changed, so a caller holding path-keyed data can follow it.
///
/// A node's [`PathKey`] *is* its path, so moving it changes the key of the node
/// and of every descendant. Everything addressed by path key — pinned detail,
/// the ledger, the chemistry and environment tables — has to be carried across,
/// and this is the list to carry it by. Parents come before children, so
/// applying it in order never orphans anything.
#[derive(Debug, Clone)]
pub struct Rehomed {
    pub moved: NodeIdx,
    pub from: NodeIdx,
    pub to: NodeIdx,
    /// `(old, new)` for the moved node and every node beneath it.
    pub keys: Vec<(PathKey, PathKey)>,
}

impl Node {
    pub fn is_materialised(&self) -> bool {
        !self.bodies.is_empty()
    }

    /// Which of this node's bodies are parts of a structure, parallel to
    /// `bodies`. `None` when the node holds no ordered matter at all.
    ///
    /// The discriminator is the joint's radius, not `SampleReport`'s
    /// `structural_parts`: the report is a diagnostic of the last
    /// materialisation and is deliberately not persisted, while the topology
    /// is. A node reloaded from a file has to give the same answer as the one
    /// that wrote it, or `docs/PLAY.md` §3.3's dispatch changes across a save.
    ///
    /// Measured on a materialised tree: 510 joints carry a radius and 90 do
    /// not, against a reported `structural_parts` of 510 — the same split, from
    /// the half that survives a round trip.
    pub fn structural_mask(&self) -> Option<Vec<bool>> {
        let t = self.topology.as_ref()?;
        let mut any = false;
        let mask: Vec<bool> = (0..self.bodies.len())
            .map(|i| {
                let ordered = t.joints.get(i).map(|j| j.radius > 0.0).unwrap_or(false);
                any |= ordered;
                ordered
            })
            .collect();
        any.then_some(mask)
    }

    /// Re-derive the angular *velocity* from the angular *momentum* the node
    /// actually holds.
    ///
    /// A node carries two spin quantities and they are not the same thing:
    /// `matter.spin` is angular momentum in kg m^2/s and `motion.spin_rate` is
    /// angular velocity in rad/s, related by the moment of inertia. Until this
    /// existed, `spin_rate` was derived at node *construction* and nowhere
    /// else — two creation paths, and then never again — so a collision added
    /// to `matter.spin` through `apply_contact` and the node's `orientation`
    /// never heard about it.
    ///
    /// Measured before the fix, on `a_ball_loose_in_a_box...`: the box banks
    /// `matter.spin` of 7.4865 kg m^2/s over five off-centre strikes while
    /// remaining, as far as `motion` is concerned, perfectly still.
    ///
    /// Called from where the momentum, the mass or the radius can have moved
    /// — the two of those are `Matter::moment_of_inertia`'s inputs — which is
    /// after a solve, after a contact, and after a child's evolved state is
    /// folded back. Not from coasting: coasting asserts that nothing changed,
    /// and re-deriving there would be either a no-op or a lie about which.
    ///
    /// **For stating a node's turning, not for keeping it.** Since Phase 5 a
    /// node's `spin_rate` is an offset from its parent's turning (the owner's
    /// decision), so `L / I` is its turning only for a node whose parent does
    /// not turn — a root, which is what a scene states. What moves it after
    /// that is angular momentum arriving, by `Tree::turn_by`.
    pub fn sync_spin_rate(&mut self) {
        self.motion.spin_rate = self.matter.angular_velocity();
    }

    /// What this node presents to a contact, in its own frame.
    ///
    /// A node is the engine's rigid body — it has one velocity and one spin,
    /// and `apply_contact` has always put an impulse straight onto them. What
    /// it did not have was a *shape*: it presented `matter.radius`, so a
    /// timber-framed building and a boulder of the same size were the same
    /// marble to anything that hit them.
    ///
    /// The partition is §3.3's, and deliberately the same one the solver
    /// dispatch uses rather than a second opinion about what a node is:
    ///
    /// - **Ordered contents have a shape.** Each structural member is a capsule
    ///   from `base` to `tip` at the joint's cross-section radius — three
    ///   numbers `Topology` has carried all along while contact measured the
    ///   midpoint bead instead.
    /// - **Disordered contents do not.** A parcel of gas, a cloud of rubble or
    ///   a star cluster has no surface to speak of, and the node falls back to
    ///   the sphere of its own radius. That is not a special case for
    ///   unstructured matter; it is the honest shape of a thing whose contents
    ///   are a distribution.
    ///
    /// Returns pieces that are each convex. A structure is many of them — a box
    /// is six walls around a cavity, and one hull over the whole thing would
    /// enclose its own contents.
    ///
    /// **In the node's own frame**, which is where a shape belongs: a shape
    /// that had the node's position and orientation baked into it would be
    /// wrong on the next frame. A caller places it with [`crate::shape::Hull::placed`],
    /// giving the node's orientation and where it is — in that order, which is
    /// the order `Motion::body_to_parent` places a body-fixed point in. The
    /// contact path used to translate without rotating, and a spinning box's
    /// walls stayed in the axes they were built in.
    /// **Superseded by `World::surface_of_node`**, which is `docs/PLAY.md`
    /// D18's stored union of solid convex primitives. This is the per-frame
    /// rebuild D13 retires: it has no material per piece, no reconciliation
    /// against the node's own pools, and no cache. It survives because the
    /// renderer and the fragment path still read a bare hull list and because a
    /// `Tree` has no access to the substance registry a measured material needs.
    pub fn collision_shape(&self) -> Vec<crate::shape::Hull> {
        let Some(mask) = self.structural_mask() else {
            return vec![crate::shape::Hull::sphere(Vec3::ZERO, self.matter.radius)];
        };
        let Some(t) = self.topology.as_ref() else {
            return vec![crate::shape::Hull::sphere(Vec3::ZERO, self.matter.radius)];
        };
        let mut out = Vec::new();
        for (i, ordered) in mask.iter().enumerate() {
            if !ordered {
                continue;
            }
            let radius = t.joints.get(i).map(|j| j.radius).unwrap_or(0.0);
            let base = t.base.get(i).copied().unwrap_or(Vec3::ZERO);
            let tip = t.tip.get(i).copied().unwrap_or(Vec3::ZERO);
            // A member with no length is a block rather than a beam — coursed
            // masonry is the case — and its own body's sphere is the right
            // shape for it.
            if (tip - base).norm2() > 0.0 {
                out.push(crate::shape::Hull::capsule(base, tip, radius));
            } else if let Some(b) = self.bodies.get(i) {
                out.push(crate::shape::Hull::sphere(b.pos, b.radius.max(radius)));
            }
        }
        if out.is_empty() {
            out.push(crate::shape::Hull::sphere(Vec3::ZERO, self.matter.radius));
        }
        out
    }

    pub fn child_of(&self, slot: usize) -> NodeIdx {
        self.children.get(slot).copied().unwrap_or(NodeIdx::NONE)
    }

    /// Bytes of fine detail held by this node alone.
    pub fn detail_bytes(&self) -> usize {
        self.bodies.len() * std::mem::size_of::<Body>()
            + self.children.len() * std::mem::size_of::<NodeIdx>()
    }
}

/// Arena of nodes plus the persistent store for pinned detail.
pub struct Tree {
    pub nodes: Vec<Node>,
    free: Vec<u32>,
    pub root: NodeIdx,
    pub world_seed: u64,
    /// Detail that can no longer be regenerated because something interacted
    /// with it. This is the only part of the world that costs permanent memory,
    /// and it grows only in proportion to what users actually touch.
    pub persisted: HashMap<PathKey, Vec<Body>>,
    pub stats: TreeStats,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct TreeStats {
    pub materialisations: u64,
    pub coarsenings: u64,
    /// Layouts the engine derived for itself, from a node's own state and its
    /// own past. See `World::derive_layout`.
    pub layouts_derived: u64,
    /// Deviations that decayed to nothing and were dropped. `docs/PLAY.md`
    /// §5.1: forgetting is the deviation reaching zero, and this counts it.
    pub deviations_forgotten: u64,
    /// Nodes moved to a different parent. Each one rekeyed a whole subtree.
    pub reparents: u64,
    /// Coarsenings where the fine detail turned out to say nothing new, so the
    /// coarse state was left exactly as it was.
    pub idempotent_coarsenings: u64,
    /// Nodes carrying a developmental state.
    pub structures: u64,
    /// Growth and construction steps advanced on matter alone, without ever
    /// materialising the structures they describe.
    pub growth_steps: u64,
    /// Structures loaded to failure.
    pub damage_events: u64,
    /// Parts attached into a composite's recipe. `docs/PLAY.md` D15's forward
    /// direction; a node went away and a part appeared.
    pub joins: u64,
    /// Parts that came away and became nodes of their own. The same transform
    /// run backwards, and the count that says how much of a world is loose.
    pub detachments: u64,
    /// Energy that has crossed a node boundary inwards to drive growth, J.
    /// The world's energy is not conserved against this — it is *balanced*
    /// against it, which is what `tests/growth.rs` asserts.
    pub external_energy_absorbed: f64,
    pub bodies_created: u64,
    pub bodies_discarded: u64,
    pub promotions: u64,
    /// Nodes whose tier was re-derived after their size changed, and moved.
    /// Non-zero means something grew or was built across a regime boundary.
    pub retiers: u64,
    pub persisted_bodies: u64,
    /// Coarsenings where blending the descriptions ran out of slots.
    ///
    /// `docs/PLAY.md` §5A.4's granularity signal: a node that cannot describe
    /// what it holds is a node that should have subdivided. Nothing subdivides
    /// on it yet — that is §5A.4's own work and Phase 3's — so this is the
    /// measurement standing where the rule will go, on §3.7's precedent that a
    /// node crossed by its own ensemble says so rather than doing it quietly.
    pub over_described: u64,
    /// Worst single such loss, as a fraction of the node's mass.
    pub worst_description_lost: f64,
    /// Worst conservation error seen across every scale transition so far.
    pub worst_conservation_error: f64,
    /// Materialised nodes whose matter was brought back into step with their
    /// own detail — see [`Tree::settle`]. Non-zero on a save taken mid-solve
    /// and zero on one taken at rest, which is exactly the shape of the defect
    /// it closes.
    pub settled: u64,
    /// Settlings where the detail turned out to say nothing new, so the matter
    /// was left exactly as it was. The idempotence guarantee, counted.
    pub settled_idempotent: u64,
    /// Nodes that stopped being one neighbourhood and became two — see
    /// [`Tree::resolve_extent`]. The operation `docs/BACKLOG.md` calls
    /// sibling-from-a-subset, which nothing could do until Phase 3.
    pub splits: u64,
    /// Sibling nodes whose contents became one neighbourhood again and were
    /// folded into one. Without it two clumps that fall back together stay two
    /// nodes for ever.
    pub merges: u64,
    /// Promoted children folded back into their parent because its detail was
    /// discarded — see [`Tree::shed_children`]. Before it existed each one of
    /// these was a live node nothing could reach and nothing would ever free.
    pub shed: u64,
}

impl Tree {
    pub fn new(world_seed: u64, root_agg: Matter, tier: Tier, spec: SampleSpec) -> Tree {
        let root = Node {
            key: PathKey::ROOT,
            parent: NodeIdx::NONE,
            slot: 0,
            depth: 0,
            tier,
            matter: root_agg,
            // A node that carries angular momentum is turning, and the frame is
            // where that is recorded. Leaving it at the default meant the one
            // node in every world that nobody promotes — the root — was the one
            // node that never rotated.
            motion: Motion {
                spin_rate: root_agg.angular_velocity(),
                orientation: crate::math::Quat::IDENTITY,
                ..Motion::default()
            },
            bodies: Vec::new(),
            potential: root_agg.gravitational_binding,
            gravity: Vec3::ZERO,
            children: Vec::new(),
            spec,
            epoch: 0,
            time: 0.0,
            last_disturbed: 0.0,
            last_solved: 0.0,
            last_grown: 0.0,
            residency: Residency::Speculative,
            pinned: false,
            contains_edit: false,
            rest_density: 0.0,
            unrest: 0.0,
            ocean: None,
            carried: 0.0,
            turning: Vec3::ZERO,
            ground: None,
            bubble: 1.0,
            alive: true,
            morphology: None,
            topology: None,
            steps_taken: 0,
            surface: None,
            surface_epoch: u32::MAX,
            last_report: SampleReport::default(),
        };
        Tree {
            nodes: vec![root],
            free: Vec::new(),
            root: NodeIdx(0),
            world_seed,
            persisted: HashMap::new(),
            stats: TreeStats::default(),
        }
    }

    /// Rebuild a tree from saved parts.
    ///
    /// The free list is *derived* rather than stored: a slot is free exactly
    /// when its node is not alive, so persisting it would be storing a fact the
    /// nodes already contain and risking the two disagreeing after a reload.
    pub fn restore(
        nodes: Vec<Node>,
        root: NodeIdx,
        world_seed: u64,
        persisted: HashMap<PathKey, Vec<Body>>,
        stats: TreeStats,
    ) -> Tree {
        let free = nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| !n.alive)
            .map(|(i, _)| i as u32)
            .collect();
        Tree { nodes, free, root, world_seed, persisted, stats }
    }

    pub fn get(&self, i: NodeIdx) -> &Node {
        &self.nodes[i.get()]
    }

    pub fn get_mut(&mut self, i: NodeIdx) -> &mut Node {
        &mut self.nodes[i.get()]
    }

    pub fn live_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.alive).count()
    }

    pub fn materialised_bodies(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| n.alive)
            .map(|n| n.bodies.len())
            .sum()
    }

    pub fn detail_bytes(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| n.alive)
            .map(|n| n.detail_bytes())
            .sum::<usize>()
            + self
                .persisted
                .values()
                .map(|v| v.len() * std::mem::size_of::<Body>())
                .sum::<usize>()
    }

    fn alloc(&mut self, node: Node) -> NodeIdx {
        if let Some(i) = self.free.pop() {
            self.nodes[i as usize] = node;
            NodeIdx(i)
        } else {
            self.nodes.push(node);
            NodeIdx((self.nodes.len() - 1) as u32)
        }
    }

    // -- materialisation --------------------------------------------------

    /// Produce this node's fine detail. Idempotent, and — for an unpinned node —
    /// bit-identical every time it is called at the same epoch.
    pub fn refine(&mut self, i: NodeIdx) -> &[Body] {
        // **An ocean is its description.** A node that carries one is a
        // planet's sea — a shell a thousand kilometres round and a few deep —
        // and its cells are what it holds. Drawn as bodies it would be a ball
        // of water parcels the size of a planet, which it is not; where water
        // is wanted at a finer scale, it is drawn in a patch of shore and
        // crosses in through that patch's open edge (`open_edge`).
        if self.nodes[i.get()].is_materialised() || self.nodes[i.get()].ocean.is_some() {
            return &self.nodes[i.get()].bodies;
        }
        let key = self.nodes[i.get()].key;

        // Pinned detail was altered by an interaction, so it cannot be
        // regenerated — it comes back from the persistent store instead.
        if let Some(saved) = self.persisted.get(&key) {
            let bodies = saved.clone();
            self.reconcile_children(i, bodies.len());
            let n = &mut self.nodes[i.get()];
            n.bodies = bodies;
            n.ground = None;
            self.stats.materialisations += 1;
            return &self.nodes[i.get()].bodies;
        }

        // The field the structure is proportioned in, derived here and stored,
        // so that the geometry and the number it was built with travel
        // together. A client regenerating this node reads the stored value
        // rather than working one out from a tree it only has part of.
        let gravity = self.gravity_at(i);
        self.nodes[i.get()].gravity = gravity;
        let (matter, spec, epoch, morph) = {
            let n = &self.nodes[i.get()];
            (n.matter, n.spec, n.epoch, n.morphology.clone())
        };
        // A structure is drawn at the node's own turning, in the axes its
        // layout is stated in, and the drawing then placed at the turn it has
        // reached (`Tree::drawn_turn`) — so what is redrawn is where the thing
        // has turned to, going round at the rate it turns.
        let turn = self.drawn_turn(i);
        let turning = turn.conjugate().rotate(self.angular_velocity(i));
        let setting = crate::sampler::Setting {
            gravity,
            rest_density: self.nodes[i.get()].rest_density,
            turning: Some(turning),
        };
        let (bodies, topo, report) = match &morph {
            Some(m) => {
                // Everything the drawing is closed against, in its own axes.
                let mut drawn = matter;
                if turn != crate::math::Quat::IDENTITY {
                    drawn.spin = turn.conjugate().rotate(matter.spin);
                    drawn.com = turn.conjugate().rotate(matter.com);
                }
                let (b, t, r) = crate::sampler::sample_structured_in(
                    &drawn,
                    m,
                    spec.count,
                    self.world_seed,
                    key.0,
                    epoch,
                    setting,
                );
                let mut b = b;
                if turn != crate::math::Quat::IDENTITY {
                    for body in b.iter_mut() {
                        body.pos = turn.rotate(body.pos);
                        body.vel = turn.rotate(body.vel);
                        body.spin = turn.rotate(body.spin);
                    }
                }
                (b, Some(t), r)
            }
            None => {
                let (b, r) =
                    crate::sampler::sample_in(&matter, spec, self.world_seed, key.0, epoch, setting);
                (b, None, r)
            }
        };
        self.stats.materialisations += 1;
        self.stats.bodies_created += bodies.len() as u64;
        self.stats.worst_conservation_error = self
            .stats
            .worst_conservation_error
            .max(report.conservation_error);
        // **A sea keeps its slot through a redraw.** A planet comes back from
        // a save with its children and none of its bodies, and its sea — which
        // holds a tide's history and its own books — is past the end of any
        // fresh draw, where `reconcile_children` would fold it back. The draw
        // is of the planet's whole matter, sea included, so the sea's water
        // is taken back out of it the way it was taken out the first time
        // (`Tree::withdraw`), and the slots between are left empty.
        let mut bodies = bodies;
        let seas: Vec<(usize, NodeIdx)> = self.nodes[i.get()]
            .children
            .iter()
            .enumerate()
            .skip(bodies.len())
            .filter(|(_, c)| !c.is_none() && self.nodes[c.get()].alive && self.nodes[c.get()].ocean.is_some())
            .map(|(s, c)| (s, *c))
            .collect();
        let last = seas.iter().map(|(s, _)| *s + 1).max().unwrap_or(0);
        if last > bodies.len() {
            bodies.resize(last, Body { mass: 0.0, ..Default::default() });
        }
        self.reconcile_children(i, bodies.len());
        // Contents held in a field by a structure are drawn near balance, not
        // in it — measured, a bucket's water drawn against its walls peaked at
        // 0.72 m/s settling — so they are solved until they say otherwise.
        let loose = topo.as_ref().map(|t| {
            (0..bodies.len()).any(|k| t.joints.get(k).map(|j| j.radius <= 0.0).unwrap_or(true))
        });
        let n = &mut self.nodes[i.get()];
        n.unrest = if gravity != Vec3::ZERO && loose == Some(true) { f64::INFINITY } else { 0.0 };
        n.bodies = bodies;
        n.topology = topo;
        n.potential = report.potential;
        n.last_report = report;
        for (slot, sea) in seas {
            let (mass, composition, radius) = {
                let m = &self.nodes[sea.get()].matter;
                (m.mass, m.composition, m.radius)
            };
            self.withdraw(i, mass, composition, radius);
            self.sync_from_child(i, slot, sea);
        }
        &self.nodes[i.get()].bodies
    }

    /// Turn one materialised body into a node of its own, one tier finer.
    ///
    /// The child's matter is *the body itself*, reinterpreted: same mass,
    /// same composition, same momentum in the parent's frame. Nothing is
    /// invented at this step — invention happens when the child is refined.
    /// State a thing that is already there: a new slot in `parent`, holding
    /// `body`, and the node it becomes.
    ///
    /// `promote` makes a node of something the parent already holds, and a
    /// parent whose contents are all members of its recipe — a patch of ground,
    /// every one of whose bodies is a cell or what is under them — has nothing
    /// to spare: making a cell into air takes the cell out of the ground. This
    /// is the other way in, the one an arrival from outside takes
    /// (`Tree::reparent`): a slot beyond the recipe's, loose among its members.
    ///
    /// **A scene's statement, not physics.** Nothing is conserved across it:
    /// the thing placed was not anywhere before.
    ///
    /// **What is placed is new matter where it is placed**, so the node it is
    /// placed in and every node holding that one hold it too: their mass grows
    /// by the body's, and the node it lands in takes its momentum. Left out, a
    /// face of an Earth holding a mass of air weighed what it did without it,
    /// and the world's books missed the part of the air's momentum that is the
    /// face carrying it — 2.7e20 kg m/s between two solves of the Earth.
    pub fn place(&mut self, parent: NodeIdx, body: Body, spec: SampleSpec) -> NodeIdx {
        self.refine(parent);
        self.nodes[parent.get()].matter.momentum += body.momentum();
        let mut up = parent;
        while !up.is_none() {
            self.nodes[up.get()].matter.mass += body.mass;
            up = self.nodes[up.get()].parent;
        }
        let slot = {
            let p = &mut self.nodes[parent.get()];
            p.bodies.push(body);
            while p.children.len() < p.bodies.len() {
                p.children.push(NodeIdx::NONE);
            }
            p.bodies.len() - 1
        };
        self.promote(parent, slot, spec)
    }

    /// Take `mass` out of a node's free bodies, each by its share of their
    /// mass, and hold it as a child of its own: what that share carried —
    /// momentum, angular momentum, internal energy — goes with it, so the
    /// books do not move. The child's slot is after everything drawn.
    ///
    /// How a planet's sea becomes a node (`World::assess_ocean`): the water
    /// was in the planet's bodies, spread through them by the draw, and it is
    /// the same water afterwards.
    pub fn hold_out(
        &mut self,
        parent: NodeIdx,
        mass: f64,
        composition: crate::state::Composition,
        radius: f64,
        spec: SampleSpec,
    ) -> NodeIdx {
        self.refine(parent);
        let Some(body) = self.withdraw(parent, mass, composition, radius) else { return NodeIdx::NONE };
        let slot = {
            let p = &mut self.nodes[parent.get()];
            p.bodies.push(body);
            while p.children.len() < p.bodies.len() {
                p.children.push(NodeIdx::NONE);
            }
            p.bodies.len() - 1
        };
        let child = self.promote(parent, slot, spec);
        // What it carried, exactly: `promote` would floor the heat at the
        // new node's own thermal energy.
        if !child.is_none() {
            self.nodes[child.get()].matter.internal_energy = body.internal_energy;
        }
        child
    }

    /// Take `mass` out of a node's free bodies by their shares of it, and
    /// return it as one body: at their shares' centre of mass, moving with
    /// their momentum, turning with their angular momentum about that centre,
    /// and holding their share of the heat. `None` where the free bodies do
    /// not hold that much.
    pub fn withdraw(
        &mut self,
        parent: NodeIdx,
        mass: f64,
        composition: crate::state::Composition,
        radius: f64,
    ) -> Option<Body> {
        let n = &mut self.nodes[parent.get()];
        let free: Vec<usize> = (0..n.bodies.len()).filter(|&s| n.child_of(s).is_none() && n.bodies[s].mass > 0.0).collect();
        let total: f64 = free.iter().map(|&s| n.bodies[s].mass).sum();
        if !(mass > 0.0) || !(total > mass) {
            return None;
        }
        let (mut p, mut x, mut l, mut u) = (Vec3::ZERO, Vec3::ZERO, Vec3::ZERO, 0.0);
        let mut temperature = 0.0;
        for &s in &free {
            let b = &mut n.bodies[s];
            let share = mass * b.mass / total;
            let f = share / b.mass;
            p += b.vel.scale(share);
            x += b.pos.scale(share);
            l += b.pos.cross(b.vel.scale(share)) + b.spin.scale(f);
            u += b.internal_energy * f;
            temperature += b.temperature * share;
            b.mass -= share;
            b.internal_energy *= 1.0 - f;
            b.spin = b.spin.scale(1.0 - f);
        }
        let (at, v) = (x.scale(1.0 / mass), p.scale(1.0 / mass));
        Some(Body {
            pos: at,
            vel: v,
            mass,
            radius,
            temperature: temperature / mass,
            composition,
            internal_energy: u,
            spin: l - at.cross(p),
            ..Default::default()
        })
    }

    pub fn promote(&mut self, i: NodeIdx, slot: usize, spec: SampleSpec) -> NodeIdx {
        self.refine(i);
        // A new thing among contents held in a field is a new balance to find.
        // See `Node::unrest`.
        if self.nodes[i.get()].gravity != Vec3::ZERO && self.nodes[i.get()].topology.is_some() {
            self.nodes[i.get()].unrest = f64::INFINITY;
        }
        {
            let n = &self.nodes[i.get()];
            if slot >= n.bodies.len() {
                return NodeIdx::NONE;
            }
            let existing = n.child_of(slot);
            if !existing.is_none() {
                return existing;
            }
        }
        let (body, key, depth, parent_tier) = {
            let n = &self.nodes[i.get()];
            (n.bodies[slot], n.key.child(slot as u64), n.depth + 1, n.tier)
        };
        // The child's tier follows from its size, not from its depth. A cloud
        // that splits into clumps is still `Stellar`; only when the pieces get
        // small enough that a different physics applies does the tier change.
        let tier = tier_for(body.radius, parent_tier);

        // The policy has to match the tier, and the caller cannot know the tier
        // until the radius is in hand. A caller extrapolating from the parent —
        // "my child will be one tier finer than me" — is guessing, and when the
        // guess is wrong the node materialises under a policy meant for a
        // different scale: eight thousand molecules packed into a node the size
        // of an atom, at three thousandths of their own interaction radius,
        // where the potential is 10^44 m/s^2 of acceleration and the solver has
        // no answer but to explode.
        //
        // So a spec meant for a *coarser* scale than the child turned out to be
        // is replaced by the tier's own policy, with the caller's count read as
        // a budget. A spec meant for a finer one is kept: asking to split an
        // atom into nucleons is a deliberate step down and not a mistake, and
        // overriding it would leave the ladder unable to reach its own bottom.
        let spec = spec_for(tier, spec);

        let mut matter = Matter::neutral(body.mass, body.radius.max(1e-30), body.temperature, body.composition);
        matter.charge = body.charge;
        matter.spin = body.spin;
        // What the parent is made of is what its contents are made of.
        // `docs/PLAY.md` D17's downward half: a `Body` carries no speciation of
        // its own — 200 bytes on a 184-byte struct is not a trade §5A.5 makes —
        // and it does not need to, because `sample` draws every body from one
        // matter. So a body becoming a node inherits that matter's mixture,
        // which is the same answer a body-level field would have given.
        //
        // Mass fractions, so nothing is scaled: a kilogram of a node that is a
        // quarter salt is a quarter salt.
        matter.mixture = self.nodes[i.get()].matter.mixture;
        // The child's own frame carries the bulk motion, so inside its frame the
        // net momentum is zero — that is what "rest frame" means. Bulk motion is
        // never double-counted.
        matter.momentum = Vec3::ZERO;
        matter.internal_energy = body.internal_energy.max(matter.thermal_energy());
        matter.luminosity = crate::state::stefan_boltzmann(matter.radius, matter.temperature);
        // **A part of a structure turns with the structure**: it was drawn
        // going round at the structure's rate, and has no turning of its own
        // relative to it. Anything else — a parcel of a cloud — turns at what
        // its angular momentum says a sphere of it does, less its parent's.
        let turning = if self.nodes[i.get()].morphology.is_some() {
            Vec3::ZERO
        } else {
            let own = matter.angular_velocity() - self.angular_velocity(i);
            self.facing(i).conjugate().rotate(own)
        };

        let child = Node {
            key,
            parent: i,
            slot: slot as u32,
            depth,
            tier,
            matter,
            motion: Motion {
                offset: body.pos,
                velocity: body.vel,
                // A body's spin becomes the node's rotation, and **the way the
                // body was facing becomes the way the node is facing**.
                //
                // This used to be `Quat::IDENTITY` with a comment explaining
                // that the child's body frame *is* how it was sampled — which
                // was true only because a `Body` had no orientation to hand
                // over. It does now, so a part promoted out of a recipe arrives
                // turned the way the recipe turned it, and a wall that comes
                // off a box is not silently squared up to its parent's axes on
                // the way out. A sampled body's orientation is identity, so
                // nothing that was right before changes.
                orientation: body.orientation,
                spin_rate: turning,
                proper_time: self.nodes[i.get()].motion.proper_time,
            },
            bodies: Vec::new(),
            potential: 0.0,
            gravity: Vec3::ZERO,
            children: Vec::new(),
            spec,
            epoch: 0,
            time: self.nodes[i.get()].time,
            last_disturbed: self.nodes[i.get()].time,
            last_solved: self.nodes[i.get()].time,
            last_grown: self.nodes[i.get()].time,
            residency: Residency::Speculative,
            pinned: false,
            contains_edit: false,
            // Drawn from the parent, so made of what the parent is made of.
            rest_density: self.nodes[i.get()].rest_density,
            unrest: 0.0,
            ocean: None,
            carried: self.nodes[i.get()].time,
            turning: Vec3::ZERO,
            ground: None,
            bubble: 1.0,
            alive: true,
            morphology: None,
            topology: None,
            steps_taken: 0,
            surface: None,
            surface_epoch: u32::MAX,
            last_report: SampleReport::default(),
        };
        let idx = self.alloc(child);
        self.nodes[i.get()].children[slot] = idx;
        self.stats.promotions += 1;
        self.inherit_recipe(i, slot, idx);
        idx
    }

    /// Hand a promoted cell the recipe its parent's recipe says it has.
    ///
    /// **This is what makes a surface a tree.** `docs/PLAY.md` D6: a planetary
    /// node "on refinement divides its sphere into patches", and a patch's
    /// children are four squares of it. The bodies a patch holds *are* those
    /// squares, so promoting one has to give the new node the piece of the
    /// parameterisation it stands for — otherwise it arrives as an anonymous
    /// lump of rock and the ladder stops at the first step.
    ///
    /// Nothing is stored per patch except its address, so the descent generates
    /// what it needs and regenerates it bit-for-bit on the way back: the
    /// address is exact integer arithmetic on a face index and two cell
    /// indices, twenty-four levels deep on an Earth.
    ///
    /// It issues no identity, which matters: which nodes a frame promotes
    /// depends on a wall-clock allowance, and `next_entity` is persisted. See
    /// `World::cross`, which promotes for the same reason and is careful the
    /// same way.
    fn inherit_recipe(&mut self, parent: NodeIdx, slot: usize, child: NodeIdx) {
        let Some((program, recipe)) = self.nodes[parent.get()].morphology.as_ref().and_then(|m| {
            let r = m.recipe.as_ref()?;
            let crate::recipe::Recipe::Tiled(t) = r else { return None };
            Some((m.program, t.child(slot)?))
        }) else {
            return;
        };
        let key = self.nodes[child.get()].key;
        let seed = self.world_seed;
        let mut m = crate::morph::Morphology::new(program, seed, key.0, 0);
        m.built = self.nodes[child.get()].matter.mass;
        m.design_mass = m.built;
        m.progress = 1.0;
        m.recipe = Some(crate::recipe::Recipe::Tiled(recipe));
        let extent = m.extent();
        let n = &mut self.nodes[child.get()];
        n.matter.radius = extent.max(1e-30);
        n.morphology = Some(m);
        self.stats.structures += 1;
        self.retier(child);
    }

    /// Fold fine detail back into the node's matter and free it.
    ///
    /// The conserved tuple is measured before and after; the difference is
    /// recorded in `stats.worst_conservation_error` and asserted on in the
    /// tests. If a solver has been sloppy, this is where it shows up.
    pub fn coarsen(&mut self, i: NodeIdx) -> f64 {
        if !self.nodes[i.get()].is_materialised() {
            return 0.0;
        }
        // Pull any promoted child's evolved state back into its body first,
        // otherwise work done at a finer tier is silently discarded.
        //
        // Speciation comes back here too, and this is the only place a blend is
        // needed at all. `docs/PLAY.md` D17's upward half: `summarise` reads
        // bodies, and a `Body` carries no mixture, so blending identical
        // descriptions would be the identity and buy nothing. What differs is a
        // *promoted child*, which is a node, has a mixture of its own, and may
        // have reacted its way somewhere the parent has not — ice that melted,
        // salt that dissolved. Those are collected before they are released and
        // blended into the parent by mass.
        let children = self.nodes[i.get()].children.clone();
        let mut speciation: Option<(crate::chem::Mixture, f64)> = None;
        let mut child_mass = 0.0;
        // Collected before the children are released, because it is measured
        // off them. See `Tree::unrepresented`.
        let mut held = 0.0;
        for (slot, c) in children.iter().enumerate() {
            if !c.is_none() {
                let (mix, mass) = {
                    let n = &self.nodes[c.get()];
                    (n.matter.mixture, n.matter.mass.max(0.0))
                };
                if !mix.is_empty() && mass > 0.0 {
                    speciation = Some(match speciation {
                        None => (mix, mass),
                        Some((acc, m)) => {
                            let (blend, lost) = crate::chem::Mixture::blend(&acc, m, &mix, mass);
                            self.note_description_lost(lost);
                            (blend, m + mass)
                        }
                    });
                }
                child_mass += mass;
                held += self.unrepresented(*c);
                self.sync_from_child(i, slot, *c);
                self.release_subtree(*c);
            }
        }

        let (before, potential, pinned, key) = {
            let n = &self.nodes[i.get()];
            (n.matter.conserved(), n.potential, n.pinned, n.key)
        };
        let bodies = std::mem::take(&mut self.nodes[i.get()].bodies);
        let (matter, err) =
            self.summarised(i, &bodies, potential, before, speciation, child_mass, held);

        if pinned {
            self.stats.persisted_bodies += bodies.len() as u64;
            self.persisted.insert(key, bodies);
        } else {
            self.stats.bodies_discarded += bodies.len() as u64;
        }

        // If the detail did not actually change the matter — the usual case
        // when a user simply pans away — keep the coarse state as the
        // authority rather than overwriting it with a summarising that differs
        // only by round-off.
        //
        // This is what makes "leave and come back" *exactly* idempotent rather
        // than merely accurate. Without it, every visit perturbs the matter
        // in the last bits, the next materialisation samples from a marginally
        // different distribution, and a region a user visits a thousand times
        // slowly drifts away from itself. With it, a region nobody has
        // disturbed is bit-for-bit the region they left.
        //
        // **The mixture is part of the matter** (D17), so it has to round-trip
        // too before the coarse state may stand. A promoted child that reacted
        // — salt that dissolved — changes what the node is made of without
        // moving a single conserved quantity, and the conserved tuple alone
        // let the early return discard the blend. It was hidden for as long
        // as a sampled condensed node's round trip happened to miss the
        // tolerance; drawing its parcels as macroscopic objects made it exact,
        // and `what_a_node_is_made_of_travels_down_and_back_up` failed.
        let unchanged = matter.mixture.same_as(&self.nodes[i.get()].matter.mixture);
        if err < IDEMPOTENT_TOLERANCE && !pinned && unchanged {
            self.stats.coarsenings += 1;
            self.stats.idempotent_coarsenings += 1;
            let n = &mut self.nodes[i.get()];
            n.children.clear();
            return err;
        }

        self.adopt(i, matter);
        self.nodes[i.get()].children.clear();
        self.stats.coarsenings += 1;
        self.stats.worst_conservation_error = self.stats.worst_conservation_error.max(err);
        err
    }

    /// Write what a node's detail has become back into its matter, and keep
    /// the detail.
    ///
    /// [`Self::coarsen`] without the destruction, and the two share every line
    /// that decides what the matter becomes. While a node is materialised its
    /// **bodies are the authority and its matter is the summary made when it
    /// last coarsened**; that is the design and it costs nothing, right up to
    /// the moment something reads the matter as though it were current. A save
    /// is that moment: it writes the matter and discards the bodies of every
    /// unpinned node, so whatever the solver did since materialisation is
    /// lost.
    ///
    /// Measured on the reference world before this existed: 2.35x10^-8 of the
    /// root's energy on the one node the scheduler was solving every frame,
    /// against 1x10^-15 or better on every node it had finished with. Zero at
    /// rest, which is why it stayed invisible — and `docs/PLAY.md` D16 makes
    /// exactly the mid-flight checkpoint routine, because a crossing is a save
    /// point in all but name.
    ///
    /// The alternative was keeping the matter in step inside the solver, which
    /// costs a `summarise` per solved node per frame; this costs one per
    /// materialised node per *save*.
    pub fn settle(&mut self, i: NodeIdx) -> f64 {
        if i.is_none() || !self.nodes[i.get()].alive || !self.nodes[i.get()].is_materialised() {
            return 0.0;
        }
        // A promoted child is the real thing and its body is a stand-in, so the
        // stand-ins have to be current before anything summarises them. This is
        // the same first step `coarsen` takes; what it does not do is release
        // the children afterwards, because nothing here is going away.
        self.sync_children(i);
        let held: f64 = self.nodes[i.get()]
            .children
            .clone()
            .iter()
            .map(|c| self.unrepresented(*c))
            .sum();
        let (before, potential, pinned) = {
            let n = &self.nodes[i.get()];
            (n.matter.conserved(), n.potential, n.pinned)
        };
        // Taken and put back rather than cloned: a node being checkpointed can
        // hold tens of thousands of bodies, and a save that allocated a second
        // copy of the whole world's detail would be a worse bargain than the
        // staleness it is fixing.
        let bodies = std::mem::take(&mut self.nodes[i.get()].bodies);
        let (matter, err) = self.summarised(i, &bodies, potential, before, None, 0.0, held);
        self.nodes[i.get()].bodies = bodies;
        // The same idempotence rule `coarsen` uses, and for the same reason: a
        // node nobody has disturbed must come back bit-for-bit, so matter that
        // differs from its own detail only by round-off is left alone.
        if err < IDEMPOTENT_TOLERANCE && !pinned {
            self.stats.settled_idempotent += 1;
            return err;
        }
        self.adopt(i, matter);
        self.stats.settled += 1;
        self.stats.worst_conservation_error = self.stats.worst_conservation_error.max(err);
        err
    }

    /// What a node's matter would be if its detail were folded back into it.
    ///
    /// The half of the scale transform that decides *what the coarse state
    /// becomes*, with no opinion about whether the detail survives. Both
    /// callers need every line of it and having had it written twice is how
    /// the two would drift apart.
    ///
    /// `speciation` is the mixture blended out of promoted children that are
    /// being released, with the mass it came from; [`Self::settle`] passes
    /// `None` because its children are staying and still speak for themselves.
    fn summarised(
        &mut self,
        i: NodeIdx,
        bodies: &[Body],
        potential: f64,
        before: crate::state::Conserved,
        speciation: Option<(crate::chem::Mixture, f64)>,
        child_mass: f64,
        held: f64,
    ) -> (Matter, f64) {
        let mut matter = summarise(bodies, potential);
        // What the stand-ins could not carry. See `Tree::unrepresented`.
        matter.internal_energy += held;
        matter.external_potential = self.nodes[i.get()].matter.external_potential;
        matter.chemical_energy = self.nodes[i.get()].matter.chemical_energy;
        // `summarise` measures where the bodies ended up, and a bond is far
        // below the scale of a body, so it reports no cohesive binding at all.
        // Reinstated here rather than after the error is measured, because the
        // error is an energy comparison and a granite block's cohesive energy
        // is the largest term in it.
        matter.cohesive_binding = self.nodes[i.get()].matter.cohesive_binding;
        // The node's own description covers whatever was not promoted; the
        // children cover the rest. With nothing promoted this is the node's own
        // mixture unchanged, which is what makes leaving and coming back
        // idempotent for chemistry as well as for matter.
        matter.mixture = match speciation {
            None => self.nodes[i.get()].matter.mixture,
            Some((child_mix, mass_from_children)) => {
                let own = self.nodes[i.get()].matter;
                let rest = (own.mass - child_mass).max(0.0);
                if own.mixture.is_empty() || rest <= 0.0 {
                    child_mix
                } else {
                    let (blend, lost) =
                        crate::chem::Mixture::blend(&own.mixture, rest, &child_mix, mass_from_children);
                    self.note_description_lost(lost);
                    blend
                }
            }
        };
        matter.entropy_exported = self.nodes[i.get()].matter.entropy_exported;
        let scales = crate::state::Scales::of(bodies);
        let err = matter.conserved().error_against(&before, &scales);
        (matter, err)
    }

    /// Take a summarised matter as the node's own.
    ///
    /// Field by field rather than wholesale, because a node is more than what
    /// `summarise` can see: its tier, its spec and its identity are not
    /// measurements of its contents, and two of the quantities that are need a
    /// rule of their own.
    fn adopt(&mut self, i: NodeIdx, matter: Matter) {
        let arrived = matter.spin - self.nodes[i.get()].matter.spin;
        let n = &mut self.nodes[i.get()];
        // Preserve the node's own frame-level bookkeeping: `summarise` measures
        // the children in the node's frame, so the node's momentum and com are
        // updated, but its tier, spec and identity are untouched.
        n.matter.mass = matter.mass;
        n.matter.com = matter.com;
        n.matter.momentum = matter.momentum;
        n.matter.spin = matter.spin;
        n.matter.internal_energy = matter.internal_energy;
        n.matter.gravitational_binding = matter.gravitational_binding;
        n.matter.cohesive_binding = matter.cohesive_binding;
        n.matter.mixture = matter.mixture;
        n.matter.radius = matter.radius;
        n.matter.temperature = matter.temperature;
        n.matter.composition = matter.composition;
        n.matter.charge = matter.charge;
        n.matter.baryon_number = matter.baryon_number;
        n.matter.lepton_number = matter.lepton_number;
        // Entropy needs two corrections that the original one-liner got wrong
        // as soon as anything in the world could become more ordered.
        //
        // First, coarse-graining itself may only *increase* entropy — you know
        // less once the detail is gone — but the quantity that is monotonic is
        // the total, local plus exported. A structure that grew since the last
        // visit has legitimately lowered its local entropy, and clamping it back
        // up would silently destroy the record of that and unbalance the books.
        //
        // Second, for a structured node `summarise` is not entitled to an
        // opinion at all: it sees an unstructured heap of parts and reports the
        // entropy of the same mass as a gas, which erases precisely the order
        // that makes the thing a structure. `Body` carries no topology, so the
        // information is not there to be recovered. The developmental state is
        // the authority.
        if n.morphology.is_none() {
            let total_stored = n.matter.total_entropy();
            let total_summarised = matter.entropy + n.matter.entropy_exported;
            if total_summarised >= total_stored {
                n.matter.entropy = matter.entropy;
            }
        }
        n.matter.luminosity = matter.luminosity;
        // The morphology owns the structure's size, for the same reason it owns
        // its entropy: `summarise` measures the parts, but what the parts add up
        // to is the program's business.
        if let Some(m) = &n.morphology {
            n.matter.radius = m.extent().max(1e-30);
        }
        // The matter this node holds has just been rewritten from its own
        // detail, and whatever angular momentum its contents gained or lost
        // since turns it by that much (`Tree::turn_by`).
        self.turn_by(i, arrived);
    }

    /// The energy a node has that the body standing in for it cannot carry.
    ///
    /// **A `Body` has a mass, a velocity and an internal energy, and that is
    /// all.** It has nowhere to put the binding holding a thing together, the
    /// grip of something that is not being refined, or the free energy a
    /// structure is storing — which is exactly what `summarise` means when it
    /// says those terms are "not knowable from the children alone". So a parent
    /// summarised from a body list in which a promoted child appears as a
    /// stand-in comes out *high* by the child's binding, and its own stand-in
    /// then carries the inflated figure one level further up.
    ///
    /// Measured on the reference world before this existed: node 2 holds
    /// -1.27x10^48 J of binding, nodes 0 and 1 each read 1.25x10^48 J high, and
    /// the world's energy moved by 1.52x10^-8 across a save — the same term
    /// arriving once per level that resamples from a summary.
    ///
    /// The fix is here, in the summary, and deliberately **not** in the body:
    /// folding a negative binding into a stand-in's `internal_energy` would
    /// hand the tier solver a parcel with less energy than its own rest mass
    /// implies, and the equation of state would price it as something that does
    /// not exist.
    ///
    /// A materialised node answers with its measured `potential`, because that
    /// is the figure [`Self::sum_conserved`] uses for it; an unmaterialised one
    /// answers with the `gravitational_binding` its matter carries.
    pub fn unrepresented(&self, i: NodeIdx) -> f64 {
        if i.is_none() || !self.nodes[i.get()].alive {
            return 0.0;
        }
        let n = &self.nodes[i.get()];
        let gravity = if n.is_materialised() {
            n.potential
        } else {
            n.matter.gravitational_binding
        };
        gravity + n.matter.cohesive_binding + n.matter.external_potential + n.matter.chemical_energy
    }

    /// Write a promoted child's evolved matter back into the parent's body.
    fn sync_from_child(&mut self, parent: NodeIdx, slot: usize, child: NodeIdx) {
        let (mass, comp, temp, charge, spin, internal, radius, frame) = {
            let c = &self.nodes[child.get()];
            (
                c.matter.mass,
                c.matter.composition,
                c.matter.temperature,
                c.matter.charge,
                c.matter.spin,
                c.matter.internal_energy,
                c.matter.radius,
                c.motion,
            )
        };
        // Where the child is at the parent's instant, which is the instant the
        // parent's contents are at: a parent that is behind the world solves
        // against where its children were then. See `Node::carried`.
        let at = self.position_at(child, self.nodes[parent.get()].time);
        let self_velocity_at = self.velocity_at(child, self.nodes[parent.get()].time);
        let _ = frame.offset;
        let p = &mut self.nodes[parent.get()];
        if let Some(b) = p.bodies.get_mut(slot) {
            b.mass = mass;
            b.composition = comp;
            b.temperature = temp;
            b.charge = charge;
            b.spin = spin;
            b.internal_energy = internal;
            b.radius = radius;
            b.pos = at;
            // And how fast it is going then, in the same frame: a held child's
            // velocity turns with the ground, and its velocity now against its
            // position then is turned 0.66 rad out of step for an Earth whose
            // contents are 9000 s behind.
            b.vel = self_velocity_at;
            // And the way it has turned, which is the other half of `promote`
            // handing its orientation down: a part that came off a box, was
            // knocked askew and then rejoined comes back askew.
            b.orientation = frame.orientation;
        }
    }

    /// Carry a node's motion to a world instant: where it is, in closed form,
    /// at constant velocity and spin. The only thing that moves a node's
    /// motion forward. See `Node::carried`.
    pub fn carry(&mut self, i: NodeIdx, instant: f64) {
        let about = self.turning_about(i, self.nodes[i.get()].carried);
        let n = &mut self.nodes[i.get()];
        let dt = instant - n.carried;
        if dt > 0.0 && dt.is_finite() {
            let (offset, velocity) = (n.motion.offset, n.motion.velocity);
            n.motion.advance(dt);
            if n.turning != Vec3::ZERO {
                let (at, v) = turning_carry_about(n.turning, about, offset, velocity, dt);
                n.motion.offset = at;
                n.motion.velocity = v;
            }
            n.carried = instant;
        }
    }

    /// Where a node is at an instant, in its parent's frame: its carried
    /// offset moved at its own velocity by the difference. Exact for the
    /// constant velocity it is carried at, in either direction — which is what
    /// a parent solving at an earlier instant than the world's needs to read.
    pub fn position_at(&self, i: NodeIdx, instant: f64) -> Vec3 {
        let n = &self.nodes[i.get()];
        let dt = instant - n.carried;
        if !dt.is_finite() {
            return n.motion.offset;
        }
        if n.turning != Vec3::ZERO {
            return turning_carry_about(n.turning, self.turning_about(i, n.carried), n.motion.offset, n.motion.velocity, dt).0;
        }
        n.motion.offset + n.motion.velocity.scale(dt)
    }

    /// Where a held node's parent's centre is at an instant, relative to the
    /// centre the node goes round, m, root-aligned: zero where the parent is
    /// the turning body itself, and otherwise the parent's own place at that
    /// instant from the same centre. A thing held by a patch of ground goes
    /// round the planet, not the patch — measured with its great circle taken
    /// about the patch's centre, air over a face of an Earth fell 955 km while
    /// the face caught up 1200 s — and where the patch is has to be read at
    /// the same instant: a copy taken at the frame's start was 21 km of turning
    /// out by its end.
    pub fn turning_about(&self, i: NodeIdx, instant: f64) -> Vec3 {
        let n = &self.nodes[i.get()];
        if n.turning == Vec3::ZERO || n.parent.is_none() {
            return Vec3::ZERO;
        }
        let p = n.parent;
        if self.nodes[p.get()].turning == Vec3::ZERO {
            return Vec3::ZERO;
        }
        self.turning_about(p, instant) + self.position_at(p, instant)
    }

    /// Where a node is at an instant relative to an ancestor: its place and
    /// each of its ancestors' up to that one, all read at the same instant —
    /// [`Tree::offset_from`] as it was then rather than as it has been carried.
    pub fn offset_at(&self, ancestor: NodeIdx, mut node: NodeIdx, instant: f64) -> Vec3 {
        let mut at = Vec3::ZERO;
        while node != ancestor && !node.is_none() {
            at += self.position_at(node, instant);
            node = self.nodes[node.get()].parent;
        }
        at
    }

    /// How fast a node is going at an instant, in its parent's frame — the
    /// companion of [`Tree::position_at`], carried the same way.
    pub fn velocity_at(&self, i: NodeIdx, instant: f64) -> Vec3 {
        let n = &self.nodes[i.get()];
        let dt = instant - n.carried;
        if dt.is_finite() && n.turning != Vec3::ZERO {
            return turning_carry_about(n.turning, self.turning_about(i, n.carried), n.motion.offset, n.motion.velocity, dt).1;
        }
        n.motion.velocity
    }

    /// Change a node's velocity by `dv`, a change that happened at `instant`
    /// — its parent's, while the node's motion has been carried to a later
    /// one. Where it was then is taken back, the change made there, and the
    /// whole carried forward again the way the node is carried: in a straight
    /// line, or turning with what holds it, which also turns the change.
    pub fn kick(&mut self, i: NodeIdx, dv: Vec3, instant: f64) {
        if !dv.is_finite() || dv == Vec3::ZERO {
            return;
        }
        let since = self.nodes[i.get()].carried - instant;
        let turning = self.nodes[i.get()].turning;
        if turning == Vec3::ZERO || !(since.abs() > 0.0) || !since.is_finite() {
            let n = &mut self.nodes[i.get()];
            n.motion.velocity = n.motion.velocity + dv;
            if since.is_finite() {
                n.motion.offset = n.motion.offset + dv.scale(since.max(0.0));
            }
            return;
        }
        let (at, v) = (self.position_at(i, instant), self.velocity_at(i, instant));
        let about = self.turning_about(i, instant);
        let (at, v) = turning_carry_about(turning, about, at, v + dv, since);
        let n = &mut self.nodes[i.get()];
        n.motion.offset = at;
        n.motion.velocity = v;
    }

    /// Free a node and everything under it.
    pub fn release_subtree(&mut self, i: NodeIdx) {
        if i.is_none() || !self.nodes[i.get()].alive {
            return;
        }
        let children = self.nodes[i.get()].children.clone();
        for c in children {
            if !c.is_none() {
                self.release_subtree(c);
            }
        }
        let n = &mut self.nodes[i.get()];
        if n.pinned && !n.bodies.is_empty() {
            let bodies = std::mem::take(&mut n.bodies);
            let key = n.key;
            self.stats.persisted_bodies += bodies.len() as u64;
            self.persisted.insert(key, bodies);
        } else {
            self.stats.bodies_discarded += n.bodies.len() as u64;
            n.bodies.clear();
        }
        let n = &mut self.nodes[i.get()];
        n.children.clear();
        n.alive = false;
        if i != self.root {
            self.free.push(i.0);
        }
    }

    /// Release every promoted child back into this node, before its detail is
    /// discarded.
    ///
    /// **The orphan.** `docs/PLAY.md` §7 Phase 3's first correction, and the
    /// phase cannot meet its done-when with it open. Every structural-change
    /// path in the engine says "my body list is stale" by clearing it — and
    /// said it by throwing away the only record of what had been promoted out
    /// of it. Measured: promote a limb out of a forty-year-old tree, load the
    /// tree to failure, and the limb's node is still `alive` with `parent`
    /// pointing at the tree while the tree's `children` array is empty.
    /// Unreachable from any walk, never freed, still holding an arena slot and
    /// still being scheduled.
    ///
    /// The answer is that a body list is not the only thing that goes stale: a
    /// child's *slot* does too, because the recipe that produced it has
    /// changed. So the child is folded back the way [`Self::coarsen`] folds one
    /// back — its evolved state written into the body that stood for it, and
    /// then released — which is collapsing rather than forgetting and conserves
    /// exactly, because the node's matter counted that child all along.
    ///
    /// **Not a crossing, deliberately.** D16's outward case is the mechanism
    /// for a child that has *left*, and it already runs over every node every
    /// frame — so by the time a structure changes, anything that had drifted
    /// out has been re-homed by the pass rather than by this. Re-homing here as
    /// well would take the child's mass out of the tree while the matter that
    /// still counts it was about to be regenerated from the recipe, which
    /// creates it twice.
    pub fn shed_children(&mut self, i: NodeIdx) -> usize {
        if i.is_none() || !self.nodes[i.get()].alive {
            return 0;
        }
        let children = self.nodes[i.get()].children.clone();
        let mut shed = 0;
        for (slot, c) in children.iter().enumerate() {
            if c.is_none() || !self.nodes[c.get()].alive {
                continue;
            }
            self.sync_from_child(i, slot, *c);
            self.release_subtree(*c);
            shed += 1;
        }
        self.nodes[i.get()].children.clear();
        self.stats.shed += shed as u64;
        shed
    }

    /// Make the children list fit a body list of `len`, without dropping any.
    ///
    /// [`Self::refine`] used to assign `vec![NONE; len]` outright, which is
    /// correct for the case it was written for — a node with no detail has
    /// nothing promoted out of it — and wrong for the case a save creates. The
    /// wire format writes `children` and, for an unpinned node, no bodies at
    /// all, so a reloaded node comes back unmaterialised *with* promoted
    /// children, and the next refinement wiped them. The same orphan as the
    /// structural paths, by a different road.
    fn reconcile_children(&mut self, i: NodeIdx, len: usize) {
        let children = self.nodes[i.get()].children.clone();
        for (slot, c) in children.iter().enumerate().skip(len) {
            if c.is_none() || !self.nodes[c.get()].alive {
                continue;
            }
            // Past the end of the new body list there is no slot to stand in,
            // so this one is folded back like any other stale slot.
            self.sync_from_child(i, slot, *c);
            self.release_subtree(*c);
            self.stats.shed += 1;
        }
        let n = &mut self.nodes[i.get()];
        n.children.resize(len, NodeIdx::NONE);
    }

    /// Pull every promoted child's evolved state into the body that stands
    /// for it.
    ///
    /// `children` runs parallel to `bodies`, and a promoted child is the real
    /// thing while its body is a stand-in. Until this ran every frame the two
    /// diverged silently: `promote` set the child's `Motion` from the body and
    /// nothing ever wrote `motion.velocity` again, so a promoted node moved
    /// ballistically in its parent's frame for the rest of its life while the
    /// parent's solver went on integrating the body it came from. Measured, in
    /// `docs/BACKLOG.md`: 0.79 of the child's own radius after forty frames.
    ///
    /// Returns the slots that hold a promoted child, so a caller that is about
    /// to run a solver knows which bodies are stand-ins.
    pub fn sync_children(&mut self, parent: NodeIdx) -> Vec<usize> {
        let children = self.nodes[parent.get()].children.clone();
        let mut slots = Vec::new();
        for (slot, c) in children.iter().enumerate() {
            if c.is_none() || !self.nodes[c.get()].alive {
                continue;
            }
            self.sync_from_child(parent, slot, *c);
            slots.push(slot);
        }
        slots
    }

    /// Hand each promoted child the velocity its body was just given.
    ///
    /// This is the half that did not exist at all. The parent's solver computes
    /// a force on every body it holds, including the stand-ins; without this
    /// the force lands on the stand-in and is thrown away on the next
    /// [`Self::sync_children`], so two promoted things in one frame could not
    /// attract, collide or perturb each other. With it, the *change* the solver
    /// made is the force, and the child integrates it through its own `Motion`.
    ///
    /// Position is deliberately not copied back. The child owns where it is —
    /// that is what "the child is the real thing" means — and its offset is
    /// carried forward by `Motion::advance` on the world's own clock. Taking
    /// the solver's position too would integrate the same motion twice.
    ///
    /// So the two do not agree exactly between syncs, and that is expected
    /// rather than a leftover of the old defect. The parent's solver integrates
    /// with leapfrog — `x + v·dt + ½a·dt²` — and the child coasts linearly, so
    /// they differ by the acceleration term until the next
    /// [`Self::sync_children`] puts the stand-in back where the child is. What
    /// changed is that the difference is now *reset* rather than accumulated:
    /// measured over four hundred frames it oscillates between 0.13 and 0.29 of
    /// the child's radius, where before it passed 0.79 in forty and kept
    /// climbing.
    pub fn apply_body_forces(&mut self, parent: NodeIdx, before: &[(usize, crate::math::Vec3)]) {
        self.hand_back(parent, before, None);
    }

    /// Hand what a parent's solve did to its stand-ins back to the children
    /// they stand for.
    ///
    /// With `adopt` of `None`, the change in each stand-in's velocity, made at
    /// the parent's instant (`Tree::kick`) — the rule for a solver whose
    /// stand-ins are points it pushes. With `Some(end)`, the stand-in's whole
    /// state, which the solve integrated to the instant `end`, carried from
    /// there to wherever the child's motion is: **the solve is the account of
    /// the child's motion over the step**, so nothing integrates it twice. A
    /// ground solve is that kind. Handing only the velocity back and letting
    /// the child's own carry move it over the same step counted the
    /// centripetal of a turning planet twice — once in the ground's support and
    /// once in the carry — and a face of a turning Earth went from 4.85e6 m to
    /// 2.9e8 m in a day.
    pub fn hand_back(&mut self, parent: NodeIdx, before: &[(usize, crate::math::Vec3)], adopt: Option<f64>) {
        for (slot, was) in before {
            let Some(child) = self.nodes[parent.get()].children.get(*slot).copied() else {
                continue;
            };
            if child.is_none() || !self.nodes[child.get()].alive {
                continue;
            }
            let Some(now) = self.nodes[parent.get()].bodies.get(*slot).map(|b| b.vel) else {
                continue;
            };
            let dv = now - *was;
            if !dv.is_finite() {
                continue;
            }
            if let Some(end) = adopt {
                let pos = self.nodes[parent.get()].bodies[*slot].pos;
                self.put_at(child, pos, now, end);
                continue;
            }
            // The change happened at the parent's instant, and the child's
            // motion is carried to a later one. See `Tree::kick`.
            let at = self.nodes[parent.get()].time;
            self.kick(child, dv, at);
        }
    }

    /// Put a node where it is at `instant` — `pos` and `vel` in its parent's
    /// frame — and carry that to wherever its motion is, the way it is
    /// carried. A node whose motion is behind `instant` is brought up to it.
    pub fn put_at(&mut self, i: NodeIdx, pos: Vec3, vel: Vec3, instant: f64) {
        if !pos.is_finite() || !vel.is_finite() {
            return;
        }
        let about = self.turning_about(i, instant);
        let n = &mut self.nodes[i.get()];
        let since = n.carried - instant;
        if !(since > 0.0) || !since.is_finite() {
            n.motion.offset = pos;
            n.motion.velocity = vel;
            n.carried = n.carried.max(instant);
            return;
        }
        if n.turning == Vec3::ZERO {
            n.motion.offset = pos + vel.scale(since);
            n.motion.velocity = vel;
        } else {
            let (at, v) = turning_carry_about(n.turning, about, pos, vel, since);
            n.motion.offset = at;
            n.motion.velocity = v;
        }
    }

    /// The velocities of the bodies standing in for promoted children, so the
    /// change across a solve can be measured.
    pub fn stand_in_velocities(&self, parent: NodeIdx, slots: &[usize]) -> Vec<(usize, crate::math::Vec3)> {
        let n = &self.nodes[parent.get()];
        slots
            .iter()
            .filter_map(|s| n.bodies.get(*s).map(|b| (*s, b.vel)))
            .collect()
    }

    /// Re-derive a node's tier from the size it is *now*, and its refinement
    /// policy with it. Returns the tier it left, if it moved.
    ///
    /// Called wherever a node's size changes — `plant`, `emplace`, the growth
    /// step, a severing, and an authored radius. Deliberately *not* called from
    /// `coarsen`: there the radius is being restored rather than changed, and a
    /// node's tier, spec and identity are preserved across a round trip on
    /// purpose.
    ///
    /// # Why this is not simply done everywhere, every frame
    ///
    /// A tier is a physics regime. Changing one changes the node's solver, its
    /// timestep and its scheduling cadence, so a node re-tiered continuously
    /// could chatter across a boundary and take its solver with it. Tying it to
    /// the operations that change size means it moves when something *happened*
    /// rather than whenever a radius drifts, and those operations are exactly
    /// the ones the tier was wrong after.
    pub fn retier(&mut self, i: NodeIdx) -> Option<Tier> {
        if i.is_none() || !self.nodes[i.get()].alive {
            return None;
        }
        let parent_tier = {
            let p = self.nodes[i.get()].parent;
            if p.is_none() {
                Tier::Galactic
            } else {
                self.nodes[p.get()].tier
            }
        };
        let n = &self.nodes[i.get()];
        let want = tier_for(n.matter.radius, parent_tier);
        if want == n.tier {
            return None;
        }
        let was = n.tier;
        let spec = spec_for(want, n.spec);
        let n = &mut self.nodes[i.get()];
        n.tier = want;
        n.spec = spec;
        self.stats.retiers += 1;
        Some(was)
    }

    /// The gravitational field a node actually sits in, in its own frame.
    ///
    /// Derived, from the masses and radii of the things it is inside. Nothing
    /// here knows what a planet is: it walks the ancestor chain and adds what
    /// each one pulls with, which gives a surface `g` on a rocky planet, a
    /// different one on a moon, and very nearly nothing in interstellar space,
    /// for the same reason and by the same arithmetic.
    ///
    /// It replaces `solvers::structure::G_EARTH` on the path that decides how
    /// debris falls — `docs/PLAY.md` D6, and a `BACKLOG.md` entry that read
    /// "debris therefore falls at Earth gravity along its own structure's
    /// negative z wherever the node actually is — on a ship under thrust, in
    /// orbit, on a body of any other mass. The engine computes real
    /// gravitational fields at every other tier and then ignores them here."
    ///
    /// # The shell term, which is not a guard
    ///
    /// An ancestor pulls with its whole mass only from outside it. A node
    /// *within* one feels the mass enclosed below it, which for a uniform
    /// sphere is `M (d/R)^3` and gives `g = G M d / R^3` — linear in `d`, and
    /// zero at the centre. That is Newton's shell theorem rather than a
    /// singularity guard, and it is why this needs no epsilon: the field falls
    /// to nothing where the naive inverse square would blow up.
    ///
    /// An ancestor's mass includes this node's own, because a promoted child's
    /// stand-in body stays in its parent's list. It is subtracted: a thing does
    /// not pull on itself, and for a node that is most of what contains it the
    /// difference is the whole answer.
    ///
    /// # The answer is in the node's own axes
    ///
    /// Each ancestor pulls along a direction expressed in *that ancestor's*
    /// axes, because [`Tree::offset_from`] walks positions and positions are
    /// offsets in the parent's frame. Summing those directly adds vectors from
    /// different frames, and the sum is only meaningful while every frame in
    /// the chain is the same one — which was true for as long as nothing
    /// oriented a node against the body it sits on.
    ///
    /// A patch of ground is oriented by definition: that is what a patch *is*,
    /// a square of surface with a local up. So each contribution is rotated
    /// into the node's own axes on the way in, with [`Tree::axes_from`], and a
    /// node on the `+x` side of a planet carries its weight along its own `-z`
    /// rather than being told the field points along `-x`. For an unoriented
    /// chain every rotation is the identity and the arithmetic is unchanged,
    /// which is why this closed a `BACKLOG.md` entry without moving a number in
    /// the suite.
    pub fn gravity_at(&self, idx: NodeIdx) -> Vec3 {
        if idx.is_none() || !self.nodes[idx.get()].alive {
            return Vec3::ZERO;
        }
        self.facing(idx).conjugate().rotate(self.gravity_at_point(idx, Vec3::ZERO))
    }

    /// The field at a point `local` from a node's centre, from everything the
    /// node is inside, in root-aligned axes — [`Tree::gravity_at`] anywhere in
    /// the node rather than at its middle. A patch of a planet a continent
    /// across is pulled harder at its deep side than its shallow one: a face of
    /// an Earth has its centre of mass 1.5x10^6 m below the air over it, and the
    /// same law gives 15.4 m/s^2 there against 8.9 at the air.
    pub fn gravity_at_point(&self, idx: NodeIdx, local: Vec3) -> Vec3 {
        self.field_at(idx, local, false)
    }

    /// The field at a point inside a node, the node's own mass included as
    /// part of what it is a piece of — the smooth field a loose thing inside a
    /// patch of ground stands in, rather than the pull of the patch's pieces as
    /// points. A ball includes itself by its own interior law.
    pub fn field_within(&self, idx: NodeIdx, local: Vec3) -> Vec3 {
        let mut g = self.field_at(idx, local, true);
        let n = &self.nodes[idx.get()];
        if n.parent.is_none() || matches!(n.morphology.as_ref().and_then(|m| m.recipe.as_ref()), Some(crate::recipe::Recipe::Tiled(t)) if t.is_ball()) {
            let d = local.norm();
            let (m, radius) = (n.matter.mass.max(0.0), n.matter.radius);
            if d > 0.0 && radius > 0.0 {
                let enclosed = if d >= radius { m } else { m * (d / radius).powi(3) };
                g += local.scale(-crate::units::G * enclosed / (d * d * d));
            }
        }
        g
    }

    fn field_at(&self, idx: NodeIdx, local: Vec3, whole: bool) -> Vec3 {
        if idx.is_none() || !self.nodes[idx.get()].alive {
            return Vec3::ZERO;
        }
        let mut g = Vec3::ZERO;
        let mut inner = idx;
        let mut anc = self.nodes[idx.get()].parent;
        let mut first = true;
        while !anc.is_none() {
            let a = &self.nodes[anc.get()];
            let r = self.offset_from(anc, idx, local).value;
            let d = r.norm();
            let (m, radius) = (a.matter.mass.max(0.0), a.matter.radius);
            // **What this ancestor holds beyond the one below it.** An
            // ancestor's mass includes every descendant's, so counting each
            // ancestor's whole mass counts the chain once per level. The right
            // decomposition is nested: the galaxy contributes what it holds
            // that the cloud does not, the cloud what it holds that the parcel
            // does not, and so on, which is the shell theorem applied at every
            // level rather than only at the last one.
            //
            // It did not bite while the ladder was clouds inside clouds, where
            // a child is a thousandth of its parent and the double count is a
            // rounding. It bites the moment the ladder is a *surface*: a face
            // of a planet is a tenth of it, standing at nine tenths of its
            // radius, and counting it twice put 0.65 m/s^2 of sideways pull on
            // everything standing on it.
            // Its own mass left out, unless the point is being measured as
            // standing in the whole of the piece it is in.
            let own = if whole && first { 0.0 } else { self.nodes[inner.get()].matter.mass.max(0.0) };
            first = false;
            // **A piece of a sphere is not a sphere.** A patch of ground holds
            // its mass in a curved shell about the *planet's* centre, not in a
            // ball about its own, and the shell theorem is what says what that
            // field is: everything above a point contributes nothing and
            // everything below it contributes as if it were at the centre. So
            // for an ancestor that is a piece of a surface the radius vector is
            // measured from the body it is a piece of.
            //
            // The difference is not small and it is not only a magnitude.
            // Modelling a face of a planet as a ball at its own centre of mass
            // — 1.9 million metres underground — left a standing thing at
            // 8.96 m/s^2 tilted 2.6 degrees off its own down, because a ninth
            // of the planet was pulling sideways. With the shell the levels
            // telescope: each contributes what it holds beyond the one below,
            // all along the same radius, and the sum is the field the planet's
            // whole mass makes.
            let shell = self.nodes[anc.get()].morphology.as_ref().and_then(|mo| {
                let Some(crate::recipe::Recipe::Tiled(t)) = mo.recipe.as_ref() else {
                    return None;
                };
                if t.is_ball() {
                    return None;
                }
                Some(t.centre_of_mass_from_planet())
            });
            let (r, d, radius) = match shell {
                Some(com) => {
                    let from_centre = r + com;
                    (from_centre, from_centre.norm(), self.nodes[anc.get()].matter.radius)
                }
                None => (r, d, radius),
            };
            let _ = radius;
            if d > 0.0 && m > own {
                let source = m - own;
                // A point below the shell feels nothing of it; a point above
                // it feels all of it. For a ball the same expression is the
                // ordinary interior law, because its "shell" is itself.
                let enclosed = match shell {
                    Some(_) => source,
                    None => {
                        let radius = self.nodes[anc.get()].matter.radius;
                        if radius <= 0.0 {
                            continue;
                        }
                        if d >= radius {
                            source
                        } else {
                            source * (d / radius).powi(3)
                        }
                    }
                };
                // **Turned by the node's own absolute facing, not by its
                // facing relative to this ancestor.** An offset is a position
                // in root-aligned axes — nothing in the tree rotates one on the
                // way up — so a pull derived from offsets is in those axes too,
                // and what takes it into the node's frame is the whole
                // composition from the node to the root. Rotating by the
                // relative part alone is exact only while every frame between
                // here and the root is the same one, which is true until a
                // surface exists and false immediately afterwards.
                let pull = r.scale(-crate::units::G * enclosed / (d * d * d));
                g += pull;
            }
            inner = anc;
            anc = a.parent;
        }
        g
    }

    /// The rotation taking a vector in `node`'s axes to `ancestor`'s.
    ///
    /// `Motion::orientation` is which way a node is pointing *in its parent's
    /// frame*, and `Motion::compose` has always composed it correctly for one
    /// step. What was missing is the walk: nothing in the tree composed
    /// orientation across more than one level, so every vector carried between
    /// two frames — a field, a velocity, a wind direction, an impulse — was
    /// implicitly assuming every frame in the chain shared its axes.
    ///
    /// **Positions are deliberately not rotated by this.** A node's `offset` is
    /// a dynamical position in its parent's frame, and the solvers integrate it
    /// there; a node's `orientation` says which way the node itself is facing.
    /// The two are independent, and conflating them would drag every child of a
    /// spinning node around with it without any of the fictitious forces that
    /// would make that an honest rotating frame. `offset_from` is therefore
    /// unchanged, and this is the companion for everything that is not a
    /// position.
    pub fn axes_from(&self, ancestor: NodeIdx, mut node: NodeIdx) -> crate::math::Quat {
        // `a.then(b)` is the Hamilton product `a b`, which applies `b` first,
        // so each ancestor's rotation goes on the left of what is below it.
        // Written the other way round, a node's facing was applied after its
        // parent's — right only while the two turn about the same axis.
        let mut q = crate::math::Quat::IDENTITY;
        while node != ancestor && !node.is_none() {
            let n = &self.nodes[node.get()];
            q = n.motion.orientation.then(q);
            node = n.parent;
        }
        q
    }

    /// Which way a node is facing, as the rotation taking a vector fixed in
    /// its body to the root-aligned axes every position is written in: its
    /// own facing after its parent's, all the way up, the root's included.
    ///
    /// **A node's facing and its turning are offsets from its parent's** —
    /// `Motion::compose`'s reading, and the owner's decision for Phase 5.
    /// Positions and velocities are not: they stay root-aligned at every level
    /// (`axes_from`), and this is only for what is drawn *on* a body — an
    /// ocean's cells, a planet's tiles, which way a structure's layout faces.
    ///
    /// The root's own facing is part of it. `axes_from` stops short of the
    /// root, and an Earth at the root of its own world had an ocean standing
    /// still while its ground went round: measured, air held over the equator
    /// swept 428 of 1536 cells in a day and raised 0.95 m of sea on the far
    /// side of the planet.
    pub fn facing(&self, node: NodeIdx) -> crate::math::Quat {
        let n = &self.nodes[node.get()];
        if n.parent.is_none() {
            n.motion.orientation
        } else {
            self.facing(n.parent).then(n.motion.orientation)
        }
    }

    /// How fast a node is turning, rad/s, in root-aligned axes: its parent's
    /// turning and its own offset from it (`Tree::facing`).
    ///
    /// A face of a turning planet has no turning of its own — it is a piece of
    /// the planet — and so turns with it exactly, whatever shape it is. Read
    /// as its angular momentum over a uniform sphere of its radius instead, a
    /// face promoted from its planet turned at 23x the planet's rate, and the
    /// pieces drawn inside it at 0.63x.
    pub fn angular_velocity(&self, node: NodeIdx) -> Vec3 {
        let n = &self.nodes[node.get()];
        if n.parent.is_none() {
            n.motion.spin_rate
        } else {
            self.angular_velocity(n.parent) + self.facing(n.parent).rotate(n.motion.spin_rate)
        }
    }

    /// Angular momentum `dl` has arrived at a node, kg m^2/s in root-aligned
    /// axes — a contact, a solve, its contents folded back — and turns it by
    /// `dl / I` more relative to its parent.
    ///
    /// **`I` is the uniform sphere of the node's radius**, which is right for
    /// something round and an approximation for anything else: a slab takes
    /// a torque about its normal more easily than this says. It is used here
    /// and only here, for the change, and never to say what the whole turning
    /// is — which is what put a planet's face at 23x its rate. What it does not
    /// see either is a change of `I` itself: a node whose radius grows while it
    /// turns keeps its rate rather than slowing, until something measures it.
    pub fn turn_by(&mut self, node: NodeIdx, dl: Vec3) {
        if dl == Vec3::ZERO || !dl.is_finite() {
            return;
        }
        let i = self.nodes[node.get()].matter.moment_of_inertia();
        if !(i > 0.0) {
            return;
        }
        let parent = self.nodes[node.get()].parent;
        let dw = dl.scale(1.0 / i);
        let dw = if parent.is_none() { dw } else { self.facing(parent).conjugate().rotate(dw) };
        self.nodes[node.get()].motion.spin_rate += dw;
    }

    /// The turn a structure's drawing is placed at: the facing of the thing
    /// its layout is stated against. A tiled recipe states every level's cells
    /// in its *planet's* axes (`Tiled::render_on_sphere`), so a patch is drawn
    /// turned by its planet's facing; anything else is drawn in its own.
    ///
    /// Drawn unturned, a face of a turning Earth that nothing pins and is
    /// redrawn each frame kept its layout where it was at the start: 0.0000 rad
    /// turned in six hours against the ground's 1.57.
    pub fn drawn_turn(&self, node: NodeIdx) -> crate::math::Quat {
        self.drawn_turn_at(node, self.nodes[node.get()].time)
    }

    /// [`Tree::drawn_turn`] at an instant: the facing of what the layout is
    /// stated against, turned back from the instant that is carried to by its
    /// own turning. **At the instant asked for, not at whatever instant the
    /// facing's owner has reached**: a planet solved earlier in the same frame
    /// is carried to the frame's end before a patch of it is solved from the
    /// frame's start, and read there, a face's supports and anchors turned a
    /// minute ahead of its pieces every time its Earth was solved — 1.5 km at
    /// the corners, and 1.77 m/s^2 of pull on them.
    pub fn drawn_turn_at(&self, node: NodeIdx, instant: f64) -> crate::math::Quat {
        let mut at = node;
        loop {
            let n = &self.nodes[at.get()];
            match n.morphology.as_ref().and_then(|m| m.recipe.as_ref()) {
                Some(crate::recipe::Recipe::Tiled(t)) if t.on_sphere() && !t.is_ball() && !n.parent.is_none() => {
                    at = n.parent;
                }
                _ => break,
            }
        }
        let since = instant - self.nodes[at.get()].carried;
        if since == 0.0 || !since.is_finite() {
            return self.facing(at);
        }
        crate::math::Quat::from_rate(self.angular_velocity(at), since).then(self.facing(at)).unit()
    }

    /// Carry a vector expressed in `ancestor`'s axes into `node`'s own.
    ///
    /// The inverse direction of [`Tree::axes_from`], which is the one nearly
    /// every consumer wants: a field, a flow or an impulse is measured where it
    /// comes from and has to be applied where it lands.
    pub fn into_axes_of(&self, ancestor: NodeIdx, node: NodeIdx, v: Vec3) -> Vec3 {
        let q = self.axes_from(ancestor, node);
        if q == crate::math::Quat::IDENTITY {
            return v;
        }
        q.conjugate().rotate(v)
    }

    /// The region this node owns — `docs/PLAY.md` D16, and the length every
    /// boundary crossing is measured against.
    ///
    /// The larger of the volume its matter occupies and the distance at which
    /// its parent's gravity takes over. See `crossing`'s module documentation
    /// for why it is not simply the radius, and for the measurements that
    /// settled it; the short version is that a surface is exactly where a
    /// sphere's boundary is, so a geometric rule alone ejects anything standing
    /// on a planet into interplanetary space.
    ///
    /// The root owns its own radius and nothing more, which costs nothing: the
    /// universe has no outside to be re-homed into, so no crossing is ever
    /// measured against it.
    pub fn domain(&self, i: NodeIdx) -> f64 {
        if i.is_none() || !self.nodes[i.get()].alive {
            return 0.0;
        }
        let n = &self.nodes[i.get()];
        let own = n.matter.radius.max(0.0);
        if n.parent.is_none() {
            return own;
        }
        let parent_mass = self.nodes[n.parent.get()].matter.mass;
        own.max(crate::crossing::hill_radius(
            n.matter.mass,
            parent_mass,
            n.motion.offset.norm(),
        ))
    }

    /// How far this node's contents reach from its centre, ignoring one of
    /// them.
    ///
    /// The third term of what a node owns, and the one that stops a crossing
    /// firing on a distribution's own tail. `matter.radius` is the *equivalent
    /// uniform sphere*, so a centrally-concentrated draw legitimately puts
    /// bodies well outside it — `docs/BACKLOG.md` measures a Plummer sphere's
    /// tail at three to four radii and the scenario shelf at 1.5 to 4.0.
    /// Measured on the ladder `drill_to` builds: seven promoted parcels sat at
    /// 1.0 to 3.6 of their parent's radius and at **0.28 to 0.93 of what the
    /// parent's other contents reach**. They had not left anything. A rule that
    /// re-homed them would flatten every ladder in the engine, and did.
    ///
    /// So the question is not "is it outside the sphere the node claims" but
    /// "is it beyond everything else the node holds", which is a measurement
    /// rather than a radius. `ignoring` is what keeps it from being circular:
    /// the escapee is usually the furthest thing there is, and a boundary that
    /// chases it can never be crossed.
    ///
    /// Distance from the node's *origin* rather than from its contents' centre
    /// of mass, because that is the frame a child's offset is expressed in.
    /// O(n) in the node's contents, and [`crate::engine::World::cross`] pays it
    /// only for a node the cheap bound has already called escaped.
    pub fn contents_reach(&self, i: NodeIdx, ignoring: NodeIdx) -> f64 {
        let (furthest, whose, second) = self.reach_pair(i);
        if !ignoring.is_none() && ignoring == whose {
            second
        } else {
            furthest
        }
    }

    /// How far the two furthest of a node's contents reach, and which occupant
    /// the furthest is.
    ///
    /// The form [`Self::contents_reach`] is really asking for, and the reason
    /// it is separate: only one occupant's exclusion can change the answer, so
    /// a caller testing every child of one parent can measure the parent *once*
    /// instead of once per child. Measured before it did: the crossing pass
    /// cost 144 ms on a world of 8193 nodes, because every child of a galaxy
    /// sitting in its tail paid a full sweep of the galaxy's contents. It is
    /// linear again with this.
    ///
    /// The second value is the furthest reach *excluding* the furthest
    /// occupant, which is all an exclusion can ever need.
    pub fn reach_pair(&self, i: NodeIdx) -> (f64, NodeIdx, f64) {
        if i.is_none() || !self.nodes[i.get()].alive {
            return (0.0, NodeIdx::NONE, 0.0);
        }
        let n = &self.nodes[i.get()];
        let slots = n.bodies.len().max(n.children.len());
        let (mut furthest, mut whose, mut second) = (0.0f64, NodeIdx::NONE, 0.0f64);
        for slot in 0..slots {
            let child = n.child_of(slot);
            let reach = if !child.is_none()
                && self.nodes.get(child.get()).is_some_and(|c| c.alive)
            {
                let c = &self.nodes[child.get()];
                c.motion.offset.norm() + c.matter.radius.max(0.0)
            } else if let Some(b) = n.bodies.get(slot) {
                b.pos.norm() + b.radius.max(0.0)
            } else {
                continue;
            };
            if reach > furthest {
                second = furthest;
                furthest = reach;
                whose = child;
            } else if reach > second {
                second = reach;
            }
        }
        (furthest, whose, second)
    }

    /// The same question asked of one of a node's *occupants*, which may be a
    /// promoted child or may still be a body.
    ///
    /// One rule for both, because `Neighbourhood` already treats them as one
    /// index: a promoted child *is* one of its parent's bodies seen one level
    /// down, and a crossing that could only land in an already-promoted sibling
    /// would be waiting for somebody else to have visited the place first.
    /// Returns where the occupant is and how far its claim reaches.
    pub fn occupant_claim(&self, parent: NodeIdx, slot: usize) -> Option<(Vec3, f64)> {
        let n = self.nodes.get(parent.get())?;
        let child = n.child_of(slot);
        if !child.is_none() {
            let c = self.nodes.get(child.get())?;
            if c.alive {
                return Some((c.motion.offset, self.domain(child)));
            }
        }
        let b = n.bodies.get(slot)?;
        let claim = b.radius.max(0.0).max(crate::crossing::hill_radius(
            b.mass,
            n.matter.mass,
            b.pos.norm(),
        ));
        Some((b.pos, claim))
    }

    /// How far this node's contents actually extend from their own centre.
    ///
    /// Bodies **and** promoted children, the same union
    /// [`crate::neighbourhood::Neighbourhood`] indexes, because a promoted
    /// child is contents at a finer resolution rather than a different kind of
    /// thing — and a child that has drifted out of its parent is exactly the
    /// case `docs/BACKLOG.md`'s fragment entry is about. A promoted child
    /// stands in for the body in its slot, so that slot is counted once.
    ///
    /// Returns an empty [`Spread`] for a node with no contents, which is
    /// honest: nothing has no extent.
    pub fn spread(&self, i: NodeIdx) -> crate::state::Spread {
        if i.is_none() || !self.nodes[i.get()].alive {
            return crate::state::Spread::default();
        }
        let n = &self.nodes[i.get()];
        let slots = n.bodies.len().max(n.children.len());
        let mut parts = Vec::with_capacity(slots);
        for slot in 0..slots {
            let promoted = n
                .children
                .get(slot)
                .copied()
                .filter(|c| !c.is_none())
                .and_then(|c| self.nodes.get(c.get()))
                .filter(|c| c.alive);
            match promoted {
                Some(c) => parts.push((c.motion.offset, c.matter.mass, c.matter.radius)),
                None => {
                    if let Some(b) = n.bodies.get(slot) {
                        parts.push((b.pos, b.mass, b.radius));
                    }
                }
            }
        }
        crate::state::Spread::of(parts)
    }

    /// Which of a node's contents are still **one neighbourhood**.
    ///
    /// `docs/PLAY.md` §7 Phase 3: the same measurement that drives a crossing
    /// drives node splitting, which Phase 1 measured and connected to nothing.
    /// A node is a region of space, and the claim a region makes is that what
    /// is in it is *near* what else is in it. When that stops being true the
    /// node is describing two places at once, and every length derived from its
    /// radius — the smoothing length, the gravity softening, the LOD's angular
    /// size, the neighbour grid's spacing — is wrong by the same factor.
    ///
    /// **Connected components under the node's own resolution**, and that is
    /// the whole criterion. No threshold is chosen: two occupants are together
    /// if they are adjacent in the sense D3's [`crate::neighbourhood`] already
    /// defines, and a component is the transitive closure of that. This is why
    /// it does not fire on a distribution's own tail — `docs/BACKLOG.md`
    /// measures a Plummer sphere's at three to four radii and warns that
    /// "outgrew its radius" is not the fault signal — because a tail is
    /// *connected* to the body it is the tail of, however far out it reaches.
    ///
    /// Returns a label per occupant and the number of distinct labels, in the
    /// neighbourhood's own occupant order. An unmaterialised node has no
    /// contents and therefore no components.
    pub fn components(&self, i: NodeIdx) -> (crate::neighbourhood::Neighbourhood, Vec<usize>, usize) {
        let nb = self.neighbourhood(i);
        let n = nb.len();
        let mut label: Vec<usize> = (0..n).collect();
        if n == 0 {
            return (nb, label, 0);
        }
        // Union-find, iterative path compression. The sets are tiny and the
        // ordering is fixed by `pairs`, so this is deterministic.
        fn find(label: &mut [usize], mut x: usize) -> usize {
            while label[x] != x {
                label[x] = label[label[x]];
                x = label[x];
            }
            x
        }
        // The widest query the index can answer, which is the node's own
        // resolution unless one oversized occupant has forced the spacing up.
        let within = nb.resolution().min(nb.reach());
        if let Some(pairs) = nb.pairs(within) {
            for (a, b) in pairs {
                let (ra, rb) = (find(&mut label, a), find(&mut label, b));
                if ra != rb {
                    label[ra.max(rb)] = ra.min(rb);
                }
            }
        }
        for x in 0..n {
            label[x] = find(&mut label, x);
        }
        let mut seen: Vec<usize> = label.clone();
        seen.sort_unstable();
        seen.dedup();
        (nb, label, seen.len())
    }

    /// What is next to what, inside this node.
    ///
    /// Built on demand rather than cached. Whether it should be cached is a
    /// question for a measurement rather than a guess: the build is O(n) over
    /// the node's contents, and until something asks for it often enough to
    /// matter, a cache is a second thing that can go stale. `Neighbourhood`
    /// carries the epoch it was built at so that a caller who does keep one can
    /// tell.
    ///
    /// Returns an empty neighbourhood for a node that is not materialised,
    /// which is honest: a node with no contents has nothing next to anything.
    pub fn neighbourhood(&self, i: crate::ids::NodeIdx) -> crate::neighbourhood::Neighbourhood {
        let n = &self.nodes[i.get()];
        let parts = n.bodies.len();
        // The node's own resolution, and deliberately the same expression
        // `World::node_resolution` and `Matter::signal_crossing` use. Three
        // definitions of one length is how they drift apart.
        let resolution = if parts > 1 {
            n.matter.radius / (parts as f64).cbrt()
        } else {
            n.matter.radius
        };
        crate::neighbourhood::Neighbourhood::build(
            &n.bodies,
            &n.children,
            resolution,
            n.epoch,
            |c| {
                let child = self.nodes.get(c.get())?;
                // Where it is at this node's instant: see `Node::carried`.
                child.alive.then(|| (self.position_at(c, n.time), child.matter.radius))
            },
        )
    }

    /// Give a node a developmental state, turning it from a statistical
    /// population into a structure with a history.
    ///
    /// The matter's radius, chemical energy and entropy are taken over by
    /// the morphology from this point on; the conserved tuple is untouched, so
    /// nothing about the surrounding world changes.
    pub fn plant(&mut self, i: NodeIdx, program: crate::morph::Program) -> &mut crate::morph::Morphology {
        let key = self.nodes[i.get()].key;
        let seed = self.world_seed;
        let epoch = self.nodes[i.get()].epoch;
        let mut m = crate::morph::Morphology::new(program, seed, key.0, epoch);
        // A seed, not a finished structure. The node's remaining mass is the
        // feedstock the thing grows out of — soil, air, water — so planting
        // neither creates nor destroys anything, and growth is bounded by what
        // is actually there.
        m.built = (self.nodes[i.get()].matter.mass * 1e-3).clamp(1e-6, 1.0);
        // Write down the rule. `Tree` has no registry and no observers, so it
        // generates against the conditions it can see; `World::plant` measures
        // the real ones and writes it again. Doing it here as well is what
        // keeps a `Tree` on its own — which is what most of the suite is — a
        // thing that can draw itself.
        m.regenerate(&crate::morph::Environment::default(), &m.program.material());
        // Whatever was promoted out of this node is folded back before its
        // detail goes: a slot in a body list that is about to be replaced is
        // not a place anything can live. See `Tree::shed_children`.
        self.shed_children(i);
        let n = &mut self.nodes[i.get()];
        n.matter.radius = m.extent().max(n.matter.radius.min(1e-3)).max(1e-30);
        n.matter.chemical_energy = m.stored_energy();
        n.bodies.clear();
        n.morphology = Some(m);
        self.stats.structures += 1;
        // A structure takes its size from its program the instant it has one,
        // and a seed is not the size of the tree it becomes.
        self.retier(i);
        self.nodes[i.get()].morphology.as_mut().unwrap()
    }

    /// Give a node a structure that is *already there*, at a stated mass.
    ///
    /// [`Self::plant`] seeds: it gives the node a thousandth of its mass and
    /// lets the program build the rest, which is what a tree or a coral does
    /// and what a construction site does. Terrain does neither. A hillside was
    /// not grown and nobody built it; it is simply the shape the ground is, and
    /// starting it as a one-kilogram seed that has to accumulate a mountain is
    /// a description of geology nobody wants to wait for.
    ///
    /// So this is the other verb: the structure exists, at this mass, complete.
    /// It is equally the right one for a town that was already standing when
    /// the player arrived — the distinction is not living against built, it is
    /// whether the world is watching it happen.
    pub fn emplace(
        &mut self,
        i: NodeIdx,
        program: crate::morph::Program,
        built: f64,
    ) -> &mut crate::morph::Morphology {
        let key = self.nodes[i.get()].key;
        let seed = self.world_seed;
        let epoch = self.nodes[i.get()].epoch;
        let mut m = crate::morph::Morphology::new(program, seed, key.0, epoch);
        // Bounded by what is actually in the node, for the same reason planting
        // is: a structure cannot be made of more than the matter available.
        m.built = built.clamp(0.0, self.nodes[i.get()].matter.mass);
        if program.is_planned() {
            // A planned program reads `progress`, not `built`, to decide how
            // much of itself to draw. Finished means finished.
            m.design_mass = m.built;
            m.progress = 1.0;
        }
        // The rule, written against the conditions a `Tree` can see. See
        // `Tree::plant`.
        m.regenerate(&crate::morph::Environment::default(), &m.program.material());
        self.shed_children(i);
        let n = &mut self.nodes[i.get()];
        n.matter.radius = m.extent().max(1e-30);
        n.matter.chemical_energy = m.stored_energy();
        n.bodies.clear();
        n.morphology = Some(m);
        self.stats.structures += 1;
        // A structure takes its size from its program the instant it has one,
        // and a seed is not the size of the tree it becomes.
        self.retier(i);
        self.nodes[i.get()].morphology.as_mut().unwrap()
    }

    /// State a composite: one node whose recipe is the parts it is made of.
    ///
    /// The third verb beside `plant` and `emplace`, and `docs/PLAY.md` D15's
    /// forward direction. `plant` seeds something that grows; `emplace` states
    /// something already grown; `assemble` states something that was *made* —
    /// six planks attached into a box — and the difference is where the
    /// geometry comes from. A grown thing derives its size from the mass it
    /// accumulated. An assembled thing is the size its parts are, and its mass
    /// and its radius both follow from them.
    ///
    /// `program` is provenance rather than species: it says what the parts are
    /// made of, so a box of oak planks weighs and burns like oak, and no
    /// `Program` variant has to exist for "box". That is D11 honoured rather
    /// than worked around — a crate and a cathedral are both parts lists.
    ///
    /// The node's mass is **not** overwritten. A box of six 40 kg planks in a
    /// node holding 300 kg is a box with 60 kg of something else in it, which
    /// is a sampled remainder exactly as litter under a tree is; `sample_structured`
    /// already splits the two. A node holding *less* than its parts claim is a
    /// caller error and the parts are scaled down to fit rather than
    /// conjuring mass.
    pub fn assemble(
        &mut self,
        i: NodeIdx,
        program: crate::morph::Program,
        parts: crate::assembly::Assembly,
    ) -> &mut crate::morph::Morphology {
        let key = self.nodes[i.get()].key;
        let seed = self.world_seed;
        let available = self.nodes[i.get()].matter.mass;
        let mut parts = parts;
        let claimed = parts.mass();
        if claimed > available && available > 0.0 && claimed > 0.0 {
            let f = available / claimed;
            for p in parts.parts.iter_mut() {
                p.mass *= f;
            }
        }
        let m = crate::morph::Morphology::assembled(program, parts, seed, key.0);
        self.shed_children(i);
        let n = &mut self.nodes[i.get()];
        n.matter.radius = m.extent().max(1e-30);
        n.matter.chemical_energy = m.stored_energy();
        n.bodies.clear();
        n.topology = None;
        n.morphology = Some(m);
        self.stats.structures += 1;
        self.retier(i);
        self.nodes[i.get()].morphology.as_mut().unwrap()
    }

    /// Record that a blend could not describe everything it was given.
    fn note_description_lost(&mut self, lost: f64) {
        if lost > 0.0 {
            self.stats.over_described += 1;
            self.stats.worst_description_lost = self.stats.worst_description_lost.max(lost);
        }
    }

    /// Mark a node as holding detail nothing can regenerate, and its ancestry
    /// as containing an edit.
    ///
    /// **The two are different claims and `docs/PLAY.md` D19 separates them.**
    /// This used to set `pinned` on the whole ancestry, which is one-way and
    /// walks to the root, and the scheduler refuses to coarsen a pinned node at
    /// all — so felling one tree made a whole planet permanently
    /// un-coarsenable, then its star, then its galaxy. Measured: pinning one
    /// leaf of a six-deep ladder left six nodes that could never be collapsed
    /// again, against one now.
    ///
    /// An ancestor that merely *contains* an edit can still collapse, because
    /// its recipe plus its descendants' edits is all it needs. What it may not
    /// do is forget — see [`Node::contains_edit`].
    pub fn pin(&mut self, i: NodeIdx) {
        if i.is_none() {
            return;
        }
        {
            let n = &mut self.nodes[i.get()];
            n.pinned = true;
            n.residency = Residency::Pinned;
        }
        self.note_edit_above(i);
    }

    /// Record a change the recipe *can* express: no pinning anywhere.
    ///
    /// `docs/PLAY.md` D19's other half. A box that lost a wall is a box with a
    /// five-wall recipe and a break, and that description regenerates — so it
    /// must not fall back to a stored body list, or the phase's own done-when
    /// ("returns to ~100 bytes when nobody is watching") is unreachable.
    ///
    /// The change itself lives where the recipe lives: `Morphology::events` for
    /// a grown or built thing, which already exists and is already replayed
    /// rather than discarded. This is the *marking* — the node and its ancestry
    /// stop being candidates for a fresh draw from the ensemble, and stay
    /// candidates for collapsing to a description.
    pub fn record_edit(&mut self, i: NodeIdx) {
        if i.is_none() {
            return;
        }
        self.nodes[i.get()].contains_edit = true;
        self.note_edit_above(i);
    }

    /// Walk the ancestry marking each node as containing an edit.
    fn note_edit_above(&mut self, i: NodeIdx) {
        let mut cur = self.nodes[i.get()].parent;
        while !cur.is_none() {
            let n = &mut self.nodes[cur.get()];
            if n.contains_edit {
                // Already marked, and so is everything above it: the walk is
                // idempotent and the flag is one-way, so there is nothing
                // further up that this call could add.
                break;
            }
            n.contains_edit = true;
            cur = n.parent;
        }
    }

    /// Reconcile a node with what it is actually holding.
    ///
    /// `docs/BACKLOG.md`'s three outcomes, and the entry that has been waiting
    /// for a caller since Phase 1 built its measurement: *within the radius, do
    /// nothing; larger than the radius but still one clump, grow the radius and
    /// re-derive everything that depends on it; genuinely bimodal, split.*
    ///
    /// The middle outcome is [`Self::settle`] — the radius is measured from the
    /// contents by `summarise`, which is the same re-derivation `coarsen` does,
    /// followed by [`Self::retier`] because a node that changed size may not be
    /// the size of thing it was. There is deliberately no separate "grow the
    /// radius" path: a radius set by anything other than a measurement is a
    /// number nothing can check.
    ///
    /// Returns the node that was split off, if one was.
    pub fn resolve_extent(&mut self, i: NodeIdx) -> Option<NodeIdx> {
        if i.is_none() || !self.nodes[i.get()].alive || !self.nodes[i.get()].is_materialised() {
            return None;
        }
        let (nb, label, count) = self.components(i);
        if count <= 1 {
            // One neighbourhood, however far it reaches, so the only question
            // left is whether the radius still describes it.
            self.follow_contents(i);
            return None;
        }
        // Which component keeps the node. The heaviest, because the node's
        // identity, its address and its pinned detail all stay with it, and
        // moving the bulk of the mass to a new address for the sake of a
        // fragment is the expensive way round.
        let children = self.nodes[i.get()].children.clone();
        let mut mass: std::collections::BTreeMap<usize, f64> = Default::default();
        let mut plain: std::collections::BTreeMap<usize, usize> = Default::default();
        for (k, l) in label.iter().enumerate() {
            let Some((occ, _, _)) = nb.at(k) else { continue };
            let Some(slot) = occ.slot(&children) else { continue };
            let m = match occ {
                crate::neighbourhood::Occupant::Child(c) => self.nodes[c.get()].matter.mass,
                crate::neighbourhood::Occupant::Body(_) => {
                    *plain.entry(*l).or_insert(0) += 1;
                    self.nodes[i.get()].bodies.get(slot).map(|b| b.mass).unwrap_or(0.0)
                }
            };
            *mass.entry(*l).or_insert(0.0) += m.max(0.0);
        }
        let total: f64 = mass.values().sum();
        let keeps = mass.iter().max_by(|a, b| a.1.total_cmp(b.1)).map(|(l, _)| *l)?;

        // **A node describes one region when one component holds the bulk of
        // it**, and the rest are that region's tail. This is the gate, and it
        // is here because connectedness on its own does not survive contact
        // with a real draw: the linking length is the mean spacing, a centrally
        // concentrated profile's outskirts are sparser than the mean, and so
        // isolated bodies fall out of the main component on every node in the
        // engine. Measured, all on worlds nobody had touched —
        //
        //   a rocky planet, 64 bodies     9 components, largest 86.5%
        //   a rocky planet, 512 bodies   84 components, largest 77.8%
        //   a planetary node, 4000        31 components, largest 99.05%
        //   a cloud pulled into two       20 components, largest 48.4%
        //
        // — and the separations do not tell them apart either: the planet's
        // stray components sit 0.9 to 2.0 radii out and the genuinely bimodal
        // pair sits at 1.79. A mass-share floor does not either; the planet's
        // tails are 4% and the bimodal pair's own noise is 0.39%. What
        // separates them by a factor of twenty is whether anything holds a
        // majority: 86.5% and 77.8% against 48.4%.
        //
        // **The known miss is a 60/40 separation**, which stays one node until
        // the shares even out or it drifts far enough to cross. That is
        // recorded on the plan rather than here, because catching it needs a
        // ratio the plan does not state.
        if total > 0.0 && mass[&keeps] / total > 0.5 {
            self.follow_contents(i);
            return None;
        }
        // **And the same question asked of the pair.** No majority does not by
        // itself mean two places: it also describes a node that has come apart
        // into *many*, where splitting one piece off achieves nothing and would
        // do it again next frame. Measured, on a draw whose bodies were spread
        // until none of them reached its neighbours: 4000 components of one
        // body each, no majority anywhere, and a split for every frame for
        // ever. That node is **dispersed**, and `docs/BACKLOG.md`'s answer to
        // dispersal is the radius following its contents rather than a new
        // node.
        //
        // So the rule is the majority rule twice over: one component holding
        // the bulk is one region; two holding it between them are two regions;
        // and no small set holding it at all is one region that has spread.
        let second = mass
            .iter()
            .filter(|(l, _)| **l != keeps)
            .map(|(_, m)| *m)
            .fold(0.0f64, f64::max);
        if total > 0.0 && (mass[&keeps] + second) / total <= 0.5 {
            self.follow_contents(i);
            return None;
        }

        let mut candidates: Vec<usize> = mass
            .iter()
            .filter(|(l, _)| **l != keeps && plain.get(l).copied().unwrap_or(0) > 0)
            .map(|(l, _)| *l)
            .collect();
        candidates.sort_by(|a, b| mass[b].total_cmp(&mass[a]));
        for leaves in candidates {
            let slots: Vec<usize> = label
                .iter()
                .enumerate()
                .filter(|(_, l)| **l == leaves)
                .filter_map(|(k, _)| nb.at(k).and_then(|(occ, _, _)| occ.slot(&children)))
                .collect();
            if slots.is_empty() {
                continue;
            }
            if let Some(new) = self.split_off(i, &slots) {
                return Some(new);
            }
        }
        // Nothing could be made into a node of its own — a component of nothing
        // but promoted children, which are nodes already, and which D16's
        // crossing re-homes when they leave.
        self.follow_contents(i);
        None
    }

    /// Let the radius follow the contents — **where the contents are the
    /// authority**.
    ///
    /// `docs/BACKLOG.md`'s second outcome: contents that are still one clump
    /// but have outgrown the radius should grow it rather than be split. The
    /// question it leaves open is *whose* answer the radius is, and the engine
    /// already draws that line.
    ///
    /// For a node whose detail is **regenerable**, the matter is the authority
    /// and the bodies are a drawing of it — `sample` scales what it draws until
    /// `1.291 x rms` comes back as the radius it was given, so the radius is an
    /// *input* to the detail, not a measurement of it. Letting it follow the
    /// bodies inverts that, and the inversion feeds back. Measured, on the
    /// biome world: a star's eight sampled parcels disperse under the hydro
    /// solver, the radius follows them from 7x10^8 m to 2.9x10^11 in four
    /// hundred frames, the radiating area grows with it, the star cools from
    /// 5800 K to 719 K, and `retier` moves it out of `Stellar` — where
    /// temperature means a velocity dispersion — into `Planetary`, where it
    /// means heat. A test about whether a patch of ground freezes in winter
    /// failed because its sun had quietly turned into something else.
    ///
    /// For a node whose detail has been **touched** — pinned, or carrying an
    /// edit below it — the bodies are the authority, because that is what
    /// pinning means, and the radius has to keep up with them. A node that has
    /// just been split is pinned by [`Self::split_off`] for exactly this
    /// reason.
    ///
    /// The regenerable case is not silently dropped: `Stats::worst_occupancy`
    /// goes on reporting how far a node's contents have outgrown what it
    /// claims, which is the measurement `docs/BACKLOG.md` asks for and the one
    /// that says when the consumers of a node radius — the smoothing length,
    /// the softening, the LOD's angular size — are being lied to.
    fn follow_contents(&mut self, i: NodeIdx) {
        let n = &self.nodes[i.get()];
        if !(n.pinned || n.contains_edit) {
            return;
        }
        self.remeasure_radius(i);
        self.retier(i);
    }

    /// Put a node's radius back in step with the contents it is holding.
    ///
    /// The *equivalent uniform sphere* of the present configuration, which is
    /// the same expression `summarise` uses and must stay the same one: a
    /// second definition of this length is how `matter.radius` and what the
    /// sampler draws come apart. Deliberately **not** the bounding radius — see
    /// `sampler::radius_scale` — so a centrally concentrated node keeps a
    /// radius its tail reaches past, exactly as a freshly sampled one does.
    ///
    /// Written only when it moved, so a node whose contents are holding still
    /// is left bit-for-bit as it was.
    pub fn remeasure_radius(&mut self, i: NodeIdx) -> f64 {
        let n = &self.nodes[i.get()];
        let count = n.bodies.len();
        if count == 0 {
            return n.matter.radius;
        }
        let mass = crate::math::det_sum_by(count, &|k| n.bodies[k].mass);
        if !(mass > 0.0) {
            return n.matter.radius;
        }
        let com = crate::math::det_sum_v3_by(count, &|k| n.bodies[k].pos.scale(n.bodies[k].mass))
            .scale(1.0 / mass);
        let r2 = crate::math::det_sum_by(count, &|k| {
            n.bodies[k].mass * (n.bodies[k].pos - com).norm2()
        }) / mass;
        let radius = (r2.max(0.0).sqrt() * crate::state::RMS_TO_RADIUS).max(1e-30);
        if radius != n.matter.radius {
            self.nodes[i.get()].matter.radius = radius;
        }
        radius
    }

    /// Make a sibling out of a subset of a node's contents.
    ///
    /// The operation `docs/BACKLOG.md` calls "sibling-from-a-subset", which it
    /// names twice: node splitting is this, and so is the promoted fragment —
    /// "`promote` takes a single slot, while a fragment is a *set* of members".
    ///
    /// The new node is a child of the same parent, holding the departing
    /// bodies re-expressed about their own centre of mass, with its matter
    /// summarised from them. Promoted children in the subset are re-homed into
    /// it by [`Self::reparent`], which is the one path that moves a node
    /// between frames and the one that carries everything keyed by its address.
    ///
    /// **Both ends are pinned.** Neither is what `sample` would draw from its
    /// parent's matter any more — the old node's body list now has holes in it
    /// and the new one was never drawn at all — so neither is regenerable, and
    /// saying so is what stops a later refinement quietly mending the split.
    ///
    /// Refused for the root, which has no parent to be a sibling in, and for a
    /// subset that is everything or nothing, which would be a rename.
    pub fn split_off(&mut self, i: NodeIdx, slots: &[usize]) -> Option<NodeIdx> {
        if i.is_none() || !self.nodes[i.get()].alive || slots.is_empty() {
            return None;
        }
        let parent = self.nodes[i.get()].parent;
        if parent.is_none() {
            return None;
        }
        let occupied = self.nodes[i.get()]
            .bodies
            .iter()
            .enumerate()
            .filter(|(slot, b)| b.mass > 0.0 || !self.nodes[i.get()].child_of(*slot).is_none())
            .count();
        if slots.len() >= occupied {
            return None;
        }

        // What is leaving, measured before anything moves. A promoted child's
        // stand-in is current — `sync_children` runs at the head of every solve
        // — so the summary counts it exactly once, at the mass the child says.
        let leaving: Vec<Body> = slots
            .iter()
            .filter_map(|s| self.nodes[i.get()].bodies.get(*s).copied())
            .collect();
        if leaving.is_empty() {
            return None;
        }
        let matter = summarise(&leaving, 0.0);
        if !(matter.mass > 0.0) {
            return None;
        }
        let com = matter.com;
        let bulk = if matter.mass > 0.0 {
            matter.momentum.scale(1.0 / matter.mass)
        } else {
            Vec3::ZERO
        };

        // The departing detail, about its own centre and in its own frame.
        let mut bodies: Vec<Body> = Vec::with_capacity(leaving.len());
        let mut movers: Vec<NodeIdx> = Vec::new();
        for s in slots {
            let child = self.nodes[i.get()].child_of(*s);
            if !child.is_none() {
                movers.push(child);
                continue;
            }
            if let Some(b) = self.nodes[i.get()].bodies.get(*s) {
                let mut b = *b;
                b.pos -= com;
                b.vel = crate::coords::velocity_add(-bulk, b.vel);
                bodies.push(b);
            }
        }
        if bodies.is_empty() {
            return None;
        }

        let (offset, velocity, spec, time) = {
            let n = &self.nodes[i.get()];
            (
                n.motion.offset + com,
                crate::coords::velocity_add(n.motion.velocity, bulk),
                n.spec,
                n.time,
            )
        };
        let radius = matter.radius.max(1e-30);
        let tier = tier_for(radius, self.nodes[parent.get()].tier);
        let spec = spec_for(tier, spec);

        // Take a slot in the parent, exactly as `reparent` does and for the
        // same reasons: pushed rather than inserted, so no sibling is
        // renumbered and no address changes but this one.
        let new_slot = {
            let p = &mut self.nodes[parent.get()];
            p.bodies.push(Body {
                pos: offset,
                vel: velocity,
                mass: matter.mass,
                radius,
                charge: matter.charge,
                internal_energy: matter.internal_energy,
                spin: matter.spin,
                temperature: matter.temperature,
                composition: matter.composition,
                kind: leaving.first().map(|b| b.kind).unwrap_or(crate::state::BodyKind::Grain),
                ..Default::default()
            });
            while p.children.len() < p.bodies.len() {
                p.children.push(NodeIdx::NONE);
            }
            p.bodies.len() - 1
        };
        let new_key = self.nodes[parent.get()].key.child(new_slot as u64);
        let new_depth = self.nodes[parent.get()].depth + 1;
        let mut matter = matter;
        matter.mixture = self.nodes[i.get()].matter.mixture;
        matter.momentum = Vec3::ZERO;
        let sibling = Node {
            key: new_key,
            parent,
            slot: new_slot as u32,
            depth: new_depth,
            tier,
            matter,
            motion: Motion {
                offset,
                velocity,
                orientation: self.nodes[i.get()].motion.orientation,
                // A piece of the node it split from, turning as it did.
                spin_rate: self.nodes[i.get()].motion.spin_rate,
                proper_time: self.nodes[i.get()].motion.proper_time,
            },
            bodies: Vec::new(),
            potential: 0.0,
            gravity: Vec3::ZERO,
            children: Vec::new(),
            spec,
            epoch: 0,
            time,
            last_disturbed: time,
            last_solved: time,
            last_grown: time,
            residency: self.nodes[i.get()].residency,
            pinned: true,
            contains_edit: false,
            rest_density: self.nodes[i.get()].rest_density,
            unrest: 0.0,
            ocean: None,
            carried: time,
            turning: Vec3::ZERO,
            ground: None,
            bubble: self.nodes[i.get()].bubble,
            alive: true,
            morphology: None,
            topology: None,
            steps_taken: 0,
            surface: None,
            surface_epoch: u32::MAX,
            last_report: SampleReport::default(),
        };
        let new = self.alloc(sibling);
        self.nodes[parent.get()].children[new_slot] = new;
        {
            let n = &mut self.nodes[new.get()];
            n.children = vec![NodeIdx::NONE; bodies.len()];
            n.bodies = bodies;
        }

        // Vacate what left. Zeroed rather than removed, because a sibling's
        // address is derived from its slot index.
        for s in slots {
            if self.nodes[i.get()].child_of(*s).is_none() {
                if let Some(b) = self.nodes[i.get()].bodies.get_mut(*s) {
                    *b = Body::default();
                }
            }
        }
        // And anything that had a node of its own goes through the one path
        // that moves a node between frames.
        for c in movers {
            self.reparent(c, new);
        }

        // **Re-centre what is left.** A node's frame origin is where its
        // contents are, and after a piece leaves from one side they are not
        // there any more — which every cheap test in the engine then gets
        // wrong, starting with the overlap check in [`Self::would_merge`],
        // which compares declared centres. The positions in the *parent's*
        // frame do not move: the offset gains exactly what the contents lose.
        let com = {
            let n = &self.nodes[i.get()];
            let mass = crate::math::det_sum_by(n.bodies.len(), &|k| n.bodies[k].mass);
            if mass > 0.0 {
                crate::math::det_sum_v3_by(n.bodies.len(), &|k| {
                    n.bodies[k].pos.scale(n.bodies[k].mass)
                })
                .scale(1.0 / mass)
            } else {
                Vec3::ZERO
            }
        };
        if com.norm() > 0.0 {
            let kids = self.nodes[i.get()].children.clone();
            {
                let n = &mut self.nodes[i.get()];
                n.motion.offset += com;
                for b in n.bodies.iter_mut() {
                    if b.mass > 0.0 {
                        b.pos -= com;
                    }
                }
            }
            for c in kids {
                if !c.is_none() && self.nodes[c.get()].alive {
                    self.nodes[c.get()].motion.offset -= com;
                }
            }
        }

        self.pin(i);
        self.pin(new);

        // **A split is a partition, so the extensive quantities divide and the
        // intensive ones do not.** `settle` was tried here first and it is the
        // wrong tool: it re-summarises the *whole* matter from the bodies, so a
        // node whose temperature somebody authored comes back at whatever
        // temperature its sampled bodies happen to carry — measured, a 50 K
        // node became 3961 K the frame it first split, because its bodies were
        // drawn from a planet at 2000 K and the summary believed them.
        //
        // The four energies that a body list cannot carry (this node's own
        // potential, and the cohesive, external and chemical terms
        // `Tree::sum_conserved` adds for it) are divided by mass share, which
        // is not their own law but is the only division that conserves exactly.
        // Both ends re-derive them from their own configuration the next time
        // they are drawn.
        let keep = {
            let n = &self.nodes[i.get()];
            if n.matter.mass > 0.0 {
                ((n.matter.mass - matter.mass) / n.matter.mass).clamp(0.0, 1.0)
            } else {
                1.0
            }
        };
        let gone = 1.0 - keep;
        {
            let n = &mut self.nodes[i.get()];
            n.matter.mass = (n.matter.mass - matter.mass).max(0.0);
            n.matter.momentum -= matter.momentum;
            n.matter.internal_energy *= keep;
            n.matter.baryon_number *= keep;
            n.matter.lepton_number *= keep;
            n.matter.charge *= keep;
            n.matter.chemical_energy *= keep;
            n.matter.cohesive_binding *= keep;
            n.matter.external_potential *= keep;
            n.matter.gravitational_binding *= keep;
            n.matter.entropy *= keep;
            n.matter.luminosity *= keep;
            n.potential *= keep;
        }
        {
            let (potential, cohesive, external, chemical) = {
                let n = &self.nodes[i.get()];
                (n.potential, n.matter.cohesive_binding, n.matter.external_potential, n.matter.chemical_energy)
            };
            let share = if keep > 0.0 { gone / keep } else { 0.0 };
            let n = &mut self.nodes[new.get()];
            n.potential = potential * share;
            n.matter.cohesive_binding = cohesive * share;
            n.matter.external_potential = external * share;
            n.matter.chemical_energy = chemical * share;
            // Intensive, so inherited rather than re-measured.
            n.matter.temperature = matter.temperature;
        }
        // What is left of the old node is a different shape and possibly a
        // different size of thing.
        self.remeasure_radius(i);
        self.retier(i);
        self.stats.splits += 1;
        Some(new)
    }

    /// Whether two siblings have become one neighbourhood again.
    ///
    /// The exact inverse of the split's own question, asked of the pair's
    /// contents in the parent's frame — **at the finer of the two resolutions**,
    /// and that asymmetry is the hysteresis. Splitting asks whether a node's
    /// contents are still connected at *its* resolution; merging asks the same
    /// of the pair at the *smaller* one, so two clumps that drift apart and
    /// back again settle rather than flipping between one node and two on
    /// alternate frames.
    ///
    /// Refused for anything with a morphology: two structures coming together
    /// is `docs/PLAY.md` D15's join, which exists, knows about seams, and would
    /// be wrong to bypass by pouring one recipe's parts into another's.
    pub fn would_merge(&self, a: NodeIdx, b: NodeIdx) -> bool {
        if a.is_none() || b.is_none() || a == b {
            return false;
        }
        let (na, nb_node) = (&self.nodes[a.get()], &self.nodes[b.get()]);
        if !na.alive || !nb_node.alive || na.parent != nb_node.parent || na.parent.is_none() {
            return false;
        }
        if na.morphology.is_some() || nb_node.morphology.is_some() {
            return false;
        }
        if !na.is_materialised() || !nb_node.is_materialised() {
            return false;
        }
        // Cheap first: two things that do not even overlap are not one
        // neighbourhood, and this is a subtraction against a grid build.
        let gap = (na.motion.offset - nb_node.motion.offset).norm()
            - na.matter.radius
            - nb_node.matter.radius;
        if gap > 0.0 {
            return false;
        }
        // The pair's contents, in the parent's frame, as one set.
        //
        // Vacated slots are skipped. A node that has been split, or had
        // something re-homed out of it, keeps the slot as a zeroed body —
        // removing it would renumber every sibling after it — and those all sit
        // at the origin with no mass and no radius. Left in, they are a
        // component of their own that nothing else ever reaches.
        let shift = nb_node.motion.offset - na.motion.offset;
        let points: Vec<(Vec3, f64, f64)> = na
            .bodies
            .iter()
            .filter(|x| x.mass > 0.0)
            .map(|x| (x.pos, x.radius, x.mass))
            .chain(
                nb_node
                    .bodies
                    .iter()
                    .filter(|x| x.mass > 0.0)
                    .map(|x| (x.pos + shift, x.radius, x.mass)),
            )
            .collect();
        let n = points.len();
        if n == 0 {
            return false;
        }
        // **At the resolution the merged node would have**, which makes this
        // the exact inverse of the split's own question rather than a second
        // rule. The first attempt used the finer of the two nodes' own
        // resolutions, on the grounds that a stricter test going back than
        // coming apart is hysteresis — and it is so strict that a node's own
        // contents are not connected at it, so nothing could ever merge.
        //
        // It does not oscillate: a pair merges only when the node they would
        // form answers `resolve_extent`'s question with "one region", which is
        // the same answer that node then goes on giving.
        let centre = points.iter().fold(Vec3::ZERO, |acc, p| acc + p.0).scale(1.0 / n as f64);
        let rms = (points.iter().map(|p| (p.0 - centre).norm2()).sum::<f64>() / n as f64)
            .max(0.0)
            .sqrt();
        let within = rms * crate::state::RMS_TO_RADIUS / (n as f64).cbrt();
        let spacing = within.max(2.0 * points.iter().fold(0.0f64, |m, p| m.max(p.1)));
        let grid = crate::neighbourhood::NeighbourGrid::of_points(
            points.iter().map(|p| p.0),
            spacing,
        );
        let mut label: Vec<usize> = (0..n).collect();
        fn find(label: &mut [usize], mut x: usize) -> usize {
            while label[x] != x {
                label[x] = label[label[x]];
                x = label[x];
            }
            x
        }
        let mut candidates = Vec::new();
        for i in 0..n {
            grid.neighbours(points[i].0, &mut candidates);
            for j in candidates.iter().map(|c| *c as usize) {
                if j <= i {
                    continue;
                }
                let d = (points[j].0 - points[i].0).norm() - points[i].1 - points[j].1;
                if d <= within {
                    let (x, y) = (find(&mut label, i), find(&mut label, j));
                    if x != y {
                        label[x.max(y)] = x.min(y);
                    }
                }
            }
        }
        // And the same majority rule, because "is everything connected" is the
        // wrong question here for exactly the reason it is wrong there: one
        // outlier in the tail would keep two clumps sitting on top of each
        // other from ever being one node again. Measured on that case — 512
        // points, 90% of them overlapping their nearest neighbour, and a single
        // straggler 1.65x10^16 m out.
        let mut share: std::collections::BTreeMap<usize, f64> = Default::default();
        let mut total = 0.0;
        for x in 0..n {
            let l = find(&mut label, x);
            *share.entry(l).or_insert(0.0) += points[x].2;
            total += points[x].2;
        }
        total > 0.0 && share.values().fold(0.0f64, |m, v| m.max(*v)) / total > 0.5
    }

    /// Fold one sibling into another.
    ///
    /// The inverse of [`Self::split_off`], and the reason `docs/BACKLOG.md`
    /// says the merge "belongs with it": without one, two clumps that fall back
    /// together stay two nodes for ever, and the tree accumulates a node per
    /// event that ever separated anything.
    ///
    /// `b`'s detail is re-expressed in `a`'s frame and appended; anything
    /// promoted out of `b` goes through [`Self::reparent`] into `a`; `b`'s slot
    /// in the parent is vacated the way a departed body's is, and `b` is freed.
    /// `a` is then re-measured against what it now holds.
    pub fn merge(&mut self, a: NodeIdx, b: NodeIdx) -> bool {
        if !self.would_merge(a, b) {
            return false;
        }
        let parent = self.nodes[a.get()].parent;
        let shift = self.nodes[b.get()].motion.offset - self.nodes[a.get()].motion.offset;
        let relative = crate::coords::velocity_add(
            -self.nodes[a.get()].motion.velocity,
            self.nodes[b.get()].motion.velocity,
        );
        // Promoted children first, while `b` still has the slots they sit in.
        for c in self.nodes[b.get()].children.clone() {
            if !c.is_none() && self.nodes[c.get()].alive {
                self.reparent(c, a);
            }
        }
        let incoming: Vec<Body> = self.nodes[b.get()]
            .bodies
            .iter()
            .filter(|x| x.mass > 0.0)
            .map(|x| {
                let mut x = *x;
                x.pos += shift;
                x.vel = crate::coords::velocity_add(relative, x.vel);
                x
            })
            .collect();
        {
            let n = &mut self.nodes[a.get()];
            n.bodies.extend(incoming);
            n.children.resize(n.bodies.len(), NodeIdx::NONE);
        }
        // Vacate `b`'s slot, then free it. `release_subtree` would take its
        // children with it, which is why they moved first.
        let slot = self.nodes[b.get()].slot as usize;
        {
            let p = &mut self.nodes[parent.get()];
            if slot < p.bodies.len() {
                p.bodies[slot] = Body::default();
            }
            if slot < p.children.len() {
                p.children[slot] = NodeIdx::NONE;
            }
        }
        self.nodes[b.get()].bodies.clear();
        self.release_subtree(b);
        self.pin(a);
        self.pin(parent);
        self.settle(a);
        self.retier(a);
        self.stats.merges += 1;
        true
    }

    /// Move a node under a different parent, re-expressing it in the new frame.
    ///
    /// This is what makes the tree a spatial index rather than a record of what
    /// owns what. Until it existed a node's place was fixed at creation for
    /// life, which is fine for containment that never changes — a star does not
    /// leave its cluster — and wrong for everything at play scale, where a
    /// thing picked up enters your frame and a thing thrown enters the ground's.
    ///
    /// Three things have to happen together, and doing any one without the
    /// others corrupts the tree:
    ///
    /// * **The frame changes.** `motion` is relative to the parent, so it is
    ///   recomputed against the new one — position through the common ancestor,
    ///   velocity by the relativistic composition the rest of the engine uses.
    ///   Computed before anything is mutated, because it reads the old chain.
    /// * **The old slot is vacated.** A promoted body is a *stand-in*:
    ///   `sum_conserved` counts the child and skips the body wherever a slot is
    ///   promoted. Leave the body behind and the old parent silently reclaims
    ///   the mass that just left it.
    /// * **The subtree is rekeyed.** A [`PathKey`] is the path, so the node and
    ///   everything under it get new ones, and anything addressed by key has to
    ///   follow. The returned [`Rehomed`] is that list; `Tree` migrates the
    ///   pinned detail it owns, and the caller migrates the rest.
    ///
    /// Refused, returning `None`, when the move is not a move: the root (the
    /// universe has no outside), a node into itself, a node into its own
    /// descendant (which would make a cycle, and every parent walk in the
    /// engine is a `while` loop that would never end), or a node into the
    /// parent it already has.
    pub fn reparent(&mut self, node: NodeIdx, new_parent: NodeIdx) -> Option<Rehomed> {
        if node.is_none() || new_parent.is_none() || node == new_parent || node == self.root {
            return None;
        }
        if !self.nodes[node.get()].alive || !self.nodes[new_parent.get()].alive {
            return None;
        }
        let old_parent = self.nodes[node.get()].parent;
        if old_parent == new_parent || old_parent.is_none() {
            return None;
        }
        // A cycle would not merely be wrong, it would hang: `lca`, `offset_from`
        // and `disturb` all walk parents with `while !cur.is_none()`.
        if self.lca(node, new_parent) == node {
            return None;
        }

        // Read the old chain before touching anything.
        let offset = self.separation(new_parent, Vec3::ZERO, node, Vec3::ZERO).value;
        let anc = self.lca(node, new_parent);
        let v_node = self.velocity_from(anc, node);
        let v_parent = self.velocity_from(anc, new_parent);
        // The node's velocity as the new parent sees it: boost by minus the
        // parent's own. `velocity_add` is what composes velocities everywhere
        // else in the engine, so the inverse uses it too rather than inventing
        // a second convention.
        let velocity = crate::coords::velocity_add(-v_parent, v_node);

        // The new parent needs a body list to hold a slot in.
        self.refine(new_parent);

        // Vacate. Zeroing rather than removing: `children` is parallel to
        // `bodies` and a sibling's `PathKey` is derived from its slot index, so
        // removing an element would renumber every sibling after it and change
        // the identity of each one.
        let old_slot = self.nodes[node.get()].slot as usize;
        // The kind travels with the object: it is the same thing, so whatever
        // its stand-in was in the old parent's list is what it is in the new
        // one. Read before the slot is cleared.
        let kind = self
            .nodes[old_parent.get()]
            .bodies
            .get(old_slot)
            .map(|b| b.kind)
            .unwrap_or(crate::state::BodyKind::Grain);
        {
            let p = &mut self.nodes[old_parent.get()];
            if old_slot < p.bodies.len() {
                p.bodies[old_slot] = Body::default();
            }
            if old_slot < p.children.len() {
                p.children[old_slot] = NodeIdx::NONE;
            }
        }

        // Take the new slot, with a stand-in body describing what arrived.
        let stand_in = {
            let n = &self.nodes[node.get()];
            Body {
                pos: offset,
                vel: velocity,
                mass: n.matter.mass,
                radius: n.matter.radius,
                charge: n.matter.charge,
                internal_energy: n.matter.internal_energy,
                spin: n.matter.spin,
                temperature: n.matter.temperature,
                composition: n.matter.composition,
                kind,
                ..Default::default()
            }
        };
        let new_slot = {
            let p = &mut self.nodes[new_parent.get()];
            p.bodies.push(stand_in);
            while p.children.len() < p.bodies.len() {
                p.children.push(NodeIdx::NONE);
            }
            let slot = p.bodies.len() - 1;
            p.children[slot] = node;
            slot
        };

        let new_depth = self.nodes[new_parent.get()].depth + 1;
        let new_key = self.nodes[new_parent.get()].key.child(new_slot as u64);
        {
            let n = &mut self.nodes[node.get()];
            n.parent = new_parent;
            n.slot = new_slot as u32;
            n.motion.offset = offset;
            n.motion.velocity = velocity;
        }

        let mut keys = Vec::new();
        self.rekey_subtree(node, new_key, new_depth, &mut keys);

        // Pinned detail is `Tree`'s own path-keyed table, so it moves here.
        // Collected first and reinserted after, because an old key and a new
        // key can belong to different nodes in the same batch.
        let mut moved_detail: Vec<(PathKey, Vec<Body>)> = Vec::new();
        for (old, _) in &keys {
            if let Some(bodies) = self.persisted.remove(old) {
                moved_detail.push((*old, bodies));
            }
        }
        for ((_, new), (_, bodies)) in keys.iter().zip(moved_detail.into_iter()) {
            self.persisted.insert(*new, bodies);
        }

        // Both ends were changed by hand, so neither is what `sample` would
        // produce any more.
        self.pin(old_parent);
        self.pin(node);
        self.stats.reparents += 1;
        Some(Rehomed { moved: node, from: old_parent, to: new_parent, keys })
    }

    /// Give a node and everything under it keys and depths for their new place.
    fn rekey_subtree(
        &mut self,
        node: NodeIdx,
        key: PathKey,
        depth: u32,
        out: &mut Vec<(PathKey, PathKey)>,
    ) {
        let old = self.nodes[node.get()].key;
        {
            let n = &mut self.nodes[node.get()];
            n.key = key;
            n.depth = depth;
        }
        out.push((old, key));
        let children = self.nodes[node.get()].children.clone();
        for (slot, c) in children.iter().enumerate() {
            if !c.is_none() && self.nodes[c.get()].alive {
                self.rekey_subtree(*c, key.child(slot as u64), depth + 1, out);
            }
        }
    }

    /// Advance the epoch of a node, invalidating its procedural detail. Called
    /// when an interaction changes the node's matter enough that the old
    /// sample is no longer a valid representative of it.
    pub fn bump_epoch(&mut self, i: NodeIdx) {
        // The detail is about to be redrawn from a different epoch, so every
        // slot in it is about to mean something else. Anything promoted out of
        // one is folded back rather than left pointing at a list that no
        // longer exists.
        self.shed_children(i);
        let n = &mut self.nodes[i.get()];
        n.epoch = n.epoch.wrapping_add(1);
        n.bodies.clear();
    }

    // -- geometry ---------------------------------------------------------

    pub fn path_to_root(&self, mut i: NodeIdx) -> Vec<NodeIdx> {
        let mut v = Vec::new();
        while !i.is_none() {
            v.push(i);
            i = self.nodes[i.get()].parent;
        }
        v
    }

    /// Is `a` an ancestor of `b` (or the same node)?
    pub fn is_ancestor(&self, a: NodeIdx, b: NodeIdx) -> bool {
        let mut cur = b;
        while !cur.is_none() {
            if cur == a {
                return true;
            }
            cur = self.nodes[cur.get()].parent;
        }
        false
    }

    /// Distance from each child node to its nearest *sibling*, which is the
    /// lookahead the scheduler is entitled to use.
    ///
    /// Ancestors are excluded on purpose. A child node is not a separate system
    /// sitting zero metres from its parent — it *is* part of its parent, and
    /// the parent's matter already accounts for it. Applying the light-speed
    /// constraint between a node and its own ancestor would force a nucleus and
    /// the galaxy containing it into lockstep, at the nucleus's zeptosecond
    /// timestep, which is precisely the catastrophe the multi-rate scheme
    /// exists to avoid. The constraint belongs between *disjoint* regions.
    pub fn sibling_separations(&self, parent: NodeIdx) -> Vec<(NodeIdx, f64)> {
        let p = &self.nodes[parent.get()];
        let kids: Vec<NodeIdx> = p.children.iter().copied().filter(|c| !c.is_none()).collect();
        let mut out = Vec::with_capacity(kids.len());
        for &a in &kids {
            let mut best = f64::INFINITY;
            for &b in &kids {
                if a == b {
                    continue;
                }
                let na = &self.nodes[a.get()];
                let nb = &self.nodes[b.get()];
                // Surface-to-surface: influence has to cross the gap, not the
                // distance between centres.
                let gap = (na.motion.offset - nb.motion.offset).norm()
                    - na.matter.radius
                    - nb.matter.radius;
                best = best.min(gap.max(0.0));
            }
            out.push((a, best));
        }
        out
    }

    /// Lowest common ancestor of two nodes.
    pub fn lca(&self, a: NodeIdx, b: NodeIdx) -> NodeIdx {
        let (mut x, mut y) = (a, b);
        let (mut dx, mut dy) = (self.nodes[a.get()].depth, self.nodes[b.get()].depth);
        while dx > dy {
            x = self.nodes[x.get()].parent;
            dx -= 1;
        }
        while dy > dx {
            y = self.nodes[y.get()].parent;
            dy -= 1;
        }
        while x != y && !x.is_none() && !y.is_none() {
            x = self.nodes[x.get()].parent;
            y = self.nodes[y.get()].parent;
        }
        x
    }

    /// Offset of `(node, local)` from `ancestor`'s origin, accumulating the
    /// round-off honestly. See `coords::Bounded` for why the error bound is
    /// carried rather than assumed negligible.
    pub fn offset_from(&self, ancestor: NodeIdx, mut node: NodeIdx, local: Vec3) -> Bounded {
        let mut acc = Bounded::exact(local);
        while node != ancestor && !node.is_none() {
            let n = &self.nodes[node.get()];
            acc = acc.add(Bounded::exact(n.motion.offset));
            node = n.parent;
        }
        acc
    }

    /// Separation between two points anywhere in the tree.
    ///
    /// The precision of the answer degrades with tree distance, which is the
    /// physically correct behaviour: two nucleons in one nucleus are located
    /// relative to each other to ~10^-31 m, while a nucleon and a star on the
    /// far side of the galaxy are located to ~10^5 m — and nothing couples them
    /// more tightly than that.
    pub fn separation(&self, a: NodeIdx, a_local: Vec3, b: NodeIdx, b_local: Vec3) -> Bounded {
        let anc = self.lca(a, b);
        let pa = self.offset_from(anc, a, a_local);
        let pb = self.offset_from(anc, b, b_local);
        pb.sub(pa)
    }

    /// Velocity of `node` relative to `ancestor`, composed relativistically.
    pub fn velocity_from(&self, ancestor: NodeIdx, mut node: NodeIdx) -> Vec3 {
        let mut chain = Vec::new();
        while node != ancestor && !node.is_none() {
            chain.push(self.nodes[node.get()].motion.velocity);
            node = self.nodes[node.get()].parent;
        }
        let mut v = Vec3::ZERO;
        for u in chain.iter().rev() {
            v = crate::coords::velocity_add(v, *u);
        }
        v
    }

    /// Depth-first walk over live nodes.
    pub fn walk<F: FnMut(NodeIdx, &Node)>(&self, start: NodeIdx, f: &mut F) {
        if start.is_none() || !self.nodes[start.get()].alive {
            return;
        }
        f(start, &self.nodes[start.get()]);
        let kids = self.nodes[start.get()].children.clone();
        for c in kids {
            if !c.is_none() {
                self.walk(c, f);
            }
        }
    }

    /// Total conserved quantities over the whole live tree, counting each
    /// region exactly once: a materialised node is represented by its bodies,
    /// except where a body has been promoted, in which case the child node
    /// speaks for it.
    pub fn total_conserved(&self) -> crate::state::Conserved {
        self.sum_conserved(self.root)
    }

    /// A child's totals, carried from its own frame into its parent's: what
    /// it holds moves with it, so the parent sees it with the child's bulk
    /// momentum added, the bulk's kinetic energy and its cross term with what
    /// the child holds, and the angular momentum all of that has about the
    /// parent's centre. Read where the child is at the parent's instant, which
    /// is the instant the parent's own bodies are at.
    ///
    /// Without it the world's books left every promoted child's bulk motion
    /// out — the stand-in body that carries it is not summed, because the child
    /// speaks for it, and the child spoke only for its contents. Measured on a
    /// turning Earth with one face promoted: the world's momentum read 1.9e26
    /// kg m/s, exactly one face's going round, and moved by 1.2e26 each time
    /// the Earth was solved.
    fn in_parents_frame(&self, c: NodeIdx, t: crate::state::Conserved) -> crate::state::Conserved {
        let n = &self.nodes[c.get()];
        let instant = self.nodes[n.parent.get()].time;
        let at = self.position_at(c, instant);
        let v = self.velocity_at(c, instant);
        let m = n.matter.mass.max(0.0);
        let gamma = crate::coords::gamma(v);
        let bulk = v.scale(m * gamma);
        let momentum = t.momentum + bulk;
        crate::state::Conserved {
            energy: t.energy + (gamma - 1.0) * m * crate::units::C2 + v.dot(t.momentum),
            momentum,
            angular_momentum: t.angular_momentum + n.matter.com.scale(m).cross(v) + at.cross(momentum),
            ..t
        }
    }

    fn sum_conserved(&self, i: NodeIdx) -> crate::state::Conserved {
        let n = &self.nodes[i.get()];
        if !n.is_materialised() {
            // An ocean's matter cannot hold what it has given and taken, and
            // its own books do: `ocean::Account`.
            return match n.ocean.as_ref() {
                Some(o) => n.matter.conserved().add(o.account.conserved()),
                None => n.matter.conserved(),
            };
        }
        let count = n.bodies.len();
        // A slot whose child speaks for it contributes nothing here; the child's
        // own total is added below.
        let stood_for = |slot: usize| {
            let c = n.child_of(slot);
            !c.is_none() && self.nodes[c.get()].alive
        };
        // **Rest mass is summed apart from everything else**, and pairwise, for
        // the reason `Matter::non_rest_energy` documents at length: rest energy
        // exceeds every other term by around 10^16 for ordinary matter, so a
        // running total that mixes the two spends its significant digits on
        // `mc^2` and reports the interesting part as noise. Sequentially
        // summing 20 000 bodies of 10^52 J apiece cost 1.6x10^45 J on the
        // reference galaxy — 10^-11 of the total, and a hundred times the
        // staleness it was being compared against. `summarise` has always
        // grouped its terms this way; this is the other account learning the
        // same lesson.
        let rest = crate::math::det_sum_by(count, &|s| {
            if stood_for(s) { 0.0 } else { n.bodies[s].mass }
        });
        let non_rest = crate::math::det_sum_by(count, &|s| {
            if stood_for(s) {
                0.0
            } else {
                let b = &n.bodies[s];
                (crate::coords::gamma(b.vel) - 1.0) * b.mass * crate::units::C2
                    + b.internal_energy
            }
        });
        let momentum = crate::math::det_sum_v3_by(count, &|s| {
            if stood_for(s) { Vec3::ZERO } else { n.bodies[s].momentum() }
        });
        let angular = crate::math::det_sum_v3_by(count, &|s| {
            if stood_for(s) {
                Vec3::ZERO
            } else {
                let b = &n.bodies[s];
                b.pos.cross(b.momentum()) + b.spin
            }
        });
        let charge = crate::math::det_sum_by(count, &|s| {
            if stood_for(s) { 0.0 } else { n.bodies[s].charge }
        });
        let baryon = crate::math::det_sum_by(count, &|s| {
            if stood_for(s) {
                0.0
            } else {
                let b = &n.bodies[s];
                b.mass * b.composition.nucleons_per_kg()
            }
        });
        let lepton = crate::math::det_sum_by(count, &|s| {
            if stood_for(s) {
                0.0
            } else {
                let b = &n.bodies[s];
                b.mass * b.composition.nucleons_per_kg() * b.composition.electrons_per_nucleon()
                    - b.charge / crate::units::E_CHARGE
            }
        });
        let mut total = crate::state::Conserved {
            energy: rest * crate::units::C2 + non_rest,
            momentum,
            angular_momentum: angular,
            charge,
            baryon,
            lepton,
        };
        for slot in 0..count {
            let c = n.child_of(slot);
            if !c.is_none() && self.nodes[c.get()].alive {
                total = total.add(self.in_parents_frame(c, self.sum_conserved(c)));
            }
        }
        total.energy += n.potential;
        // **Three terms a body list cannot carry.** `summarise` says so in as
        // many words — "not knowable from the children alone; the caller
        // reinstates these" — and this is the other side of that sentence. A
        // `Body` has a mass, a velocity and an internal energy; it has nowhere
        // to put the grip of a dark halo that is not being refined, the bonds
        // holding a solid together below the scale of any body, or the free
        // energy a structure is storing. Leave them out here and a node's
        // energy *changes when it materialises*, which is the one thing a scale
        // transform may never do.
        //
        // Measured, and it is how this was found: the reference world's root is
        // a galaxy whose halo contributes -3.861x10^48 J of `external_potential`
        // against a total of 1.644x10^56, and the tree's total jumped by exactly
        // that — 2.35x10^-8 relative — the moment the root was materialised.
        // `docs/BACKLOG.md` had recorded that number as the staleness of a
        // solved node's matter, because in that world the root is both the only
        // node carrying an external potential and the only node solved every
        // frame. The staleness is real and is what [`Tree::settle`] closes, but
        // it is 10^-11, not 2.35x10^-8.
        total.energy +=
            n.matter.cohesive_binding + n.matter.external_potential + n.matter.chemical_energy;
        total
    }
}

/// Carry a point held in a frame turning at `w` for `dt`, in either direction:
/// where it is and how fast it is going afterwards, in the same axes.
///
/// Held means supported at its height: what it does *relative to the turning
/// frame* — `u = v - w x r` — is to go round the centre along the surface at
/// its tangential part and to change its height at its radial part, keeping
/// both. So the relative motion is a great circle, `|u_t| dt / |r|` of arc
/// about `r x u_t`, and the whole is then turned by `w dt`. A straight line in
/// the turning frame would leave the curve by `(u t)^2 / 2R`: ten metres a
/// second along the surface of an Earth gains a metre of height in six minutes
/// and a kilometre in three hours.
pub fn turning_carry(w: Vec3, r: Vec3, v: Vec3, dt: f64) -> (Vec3, Vec3) {
    let u = v - w.cross(r);
    let d = r.norm();
    let (at, rel) = if d > 0.0 {
        let up = r.scale(1.0 / d);
        let radial = u.dot(up);
        let along = u - up.scale(radial);
        let height = d + radial * dt;
        let s = along.norm();
        let (dir, tangent) = if s > 0.0 {
            let axis = up.cross(along).unit();
            let arc = crate::math::Quat::from_rate(axis.scale(s / d), dt);
            (arc.rotate(up), arc.rotate(along))
        } else {
            (up, along)
        };
        (dir.scale(height), tangent + dir.scale(radial))
    } else {
        (r + u.scale(dt), u)
    };
    let turn = crate::math::Quat::from_rate(w, dt);
    let at = turn.rotate(at);
    (at, w.cross(at) + turn.rotate(rel))
}

/// [`turning_carry`] for a thing whose parent's centre is at `about` from the
/// centre it goes round: carried about that centre, and handed back relative to
/// where its parent's centre has gone round to. The parent is taken to be going
/// round rigidly, which is what a held parent is.
pub fn turning_carry_about(w: Vec3, about: Vec3, r: Vec3, v: Vec3, dt: f64) -> (Vec3, Vec3) {
    if about == Vec3::ZERO {
        return turning_carry(w, r, v, dt);
    }
    let (at, vel) = turning_carry(w, about + r, w.cross(about) + v, dt);
    let about = crate::math::Quat::from_rate(w, dt).rotate(about);
    (at - about, vel - w.cross(about))
}

// ---------------------------------------------------------------------------
// Tier, and keeping it true
// ---------------------------------------------------------------------------

/// The tier something of this size belongs to, inside a parent at `parent_tier`.
///
/// One rule with two callers, and having had only one of them is the whole of
/// `docs/BACKLOG.md`'s "a node's tier is decided once and never revisited":
/// [`Tree::promote`] asks it of a body it is about to turn into a node, and
/// [`Tree::retier`] asks it of a node whose size has changed since. Before
/// `retier` existed a terrain patch promoted out of a moon inherited a body
/// radius in the planetary band, was emplaced 1.4 km across, and stayed
/// `Planetary` — two tiers from what its own radius said.
///
/// The clamp against the parent is structural rather than physical: a node is
/// inside its parent, so it cannot be larger than one, and a rounding that said
/// otherwise would put a child on a coarser solver than the thing containing
/// it.
pub fn tier_for(radius: f64, parent_tier: Tier) -> Tier {
    Tier::containing(radius).max(parent_tier)
}

/// Reconcile a refinement policy with the tier it is about to be used at.
///
/// A spec meant for a *coarser* scale than the node turned out to be is
/// replaced by the tier's own policy, with the caller's count read as a budget.
/// A spec meant for a finer one is kept: asking to split an atom into nucleons
/// is a deliberate step down and not a mistake, and overriding it would leave
/// the ladder unable to reach its own bottom.
///
/// This travels with [`tier_for`] and is not optional alongside it. A tier that
/// moves without its spec is the failure `promote` documents at length —
/// materialising under a policy meant for a different scale, which is how eight
/// thousand molecules ended up inside a node the size of an atom — and moving
/// the tier at a *later* point than promotion reopens exactly that hole.
pub fn spec_for(tier: Tier, spec: SampleSpec) -> SampleSpec {
    if crate::sampler::tier_of(spec.kind) >= tier {
        spec
    } else {
        crate::sampler::budgeted_spec(tier, spec.count)
    }
}
