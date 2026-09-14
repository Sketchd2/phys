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
use crate::sampler::{sample, SampleReport, SampleSpec};
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
    pub pinned: bool,
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
    /// Worst conservation error seen across every scale transition so far.
    pub worst_conservation_error: f64,
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
            potential: root_agg.binding_energy,
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
            bubble: 1.0,
            alive: true,
            morphology: None,
            topology: None,
            steps_taken: 0,
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
        if self.nodes[i.get()].is_materialised() {
            return &self.nodes[i.get()].bodies;
        }
        let key = self.nodes[i.get()].key;

        // Pinned detail was altered by an interaction, so it cannot be
        // regenerated — it comes back from the persistent store instead.
        if let Some(saved) = self.persisted.get(&key) {
            let bodies = saved.clone();
            let n = &mut self.nodes[i.get()];
            n.children = vec![NodeIdx::NONE; bodies.len()];
            n.bodies = bodies;
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
        let (bodies, topo, report) = match &morph {
            Some(m) => {
                let (b, t, r) = crate::sampler::sample_structured(
                    &matter,
                    m,
                    spec.count,
                    self.world_seed,
                    key.0,
                    epoch,
                    gravity,
                );
                (b, Some(t), r)
            }
            None => {
                let (b, r) = sample(&matter, spec, self.world_seed, key.0, epoch);
                (b, None, r)
            }
        };
        self.stats.materialisations += 1;
        self.stats.bodies_created += bodies.len() as u64;
        self.stats.worst_conservation_error = self
            .stats
            .worst_conservation_error
            .max(report.conservation_error);
        let n = &mut self.nodes[i.get()];
        n.children = vec![NodeIdx::NONE; bodies.len()];
        n.bodies = bodies;
        n.topology = topo;
        n.potential = report.potential;
        n.last_report = report;
        &self.nodes[i.get()].bodies
    }

    /// Turn one materialised body into a node of its own, one tier finer.
    ///
    /// The child's matter is *the body itself*, reinterpreted: same mass,
    /// same composition, same momentum in the parent's frame. Nothing is
    /// invented at this step — invention happens when the child is refined.
    pub fn promote(&mut self, i: NodeIdx, slot: usize, spec: SampleSpec) -> NodeIdx {
        self.refine(i);
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
        // The child's own frame carries the bulk motion, so inside its frame the
        // net momentum is zero — that is what "rest frame" means. Bulk motion is
        // never double-counted.
        matter.momentum = Vec3::ZERO;
        matter.internal_energy = body.internal_energy.max(matter.thermal_energy());
        matter.luminosity = crate::state::stefan_boltzmann(matter.radius, matter.temperature);

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
                // A body's spin becomes the node's rotation. Orientation starts
                // at identity: the child's body frame *is* how it was sampled,
                // and everything after is what the rotation did to it.
                orientation: crate::math::Quat::IDENTITY,
                spin_rate: matter.angular_velocity(),
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
            bubble: 1.0,
            alive: true,
            morphology: None,
            topology: None,
            steps_taken: 0,
            last_report: SampleReport::default(),
        };
        let idx = self.alloc(child);
        self.nodes[i.get()].children[slot] = idx;
        self.stats.promotions += 1;
        idx
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
        let children = self.nodes[i.get()].children.clone();
        for (slot, c) in children.iter().enumerate() {
            if !c.is_none() {
                self.sync_from_child(i, slot, *c);
                self.release_subtree(*c);
            }
        }

        let (before, potential, pinned, key) = {
            let n = &self.nodes[i.get()];
            (n.matter.conserved(), n.potential, n.pinned, n.key)
        };
        let bodies = std::mem::take(&mut self.nodes[i.get()].bodies);
        let mut matter = summarise(&bodies, potential);
        matter.external_potential = self.nodes[i.get()].matter.external_potential;
        matter.chemical_energy = self.nodes[i.get()].matter.chemical_energy;
        matter.entropy_exported = self.nodes[i.get()].matter.entropy_exported;
        let scales = crate::state::Scales::of(&bodies);
        let err = matter.conserved().error_against(&before, &scales);

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
        if err < IDEMPOTENT_TOLERANCE && !pinned {
            self.stats.coarsenings += 1;
            self.stats.idempotent_coarsenings += 1;
            let n = &mut self.nodes[i.get()];
            n.children.clear();
            return err;
        }

        let n = &mut self.nodes[i.get()];
        // Preserve the node's own frame-level bookkeeping: `summarise` measures
        // the children in the node's frame, so the node's momentum and com are
        // updated, but its tier, spec and identity are untouched.
        n.matter.mass = matter.mass;
        n.matter.com = matter.com;
        n.matter.momentum = matter.momentum;
        n.matter.spin = matter.spin;
        n.matter.internal_energy = matter.internal_energy;
        n.matter.binding_energy = matter.binding_energy;
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
        n.children.clear();
        self.stats.coarsenings += 1;
        self.stats.worst_conservation_error = self.stats.worst_conservation_error.max(err);
        err
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
        let p = &mut self.nodes[parent.get()];
        if let Some(b) = p.bodies.get_mut(slot) {
            b.mass = mass;
            b.composition = comp;
            b.temperature = temp;
            b.charge = charge;
            b.spin = spin;
            b.internal_energy = internal;
            b.radius = radius;
            b.pos = frame.offset;
            b.vel = frame.velocity;
        }
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
            self.nodes[child.get()].motion.velocity = self.nodes[child.get()].motion.velocity + dv;
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
    /// # What it does not do yet: orientation
    ///
    /// The answer is in the node's *frame axes*, which today are its parent's,
    /// because [`Tree::offset_from`] composes offsets and not rotations. A
    /// `Motion` carries an `orientation` and this does not consult it. So a node
    /// on the `+x` side of a planet is told gravity points along `-x`, which is
    /// true in the parent's axes and is not what a structure generated with
    /// `+z` up expects to carry its weight along.
    ///
    /// It does not bite yet because nothing orients a node against the body it
    /// sits on — there is no terrain, which is the same `PLAY.md` D6 work this
    /// field was built for. It will the moment there is: a patch on a sphere is
    /// oriented by definition. Recorded in `BACKLOG.md`.
    pub fn gravity_at(&self, idx: NodeIdx) -> Vec3 {
        if idx.is_none() || !self.nodes[idx.get()].alive {
            return Vec3::ZERO;
        }
        let own = self.nodes[idx.get()].matter.mass.max(0.0);
        let mut g = Vec3::ZERO;
        let mut anc = self.nodes[idx.get()].parent;
        while !anc.is_none() {
            let a = &self.nodes[anc.get()];
            let r = self.offset_from(anc, idx, Vec3::ZERO).value;
            let d = r.norm();
            let (m, radius) = (a.matter.mass.max(0.0), a.matter.radius);
            if d > 0.0 && radius > 0.0 && m > own {
                let source = m - own;
                let enclosed = if d >= radius {
                    source
                } else {
                    source * (d / radius).powi(3)
                };
                g += r.scale(-crate::units::G * enclosed / (d * d * d));
            }
            anc = a.parent;
        }
        g
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
                child.alive.then(|| (child.motion.offset, child.matter.radius))
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
        let n = &mut self.nodes[i.get()];
        n.matter.radius = m.extent().max(n.matter.radius.min(1e-3)).max(1e-30);
        n.matter.chemical_energy = m.stored_energy();
        n.bodies.clear();
        n.children.clear();
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
        let n = &mut self.nodes[i.get()];
        n.matter.radius = m.extent().max(1e-30);
        n.matter.chemical_energy = m.stored_energy();
        n.bodies.clear();
        n.children.clear();
        n.morphology = Some(m);
        self.stats.structures += 1;
        // A structure takes its size from its program the instant it has one,
        // and a seed is not the size of the tree it becomes.
        self.retier(i);
        self.nodes[i.get()].morphology.as_mut().unwrap()
    }

    /// Mark a node — and its whole ancestry — as holding non-derivable detail.
    /// Ancestors must be pinned too: a changed child means the parent's
    /// materialisation no longer matches what `sample` would produce.
    pub fn pin(&mut self, i: NodeIdx) {
        let mut cur = i;
        while !cur.is_none() {
            let n = &mut self.nodes[cur.get()];
            n.pinned = true;
            n.residency = Residency::Pinned;
            cur = n.parent;
        }
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
        let n = &mut self.nodes[i.get()];
        n.epoch = n.epoch.wrapping_add(1);
        n.bodies.clear();
        n.children.clear();
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

    fn sum_conserved(&self, i: NodeIdx) -> crate::state::Conserved {
        let n = &self.nodes[i.get()];
        if !n.is_materialised() {
            return n.matter.conserved();
        }
        let mut total = crate::state::Conserved::zero();
        for (slot, b) in n.bodies.iter().enumerate() {
            let c = n.child_of(slot);
            if !c.is_none() && self.nodes[c.get()].alive {
                total = total.add(self.sum_conserved(c));
            } else {
                total = total.add(crate::state::Conserved {
                    energy: (crate::coords::gamma(b.vel)) * b.mass * crate::units::C2
                        + b.internal_energy,
                    momentum: b.momentum(),
                    angular_momentum: b.pos.cross(b.momentum()) + b.spin,
                    charge: b.charge,
                    baryon: b.mass * b.composition.nucleons_per_kg(),
                    lepton: b.mass * b.composition.nucleons_per_kg()
                        * b.composition.electrons_per_nucleon()
                        - b.charge / crate::units::E_CHARGE,
                });
            }
        }
        total.energy += n.potential;
        total
    }
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
