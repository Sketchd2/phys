//! The scale tree: the object that is 10^68 particles without storing them.
//!
//! A `Node` is a region of space at a tier, holding a bulk `Aggregate`. It has
//! two optional finer representations:
//!
//! * **materialised bodies** — a `Vec<Body>` produced by `prolong`. Cheap to
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
//! it safe is stated in `prolong.rs` and enforced in `tests/consistency.rs`:
//! restriction after prolongation returns the same conserved tuple. Detail that
//! has been *touched* — measured, or hit by something — is different, and is
//! pinned (see `Node::pinned` and `observe::Ledger`).

use crate::coords::{Frame, Located};
use crate::ids::{NodeIdx, PathKey};
use crate::math::Vec3;
use crate::prolong::{prolong, ProlongReport, ProlongSpec};
use crate::state::{restrict, Aggregate, Body};
use crate::units::Tier;
use std::collections::HashMap;

/// How closely a restriction must agree with the stored aggregate before the
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

    /// Bulk state. Always present — this is what a node *is*.
    pub agg: Aggregate,
    /// Position and velocity relative to the parent node's frame.
    pub frame: Frame,

    /// Fine detail, if currently materialised.
    pub bodies: Vec<Body>,
    /// Self-potential the materialisation was built against. Restriction must
    /// use this same number (see `ProlongReport::potential`).
    pub potential: f64,
    /// Children promoted from `bodies`; `NodeIdx::NONE` where not promoted.
    /// Parallel to `bodies`, and empty when nothing is promoted.
    pub children: Vec<NodeIdx>,

    /// How this node is to be split when refined.
    pub spec: ProlongSpec,
    /// Bumped whenever a recorded interaction changes the node's contents.
    /// Detail regenerated at the same epoch is identical; a new epoch means the
    /// old detail is gone for good.
    pub epoch: u32,
    /// Local coordinate time reached by this node's integrator.
    pub time: f64,
    pub residency: Residency,
    /// Set when the node's detail has been altered away from what `prolong`
    /// would produce, so it must be stored rather than regenerated.
    pub pinned: bool,
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
    pub last_report: ProlongReport,
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
    /// Coarsenings where the fine detail turned out to say nothing new, so the
    /// coarse state was left exactly as it was.
    pub idempotent_coarsenings: u64,
    /// Nodes carrying a developmental state.
    pub structures: u64,
    /// Growth and construction steps advanced on aggregates, without ever
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
    pub persisted_bodies: u64,
    /// Worst conservation error seen across every scale transition so far.
    pub worst_conservation_error: f64,
}

impl Tree {
    pub fn new(world_seed: u64, root_agg: Aggregate, tier: Tier, spec: ProlongSpec) -> Tree {
        let root = Node {
            key: PathKey::ROOT,
            parent: NodeIdx::NONE,
            slot: 0,
            depth: 0,
            tier,
            agg: root_agg,
            frame: Frame::default(),
            bodies: Vec::new(),
            potential: root_agg.binding_energy,
            children: Vec::new(),
            spec,
            epoch: 0,
            time: 0.0,
            residency: Residency::Speculative,
            pinned: false,
            alive: true,
            morphology: None,
            topology: None,
            steps_taken: 0,
            last_report: ProlongReport::default(),
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

        let (agg, spec, epoch, morph) = {
            let n = &self.nodes[i.get()];
            (n.agg, n.spec, n.epoch, n.morphology.clone())
        };
        let (bodies, topo, report) = match &morph {
            Some(m) => {
                let (b, t, r) = crate::prolong::prolong_structured(
                    &agg,
                    m,
                    spec.count,
                    self.world_seed,
                    key.0,
                    epoch,
                );
                (b, Some(t), r)
            }
            None => {
                let (b, r) = prolong(&agg, spec, self.world_seed, key.0, epoch);
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
    /// The child's aggregate is *the body itself*, reinterpreted: same mass,
    /// same composition, same momentum in the parent's frame. Nothing is
    /// invented at this step — invention happens when the child is refined.
    pub fn promote(&mut self, i: NodeIdx, slot: usize, spec: ProlongSpec) -> NodeIdx {
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
        let tier = Tier::containing(body.radius).max(parent_tier);

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
        let spec = if crate::prolong::tier_of(spec.kind) >= tier {
            spec
        } else {
            crate::prolong::budgeted_spec(tier, spec.count)
        };

        let mut agg = Aggregate::neutral(body.mass, body.radius.max(1e-30), body.temperature, body.composition);
        agg.charge = body.charge;
        agg.spin = body.spin;
        // The child's own frame carries the bulk motion, so inside its frame the
        // net momentum is zero — that is what "rest frame" means. Bulk motion is
        // never double-counted.
        agg.momentum = Vec3::ZERO;
        agg.internal_energy = body.internal_energy.max(agg.thermal_energy());
        agg.luminosity = crate::state::stefan_boltzmann(agg.radius, agg.temperature);

        let child = Node {
            key,
            parent: i,
            slot: slot as u32,
            depth,
            tier,
            agg,
            frame: Frame {
                offset: body.pos,
                velocity: body.vel,
                proper_time: self.nodes[i.get()].frame.proper_time,
            },
            bodies: Vec::new(),
            potential: 0.0,
            children: Vec::new(),
            spec,
            epoch: 0,
            time: self.nodes[i.get()].time,
            residency: Residency::Speculative,
            pinned: false,
            alive: true,
            morphology: None,
            topology: None,
            steps_taken: 0,
            last_report: ProlongReport::default(),
        };
        let idx = self.alloc(child);
        self.nodes[i.get()].children[slot] = idx;
        self.stats.promotions += 1;
        idx
    }

    /// Fold fine detail back into the bulk state and free it.
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
            (n.agg.conserved(), n.potential, n.pinned, n.key)
        };
        let bodies = std::mem::take(&mut self.nodes[i.get()].bodies);
        let mut agg = restrict(&bodies, potential);
        agg.external_potential = self.nodes[i.get()].agg.external_potential;
        agg.chemical_energy = self.nodes[i.get()].agg.chemical_energy;
        agg.entropy_exported = self.nodes[i.get()].agg.entropy_exported;
        let scales = crate::state::Scales::of(&bodies);
        let err = agg.conserved().error_against(&before, &scales);

        if pinned {
            self.stats.persisted_bodies += bodies.len() as u64;
            self.persisted.insert(key, bodies);
        } else {
            self.stats.bodies_discarded += bodies.len() as u64;
        }

        // If the detail did not actually change the bulk state — the usual case
        // when a user simply pans away — keep the coarse state as the
        // authority rather than overwriting it with a restriction that differs
        // only by round-off.
        //
        // This is what makes "leave and come back" *exactly* idempotent rather
        // than merely accurate. Without it, every visit perturbs the aggregate
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
        // Preserve the node's own frame-level bookkeeping: `restrict` measures
        // the children in the node's frame, so the node's momentum and com are
        // updated, but its tier, spec and identity are untouched.
        n.agg.mass = agg.mass;
        n.agg.com = agg.com;
        n.agg.momentum = agg.momentum;
        n.agg.spin = agg.spin;
        n.agg.internal_energy = agg.internal_energy;
        n.agg.binding_energy = agg.binding_energy;
        n.agg.radius = agg.radius;
        n.agg.temperature = agg.temperature;
        n.agg.composition = agg.composition;
        n.agg.charge = agg.charge;
        n.agg.baryon_number = agg.baryon_number;
        n.agg.lepton_number = agg.lepton_number;
        // Entropy needs two corrections that the original one-liner got wrong
        // as soon as anything in the world could become more ordered.
        //
        // First, coarse-graining itself may only *increase* entropy — you know
        // less once the detail is gone — but the quantity that is monotonic is
        // the total, local plus exported. A structure that grew since the last
        // visit has legitimately lowered its local entropy, and clamping it back
        // up would silently destroy the record of that and unbalance the books.
        //
        // Second, for a structured node `restrict` is not entitled to an
        // opinion at all: it sees an unstructured heap of parts and reports the
        // entropy of the same mass as a gas, which erases precisely the order
        // that makes the thing a structure. `Body` carries no topology, so the
        // information is not there to be recovered. The developmental state is
        // the authority.
        if n.morphology.is_none() {
            let total_stored = n.agg.total_entropy();
            let total_restricted = agg.entropy + n.agg.entropy_exported;
            if total_restricted >= total_stored {
                n.agg.entropy = agg.entropy;
            }
        }
        n.agg.luminosity = agg.luminosity;
        // The morphology owns the structure's size, for the same reason it owns
        // its entropy: `restrict` measures the parts, but what the parts add up
        // to is the program's business.
        if let Some(m) = &n.morphology {
            n.agg.radius = m.extent().max(1e-30);
        }
        n.children.clear();
        self.stats.coarsenings += 1;
        self.stats.worst_conservation_error = self.stats.worst_conservation_error.max(err);
        err
    }

    /// Write a promoted child's evolved bulk state back into the parent's body.
    fn sync_from_child(&mut self, parent: NodeIdx, slot: usize, child: NodeIdx) {
        let (mass, comp, temp, charge, spin, internal, radius, frame) = {
            let c = &self.nodes[child.get()];
            (
                c.agg.mass,
                c.agg.composition,
                c.agg.temperature,
                c.agg.charge,
                c.agg.spin,
                c.agg.internal_energy,
                c.agg.radius,
                c.frame,
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

    /// Give a node a developmental state, turning it from a statistical
    /// population into a structure with a history.
    ///
    /// The aggregate's radius, chemical energy and entropy are taken over by
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
        m.built = (self.nodes[i.get()].agg.mass * 1e-3).clamp(1e-6, 1.0);
        let n = &mut self.nodes[i.get()];
        n.agg.radius = m.extent().max(n.agg.radius.min(1e-3)).max(1e-30);
        n.agg.chemical_energy = m.stored_energy();
        n.bodies.clear();
        n.children.clear();
        n.morphology = Some(m);
        self.stats.structures += 1;
        self.nodes[i.get()].morphology.as_mut().unwrap()
    }

    /// Mark a node — and its whole ancestry — as holding non-derivable detail.
    /// Ancestors must be pinned too: a changed child means the parent's
    /// materialisation no longer matches what `prolong` would produce.
    pub fn pin(&mut self, i: NodeIdx) {
        let mut cur = i;
        while !cur.is_none() {
            let n = &mut self.nodes[cur.get()];
            n.pinned = true;
            n.residency = Residency::Pinned;
            cur = n.parent;
        }
    }

    /// Advance the epoch of a node, invalidating its procedural detail. Called
    /// when an interaction changes the node's bulk state enough that the old
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
    /// the parent's aggregate already accounts for it. Applying the light-speed
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
                let gap = (na.frame.offset - nb.frame.offset).norm()
                    - na.agg.radius
                    - nb.agg.radius;
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
    /// round-off honestly. See `coords::Located` for why the error bound is
    /// carried rather than assumed negligible.
    pub fn offset_from(&self, ancestor: NodeIdx, mut node: NodeIdx, local: Vec3) -> Located {
        let mut acc = Located::exact(local);
        while node != ancestor && !node.is_none() {
            let n = &self.nodes[node.get()];
            acc = acc.add(Located::exact(n.frame.offset));
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
    pub fn separation(&self, a: NodeIdx, a_local: Vec3, b: NodeIdx, b_local: Vec3) -> Located {
        let anc = self.lca(a, b);
        let pa = self.offset_from(anc, a, a_local);
        let pb = self.offset_from(anc, b, b_local);
        pb.sub(pa)
    }

    /// Velocity of `node` relative to `ancestor`, composed relativistically.
    pub fn velocity_from(&self, ancestor: NodeIdx, mut node: NodeIdx) -> Vec3 {
        let mut chain = Vec::new();
        while node != ancestor && !node.is_none() {
            chain.push(self.nodes[node.get()].frame.velocity);
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
            return n.agg.conserved();
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
