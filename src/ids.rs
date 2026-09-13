//! Stable identity for nodes across materialisation cycles.
//!
//! A node's *arena index* is an implementation detail that changes every time
//! detail is discarded and rebuilt. Its *path key* is not: it is a hash of the
//! child-index path from the root, so the third child of the seventh child of
//! the root has the same key today, after a coarsen, and after a reload. All
//! randomness and all ledger entries are addressed by path key.

use crate::rng::mix2;

/// Arena handle. Cheap, dense, invalidated by coarsening.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeIdx(pub u32);

impl NodeIdx {
    pub const NONE: NodeIdx = NodeIdx(u32::MAX);
    #[inline]
    pub fn is_none(self) -> bool {
        self.0 == u32::MAX
    }
    #[inline]
    pub fn get(self) -> usize {
        self.0 as usize
    }
}

/// An issued name for a thing, as distinct from where that thing currently is.
///
/// [`PathKey`] does two jobs and is good at one of them. As an *address* it is
/// perfect: derived from the child-index path, so it survives a node being
/// discarded and rebuilt, and it seeds the sampler. As an *identity* it fails
/// the moment anything moves — `reparent` changes a node's address, and every
/// side table keyed by the old one is then pointing at a stranger.
///
/// So identity is issued rather than derived. An `EntityId` is handed out once,
/// never reused, and does not change for any reason: not a move, not a
/// coarsen, not a reload. What a thing *is* stops being a function of where it
/// happens to be.
///
/// # Why this is not just another key
///
/// The trap it removes is that a node's identity was spread across side tables
/// on `World` — chemistry, environment, clocks, history — each keyed by address,
/// and `World::reparent` held the only enumeration of them. Adding a table
/// meant remembering to add a line there, and forgetting meant a moved object
/// arrived without its chemistry, silently. Keyed by `EntityId` the tables do
/// not move when the thing does, so there is nothing to remember.
///
/// One index does still have to be migrated, and it is worth being honest that
/// this is not quite "by construction": a node discarded and rebuilt has to
/// recover its identity from *somewhere*, and the only somewhere is a map from
/// address to identity. `reparent` migrates that one map. The difference is
/// that it is one line in one place for ever, rather than one line per table
/// for ever.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EntityId(pub u64);

impl EntityId {
    /// Not a thing. Distinct from any issued id, which start at one.
    pub const NONE: EntityId = EntityId(0);

    #[inline]
    pub fn is_none(self) -> bool {
        self.0 == 0
    }
}

impl std::fmt::Display for EntityId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_none() {
            f.write_str("@none")
        } else {
            write!(f, "@{}", self.0)
        }
    }
}

/// Persistent identity. 128 bits: with ~10^12 live nodes the collision
/// probability over the life of a simulation is ~10^-15.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PathKey(pub u128);

impl PathKey {
    pub const ROOT: PathKey = PathKey(0x5EED_0000_0000_0001_0000_0000_0000_0001);

    /// Derive a child's key. Deliberately *not* a plain concatenation: paths
    /// can be thousands deep (galaxy → … → nucleus is 7 tiers but each tier
    /// may nest many levels), and a rolling hash keeps the key fixed-width.
    pub fn child(self, index: u64) -> PathKey {
        let lo = self.0 as u64;
        let hi = (self.0 >> 64) as u64;
        let nlo = mix2(lo, index);
        let nhi = mix2(hi, nlo ^ index.rotate_left(17));
        PathKey(((nhi as u128) << 64) | nlo as u128)
    }

    /// Short human-readable form for logs and the ledger UI.
    pub fn short(self) -> String {
        format!("{:016x}", (self.0 >> 64) as u64)
    }

    #[inline]
    pub fn lo(self) -> u64 {
        self.0 as u64
    }
}

impl std::fmt::Display for PathKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "#{}", self.short())
    }
}
