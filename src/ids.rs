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
