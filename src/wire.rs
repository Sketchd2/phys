//! A binary wire format for the persistent world.
//!
//! # Why this is hand-written
//!
//! The crate has no dependencies, and this file is the one place that was most
//! tempting to add some. A derive macro would write the encoders for free — and
//! would also decide the format, silently, from the order the fields happen to
//! be declared in. A world file outlives the code that wrote it: it has to be
//! readable by a build six months and four refactors later, which means the
//! layout is a decision to be made deliberately once rather than a consequence
//! of `struct` field order. So the format is explicit, versioned at the header,
//! and every field is named at the call site in `persist.rs`.
//!
//! The cost is that adding a field means editing two functions instead of none,
//! and forgetting the second is a silent data-loss bug. `tests/persistence.rs`
//! answers that with a round-trip over a fully-populated value of every type:
//! a field that is written and not read, or read and not written, fails it.
//!
//! # Reading is hostile-input code
//!
//! A world file arrives from disk, and later from a network. Every read is
//! bounds-checked, no read panics, and a length prefix is validated against the
//! bytes actually remaining *before* anything is allocated — otherwise a
//! four-byte claim of four billion elements is an out-of-memory abort rather
//! than a parse error.
//!
//! # Floats round-trip by bits
//!
//! `f64` is stored as `to_bits`, not as a decimal. The engine's central claim
//! is that leaving a region and coming back returns *bit-for-bit* the same
//! state (`tree::IDEMPOTENT_TOLERANCE` and the idempotent-coarsening path), and
//! a save that quietly renormalised `-0.0` or flattened a NaN payload would
//! break that in a way no conservation check would notice.

use crate::math::{Quat, Vec3};

/// Every world file starts with these four bytes.
pub const MAGIC: [u8; 4] = *b"PHYS";

/// The format version. Bump it whenever the layout changes in a way an older
/// reader would misinterpret; add a migration in `persist::load` when you do.
///
/// * **3** — the world carries its substance catalogue and every node's
///   mixture. A node's mixture names substances by position in that catalogue,
///   so the two are one change: a file with mixtures and no registry would
///   point at nothing.
/// * **2** — a node carries its `bubble`, the administrative multiplier on how
///   fast its interior runs. A version 1 file has no such field and every node
///   in it was implicitly at 1.0, but the field sits in the middle of the node
///   payload rather than at the end, so a v1 reader and a v2 file disagree
///   about every byte after it. Nothing has shipped that writes v1, so there is
///   no migration: the version simply refuses the old layout instead of
///   misreading it.
pub const FORMAT_VERSION: u16 = 3;

/// A hard ceiling on any single length prefix, independent of the bytes
/// available. Nothing legitimate in this engine has a billion of anything in
/// one vector, and the cap keeps a corrupt length from being merely *plausible*
/// enough to survive the remaining-bytes check on a very large file.
pub const MAX_SEQUENCE: usize = 1 << 28;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WireError {
    /// The read ran off the end of the buffer.
    Truncated { what: &'static str, need: usize, have: usize },
    /// The first four bytes were not `PHYS`.
    BadMagic,
    /// A version this build does not know how to read.
    UnsupportedVersion { found: u16, supported: u16 },
    /// An enum tag outside the set this build defines.
    BadTag { what: &'static str, tag: u64 },
    /// A length prefix that could not possibly be satisfied.
    TooLong { what: &'static str, len: u64, remaining: usize },
    /// Bytes left over after the value was fully decoded, which means the
    /// reader and writer disagree about the layout.
    TrailingBytes { count: usize },
    /// Text that was not valid UTF-8.
    BadText,
    /// A recipe did not check out: the blob was mangled, or the sampler in
    /// this build does not reproduce what the one that wrote it produced.
    /// Either way the client should ask again with `allow_recipes: false`.
    BadRecipe { what: &'static str },
    Io(String),
}

impl std::fmt::Display for WireError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WireError::Truncated { what, need, have } => {
                write!(f, "truncated reading {what}: needed {need} bytes, {have} left")
            }
            WireError::BadMagic => write!(f, "not a world file (bad magic)"),
            WireError::UnsupportedVersion { found, supported } => {
                write!(f, "world file is format version {found}, this build reads {supported}")
            }
            WireError::BadTag { what, tag } => write!(f, "unknown {what} tag {tag}"),
            WireError::TooLong { what, len, remaining } => {
                write!(f, "{what} claims {len} elements with {remaining} bytes left")
            }
            WireError::TrailingBytes { count } => {
                write!(f, "{count} bytes left over after decoding")
            }
            WireError::BadText => write!(f, "invalid UTF-8"),
            WireError::BadRecipe { what } => write!(f, "recipe failed its {what} check"),
            WireError::Io(e) => write!(f, "io: {e}"),
        }
    }
}

impl std::error::Error for WireError {}

pub type Result<T> = std::result::Result<T, WireError>;

// ---------------------------------------------------------------------------
// writing
// ---------------------------------------------------------------------------

#[derive(Default)]
pub struct Writer {
    buf: Vec<u8>,
}

impl Writer {
    pub fn new() -> Writer {
        Writer { buf: Vec::new() }
    }

    /// Magic and version. Call once, first.
    pub fn header(&mut self) {
        self.buf.extend_from_slice(&MAGIC);
        self.u16(FORMAT_VERSION);
    }

    #[inline]
    pub fn u8(&mut self, v: u8) {
        self.buf.push(v);
    }

    #[inline]
    pub fn bool(&mut self, v: bool) {
        self.buf.push(v as u8);
    }

    #[inline]
    pub fn u16(&mut self, v: u16) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    #[inline]
    pub fn u32(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    #[inline]
    pub fn u64(&mut self, v: u64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    #[inline]
    pub fn u128(&mut self, v: u128) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    /// By bits, so `-0.0` stays `-0.0` and a NaN keeps its payload.
    #[inline]
    pub fn f64(&mut self, v: f64) {
        self.buf.extend_from_slice(&v.to_bits().to_le_bytes());
    }

    #[inline]
    pub fn vec3(&mut self, v: Vec3) {
        self.f64(v.x);
        self.f64(v.y);
        self.f64(v.z);
    }

    #[inline]
    pub fn quat(&mut self, q: Quat) {
        self.f64(q.w);
        self.vec3(q.v);
    }

    /// A length prefix. Every variable-length thing is written as `seq(n)`
    /// followed by `n` items, so the reader can bound its allocation.
    #[inline]
    pub fn seq(&mut self, n: usize) {
        self.u32(n as u32);
    }

    pub fn str(&mut self, s: &str) {
        self.seq(s.len());
        self.buf.extend_from_slice(s.as_bytes());
    }

    pub fn bytes(&mut self, b: &[u8]) {
        self.seq(b.len());
        self.buf.extend_from_slice(b);
    }

    pub fn len(&self) -> usize {
        self.buf.len()
    }

    pub fn is_empty(&self) -> bool {
        self.buf.is_empty()
    }

    pub fn finish(self) -> Vec<u8> {
        self.buf
    }
}

// ---------------------------------------------------------------------------
// reading
// ---------------------------------------------------------------------------

pub struct Reader<'a> {
    buf: &'a [u8],
    at: usize,
}

impl<'a> Reader<'a> {
    pub fn new(buf: &'a [u8]) -> Reader<'a> {
        Reader { buf, at: 0 }
    }

    #[inline]
    pub fn remaining(&self) -> usize {
        self.buf.len() - self.at
    }

    #[inline]
    fn take(&mut self, what: &'static str, n: usize) -> Result<&'a [u8]> {
        if self.remaining() < n {
            return Err(WireError::Truncated { what, need: n, have: self.remaining() });
        }
        let s = &self.buf[self.at..self.at + n];
        self.at += n;
        Ok(s)
    }

    /// Magic and version. Returns the version so a caller can migrate.
    pub fn header(&mut self) -> Result<u16> {
        let m = self.take("magic", 4)?;
        if m != MAGIC {
            return Err(WireError::BadMagic);
        }
        let v = self.u16()?;
        if v != FORMAT_VERSION {
            return Err(WireError::UnsupportedVersion { found: v, supported: FORMAT_VERSION });
        }
        Ok(v)
    }

    #[inline]
    pub fn u8(&mut self) -> Result<u8> {
        Ok(self.take("u8", 1)?[0])
    }

    #[inline]
    pub fn bool(&mut self) -> Result<bool> {
        Ok(self.u8()? != 0)
    }

    #[inline]
    pub fn u16(&mut self) -> Result<u16> {
        let b = self.take("u16", 2)?;
        Ok(u16::from_le_bytes([b[0], b[1]]))
    }

    #[inline]
    pub fn u32(&mut self) -> Result<u32> {
        let b = self.take("u32", 4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    #[inline]
    pub fn u64(&mut self) -> Result<u64> {
        let b = self.take("u64", 8)?;
        let mut a = [0u8; 8];
        a.copy_from_slice(b);
        Ok(u64::from_le_bytes(a))
    }

    #[inline]
    pub fn u128(&mut self) -> Result<u128> {
        let b = self.take("u128", 16)?;
        let mut a = [0u8; 16];
        a.copy_from_slice(b);
        Ok(u128::from_le_bytes(a))
    }

    #[inline]
    pub fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_bits(self.u64()?))
    }

    #[inline]
    pub fn vec3(&mut self) -> Result<Vec3> {
        Ok(Vec3 { x: self.f64()?, y: self.f64()?, z: self.f64()? })
    }

    #[inline]
    pub fn quat(&mut self) -> Result<Quat> {
        Ok(Quat { w: self.f64()?, v: self.vec3()? })
    }

    /// A length prefix, validated before anything is allocated.
    ///
    /// `min_item_bytes` is the smallest number of bytes one element can
    /// possibly occupy. A claim that cannot be satisfied by the bytes actually
    /// remaining is rejected here rather than after a several-gigabyte
    /// allocation — which is the difference between a parse error and an abort.
    pub fn seq(&mut self, what: &'static str, min_item_bytes: usize) -> Result<usize> {
        let n = self.u32()? as usize;
        if n > MAX_SEQUENCE {
            return Err(WireError::TooLong {
                what,
                len: n as u64,
                remaining: self.remaining(),
            });
        }
        let need = n.saturating_mul(min_item_bytes.max(1));
        if need > self.remaining() {
            return Err(WireError::TooLong { what, len: n as u64, remaining: self.remaining() });
        }
        Ok(n)
    }

    pub fn str(&mut self) -> Result<String> {
        let n = self.seq("string", 1)?;
        let b = self.take("string body", n)?;
        std::str::from_utf8(b).map(|s| s.to_string()).map_err(|_| WireError::BadText)
    }

    pub fn bytes(&mut self) -> Result<Vec<u8>> {
        let n = self.seq("bytes", 1)?;
        Ok(self.take("bytes body", n)?.to_vec())
    }

    /// Decode an enum from a `u8` tag, with the valid range checked.
    pub fn tag(&mut self, what: &'static str, count: u8) -> Result<u8> {
        let t = self.u8()?;
        if t >= count {
            return Err(WireError::BadTag { what, tag: t as u64 });
        }
        Ok(t)
    }

    /// Assert the buffer was fully consumed. A disagreement between reader and
    /// writer usually shows up here rather than as bad data.
    pub fn finish(&self) -> Result<()> {
        if self.remaining() != 0 {
            return Err(WireError::TrailingBytes { count: self.remaining() });
        }
        Ok(())
    }
}
