//! Chemistry: what matter *is*, as opposed to what elements it contains.
//!
//! See `elements` for why this coexists with `units::CoarseElement` rather than
//! replacing it, `arrange` for how a substance is described, `analyse` for
//! how its properties are derived rather than looked up, and `registry` for
//! where the derivation is kept.

pub mod analyse;
pub mod arrange;
pub mod elements;
pub mod geometry;
pub mod react;
pub mod registry;

pub use analyse::{analyse, Confidence, Illegal, Properties};
pub use arrange::{Arrangement, Bond, Formula, Lattice, Order};
pub use elements::Element;
pub use react::{equilibrate, react, ReactionReport};
pub use registry::{Mixture, Phase, Pool, Provenance, Registry, Substance, SubstanceId};
