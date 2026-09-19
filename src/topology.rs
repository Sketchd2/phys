//! Joints: what makes a heap of parts into a structure.
//!
//! # The gap this closes
//!
//! Until now a generated tree was geometry that happened to hold still.
//! Materialise it and hand the parts to the molecular dynamics solver and the
//! trunk falls apart, because nothing in a `Body` says that this segment is
//! attached to that one. The conserved tuple pins mass and momentum; the
//! developmental state pins shape; neither pins *cohesion*.
//!
//! A `Topology` is the missing piece: an explicit list of joints, each with a
//! cross-section, a material and a remaining integrity. It costs one array
//! alongside the bodies and it is regenerated from the program exactly as the
//! geometry is, so it is as free as everything else the engine chooses not to
//! store.
//!
//! # Why the support graph is a tree, and why that matters enormously
//!
//! Every structure the engine generates has a *load path* that is a tree: a
//! branch hangs off a branch, a floor stands on the floor below, a course of
//! bricks rests on the course beneath. A real framed building is of course a
//! redundant lattice, and getting its true internal forces needs a stiffness
//! matrix and a sparse solve — thousands of unknowns, iterative, and awkward to
//! budget for at 20 frames a second.
//!
//! On a tree it is exact in one pass. Accumulate force and moment from the
//! leaves inward, and every joint's internal load is known in O(n) with no
//! solve at all. That is the difference between structural failure being an
//! occasional expensive event and it being something the engine can afford to
//! check on every structure, every frame.
//!
//! The approximation is stated rather than hidden: for a redundant structure
//! this is the load path under gravity with no force sharing between
//! alternative routes, which is conservative — it over-predicts the load on the
//! nominal path and therefore fails early rather than late.

use crate::math::Vec3;
use crate::morph::{Skeleton, NO_SUPPORT};

/// Re-exported so `topology::Material` keeps working. The type moved to its own
/// module when it stopped being a property of structures and became something
/// measured off any matter at all — see `docs/PLAY.md` D13 and `material.rs`.
pub use crate::material::Material;

/// One joint between two parts.
#[derive(Debug, Clone, Copy)]
pub struct Joint {
    /// The supported part.
    pub child: u32,
    /// The supporting part, or [`NO_SUPPORT`] for a ground anchor.
    pub parent: u32,
    /// Where the joint is. Bending stress peaks here, so this is where it goes.
    pub at: Vec3,
    /// Cross-sectional radius of the joint, metres.
    pub radius: f64,
    /// Fraction of nominal strength remaining, 0..1. Reduced by heat, decay and
    /// previous damage; a joint at zero has already failed.
    pub integrity: f64,
}

impl Joint {
    /// Section modulus `I/c = pi r^3 / 4`, in m^3. Divide a bending moment by
    /// this to get the peak fibre stress.
    #[inline]
    pub fn section_modulus(&self) -> f64 {
        std::f64::consts::PI * self.radius.powi(3) / 4.0
    }

    #[inline]
    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.radius * self.radius
    }
}

/// The joints of one structure, in the same index space as its bodies.
#[derive(Debug, Clone, Default)]
pub struct Topology {
    pub joints: Vec<Joint>,
    /// Supporting part per body, parallel to the body list. `NO_SUPPORT` means
    /// the part is anchored, or is loose matter with no structural role.
    pub support: Vec<u32>,
    /// Program-stable name per body, for naming a failure in an event.
    pub site: Vec<u32>,
    pub base: Vec<Vec3>,
    pub tip: Vec<Vec3>,
    pub material: Material,
    /// Connections *beyond* the support forest: bracing, ties, redundant load
    /// paths. Their presence is what decides whether the structure is
    /// statically determinate, and therefore which solver applies.
    pub ties: Vec<Tie>,
}

/// One member of an explicitly-described structure.
#[derive(Debug, Clone, Copy)]
pub struct Member {
    pub base: Vec3,
    pub tip: Vec3,
    pub radius: f64,
    /// Index of the member that supports this one, or [`NO_SUPPORT`] for a
    /// ground anchor.
    pub support: u32,
}

impl Member {
    pub fn new(base: Vec3, tip: Vec3, radius: f64, support: u32) -> Member {
        Member { base, tip, radius, support }
    }
    pub fn anchored(base: Vec3, tip: Vec3, radius: f64) -> Member {
        Member { base, tip, radius, support: NO_SUPPORT }
    }
}

/// A redundant connection between two parts.
#[derive(Debug, Clone, Copy)]
pub struct Tie {
    pub a: u32,
    pub b: u32,
    /// Cross-sectional area of the tie, m^2.
    pub area: f64,
    pub integrity: f64,
}

impl Topology {
    /// Build the joint list from a generated skeleton.
    ///
    /// `scale` converts the skeleton's normalised units into metres — the same
    /// factor `sample_structured` applies to the positions, so the joints land
    /// exactly on the parts.
    /// `shift` and `scale` are the same centring and scaling `sample_structured`
    /// applied to the positions, and `radii` the density-corrected member radii.
    /// Passing anything else puts the joints somewhere the parts are not.
    pub fn from_skeleton(
        skel: &Skeleton,
        material: Material,
        shift: Vec3,
        scale: f64,
        radii: &[f64],
        parts: usize,
    ) -> Topology {
        let n = skel.len();
        let place = |p: Vec3| (p - shift).scale(scale);
        let mut joints = Vec::with_capacity(n);
        for i in 0..n {
            joints.push(Joint {
                child: i as u32,
                parent: skel.support[i],
                at: place(skel.base[i]),
                // The member's own section, unless the generator stated a
                // narrower joint — `docs/PLAY.md` D15's glued seam. Scaled
                // alongside everything else, because it is stated in the same
                // units the skeleton is.
                radius: match skel.joint_radius.get(i).copied().unwrap_or(0.0) {
                    r if r > 0.0 => r * scale,
                    _ => radii.get(i).copied().unwrap_or(skel.radius[i] * scale),
                },
                integrity: 1.0,
            });
        }
        // The body list is longer than the skeleton whenever the node also
        // holds unstructured matter, so *every* parallel array is padded to the
        // same length — including the joints. Leaving `joints` short while the
        // others were padded meant any index derived from a body number could
        // walk off the end of it, which is exactly what an insult entering at a
        // litter particle did.
        joints.resize(
            parts,
            Joint {
                child: 0,
                parent: NO_SUPPORT,
                at: Vec3::ZERO,
                // Zero radius marks a part with no structural role, which is
                // also how the renderer and the failure check identify litter.
                radius: 0.0,
                integrity: 0.0,
            },
        );
        let mut support = skel.support.clone();
        let mut site = skel.site.clone();
        let mut base: Vec<Vec3> = skel.base.iter().map(|p| place(*p)).collect();
        let mut tip: Vec<Vec3> = skel.tip.iter().map(|p| place(*p)).collect();
        support.resize(parts, NO_SUPPORT);
        site.resize(parts, NO_SUPPORT);
        base.resize(parts, Vec3::ZERO);
        tip.resize(parts, Vec3::ZERO);
        Topology {
            joints,
            support,
            site,
            base,
            tip,
            material,
            ties: skel
                .ties
                .iter()
                .map(|&(a, b, fraction)| {
                    // Sized from the members actually being joined, using the
                    // density-corrected radii, so a tie is always in proportion
                    // to its neighbours however the structure has been scaled.
                    let ra = radii.get(a as usize).copied().unwrap_or(0.0);
                    let rb = radii.get(b as usize).copied().unwrap_or(0.0);
                    let r = ra.min(rb);
                    Tie {
                        a,
                        b,
                        area: fraction * std::f64::consts::PI * r * r,
                        integrity: 1.0,
                    }
                })
                .collect(),
        }
    }

    /// Build a topology directly from an explicit list of members and ties.
    ///
    /// The generated structures in `morph.rs` are one source of geometry, not
    /// the only one. Anything that can describe itself as members with a
    /// support relation — a truss, a bridge, a scaffold, a machine, a skeleton
    /// imported from elsewhere — gets the same analysis.
    pub fn from_parts(
        members: &[Member],
        ties: &[(u32, u32, f64)],
        material: Material,
    ) -> Topology {
        let n = members.len();
        let mut joints = Vec::with_capacity(n);
        let (mut support, mut site, mut base, mut tip) = (
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
        );
        for (i, m) in members.iter().enumerate() {
            joints.push(Joint {
                child: i as u32,
                parent: m.support,
                at: m.base,
                radius: m.radius,
                integrity: 1.0,
            });
            support.push(m.support);
            site.push(i as u32);
            base.push(m.base);
            tip.push(m.tip);
        }
        Topology {
            joints,
            support,
            site,
            base,
            tip,
            material,
            ties: ties
                .iter()
                .map(|&(a, b, area)| Tie { a, b, area, integrity: 1.0 })
                .collect(),
        }
    }

    /// Is the load path a forest — every part supported by at most one other,
    /// with no alternative routes?
    ///
    /// This is the question that decides which solver runs. A forest is exactly
    /// solvable in one pass; anything else is statically indeterminate and the
    /// internal forces depend on relative stiffness, which needs a solve.
    pub fn is_determinate(&self) -> bool {
        self.ties.iter().all(|t| t.integrity <= 0.0)
    }

    pub fn is_empty(&self) -> bool {
        self.joints.is_empty()
    }

    /// Parts that are held by nothing — the ones that fall.
    pub fn loose(&self) -> impl Iterator<Item = u32> + '_ {
        self.support
            .iter()
            .enumerate()
            .filter(|(_, s)| **s == NO_SUPPORT)
            .map(|(i, _)| i as u32)
    }

    /// Bytes held. Compared against the geometry it makes coherent.
    pub fn bytes(&self) -> usize {
        self.joints.len() * std::mem::size_of::<Joint>()
            + (self.support.len() + self.site.len()) * 4
            + (self.base.len() + self.tip.len()) * std::mem::size_of::<Vec3>()
    }
}
