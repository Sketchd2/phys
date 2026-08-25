//! Bonds: what makes a heap of parts into a structure.
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

/// What a joint is made of. Strength here is the modulus of rupture: the
/// bending stress at which the section fails.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Material {
    #[default]
    /// Living wood, wet. Strong in bending, and the reason a branch bends a
    /// long way before it goes.
    GreenWood,
    /// Coral skeleton. Stiff and brittle.
    Aragonite,
    /// Reinforced concrete and steel.
    ReinforcedFrame,
    /// Mortared brick or stone: strong in compression, nearly useless in
    /// tension, which is why walls fall over rather than snapping.
    Masonry,
}

impl Material {
    /// Modulus of rupture at room temperature, pascals.
    pub fn strength(self) -> f64 {
        match self {
            Material::GreenWood => 45.0e6,
            Material::Aragonite => 12.0e6,
            Material::ReinforcedFrame => 180.0e6,
            Material::Masonry => 2.0e6,
        }
    }

    /// Ratio of tensile to compressive strength. Masonry's asymmetry is the
    /// whole reason its failure mode differs from wood's.
    pub fn tensile_ratio(self) -> f64 {
        match self {
            Material::GreenWood => 1.0,
            Material::Aragonite => 0.4,
            Material::ReinforcedFrame => 0.9,
            Material::Masonry => 0.05,
        }
    }

    /// Young's modulus, pascals. Used for deflection, not for failure.
    pub fn stiffness(self) -> f64 {
        match self {
            Material::GreenWood => 10.0e9,
            Material::Aragonite => 60.0e9,
            Material::ReinforcedFrame => 30.0e9,
            Material::Masonry => 15.0e9,
        }
    }

    /// Temperature at which the material starts losing strength, and the
    /// temperature at which it has none left.
    pub fn thermal_limits(self) -> (f64, f64) {
        match self {
            // Wood holds its strength until water is driven off, then pyrolyses.
            Material::GreenWood => (373.0, 575.0),
            // Calcium carbonate calcines.
            Material::Aragonite => (600.0, 1100.0),
            // Steel loses half its yield by 800 K; concrete spalls.
            Material::ReinforcedFrame => (600.0, 1000.0),
            Material::Masonry => (900.0, 1500.0),
        }
    }

    /// Fraction of nominal strength remaining at a given temperature.
    pub fn strength_at(self, temperature: f64) -> f64 {
        let (onset, gone) = self.thermal_limits();
        if temperature <= onset {
            1.0
        } else if temperature >= gone {
            0.0
        } else {
            1.0 - (temperature - onset) / (gone - onset)
        }
    }

    /// Enthalpy needed to destroy a kilogram of this material outright —
    /// vaporising the water in it and pyrolysing what is left. This is what a
    /// lightning channel has to supply to blow a tree apart.
    pub fn destruction_enthalpy(self) -> f64 {
        match self {
            // Mostly the latent heat of the sap: ~40% water at 2.26 MJ/kg,
            // plus heating and pyrolysis.
            Material::GreenWood => 1.4e6,
            Material::Aragonite => 1.8e6,
            Material::ReinforcedFrame => 1.2e6,
            Material::Masonry => 1.0e6,
        }
    }
}

/// One joint between two parts.
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    /// The supported part.
    pub child: u32,
    /// The supporting part, or [`NO_SUPPORT`] for a ground anchor.
    pub parent: u32,
    /// Where the joint is. Bending stress peaks here, so this is where it goes.
    pub at: Vec3,
    /// Cross-sectional radius of the joint, metres.
    pub radius: f64,
    /// Fraction of nominal strength remaining, 0..1. Reduced by heat, decay and
    /// previous damage; a bond at zero has already failed.
    pub integrity: f64,
}

impl Bond {
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
    pub bonds: Vec<Bond>,
    /// Supporting part per body, parallel to the body list. `NO_SUPPORT` means
    /// the part is anchored, or is loose matter with no structural role.
    pub support: Vec<u32>,
    /// Program-stable name per body, for naming a failure in an event.
    pub site: Vec<u32>,
    pub base: Vec<Vec3>,
    pub tip: Vec<Vec3>,
    pub material: Material,
}

impl Topology {
    /// Build the joint list from a generated skeleton.
    ///
    /// `scale` converts the skeleton's normalised units into metres — the same
    /// factor `prolong_structured` applies to the positions, so the joints land
    /// exactly on the parts.
    /// `shift` and `scale` are the same centring and scaling `prolong_structured`
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
        let mut bonds = Vec::with_capacity(n);
        for i in 0..n {
            bonds.push(Bond {
                child: i as u32,
                parent: skel.support[i],
                at: place(skel.base[i]),
                radius: radii.get(i).copied().unwrap_or(skel.radius[i] * scale),
                integrity: 1.0,
            });
        }
        // The body list is longer than the skeleton whenever the node also
        // holds unstructured matter, so *every* parallel array is padded to the
        // same length — including the bonds. Leaving `bonds` short while the
        // others were padded meant any index derived from a body number could
        // walk off the end of it, which is exactly what an insult entering at a
        // litter particle did.
        bonds.resize(
            parts,
            Bond {
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
            bonds,
            support,
            site,
            base,
            tip,
            material,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty()
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
        self.bonds.len() * std::mem::size_of::<Bond>()
            + (self.support.len() + self.site.len()) * 4
            + (self.base.len() + self.tip.len()) * std::mem::size_of::<Vec3>()
    }
}
