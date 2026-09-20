//! A boundary crossing is the event — `docs/PLAY.md` D16.
//!
//! One mechanism answers three questions the tree had no answer to. **Inward**,
//! a crossing generates the detail about to be met. **Outward**, it re-homes to
//! the node above. **Into a sibling**, it re-homes sideways, with the parent
//! arbitrating. The rule is "re-home to whatever now contains you", and the
//! parent's own index — which holds its bodies *and* its promoted children in
//! one list — is what answers it.
//!
//! Until this existed, `reparent` was only ever called by hand, through
//! `Interaction::Rehome`. Measured before any of it: a rocket at escape
//! velocity leaves a 1 km forest node in 0.100 s and nothing notices; a hundred
//! seconds later it is 1.12x10^6 m away, still a child of a node that claims a
//! kilometre, and the spatial index has clamped it into a corner cell where it
//! can neither collide with nor exchange heat with anything.
//!
//! # What a node owns, and why it is not simply its radius
//!
//! The obvious reading of "leaving the patch's volume" is geometric: a node has
//! left when its bulk is wholly outside its parent's claimed sphere. Measured,
//! that cannot work, and the case that breaks it is every case that matters.
//! **A surface is exactly where a sphere's boundary is.** The rocket leaves the
//! forest at 1.1 km altitude, re-homes to the planet — and is immediately
//! 6372 km from a planet claiming 6371, so it is already wholly outside its new
//! parent and re-homes to the star on the next frame. `Tree::gravity_at` walks
//! ancestors, by design, so the planet becomes a sibling and its gravity
//! vanishes while the rocket is a kilometre above the ground.
//!
//! So a node owns the region where *it* is the thing that holds you, which is
//! the larger of two lengths:
//!
//! * the volume its matter occupies — `matter.radius`, the equivalent uniform
//!   sphere, which is exactly the radius `sample` draws inside; and
//! * its **Hill radius** about its own parent, `d·(m/3M)^(1/3)`, the distance
//!   at which its parent's gravity takes over.
//!
//! Nothing is supplied: the two masses and the separation are already on the
//! nodes. The atmosphere was considered first and is *inverted* for this
//! purpose — scale height is `kT/(μg)`, so a 500 m asteroid's is 2.15x10^5 km,
//! 430,000 times its own radius, precisely because its gravity is too weak to
//! hold gas at all, while Earth's is 8.5 km, 0.13% of its radius. It hands the
//! large multiple to the wrong body. The Hill radius has the right shape
//! measured on the same pair: **495 radii for the asteroid, 235 for Earth**,
//! and it stays out of the way where it must — 24.3 m for a 1 km forest patch
//! and 0.235 m for a 900 kg tree, so D6's walk-off-a-patch is decided by
//! geometry exactly as D6 says.
//!
//! An atmosphere needs no term of its own. If a planet ever holds a modelled
//! envelope, that envelope is matter in the node and `matter.radius` grows with
//! it, so the domain does too.
//!
//! # The band, and why leaving and arriving are not the same test
//!
//! Leaving asks whether the node is *wholly* outside; arriving asks whether it
//! is *wholly* inside. Between the two is a band — the thing straddles the
//! boundary — and in that band nothing happens at all. That asymmetry is the
//! whole of the hysteresis: a node sitting on a boundary neither leaves nor
//! arrives, so it cannot flip between two parents frame after frame, being
//! rekeyed and re-indexed each time. It is also what lets something *stand* on
//! a planet: resting on the surface, its bulk crosses the boundary, so it
//! straddles and stays.

use crate::ids::NodeIdx;
use crate::math::Vec3;

/// Where something of `radius`, sitting `offset` from the centre of a volume of
/// `bound`, stands with respect to that volume.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Standing {
    /// Wholly within it. The only verdict that lets something *arrive*.
    Inside,
    /// Across its surface — neither wholly in nor wholly out. Nothing happens.
    Straddling,
    /// Wholly beyond it. The only verdict that makes something *leave*.
    Outside,
}

/// Measure one thing against one volume.
///
/// Radii are taken into account on the thing's side, the same convention
/// `state::Spread` uses for its `furthest`: a body is not inside a volume its
/// own bulk sticks out of, and it has not left one it is still poking into.
pub fn standing(offset: Vec3, radius: f64, bound: f64) -> Standing {
    let d = offset.norm();
    let r = radius.max(0.0);
    if !(bound > 0.0) || !d.is_finite() {
        // A volume that claims nothing contains nothing. Reported as outside
        // rather than straddling, because something has to decide, and a node
        // with no extent is not a place to live.
        return Standing::Outside;
    }
    if d + r <= bound {
        Standing::Inside
    } else if d - r > bound {
        Standing::Outside
    } else {
        Standing::Straddling
    }
}

/// The distance at which a parent's gravity takes over from a child's.
///
/// `d·(m/3M)^(1/3)`, with `M` the mass of the parent *excluding* the child —
/// the same subtraction `Tree::gravity_at` makes, and for the same reason: a
/// node's matter counts its promoted children, so a child attracted by its own
/// mass would be attracted by itself.
///
/// Zero at the centre, which is right: a node sitting exactly on its parent's
/// centre of mass has no region where it wins. A parent lighter than its own
/// child is answered with the whole separation rather than with infinity —
/// nothing in the engine is supposed to reach that state, and returning a
/// number that poisons every comparison downstream is a worse answer than a
/// bounded one.
pub fn hill_radius(mass: f64, parent_mass: f64, separation: f64) -> f64 {
    let d = separation.max(0.0);
    if !(mass > 0.0) || !d.is_finite() {
        return 0.0;
    }
    let source = parent_mass - mass;
    if source <= 0.0 {
        return d;
    }
    d * (mass / (3.0 * source)).cbrt()
}

/// Which way a crossing went.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Out of the parent and into the node above it, because nothing the
    /// arbiter holds contains the crosser any more.
    Outward,
    /// Into something the arbiter holds — the creature walking from the forest
    /// into the desert, which should never become a direct child of the planet.
    Sideways,
}

/// One crossing, as it happened.
#[derive(Debug, Clone, Copy)]
pub struct Crossed {
    pub node: NodeIdx,
    /// The parent it left.
    pub from: NodeIdx,
    /// The parent it arrived at.
    pub to: NodeIdx,
    pub direction: Direction,
    /// Whether the destination had to be *generated* to be arrived at — the
    /// inward row of D16's table. A crossing into a place that was still only
    /// one of its parent's bodies promotes that body first.
    pub generated: bool,
}
