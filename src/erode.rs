//! What a deviation's own physics does to it.
//!
//! `docs/PLAY.md` §5. "Detail exists where something is happening" is
//! implemented as a binary — anything touched is pinned and persisted, for
//! ever — and §5.1 says that at play resolution the binary is wrong, because it
//! makes a footprint as permanent as a felled trunk:
//!
//! > **The fix is to make the exemption a decay.** A deviation carries an
//! > amplitude, and the amplitude falls at a rate the physics sets. When it
//! > drops below the resolution anything could observe it at, the deviation is
//! > dropped and the node goes back to being purely regenerable.
//!
//! So forgetting stops being the deletion of a memory and becomes the deviation
//! reaching zero.
//!
//! # The rate is one expression
//!
//! §5.2 sets the test this module has to pass: *can one expression produce the
//! granite case and the wet-sand case with only material and flux differing? If
//! it needs a per-material correction, it is a table wearing a derivation's
//! clothes and it has failed.*
//!
//! It is one expression, and it is three lines:
//!
//! ```text
//!     tau = 1/2 rho v^2               what the flux presses on the surface, Pa
//!     e   = (tau - sigma_c)+ v / sigma_c    how fast it strips it, m/s
//!     k   = e / w                     and a feature of half-width w relaxes so
//! ```
//!
//! **The first line is a stress and the second is a threshold**, and the
//! threshold is the whole reason a scar in granite is not gone by lunchtime.
//! A flux presses on a surface with its own dynamic pressure; a grain is held
//! by `sigma_c`, the stress it takes to detach one. Below that nothing moves at
//! all — which is not a modelling convenience but the observed behaviour of
//! every real surface, and it is why wind does not erode rock. Above it, the
//! excess drives transport, and the stripping velocity is the excess as a
//! fraction of what holds the material, carried past at the speed of the flow.
//!
//! The third line is the geometry, and it is what makes a notch outlast a
//! scratch rather than the other way round. Nothing is stripped from a flat
//! plain, because transport is downslope and a plain has no slope; what a
//! feature loses is proportional to its own gradient, which for an amplitude
//! `h` over a half-width `w` is `h/w`. So `dh/dt = -e h / w`, the amplitude
//! decays exponentially, and the rate has no `h` in it — a deep scratch and a
//! shallow one of the same width take the same time to go, which is the right
//! behaviour and was not put in by hand.
//!
//! **Only material and flux differ**, which is the honesty test. `sigma_c` is
//! [`crate::material::Material::strength_of`] at the size of the material's own
//! worst flaw — the stress it takes to detach one grain rather than to break
//! the whole piece — and `rho` and `v` are the fluid and the flow the node is
//! actually in, both of which `morph::Environment` already measures rather than
//! being told.
//!
//! # What this cannot yet do, measured
//!
//! §5.2 promises wet sand at the tide line in hours and a scar in granite in
//! millennia from these three lines. **Granite it gives; sand it does not**, and
//! the reason is not the expression:
//!
//! ```text
//!   silica at 1600 kg/m^3, deposited     grain-scale strength 3.26e7 Pa
//!   silica at 1400 kg/m^3, deposited                          2.85e7 Pa
//!   granite, frozen at 1e-4 K/s                               9.72e7 Pa
//!   water at 3 m/s                       presses               4.5e3 Pa
//!   air at 25 m/s                                              3.8e2 Pa
//! ```
//!
//! A poured pile of dry sand comes out of `Material::measured` at 33 MPa in
//! tension, four orders above anything a flux can press with, because the
//! engine has **no representation for an uncemented aggregate**: every solid it
//! measures is a bonded one, and `Formation::Deposited` scales strength by the
//! packing squared — Gibson and Ashby's cellular solid, which assumes the cell
//! walls are joined. A sandpile's grains are not joined; what holds a surface
//! grain there is its own weight, and what holds a sandstone's is a cemented
//! neck whose area nothing in the engine states.
//!
//! So the threshold is right and the number going into it is wrong, for one
//! class of material. Finishing it needed a derived strength for a granular
//! aggregate — zero at zero cementation, rising from the jamming point — which
//! was a `PHYSICS.md` decision. Phase 5 took it: see [`loose_grain_threshold`],
//! [`capillary_cohesion`], and `World::is_loose` for which ground is loose.
//!
//! # What may decay, and what may not
//!
//! §5.8: **a deviation may be dropped only if dropping it leaves the conserved
//! tuple unchanged.** Forgetting is not permitted to be a source or a sink. So a
//! deviation records the mass that left the node when it was made, and one that
//! took mass away cannot decay at all — dropping it would put the mass back.
//! That is [`Deviation::conservative`], it is a measurement over the conserved
//! set rather than a judgement about which edits matter, and it is what
//! separates a squiggle drawn in sand from a channel cut through it without
//! either of them having been tagged when it was made.

use crate::math::Vec3;

/// Something that happened to a surface and is not in its derived baseline.
///
/// **A field rather than an edit list**, which is §5.7's decision for terrain:
/// ground is never eventful one grain at a time, and a field superposes for
/// free, so a thousand footprints are a thousand of these summed and cost one
/// evaluation each rather than a thousand stored heightmaps.
///
/// Position is in the patch's own normalised coordinates, which run from -1 to
/// 1 across it, so a deviation means the same place when the patch is drawn at
/// any level of detail. Everything else is in metres and kilograms.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Deviation {
    /// Where on the patch, -1 to 1 across it.
    pub x: f64,
    pub y: f64,
    /// Half-width of the feature, metres. **Its own geometry**, and the third
    /// input to the rate — a sharp narrow notch has a steeper gradient than a
    /// broad shallow dish and goes faster.
    pub span: f64,
    /// How far it moves the surface, metres. Negative is a cut.
    pub amplitude: f64,
    /// How much mass left the node when this was made, kg.
    ///
    /// Zero for a reshaping — a finger drawn through sand pushes grains aside
    /// and the node still holds them. Non-zero for a cut that carried material
    /// away, and that is what makes it permanent: see [`Self::conservative`].
    pub moved: f64,
}

impl Deviation {
    /// A reshaping: the surface moved and the node kept everything it had.
    pub fn reshaped(x: f64, y: f64, span: f64, amplitude: f64) -> Deviation {
        Deviation { x, y, span: span.max(1e-30), amplitude, moved: 0.0 }
    }

    /// A cut: this much mass left the node when the feature was made.
    pub fn cut(x: f64, y: f64, span: f64, amplitude: f64, moved: f64) -> Deviation {
        Deviation { x, y, span: span.max(1e-30), amplitude, moved: moved.abs() }
    }

    /// Whether dropping this would leave the node's conserved tuple unchanged.
    ///
    /// §5.8's rule, and the reason the exploit it describes does not work: a
    /// deviation that walked off with mass may not be forgotten however long
    /// nobody looks at it, because forgetting it would put the mass back while
    /// it is still somewhere else.
    ///
    /// The threshold is the coarse level's own representational resolution,
    /// which §5.3 says it must be and which is the shape
    /// `tree::IDEMPOTENT_TOLERANCE` already has, rather than a number chosen
    /// here.
    pub fn conservative(&self, node_mass: f64) -> bool {
        self.moved <= crate::tree::IDEMPOTENT_TOLERANCE * node_mass.abs().max(1e-300)
    }

    /// How much this raises the surface at a point, metres.
    ///
    /// A Gaussian of its own half-width, which is what makes the field
    /// superpose: a sum of these is a surface, and there is no seam where two
    /// of them overlap.
    pub fn height_at(&self, x: f64, y: f64, half_side: f64) -> f64 {
        let w = (self.span / half_side.max(1e-30)).max(1e-12);
        let dx = (x - self.x) / w;
        let dy = (y - self.y) / w;
        let r2 = dx * dx + dy * dy;
        if r2 > 32.0 {
            return 0.0;
        }
        self.amplitude * (-r2).exp()
    }

    /// Volume of the dimple this makes, m^3. A Gaussian's integral.
    pub fn volume(&self) -> f64 {
        std::f64::consts::PI * self.span * self.span * self.amplitude
    }
}

/// How fast a surface feature of this half-width relaxes, per second.
///
/// The three lines of this module's documentation, in that order.
///
/// * `cohesion` — the stress it takes to detach a grain of the material, Pa.
/// * `fluid_density`, `flow_speed` — what is going past, measured.
/// * `span` — the feature's own half-width, metres.
///
/// Zero when the flux cannot reach the threshold, when there is no flux, or
/// when the material has no strength left to have a threshold: a feature in
/// vacuum lasts, and so does one in something nothing can lift.
///
/// A material with **no** strength at all — a surface at or past its own
/// softening point, which `Material::strength_at` already reports — is not a
/// division by zero but a surface that flows, and the rate is capped at the
/// flux's own speed over the feature, because nothing relaxes faster than the
/// thing moving it goes.
pub fn relaxation_rate(cohesion: f64, fluid_density: f64, flow_speed: f64, span: f64) -> f64 {
    if !(span > 0.0) || !(fluid_density > 0.0) || !(flow_speed > 0.0) {
        return 0.0;
    }
    let pressing = 0.5 * fluid_density * flow_speed * flow_speed;
    let ceiling = flow_speed / span;
    if !(cohesion > 0.0) {
        return ceiling;
    }
    let excess = pressing - cohesion;
    if excess <= 0.0 {
        return 0.0;
    }
    let stripping = excess * flow_speed / cohesion;
    let rate = (stripping / span).min(ceiling);
    if rate.is_finite() { rate.max(0.0) } else { 0.0 }
}

/// Where a point on a sphere lands in a patch's normalised coordinates.
///
/// A convenience for the callers that have a direction rather than a pair of
/// numbers — an actor standing somewhere and drawing in the sand has one.
pub fn local_of(patch: &crate::recipe::Tiled, dir: Vec3) -> Option<(f64, f64)> {
    patch.local_of_direction(dir)
}

/// Close packing of equal spheres, the densest there is: `pi / (3 sqrt 2)`.
/// A universal of geometry, and the reference a looser packing is measured
/// from.
pub const CLOSE_PACKING: f64 = 0.740_480_489_693_061;

/// How steep a pocket a grain sits in at the surface of a loose pile: the
/// tangent of the angle it has to be tipped through to roll out.
///
/// A grain resting on three others that touch each other — a close-packed
/// pocket — tips out over the edge between two of its supports at
/// `tan a = 1 / (2 sqrt 2)`, 19.5 degrees, which is geometry and nothing else.
/// A loose pile's supports do not touch: at packing `phi` the spacing between
/// them is the close-packed one opened by `(phi_close / phi)^(1/3)`, the grain
/// sits deeper between them, and the pocket is steeper. At random loose
/// packing that gives 22.5 degrees.
///
/// This is Phase 5's choice for what holds an uncemented grain, taken with
/// three alternatives on the table: the pocket of a random packing rather than
/// the Bowden–Tabor surface friction the contact code derives, because a
/// grain in a pile is held by the grains around it, not by rubbing.
pub fn pocket_friction(packing: f64) -> f64 {
    let open = (CLOSE_PACKING / packing.clamp(1e-3, CLOSE_PACKING)).cbrt();
    let height = (1.0 - open * open / 3.0).max(1e-9).sqrt();
    open / (2.0 * 3.0f64.sqrt() * height)
}

/// The stress a flux has to press with to lift a grain of an uncemented
/// aggregate, Pa, in the terms this module compares it in.
///
/// **What holds a loose grain is its own weight**, not a bond. A sphere of
/// diameter `D` presents `pi D^2 / 4` to a flow and weighs
/// `(rho_s - rho_f) g pi D^3 / 6` in the fluid it is in, and it leaves when the
/// flow's push tips it out of its pocket:
///
/// ```text
///     sigma_c = tan(a) (rho_s - rho_f) g (2 D / 3)
/// ```
///
/// Shields' threshold is the same balance written against a bed shear stress
/// rather than a dynamic pressure; in those terms this is a Shields number of
/// `(2/3) tan a`, 0.28 at random loose packing, where Shields' own curve is
/// quoted in shear-stress units and sits at 0.03 to 0.06. The two differ by
/// the ratio of the pressure at a grain to the stress on the bed, which is a
/// boundary layer this module does not model and which [`relaxation_rate`]
/// already takes the dynamic pressure for.
pub fn loose_grain_threshold(grain_density: f64, fluid_density: f64, g: f64, grain: f64, packing: f64) -> f64 {
    let buoyant = (grain_density - fluid_density).max(0.0);
    (pocket_friction(packing) * buoyant * g.abs() * 2.0 / 3.0 * grain.max(0.0)).max(0.0)
}

/// What a liquid's menisci hold a damp aggregate together with, Pa.
///
/// Rumpf's relation for a packing held by bridges between its grains,
/// `sigma = (9/8) (phi / (1 - phi)) F / D^2`, with the bridge force between two
/// wetted spheres of radius `D/2` taken as `2 pi gamma (D/2)`:
///
/// ```text
///     sigma = (9/8) (phi / (1 - phi)) pi gamma / D
/// ```
///
/// Why a sandcastle stands and a heap of dry sand does not. It holds only
/// where the pores hold both liquid and air: fully submerged there are no
/// menisci and nothing here applies, which is why the same squiggle that
/// holds its shape on a damp beach at low tide is gone when the water is over
/// it. The caller decides which; see `World::weather`.
pub fn capillary_cohesion(surface_tension: f64, grain: f64, packing: f64) -> f64 {
    if !(grain > 0.0) || !(surface_tension > 0.0) {
        return 0.0;
    }
    let phi = packing.clamp(1e-3, 0.999);
    9.0 / 8.0 * phi / (1.0 - phi) * std::f64::consts::PI * surface_tension / grain
}

/// A liquid's surface tension, J/m^2, from what holds it together.
///
/// The rule `Material::of` uses for a solid's surface energy — a quarter of the
/// energy per unit spread over the area one unit presents — applied to the
/// energy that actually holds a liquid, which is its vaporisation enthalpy
/// rather than the bonds inside its molecules (see `eos.rs`, where the same
/// distinction sets a liquid's stiffness).
pub fn surface_tension(props: &crate::chem::Properties) -> f64 {
    let l = crate::chem::react::heat_of_vaporisation(props);
    let m = props.unit_mass;
    let rho = props.density;
    if !(l > 0.0) || !(m > 0.0) || !(rho > 0.0) {
        return 0.0;
    }
    let area = (m / rho).powf(2.0 / 3.0);
    l * m / (4.0 * area)
}
