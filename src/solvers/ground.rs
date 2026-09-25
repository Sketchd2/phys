//! Ground: pieces of rock held together by gravity and kept apart by their own
//! stiffness.
//!
//! # Why ground is solved at all
//!
//! Phase 4 made ground a *recipe*: a tiled planet's faces, a patch's cells and
//! the substrate under them are members of a structure, and a structure's
//! members were not the solver's to move. That held for as long as nothing
//! asked a piece of ground to be anything but where the recipe drew it. Phase 5
//! promoted one — a face of a turning Earth, to hold the air over its sea — and
//! the face, now its own node, met the planet's interior as a separate object
//! and was thrown out at 76 m/s on the first contact. Nothing had ever said
//! what holds ground together, and the answer the owner gave for this phase is
//! the physical one: **its gravity, against its own stiffness**.
//!
//! # The law
//!
//! Every body feels gravity: the field from outside the node where that body
//! is (`Tree::gravity_at_point`), and the pull of the node's own contents on
//! each other. Neighbouring pieces are
//! joined by the rock's stiffness, a column of it between their centres:
//!
//! ```text
//!   k = E A / L
//! ```
//!
//! with `E` the Young's modulus `Material` derives for what the ground is made
//! of, `A` the face the two pieces share and `L` the distance between their
//! centres. That is the flat-faced form of the Hertz contact that was named for
//! this: two slabs sharing a face are a column in compression, and Hertz's
//! `sqrt(R) d^(3/2)` is what the same modulus gives for two curved bodies
//! touching at a point, which tiles are not.
//!
//! # Ground as drawn is at rest
//!
//! **Each piece carries the load its own weight puts on it where the recipe
//! draws it.** A planet's rock is under the lithostatic pressure of everything
//! above it; that stress is what holds a piece up against its gravity, and it
//! is present whether or not anything has moved. So the springs are unstretched
//! at the recipe's positions, and each piece is given, once, the support that
//! exactly balances the gravity it feels there (`preload`). The configuration
//! the recipe draws is then an equilibrium: ground nobody disturbs stays where
//! it is drawn and regenerates bit for bit, and a disturbance rings through the
//! springs at the rock's own sound speed. Nothing damps it, which is what rock
//! does over the time anything here watches it.

use crate::math::Vec3;
use crate::solvers::gravity::{GravityParams, Octree};
use crate::solvers::SolveReport;
use crate::state::Body;

/// A column of rock between two pieces of ground.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spring {
    pub a: u32,
    pub b: u32,
    /// Stiffness, N/m: `E A / L`.
    pub k: f64,
    /// Length at which it carries nothing, m: the pieces' distance as drawn.
    pub rest: f64,
}

/// What a node's ground is, derived on its first solve after it is drawn and
/// kept until it is drawn again: the springs, with the lengths they rest at,
/// and each piece's support in the node's own axes, so that it turns with it.
/// Not persisted — it is derived from the drawing, and a reload draws again.
#[derive(Debug, Clone, Default)]
pub struct Cache {
    pub springs: Vec<Spring>,
    /// Support per body, m/s^2, in the node's own axes; zero past the pieces.
    pub support: Vec<Vec3>,
    /// How many of the node's bodies are pieces of its ground.
    pub pieces: usize,
}

/// Everything a ground solve needs besides the bodies.
#[derive(Debug, Clone)]
pub struct Ground {
    pub springs: Vec<Spring>,
    /// The field from outside the node at each body, m/s^2, in the bodies'
    /// axes — at each, because a patch of a planet a continent across is not
    /// in one field.
    pub field: Vec<Vec3>,
    /// The support each body is given against its weight as drawn, m/s^2 —
    /// zero for anything that is not a piece of the ground.
    pub preload: Vec<Vec3>,
    /// Self-gravity among the pieces — the first `pieces` bodies — which pull
    /// each other as the points they are drawn as. Anything loose among them
    /// is in the smooth field `field` gives it instead: a point standing in for
    /// a continent of rock pulls what is close to it far too hard, 11.9 m/s^2
    /// on air over a face of an Earth where the whole planet gives 9.8.
    pub gravity: Option<GravityParams>,
    pub pieces: usize,
    /// The share of its weight the fluid around each body carries, within the
    /// step — the ambient density over its own — so that something at the
    /// density around it is weightless while it moves rather than only after.
    /// Pushed back after a step instead, a parcel of air ten metres over a sea
    /// fell 24 m in a 2.2 s step, found itself under water, and was pushed up
    /// by 1450 times its own weight. Zero for the pieces.
    pub lift: Vec<f64>,
    /// What else moves what is loose among the pieces, for a node inside a
    /// fluid that goes round with its planet — `None` where nothing does.
    pub swirl: Option<Swirl>,
    /// The ground's angular velocity, rad/s, in the bodies' axes: the support
    /// is a stress fixed in the ground, so within a step it turns with it.
    pub turn: Vec3,
}

/// The fluid around a node going round with the planet that describes it, and
/// the node's own frame going round in it.
///
/// **A node that goes round with its planet is not an inertial frame.** Its
/// centre accelerates toward the planet's axis, and what is loose in it — not
/// its pieces, whose support takes it in — is measured against that. And the
/// fluid around a loose thing is going round too, so its push is Archimedes in
/// an accelerating fluid, `rho V (a - g)`: the `a` part here, the `g` part in
/// `lift`. Air at the density around it then goes round with the ground
/// exactly.
///
/// **Both are laws of where the body is and when**, not numbers read once at
/// the start of a step. Held for a step, the pull toward the axis is an Euler
/// step on circular motion, which adds `(a dt)^2 / 2v` to the speed every
/// step: measured, air at rest over the sea of a turning Earth went east at
/// 0.0135 m/s an hour with nothing pushing it.
#[derive(Debug, Clone, Copy)]
pub struct Swirl {
    /// The fluid's angular velocity, rad/s, in the bodies' axes.
    pub spin: Vec3,
    /// The node's centre from the describing body's, m, at the step's start.
    pub centre: Vec3,
    /// The acceleration of the node's own frame, m/s^2, at the step's start.
    pub frame: Vec3,
    /// The node's own turning, rad/s, which carries `centre` and `frame` round
    /// within the step — zero for a node that does not go round.
    pub round: Vec3,
}

/// The acceleration of every body: gravity, the springs, and the support.
pub fn accelerations(bodies: &[Body], ground: &Ground) -> Vec<Vec3> {
    accelerations_at(bodies, ground, 0.0)
}

/// The accelerations `elapsed` seconds into a step, with the support and the
/// node's frame turned that far with the ground. Holding the support still
/// across a step while the pieces turn under it is a sideways push of the
/// support times `w dt`: measured, a turning Earth's equatorial faces rang
/// +-490 m about where they are drawn.
fn accelerations_at(bodies: &[Body], ground: &Ground, elapsed: f64) -> Vec<Vec3> {
    let turned = |w: Vec3, v: Vec3| {
        if w == Vec3::ZERO || elapsed == 0.0 { v } else { crate::math::Quat::from_rate(w, elapsed).rotate(v) }
    };
    let swirl = ground.swirl.map(|s| (s.spin, turned(s.round, s.centre), turned(s.round, s.frame)));
    let mut acc: Vec<Vec3> = (0..bodies.len())
        .map(|i| {
            let lift = ground.lift.get(i).copied().unwrap_or(0.0);
            let mut a = ground.field.get(i).copied().unwrap_or(Vec3::ZERO).scale(1.0 - lift);
            if let (Some((w, centre, frame)), true) = (swirl, i >= ground.pieces) {
                let at = centre + bodies[i].pos;
                a += w.cross(w.cross(at)).scale(lift) - frame;
            }
            a
        })
        .collect();
    if let Some(params) = ground.gravity {
        let pieces = &bodies[..ground.pieces.min(bodies.len())];
        if !pieces.is_empty() {
            let mut tree = Octree::build(pieces, params);
            for (a, g) in acc.iter_mut().zip(tree.accelerate_all(pieces)) {
                *a += g;
            }
        }
    }
    for s in &ground.springs {
        let (i, j) = (s.a as usize, s.b as usize);
        let (Some(bi), Some(bj)) = (bodies.get(i), bodies.get(j)) else { continue };
        let d = bj.pos - bi.pos;
        let len = d.norm();
        if !(len > 0.0) {
            continue;
        }
        // Positive when stretched: it pulls the two together.
        let f = d.scale(s.k * (len - s.rest) / len);
        if bi.mass > 0.0 {
            acc[i] += f.scale(1.0 / bi.mass);
        }
        if bj.mass > 0.0 {
            acc[j] -= f.scale(1.0 / bj.mass);
        }
    }
    for (a, p) in acc.iter_mut().zip(&ground.preload) {
        *a += turned(ground.turn, *p);
    }
    acc
}

/// The elastic energy the springs hold, J.
pub fn elastic_energy(bodies: &[Body], springs: &[Spring]) -> f64 {
    springs
        .iter()
        .filter_map(|s| {
            let d = bodies.get(s.b as usize)?.pos - bodies.get(s.a as usize)?.pos;
            let x = d.norm() - s.rest;
            Some(0.5 * s.k * x * x)
        })
        .sum()
}

/// The support that holds each body where it is: minus the acceleration
/// gravity and the springs give it there. Evaluated once, on the ground as
/// drawn, for the pieces; zero for everything else.
pub fn support(bodies: &[Body], ground: &Ground, pieces: &[bool]) -> Vec<Vec3> {
    let free = Ground { preload: vec![Vec3::ZERO; bodies.len()], turn: Vec3::ZERO, lift: Vec::new(), swirl: None, ..ground.clone() };
    accelerations(bodies, &free)
        .into_iter()
        .zip(pieces)
        .map(|(a, &p)| if p { Vec3::ZERO - a } else { Vec3::ZERO })
        .collect()
}

/// One leapfrog step of `dt`.
pub fn step(bodies: &mut [Body], dt: f64, ground: &Ground) -> SolveReport {
    let before = crate::solvers::measure(bodies, 0.0);
    if bodies.is_empty() || !(dt > 0.0) {
        return SolveReport { before, after: before, dt_used: dt, ..Default::default() };
    }
    let half = 0.5 * dt;
    let acc = accelerations_at(bodies, ground, 0.0);
    for (b, a) in bodies.iter_mut().zip(&acc) {
        b.vel += a.scale(half);
        b.pos += b.vel.scale(dt);
    }
    let acc = accelerations_at(bodies, ground, dt);
    let mut unrest: f64 = 0.0;
    for (b, a) in bodies.iter_mut().zip(&acc) {
        b.vel += a.scale(half);
        unrest = unrest.max(a.norm());
    }
    let after = crate::solvers::measure(bodies, 0.0);
    SolveReport { steps: 1, interactions: 0, dt_used: dt, before, after, non_mechanical_energy: 0.0, unrest }
}

/// The longest step the springs are stable at, s: a fraction of the fastest
/// piece's own period, `sqrt(m / k)`.
pub fn stable_step(bodies: &[Body], springs: &[Spring]) -> f64 {
    let mut stiffest = vec![0.0f64; bodies.len()];
    for s in springs {
        for i in [s.a as usize, s.b as usize] {
            if let Some(k) = stiffest.get_mut(i) {
                *k += s.k;
            }
        }
    }
    bodies
        .iter()
        .zip(&stiffest)
        .filter(|(_, &k)| k > 0.0)
        .map(|(b, &k)| 0.2 * (b.mass / k).sqrt())
        .fold(f64::INFINITY, f64::min)
}
