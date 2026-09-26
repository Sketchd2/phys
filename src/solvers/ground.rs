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
//! does over the time anything here watches it. Measured, a face of a turning
//! Earth solved under a mass of air for a day rings at 8.1 m/s at its worst
//! piece, and does not grow over two — once the pieces pull each other pair by
//! pair; through a Barnes-Hut tree the error changed as they turned through its
//! cells, and drove the ringing to 27 m/s in two days.

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
    /// Length at which it carries no more than `tension`, m: the pieces'
    /// distance as drawn.
    pub rest: f64,
    /// What it carries, N — positive pulling the two together, negative
    /// pushing them apart: the rock's stress as drawn (`stress`).
    pub tension: f64,
    /// The line from `a` to `b` as drawn, unit, which that stress acts along.
    /// Kept in the node's own axes in a `Cache` and in the bodies' within a
    /// step, where it turns with the ground (`Ground::turn`) — not with the
    /// two pieces. Along their current line instead, a compression is also a
    /// sideways push the moment the line bends, and a face of a turning Earth
    /// carrying its weight that way buckled: its pieces slipped over the
    /// ground at 124 m/s by the second hour and climbing.
    pub along: Vec3,
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
    /// The share of what the pieces need that the springs could not carry
    /// and is left on each piece instead, rms against the need (`stress`).
    pub unresolved: f64,
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
    /// Whether something outside holds these pieces — a patch of a planet,
    /// held by the rest of it — so that what holds them is exactly what keeps
    /// their centre of mass with the node: the pieces' mean acceleration is
    /// taken off each of them. Measured without it, a face of a turning Earth
    /// held by a force derived where it was drawn sloshed inside itself at 2.5
    /// m/s against the face it is, as the face's real motion drifted from the
    /// drawing's. `false` for a body nothing holds, whose pieces' forces
    /// already sum to nothing.
    pub held: bool,
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
        // Positive when stretched: it pulls the two together. The stress it
        // was drawn carrying acts along the line as drawn, turned with the
        // ground.
        let f = d.scale(s.k * (len - s.rest) / len) + turned(ground.turn, s.along).scale(s.tension);
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
    if ground.held {
        let n = ground.pieces.min(bodies.len());
        let mass: f64 = bodies[..n].iter().map(|b| b.mass).sum();
        if mass > 0.0 {
            let mean = (0..n).fold(Vec3::ZERO, |m, i| m + acc[i].scale(bodies[i].mass)).scale(1.0 / mass);
            for a in acc[..n].iter_mut() {
                *a -= mean;
            }
            // And turning with it: the net torque about the pieces' centre is
            // what going round at the ground's rate needs, and no more.
            // Measured without it, the same face twisted inside itself by
            // 3e29 kg m^2/s against the ground's turning in forty minutes.
            let centre = (0..n).fold(Vec3::ZERO, |c, i| c + bodies[i].pos.scale(bodies[i].mass)).scale(1.0 / mass);
            let w = ground.turn;
            let mut excess = Vec3::ZERO;
            let mut inertia = crate::math::Mat3::zero();
            for i in 0..n {
                let (m, d) = (bodies[i].mass, bodies[i].pos - centre);
                excess += d.cross(acc[i] - w.cross(w.cross(d))).scale(m);
                let o = d.outer(d);
                for r in 0..3 {
                    for c in 0..3 {
                        let delta = if r == c { d.norm2() } else { 0.0 };
                        inertia.0[r][c] += m * (delta - o.0[r][c]);
                    }
                }
            }
            if let Some(beta) = inertia.solve(excess) {
                for i in 0..n {
                    acc[i] -= beta.cross(bodies[i].pos - centre);
                }
            }
        }
    }
    acc
}

/// The elastic energy the springs hold beyond their stress as drawn, J.
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
    let free = Ground { preload: vec![Vec3::ZERO; bodies.len()], turn: Vec3::ZERO, lift: Vec::new(), swirl: None, held: false, ..ground.clone() };
    accelerations(bodies, &free)
        .into_iter()
        .zip(pieces)
        .map(|(a, &p)| if p { Vec3::ZERO - a } else { Vec3::ZERO })
        .collect()
}

/// The rock's stress as drawn: how each piece's weight is carried.
///
/// `need` is the force each piece needs to stay where it is drawn, N —
/// against gravity and the springs, and toward the axis it goes round. Part of
/// it is what holds the whole: a patch of a planet is held up by the rest of
/// the planet, and that is one rigid acceleration of all its pieces with the
/// same net force and torque as `need`. That part comes back as a force on
/// each piece (`Vec3`s, N). **The rest sums to no force and no torque, and the
/// springs carry it** — a tension or a compression along each, found by least
/// squares — so it is equal and opposite between two pieces along the line
/// joining them, whatever they then do. What the springs cannot carry is added
/// to the per-piece forces, and its size returned, so that it is measured
/// rather than assumed away.
///
/// The owner's decision for Phase 5. Given as a fixed force per piece, the
/// support stopped summing to no torque the moment the pieces rang off where
/// they were drawn: a turning Earth with one face promoted moved its angular
/// momentum by 4e-7 of itself in six hours, in steps at each ground solve.
pub fn stress(bodies: &[Body], springs: &[Spring], pieces: usize, need: &[Vec3]) -> (Vec<f64>, Vec<Vec3>, f64) {
    let n = pieces.min(bodies.len()).min(need.len());
    let mut outside = vec![Vec3::ZERO; need.len()];
    if n == 0 {
        return (vec![0.0; springs.len()], outside, 0.0);
    }
    // The rigid part: net force and torque about the pieces' centre of mass.
    let mass: f64 = bodies[..n].iter().map(|b| b.mass).sum();
    if !(mass > 0.0) {
        return (vec![0.0; springs.len()], need.to_vec(), 0.0);
    }
    let centre = bodies[..n].iter().fold(Vec3::ZERO, |a, b| a + b.pos.scale(b.mass)).scale(1.0 / mass);
    let net: Vec3 = need[..n].iter().fold(Vec3::ZERO, |a, f| a + *f);
    let torque: Vec3 = (0..n).fold(Vec3::ZERO, |a, i| a + (bodies[i].pos - centre).cross(need[i]));
    let mut inertia = crate::math::Mat3::zero();
    for b in &bodies[..n] {
        let d = b.pos - centre;
        let o = d.outer(d);
        for r in 0..3 {
            for c in 0..3 {
                let delta = if r == c { d.norm2() } else { 0.0 };
                inertia.0[r][c] += b.mass * (delta - o.0[r][c]);
            }
        }
    }
    let linear = net.scale(1.0 / mass);
    let angular = inertia.solve(torque).unwrap_or(Vec3::ZERO);
    for i in 0..n {
        outside[i] = (linear + angular.cross(bodies[i].pos - centre)).scale(bodies[i].mass);
    }
    // The rest, carried along the springs: minimise |A t - g|^2 by conjugate
    // gradients on the normal equations, where `A` puts a spring's tension on
    // its two ends along the line between them.
    let g: Vec<Vec3> = (0..n).map(|i| need[i] - outside[i]).collect();
    let usable: Vec<Option<(usize, usize, Vec3)>> = springs
        .iter()
        .map(|s| {
            let (a, b) = (s.a as usize, s.b as usize);
            if a >= n || b >= n {
                return None;
            }
            let d = bodies[b].pos - bodies[a].pos;
            let len = d.norm();
            (len > 0.0).then(|| (a, b, d.scale(1.0 / len)))
        })
        .collect();
    let apply = |t: &[f64]| -> Vec<Vec3> {
        let mut f = vec![Vec3::ZERO; n];
        for (e, u) in usable.iter().enumerate() {
            if let Some((a, b, u)) = u {
                f[*a] += u.scale(t[e]);
                f[*b] -= u.scale(t[e]);
            }
        }
        f
    };
    let apply_t = |v: &[Vec3]| -> Vec<f64> {
        usable.iter().map(|u| u.map(|(a, b, u)| (v[a] - v[b]).dot(u)).unwrap_or(0.0)).collect()
    };
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(a, b)| a * b).sum::<f64>();
    let mut t = vec![0.0; springs.len()];
    let mut r = apply_t(&g);
    let mut p = r.clone();
    let mut rr = dot(&r, &r);
    let scale = rr;
    for _ in 0..(4 * springs.len().max(1)) {
        if !(rr > 1e-28 * scale) {
            break;
        }
        let ap = apply_t(&apply(&p));
        let pap = dot(&p, &ap);
        if !(pap > 0.0) {
            break;
        }
        let alpha = rr / pap;
        for e in 0..t.len() {
            t[e] += alpha * p[e];
            r[e] -= alpha * ap[e];
        }
        let next = dot(&r, &r);
        let beta = next / rr;
        rr = next;
        for e in 0..p.len() {
            p[e] = r[e] + beta * p[e];
        }
    }
    let carried = apply(&t);
    let mut left: f64 = 0.0;
    let mut total: f64 = 0.0;
    for i in 0..n {
        let rest = g[i] - carried[i];
        outside[i] += rest;
        left += rest.norm2();
        total += need[i].norm2();
    }
    let unresolved = if total > 0.0 { (left / total).sqrt() } else { 0.0 };
    (t, outside, unresolved)
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
