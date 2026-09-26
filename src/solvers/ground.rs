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
//! Within a layer, the diagonals carry shear (`recipe::Touch::Across`). A
//! patch is two layers — its cells, and a column of rock under each — so each
//! cell rests on what is directly beneath it.
//!
//! # Ground as drawn is at rest
//!
//! **The rock carries its weight as a stress between neighbours** — the
//! owner's decision for Phase 5. A planet's rock is under the lithostatic
//! pressure of everything above it, and that stress is present whether or not
//! anything has moved. So each spring is given, once, the tension or
//! compression that holds the pieces where the recipe draws them (`stress`):
//! the one with the least strain energy, which is how the rock would carry it,
//! and along the pieces' current line, so it is central and conserves angular
//! momentum exactly. What is outside the ground holds the rest: a patch's floor
//! bears its stack and rests on the rest of its planet (`Ground::foundation`),
//! and whatever net force and torque are left keep the patch with its node
//! (`Ground::held`). The configuration the recipe draws is then an equilibrium:
//! ground nobody disturbs stays where it is drawn and regenerates bit for bit,
//! and a disturbance rings through the springs at the rock's own sound speed.
//!
//! Measured on a turning Earth with one face promoted: its angular momentum
//! holds to 1.8e-15 over six hours; the face, solved every frame under a mass
//! of air for a day, has its pieces slip over the turning ground at 2.4e-3 m/s
//! at worst. Each step there was measured: a fixed support per piece moved the
//! world's angular momentum by 4e-7; a stress along the line as drawn by as
//! much; along the current line but spread by least squares, the face buckled
//! like a shell; with one substrate body under a face its corners buckled; and
//! pulling its pieces through a Barnes-Hut tree, whose error changed as they
//! turned through its cells, drove its ringing to 27 m/s in two days.

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
    /// The line from `a` to `b` as drawn, unit, in the node's own axes: what
    /// `stress` solved the tensions along. The stress then acts along the
    /// pieces' current line — central, so it conserves angular momentum
    /// exactly — which is the owner's decision for Phase 5.
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
    /// Where each floor piece rests on what is under the ground, in the
    /// node's own axes, and how stiffly (`Ground::foundation`).
    pub foundation: Vec<(usize, Vec3, f64)>,
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
    /// The first of the pieces that rest on what is outside the ground — its
    /// floor (`recipe::Tiled::floor`), which is what `foundation` anchors.
    pub floor: usize,
    /// **What a patch's floor rests on**: the rest of its planet, as an
    /// elastic anchor at each floor piece's place as drawn — in the bodies'
    /// axes at the step's start, turning with the ground — of `E A / L` over
    /// the depth beneath it: `(piece, place, stiffness in N/m)`. Resting on a
    /// fixed support alone, nothing pushed back on a face's corner column
    /// pushed down but its neighbours, and its corners grew at 7.1e-7 s^-2.
    /// What it does to the patch moving as a whole is taken off with the rest
    /// of what holds it (`held`).
    pub foundation: Vec<(usize, Vec3, f64)>,
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
        // Positive when stretched: it pulls the two together, and the stress
        // it was drawn carrying acts along the same line.
        let f = d.scale((s.k * (len - s.rest) + s.tension) / len);
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
    for &(i, at, k) in &ground.foundation {
        if let Some(b) = bodies.get(i) {
            if b.mass > 0.0 {
                acc[i] += (turned(ground.turn, at) - b.pos).scale(k / b.mass);
            }
        }
    }
    if ground.held {
        // What holds the ground is exactly what keeps its centre of mass with
        // the node and its turning with the ground: whatever net force and
        // torque the pieces' own forces leave over that. Measured without it, a face of a turning Earth sloshed inside
        // itself at 2.5 m/s and then twisted by 3e29 kg m^2/s in forty minutes.
        let n = ground.pieces.min(bodies.len());
        let mass: f64 = bodies[..n].iter().map(|b| b.mass).sum();
        if mass > 0.0 {
            let centre = (0..n).fold(Vec3::ZERO, |c, i| c + bodies[i].pos.scale(bodies[i].mass)).scale(1.0 / mass);
            let w = ground.turn;
            let force = (0..n).fold(Vec3::ZERO, |f, i| f + acc[i].scale(bodies[i].mass));
            let torque = (0..n).fold(Vec3::ZERO, |t, i| {
                let d = bodies[i].pos - centre;
                t + d.cross(acc[i] - w.cross(w.cross(d))).scale(bodies[i].mass)
            });
            // Taken off every piece alike, by mass — watching the ground from
            // its own centre of mass, turning with it. Taken off the floor
            // alone it was an oblique projection, and it pumped a face under
            // air from 30 to 309 m/s over twelve hours; where the ground's
            // weight goes is its stress's business (`stress`), not this.
            let all: Vec<usize> = (0..n).collect();
            for (i, f) in rigid_share(bodies, &all, centre, force, torque) {
                acc[i] -= f.scale(1.0 / bodies[i].mass);
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
    let free = Ground { preload: vec![Vec3::ZERO; bodies.len()], turn: Vec3::ZERO, lift: Vec::new(), swirl: None, held: false, floor: 0, foundation: Vec::new(), ..ground.clone() };
    accelerations(bodies, &free)
        .into_iter()
        .zip(pieces)
        .map(|(a, &p)| if p { Vec3::ZERO - a } else { Vec3::ZERO })
        .collect()
}

/// Solve `A x = b` for a symmetric positive semi-definite `A` (row-major,
/// `dim` square, overwritten) by Cholesky, with the null space a ground's
/// free rigid motions leave regularised away at a part in 10^12 of its
/// largest diagonal. `b` is in the range whenever it is a load the pieces'
/// own springs can carry, which is what `stress` hands it.
fn solve_spd(a: &mut [f64], dim: usize, b: &[f64]) -> Vec<f64> {
    let top = (0..dim).map(|i| a[i * dim + i]).fold(0.0f64, f64::max);
    let eps = 1e-12 * top.max(1e-300);
    for i in 0..dim {
        a[i * dim + i] += eps;
    }
    for j in 0..dim {
        let mut d = a[j * dim + j];
        for k in 0..j {
            d -= a[j * dim + k] * a[j * dim + k];
        }
        let d = d.max(eps).sqrt();
        a[j * dim + j] = d;
        for i in (j + 1)..dim {
            let mut v = a[i * dim + j];
            for k in 0..j {
                v -= a[i * dim + k] * a[j * dim + k];
            }
            a[i * dim + j] = v / d;
        }
    }
    let mut z = b.to_vec();
    for i in 0..dim {
        for k in 0..i {
            z[i] -= a[i * dim + k] * z[k];
        }
        z[i] /= a[i * dim + i];
    }
    for i in (0..dim).rev() {
        for k in (i + 1)..dim {
            z[i] -= a[k * dim + i] * z[k];
        }
        z[i] /= a[i * dim + i];
    }
    z
}

/// A net force and a torque about `centre`, shared among the pieces `on` as
/// one rigid acceleration of them — `m (A + B x (r - c))` about their own
/// centre of mass `c` — so that the shares add up to exactly that force and
/// that torque. Returns each piece's share, N.
fn rigid_share(bodies: &[Body], on: &[usize], centre: Vec3, force: Vec3, torque: Vec3) -> Vec<(usize, Vec3)> {
    let mass: f64 = on.iter().map(|&i| bodies[i].mass).sum();
    if !(mass > 0.0) {
        return Vec::new();
    }
    let own = on.iter().fold(Vec3::ZERO, |c, &i| c + bodies[i].pos.scale(bodies[i].mass)).scale(1.0 / mass);
    let linear = force.scale(1.0 / mass);
    let mut inertia = crate::math::Mat3::zero();
    for &i in on {
        let (m, d) = (bodies[i].mass, bodies[i].pos - own);
        let o = d.outer(d);
        for r in 0..3 {
            for c in 0..3 {
                let delta = if r == c { d.norm2() } else { 0.0 };
                inertia.0[r][c] += m * (delta - o.0[r][c]);
            }
        }
    }
    // The force applied at their own centre is a torque about `centre` too.
    let angular = inertia.solve(torque - (own - centre).cross(force)).unwrap_or(Vec3::ZERO);
    on.iter().map(|&i| (i, (linear + angular.cross(bodies[i].pos - own)).scale(bodies[i].mass))).collect()
}

/// The rock's stress as drawn: how each piece's weight is carried.
///
/// `need` is the force each piece needs to stay where it is drawn, N —
/// against gravity and the springs, and toward the axis it goes round. Part of
/// it is what holds the whole. A patch of a planet is held up by the rest of
/// the planet under its floor — the pieces from `floor` on — and each of them
/// bears what rests on it, `bears[i]` naming the floor piece piece `i`'s load
/// goes to. A body resting on nothing is held by no more than one rigid
/// acceleration of it with the same net force and torque as `need`. That part comes back as a force on
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
pub fn stress(bodies: &[Body], springs: &[Spring], pieces: usize, floor: usize, bears: &[usize], need: &[Vec3]) -> (Vec<f64>, Vec<Vec3>, f64) {
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
    if floor > 0 && bears.len() >= n {
        // **Each piece of the floor bears what rests on it**: the pressure
        // under a column is along its own vertical and is the weight of its
        // stack, which is what holds a curved patch up. Held instead as one
        // rigid acceleration of the whole floor, a face's edge columns — whose
        // gravity is 45 degrees off the face's middle — pulled sideways on the
        // lattice at fifteen times what its rock can carry.
        for i in 0..n {
            outside[bears[i].min(n - 1)] += need[i];
        }
    } else {
        let on: Vec<usize> = (floor.min(n)..n).collect();
        for (i, f) in rigid_share(bodies, &on, centre, net, torque) {
            outside[i] = f;
        }
    }
    // The rest, carried along the springs — **as the rock carries it**: of
    // every set of tensions that holds the pieces, the one with the least
    // strain energy, `sum T^2 / 2k`, which is Menabrea's theorem for a
    // structure with more members than it needs. Written as `T = sqrt(k) s`,
    // that is the least-norm `s` with `A sqrt(k) s = g`, which conjugate
    // gradients on the normal equations find from zero. The least-norm `T`
    // instead spread each cell's weight across its neighbours as much as
    // down to what is under it, put a face's surface into compression in its
    // own plane, and it buckled like a shell: two modes of the face grew at
    // 6.4e-6 s^-2, a tenfold every fifteen minutes.
    let g: Vec<Vec3> = (0..n).map(|i| need[i] - outside[i]).collect();
    let usable: Vec<Option<(usize, usize, Vec3)>> = springs
        .iter()
        .map(|s| {
            let (a, b) = (s.a as usize, s.b as usize);
            if a >= n || b >= n || !(s.k > 0.0) {
                return None;
            }
            let d = bodies[b].pos - bodies[a].pos;
            let len = d.norm();
            (len > 0.0).then(|| (a, b, d.scale(s.k.sqrt() / len)))
        })
        .collect();
    // Exactly, by the dual: `B B^T y = g` over the pieces' coordinates, with
    // `B` putting `sqrt(k) s` on a spring's two ends, and then `s = B^T y`.
    // Solving for the tensions directly by conjugate gradients instead let in
    // self-stress — tensions that cancel at every piece and so move nothing
    // but soften every compressed line sideways — at 2000 times what a spring
    // carries at its own length, and a face buckled at 8e-2 s^-2.
    let dim = 3 * n;
    let mut normal = vec![0.0f64; dim * dim];
    for u in usable.iter().flatten() {
        let (a, b, w) = *u;
        let w = [w.x, w.y, w.z];
        for r in 0..3 {
            for c in 0..3 {
                let x = w[r] * w[c];
                normal[(3 * a + r) * dim + 3 * a + c] += x;
                normal[(3 * b + r) * dim + 3 * b + c] += x;
                normal[(3 * a + r) * dim + 3 * b + c] -= x;
                normal[(3 * b + r) * dim + 3 * a + c] -= x;
            }
        }
    }
    let rhs: Vec<f64> = g.iter().flat_map(|v| [v.x, v.y, v.z]).collect();
    let y = solve_spd(&mut normal, dim, &rhs);
    let mut t: Vec<f64> = usable
        .iter()
        .map(|u| match u {
            Some((a, b, w)) => w.dot(crate::math::v3(y[3 * a] - y[3 * b], y[3 * a + 1] - y[3 * b + 1], y[3 * a + 2] - y[3 * b + 2])),
            None => 0.0,
        })
        .collect();
    for (e, s) in springs.iter().enumerate() {
        t[e] *= s.k.max(0.0).sqrt();
    }
    let carried = {
        let mut f = vec![Vec3::ZERO; n];
        for (e, u) in usable.iter().enumerate() {
            if let Some((a, b, u)) = u {
                let dir = u.scale(1.0 / springs[e].k.sqrt());
                f[*a] += dir.scale(t[e]);
                f[*b] -= dir.scale(t[e]);
            }
        }
        f
    };
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

/// The longest step the ground is stable at, s: a fraction of the fastest
/// piece's own period, `sqrt(m / k)` over every stiffness it has — its springs
/// and what its floor rests on. Counting the springs alone, a face's corner
/// columns, shallow and so stiffly anchored, were stepped past their own
/// period and grew to 500 m/s in twelve hours; at a quarter of the step the
/// same face slipped at 5e-4 m/s.
pub fn stable_step(bodies: &[Body], springs: &[Spring], foundation: &[(usize, Vec3, f64)]) -> f64 {
    let mut stiffest = vec![0.0f64; bodies.len()];
    for s in springs {
        for i in [s.a as usize, s.b as usize] {
            if let Some(k) = stiffest.get_mut(i) {
                *k += s.k;
            }
        }
    }
    for &(i, _, k) in foundation {
        if let Some(s) = stiffest.get_mut(i) {
            *s += k;
        }
    }
    bodies
        .iter()
        .zip(&stiffest)
        .filter(|(_, &k)| k > 0.0)
        .map(|(b, &k)| 0.2 * (b.mass / k).sqrt())
        .fold(f64::INFINITY, f64::min)
}
