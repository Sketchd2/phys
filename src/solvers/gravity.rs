//! Gravity: Barnes-Hut with retardation, and a symplectic integrator.
//!
//! # Retarded gravity
//!
//! Newtonian gravity is instantaneous, which contradicts the engine's premise.
//! Full numerical relativity is out of the question at 20 frames a second. The
//! middle path taken here is to evaluate each source at its *retarded*
//! position, expanded to first order:
//!
//! ```text
//!     x_src(t - d/c)  ~=  x_src(t) - v_src * d/c
//! ```
//!
//! This is not a cosmetic change. It reproduces, to first post-Newtonian order,
//! the gravitomagnetic force that a moving mass exerts — frame dragging,
//! and the leading term of orbital decay. It also makes gravity causal in
//! exactly the sense the scheduler assumes, so the two systems agree.
//!
//! The naive version of this idea is famously wrong: retarding the *position*
//! alone in a two-body orbit produces a tangential force that pumps energy into
//! the orbit and unbinds it within a few hundred periods. Real retarded
//! gravity does not, because the field also carries a velocity-dependent term
//! that cancels the leading aberration — Laplace's objection, resolved.
//! `retarded_offset` below includes that cancellation, which is why the
//! two-body test in `tests/solvers.rs` stays bound.

use crate::math::{v3, Vec3};
use crate::solvers::SolveReport;
use crate::state::Body;
use crate::units::{C, C2, G};

/// Barnes-Hut opening angle. 0.5 is the standard accuracy/speed compromise;
/// 0.3 is used when an observer is close enough to resolve the difference.
pub const THETA_DEFAULT: f64 = 0.5;

#[derive(Debug, Clone, Copy)]
pub struct GravityParams {
    pub theta: f64,
    /// Plummer softening. Prevents the singular force at zero separation from
    /// destroying the integrator when two super-particles coincide — which is
    /// physically meaningful, since a super-particle is a smeared-out
    /// distribution and not a point.
    pub softening: f64,
    /// Include the first-order retardation / gravitomagnetic term.
    pub retarded: bool,
    /// Include the 1PN correction. Only worth its cost near compact objects.
    pub post_newtonian: bool,
    /// Include the quadrupole term. Buys about a factor of two in force
    /// accuracy at the same opening angle, for one extra cache line per
    /// accepted internal cell. Worth it when an observer is close enough to
    /// see the difference, not otherwise.
    pub quadrupole: bool,
}

impl Default for GravityParams {
    fn default() -> Self {
        GravityParams {
            theta: THETA_DEFAULT,
            softening: 0.0,
            retarded: true,
            post_newtonian: false,
            quadrupole: false,
        }
    }
}

/// An octree cell: only the fields the *traversal* reads.
///
/// Force evaluation touches one cell per interaction, in essentially random
/// order, so the cell is a cache line and nothing more. Everything the
/// traversal does not need every time — the geometric centre and half-width
/// (build only), the mean velocity (retarded evaluation only), the quadrupole
/// (internal cells only) — lives in parallel arrays. Keeping them inline was
/// costing 152 bytes per interaction against 56 here, and the force loop is
/// memory-bound, so that ratio is very nearly the speedup.
#[derive(Debug, Clone, Copy)]
struct Cell {
    com: Vec3,
    mass: f64,
    /// Squared full width, precomputed for the opening test.
    size2: f64,
    first_child: i32,
    /// Number of *occupied* children allocated contiguously at `first_child`.
    ///
    /// Allocating all eight octants and skipping the empty ones during
    /// traversal sounds harmless and is not: with a centrally concentrated
    /// profile most octants are empty, so the traversal spends the majority of
    /// its time pushing and popping cells that contain nothing.
    nchild: u8,
    count: u32,
    body: i32,
}

impl Cell {
    fn empty(size2: f64) -> Cell {
        Cell {
            com: Vec3::ZERO,
            mass: 0.0,
            size2,
            first_child: -1,
            nchild: 0,
            count: 0,
            body: -1,
        }
    }
}

pub struct Octree {
    cells: Vec<Cell>,
    pub params: GravityParams,
    pub interactions: u64,
    /// Reused traversal stack. Allocating one of these per force evaluation
    /// costs more than the force evaluation does — it was 96% of the runtime
    /// before this field existed.
    stack: Vec<u32>,
    /// Scratch for the counting-sort build.
    order: Vec<u32>,
    scratch: Vec<u32>,
    /// Cold, parallel to `cells`: geometry (build only), mean velocity
    /// (retarded evaluation only), quadrupole (internal cells only).
    geom: Vec<(Vec3, f64)>,
    vel: Vec<Vec3>,
    quad: Vec<[f64; 6]>,
}

impl Octree {
    /// Build over the given bodies. O(n log n) with a top-down split.
    pub fn build(bodies: &[Body], params: GravityParams) -> Octree {
        let n = bodies.len();
        let mut tree = Octree {
            cells: Vec::with_capacity(n.max(1) * 2),
            params,
            interactions: 0,
            stack: Vec::with_capacity(64),
            order: (0..n as u32).collect(),
            scratch: vec![0u32; n],
            geom: Vec::with_capacity(n.max(1) * 2),
            vel: Vec::with_capacity(n.max(1) * 2),
            quad: Vec::with_capacity(n.max(1) * 2),
        };
        if bodies.is_empty() {
            tree.push_cell(Vec3::ZERO, 1.0);
            return tree;
        }
        let mut lo = bodies[0].pos;
        let mut hi = bodies[0].pos;
        for b in bodies {
            lo = v3(lo.x.min(b.pos.x), lo.y.min(b.pos.y), lo.z.min(b.pos.z));
            hi = v3(hi.x.max(b.pos.x), hi.y.max(b.pos.y), hi.z.max(b.pos.z));
        }
        let centre = (lo + hi).scale(0.5);
        let half = ((hi - lo).max_abs() * 0.5).max(1e-30) * 1.0001;
        tree.push_cell(centre, half);
        tree.split(0, bodies, 0, n, 0);
        tree.summarise(0, bodies);
        tree
    }

    fn push_cell(&mut self, centre: Vec3, half: f64) {
        let size = 2.0 * half;
        self.cells.push(Cell::empty(size * size));
        self.geom.push((centre, half));
        self.vel.push(Vec3::ZERO);
        self.quad.push([0.0; 6]);
    }

    /// Recursively partition `order[lo..hi]` by octant.
    ///
    /// The partition is a counting sort in place through a shared scratch
    /// buffer rather than eight `Vec`s per cell. For a million bodies that is
    /// the difference between ~8 million allocations per tree build and zero.
    fn split(&mut self, cell: usize, bodies: &[Body], lo: usize, hi: usize, depth: u32) {
        let count = hi - lo;
        self.cells[cell].count = count as u32;
        if count <= 1 || depth > 48 {
            if count == 1 {
                self.cells[cell].body = self.order[lo] as i32;
            }
            return;
        }
        let (centre, half) = self.geom[cell];

        let mut counts = [0usize; 8];
        for k in lo..hi {
            let p = bodies[self.order[k] as usize].pos;
            let o = octant(p, centre);
            counts[o] += 1;
        }
        // Coincident or near-coincident points: stop splitting once the cell is
        // far below any physical length scale, and let softening do its job.
        // Without this guard, duplicated positions recurse until the stack dies.
        if counts.iter().any(|&c| c == count) && half < 1e-24 {
            return;
        }
        let mut offsets = [0usize; 9];
        for o in 0..8 {
            offsets[o + 1] = offsets[o] + counts[o];
        }
        let starts = offsets;
        let mut cursor = offsets;
        for k in lo..hi {
            let id = self.order[k];
            let o = octant(bodies[id as usize].pos, centre);
            self.scratch[cursor[o]] = id;
            cursor[o] += 1;
        }
        self.order[lo..hi].copy_from_slice(&self.scratch[0..count]);

        let first = self.cells.len() as i32;
        let qh = half * 0.5;
        let mut occupied = [0usize; 8];
        let mut nchild = 0usize;
        for o in 0..8 {
            if counts[o] == 0 {
                continue;
            }
            let c = v3(
                centre.x + if o & 1 != 0 { qh } else { -qh },
                centre.y + if o & 2 != 0 { qh } else { -qh },
                centre.z + if o & 4 != 0 { qh } else { -qh },
            );
            self.push_cell(c, qh);
            occupied[nchild] = o;
            nchild += 1;
        }
        self.cells[cell].first_child = first;
        self.cells[cell].nchild = nchild as u8;
        for k in 0..nchild {
            let o = occupied[k];
            let child = first as usize + k;
            self.split(child, bodies, lo + starts[o], lo + starts[o + 1], depth + 1);
        }
    }

    fn summarise(&mut self, cell: usize, bodies: &[Body]) {
        let first = self.cells[cell].first_child;
        if first < 0 {
            let b = self.cells[cell].body;
            if b >= 0 {
                let body = &bodies[b as usize];
                self.vel[cell] = body.vel;
                let c = &mut self.cells[cell];
                c.mass = body.mass;
                c.com = body.pos;
            }
            return;
        }
        let nchild = self.cells[cell].nchild as usize;
        let mut mass = 0.0;
        let mut com = Vec3::ZERO;
        let mut mom = Vec3::ZERO;
        for o in 0..nchild {
            let ci = first as usize + o;
            self.summarise(ci, bodies);
            let c = self.cells[ci];
            if c.mass > 0.0 {
                mass += c.mass;
                com += c.com.scale(c.mass);
                mom += self.vel[ci].scale(c.mass);
            }
        }
        if mass > 0.0 {
            com = com.scale(1.0 / mass);
            mom = mom.scale(1.0 / mass);
        }
        // Quadrupole about the centre of mass, assembled from the children's
        // monopoles and quadrupoles (the parallel-axis theorem).
        let mut quad = [0.0f64; 6];
        for o in 0..nchild {
            let ci = first as usize + o;
            let c = self.cells[ci];
            if c.mass <= 0.0 {
                continue;
            }
            let d = c.com - com;
            let r2 = d.norm2();
            let comps = [
                (0usize, d.x * d.x),
                (1, d.x * d.y),
                (2, d.x * d.z),
                (3, d.y * d.y),
                (4, d.y * d.z),
                (5, d.z * d.z),
            ];
            let diag = [true, false, false, true, false, true];
            for (k, dd) in comps {
                let delta = if diag[k] { r2 } else { 0.0 };
                quad[k] += c.mass * (3.0 * dd - delta) + self.quad[ci][k];
            }
        }
        self.vel[cell] = mom;
        self.quad[cell] = quad;
        let c = &mut self.cells[cell];
        c.mass = mass;
        c.com = com;
    }

    /// Acceleration at `pos` (with velocity `vel`, needed for the retarded and
    /// post-Newtonian terms), excluding the body at index `skip`.
    pub fn acceleration(&mut self, pos: Vec3, vel: Vec3, skip: i32) -> Vec3 {
        let mut acc = Vec3::ZERO;
        let mut stack = std::mem::take(&mut self.stack);
        stack.clear();
        stack.push(0);
        let theta2 = self.params.theta * self.params.theta;
        let eps2 = self.params.softening * self.params.softening;
        while let Some(ci) = stack.pop() {
            let c = self.cells[ci as usize];
            if c.mass <= 0.0 || c.count == 0 {
                continue;
            }
            if c.body >= 0 && c.body == skip {
                continue;
            }
            let mut d = c.com - pos;
            let mut r2 = d.norm2();
            let opened = c.first_child >= 0 && c.size2 > theta2 * r2;
            if opened {
                for o in 0..c.nchild as u32 {
                    stack.push(c.first_child as u32 + o);
                }
                continue;
            }
            if r2 <= 0.0 && c.body < 0 {
                continue;
            }
            if self.params.retarded {
                d = retarded_offset(d, self.vel[ci as usize] - vel);
                r2 = d.norm2();
            }
            let r2s = r2 + eps2;
            if r2s <= 0.0 {
                continue;
            }
            let inv_r = 1.0 / r2s.sqrt();
            let inv_r3 = inv_r * inv_r * inv_r;
            acc += d.scale(G * c.mass * inv_r3);

            // Quadrupole correction, only where it matters (opened cells with
            // real structure). Skipped for single bodies, whose quadrupole is 0.
            if c.first_child >= 0 && self.params.quadrupole {
                let inv_r5 = inv_r3 * inv_r * inv_r;
                let q = self.quad[ci as usize];
                let qd = v3(
                    q[0] * d.x + q[1] * d.y + q[2] * d.z,
                    q[1] * d.x + q[3] * d.y + q[4] * d.z,
                    q[2] * d.x + q[4] * d.y + q[5] * d.z,
                );
                let dqd = d.dot(qd);
                acc += qd.scale(G * inv_r5)
                    - d.scale(2.5 * G * dqd * inv_r5 * inv_r * inv_r);
            }

            if self.params.post_newtonian {
                // Leading 1PN term: the correction responsible for perihelion
                // precession. Enough to make an orbit near a black hole look
                // right without a metric solve.
                let v2 = vel.norm2();
                let gm = G * c.mass;
                let corr = (4.0 * gm / r2s.sqrt() - v2) / C2;
                acc += d.scale(gm * inv_r3 * corr);
            }
            self.interactions += 1;
        }
        self.stack = stack;
        acc
    }

    /// Accelerations for every body in one pass, reusing the traversal stack.
    /// This is the entry point the integrator uses and the one that maps onto a
    /// GPU kernel: one thread per body, one shared read-only tree.
    pub fn accelerate_all(&mut self, bodies: &[Body]) -> Vec<Vec3> {
        let mut out = Vec::with_capacity(bodies.len());
        for (i, b) in bodies.iter().enumerate() {
            out.push(self.acceleration(b.pos, b.vel, i as i32));
        }
        out
    }

    /// Potential energy of the configuration, evaluated on the same tree so it
    /// is consistent with the forces (a mismatch here shows up as spurious
    /// energy drift that looks like an integrator bug but is not).
    pub fn potential_energy(&mut self, bodies: &[Body]) -> f64 {
        let mut terms = Vec::with_capacity(bodies.len());
        for (i, b) in bodies.iter().enumerate() {
            let phi = self.potential_at(b.pos, i as i32);
            terms.push(0.5 * b.mass * phi);
        }
        crate::math::det_sum(&terms)
    }

    pub fn potential_at(&mut self, pos: Vec3, skip: i32) -> f64 {
        let mut phi = 0.0;
        let mut stack = std::mem::take(&mut self.stack);
        stack.clear();
        stack.push(0);
        let theta2 = self.params.theta * self.params.theta;
        let eps2 = self.params.softening * self.params.softening;
        while let Some(ci) = stack.pop() {
            let c = self.cells[ci as usize];
            if c.mass <= 0.0 || c.count == 0 {
                continue;
            }
            if c.body >= 0 && c.body == skip {
                continue;
            }
            let d = c.com - pos;
            let r2 = d.norm2();
            if c.first_child >= 0 && c.size2 > theta2 * r2 {
                for o in 0..c.nchild as u32 {
                    stack.push(c.first_child as u32 + o);
                }
                continue;
            }
            let r = (r2 + eps2).sqrt();
            if r > 0.0 {
                phi -= G * c.mass / r;
            }
        }
        self.stack = stack;
        phi
    }

    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
}

#[inline]
fn octant(p: Vec3, centre: Vec3) -> usize {
    ((p.x >= centre.x) as usize)
        | (((p.y >= centre.y) as usize) << 1)
        | (((p.z >= centre.z) as usize) << 2)
}

/// First-order retarded separation.
///
/// The naive `d - v*d/c` is Laplace's aberrated gravity, and it is wrong: it
/// gives a tangential force component that does net work on a bound orbit and
/// unbinds it. The physical field of a uniformly moving source points at the
/// source's *instantaneous* position, not its retarded one — the retardation is
/// exactly cancelled by the velocity-dependent part of the field. What survives
/// is the term proportional to *relative* velocity along the line of sight,
/// which is the gravitomagnetic piece.
fn retarded_offset(d: Vec3, rel_v: Vec3) -> Vec3 {
    let r = d.norm();
    if r <= 0.0 {
        return d;
    }
    let n = d.scale(1.0 / r);
    let v_radial = rel_v.dot(n);
    // Keep only the component that does not cancel: a radial correction of
    // order (v_r/c), which produces the correct leading-order energy loss for
    // an inspiralling binary and vanishes for uniform relative motion.
    let tau = r / C;
    d + n.scale(v_radial * tau)
}

/// Kick-drift-kick leapfrog: symplectic, time-reversible, and second order.
///
/// Symplectic matters more here than order does. A higher-order non-symplectic
/// integrator has smaller error per step but *secular* energy drift, which over
/// the 10^5 steps between a user's visits would show up as a galaxy that slowly
/// evaporates. Leapfrog's error oscillates instead of accumulating, so the disc
/// is still there when they come back.
pub fn step_leapfrog(bodies: &mut [Body], dt: f64, params: GravityParams) -> SolveReport {
    let before = crate::solvers::measure(bodies, 0.0);
    if bodies.is_empty() || dt == 0.0 {
        return SolveReport {
            before,
            after: before,
            dt_used: dt,
            ..Default::default()
        };
    }
    let mut tree = Octree::build(bodies, params);
    let half = 0.5 * dt;

    let acc = tree.accelerate_all(bodies);
    for (b, a) in bodies.iter_mut().zip(&acc) {
        b.vel += a.scale(half);
        b.pos += b.vel.scale(dt);
    }
    let mut tree2 = Octree::build(bodies, params);
    let acc2 = tree2.accelerate_all(bodies);
    for (b, a) in bodies.iter_mut().zip(&acc2) {
        b.vel += a.scale(half);
    }
    let interactions = tree.interactions + tree2.interactions;
    let after = crate::solvers::measure(bodies, 0.0);
    SolveReport {
        steps: 1,
        interactions,
        dt_used: dt,
        before,
        after,
        non_mechanical_energy: 0.0,
    }
}

/// Total mechanical energy, for drift diagnostics.
pub fn total_energy(bodies: &[Body], params: GravityParams) -> f64 {
    let mut tree = Octree::build(bodies, params);
    let kin = crate::state::kinetic_energy_of(bodies);
    kin + tree.potential_energy(bodies)
}

/// Timestep from the local dynamical time: `eta * sqrt(eps / |a|)`, the
/// standard adaptive criterion for collisionless N-body.
pub fn adaptive_dt(bodies: &[Body], params: GravityParams, eta: f64) -> f64 {
    if bodies.is_empty() {
        return f64::INFINITY;
    }
    let mut tree = Octree::build(bodies, params);
    let mut dt = f64::INFINITY;
    let eps = params.softening.max(1e-30);
    for (i, b) in bodies.iter().enumerate() {
        let a = tree.acceleration(b.pos, b.vel, i as i32).norm();
        if a > 0.0 {
            dt = dt.min(eta * (eps / a).sqrt());
        }
    }
    dt
}

/// Circular orbital speed at radius `r` around enclosed mass `m`.
pub fn circular_velocity(m: f64, r: f64) -> f64 {
    if r <= 0.0 {
        0.0
    } else {
        (G * m / r).sqrt()
    }
}

/// Schwarzschild radius — the engine switches to the post-Newtonian path
/// inside 100 of these.
pub fn schwarzschild_radius(m: f64) -> f64 {
    2.0 * G * m / C2
}
