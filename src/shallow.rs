//! Water over a patch of ground, as a sheet: the shallow-water equations on
//! the patch's own columns — the owner's decision for Phase 5, over SPH.
//!
//! # Why a sheet
//!
//! `docs/PLAY.md` §4's beach wants waves flowing through a 5 cm channel in
//! real time. Drawn as parcels, a patch of shore could not hold still under a
//! still sea — an open boundary in weakly compressible SPH has nothing to
//! absorb what the region sends back — and cost 0.5 s of computing for 0.02 s
//! of world at 12.6 cm, where the beach wants a centimetre or so. Water a few
//! tens of centimetres deep over ground metres across is a sheet, and the
//! equations of a sheet are the ones the ocean already uses: depth and
//! depth-averaged velocity, at the patch's resolution instead of the planet's.
//! SPH stays for water that is not a sheet, such as a bucket or a splash.
//!
//! # What it is made of
//!
//! Nothing stored is a choice. The bed is the patch's own floor, measured off
//! its ordered members column by column down the field (`Floor`), at the
//! floor's own resolution; a channel carved into the ground is a change in
//! the bed and nothing else. The water is lent by the sea the patch stands
//! in (`ocean::Account`), and the edge is where the floor ends.
//!
//! # The scheme
//!
//! Finite volumes on the columns, first order:
//!
//! - **Hydrostatic reconstruction** (Audusse et al., 2004) at every face, so
//!   that still water over any bed stays still to rounding — a lake at rest
//!   is a steady state of the scheme, not an approximation to one — and a
//!   depth never goes negative, so water runs up a beach and down a channel
//!   with cells going dry and wet.
//! - **HLL** fluxes for mass and the normal momentum, the tangential momentum
//!   carried with the mass.
//! - **The seabed's log-law drag** against the ground's own grain, implicit,
//!   so friction cannot reverse a flow.
//! - **The edge is a characteristic boundary.** What leaves travels out on
//!   the invariant the sheet carries and what arrives is the sea's own — its
//!   level, current and wave train there — so a wave reflected from the beach
//!   leaves rather than ringing, and the sea's wave comes in.

use crate::math::Vec3;
use crate::ocean::{Sea, Train};
use crate::solvers::hydro::Wall;
use std::collections::HashMap;

/// A column index on the plane across the field.
pub type Column = (i64, i64);

/// Below this depth a column is dry, m: its velocity is not defined and is
/// held at zero.
pub const DRY: f64 = 1e-6;

/// The floor under a region's water, as columns `spacing` wide: the height of
/// the top of whatever solid each column meets, looking down the field.
#[derive(Debug, Clone)]
pub struct Floor {
    pub up: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub spacing: f64,
    /// A point of the lattice, on the plane through the region's centre.
    pub origin: Vec3,
    tops: HashMap<Column, f64>,
}

impl Floor {
    /// The floor `walls` make, down `gravity`, in columns `spacing` wide with
    /// one of them centred on `origin`.
    pub fn of(walls: &[Wall], gravity: Vec3, spacing: f64, origin: Vec3) -> Floor {
        let g = gravity.norm();
        let down = if g > 0.0 { gravity.scale(1.0 / g) } else { crate::math::v3(0.0, 0.0, -1.0) };
        let up = Vec3::ZERO - down;
        // The same basis the sampler's own floor uses.
        let seed = if down.x.abs() < 0.9 { crate::math::v3(1.0, 0.0, 0.0) } else { crate::math::v3(0.0, 1.0, 0.0) };
        let w = seed.cross(down);
        let u = w.scale(1.0 / w.norm());
        let v = down.cross(u);
        let origin = origin - up.scale(origin.dot(up));
        let mut tops: HashMap<Column, f64> = HashMap::new();
        for wall in walls {
            let reach = wall.reach();
            let (cu, cv) = ((wall.centre - origin).dot(u), (wall.centre - origin).dot(v));
            let (i0, i1) = (((cu - reach) / spacing).floor() as i64, ((cu + reach) / spacing).ceil() as i64);
            let (j0, j1) = (((cv - reach) / spacing).floor() as i64, ((cv + reach) / spacing).ceil() as i64);
            let (high, low) = (wall.centre.dot(up) + reach, wall.centre.dot(up) - reach);
            for i in i0..=i1 {
                for j in j0..=j1 {
                    let foot = origin + u.scale(i as f64 * spacing) + v.scale(j as f64 * spacing);
                    // Down the column from above the wall until it is met:
                    // each step the distance to it, which cannot overshoot.
                    let mut z = high;
                    let mut met = None;
                    for _ in 0..64 {
                        let (d, _) = wall.distance(foot + up.scale(z));
                        if d <= 1e-9 * spacing {
                            met = Some(z);
                            break;
                        }
                        z -= d;
                        if z < low {
                            break;
                        }
                    }
                    if let Some(top) = met {
                        let e = tops.entry((i, j)).or_insert(f64::NEG_INFINITY);
                        *e = e.max(top);
                    }
                }
            }
        }
        Floor { up, u, v, spacing, origin, tops }
    }

    /// The top of the floor in a column, or `None` where it has none.
    pub fn top(&self, c: Column) -> Option<f64> {
        self.tops.get(&c).copied()
    }

    /// The columns the floor covers, as a box of indices: `None` for a floor
    /// with none.
    pub fn extent(&self) -> Option<(Column, Column)> {
        let mut it = self.tops.keys();
        let first = *it.next()?;
        let (mut lo, mut hi) = (first, first);
        for &(i, j) in it {
            lo = (lo.0.min(i), lo.1.min(j));
            hi = (hi.0.max(i), hi.1.max(j));
        }
        Some((lo, hi))
    }
}

/// What crossed a sheet's edge, or left it by the bed, in one step or many.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Crossing {
    /// Mass that came in from the sea, kg (negative for what went out).
    pub mass: f64,
    /// Momentum that came in with it and with the sea's pressure, kg m/s, in
    /// the region's axes.
    pub momentum: Vec3,
    /// Its moment about the region's centre, kg m^2/s.
    pub moment: Vec3,
    /// Momentum the bed took, by its slope's push and its drag, kg m/s —
    /// from outside the water, as a wall's is (`SolveReport::outside`).
    pub bed: Vec3,
    /// Kinetic energy the bed's drag turned to heat, J. What the scheme
    /// dissipates at a bore — a breaking wave — is not counted here, and
    /// neither is angular momentum: a first-order finite-volume scheme
    /// carries neither exactly, and a sheet's books close on mass and
    /// momentum.
    pub heat: f64,
}

/// The water over a patch of ground.
#[derive(Debug, Clone)]
pub struct Sheet {
    /// Columns along `u` and `v`.
    pub nx: usize,
    pub ny: usize,
    /// Their width, m.
    pub dx: f64,
    /// The region's axes: across the field and up it.
    pub u: Vec3,
    pub v: Vec3,
    pub up: Vec3,
    /// The centre of column `(0, 0)`, on the plane through the region's
    /// centre.
    pub corner: Vec3,
    /// The bed's height in each column, m up the field from the region's
    /// centre; NaN where the floor has none, which is the sea beyond the edge.
    pub bed: Vec<f64>,
    /// Depth of water in each column, m.
    pub depth: Vec<f64>,
    /// Depth times velocity along `u` and `v`, m^2/s.
    pub flow: Vec<[f64; 2]>,
    /// The water's rest density, kg/m^3.
    pub density: f64,
    /// The seabed's roughness: the ground's own grain, m.
    pub grain: f64,
    /// Heat the bed's drag has made in the water since it was laid, J.
    pub heat: f64,
    /// The arrangement of the patch its bed was measured from (`Node::epoch`):
    /// a carve redraws the ground, and the bed follows it.
    pub epoch: u32,
}

/// One face's state on each side, reconstructed.
#[derive(Clone, Copy)]
struct Side {
    h: f64,
    /// Velocity normal to the face (positive from left to right) and along it.
    un: f64,
    ut: f64,
}

impl Sheet {
    /// A sheet on a floor, filled to `surface` m up the field from the
    /// region's centre wherever the bed is below it.
    pub fn on(floor: &Floor, surface: f64, density: f64, grain: f64) -> Option<Sheet> {
        let (lo, hi) = floor.extent()?;
        let (nx, ny) = ((hi.0 - lo.0 + 1) as usize, (hi.1 - lo.1 + 1) as usize);
        let mut bed = vec![f64::NAN; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                if let Some(t) = floor.top((lo.0 + i as i64, lo.1 + j as i64)) {
                    bed[j * nx + i] = t;
                }
            }
        }
        let depth = bed.iter().map(|b| if b.is_finite() { (surface - b).max(0.0) } else { 0.0 }).collect();
        let corner = floor.origin + floor.u.scale(lo.0 as f64 * floor.spacing) + floor.v.scale(lo.1 as f64 * floor.spacing);
        Some(Sheet {
            nx,
            ny,
            dx: floor.spacing,
            u: floor.u,
            v: floor.v,
            up: floor.up,
            corner,
            bed,
            depth,
            flow: vec![[0.0; 2]; nx * ny],
            density,
            grain,
            heat: 0.0,
            epoch: 0,
        })
    }

    /// Lay the sheet on a floor that has changed under it — the ground was
    /// carved, or filled — keeping the water each column holds: a channel
    /// cut into a flooded flat is empty the moment it is cut, and the water
    /// round it flows into it. A column the floor no longer covers gives its
    /// water back, which is what this returns, kg; the grid is the same
    /// floor's, so that is only ever a column at its edge.
    pub fn rebed(&mut self, floor: &Floor) -> f64 {
        let mut lost = 0.0;
        let a = self.density * self.dx * self.dx;
        for k in 0..self.bed.len() {
            let foot = self.foot(k) - floor.origin;
            let c = ((foot.dot(floor.u) / floor.spacing).round() as i64, (foot.dot(floor.v) / floor.spacing).round() as i64);
            match floor.top(c) {
                Some(t) => self.bed[k] = t,
                None => {
                    if self.holds(k) {
                        lost += a * self.depth[k];
                    }
                    self.bed[k] = f64::NAN;
                    self.depth[k] = 0.0;
                    self.flow[k] = [0.0, 0.0];
                }
            }
        }
        lost
    }

    fn index(&self, i: i64, j: i64) -> Option<usize> {
        (i >= 0 && j >= 0 && (i as usize) < self.nx && (j as usize) < self.ny).then(|| j as usize * self.nx + i as usize)
    }

    /// Whether a column is water's to hold: the floor covers it.
    pub fn holds(&self, k: usize) -> bool {
        self.bed[k].is_finite()
    }

    /// The centre of a column's footprint, in the region's axes.
    pub fn foot(&self, k: usize) -> Vec3 {
        let (i, j) = (k % self.nx, k / self.nx);
        self.corner + self.u.scale(i as f64 * self.dx) + self.v.scale(j as f64 * self.dx)
    }

    /// The column a point is over, if the sheet has one there.
    pub fn column_of(&self, p: Vec3) -> Option<usize> {
        let d = p - self.corner;
        let i = (d.dot(self.u) / self.dx).round() as i64;
        let j = (d.dot(self.v) / self.dx).round() as i64;
        self.index(i, j).filter(|&k| self.holds(k))
    }

    /// The water surface in a column, m up the field, where it is wet.
    pub fn surface(&self, k: usize) -> Option<f64> {
        (self.holds(k) && self.depth[k] > DRY).then(|| self.bed[k] + self.depth[k])
    }

    /// The water's depth-averaged velocity in a column, m/s, in the region's
    /// axes.
    pub fn velocity(&self, k: usize) -> Vec3 {
        let h = self.depth[k];
        if !(h > DRY) {
            return Vec3::ZERO;
        }
        self.u.scale(self.flow[k][0] / h) + self.v.scale(self.flow[k][1] / h)
    }

    /// Mass, kg.
    pub fn mass(&self) -> f64 {
        let a = self.dx * self.dx;
        (0..self.depth.len()).filter(|&k| self.holds(k)).map(|k| self.density * a * self.depth[k]).sum()
    }

    /// Momentum, kg m/s, in the region's axes.
    pub fn momentum(&self) -> Vec3 {
        let a = self.density * self.dx * self.dx;
        (0..self.depth.len())
            .filter(|&k| self.holds(k))
            .fold(Vec3::ZERO, |p, k| p + (self.u.scale(self.flow[k][0]) + self.v.scale(self.flow[k][1])).scale(a))
    }

    /// Angular momentum about the region's centre, kg m^2/s: each column's
    /// water at its own height.
    pub fn angular_momentum(&self) -> Vec3 {
        let a = self.density * self.dx * self.dx;
        (0..self.depth.len()).filter(|&k| self.holds(k)).fold(Vec3::ZERO, |l, k| {
            let at = self.foot(k) + self.up.scale(self.bed[k] + 0.5 * self.depth[k]);
            let p = (self.u.scale(self.flow[k][0]) + self.v.scale(self.flow[k][1])).scale(a);
            l + at.cross(p)
        })
    }

    /// Centre of mass, in the region's axes.
    pub fn centre_of_mass(&self) -> Vec3 {
        let a = self.density * self.dx * self.dx;
        let (mut m, mut x) = (0.0, Vec3::ZERO);
        for k in (0..self.depth.len()).filter(|&k| self.holds(k)) {
            let mk = a * self.depth[k];
            m += mk;
            x += (self.foot(k) + self.up.scale(self.bed[k] + 0.5 * self.depth[k])).scale(mk);
        }
        if m > 0.0 { x.scale(1.0 / m) } else { Vec3::ZERO }
    }

    /// Kinetic energy, J.
    pub fn kinetic_energy(&self) -> f64 {
        let a = self.density * self.dx * self.dx;
        (0..self.depth.len())
            .filter(|&k| self.holds(k) && self.depth[k] > DRY)
            .map(|k| 0.5 * a * (self.flow[k][0].powi(2) + self.flow[k][1].powi(2)) / self.depth[k])
            .sum()
    }

    /// The longest step the scheme is stable at, s.
    pub fn stable_step(&self, g: f64) -> f64 {
        let mut fastest: f64 = 0.0;
        for k in 0..self.depth.len() {
            if self.holds(k) && self.depth[k] > DRY {
                let h = self.depth[k];
                let s = (self.flow[k][0].abs().max(self.flow[k][1].abs())) / h + (g * h).sqrt();
                fastest = fastest.max(s);
            }
        }
        if fastest > 0.0 { 0.45 * self.dx / fastest } else { f64::INFINITY }
    }

    /// The state of the sea beyond an open face: its depth over the bed
    /// carried on, and its velocity, in the sheet's own `(u, v)`.
    fn beyond(&self, sea: &SeaAtEdge, at: Vec3, bed: f64, t: f64) -> (f64, [f64; 2]) {
        let place = sea.centre + at;
        let train = sea.train_at(sea.sea.level - (sea.height + bed));
        let surface = sea.sea.surface(train.as_ref(), place, t) - sea.height;
        let h = (surface - bed).max(0.0);
        let velocity = sea.sea.current + train.map(|w| w.mean_velocity(place, sea.sea.up, t)).unwrap_or(Vec3::ZERO);
        (h, [velocity.dot(self.u), velocity.dot(self.v)])
    }

    /// Advance by `dt` in a field of `g`, with the sea beyond the edge; `t`
    /// is the time the sea is read at. Returns what crossed.
    pub fn step(&mut self, dt: f64, g: f64, sea: &SeaAtEdge, t: f64) -> Crossing {
        let (nx, ny) = (self.nx as i64, self.ny as i64);
        let a = self.dx * self.dx;
        let rho = self.density;
        let mut out = Crossing::default();
        let mut dh = vec![0.0; self.depth.len()];
        let mut dq = vec![[0.0f64; 2]; self.depth.len()];
        // Along u (axis 0) then v (axis 1): each face once.
        for axis in 0..2usize {
            for j in 0..ny + (axis == 1) as i64 {
                for i in 0..nx + (axis == 0) as i64 {
                    let (li, lj) = if axis == 0 { (i - 1, j) } else { (i, j - 1) };
                    let (l, r) = (self.index(li, lj).filter(|&k| self.holds(k)), self.index(i, j).filter(|&k| self.holds(k)));
                    if l.is_none() && r.is_none() {
                        continue;
                    }
                    let normal = if axis == 0 { self.u } else { self.v };
                    let face = self.corner
                        + self.u.scale((li as f64 + if axis == 0 { 0.5 } else { 0.0 }) * self.dx)
                        + self.v.scale((lj as f64 + if axis == 1 { 0.5 } else { 0.0 }) * self.dx);
                    let side = |k: usize| -> (f64, f64, [f64; 2]) { (self.bed[k], self.depth[k], self.flow[k]) };
                    // A missing side is the sea beyond the edge, over the bed
                    // carried on from the side that is there.
                    let (zl, hl, ql) = match l {
                        Some(k) => side(k),
                        None => {
                            let z = self.bed[r.unwrap()];
                            let (h, v) = self.beyond(sea, face, z, t);
                            (z, h, [v[0] * h, v[1] * h])
                        }
                    };
                    let (zr, hr, qr) = match r {
                        Some(k) => side(k),
                        None => {
                            let z = self.bed[l.unwrap()];
                            let (h, v) = self.beyond(sea, face, z, t);
                            (z, h, [v[0] * h, v[1] * h])
                        }
                    };
                    let (mut hl, mut hr, mut ql, mut qr) = (hl, hr, ql, qr);
                    // At an open face, what the sea offers meets what the
                    // sheet sends out on the invariants each carries.
                    if l.is_none() || r.is_none() {
                        let inside_left = r.is_none();
                        let (hi, qi, ho, qo) = if inside_left { (hl, ql, hr, qr) } else { (hr, qr, hl, ql) };
                        let n = if inside_left { 1.0 } else { -1.0 };
                        let un_i = if hi > DRY { n * qi[axis] / hi } else { 0.0 };
                        let un_o = if ho > DRY { n * qo[axis] / ho } else { 0.0 };
                        let (ci, co) = ((g * hi).sqrt(), (g * ho).sqrt());
                        let out_inv = un_i + 2.0 * ci;
                        let in_inv = un_o - 2.0 * co;
                        let un = 0.5 * (out_inv + in_inv);
                        let c = (0.25 * (out_inv - in_inv)).max(0.0);
                        let hg = c * c / g;
                        let ut = if un > 0.0 {
                            if hi > DRY { qi[1 - axis] / hi } else { 0.0 }
                        } else if ho > DRY {
                            qo[1 - axis] / ho
                        } else {
                            0.0
                        };
                        let mut qg = [0.0; 2];
                        qg[axis] = n * un * hg;
                        qg[1 - axis] = ut * hg;
                        if inside_left {
                            hr = hg;
                            qr = qg;
                        } else {
                            hl = hg;
                            ql = qg;
                        }
                    }
                    // Hydrostatic reconstruction.
                    let zs = zl.max(zr);
                    let hl_s = (hl + zl - zs).max(0.0);
                    let hr_s = (hr + zr - zs).max(0.0);
                    let vel = |h: f64, q: [f64; 2]| if h > DRY { (q[axis] / h, q[1 - axis] / h) } else { (0.0, 0.0) };
                    let (unl, utl) = vel(hl, ql);
                    let (unr, utr) = vel(hr, qr);
                    let left = Side { h: hl_s, un: unl, ut: utl };
                    let right = Side { h: hr_s, un: unr, ut: utr };
                    let (fm, fn_, ft) = hll(left, right, g);
                    // The pressure the reconstruction took off each side's
                    // own column: what the bed's step pushes with.
                    let pl = 0.5 * g * (hl * hl - hl_s * hl_s);
                    let pr = 0.5 * g * (hr * hr - hr_s * hr_s);
                    let per = dt / self.dx;
                    if let Some(k) = l {
                        dh[k] -= per * fm;
                        dq[k][axis] -= per * (fn_ + pl);
                        dq[k][1 - axis] -= per * ft;
                        // The step's push on the left column is the bed's.
                        out.bed += normal.scale(-rho * a * per * pl);
                    }
                    if let Some(k) = r {
                        dh[k] += per * fm;
                        dq[k][axis] += per * (fn_ + pr);
                        dq[k][1 - axis] += per * ft;
                        out.bed += normal.scale(rho * a * per * pr);
                    }
                    // Across an open face, the sea's side of the books.
                    if l.is_none() || r.is_none() {
                        let n = if l.is_none() { 1.0 } else { -1.0 };
                        let mass = n * rho * self.dx * dt * fm;
                        let across = if axis == 0 { self.v } else { self.u };
                        let p = (normal.scale(fn_ + if l.is_none() { pr } else { pl }) + across.scale(ft)).scale(n * rho * self.dx * dt);
                        // The bed's step at the open face is the sea's: it is
                        // the sea's column the reconstruction is against.
                        let bed_there = normal.scale(if l.is_none() { rho * a * per * pr } else { -rho * a * per * pl });
                        out.bed -= bed_there;
                        out.mass += mass;
                        out.momentum += p;
                        out.moment += (face + self.up.scale(zs)).cross(p);
                    }
                }
            }
        }
        // The drag of the bed, implicit: the slope the bed makes is already in
        // the faces' reconstruction.
        for k in 0..self.depth.len() {
            if !self.holds(k) {
                continue;
            }
            let h = (self.depth[k] + dh[k]).max(0.0);
            let mut q = [self.flow[k][0] + dq[k][0], self.flow[k][1] + dq[k][1]];
            if h > DRY {
                let speed = (q[0] * q[0] + q[1] * q[1]).sqrt() / h;
                let cd = crate::ocean::log_law_drag(h, self.grain);
                let keep = 1.0 / (1.0 + dt * cd * speed / h);
                let lost = [q[0] * (1.0 - keep), q[1] * (1.0 - keep)];
                out.bed -= (self.u.scale(lost[0]) + self.v.scale(lost[1])).scale(rho * a);
                out.heat += 0.5 * rho * a * (q[0] * q[0] + q[1] * q[1]) * (1.0 - keep * keep) / h;
                q = [q[0] * keep, q[1] * keep];
            } else {
                q = [0.0, 0.0];
            }
            self.depth[k] = h;
            self.flow[k] = q;
        }
        self.heat += out.heat;
        out
    }
}

/// The HLL flux across a face from `l` to `r`: mass, normal momentum and
/// tangential momentum, per unit width and density.
fn hll(l: Side, r: Side, g: f64) -> (f64, f64, f64) {
    let flux = |s: Side| (s.h * s.un, s.h * s.un * s.un + 0.5 * g * s.h * s.h, s.h * s.un * s.ut);
    if !(l.h > DRY) && !(r.h > DRY) {
        return (0.0, 0.0, 0.0);
    }
    let (cl, cr) = ((g * l.h).sqrt(), (g * r.h).sqrt());
    let sl = (l.un - cl).min(r.un - cr);
    let sr = (l.un + cl).max(r.un + cr);
    let (fl, fr) = (flux(l), flux(r));
    if sl >= 0.0 {
        return fl;
    }
    if sr <= 0.0 {
        return fr;
    }
    let (ul, ur) = ((l.h, l.h * l.un, l.h * l.ut), (r.h, r.h * r.un, r.h * r.ut));
    let w = 1.0 / (sr - sl);
    let mass = (sr * fl.0 - sl * fr.0 + sl * sr * (ur.0 - ul.0)) * w;
    let normal = (sr * fl.1 - sl * fr.1 + sl * sr * (ur.1 - ul.1)) * w;
    // The tangential momentum goes with the mass, upwind.
    let tangential = if mass >= 0.0 { mass * l.ut } else { mass * r.ut };
    (mass, normal, tangential)
}

/// The sea round a region, as its sheet's edge reads it: the sea there
/// (`ocean::Sea`), where the region's centre is from the planet's, and how
/// far that centre stands above the sea's mean surface — all in the region's
/// axes. The train is carried into the depth at each face.
#[derive(Debug, Clone, Copy)]
pub struct SeaAtEdge {
    pub sea: Sea,
    pub centre: Vec3,
    pub height: f64,
}

impl SeaAtEdge {
    fn train_at(&self, depth: f64) -> Option<Train> {
        if depth > DRY { self.sea.train_in(depth) } else { None }
    }
}
