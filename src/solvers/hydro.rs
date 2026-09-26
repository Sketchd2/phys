//! Smoothed-particle hydrodynamics for the gas and continuum tiers.
//!
//! SPH rather than a grid, for one architectural reason: it is *meshless*, so
//! it composes with the tiers on either side of it without a remeshing step. A
//! gas parcel that condenses into a protostar is promoted to a `Planetary` node
//! by changing its interpretation, not by interpolating it onto a new grid —
//! and interpolation between grids is exactly where conservation goes to die.
//!
//! The formulation is the standard density-entropy one with Monaghan
//! artificial viscosity, which conserves momentum and angular momentum exactly
//! (forces are pairwise and antisymmetric) and energy to integrator accuracy.

use crate::eos::Eos;
use crate::math::Vec3;
use crate::neighbourhood::NeighbourGrid;
use crate::solvers::SolveReport;
use crate::state::Body;
use crate::units::*;

/// Cubic spline kernel (Monaghan & Lattanzio 1985), normalised in 3D.
#[inline]
pub fn kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        sigma * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}

/// Radial derivative of the kernel, `dW/dr` (negative inside the support).
#[inline]
pub fn kernel_grad(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        sigma * -0.75 * (2.0 - q).powi(2) / h
    } else {
        0.0
    }
}

#[derive(Debug, Clone, Copy)]
pub struct HydroParams {
    /// Smoothing length. Set from the mean interparticle spacing so that each
    /// particle has ~50 neighbours, the standard compromise between noise and
    /// resolution.
    pub h: f64,
    pub gamma: f64,
    /// Monaghan viscosity coefficients. Without these, shocks are not captured
    /// and supernova blast waves simply pass through each other.
    pub alpha: f64,
    pub beta: f64,
    /// Enable optically-thin radiative cooling. This is what lets gas collapse:
    /// without a cooling channel, compression heats gas until pressure stops
    /// it, and no star ever forms.
    pub cooling: bool,
    /// A uniform field the contents are held in, m/s^2, in the node's axes.
    ///
    /// Zero for a fluid region that nothing inside its own node holds up —
    /// its weight is carried by the fluid around it, which is outside the
    /// node — and the node's own `gravity` for contents resting on something
    /// the node holds. See `World::advance_node`.
    pub gravity: Vec3,
}

impl Default for HydroParams {
    fn default() -> Self {
        HydroParams {
            h: 1.0,
            gamma: 5.0 / 3.0,
            alpha: 1.0,
            beta: 2.0,
            cooling: true,
            gravity: Vec3::ZERO,
        }
    }
}

/// Something solid the fluid cannot enter, in the node's frame: a sphere of
/// `radius`, a capsule of `radius` about the segment `centre +- axis`, or a box
/// of `half`-extents turned by `orientation`. What a node's ordered members,
/// and the stand-ins for its solid children, present to its loose contents.
///
/// `owner` is the index of the body it moves with, among the bodies being
/// solved, for a wall that is itself one of them: a stand-in for a promoted
/// child. Its position follows that body every substep, the body is taken out
/// of the pairwise sums, and **every push the wall gives a parcel is given back
/// to its owner** — so a thing standing in water is held up by exactly the
/// pressure the water puts on it, which is what buoyancy is, and slowed by
/// exactly the drag, with nothing written for either.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wall {
    pub centre: Vec3,
    pub radius: f64,
    pub half: Vec3,
    pub orientation: crate::math::Quat,
    pub axis: Vec3,
    pub owner: Option<u32>,
}

impl Wall {
    /// A body as a wall: its box if it states one, its sphere otherwise.
    pub fn of(b: &Body) -> Wall {
        Wall {
            centre: b.pos,
            radius: b.radius,
            half: b.half,
            orientation: b.orientation,
            axis: Vec3::ZERO,
            owner: None,
        }
    }

    /// A member between two points: a structure's beam, which a body alone
    /// cannot state because a body is a sphere or a box.
    pub fn capsule(base: Vec3, tip: Vec3, radius: f64) -> Wall {
        Wall {
            centre: (base + tip).scale(0.5),
            radius,
            half: Vec3::ZERO,
            orientation: crate::math::Quat::IDENTITY,
            axis: (tip - base).scale(0.5),
            owner: None,
        }
    }

    /// How far from its centre anything of it reaches.
    pub fn reach(&self) -> f64 {
        self.radius.max(self.half.norm()) + self.axis.norm()
    }

    /// Signed distance from a point to the surface, and the outward normal.
    pub fn distance(&self, p: Vec3) -> (f64, Vec3) {
        let d = p - self.centre;
        if self.half == Vec3::ZERO {
            // A sphere, or a capsule: the distance to the nearest point of
            // its segment, less the radius.
            let a2 = self.axis.norm2();
            let along = if a2 > 0.0 { (d.dot(self.axis) / a2).clamp(-1.0, 1.0) } else { 0.0 };
            let off = d - self.axis.scale(along);
            let r = off.norm();
            let n = if r > 0.0 { off.scale(1.0 / r) } else { crate::math::v3(0.0, 0.0, 1.0) };
            return (r - self.radius, n);
        }
        let q = self.orientation.conjugate().rotate(d);
        let h = self.half;
        let e = crate::math::v3(q.x.abs() - h.x, q.y.abs() - h.y, q.z.abs() - h.z);
        let out = crate::math::v3(e.x.max(0.0), e.y.max(0.0), e.z.max(0.0));
        let outside = out.norm();
        let (dist, local) = if outside > 0.0 {
            (
                outside,
                crate::math::v3(out.x * q.x.signum(), out.y * q.y.signum(), out.z * q.z.signum())
                    .scale(1.0 / outside),
            )
        } else {
            // Inside: the nearest face.
            let m = e.x.max(e.y).max(e.z);
            let n = if m == e.x {
                crate::math::v3(q.x.signum(), 0.0, 0.0)
            } else if m == e.y {
                crate::math::v3(0.0, q.y.signum(), 0.0)
            } else {
                crate::math::v3(0.0, 0.0, q.z.signum())
            };
            (m, n)
        };
        (dist, self.orientation.rotate(local))
    }
}

/// One solve's walls, found by where they are.
///
/// **Indexed once a solve, not searched once a parcel.** Ground is thousands
/// of columns, and every parcel asked every one of them how far away it was,
/// twice a substep — and integrating the kernel over the solid asked each
/// point of the kernel the same of every wall in reach. Measured on a patch
/// of shore holding 450 parcels over 2500 columns: four seconds a solve. The
/// walls are found through a grid of their centres a cell wide enough that
/// anything reaching within a kernel's support of a point is in the 27 cells
/// round it, and whether a small cube of space is solid is worked out the
/// first time a kernel asks and remembered for the rest of the solve.
pub struct Solid {
    /// The walls nothing moves, by the cells their bounding boxes touch, a
    /// kernel's support wide: what is near a parcel.
    coarse: Cells,
    /// The walls that move with a body, which are few and are asked
    /// directly.
    moving: Vec<usize>,
    /// The kernel's points (`kernel_points`).
    points: Vec<(Vec3, f64)>,
    /// How far into or out of the walls nothing moves each point of a
    /// lattice round them is — negative inside — on the kernel points' own
    /// spacing: `step` apart, from the point at `corner`, `dims` of them each
    /// way.
    step: f64,
    corner: (i64, i64, i64),
    dims: (usize, usize, usize),
    distance: Vec<f64>,
}

/// Indices filed under every cell a box touches.
struct Cells {
    side: f64,
    map: std::collections::HashMap<(i64, i64, i64), Vec<u32>>,
}

impl Cells {
    fn key(&self, p: Vec3) -> (i64, i64, i64) {
        ((p.x / self.side).floor() as i64, (p.y / self.side).floor() as i64, (p.z / self.side).floor() as i64)
    }

    fn of(boxes: &[(usize, Vec3, Vec3)], side: f64) -> Cells {
        let mut cells = Cells { side, map: Default::default() };
        for (k, lo, hi) in boxes {
            let (a, b) = (cells.key(*lo), cells.key(*hi));
            for x in a.0..=b.0 {
                for y in a.1..=b.1 {
                    for z in a.2..=b.2 {
                        cells.map.entry((x, y, z)).or_default().push(*k as u32);
                    }
                }
            }
        }
        cells
    }

    /// Everything filed under a cell that the box `lo..hi` touches, once each.
    fn within(&self, lo: Vec3, hi: Vec3, out: &mut Vec<usize>) {
        let (a, b) = (self.key(lo), self.key(hi));
        for x in a.0..=b.0 {
            for y in a.1..=b.1 {
                for z in a.2..=b.2 {
                    if let Some(list) = self.map.get(&(x, y, z)) {
                        out.extend(list.iter().map(|&k| k as usize));
                    }
                }
            }
        }
        out.sort_unstable();
        out.dedup();
    }
}

impl Wall {
    /// The box that holds all of it, in the node's axes.
    fn bounds(&self) -> (Vec3, Vec3) {
        let e = if self.half == Vec3::ZERO {
            let r = self.radius;
            crate::math::v3(self.axis.x.abs() + r, self.axis.y.abs() + r, self.axis.z.abs() + r)
        } else {
            let q = self.orientation;
            let (a, b, c) = (
                q.rotate(crate::math::v3(self.half.x, 0.0, 0.0)),
                q.rotate(crate::math::v3(0.0, self.half.y, 0.0)),
                q.rotate(crate::math::v3(0.0, 0.0, self.half.z)),
            );
            crate::math::v3(
                a.x.abs() + b.x.abs() + c.x.abs(),
                a.y.abs() + b.y.abs() + c.y.abs(),
                a.z.abs() + b.z.abs() + c.z.abs(),
            )
        };
        (self.centre - e, self.centre + e)
    }
}

impl Solid {
    /// The walls of a solve at smoothing length `h`, indexed. **Once for all
    /// of a node's substeps**, not once a substep: the walls nothing moves do
    /// not move between them.
    pub fn new(walls: &[Wall], h: f64) -> Solid {
        let boxes: Vec<(usize, Vec3, Vec3)> = walls
            .iter()
            .enumerate()
            .filter(|(_, w)| w.owner.is_none())
            .map(|(k, w)| {
                let (lo, hi) = w.bounds();
                (k, lo, hi)
            })
            .collect();
        let moving = walls.iter().enumerate().filter(|(_, w)| w.owner.is_some()).map(|(k, _)| k).collect();
        let points = kernel_points(h);
        let step = 0.5 * h;
        // Every cube a kernel touching a fixed wall could ask about: the
        // walls' own bounds and a kernel's support round them.
        let (mut lo, mut hi) = (Vec3::ZERO, Vec3::ZERO);
        for (n, (_, a, b)) in boxes.iter().enumerate() {
            if n == 0 {
                (lo, hi) = (*a, *b);
            } else {
                lo = crate::math::v3(lo.x.min(a.x), lo.y.min(a.y), lo.z.min(a.z));
                hi = crate::math::v3(hi.x.max(b.x), hi.y.max(b.y), hi.z.max(b.z));
            }
        }
        let cube = |x: f64| (x / step).round() as i64;
        let (corner, dims, distance) = if boxes.is_empty() {
            ((0, 0, 0), (0, 0, 0), Vec::new())
        } else {
            let pad = 2.0 * h + step;
            let a = (cube(lo.x - pad), cube(lo.y - pad), cube(lo.z - pad));
            let b = (cube(hi.x + pad), cube(hi.y + pad), cube(hi.z + pad));
            let dims = ((b.0 - a.0 + 1) as usize, (b.1 - a.1 + 1) as usize, (b.2 - a.2 + 1) as usize);
            // Only the nearby walls matter to a point, and nothing past two
            // steps from a surface needs its distance exactly.
            let fine = Cells::of(&boxes, 2.0 * step);
            let far = 2.0 * step;
            let r = crate::math::v3(far, far, far);
            let mut distance = vec![far; dims.0 * dims.1 * dims.2];
            let mut near = Vec::new();
            for z in 0..dims.2 {
                for y in 0..dims.1 {
                    for x in 0..dims.0 {
                        let c = crate::math::v3((a.0 + x as i64) as f64, (a.1 + y as i64) as f64, (a.2 + z as i64) as f64).scale(step);
                        near.clear();
                        fine.within(c - r, c + r, &mut near);
                        let d = near.iter().map(|&k| walls[k].distance(c).0).fold(far, f64::min);
                        distance[(z * dims.1 + y) * dims.0 + x] = d.max(-far);
                    }
                }
            }
            (a, dims, distance)
        };
        Solid { coarse: Cells::of(&boxes, 2.0 * h), moving, points, step, corner, dims, distance }
    }

    /// The walls whose surfaces could be within `reach` of `p`.
    fn near(&self, p: Vec3, reach: f64) -> Vec<usize> {
        let r = crate::math::v3(reach, reach, reach);
        let mut out = self.moving.clone();
        self.coarse.within(p - r, p + r, &mut out);
        out
    }

    /// The share of a kernel at `p` inside the walls nothing moves, by
    /// integrating over them.
    fn integrated(&self, p: Vec3) -> f64 {
        self.points.iter().map(|(q, wt)| wt * self.solid_at(p + *q)).sum()
    }

    /// How much of the small cube of space round `p` is inside a wall nothing
    /// moves: the distance to them, interpolated between the lattice's
    /// points, read as a surface smoothed over one step. Read as inside or
    /// out by the nearest point instead, a floor came out a fifth of a kernel
    /// off the plane's closed form.
    fn solid_at(&self, p: Vec3) -> f64 {
        let far = 2.0 * self.step;
        let at = |x: i64, y: i64, z: i64| -> f64 {
            let (x, y, z) = (x - self.corner.0, y - self.corner.1, z - self.corner.2);
            if x < 0 || y < 0 || z < 0 || x as usize >= self.dims.0 || y as usize >= self.dims.1 || z as usize >= self.dims.2 {
                return far;
            }
            self.distance[(z as usize * self.dims.1 + y as usize) * self.dims.0 + x as usize]
        };
        let (fx, fy, fz) = (p.x / self.step, p.y / self.step, p.z / self.step);
        let (x0, y0, z0) = (fx.floor(), fy.floor(), fz.floor());
        let (tx, ty, tz) = (fx - x0, fy - y0, fz - z0);
        let (x0, y0, z0) = (x0 as i64, y0 as i64, z0 as i64);
        let mut d = 0.0;
        for (dx, wx) in [(0, 1.0 - tx), (1, tx)] {
            for (dy, wy) in [(0, 1.0 - ty), (1, ty)] {
                for (dz, wz) in [(0, 1.0 - tz), (1, tz)] {
                    d += wx * wy * wz * at(x0 + dx, y0 + dy, z0 + dz);
                }
            }
        }
        (0.5 - d / self.step).clamp(0.0, 1.0)
    }
}

/// The share of a kernel at `p` that lies inside the walls: what a wall
/// stands in place of in a density sum.
///
/// **One face, exactly; more than one, by integrating over the solid.** A
/// single face within reach is a plane to the kernel and `beyond_plane` is
/// its closed form — and exactness matters here more than anywhere: weakly
/// compressible, a liquid carries its weight on a per cent of its density, so
/// an error of a per cent in what a wall stands in for is an error the size
/// of the whole column's weight in the pressure. Where more than one face is
/// near, adding a plane for each counts every one of them as reaching away to
/// infinity, and a step two centimetres high in a floor of generated ground
/// was priced as a whole wall: measured, 5500 m/s^2 on a parcel standing in a
/// dip. There the kernel is integrated over the union itself, on
/// `kernel_points` against the distance to it (`Solid`), which is within
/// 0.0153 of a kernel of the closed form on a floor
/// (`a_floor_integrated_is_the_plane`) and which the closed form is kept
/// ahead of wherever it applies.
fn inside_walls(solid: &Solid, walls: &[Wall], p: Vec3, i: usize, h: f64) -> f64 {
    let near = facing(solid, walls, p, i, 2.0 * h);
    // A wall that moves with a body is that body's own surface, and is the
    // plane it always was: a thing floating in water is held up by exactly
    // what this gives it, and `one_ball_floats_half_under` is measured on it.
    let moving: f64 = near.iter().filter(|(k, _, _)| walls[*k].owner.is_some()).map(|(_, d, _)| beyond_plane(*d, h)).sum();
    let fixed: Vec<&(usize, f64, Vec3)> = near.iter().filter(|(k, _, _)| walls[*k].owner.is_none()).collect();
    moving
        + match fixed.len() {
            0 => 0.0,
            1 => beyond_plane(fixed[0].1, h),
            _ => solid.integrated(p),
        }
}

/// A kernel's support as points on a cubic lattice at a quarter of its
/// radius, each weighted by the kernel there, normalised so that the weights
/// sum to one — the kernel's own integral.
fn kernel_points(h: f64) -> Vec<(Vec3, f64)> {
    let step = 0.5 * h;
    let n = 4i32;
    let mut out = Vec::new();
    for x in -n..=n {
        for y in -n..=n {
            for z in -n..=n {
                let q = crate::math::v3(x as f64, y as f64, z as f64).scale(step);
                let w = kernel(q.norm(), h);
                if w > 0.0 {
                    out.push((q, w));
                }
            }
        }
    }
    let total: f64 = out.iter().map(|p| p.1).sum();
    out.into_iter().map(|(q, w)| (q, w / total)).collect()
}

/// The walls a body at `p` is within `reach` of, one per way the solid
/// faces it: `(wall, distance, outward normal)`, nearest first.
///
/// **A kernel sees each face of the solid once**, however many walls the
/// solid is built of. A floor of generated ground is thousands of columns a
/// few centimetres across, and a parcel standing on it is within reach of
/// twenty of them, every one facing it the same way; counted wall by wall,
/// that is twenty floors' worth of missing kernel in its density and twenty
/// springs under it. Measured over a patch of ground with water on it: a
/// density several times rest, which Tait prices at 10^8 Pa, and a parcel
/// flung at 1.2x10^5 m/s in the first substep. A wall that faces the body
/// within 45 degrees of a nearer one is the same face of the union and is
/// left out; a bucket's floor and its side, at right angles, are two faces
/// and are both counted, as they always were. A wall that moves with a body
/// is that body's own surface and always counts, so every push it gives
/// still comes back to it.
fn facing(solid: &Solid, walls: &[Wall], p: Vec3, i: usize, reach: f64) -> Vec<(usize, f64, Vec3)> {
    let mut near: Vec<(usize, f64, Vec3)> = solid
        .near(p, reach)
        .into_iter()
        .filter(|&k| walls[k].owner != Some(i as u32))
        .map(|k| {
            let (d, n) = walls[k].distance(p);
            (k, d, n)
        })
        .filter(|(_, d, _)| *d < reach)
        .collect();
    near.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let mut kept: Vec<(usize, f64, Vec3)> = Vec::with_capacity(near.len().min(6));
    for (k, d, n) in near {
        let moving = walls[k].owner.is_some();
        if moving || !kept.iter().any(|&(j, _, m)| walls[j].owner.is_none() && m.dot(n) > std::f64::consts::FRAC_1_SQRT_2) {
            kept.push((k, d, n));
        }
    }
    kept
}

/// The fraction of the kernel's mass that lies beyond a plane `d` from its
/// centre: `integral_{z > d} W dV`, for the cubic spline, as a function of
/// `d / h` alone.
///
/// **What a wall is to a density sum.** A parcel beside a wall has no
/// neighbours on the wall's side, so the sum over its neighbours comes out
/// short by exactly this much of a full kernel, and Tait's floor turns that
/// deficit into zero pressure. Measured, in a bucket two layers deep: the
/// bottom layer read below rest, pressed with nothing, and the top layer sank
/// into it until all 120 parcels lay in one sheet on the floor. A wall stands
/// where fluid would otherwise be, and counting it as the rest density it
/// displaces is the semi-analytic wall of Kulasegaram and Ferrand — derived
/// from the kernel, which is the only number in it. Half at contact, zero past
/// two smoothing lengths. The free surface is deliberately *not* given one,
/// because a deficit there is what makes it free.
pub fn beyond_plane(d: f64, h: f64) -> f64 {
    static TABLE: std::sync::OnceLock<Vec<f64>> = std::sync::OnceLock::new();
    const N: usize = 256;
    let table = TABLE.get_or_init(|| {
        // A(z) = 2 pi integral_{|z|}^{2} W(r, 1) r dr, then the tail of A.
        let steps = 2000;
        let dz = 2.0 / steps as f64;
        let area: Vec<f64> = (0..=steps)
            .map(|k| {
                let z = k as f64 * dz;
                let m = 400;
                let dr = (2.0 - z) / m as f64;
                let mut a = 0.0;
                for j in 0..m {
                    let r = z + (j as f64 + 0.5) * dr;
                    a += kernel(r, 1.0) * r * dr;
                }
                2.0 * std::f64::consts::PI * a
            })
            .collect();
        // tail[k] = integral_{z_k}^{2} A(z) dz, trapezoid from the far end.
        let mut tail = vec![0.0; steps + 1];
        for k in (0..steps).rev() {
            tail[k] = tail[k + 1] + 0.5 * (area[k] + area[k + 1]) * dz;
        }
        (0..=N).map(|i| tail[(i * steps) / N]).collect()
    });
    if !(h > 0.0) {
        return 0.0;
    }
    let q = (d / h).max(0.0);
    if q >= 2.0 {
        return 0.0;
    }
    let x = q / 2.0 * N as f64;
    let i = (x.floor() as usize).min(N - 1);
    let f = x - i as f64;
    table[i] * (1.0 - f) + table[i + 1] * f
}

/// How far the discrete kernel sum on a cubic lattice of `spacing` departs from
/// the density it is summing: `spacing^3 sum_n W(|n| spacing, h)`.
///
/// Derived, not tuned. A kernel integrates to one and a lattice sum of it does
/// not quite, and the difference is a pure function of `h / spacing`; for a
/// gas it is noise in a density that is free to be anything, and for a liquid
/// at rest it is a pressure, seven times over.
///
/// Measured, for the cubic spline, against `h / spacing`:
///
/// ```text
///   1.0  0.999972     1.5  1.001795     4.0  1.000008
///   1.2  1.000810     2.0  1.000079     6.0  1.000000
///   1.3  0.997262     3.0  0.999979
/// ```
///
/// Past four spacings it is the integral to eight parts in a million, and it
/// is taken as one there rather than summed: a stand-in for a 10 kg ball in a
/// 12 m room has a smoothing length of 110 of its own spacings, and summing
/// that lattice cost 221^3 kernel evaluations per body per substep.
pub fn lattice_sum(h: f64, spacing: f64) -> f64 {
    if !(h > 0.0) || !(spacing > 0.0) || h >= 4.0 * spacing {
        return 1.0;
    }
    let k = (2.0 * h / spacing).ceil() as i64;
    let mut sum = 0.0;
    for ix in -k..=k {
        for iy in -k..=k {
            for iz in -k..=k {
                let r = spacing * ((ix * ix + iy * iy + iz * iz) as f64).sqrt();
                sum += kernel(r, h);
            }
        }
    }
    sum * spacing.powi(3)
}

/// Densities by kernel summation.
pub fn densities(bodies: &[Body], params: HydroParams) -> Vec<f64> {
    let grid = NeighbourGrid::build(bodies, 2.0 * params.h);
    let mut out = vec![0.0; bodies.len()];
    let mut nb = Vec::with_capacity(128);
    for (i, b) in bodies.iter().enumerate() {
        grid.neighbours(b.pos, &mut nb);
        let mut rho = 0.0;
        for &j in nb.iter() {
            let o = &bodies[j as usize];
            let r = (b.pos - o.pos).norm();
            rho += o.mass * kernel(r, params.h);
        }
        out[i] = rho;
    }
    out
}

/// One SPH step: densities, pressures, forces, then a velocity-Verlet update.
///
/// Returns the energy explicitly radiated away, so the caller can subtract it
/// before checking conservation instead of quietly tolerating a drift.
pub fn step(bodies: &mut [Body], dt: f64, params: HydroParams) -> SolveReport {
    step_with(bodies, dt, params, &[], &[])
}

/// The equation of state body `i` answers to: `eos[i]`, or the gas law where
/// the slice does not reach — so an empty slice is every body a gas, which is
/// what [`step`] has always been.
#[inline]
fn eos_at(eos: &[Eos], i: usize) -> Eos {
    eos.get(i).copied().unwrap_or(Eos::Gas)
}

/// Pressure and sound speed of one body at a density.
fn state_of(b: &Body, rho: f64, eos: Eos, gamma: f64) -> (f64, f64) {
    match eos {
        Eos::Condensed(c) => (c.pressure(rho), c.sound_speed(rho)),
        Eos::Gas => {
            let mu = b.composition.mean_molecular_mass(b.temperature);
            let p = if mu > 0.0 { rho * K_B * b.temperature / mu } else { 0.0 }
                + A_RAD * b.temperature.powi(4) / 3.0;
            let cs = if rho > 0.0 { (gamma * p / rho).sqrt() } else { 0.0 };
            (p, cs)
        }
    }
}

/// [`step`], with an equation of state per body.
///
/// `docs/PLAY.md` Phase 5's second piece. A body whose entry is
/// [`Eos::Condensed`] is priced by Tait or Murnaghan rather than by the gas
/// law, and two more things change with it, both because they are gas laws too:
///
/// - **The compression work is carried.** A barotropic pressure does work that
///   the gas branch's temperature-based pressure never booked, and the books
///   close only if it goes somewhere: the standard SPH energy equation,
///   `du_i = (p_i / rho_i^2) sum_j m_j v_ij . grad W_ij`, puts it in the body's
///   internal energy. It is elastic, not heat, so it does not move the
///   temperature.
/// - **No optically thin cooling.** `cooling_rate` is a law for a gas thin
///   enough that every photon escapes. Asked about water it answers 10^20
///   W/m^3, and the half-the-energy cap below turned that into a liquid losing
///   half its internal energy every substep. Condensed matter is optically
///   thick and radiates from its surface, which `evolve_matter` already does.
///
/// And three things a liquid needs that a gas does not, which together are
/// weakly-compressible SPH (`docs/PLAY.md` §4.2) — the stiffness itself is the
/// caller's, since it is an equation of state:
///
/// - **A condensed body's density is summed over condensed neighbours**, and
///   divided by [`lattice_sum`] so that a liquid at rest reads its rest
///   density. Counting the air above it would put the free surface at the
///   wrong density, and counting the water below a gas parcel would give the
///   air a thousand times the pressure it has.
/// - **`params.gravity`** acts on every body.
/// - **`walls`** push a body out along their normal once it is within half a
///   parcel spacing (`h / 1.3`) of them — where a cell of liquid resting on a
///   floor has its centre, and where `sampler::packed_positions` lays it — at
///   the stiffness its own sound speed gives: `c^2 (g - d) / g^2` for a gap
///   `g`. Reaching out a whole spacing instead, which was the first version,
///   put every bottom parcel of a fresh draw half a spacing into a spring and
///   launched a lone one 2.8 m into the air.
///
/// Both of the last two do work on the contents from outside them, and that is
/// booked in `non_mechanical_energy` rather than left as drift.
pub fn step_with(
    bodies: &mut [Body],
    dt: f64,
    params: HydroParams,
    eos: &[Eos],
    walls: &[Wall],
) -> SolveReport {
    step_indexed(bodies, dt, params, eos, walls, None)
}

/// [`step_with`], with the walls indexed once for every substep of a solve
/// (see [`Solid`]; `None` indexes them for this step alone). The index must
/// have been made from these walls at this smoothing length.
///
/// A body with no mass is an empty slot — one a child was re-homed out of
/// (`Tree::reparent`) or that a sea's stand-in is left out through
/// (`World::advance_node`) — and takes no part.
pub fn step_indexed(
    bodies: &mut [Body],
    dt: f64,
    params: HydroParams,
    eos: &[Eos],
    walls: &[Wall],
    solid: Option<&Solid>,
) -> SolveReport {
    let before = crate::solvers::measure(bodies, 0.0);
    let n = bodies.len();
    if n == 0 || dt == 0.0 {
        return SolveReport {
            before,
            after: before,
            dt_used: dt,
            ..Default::default()
        };
    }
    let absent: Vec<bool> = bodies.iter().map(|b| !(b.mass > 0.0)).collect();

    // A body that is a wall is taken out of the pairwise sums: the fluid meets
    // it as a surface, not as a kernel's worth of its mass. See `Wall`.
    let mut walls: Vec<Wall> = walls.to_vec();
    let mut is_wall = vec![false; n];
    for w in walls.iter_mut() {
        if let Some(o) = w.owner {
            match bodies.get(o as usize) {
                // It goes where its body has gone.
                Some(b) => {
                    w.centre = b.pos;
                    is_wall[o as usize] = true;
                }
                None => w.owner = None,
            }
        }
    }
    let condensed: Vec<bool> = (0..n).map(|i| matches!(eos_at(eos, i), Eos::Condensed(_))).collect();
    let mut rho = if condensed.iter().any(|c| *c) {
        let grid = NeighbourGrid::build(bodies, 2.0 * params.h);
        let mut out = vec![0.0; n];
        let mut nb = Vec::with_capacity(128);
        for (i, b) in bodies.iter().enumerate() {
            if absent[i] {
                continue;
            }
            grid.neighbours(b.pos, &mut nb);
            out[i] = nb
                .iter()
                .map(|&j| j as usize)
                .filter(|&j| condensed[j] == condensed[i] && !is_wall[j])
                .map(|j| bodies[j].mass * kernel((b.pos - bodies[j].pos).norm(), params.h))
                .sum();
        }
        out
    } else {
        densities(bodies, params)
    };
    // The kernel's own bias on a lattice at the spacing the rest density
    // gives, per condensed body, and then the walls' share of its sum: see
    // `beyond_plane`.
    // Every parcel of one liquid shares a spacing, so the lattice sum is
    // worked out once per spacing rather than once per parcel per substep.
    let mut sums: Vec<(u64, f64)> = Vec::new();
    let own;
    let solid = match solid {
        Some(s) => s,
        None => {
            own = Solid::new(&walls, params.h);
            &own
        }
    };
    for i in 0..n {
        if absent[i] {
            continue;
        }
        if let Eos::Condensed(c) = eos_at(eos, i) {
            let spacing = (bodies[i].mass / c.rest_density).cbrt();
            let bits = spacing.to_bits();
            let sum = match sums.iter().find(|(b, _)| *b == bits) {
                Some((_, v)) => *v,
                None => {
                    let v = lattice_sum(params.h, spacing);
                    sums.push((bits, v));
                    v
                }
            };
            rho[i] /= sum;
            rho[i] += c.rest_density * inside_walls(solid, &walls, bodies[i].pos, i, params.h);
        }
    }
    let mut pressure = vec![0.0; n];
    let mut cs = vec![0.0; n];
    for i in 0..n {
        let (p, c) = state_of(&bodies[i], rho[i], eos_at(eos, i), params.gamma);
        pressure[i] = p;
        cs[i] = c;
    }

    let grid = NeighbourGrid::build(bodies, 2.0 * params.h);
    let mut acc = vec![Vec3::ZERO; n];
    let mut du = vec![0.0; n];
    // Compression work on condensed bodies, J/kg/s. See `step_with`.
    let mut dw = vec![0.0; n];
    let mut nb = Vec::with_capacity(128);
    let mut interactions = 0u64;

    for i in 0..n {
        grid.neighbours(bodies[i].pos, &mut nb);
        let bi = bodies[i];
        if rho[i] <= 0.0 || is_wall[i] || absent[i] {
            continue;
        }
        for &jj in nb.iter() {
            let j = jj as usize;
            if j == i || rho[j] <= 0.0 || is_wall[j] || absent[j] {
                continue;
            }
            let bj = bodies[j];
            let d = bi.pos - bj.pos;
            let r = d.norm();
            if r <= 0.0 || r >= 2.0 * params.h {
                continue;
            }
            let grad = kernel_grad(r, params.h);
            let dir = d.scale(1.0 / r);

            // Symmetric pressure form: the force on i from j is exactly minus
            // the force on j from i, so momentum and angular momentum are
            // conserved to machine precision rather than to truncation.
            let term = pressure[i] / (rho[i] * rho[i]) + pressure[j] / (rho[j] * rho[j]);

            // Monaghan artificial viscosity, active only in compression.
            let v_ij = bi.vel - bj.vel;
            let vr = v_ij.dot(d);
            let visc = if vr < 0.0 {
                let h = params.h;
                let mu_ij = h * vr / (r * r + 0.01 * h * h);
                let c_bar = 0.5 * (cs[i] + cs[j]);
                let rho_bar = 0.5 * (rho[i] + rho[j]);
                (-params.alpha * c_bar * mu_ij + params.beta * mu_ij * mu_ij) / rho_bar
            } else {
                0.0
            };

            let f = bj.mass * (term + visc) * grad;
            acc[i] += dir.scale(-f);
            // Viscous heating: the energy the viscosity removes from bulk
            // motion reappears as heat. Dropping this term is the classic way
            // to lose 10% of a shock's energy.
            du[i] += 0.5 * bj.mass * visc * grad * v_ij.dot(dir);
            if matches!(eos_at(eos, i), Eos::Condensed(_)) {
                dw[i] += pressure[i] / (rho[i] * rho[i]) * bj.mass * grad * v_ij.dot(dir);
            }
            interactions += 1;
        }
    }

    // The field, and the walls, as accelerations from outside — or, for a
    // wall that moves with a body, from that body, which takes the reaction.
    let gap = 0.5 * params.h / 1.3;
    let mut pushed = vec![Vec3::ZERO; n];
    // The dashpot's share, kept apart: its work is heat, booked in `du`, but
    // the momentum it takes is still from outside.
    let mut damped = vec![Vec3::ZERO; n];
    for i in 0..n {
        if absent[i] {
            continue;
        }
        acc[i] += params.gravity;
        if walls.is_empty() || is_wall[i] {
            continue;
        }
        let p = bodies[i].pos;
        for (k, d, normal) in facing(solid, &walls, p, i, gap) {
            let w = &walls[k];
            if d < gap {
                let c = cs[i].max(1e-30);
                let spring = normal.scale(c * c * (gap - d) / (gap * gap));
                // And the wall's share of the artificial viscosity. A spring
                // with nothing to damp it rings for ever: measured, a bucket's
                // bottom layer bounced at +-0.5 m/s for the whole of a second
                // while the parcels above it, which the pairwise viscosity
                // does reach, were still. A dashpot at the spring's own
                // frequency `c / g` — damping ratio a half — against the
                // velocity *relative to the wall*, and what it takes out of the
                // motion goes into the parcel as heat, exactly as the pairwise
                // term's does.
                let moving = w.owner.map(|o| bodies[o as usize].vel).unwrap_or(Vec3::ZERO);
                let vn = (bodies[i].vel - moving).dot(normal);
                let damp = normal.scale(-(c / gap) * vn);
                acc[i] += spring + damp;
                if w.owner.is_none() {
                    pushed[i] += spring;
                    damped[i] += damp;
                }
                du[i] += (c / gap) * vn * vn;
                if let Some(o) = w.owner {
                    let o = o as usize;
                    if bodies[o].mass > 0.0 {
                        acc[o] -= (spring + damp).scale(bodies[i].mass / bodies[o].mass);
                    }
                }
            }
        }
    }

    let mut radiated = 0.0;
    let mut external = 0.0;
    let mut impulse = Vec3::ZERO;
    for i in 0..n {
        if absent[i] {
            continue;
        }
        let b = &mut bodies[i];
        // Work done from outside: the field and the walls, less the pairwise
        // part, which is internal and conserves. Measured on the velocity the
        // step actually moves the body with.
        // Work done on the contents from outside them: the field, and the
        // push of a wall nothing in the solve owns. A moving wall's push is
        // internal — it is given back — and every dashpot's work is booked as
        // heat in `du`, which the conservation check already sees.
        let outside = params.gravity + pushed[i];
        b.vel += acc[i].scale(dt);
        external += b.mass * outside.dot(b.vel) * dt;
        impulse += (outside + damped[i]).scale(b.mass * dt);
        b.pos += b.vel.scale(dt);
        let heat = du[i] * dt * b.mass;
        b.internal_energy += heat + dw[i] * dt * b.mass;
        let mu = b.composition.mean_molecular_mass(b.temperature);
        let particles = if mu > 0.0 { b.mass / mu } else { 0.0 };
        if particles > 0.0 {
            b.temperature = (b.temperature + heat / (1.5 * particles * K_B)).max(2.725);
        }
        if params.cooling && rho[i] > 0.0 && !matches!(eos_at(eos, i), Eos::Condensed(_)) {
            let loss = cooling_rate(b.temperature, rho[i], b.composition.metallicity())
                * (b.mass / rho[i])
                * dt;
            let capped = loss.min(b.internal_energy.max(0.0) * 0.5);
            b.internal_energy -= capped;
            radiated += capped;
            if particles > 0.0 {
                b.temperature = (b.temperature - capped / (1.5 * particles * K_B)).max(2.725);
            }
        }
    }

    let after = crate::solvers::measure(bodies, 0.0);
    SolveReport {
        steps: 1,
        interactions,
        dt_used: dt,
        before,
        after,
        non_mechanical_energy: external - radiated,
        outside: impulse,
        // Including anything that is a wall: a ball put in water is out of
        // balance until the water holds it, and that is precisely the thing a
        // cadence has to see.
        unrest: (0..n).map(|i| acc[i].norm()).fold(0.0, f64::max),
    }
}

/// Optically-thin cooling rate in W/m^3.
///
/// A three-regime fit to the standard collisional-ionisation-equilibrium curve:
/// molecular/atomic line cooling below 10^4 K, the Lyman-alpha peak around
/// 10^4-10^5 K, and bremsstrahlung above 10^7 K. Metallicity scales the line
/// cooling, which is why the first generation of stars formed differently from
/// later ones — a difference this engine reproduces for free.
pub fn cooling_rate(t: f64, rho: f64, metallicity: f64) -> f64 {
    if t <= 10.0 {
        return 0.0;
    }
    let n = rho / M_PROTON; // number density, m^-3
    let lambda = if t < 1e4 {
        // Molecular/fine-structure cooling, strongly metallicity-dependent.
        1e-40 * (t / 100.0).powf(2.0) * (0.01 + metallicity * 30.0)
    } else if t < 1e7 {
        // Line cooling: peak near 10^5 K.
        let x = (t / 1e5).ln();
        1e-35 * (-x * x * 0.5).exp() * (0.1 + metallicity * 30.0)
    } else {
        // Free-free.
        2.3e-40 * t.sqrt()
    };
    lambda * n * n
}

/// Courant condition, including the viscous signal speed.
pub fn courant_dt(bodies: &[Body], params: HydroParams, cfl: f64) -> f64 {
    courant_dt_with(bodies, params, cfl, &[], &[])
}

/// [`courant_dt`], with an equation of state per body and the walls — as
/// [`step_with`].
///
/// **A body that is a wall is not a fluid parcel**, and its own sound speed is
/// a fact about its inside, which is a different node's business. Counting it
/// held a bucket's water to the step a wooden ball's 11 km/s asks for: 256
/// substeps a millisecond, and a second of floating that did not finish in ten
/// minutes. What a wall does ask of the step is its spring, `c / g` at the
/// parcels' own speed, which the parcels' term already carries.
pub fn courant_dt_with(bodies: &[Body], params: HydroParams, cfl: f64, eos: &[Eos], walls: &[Wall]) -> f64 {
    // A body in a field accelerates across its own smoothing length in
    // `sqrt(h / g)`, and a step longer than that lets it fall through
    // whatever holds it up.
    let g = params.gravity.norm();
    let mut dt = if g > 0.0 { cfl * (params.h / g).sqrt() } else { f64::INFINITY };
    for (i, b) in bodies.iter().enumerate() {
        if walls.iter().any(|w| w.owner == Some(i as u32)) {
            continue;
        }
        let rho = b.mass / (4.0 / 3.0 * std::f64::consts::PI * params.h.powi(3));
        let c = match eos_at(eos, i) {
            Eos::Condensed(c) => c.sound_speed(rho),
            Eos::Gas => {
                let mu = b.composition.mean_molecular_mass(b.temperature);
                let p = if mu > 0.0 { rho * K_B * b.temperature / mu } else { 0.0 };
                if rho > 0.0 { (params.gamma * p / rho).sqrt() } else { 0.0 }
            }
        };
        let v = b.vel.norm();
        let signal = c + v + 1e-30;
        dt = dt.min(cfl * params.h / signal);
    }
    dt
}

/// Jeans criterion: does this parcel have to be refined, or may it stay coarse?
///
/// A node whose Jeans length is unresolved is on the verge of collapsing into
/// structure the engine cannot see, and refusing to refine it is how a
/// simulation quietly produces a galaxy with no stars in it.
pub fn needs_refinement(rho: f64, temperature: f64, mu: f64, h: f64) -> bool {
    if rho <= 0.0 {
        return false;
    }
    let cs = (1.6667 * K_B * temperature / mu).sqrt();
    let jeans = cs * (std::f64::consts::PI / (G * rho)).sqrt();
    jeans < 4.0 * h
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The kernel integrated over a solid agrees with the plane's closed
    /// form**, for a floor, wherever a parcel can stand over it: what
    /// `Solid::integrated` gives a flat floor against `beyond_plane`, at
    /// heights from touching it to a kernel's support above it. Read as
    /// inside or out by the nearest point of the lattice instead of by the
    /// distance, it was 0.218 of a kernel off.
    #[test]
    fn a_floor_integrated_is_the_plane() {
        let h = 0.13;
        let floor = Wall {
            centre: crate::math::v3(0.0, 0.0, -1.0),
            radius: 0.0,
            half: crate::math::v3(3.0, 3.0, 1.0),
            orientation: crate::math::Quat::IDENTITY,
            axis: Vec3::ZERO,
            owner: None,
        };
        let walls = [floor];
        let solid = Solid::new(&walls, h);
        let mut worst: f64 = 0.0;
        for k in 0..=40 {
            let d = 2.0 * h * k as f64 / 40.0;
            // Off the cubes' own lattice, where the rounding is worst.
            let p = crate::math::v3(0.013, -0.029, d);
            let integrated = solid.integrated(p);
            let exact = beyond_plane(d, h);
            worst = worst.max((integrated - exact).abs());
        }
        println!("  a floor integrated on the kernel's points: worst {worst:.4} of a whole kernel off the plane's closed form");
        assert!(worst < 0.03, "{worst}");
    }
}
