//! Where a region whose water is drawn as parcels meets the sea that describes
//! it: `docs/PLAY.md` §4.3's fifth piece, multi-resolution transport.
//!
//! # What the edge is
//!
//! A patch of shore holds its water as parcels, standing on its floor — the
//! solids its recipe states. Beyond the floor is the rest of the sea, which is
//! not drawn at all: it is the ocean on the planet, a level, a current and one
//! wave train (`ocean::Sea`). **The edge is where the floor ends**, measured
//! off the floor itself, column by column down the field: water over a column
//! the floor covers is the region's, and water anywhere else has left it.
//! Nothing states which side faces the sea. A side whose ground stands above
//! the sea's surface has no water beyond it and is an edge to nothing, and a
//! side under the sea is open.
//!
//! # What crosses it
//!
//! Three things, and each is paid for by the body that describes the ocean,
//! so the world's books close across the scale change:
//!
//! - **The sea's push.** Beyond the edge the sea is drawn as parcels for the
//!   region's parcels to feel (`hydro::Ghosts`): on the region's own lattice,
//!   from the bed carried on past the floor up to the sea's surface there,
//!   moving with the train and pressing with its pressure. They are drawn
//!   afresh each substep and never integrated, because the ocean is far larger
//!   than anything a patch of shore does to it.
//! - **Water leaving.** A parcel over a column the floor does not cover has
//!   gone back into the sea, with its mass, momentum and energy.
//! - **Water arriving.** A site on the region's lattice just inside the edge,
//!   under the sea's surface there and empty, is filled from the sea: a
//!   parcel of the region's own water, moving as the sea does there.

use crate::math::Vec3;
use crate::ocean::{Sea, Train};
use crate::solvers::hydro::Wall;
use crate::state::Body;
use std::collections::HashMap;

/// A column index on the plane across the field.
pub type Column = (i64, i64);

/// The floor under a region's water, as columns one parcel spacing wide: the
/// height of the top of whatever solid each column meets, looking down the
/// field.
#[derive(Debug, Clone)]
pub struct Floor {
    pub up: Vec3,
    u: Vec3,
    v: Vec3,
    pub spacing: f64,
    /// A point of the lattice the water lies on — one of its own parcels —
    /// on the plane through the centre, and the height it stood at.
    origin: Vec3,
    height: f64,
    tops: HashMap<Column, f64>,
}

impl Floor {
    /// The floor `walls` make, down `gravity`, in columns `spacing` wide on
    /// the lattice the water lies on, which `origin` — any one of its parcels
    /// — fixes.
    ///
    /// **The water's own lattice, not one this assumes.** The draw lays a
    /// liquid on multiples of its spacing and then moves the whole draw to
    /// put its centre of mass at the node's centre, so the lattice is off the
    /// multiples by whatever that took: measured, half a spacing across. On
    /// columns assumed at the multiples, 74 of 450 parcels at the edge of a
    /// patch of shore read as standing over nothing, and the sites drawn in
    /// to replace them sat seven tenths of a spacing from the parcels already
    /// there, which Tait prices at thousands of metres a second.
    pub fn of(walls: &[Wall], gravity: Vec3, spacing: f64, origin: Vec3) -> Floor {
        let g = gravity.norm();
        let down = if g > 0.0 { gravity.scale(1.0 / g) } else { crate::math::v3(0.0, 0.0, -1.0) };
        let up = Vec3::ZERO - down;
        // The same basis the sampler's own floor uses, so a column here is a
        // column of the draw there.
        let seed = if down.x.abs() < 0.9 { crate::math::v3(1.0, 0.0, 0.0) } else { crate::math::v3(0.0, 1.0, 0.0) };
        let w = seed.cross(down);
        let u = w.scale(1.0 / w.norm());
        let v = down.cross(u);
        let height = origin.dot(up);
        let origin = origin - up.scale(height);
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
        Floor { up, u, v, spacing, origin, height, tops }
    }

    /// The height of one of the water lattice's layers.
    pub fn origin_height(&self) -> f64 {
        self.height
    }

    /// The column a point is over.
    pub fn column(&self, p: Vec3) -> Column {
        let d = p - self.origin;
        ((d.dot(self.u) / self.spacing).round() as i64, (d.dot(self.v) / self.spacing).round() as i64)
    }

    /// Where a column's axis crosses the plane through the centre.
    pub fn foot(&self, c: Column) -> Vec3 {
        self.origin + self.u.scale(c.0 as f64 * self.spacing) + self.v.scale(c.1 as f64 * self.spacing)
    }

    /// The heights of the lattice's sites in a column above a bed: the lowest
    /// a quarter of a spacing clear of it, and up from there. The water's own
    /// sites stand half a spacing clear of a flat floor, and a quarter is
    /// what leaves room for a floor that is not flat.
    fn first_above(&self, bed: f64, lattice_height: f64) -> f64 {
        let s = self.spacing;
        lattice_height + ((bed + 0.25 * s - lattice_height) / s).ceil() * s
    }

    /// The top of the floor in a column, or `None` where it has none.
    pub fn top(&self, c: Column) -> Option<f64> {
        self.tops.get(&c).copied()
    }

    /// Whether the floor covers a column.
    pub fn covers(&self, c: Column) -> bool {
        self.tops.contains_key(&c)
    }

    /// The columns the floor does not cover within `reach` of one it does,
    /// each with the bed carried on past the edge: the top of the nearest
    /// column it covers. Deterministic order.
    pub fn beyond(&self, reach: f64) -> Vec<(Column, f64)> {
        let r = (reach / self.spacing).ceil() as i64;
        let mut out: HashMap<Column, (f64, f64)> = HashMap::new();
        for (&(i, j), &top) in &self.tops {
            for a in -r..=r {
                for b in -r..=r {
                    let c = (i + a, j + b);
                    if self.covers(c) {
                        continue;
                    }
                    let d = ((a * a + b * b) as f64).sqrt() * self.spacing;
                    if d > reach {
                        continue;
                    }
                    let e = out.entry(c).or_insert((f64::INFINITY, top));
                    if d < e.0 || (d == e.0 && top > e.1) {
                        *e = (d, top);
                    }
                }
            }
        }
        let mut list: Vec<(Column, f64)> = out.into_iter().map(|(c, (_, top))| (c, top)).collect();
        list.sort_by(|a, b| a.0.cmp(&b.0));
        list
    }

    /// The columns the floor covers that have one it does not beside them:
    /// the rim. Deterministic order.
    pub fn rim(&self) -> Vec<(Column, f64)> {
        let mut list: Vec<(Column, f64)> = self
            .tops
            .iter()
            .filter(|(&(i, j), _)| {
                (-1..=1).any(|a| (-1..=1).any(|b| (a, b) != (0, 0) && !self.covers((i + a, j + b))))
            })
            .map(|(&c, &t)| (c, t))
            .collect();
        list.sort_by(|a, b| a.0.cmp(&b.0));
        list
    }
}

/// The open edge of one region, for one solve: its floor, the sea beyond it,
/// and the columns on either side of the edge with the train each one's
/// water carries.
#[derive(Debug, Clone)]
pub struct Edge {
    pub floor: Floor,
    pub sea: Sea,
    /// How far the region's centre stands above the sea's mean surface, m.
    pub centre_height: f64,
    /// The region's centre from the planet's centre, in the region's axes,
    /// which is where the train's phase is read from.
    pub centre: Vec3,
    /// Columns beyond the edge, each with its bed and its train.
    pub beyond: Vec<(Column, f64, Option<Train>)>,
    /// Columns just inside it, likewise.
    pub rim: Vec<(Column, f64, Option<Train>)>,
    /// The height of one of the water lattice's layers.
    lattice_height: f64,
}

/// One substep's sea beyond an edge, as `hydro::Ghosts` wants it.
#[derive(Debug, Clone, Default)]
pub struct Beyond {
    pub bodies: Vec<Body>,
    pub pressure: Vec<f64>,
    pub density: Vec<f64>,
}

impl Edge {
    /// The edge of a region with this floor in this sea, reaching `reach`
    /// past the floor — the kernel's support, which is as far as a parcel
    /// feels.
    pub fn new(floor: Floor, sea: Sea, centre: Vec3, centre_height: f64, reach: f64) -> Edge {
        let lattice_height = floor.origin_height();
        // Each column's train is the sea's carried into that column's depth:
        // what reaches a shallow edge has shoaled and may have broken.
        let train_at = |bed: f64| sea.train_in(sea.level - (centre_height + bed));
        let beyond = floor.beyond(reach).into_iter().map(|(c, bed)| (c, bed, train_at(bed))).collect();
        let rim = floor.rim().into_iter().map(|(c, bed)| (c, bed, train_at(bed))).collect();
        Edge { floor, sea, centre_height, centre, beyond, rim, lattice_height }
    }

    /// Where a column's water is read from: its place from the planet's
    /// centre.
    fn place(&self, c: Column) -> Vec3 {
        self.centre + self.floor.foot(c)
    }

    /// The sea's surface over a column at `t`, as a height in the region's
    /// frame.
    fn surface(&self, c: Column, train: Option<&Train>, t: f64) -> f64 {
        self.sea.surface(train, self.place(c), t) - self.centre_height
    }

    /// Whether a parcel has left the region: it is over a column the floor
    /// does not cover.
    pub fn left(&self, p: Vec3) -> bool {
        !self.floor.covers(self.floor.column(p))
    }

    /// The sea beyond the edge at `t`, as parcels of `template`'s mass on the
    /// region's lattice, at the density `density_at` gives each pressure.
    pub fn beyond_at(&self, t: f64, template: &Body, density_at: &dyn Fn(f64) -> f64) -> Beyond {
        let s = self.floor.spacing;
        let mut out = Beyond::default();
        for (c, bed, train) in &self.beyond {
            let top = self.surface(*c, train.as_ref(), t);
            let foot = self.floor.foot(*c);
            let place = self.place(*c);
            let mut z = self.floor.first_above(*bed, self.lattice_height);
            while z < top {
                let height = z + self.centre_height;
                let p = self.sea.pressure(train.as_ref(), place, height, t);
                out.bodies.push(Body {
                    pos: foot + self.floor.up.scale(z),
                    vel: self.sea.velocity(train.as_ref(), place, height, t),
                    ..*template
                });
                out.pressure.push(p);
                out.density.push(density_at(p));
                z += s;
            }
        }
        out
    }

    /// Sites just inside the edge that are under the sea's surface at `t`,
    /// with the sea's velocity at each, in deterministic order. Which of them
    /// are empty is the caller's to say.
    pub fn inflow_sites(&self, t: f64) -> Vec<(Vec3, Vec3)> {
        let s = self.floor.spacing;
        let mut out = Vec::new();
        for (c, bed, train) in &self.rim {
            let top = self.surface(*c, train.as_ref(), t);
            let foot = self.floor.foot(*c);
            let place = self.place(*c);
            let mut z = self.floor.first_above(*bed, self.lattice_height);
            while z < top {
                let height = z + self.centre_height;
                out.push((foot + self.floor.up.scale(z), self.sea.velocity(train.as_ref(), place, height, t)));
                z += s;
            }
        }
        out
    }
}
