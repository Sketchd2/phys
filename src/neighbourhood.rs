//! Which things are next to which.
//!
//! The engine models what happens *inside* a node very well and what happens
//! *between* nodes barely at all, and `docs/BACKLOG.md`'s coupling audit found
//! the reason: contact, heat conduction, mass diffusion, sibling-to-sibling
//! radiation, friction and debris landing on anything but its own parent are
//! not six missing features. They are one missing primitive — *two adjacent
//! things exchange a conserved quantity across the boundary they share* — and
//! building the six first would have produced six incompatible answers.
//!
//! This module is that primitive. See `docs/PLAY.md` D3.
//!
//! # Why the grid moved here
//!
//! [`NeighbourGrid`] lived in `solvers::hydro` and was already being imported
//! by `solvers::md`, which is the usual sign that something is not what its
//! address says it is. A uniform spatial hash is not a hydrodynamics concept;
//! it is how anything finds what is near it. Nothing about it changed in the
//! move except that it now indexes points rather than only bodies, because a
//! node's contents are bodies *and* promoted children.

use crate::math::Vec3;
use crate::state::Body;
use std::collections::HashMap;

/// Uniform spatial hash for neighbour finding — O(n) build, O(1) query.
///
/// A query touches 27 cells and nothing outside them can be a neighbour, so
/// the caller picks a spacing at least as large as the reach it will ask
/// about. SPH passes `2h`, the kernel support radius; molecular dynamics
/// passes its cutoff; a [`Neighbourhood`] passes the node's own resolution.
///
/// It indexes *points*, not bodies. That is not a generalisation for its own
/// sake: a node's contents are its materialised bodies **and** its promoted
/// children, and those are the same thing seen at two resolutions rather than
/// two kinds of thing. An index that could only hold bodies would have forced
/// adjacency to be written twice.
pub struct NeighbourGrid {
    cells: HashMap<(i64, i64, i64), Vec<u32>>,
    spacing: f64,
}

impl NeighbourGrid {
    /// Index a list of points.
    pub fn of_points(points: impl IntoIterator<Item = Vec3>, spacing: f64) -> NeighbourGrid {
        let mut cells: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
        let s = spacing.max(1e-30);
        for (i, p) in points.into_iter().enumerate() {
            cells.entry(key_of(p, s)).or_default().push(i as u32);
        }
        NeighbourGrid { cells, spacing: s }
    }

    /// Index a body list by position. The indices returned by
    /// [`Self::neighbours`] are then indices into `bodies`.
    pub fn build(bodies: &[Body], spacing: f64) -> NeighbourGrid {
        NeighbourGrid::of_points(bodies.iter().map(|b| b.pos), spacing)
    }

    /// The spacing the grid was built at. A query wider than this reaches
    /// beyond the 27 cells it visits and would miss neighbours, so callers that
    /// choose a radius at query time check it against this.
    pub fn spacing(&self) -> f64 {
        self.spacing
    }

    pub fn neighbours(&self, pos: Vec3, out: &mut Vec<u32>) {
        out.clear();
        let (cx, cy, cz) = key_of(pos, self.spacing);
        // The 27 cells are visited in a fixed order and each cell's contents
        // are in increasing body index (they were pushed in that order), so the
        // result is already deterministic and needs no sort. Sorting here — the
        // obvious way to guarantee determinism against HashMap iteration order,
        // which this loop never uses — cost more than the physics did: it was
        // 60% of the molecular dynamics runtime.
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(v) = self.cells.get(&(cx + dx, cy + dy, cz + dz)) {
                        out.extend_from_slice(v);
                    }
                }
            }
        }
    }

    pub fn occupied_cells(&self) -> usize {
        self.cells.len()
    }
}

/// The largest cell index the grid will name.
///
/// A body far enough from the origin bucketed naively produces an index that
/// does not fit in an `i64`, and the neighbour walk then adds one to it. In a
/// debug build that panics; in a release build it wraps to `i64::MIN`, the
/// lookup consults an arbitrary far-away cell, and that body's forces are
/// computed against whatever happens to be in it. Silently wrong physics is
/// the worse of the two, so the index is clamped instead — a quarter of the
/// range leaves the `+/- 1` walk room on both sides.
///
/// Clamping puts everything beyond the grid into one of eight corner cells.
/// That is a real loss of resolution and it is meant to be: a body 10^20 cells
/// from its own node has no neighbourhood, and pretending otherwise is what
/// this guard exists to stop. Nothing legitimate reaches it — a node's bodies
/// are node-relative and lie within a few radii of the origin.
const KEY_LIMIT: i64 = i64::MAX / 4;

#[inline]
fn axis_key(x: f64, s: f64) -> i64 {
    // `f64 as i64` saturates on overflow, which is what makes the clamp
    // sufficient, but it also maps NaN to *zero* — a real cell in the middle
    // of the grid, where a body with a NaN coordinate would silently become
    // everybody's neighbour. Send those off-grid with the rest.
    if !x.is_finite() {
        return KEY_LIMIT;
    }
    ((x / s).floor() as i64).clamp(-KEY_LIMIT, KEY_LIMIT)
}

#[inline]
fn key_of(p: Vec3, s: f64) -> (i64, i64, i64) {
    (axis_key(p.x, s), axis_key(p.y, s), axis_key(p.z, s))
}
