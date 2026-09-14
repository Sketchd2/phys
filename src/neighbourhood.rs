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

// ---------------------------------------------------------------------------
// what a node holds
// ---------------------------------------------------------------------------

/// One thing a node contains, at the finest resolution the node has for it.
///
/// A node's `children` runs parallel to its `bodies`: a promoted child *is* one
/// of those bodies, seen one level down. So a slot is one occupant, never two —
/// indexing both would make a thing its own neighbour and double every mass
/// that crossed a boundary.
///
/// Which of the two represents the slot is settled by the same rule promotion
/// is: the child is the real thing and the body is its stand-in, so where a
/// child exists it is what the neighbourhood holds.
/// Ordered so that a caller who needs a stable order can sort for one. The
/// query does not: `NeighbourGrid` visits its 27 cells in a fixed order and
/// each cell's contents are already in increasing index, so the result is
/// deterministic without sorting — and sorting it once measured 60% of the
/// molecular dynamics runtime, which is why it is the caller's choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Occupant {
    /// Index into the node's `bodies`.
    Body(u32),
    /// A promoted child, standing for the body in the same slot.
    Child(crate::ids::NodeIdx),
}

impl Occupant {
    /// The slot this occupant sits in, which is its index in either list.
    pub fn slot(self, children: &[crate::ids::NodeIdx]) -> Option<usize> {
        match self {
            Occupant::Body(i) => Some(i as usize),
            Occupant::Child(c) => children.iter().position(|x| *x == c),
        }
    }
}

/// What is next to what, inside one node.
///
/// Positions are in the node's own frame and distances in its own units, which
/// is what makes this scale-free: the same code indexes a galaxy's arms and a
/// nucleus's nucleons, because the numbers it sees are of order one either way.
/// There is no per-tier variant and there must not be one.
///
/// # The spacing, and why it is not simply the resolution
///
/// The grid's 27-cell walk finds everything within one cell of a point and
/// nothing beyond it, so the spacing has to be at least the largest distance
/// any query will ask about. A node's own resolution — its radius over the cube
/// root of its count, the same length SPH uses for a smoothing length — is the
/// natural scale, but an occupant larger than that would have its overlaps
/// missed. So the spacing is the greater of the resolution and twice the
/// largest occupant radius, which is the smallest spacing at which two touching
/// things are guaranteed to share or neighbour a cell.
///
/// The degenerate case is real and worth knowing: one occupant nearly as large
/// as the node forces a spacing that puts everything in a handful of cells, and
/// the query goes back to O(n). That is correct, merely slow, and it is what a
/// node holding one enormous thing and a thousand small ones actually deserves.
pub struct Neighbourhood {
    grid: NeighbourGrid,
    occupants: Vec<Occupant>,
    positions: Vec<Vec3>,
    radii: Vec<f64>,
    /// The node epoch this was built against. A node whose epoch has moved has
    /// different contents, and a neighbourhood built before it is stale.
    epoch: u32,
}

impl Neighbourhood {
    /// Build over a node's contents.
    ///
    /// `resolution` is the node's own, and `child_at(slot)` supplies the
    /// position and radius of the promoted child in that slot if there is one.
    /// Taking it as a closure keeps this module free of the tree: adjacency is
    /// about geometry, and which arena a child lives in is not its business.
    pub fn build(
        bodies: &[Body],
        children: &[crate::ids::NodeIdx],
        resolution: f64,
        epoch: u32,
        mut child_at: impl FnMut(crate::ids::NodeIdx) -> Option<(Vec3, f64)>,
    ) -> Neighbourhood {
        let n = bodies.len().max(children.len());
        let mut occupants = Vec::with_capacity(n);
        let mut positions = Vec::with_capacity(n);
        let mut radii = Vec::with_capacity(n);

        for slot in 0..n {
            let promoted = children
                .get(slot)
                .copied()
                .filter(|c| !c.is_none())
                .and_then(|c| child_at(c).map(|(p, r)| (c, p, r)));
            match promoted {
                Some((c, p, r)) => {
                    occupants.push(Occupant::Child(c));
                    positions.push(p);
                    radii.push(r.max(0.0));
                }
                None => {
                    let Some(b) = bodies.get(slot) else { continue };
                    occupants.push(Occupant::Body(slot as u32));
                    positions.push(b.pos);
                    radii.push(b.radius.max(0.0));
                }
            }
        }

        let widest = radii.iter().copied().fold(0.0f64, f64::max);
        let spacing = resolution.max(2.0 * widest).max(1e-30);
        let grid = NeighbourGrid::of_points(positions.iter().copied(), spacing);
        Neighbourhood { grid, occupants, positions, radii, epoch }
    }

    /// Whether this was built against a node in its present state.
    pub fn is_current(&self, epoch: u32) -> bool {
        self.epoch == epoch
    }

    pub fn len(&self) -> usize {
        self.occupants.len()
    }

    pub fn is_empty(&self) -> bool {
        self.occupants.is_empty()
    }

    /// The widest query this index can answer. Beyond it the 27-cell walk
    /// would miss neighbours, so a caller asking for more is told rather than
    /// quietly given a short answer.
    pub fn reach(&self) -> f64 {
        self.grid.spacing()
    }

    /// Everything whose *surface* lies within `within` of `point`.
    ///
    /// Radii are taken into account on the occupant's side, so a large thing is
    /// found by a query that would have missed its centre. Returns `None` if
    /// the query is wider than [`Self::reach`], because a short answer to a
    /// question the index cannot answer is worse than no answer.
    pub fn near(&self, point: Vec3, within: f64) -> Option<Vec<Occupant>> {
        if within > self.grid.spacing() {
            return None;
        }
        let mut candidates = Vec::new();
        self.grid.neighbours(point, &mut candidates);
        let mut out = Vec::new();
        for i in candidates {
            let i = i as usize;
            let gap = (self.positions[i] - point).norm() - self.radii[i];
            if gap <= within {
                out.push(self.occupants[i]);
            }
        }
        Some(out)
    }

    /// Everything overlapping the occupant at `index` — the impulsive case,
    /// where two things are not merely near but interpenetrating.
    pub fn touching(&self, index: usize) -> Vec<Occupant> {
        let Some(&me) = self.occupants.get(index) else { return Vec::new() };
        let (p, r) = (self.positions[index], self.radii[index]);
        let mut candidates = Vec::new();
        self.grid.neighbours(p, &mut candidates);
        let mut out = Vec::new();
        for j in candidates {
            let j = j as usize;
            if self.occupants[j] == me {
                continue;
            }
            if (self.positions[j] - p).norm() < r + self.radii[j] {
                out.push(self.occupants[j]);
            }
        }
        out
    }

    /// How many cells the contents actually landed in.
    ///
    /// The grid is a performance structure, not a correctness one: `near`
    /// filters by true distance afterwards, so a badly chosen spacing gives the
    /// right answer slowly rather than the wrong answer. That makes the failure
    /// invisible to a test that only checks results, which is why this is
    /// exposed — a spacing with a length baked into it collapses everything
    /// into one cell at some scales, and that is the symptom to assert on.
    pub fn cells(&self) -> usize {
        self.grid.occupied_cells()
    }

    /// The occupants, in slot order.
    pub fn occupants(&self) -> &[Occupant] {
        &self.occupants
    }
}
