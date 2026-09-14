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
    /// The length below which the node has no structure. Kept because callers
    /// need it and deriving it a second time from the node is how two
    /// definitions of one length drift apart.
    resolution: f64,
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
        Neighbourhood { grid, occupants, positions, radii, resolution, epoch }
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

    /// The node's own resolution — the length below which it has no structure.
    ///
    /// Distinct from [`Self::reach`], and the distinction matters to anything
    /// choosing a cutoff. `reach` is how far the *index* can see, which one
    /// oversized occupant can inflate until every pair is a candidate and the
    /// query is O(n^2). `resolution` is how far it is *meaningful* to look:
    /// below it the node has nothing to say, and above it the structure being
    /// described belongs to the parent.
    pub fn resolution(&self) -> f64 {
        self.resolution
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

    /// The occupant at an index, with where it is and how big it is.
    ///
    /// The index is into this neighbourhood, not into the node's bodies or its
    /// children — [`Occupant`] carries which of those it is. Keeping the two
    /// apart is what lets a promoted child and a plain body be handled by one
    /// piece of code.
    pub fn at(&self, i: usize) -> Option<(Occupant, Vec3, f64)> {
        Some((*self.occupants.get(i)?, self.positions[i], self.radii[i]))
    }

    /// Every pair whose surfaces lie within `within` of each other, each pair
    /// once, in a fixed order.
    ///
    /// The caller wants pairs, not a neighbour list per occupant: an exchange
    /// across a boundary has to be applied once or it moves twice as much of
    /// the quantity as it should, and deduplicating a per-occupant list at the
    /// call site is the kind of thing that is got right in one of the six
    /// callers D3 exists to unify.
    ///
    /// Ordered by `(i, j)` with `i < j`, which is deterministic without a sort
    /// for the same reason [`NeighbourGrid::neighbours`] is. Returns `None` if
    /// the query is wider than [`Self::reach`], on the same grounds as
    /// [`Self::near`].
    pub fn pairs(&self, within: f64) -> Option<Vec<(usize, usize)>> {
        if within > self.grid.spacing() {
            return None;
        }
        let mut out = Vec::new();
        let mut candidates = Vec::new();
        for i in 0..self.occupants.len() {
            self.grid.neighbours(self.positions[i], &mut candidates);
            for j in candidates.iter().map(|c| *c as usize) {
                if j <= i {
                    continue;
                }
                let gap =
                    (self.positions[j] - self.positions[i]).norm() - self.radii[i] - self.radii[j];
                if gap <= within {
                    out.push((i, j));
                }
            }
        }
        Some(out)
    }
}

// ---------------------------------------------------------------------------
// Exchange
// ---------------------------------------------------------------------------

/// One side of a shared boundary.
///
/// Deliberately not a `Body`, a `Matter` or a `Node`. Exchange does not care
/// what is on either side of the boundary — it cares how hard the quantity is
/// being pushed and how much of it one unit of push is worth. Heat sees a
/// temperature and a heat capacity; diffusing mass sees a concentration and a
/// volume; charge sees a voltage and a capacitance. Writing the transport three
/// times, once per caller, is exactly what `docs/BACKLOG.md`'s coupling audit
/// warned would produce three incompatible answers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reservoir {
    /// What drives the flow: temperature, concentration, potential.
    pub potential: f64,
    /// How much of the conserved quantity one unit of potential buys — heat
    /// capacity in J/K, volume for a concentration. `f64::INFINITY` is legal
    /// and means a bath: something whose potential the exchange cannot move.
    pub capacity: f64,
}

impl Reservoir {
    pub fn new(potential: f64, capacity: f64) -> Reservoir {
        Reservoir { potential, capacity }
    }

    /// A reservoir so large the exchange cannot move it — the sky, the ground,
    /// the cosmic microwave background.
    pub fn bath(potential: f64) -> Reservoir {
        Reservoir { potential, capacity: f64::INFINITY }
    }
}

/// Move a conserved quantity across a boundary. This is the one function.
///
/// Returns how much went **from `a` to `b`** over `dt`, positive when `a` was
/// the hotter/fuller side. Conduction, diffusion and radiative exchange are
/// three calls to it with different conductances; nothing else about them
/// differs, which is the whole claim of `docs/PLAY.md` D3.
///
/// # Why this is not `conductance * delta * dt`
///
/// That is the same expression, and it is what this returns in the limit of a
/// short step. But the engine's steps are not short by construction: a frame
/// covers a second, and two small things in good thermal contact equilibrate in
/// microseconds. The explicit form then moves more than the whole difference,
/// overshoots, and oscillates with growing amplitude — the classic stiff
/// failure, and it would have arrived as "a coffee cup on a table heats the
/// table to 900 K and then freezes it".
///
/// So the exact solution of the two-body problem is used instead. Two lumped
/// capacities coupled by a conductance obey `d(delta)/dt = -delta / tau` with
/// `tau = C_h / G` and `C_h = C_a C_b / (C_a + C_b)`, so the amount that
/// actually crosses in a span is `C_h * delta * (1 - exp(-dt / tau))`. That is
/// not a clamp bolted onto an unstable scheme: it is unconditionally stable
/// because it is *correct*, it conserves the quantity identically (what leaves
/// one side is what arrives at the other, by construction of a single scalar),
/// it can never carry the two sides past each other, and for small `dt` it
/// reduces to `G * delta * dt` exactly.
///
/// An infinite capacity on one side degrades gracefully to the one-body law
/// `C * delta * (1 - exp(-G dt / C))`, which is how a bath is written.
pub fn exchange(a: Reservoir, b: Reservoir, conductance: f64, dt: f64) -> f64 {
    let delta = a.potential - b.potential;
    if !(conductance > 0.0) || !(dt > 0.0) || !delta.is_finite() || delta == 0.0 {
        return 0.0;
    }
    // Harmonic mean of the two capacities: the quantity that can cross before
    // the potentials meet. Zero on either side means one of them cannot hold
    // any of the quantity, so nothing can cross.
    let (ca, cb) = (a.capacity, b.capacity);
    if !(ca > 0.0) || !(cb > 0.0) {
        return 0.0;
    }
    let harmonic = if ca.is_infinite() {
        cb
    } else if cb.is_infinite() {
        ca
    } else {
        ca * cb / (ca + cb)
    };
    if !harmonic.is_finite() || harmonic <= 0.0 {
        return 0.0;
    }
    // `-expm1(-x)` rather than `1 - exp(-x)`: for the small `x` of a brief step
    // between weakly coupled things the subtraction loses every significant
    // digit, and the transport silently stops happening at the point where it
    // is slowest — which is exactly where a slow leak matters.
    let reached = -(-(conductance * dt / harmonic)).exp_m1();
    harmonic * delta * reached
}

/// The area two spheres exchange radiation across, square metres.
///
/// Not a shared surface — nothing is shared, they are not touching. It is the
/// reciprocal area `A_a F_ab = A_b F_ba` that makes the grey-body law come out
/// symmetric, which is the thing that has to hold if the exchange is to
/// conserve energy: `pi r_a^2 r_b^2 / d^2`, the far-field limit of the
/// sphere-to-sphere view factor. Symmetric by construction rather than by
/// arithmetic luck, so neither side can be given a different answer than the
/// other about the same boundary.
///
/// The cap is the geometric bound for two spheres that do not interpenetrate:
/// at their closest the smaller one can present no more than a hemisphere to
/// the larger. It is almost never the binding term — equal spheres in contact
/// come out at an eighth of it — and that it is loose is the point. A cap that
/// bound often would be a fudge factor wearing a bound's clothes.
pub fn radiative_area(r_a: f64, r_b: f64, distance: f64) -> f64 {
    let d = distance.max(r_a + r_b).max(1e-300);
    // `(r_a r_b)^2` and not `r_a r_a r_b r_b`. The two are the same number in
    // exact arithmetic and *not* the same `f64`: the second rounds four times
    // in an order that depends on which side was named first, so a boundary
    // asked about from `a` and from `b` came back with areas differing in the
    // last bits. One conserved quantity crossing one boundary cannot have two
    // sizes, however small the difference — and a test asserting the two are
    // equal is what found this. A single product is commutative exactly,
    // because IEEE multiplication is, so squaring it is symmetric by
    // construction rather than by hoping the rounding cancels.
    let rr = r_a * r_b;
    let far = std::f64::consts::PI * rr * rr / (d * d);
    let cap = 2.0 * std::f64::consts::PI * r_a.min(r_b).powi(2);
    far.min(cap)
}

/// The conductance of a radiative boundary, watts per kelvin.
///
/// Stefan-Boltzmann is a fourth-power law and [`exchange`] wants a linear one,
/// so the usual move is to linearise about a mean temperature and accept the
/// error. There is no need: `T_a^4 - T_b^4` factors exactly as
/// `(T_a + T_b)(T_a^2 + T_b^2)(T_a - T_b)`, so dividing out the difference
/// leaves a conductance that reproduces the fourth-power law with no
/// approximation at all at the start of the step.
///
/// Emissivity is one, which is the same black body `state::stefan_boltzmann`
/// already assumes. A grey-body emissivity would have to come from somewhere,
/// and the only honest somewheres are a measurement or a table; the table is
/// forbidden and the measurement does not exist yet.
pub fn radiative_conductance(t_a: f64, t_b: f64, area: f64) -> f64 {
    if !(area > 0.0) || !t_a.is_finite() || !t_b.is_finite() {
        return 0.0;
    }
    let (a, b) = (t_a.max(0.0), t_b.max(0.0));
    crate::units::SIGMA_SB * area * (a + b) * (a * a + b * b)
}
