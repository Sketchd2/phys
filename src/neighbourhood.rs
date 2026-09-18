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

    /// Everything whose *resilience* lies within `within` of `point`.
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
/// Not a shared resilience — nothing is shared, they are not touching. It is the
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

// ---------------------------------------------------------------------------
// Contact
// ---------------------------------------------------------------------------

/// How much of an impact a material gives back, and what it takes to stop it
/// giving any back at all.
///
/// Three numbers, all of them on `Material` — which is what `docs/PLAY.md` D3
/// means by "derived from the materials `topology.rs` already carries as data",
/// and which D14 then made derived rather than tabulated. Nothing here is a
/// coefficient of restitution or a coefficient of friction: those are
/// *results*, computed below from these and from how fast the two things are
/// closing.
///
/// **This was called `Surface`.** D18 gives that word to the thing a solid
/// presents geometrically — a union of solid convex primitives, each with its
/// own material — and no two things may share a name. This is not that: it is
/// the three numbers a *contact* needs from whichever primitive was struck, and
/// resilience is what a materials engineer calls the energy a material returns
/// rather than keeps.
///
/// That distinction is the whole point. A tabulated restitution is a frozen
/// answer to a question whose answer depends on the impact speed — the same
/// two blocks bounce at a walking pace and do not bounce when dropped from a
/// roof — so a table gets one of those two right and the rest of the range
/// wrong. `drop_fragments` carried `0.15` with the comment "Wood on wood: it
/// does not bounce", which was a reasonable guess for the one speed it was
/// tuned at.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Resilience {
    /// kg/m^3.
    pub density: f64,
    /// Young's modulus, Pa.
    pub stiffness: f64,
    /// The stress at which the contact stops returning what is put into it,
    /// Pa. Yield for a material that yields; fracture for one that breaks
    /// instead, because that is where a brittle contact stops being elastic.
    pub strength: f64,
}

impl Resilience {
    /// Read a resilience off a structural material.
    ///
    /// `ductility` is documented as "yield strength as a fraction of
    /// `strength()`, or zero for a brittle material that fractures instead of
    /// yielding", so a zero is not a material with no strength — it is a
    /// material whose elastic range ends at fracture.
    pub fn of(m: &crate::topology::Material) -> Resilience {
        let strength = if m.ductility > 0.0 {
            m.strength() * m.ductility
        } else {
            m.strength()
        };
        Resilience {
            density: m.density.max(0.0),
            stiffness: m.stiffness.max(0.0),
            strength: strength.max(0.0),
        }
    }

    fn is_usable(&self) -> bool {
        self.density > 0.0 && self.stiffness > 0.0 && self.strength > 0.0
    }
}

/// The closing speed above which a contact stops being elastic, m/s.
///
/// Derived rather than cited, because the cited constants disagree and the
/// derivation is short. Two spheres in Hertzian contact at approach `d` carry
/// `F = (4/3) E* sqrt(R*) d^(3/2)` over a circle of radius `sqrt(R* d)`, so the
/// mean contact pressure is `(4 E* / 3pi) sqrt(d / R*)`. Yield begins when that
/// reaches about `1.1 Y` — Johnson's result, the maximum shear under a Hertzian
/// circle being at `p_0 = 1.6 Y` and the mean being two thirds of the peak — so
///
/// ```text
///     d_y / R* = k^2,       k = (3 pi / 4) * 1.1 * (Y / E*) = 2.592 Y / E*
/// ```
///
/// The work stored up to that point is `(8/15) E* sqrt(R*) d_y^(5/2)`, which is
/// `(8/15) E* R*^3 k^5`. Setting it equal to `(1/2) m* v^2` for two equal
/// spheres — where `R* = R/2` and `m* = (2/3) pi R^3 rho` — the `R^3` cancels
/// on both sides, which is the interesting part: **the yield velocity does not
/// depend on how big the things are.** What is left is
///
/// ```text
///     v_y = 2.73 * sqrt( Y^5 / (E*^4 rho) )
/// ```
///
/// Checked against what it should reproduce: hardened steel comes out at
/// 0.22 m/s against a literature 0.1-0.2, and green wood on green wood at
/// 0.014 m/s, which puts a 5 m/s branch-fall at a restitution of 0.23 where the
/// hand-tuned constant it replaces was 0.15.
///
/// The `2.73` is for two equal spheres and is a scale rather than a precision
/// constant; the `(1 - nu^2)` plane-strain correction is left out because
/// `Material` carries no Poisson's ratio and it is worth about 10% per side.
pub fn yield_velocity(a: &Resilience, b: &Resilience) -> f64 {
    if !a.is_usable() || !b.is_usable() {
        return 0.0;
    }
    // Series stiffness: the softer side does most of the deflecting, which is
    // also why the softer side is the one that yields first.
    let e_star = 1.0 / (1.0 / a.stiffness + 1.0 / b.stiffness);
    // And the weaker side is the one that decides when the contact stops
    // returning energy, for the same reason.
    let y = a.strength.min(b.strength);
    let rho = 0.5 * (a.density + b.density);
    let ratio = y.powi(5) / (e_star.powi(4) * rho);
    if !ratio.is_finite() || ratio <= 0.0 {
        return 0.0;
    }
    2.73 * ratio.sqrt()
}

/// How much of the closing speed comes back, for this pair at this speed.
///
/// Elastic below the yield velocity and `(v_y / v)^(1/4)` above it — Johnson's
/// elastic-plastic result, and the exponent is the one thing here worth
/// remembering: restitution falls off *slowly*, so a contact that is 10,000
/// times past yield still returns a tenth of what it was given.
///
/// Speed-dependent by construction, which no tabulated coefficient can be.
pub fn restitution(a: &Resilience, b: &Resilience, closing: f64) -> f64 {
    let v_y = yield_velocity(a, b);
    let v = closing.abs();
    if !(v > 0.0) || !(v_y > 0.0) || v <= v_y {
        return 1.0;
    }
    (v_y / v).powf(0.25).clamp(0.0, 1.0)
}

/// Coulomb friction for a pair of surfaces.
///
/// Bowden and Tabor's adhesion account, and it is worth following because the
/// answer it gives is a *result* rather than a number: the real area of contact
/// is the load over the softer side's hardness, `A = W / H`; the junctions
/// formed there shear at that same side's shear strength; so
///
/// ```text
///     mu = tau A / W = tau / H = (Y / sqrt(3)) / (3 Y) = 1 / (3 sqrt(3)) = 0.192
/// ```
///
/// using von Mises for the shear yield and the standard `H ~ 3Y` for indentation
/// hardness. **The strength cancels.** Both terms come from the softer of the
/// two materials — it is the one that flows to make the junction and the one
/// that shears to break it — so the material drops out entirely, and that is
/// not a simplification made here but the reason most dry coefficients between
/// unlubricated solids sit between 0.2 and 0.5 whatever they are made of.
///
/// What the engine cannot see is what *does* vary: resilience films, roughness,
/// and the melt layer that makes ice 0.05 rather than 0.2. None of those are
/// represented, so none of them are guessed at. The signature still takes both
/// surfaces, because the day one of those is measurable this is where it goes.
pub fn friction(_a: &Resilience, _b: &Resilience) -> f64 {
    1.0 / (3.0 * 3.0f64.sqrt())
}

/// One side of a contact, at the moment of it.
#[derive(Debug, Clone)]
pub struct Side {
    /// Where this side's mass is. The lever arms of the contact couple are
    /// measured from here, so it is the centre of mass and not the centre of
    /// the shape — for a single sphere those coincide and for a wall of panels
    /// they need not.
    pub pos: Vec3,
    pub velocity: Vec3,
    pub mass: f64,
    /// Bounding radius about `pos`. The broad phase indexes this; the narrow
    /// phase does not use it, because `shape` says where the resilience actually
    /// is.
    pub radius: f64,
    /// J/K, for the heat the contact makes. See `state::Matter::heat_capacity`.
    pub heat_capacity: f64,
    pub resilience: Resilience,
    /// What this side actually *is*, geometrically, as one or more convex
    /// pieces. A lone body presents one sphere and gets exactly the arithmetic
    /// it always got; a structure presents a capsule per member and stops being
    /// a row of beads.
    pub shape: Vec<crate::shape::Hull>,
}

impl Side {
    /// The sphere case, which is every plain body and every promoted child.
    pub fn sphere(
        pos: Vec3,
        velocity: Vec3,
        mass: f64,
        radius: f64,
        heat_capacity: f64,
        resilience: Resilience,
    ) -> Side {
        Side {
            pos,
            velocity,
            mass,
            radius,
            heat_capacity,
            resilience,
            shape: vec![crate::shape::Hull::sphere(pos, radius)],
        }
    }

    /// A side whose geometry is a hull rather than its bounding sphere.
    ///
    /// `pos` stays the centre of mass: the hull says where the surfaces meet,
    /// and the couple is still taken about the mass.
    pub fn shaped(
        pos: Vec3,
        velocity: Vec3,
        mass: f64,
        heat_capacity: f64,
        resilience: Resilience,
        shape: Vec<crate::shape::Hull>,
    ) -> Side {
        let radius = shape
            .iter()
            .map(|h| (h.centre() - pos).norm() + h.bound())
            .fold(0.0f64, f64::max);
        Side { pos, velocity, mass, radius, heat_capacity, resilience, shape }
    }
}

/// What a contact does to the pair.
///
/// Expressed as what happens to **b**; `a` gets the negative of each impulse,
/// which is what makes momentum conservation structural rather than something
/// to be checked afterwards.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Collision {
    /// N s, along the line of centres.
    pub normal: Vec3,
    /// N s, across it.
    pub friction: Vec3,
    /// Angular momentum the friction couple puts into each side.
    pub spin_a: Vec3,
    pub spin_b: Vec3,
    /// J. The kinetic energy the contact did not give back, split so that both
    /// sides rise by the same temperature.
    pub heat_a: f64,
    pub heat_b: f64,
}

/// Resolve an overlap into an impulse pair.
///
/// Returns `None` when there is nothing to resolve: the two are separating
/// already, either side has no mass, or either side has no resilience to collide
/// with. The last is not a failure — a gas parcel and a star cluster are things
/// a node holds that have no resilience, and giving them one would be the engine
/// being *told* they are solid rather than measuring it.
///
/// # What is conserved, and how
///
/// Momentum: exactly, because one impulse is applied with both signs.
///
/// Angular momentum: exactly, and this is the part that is easy to drop.
/// Friction acts at the contact point and not at either centre, so it is a
/// couple as well as a force. Applying only the linear part changes
/// `sum r x p` without changing any spin, which loses angular momentum every
/// time anything slides. Applying the torque about each centre restores it
/// identically — `(c - p_a) x (-J) + (c - p_b) x J` is exactly the `(p_a - p_b) x J`
/// that the linear terms gained.
///
/// Energy: what the restitution did not return is not discarded, it is heat,
/// and it goes back into the two sides. Split by heat capacity rather than by
/// mass, because the heat is made at the interface and an interface has one
/// temperature; the exchange pass then moves it from there like any other heat.
pub fn contact(a: &Side, b: &Side) -> Option<Collision> {
    if !(a.mass > 0.0) || !(b.mass > 0.0) {
        return None;
    }
    if !a.resilience.is_usable() || !b.resilience.is_usable() {
        return None;
    }
    // Where the two surfaces actually meet. For a pair of spheres this is the
    // line of centres and the arithmetic below is unchanged; for anything else
    // it is the difference between hitting a wall and hitting whichever of its
    // panels happened to be nearest.
    let Some(near) = crate::shape::closest_of(&a.shape, &b.shape) else { return None };
    // The broad phase screens on bounding radii, which for a hull is a sphere
    // around the whole thing. Two walls whose bounds overlap and whose surfaces
    // do not are not in contact, and only the narrow phase knows.
    if near.gap > 0.0 {
        return None;
    }
    let n = near.normal;
    if !n.is_finite() {
        return None;
    }
    // One point, used by both sides. The midpoint of the two witness points is
    // the touching point exactly when the surfaces are just touching, and stays
    // on the overlap when they are not. Which point is chosen does not affect
    // conservation — the orbital and spin terms cancel for any of them, as
    // below — but it has to be the *same* point for both sides, which is the
    // thing the sphere path documents at length and the reason this is computed
    // once here.
    let point = (near.on_a + near.on_b).scale(0.5);
    let rel = b.velocity - a.velocity;
    let closing = rel.dot(n);
    // Positive means b is moving away from a. Nothing to resolve, and resolving
    // it anyway is how two overlapping things get stuck vibrating against each
    // other for the rest of their lives.
    if closing >= 0.0 {
        return None;
    }
    let reduced = 1.0 / (1.0 / a.mass + 1.0 / b.mass);
    let e = restitution(&a.resilience, &b.resilience, closing);
    let jn = -(1.0 + e) * closing * reduced;
    if !jn.is_finite() || jn <= 0.0 {
        return None;
    }

    // Tangential: arrest the sliding if friction can afford to, and slide at
    // the Coulomb limit if it cannot.
    let tangent_v = rel - n.scale(closing);
    let slide = tangent_v.norm();
    let (friction_impulse, jt) = if slide > 0.0 {
        let t = tangent_v.scale(1.0 / slide);
        let mu = friction(&a.resilience, &b.resilience);
        let jt = (slide * reduced).min(mu * jn);
        (t.scale(-jt), jt)
    } else {
        (Vec3::ZERO, 0.0)
    };

    let total_impulse = n.scale(jn) + friction_impulse;

    // The couple, about the one contact point both sides share.
    //
    // The obvious spelling — a's radius out from a, b's radius back from b — is
    // wrong, and wrong in a way that only shows up once the two are actually
    // interpenetrating. Those are two different points whenever the surfaces
    // overlap, and the total angular momentum then comes out as
    // `(dist - r_a - r_b) (n x J)` instead of zero: a 1.3% leak at a tenth of a
    // radius of overlap, which is what the test measured before this was fixed.
    //
    // With one point the spin terms contribute `(r_b - r_a) x J` and the
    // orbital terms `(b.pos - a.pos) x J`, and those are exact negatives for
    // *any* choice of point, so conservation does not depend on picking well.
    //
    // The **whole** impulse enters the couple, not only the friction part. For
    // two spheres the contact point lies on the line of centres, the normal
    // impulse is parallel to its own lever arm, and its cross product is
    // identically zero — which is why the sphere path could leave it out and
    // stay exact. Off the line of centres, which is where a hull puts a glancing
    // blow on a wall, the normal impulse turns what it hits, and omitting it
    // would lose that torque rather than cancel it.
    let ra = point - a.pos;
    let rb = point - b.pos;
    let spin_a = ra.cross(total_impulse.scale(-1.0));
    let spin_b = rb.cross(total_impulse);

    // Energy not returned. The normal direction loses `(1 - e^2)` of the
    // approach energy; the tangential direction loses whatever the friction
    // impulse took out of the sliding.
    let normal_loss = 0.5 * reduced * closing * closing * (1.0 - e * e);
    let slide_after = (slide - jt / reduced).max(0.0);
    let tangent_loss = 0.5 * reduced * (slide * slide - slide_after * slide_after);
    let heat = (normal_loss + tangent_loss).max(0.0);
    let (ca, cb) = (a.heat_capacity.max(0.0), b.heat_capacity.max(0.0));
    let total = ca + cb;
    let (heat_a, heat_b) = if total > 0.0 {
        (heat * ca / total, heat * cb / total)
    } else {
        (0.0, 0.0)
    };

    Some(Collision {
        normal: n.scale(jn),
        friction: friction_impulse,
        spin_a,
        spin_b,
        heat_a,
        heat_b,
    })
}
