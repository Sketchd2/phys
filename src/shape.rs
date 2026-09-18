//! Convex collision geometry, baked from spheres.
//!
//! The engine had exactly one collision shape on its general path — a sphere,
//! known by a centre and a radius — and `docs/BACKLOG.md` measured what that
//! costs: a 26 m settlement beam presents a 13.5 cm bead at its midpoint, 164
//! times shorter than the member the renderer draws. A ball thrown at a timber
//! frame passes through it unless it happens to hit a bead.
//!
//! # Why a set of spheres, and not a polytope
//!
//! The obvious spelling of "convex hull" is a list of vertices and faces. This
//! is not that, and the reason is that everything the engine already holds is a
//! sphere:
//!
//! - A materialised body is a centre and a radius. That is what `sample`
//!   produces and what `Matter` summarises back.
//! - A structural member is `base`, `tip` and a cross-section radius — which is
//!   a capsule, and a capsule is the hull of exactly two spheres.
//! - A parcel of fluid, a grain, a star: one sphere each.
//!
//! So a hull over a *set of spheres* covers all three without converting any of
//! them, and it does it exactly rather than approximately: the hull of two
//! spheres **is** the capsule, not a faceted stand-in for one. A polytope would
//! have to discard the radii and then re-inflate by some margin, which is the
//! same shape with two extra sources of error.
//!
//! It also means the bake stores nothing that was not already there. `PLAY.md`
//! axiom three allows a derived shortcut to be stored; a shape that is the
//! *source* of the geometry rather than a consequence of it is what the backlog
//! entry refuses, and a sphere set cannot become one — it is a view of the
//! bodies, and it is rebuilt when they change.
//!
//! # What varies is the grouping, not the primitive
//!
//! One hull over one sphere is that sphere. One hull over a member's two end
//! spheres is that capsule. One hull over the sixteen panels of a wall is that
//! wall — a slab with rounded edges, which is what a wall made of blocks
//! actually is. There is one shape type and one narrow phase, and what a caller
//! decides is *which occupants share a hull*.
//!
//! This is the part worth being careful about, because a hull is convex and
//! most interesting things are not. A box is six convex walls around a concave
//! cavity, and one hull over the whole box would enclose its own contents —
//! anything inside reads as deeply penetrating on every frame. So a hull is
//! never "the structure"; it is one convex piece of one, and a structure
//! presents as many of them as it has convex pieces.
//!
//! # The narrow phase is a distance query, and that is all contact needs
//!
//! [`closest`] runs GJK over the two hulls' *cores* — the convex hulls of their
//! sphere centres — and converts the result back to surfaces by subtracting the
//! radii at the witness points.
//!
//! There is deliberately no EPA here, and it is worth saying why, because its
//! absence looks like an omission. `neighbourhood::contact` resolves an overlap
//! from a normal, a contact point, the two masses and the closing velocity; it
//! never reads a penetration depth, and it does not push the pair apart
//! positionally. So the expensive half of a conventional narrow phase computes
//! a number nothing in this engine consumes. What is needed is the direction
//! and the point, and a distance query gives both.
//!
//! The one case a distance query cannot answer is two *cores* interpenetrating
//! — not two surfaces, which is the ordinary contact, but the centres of one
//! inside the hull of the other. That is a ball already halfway through a wall,
//! and [`Closest::degenerate`] reports it rather than inventing a normal.

use crate::math::{det_sum_v3_by, v3, Mat3, Vec3};

/// One sphere of a hull, in the frame the hull is expressed in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    pub centre: Vec3,
    pub radius: f64,
}

impl Sphere {
    pub fn new(centre: Vec3, radius: f64) -> Sphere {
        Sphere { centre, radius: radius.max(0.0) }
    }
}

/// A convex hull over a set of spheres.
///
/// Cheap to copy is not a goal — this is built once per node epoch and read
/// many times, which is what makes the bake worth doing at all.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Hull {
    /// The spheres that actually reach the surface. Interior ones are dropped
    /// by the bake, because a sphere swallowed by another can never be the
    /// answer to a support query and costs a comparison on every one.
    spheres: Vec<Sphere>,
    /// Centroid of the surviving centres. A broad-phase handle, not a centre of
    /// mass — mass lives on the bodies and [`Hull::mass_properties`] takes it
    /// as an argument rather than assuming the spheres carry it.
    centre: Vec3,
    /// Radius about `centre` that contains every sphere, surfaces included.
    bound: f64,
}

impl Hull {
    /// The empty hull. Collides with nothing, and [`closest`] says so.
    pub fn empty() -> Hull {
        Hull::default()
    }

    /// One sphere — what every loose body presents, and the fast path.
    pub fn sphere(centre: Vec3, radius: f64) -> Hull {
        Hull::of_spheres([Sphere::new(centre, radius)])
    }

    /// A capsule, which is the hull of its two end spheres exactly.
    ///
    /// This is what a structural member is: `Topology` carries `base`, `tip`
    /// and a cross-section radius, and those three numbers *are* a capsule. The
    /// engine has been storing them all along and collided the midpoint bead
    /// instead.
    pub fn capsule(base: Vec3, tip: Vec3, radius: f64) -> Hull {
        Hull::of_spheres([Sphere::new(base, radius), Sphere::new(tip, radius)])
    }

    /// Bake a hull over a set of spheres.
    ///
    /// Spheres wholly inside another are dropped. The test is containment in a
    /// *single* other sphere, which is sufficient rather than necessary: a
    /// sphere covered only by the combination of several others survives, and
    /// costs one comparison per support query for as long as it does. The exact
    /// test is a linear program per sphere, and paying that at bake time to
    /// save a handful of dot products at query time is the wrong trade — a wall
    /// of equal panels culls nothing either way, because none of its panels
    /// contains another.
    pub fn of_spheres(spheres: impl IntoIterator<Item = Sphere>) -> Hull {
        let all: Vec<Sphere> = spheres
            .into_iter()
            .filter(|s| s.centre.is_finite() && s.radius.is_finite() && s.radius >= 0.0)
            .collect();

        // The swallow scan below is O(n^2), and it is an *optimisation*: an
        // interior sphere can never answer a support query, so dropping it
        // saves one comparison per query and keeping it costs nothing but that.
        // Paying n^2 to save n per query stops being worth it quickly, and
        // measured it was catastrophic — a node of 4,000 sampled bodies took
        // 16 million distance computations to bake, and the frame went from
        // 25 ms to 121 ms.
        //
        // So above a threshold the scan is skipped and everything is kept.
        // Nothing downstream can tell the difference except in time.
        const SCAN_LIMIT: usize = 64;
        if all.len() > SCAN_LIMIT {
            let n = all.len();
            let centre = det_sum_v3_by(n, &|i| all[i].centre).scale(1.0 / n as f64);
            let bound = all
                .iter()
                .map(|s| (s.centre - centre).norm() + s.radius)
                .fold(0.0f64, f64::max);
            return Hull { spheres: all, centre, bound };
        }

        let mut kept: Vec<Sphere> = Vec::with_capacity(all.len());
        for (i, s) in all.iter().enumerate() {
            let swallowed = all.iter().enumerate().any(|(j, o)| {
                if i == j {
                    return false;
                }
                // Contained in `o`, and on an exact tie the later index yields
                // so that two identical spheres leave exactly one behind.
                let d = (s.centre - o.centre).norm();
                let inside = d + s.radius <= o.radius;
                let identical = d == 0.0 && s.radius == o.radius;
                if identical {
                    j < i
                } else {
                    inside
                }
            });
            if !swallowed {
                kept.push(*s);
            }
        }

        if kept.is_empty() {
            return Hull::default();
        }
        let n = kept.len();
        let centre = det_sum_v3_by(n, &|i| kept[i].centre).scale(1.0 / n as f64);
        let bound = kept
            .iter()
            .map(|s| (s.centre - centre).norm() + s.radius)
            .fold(0.0f64, f64::max);
        Hull { spheres: kept, centre, bound }
    }

    /// Move the whole hull, which is what expressing it in another node's frame
    /// costs. The shape is unchanged, so the bound is too.
    pub fn translated(&self, by: Vec3) -> Hull {
        Hull {
            spheres: self.spheres.iter().map(|s| Sphere::new(s.centre + by, s.radius)).collect(),
            centre: self.centre + by,
            bound: self.bound,
        }
    }

    /// Turn the whole hull about its own frame's origin, then move it.
    ///
    /// The composition matters and is the order a node's frame is defined in:
    /// `Node::collision_shape` is expressed about the node's own origin, and
    /// `Motion::body_to_parent` places a body-fixed point at
    /// `offset + orientation * local`. Translating without rotating is what
    /// `contact_within` used to do, and it left a spinning box's walls in the
    /// axes they were built in while the box's `orientation` turned — the
    /// second of the two rigid-body defects `PLAY.md` §2A found.
    ///
    /// A rotation is rigid, so the bound about the centre is unchanged and the
    /// spheres that survived the bake still survive it.
    pub fn placed(&self, orientation: crate::math::Quat, offset: Vec3) -> Hull {
        Hull {
            spheres: self
                .spheres
                .iter()
                .map(|s| Sphere::new(orientation.rotate(s.centre) + offset, s.radius))
                .collect(),
            centre: orientation.rotate(self.centre) + offset,
            bound: self.bound,
        }
    }

    pub fn spheres(&self) -> &[Sphere] {
        &self.spheres
    }

    pub fn is_empty(&self) -> bool {
        self.spheres.is_empty()
    }

    pub fn len(&self) -> usize {
        self.spheres.len()
    }

    /// Centroid of the surviving sphere centres.
    pub fn centre(&self) -> Vec3 {
        self.centre
    }

    /// Radius about [`Hull::centre`] containing the whole shape. This is what a
    /// broad phase indexes, and it is why `Neighbourhood` needs no change to
    /// carry hulls: it already indexes points with radii.
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// The furthest point of the shape along `dir`.
    ///
    /// Exact for the sphere-swept hull: the furthest point of a union of
    /// spheres along a direction is the furthest of their individual furthest
    /// points, and a sphere's is its centre plus its radius along `dir`.
    pub fn support(&self, dir: Vec3) -> Vec3 {
        let d = dir.unit();
        let mut best = Vec3::ZERO;
        let mut best_v = f64::NEG_INFINITY;
        for s in &self.spheres {
            let v = s.centre.dot(d) + s.radius;
            if v > best_v {
                best_v = v;
                best = s.centre + d.scale(s.radius);
            }
        }
        best
    }

    /// The furthest *centre* along `dir`, and its radius — the support of the
    /// core, which is what the narrow phase runs on.
    fn core_support(&self, dir: Vec3) -> (Vec3, f64) {
        let mut best = (Vec3::ZERO, 0.0);
        let mut best_v = f64::NEG_INFINITY;
        for s in &self.spheres {
            let v = s.centre.dot(dir);
            if v > best_v {
                best_v = v;
                best = (s.centre, s.radius);
            }
        }
        best
    }

    /// Centre of mass and inertia tensor for a rigid group presenting this
    /// hull, given a mass per sphere.
    ///
    /// Each sphere is treated as a solid sphere of its own mass about its own
    /// centre, plus the parallel-axis term. That is exact when the spheres are
    /// the bodies — which is the case this exists for, a wall's panels or a
    /// member's ends — and it is *not* an attempt to integrate the swept volume,
    /// which would double-count every overlap.
    ///
    /// Returns `None` when the masses do not line up with the spheres or sum to
    /// nothing, rather than guessing at either.
    pub fn mass_properties(&self, masses: &[f64]) -> Option<MassProperties> {
        if masses.len() != self.spheres.len() || self.spheres.is_empty() {
            return None;
        }
        let total: f64 = crate::math::det_sum(masses);
        if !(total > 0.0) || !total.is_finite() {
            return None;
        }
        let n = self.spheres.len();
        let com = det_sum_v3_by(n, &|i| self.spheres[i].centre.scale(masses[i]))
            .scale(1.0 / total);

        let mut inertia = Mat3::zero();
        for (s, &m) in self.spheres.iter().zip(masses) {
            if !(m > 0.0) {
                continue;
            }
            // Solid sphere about its own centre: (2/5) m r^2 on the diagonal.
            let own = 0.4 * m * s.radius * s.radius;
            // Parallel axis: m (|d|^2 I - d d^T).
            let d = s.centre - com;
            let shift = Mat3::identity().scaled(d.norm2()).add(d.outer(d).scaled(-1.0));
            inertia = inertia.add(Mat3::identity().scaled(own)).add(shift.scaled(m));
        }
        Some(MassProperties { mass: total, centre_of_mass: com, inertia })
    }
}

/// What a rigid group weighs and how it resists being turned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassProperties {
    pub mass: f64,
    pub centre_of_mass: Vec3,
    /// kg m^2, about `centre_of_mass`, in the hull's own frame.
    pub inertia: Mat3,
}

/// The closest approach of two hulls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Closest {
    /// Surface-to-surface separation. Negative means the surfaces overlap,
    /// which is the ordinary contact case.
    pub gap: f64,
    /// The witness point on each surface.
    pub on_a: Vec3,
    pub on_b: Vec3,
    /// Unit, pointing from `a` towards `b`.
    pub normal: Vec3,
    /// The cores interpenetrate, so the normal is a fallback along the line of
    /// centroids rather than a measured one. A caller that cares — anything
    /// deciding whether to trust the direction — checks this instead of
    /// discovering it from a contact that pushes the wrong way.
    pub degenerate: bool,
}

/// Iterations before GJK gives up. A converging query on shapes this small
/// takes single digits; the cap is here so a degenerate input cannot spin.
const GJK_MAX_ITERS: usize = 32;

/// Relative tolerance for "the support found nothing further out".
const GJK_TOLERANCE: f64 = 1e-12;

/// Closest approach of two hulls, by GJK over their cores.
///
/// Returns `None` only when either hull is empty, which is a shape that is not
/// there rather than a shape that does not touch.
///
/// # Why the cores, and what the radii do afterwards
///
/// Running GJK on the sphere-swept supports directly would give the surface
/// distance in one step, and it would also collapse the moment the surfaces
/// touch — which is every case contact is called for. Running it on the cores
/// instead keeps the query well away from its own degeneracy: two surfaces in
/// contact have cores a full radius apart, and the simplex stays a simplex.
///
/// The radii come back afterwards, at the witness points. Each witness is a
/// convex combination of sphere centres, so the radius there is the same
/// combination of their radii. That is exact when the closest feature is a
/// single sphere — the usual case, and always the case for a sphere against
/// anything — and exact whenever the spheres involved share a radius, which is
/// what a member's two ends and a wall's panels do. It is an interpolation only
/// for a closest feature spanning spheres of *different* radii, where the true
/// surface is a cone tangent to both and this reads it as a chord.
pub fn closest(a: &Hull, b: &Hull) -> Option<Closest> {
    if a.is_empty() || b.is_empty() {
        return None;
    }

    // The Minkowski difference of the two cores, and the simplex over it. Each
    // vertex remembers which side produced it, so the witness points fall out
    // of the same barycentric weights that locate the closest point.
    struct Vertex {
        d: Vec3,
        pa: Vec3,
        ra: f64,
        pb: Vec3,
        rb: f64,
    }
    let difference = |dir: Vec3| -> Vertex {
        let (pa, ra) = a.core_support(dir);
        let (pb, rb) = b.core_support(dir.scale(-1.0));
        Vertex { d: pa - pb, pa, ra, pb, rb }
    };

    // Any direction will do to start; the centroid separation is the one most
    // likely to be nearly right, and costs nothing.
    let seed = {
        let s = b.centre() - a.centre();
        if s.norm2() > 0.0 {
            s.scale(-1.0)
        } else {
            v3(1.0, 0.0, 0.0)
        }
    };
    let mut simplex: Vec<Vertex> = vec![difference(seed)];
    let mut weights: Vec<f64> = vec![1.0];

    // The length this query is measured in. "The origin is inside the
    // Minkowski difference" cannot be tested as an exact zero: the origin
    // inside a triangle projects to a point whose computed norm is whatever
    // the cancellation leaves — 10^-17 of the triangle, not 0.0 — and an exact
    // test falls straight through to a search direction of length nothing and
    // a normal made of rounding error. Measured on a ball inside a three-panel
    // core: 4.6x10^-34 for `dist2` against a shape 2 m across.
    let scale = (a.bound() + b.bound() + (b.centre() - a.centre()).norm()).max(1e-300);
    let inside2 = {
        let e = 1e-10 * scale;
        e * e
    };

    for _ in 0..GJK_MAX_ITERS {
        let points: Vec<Vec3> = simplex.iter().map(|v| v.d).collect();
        let Some((closest_point, w, keep)) = closest_on_simplex(&points) else {
            break;
        };
        // Reduce to the sub-simplex that actually supports the closest point,
        // or the simplex grows without bound and goes degenerate.
        let mut reduced: Vec<Vertex> = Vec::with_capacity(keep.len());
        let mut reduced_w: Vec<f64> = Vec::with_capacity(keep.len());
        for (slot, &k) in keep.iter().enumerate() {
            let v = &simplex[k];
            reduced.push(Vertex { d: v.d, pa: v.pa, ra: v.ra, pb: v.pb, rb: v.rb });
            reduced_w.push(w[slot]);
        }
        simplex = reduced;
        weights = reduced_w;

        let dist2 = closest_point.norm2();
        if dist2 <= inside2 || simplex.len() == 4 {
            // The origin is inside the Minkowski difference: the cores
            // themselves overlap.
            return Some(degenerate_contact(a, b));
        }

        // Search away from the origin, towards the shapes.
        let dir = closest_point.scale(-1.0);
        let next = difference(dir);
        // Converged when nothing lies further along the search direction than
        // the point already found.
        let progress = closest_point.dot(dir) - next.d.dot(dir);
        if progress >= -GJK_TOLERANCE * dist2.max(1.0) {
            break;
        }
        simplex.push(next);
    }

    // Witness points, from the weights the last reduction produced.
    let n = simplex.len();
    if n == 0 {
        return Some(degenerate_contact(a, b));
    }
    let total: f64 = crate::math::det_sum(&weights);
    if !(total > 0.0) {
        return Some(degenerate_contact(a, b));
    }
    let inv = 1.0 / total;
    let core_a = det_sum_v3_by(n, &|i| simplex[i].pa.scale(weights[i])).scale(inv);
    let core_b = det_sum_v3_by(n, &|i| simplex[i].pb.scale(weights[i])).scale(inv);
    let ra: f64 = crate::math::det_sum_by(n, &|i| simplex[i].ra * weights[i]) * inv;
    let rb: f64 = crate::math::det_sum_by(n, &|i| simplex[i].rb * weights[i]) * inv;

    let between = core_b - core_a;
    let d = between.norm();
    if !(d * d > inside2) || !d.is_finite() {
        return Some(degenerate_contact(a, b));
    }
    let normal = between.scale(1.0 / d);
    Some(Closest {
        gap: d - ra - rb,
        on_a: core_a + normal.scale(ra),
        on_b: core_b - normal.scale(rb),
        normal,
        degenerate: false,
    })
}

/// What to report when the two cores interpenetrate.
///
/// The line of centroids is the only direction left that is defined, and it is
/// flagged rather than presented as a measurement. A caller resolving this as
/// an ordinary contact would be pushing along a guess.
fn degenerate_contact(a: &Hull, b: &Hull) -> Closest {
    let between = b.centre() - a.centre();
    let d = between.norm();
    let normal = if d > 0.0 { between.scale(1.0 / d) } else { v3(1.0, 0.0, 0.0) };
    Closest {
        gap: d - a.bound() - b.bound(),
        on_a: a.centre() + normal.scale(a.bound()),
        on_b: b.centre() - normal.scale(b.bound()),
        normal,
        degenerate: true,
    }
}

/// The point of a simplex's convex hull closest to the origin.
///
/// Returns the point, the barycentric weights over the *surviving* vertices,
/// and which vertices those were. Every non-empty subset is tried and the
/// best valid one wins — valid meaning all weights non-negative, which is what
/// distinguishes a projection that lands inside a face from one that does not.
///
/// Johnson's recursive sub-distance algorithm does the same thing with fewer
/// operations and considerably more ways to be subtly wrong; with at most four
/// points there are fifteen subsets, and checking all of them is both exact and
/// obviously exact. The cost is a handful of 4x4 solves per GJK iteration on
/// shapes that have single-digit iterations.
fn closest_on_simplex(points: &[Vec3]) -> Option<(Vec3, Vec<f64>, Vec<usize>)> {
    let n = points.len();
    if n == 0 || n > 4 {
        return None;
    }
    let mut best: Option<(f64, Vec3, Vec<f64>, Vec<usize>)> = None;

    for mask in 1u32..(1 << n) {
        let idx: Vec<usize> = (0..n).filter(|i| mask & (1 << i) != 0).collect();
        let Some(w) = project_origin(points, &idx) else { continue };
        if w.iter().any(|x| *x < 0.0) {
            continue;
        }
        let p = det_sum_v3_by(idx.len(), &|k| points[idx[k]].scale(w[k]));
        let d2 = p.norm2();
        let better = match &best {
            None => true,
            // Fixed tie-break by subset size then by mask, so an origin
            // equidistant from two features always picks the same one. A query
            // that answered differently run to run would put a contact normal
            // on a coin flip.
            Some((bd2, _, bw, _)) => d2 < *bd2 || (d2 == *bd2 && idx.len() < bw.len()),
        };
        if better {
            best = Some((d2, p, w, idx));
        }
    }
    best.map(|(_, p, w, idx)| (p, w, idx))
}

/// Barycentric weights of the origin's projection onto the affine hull of the
/// selected points.
///
/// Solves the constrained least-squares system
///
/// ```text
///     [ 2G  1 ] [ lambda ]   [ 0 ]
///     [ 1^T 0 ] [   mu   ] = [ 1 ]
/// ```
///
/// where `G` is the Gram matrix of the points. `None` when the points are
/// degenerate — coincident or collinear within a face — which is a subset that
/// simply does not support a projection and is skipped rather than repaired.
fn project_origin(points: &[Vec3], idx: &[usize]) -> Option<Vec<f64>> {
    let k = idx.len();
    if k == 1 {
        return Some(vec![1.0]);
    }
    let m = k + 1;
    // Row-major (k+1) x (k+2) augmented system, m <= 5.
    let mut a = [[0.0f64; 6]; 5];
    for i in 0..k {
        for j in 0..k {
            a[i][j] = 2.0 * points[idx[i]].dot(points[idx[j]]);
        }
        a[i][k] = 1.0;
        a[i][m] = 0.0;
    }
    for j in 0..k {
        a[k][j] = 1.0;
    }
    a[k][k] = 0.0;
    a[k][m] = 1.0;

    // Gaussian elimination with partial pivoting. The pivot search takes the
    // first row of maximal magnitude, so the choice is deterministic on a tie.
    for col in 0..m {
        let mut pivot = col;
        for r in (col + 1)..m {
            if a[r][col].abs() > a[pivot][col].abs() {
                pivot = r;
            }
        }
        if a[pivot][col].abs() <= 1e-300 {
            return None;
        }
        if pivot != col {
            a.swap(pivot, col);
        }
        let d = a[col][col];
        for j in col..=m {
            a[col][j] /= d;
        }
        for r in 0..m {
            if r == col {
                continue;
            }
            let f = a[r][col];
            if f == 0.0 {
                continue;
            }
            for j in col..=m {
                a[r][j] -= f * a[col][j];
            }
        }
    }
    let w: Vec<f64> = (0..k).map(|i| a[i][m]).collect();
    if w.iter().any(|x| !x.is_finite()) {
        return None;
    }
    Some(w)
}

/// Closest approach of two *compound* shapes — the nearest pair of convex
/// pieces.
///
/// A structure is many convex pieces, and a contact between two of them is a
/// contact between one piece of each. Taking the minimum-gap pair resolves one
/// contact per pair of things per frame, which is the right answer for a ball
/// bouncing around inside a box and an approximation for a ball wedged into a
/// corner: the second wall is resolved on the following frame, by which time
/// the first has already turned the ball away from it.
///
/// Resolving *every* overlapping pair at once is the alternative, and it is not
/// obviously better — two normal impulses applied independently in one step
/// over-correct, which is the standard reason a solver iterates instead. That
/// is a scheduling question of the same kind §3.3 left open, and this does not
/// pre-empt it.
pub fn closest_of(a: &[Hull], b: &[Hull]) -> Option<Closest> {
    nearest_of(a, b).map(|(c, _, _)| c)
}

/// One bounding sphere over a set of hulls, in whatever frame they are in.
///
/// Every hull already knows its own centre and the radius that contains it, so
/// this is a sum and a max rather than anything geometric.
fn bounds(hulls: &[Hull]) -> Option<(Vec3, f64)> {
    if hulls.is_empty() {
        return None;
    }
    let n = hulls.len();
    let centre = det_sum_v3_by(n, &|i| hulls[i].centre()).scale(1.0 / n as f64);
    let bound = hulls
        .iter()
        .map(|h| (h.centre() - centre).norm() + h.bound())
        .fold(0.0f64, f64::max);
    Some((centre, bound))
}

/// As [`closest_of`], and says *which* pieces were nearest.
///
/// # The level of detail the distance deserves
///
/// `PLAY.md` §7's item 8. The N x M is linear in the piece count and flat per
/// piece — `PERFORMANCE.md` measures 0.31 µs each — so a generated tree of
/// 3,400 members is 1.1 ms for a single pair. That is a fiftieth of a frame for
/// one tree touching one thing, and a recipe is entitled to emit as many pieces
/// as the thing has.
///
/// **The saving here is exact rather than approximate**, which is what makes it
/// a level of detail rather than a fudge. A hull's `bound` contains it, so a
/// piece further from the other side's bounding sphere than the two bounds
/// together *cannot* be the nearest pair, and skipping it changes no answer. A
/// tree a hundred metres away costs one sphere test; a ball resting against one
/// of its branches costs the branches within reach of the ball and not the
/// three thousand that are not.
///
/// The returned indices are what a caller reads a per-piece material from: D18
/// puts the material on the piece, and a contact has to read the one it
/// actually struck.
pub fn nearest_of(a: &[Hull], b: &[Hull]) -> Option<(Closest, usize, usize)> {
    let (ca, ra) = bounds(a)?;
    let (cb, rb) = bounds(b)?;
    let separation = (cb - ca).norm();

    // Both sides at once: nothing can touch, so nothing is tested.
    if separation > ra + rb {
        // Still an answer, and an honest one — the caller asked how far apart
        // they are, and the bounds give a lower bound on it that is enough for
        // any caller that only wants to know whether they are in contact.
        let n = cb - ca;
        let normal = if separation > 0.0 { n.scale(1.0 / separation) } else { v3(1.0, 0.0, 0.0) };
        return Some((
            Closest {
                gap: separation - ra - rb,
                on_a: ca + normal.scale(ra),
                on_b: cb - normal.scale(rb),
                normal,
                degenerate: false,
            },
            0,
            0,
        ));
    }

    let mut best: Option<(Closest, usize, usize)> = None;
    for (i, ha) in a.iter().enumerate() {
        // One side at a time: a piece out of reach of everything on the other
        // side cannot win.
        if (ha.centre() - cb).norm() > ha.bound() + rb {
            continue;
        }
        for (j, hb) in b.iter().enumerate() {
            if (hb.centre() - ca).norm() > hb.bound() + ra {
                continue;
            }
            let Some(c) = closest(ha, hb) else { continue };
            // Strictly less, so the first piece wins a tie and the answer does
            // not depend on how the pieces happen to be ordered.
            if best.is_none_or(|(x, _, _)| c.gap < x.gap) {
                best = Some((c, i, j));
            }
        }
    }
    best
}

// ---------------------------------------------------------------------------
// what a solid presents — PLAY.md D18
// ---------------------------------------------------------------------------

/// One convex piece of a solid's boundary, with what it is made of.
///
/// **A primitive is always a filled solid. Never a shell, never hollow.**
/// `docs/PLAY.md` D18 states that first and without qualification, because the
/// rest of the decision rests on it: a hollow wooden box is not one primitive,
/// it is six solid slabs generated together, and a void is not represented at
/// all — it is simply where no primitive is.
///
/// That rule is the general statement of a failure this module already
/// measured. One convex hull over a box *encloses its own cavity*, so anything
/// inside reads as deeply interpenetrating on every frame. A shape language
/// that can express "hollow" invites exactly that mistake; one that cannot,
/// cannot.
///
/// # Material per piece, not per node
///
/// A house is stone walls, an oak door and glass panes, and a contact has to
/// read the material of the piece it actually struck. `Topology` carries one
/// material for a whole structure, which was enough while a structure was one
/// substance and is not enough for anything built.
///
/// The `substance` is carried alongside so the bake can check what D18 calls an
/// invariant rather than a convention: **a node's surface materials are a
/// partition of its solid pools.** See [`Surface::reconcile`].
#[derive(Debug, Clone, PartialEq)]
pub struct Piece {
    /// The convex solid, in the node's own frame.
    pub hull: Hull,
    /// What this piece is made of, for the contact that strikes it.
    pub material: crate::material::Material,
    /// Which substance that material was measured from, or
    /// `SubstanceId::UNSPECIATED` for a piece nobody has speciated.
    pub substance: crate::chem::SubstanceId,
    /// How much of the node's mass this piece accounts for, kg.
    ///
    /// Not used by the narrow phase, which only needs the geometry. It is what
    /// makes the partition checkable: a surface whose pieces claim more steel
    /// than the node has is describing something the node is not made of.
    pub mass: f64,
}

/// What a solid presents to the world: a union of solid convex primitives.
///
/// `docs/PLAY.md` D18. **The generator emits the pieces; nothing infers a
/// decomposition.** That is the part which makes convex decomposition — the
/// hard, unsolved half of this problem everywhere else — not arise at all.
/// Inferring convex pieces from arbitrary geometry is difficult; a generator
/// never infers, because it *knows*. A wall with a doorway emits four boxes
/// around the opening, because the recipe is what put the opening there. A tree
/// emits capsules. A rock emits hulls over its sampled bodies.
///
/// # Derived once and stored
///
/// D13 read literally: "deriving *once* and storing the result", applied to
/// shape, which is the one place the engine never applied it.
/// `Node::collision_shape` rebuilt a proxy from the member list every frame and
/// persisted nothing — a derivation with no shortcut, which is the half of the
/// third axiom that saves nothing. A `Surface` is baked when it is first asked
/// for and kept until the node's `epoch` moves, which is precisely when its
/// arrangement changed.
///
/// # What was rejected
///
/// A **baked triangle mesh**: storage grows with visual complexity, and D15's
/// whole argument is that a house is hundreds of bytes. A **CSG tree with
/// subtraction**: expressive, but it needs a second narrow phase —
/// sphere-tracing rather than GJK, iterative where GJK is exact — and the
/// solid-primitive rule removes the need, because a generator that knows where
/// the door goes emits the pieces around it instead of subtracting one.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Surface {
    pieces: Vec<Piece>,
}

/// Why a surface did not reconcile with the matter it belongs to.
///
/// D18: "a node whose mixture is all water by mass while its primitives present
/// steel would be *telling* a contact something its own bulk contradicts. That
/// is the second axiom failing by way of having two answers to one question."
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Mismatch {
    /// A piece names a substance the node has no solid pool of.
    NotHeld { substance: crate::chem::SubstanceId },
    /// The pieces of one substance claim more mass than the node holds of it.
    OverClaimed { substance: crate::chem::SubstanceId, claimed: f64, held: f64 },
}

impl Surface {
    pub fn new(pieces: Vec<Piece>) -> Surface {
        Surface { pieces }
    }

    pub fn pieces(&self) -> &[Piece] {
        &self.pieces
    }

    pub fn is_empty(&self) -> bool {
        self.pieces.is_empty()
    }

    pub fn len(&self) -> usize {
        self.pieces.len()
    }

    /// The geometry alone, which is what the narrow phase consumes.
    pub fn hulls(&self) -> impl Iterator<Item = &Hull> {
        self.pieces.iter().map(|p| &p.hull)
    }

    /// Bytes this surface costs, for the detail budget.
    pub fn bytes(&self) -> usize {
        self.pieces.iter().map(|p| std::mem::size_of::<Piece>() + p.hull.len() * std::mem::size_of::<Sphere>()).sum()
    }

    /// Check the surface against the matter it belongs to. **D18's invariant.**
    ///
    /// A node's surface materials must be a *partition of its solid pools*: a
    /// cup that is 60% ceramic may present ceramic surfaces and may not present
    /// steel ones, and the mass its ceramic primitives carry has to reconcile
    /// against its ceramic pool.
    ///
    /// This is a conservation statement about material, in the same family as
    /// `summarise(sample(m)) == m`, and D18 asks for it to be asserted rather
    /// than documented — it is cheap to check at bake time and unpleasant to
    /// retrofit once both representations exist and have drifted.
    ///
    /// A piece whose substance is `UNSPECIATED` is exempt: it is a piece nobody
    /// has said anything about, which is a *missing* answer rather than a
    /// contradicting one. `tolerance` is a fraction of the node's mass.
    pub fn reconcile(
        &self,
        mixture: &crate::chem::Mixture,
        mass: f64,
        tolerance: f64,
    ) -> Result<(), Mismatch> {
        use crate::chem::{Phase, SubstanceId};
        let mut claimed: Vec<(SubstanceId, f64)> = Vec::new();
        for p in &self.pieces {
            if p.substance == SubstanceId::UNSPECIATED {
                continue;
            }
            match claimed.iter_mut().find(|(s, _)| *s == p.substance) {
                Some((_, m)) => *m += p.mass.max(0.0),
                None => claimed.push((p.substance, p.mass.max(0.0))),
            }
        }
        let slack = tolerance * mass.abs().max(1e-30);
        for (substance, want) in claimed {
            let held = mixture.pool(substance, Phase::Solid) * mass;
            if held <= 0.0 {
                return Err(Mismatch::NotHeld { substance });
            }
            if want > held + slack {
                return Err(Mismatch::OverClaimed { substance, claimed: want, held });
            }
        }
        Ok(())
    }

    /// The same surface expressed in the parent's frame: turned by the node's
    /// orientation, then moved to where the node is.
    ///
    /// The composition is [`Hull::placed`]'s and is the order a body-fixed
    /// point is placed in. A surface is stored in the node's own frame, because
    /// one with a position baked into it would be wrong on the next frame.
    pub fn placed(&self, orientation: crate::math::Quat, offset: Vec3) -> Vec<Hull> {
        self.pieces.iter().map(|p| p.hull.placed(orientation, offset)).collect()
    }

    /// The material of the piece nearest a point, which is what a contact at
    /// that point struck.
    ///
    /// D13: "Contact then reads the material *at the point of impact*, which is
    /// also what damage and the renderer need."
    pub fn material_at(&self, point: Vec3) -> Option<&crate::material::Material> {
        let mut best: Option<(f64, &Piece)> = None;
        for p in &self.pieces {
            let d = (p.hull.centre() - point).norm() - p.hull.bound();
            if best.is_none_or(|(b, _)| d < b) {
                best = Some((d, p));
            }
        }
        best.map(|(_, p)| &p.material)
    }
}
