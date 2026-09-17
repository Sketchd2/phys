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
