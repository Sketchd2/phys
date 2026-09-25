//! A planet's ocean at the scale of the planet: its dynamic tide, and the sea
//! the wind raises on it.
//!
//! # Why a planet needs this at all
//!
//! `docs/PLAY.md` Phase 5's done-when asks for a tide: "a squiggle drawn in
//! sand is gone by the next tide". A tide is the ocean's response to the
//! differential gravity of the bodies outside the planet, and it is a property
//! of the whole ocean — a point on a shore does not know when high water comes
//! without the sea that brings it. The owner's decision for this phase is a
//! **dynamic** tide rather than an equilibrium one: the ocean responds through
//! the shallow-water equations, with its own resonances, rather than standing
//! at `-potential / g` everywhere.
//!
//! # What the ocean is made of
//!
//! Nothing stored about it is a choice. The water is the liquid in the node's
//! own mixture; its depth is that water's volume spread over the node's
//! surface; `g` is the node's own surface gravity; the Coriolis term is its
//! own spin; the forcing is the exact tidal potential of every other body its
//! parent holds, at their actual positions; and the seabed's friction is the
//! log-law against the ground's own grain.
//!
//! **The seabed is flat for now, and that is the owner's call rather than
//! this module's.** Ground on a sphere carries no relief yet: it is to come
//! from the planet's root `Program`, generated offline, and until then a scene
//! that needs a shore states one by hand. So the ocean here is an aqua-planet:
//! one depth everywhere. Relief would enter as a depth per cell and nothing
//! else in the scheme changes.
//!
//! # The scheme
//!
//! Linear shallow water on a cubed sphere of `n x n` cells a face, laid out
//! on the same faces the ground uses (`recipe::face_axes`):
//!
//! ```text
//!   du/dt   = -g grad(eta - eta_eq) - f k x u - r u
//!   deta/dt = -div(H u)
//! ```
//!
//! - **Water is conserved exactly.** Each edge's flux is computed once and
//!   given to its two cells with opposite signs.
//! - **Coriolis cannot pump energy.** It is applied as the exact rotation it
//!   is, by `f dt` about the local vertical.
//! - **Friction cannot overshoot.** It is implicit.
//! - **The pressure gradient and the divergence are adjoint**, so the discrete
//!   energy is conserved by everything but friction: see `Ocean::step`.
//! - Forward-backward in time: the velocity steps with the old surface, and
//!   the surface with the new velocity, which is stable to a Courant number of
//!   about one against `sqrt(g H)`.

use crate::math::Vec3;
use crate::units::G;

/// One cell of the ocean.
#[derive(Debug, Clone, Copy)]
pub struct Cell {
    /// Unit vector from the planet's centre, in the planet's own axes.
    pub up: Vec3,
    /// Area, m^2.
    pub area: f64,
    /// Surface elevation above the rest level, m.
    pub eta: f64,
    /// Depth-averaged velocity, tangent to the sphere, m/s.
    pub velocity: Vec3,
    /// Energy of the wind-raised sea over this cell, J/m^2. See `Ocean::raise`.
    pub waves: f64,
    /// Which way that sea is travelling: a unit vector tangent to the sphere,
    /// or zero where there is none.
    pub heading: Vec3,
}

/// An edge between two cells, stored once.
#[derive(Debug, Clone, Copy)]
struct Edge {
    a: u32,
    b: u32,
    /// Unit vector from `a` towards `b`, along the chord.
    toward: Vec3,
    /// Length of the edge, m.
    length: f64,
}

/// A planet's ocean.
#[derive(Debug, Clone)]
pub struct Ocean {
    /// Cells a face is divided into, along a side.
    pub n: usize,
    /// Radius of the sea surface, m.
    pub radius: f64,
    /// Depth at rest, m. One number while the seabed is flat.
    pub depth: f64,
    /// Surface gravity, m/s^2.
    pub g: f64,
    /// Seabed drag coefficient, from the log-law against the ground's grain.
    pub drag: f64,
    /// Density of the water, kg/m^3: its liquid's rest density.
    pub density: f64,
    /// Surface tension of the water, N/m, which sets the slowest wave there is.
    pub tension: f64,
    pub cells: Vec<Cell>,
    edges: Vec<Edge>,
    /// Time the ocean has been stepped to, s, on the planet's clock.
    pub time: f64,
    /// World time owed to the ocean and not yet stepped, s: less than one
    /// stable step, carried from frame to frame so that nothing is lost to
    /// rounding a frame to whole steps.
    pub owed: f64,
}

/// Von Kármán's constant, the log-law's one number: the universal slope of a
/// turbulent boundary layer's velocity against the logarithm of height. The
/// owner's decision for this phase is that the wind's stress is the log-law's,
/// with this stated as a universal in the same family as random loose packing.
pub const VON_KARMAN: f64 = 0.41;

/// Nikuradse's result: a surface of roughness elements of size `k` has a
/// roughness length of `k / 30`. The log-law's other half, used for the seabed
/// against its grain and for the sea surface against its own waves.
pub const ROUGHNESS_PER_HEIGHT: f64 = 1.0 / 30.0;

/// The log-law's drag coefficient for a flow of depth or height `z` over a
/// surface of roughness elements of size `k`: `(kappa / ln(z / z0))^2` with
/// `z0 = k / 30`.
pub fn log_law_drag(z: f64, k: f64) -> f64 {
    let z0 = (k * ROUGHNESS_PER_HEIGHT).max(1e-12);
    let l = (z.max(z0 * std::f64::consts::E) / z0).ln();
    (VON_KARMAN / l).powi(2)
}

/// Michell's limiting steepness: a wave of height `H` and length `L` breaks
/// past `H / L = 1/7`. The owner's decision for Phase 5 is that a wind-raised
/// sea stands at this steepness while it grows, so its height sets its length
/// and its length its speed — see [`phase_speed`].
pub const BREAKING_STEEPNESS: f64 = 1.0 / 7.0;

/// McCowan's depth-limited breaking: a wave breaks once its height passes
/// `0.78` of the depth it is in, which is where a sea meets a shore.
pub const MCCOWAN: f64 = 0.78;

/// Significant height of a sea holding `energy` J/m^2: `4 sqrt(E / rho g)`,
/// which is four standard deviations of the surface.
pub fn significant_height(energy: f64, density: f64, g: f64) -> f64 {
    if !(energy > 0.0) || !(density > 0.0) || !(g > 0.0) {
        return 0.0;
    }
    4.0 * (energy / (density * g)).sqrt()
}

/// The energy of a sea of significant height `h`, J/m^2 — the inverse of
/// [`significant_height`].
pub fn energy_of_height(h: f64, density: f64, g: f64) -> f64 {
    density * g * h * h / 16.0
}

/// The slowest a wave on this water can go, m/s: the minimum of the
/// gravity–capillary dispersion relation, `(4 g sigma / rho)^(1/4)`. Shorter
/// waves are held by surface tension and go faster, longer ones by gravity and
/// go faster too; a wind slower than this raises no sea at all.
pub fn slowest_wave(g: f64, tension: f64, density: f64) -> f64 {
    if !(g > 0.0) || !(tension > 0.0) || !(density > 0.0) {
        return 0.0;
    }
    (4.0 * g * tension / density).powf(0.25)
}

/// Phase and group speed of a sea of significant height `h` standing at the
/// limiting steepness, in water `depth` deep, m/s.
///
/// Its length is `h / BREAKING_STEEPNESS` and its speed the gravity-wave
/// dispersion relation's at that length, `sqrt(g/k tanh(k d))`, with the group
/// speed `c/2 (1 + 2kd / sinh 2kd)` — half the phase speed in deep water, all
/// of it in shallow. Neither is ever slower than [`slowest_wave`], which is
/// what a sea is made of before it has any height at all.
pub fn phase_speed(h: f64, depth: f64, g: f64, tension: f64, density: f64) -> (f64, f64) {
    let floor = slowest_wave(g, tension, density);
    if !(h > 0.0) || !(g > 0.0) {
        return (floor, floor);
    }
    let k = std::f64::consts::TAU * BREAKING_STEEPNESS / h;
    let kd = k * depth.max(0.0);
    let c = (g / k * kd.tanh()).sqrt();
    if !(c > floor) {
        return (floor, floor);
    }
    let group = if kd > 20.0 { 0.5 * c } else { 0.5 * c * (1.0 + 2.0 * kd / (2.0 * kd).sinh()) };
    (c, group)
}

/// The height a wind of speed `u` raises a deep sea to: where the phase speed
/// at the limiting steepness reaches the wind's, `2 pi s u^2 / g`. The owner's
/// limit for Phase 5 — the sea grows until it outruns the wind. At 10 m/s it
/// is 9.2 m, about four times what is observed of a fully developed sea,
/// because nothing here takes energy out of a growing sea but breaking.
pub fn saturated_height(u: f64, g: f64) -> f64 {
    if !(g > 0.0) {
        return 0.0;
    }
    std::f64::consts::TAU * BREAKING_STEEPNESS * u * u / g
}

/// The height at which a sea in water `depth` deep goes as fast as a wind of
/// speed `u`, m — [`saturated_height`] with the depth's dispersion in it, by
/// bisection on [`phase_speed`]. Zero for a wind no wave is slower than.
pub fn outrun_height(u: f64, depth: f64, g: f64, tension: f64, density: f64) -> f64 {
    if !(u > slowest_wave(g, tension, density)) {
        return 0.0;
    }
    // Shallow water caps the phase speed at sqrt(g d): a wind faster than
    // that is never outrun, and only breaking stops the sea.
    if u * u >= g * depth.max(0.0) {
        return f64::INFINITY;
    }
    let (mut lo, mut hi) = (0.0, saturated_height(u, g).max(1e-6));
    while phase_speed(hi, depth, g, tension, density).0 < u {
        hi *= 2.0;
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if phase_speed(mid, depth, g, tension, density).0 < u {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// A wind over one cell of the ocean, in the planet's own axes.
#[derive(Debug, Clone, Copy)]
pub struct Wind {
    /// The cell it blows over.
    pub cell: usize,
    /// The air's velocity relative to the water's surface, m/s, tangent.
    pub velocity: Vec3,
    /// Height the velocity is measured at, m.
    pub height: f64,
    /// Density of the air, kg/m^3.
    pub air_density: f64,
    /// Area of the cell it covers, m^2.
    pub area: f64,
}

/// What a wind did to the ocean in one step.
#[derive(Debug, Clone, Copy, Default)]
pub struct Blown {
    /// Momentum the air gave the water, kg m/s, in the planet's axes.
    pub momentum: Vec3,
    /// Energy the sea took up, J.
    pub to_waves: f64,
    /// Energy the currents took up, J.
    pub to_currents: f64,
    /// Energy the sea gave up by breaking, J — heat.
    pub broken: f64,
}

/// The exact tidal potential at `x` of a mass `m` at `r`, both from the
/// planet's centre: the potential less its value and its gradient at the
/// centre, which is what moves the whole planet rather than its ocean.
pub fn tidal_potential(x: Vec3, r: Vec3, m: f64) -> f64 {
    let d = r.norm();
    if !(d > 0.0) || !(m > 0.0) {
        return 0.0;
    }
    let s = (r - x).norm().max(1e-300);
    -G * m * (1.0 / s - 1.0 / d - r.dot(x) / (d * d * d))
}

fn face_point(face: u8, alpha: f64, beta: f64) -> Vec3 {
    let (right, up, out) = crate::recipe::face_axes(face);
    (out + right.scale(alpha.tan()) + up.scale(beta.tan())).unit()
}

impl Ocean {
    /// An ocean of `depth` on a sphere of `radius`, `n` cells a face side.
    pub fn new(n: usize, radius: f64, depth: f64, g: f64, drag: f64, density: f64, tension: f64) -> Ocean {
        let n = n.max(2);
        let q = std::f64::consts::FRAC_PI_4;
        let step = 2.0 * q / n as f64;
        let index = |f: usize, i: usize, j: usize| (f * n + j) * n + i;
        let mut cells = Vec::with_capacity(6 * n * n);
        for f in 0..6u8 {
            for j in 0..n {
                for i in 0..n {
                    let (a0, b0) = (-q + i as f64 * step, -q + j as f64 * step);
                    let up = face_point(f, a0 + 0.5 * step, b0 + 0.5 * step);
                    // Area from the four corners, as two spherical triangles'
                    // flat approximations: exact enough at any n that resolves
                    // anything, and the areas sum to the sphere's within it.
                    let c = [
                        face_point(f, a0, b0),
                        face_point(f, a0 + step, b0),
                        face_point(f, a0 + step, b0 + step),
                        face_point(f, a0, b0 + step),
                    ];
                    let tri = |p: Vec3, q: Vec3, r: Vec3| 0.5 * (q - p).cross(r - p).norm();
                    let area = (tri(c[0], c[1], c[2]) + tri(c[0], c[2], c[3])) * radius * radius;
                    cells.push(Cell { up, area, eta: 0.0, velocity: Vec3::ZERO, waves: 0.0, heading: Vec3::ZERO });
                }
            }
        }
        // Scale the areas so they sum to the sphere's exactly: a flat facet
        // under-counts a curved one by the same small fraction everywhere, and
        // conserving water needs the total right.
        let total: f64 = cells.iter().map(|c| c.area).sum();
        let sphere = 4.0 * std::f64::consts::PI * radius * radius;
        for c in cells.iter_mut() {
            c.area *= sphere / total;
        }
        // Edges: each cell to its right and upper neighbour, found across a
        // face boundary by stepping past it and asking which cell the point is
        // in. Every edge is then found exactly once from one side.
        let locate = |d: Vec3| -> usize {
            let (mut best, mut bi) = (f64::NEG_INFINITY, 0usize);
            for f in 0..6u8 {
                let (_, _, out) = crate::recipe::face_axes(f);
                let o = d.dot(out);
                if o > best {
                    best = o;
                    bi = f as usize;
                }
            }
            let (right, up, out) = crate::recipe::face_axes(bi as u8);
            let o = d.dot(out);
            let a = (d.dot(right) / o).atan();
            let b = (d.dot(up) / o).atan();
            let i = (((a + q) / step).floor() as i64).clamp(0, n as i64 - 1) as usize;
            let j = (((b + q) / step).floor() as i64).clamp(0, n as i64 - 1) as usize;
            index(bi, i, j)
        };
        let mut edges = Vec::with_capacity(12 * n * n);
        let mut seen = std::collections::HashSet::new();
        for f in 0..6u8 {
            for j in 0..n {
                for i in 0..n {
                    let a = index(f as usize, i, j);
                    let (ac, bc) = (-q + (i as f64 + 0.5) * step, -q + (j as f64 + 0.5) * step);
                    for (da, db) in [(step, 0.0), (-step, 0.0), (0.0, step), (0.0, -step)] {
                        let b = locate(face_point(f, ac + da, bc + db));
                        if b == a {
                            continue;
                        }
                        let key = (a.min(b), a.max(b));
                        if !seen.insert(key) {
                            continue;
                        }
                        let (pa, pb) = (cells[a].up, cells[b].up);
                        let toward = (pb - pa).unit();
                        // The shared side: the mean of the two cells' sides.
                        let length = 0.5 * (cells[a].area.sqrt() + cells[b].area.sqrt());
                        edges.push(Edge { a: a as u32, b: b as u32, toward, length });
                    }
                }
            }
        }
        Ocean { n, radius, depth, g, drag, density, tension, cells, edges, time: 0.0, owed: 0.0 }
    }

    /// The largest step the scheme is stable at, s.
    pub fn stable_step(&self) -> f64 {
        let c = (self.g * self.depth).sqrt().max(1e-9);
        let side = self.cells.iter().map(|c| c.area.sqrt()).fold(f64::INFINITY, f64::min);
        0.5 * side / c
    }

    /// Water held above the rest level, m^3. Zero to rounding for all time,
    /// which is the conservation the scheme is built for.
    pub fn excess_volume(&self) -> f64 {
        self.cells.iter().map(|c| c.eta * c.area).sum()
    }

    /// Which cell a direction is over.
    pub fn cell_of(&self, dir: Vec3) -> usize {
        let d = dir.unit();
        let mut best = (f64::NEG_INFINITY, 0usize);
        for (i, c) in self.cells.iter().enumerate() {
            let s = c.up.dot(d);
            if s > best.0 {
                best = (s, i);
            }
        }
        best.1
    }

    /// Step by `dt` against the equilibrium elevation `eq` of each cell (the
    /// tidal potential over `-g`) with the planet spinning at `spin` rad/s,
    /// in the planet's own axes.
    pub fn step(&mut self, dt: f64, eq: &[f64], spin: Vec3) {
        if !(dt > 0.0) {
            return;
        }
        let nc = self.cells.len();
        // Momentum: the pressure gradient over each cell's edges, as the
        // *difference* of the heads across each edge, half to each side.
        //
        // That is the form whose adjoint is the continuity step below, and so
        // the one that conserves the discrete energy. The mean-of-heads form —
        // Gauss's theorem read literally — differs from it by each cell's head
        // times the sum of its own edge normals, which on a sphere is not
        // zero; measured, it grew a slow mode from a tide's 0.8 m to 239 m in
        // a week.
        let mut push = vec![Vec3::ZERO; nc];
        for e in &self.edges {
            let (a, b) = (e.a as usize, e.b as usize);
            let head = |k: usize| self.cells[k].eta - eq.get(k).copied().unwrap_or(0.0);
            let f = e.toward.scale(0.5 * (head(b) - head(a)) * e.length);
            push[a] += f;
            push[b] += f;
        }
        for (k, c) in self.cells.iter_mut().enumerate() {
            let up = c.up;
            // Gradient of the head, projected onto the surface.
            let grad = push[k].scale(1.0 / c.area);
            let grad = grad - up.scale(grad.dot(up));
            let mut u = c.velocity - grad.scale(self.g * dt);
            // Coriolis: an exact rotation by f dt about the local vertical.
            let f = 2.0 * spin.dot(up);
            let angle = -f * dt;
            let (s, co) = angle.sin_cos();
            u = u.scale(co) + up.cross(u).scale(s) + up.scale(up.dot(u) * (1.0 - co));
            // The seabed, implicitly: quadratic drag, linearised about the
            // current speed.
            let r = self.drag * u.norm() / self.depth.max(1e-9);
            u = u.scale(1.0 / (1.0 + r * dt));
            // And kept on the surface.
            c.velocity = u - up.scale(u.dot(up));
        }
        // Continuity: each edge's flux once, both signs.
        let mut dv = vec![0.0; nc];
        for e in &self.edges {
            let (a, b) = (e.a as usize, e.b as usize);
            let u = (self.cells[a].velocity + self.cells[b].velocity).scale(0.5);
            let flux = self.depth * u.dot(e.toward) * e.length * dt;
            dv[a] -= flux;
            dv[b] += flux;
        }
        for (k, c) in self.cells.iter_mut().enumerate() {
            c.eta += dv[k] / c.area;
        }
        self.time += dt;
    }

    /// The significant height of the sea over a cell, m.
    pub fn wave_height(&self, cell: usize) -> f64 {
        self.cells.get(cell).map(|c| significant_height(c.waves, self.density, self.g)).unwrap_or(0.0)
    }

    /// Let each wind work on the water under it for `dt`.
    ///
    /// **The stress is the log-law's**, `rho_air C_d U^2` with `C_d` from
    /// [`log_law_drag`] at the wind's height over roughness elements the size
    /// of the sea's own height — the owner's decision for Phase 5, which
    /// measures the roughness instead of reading Charnock's constant. The
    /// momentum it carries goes into the water column as current. **The sea
    /// takes energy at `tau c`** while it is slower than the wind, and stops
    /// taking it when its phase speed reaches the wind's, which is the owner's
    /// limit; what the air loses beyond both is turbulence, and the caller
    /// books it as heat. A sea past McCowan's depth limit breaks, and what it
    /// gives up is reported.
    pub fn raise(&mut self, dt: f64, winds: &[Wind]) -> Vec<Blown> {
        let mut out = Vec::with_capacity(winds.len());
        for w in winds {
            let mut b = Blown::default();
            let u = w.velocity.norm();
            if !(dt > 0.0) || !(u > 0.0) || !(w.area > 0.0) || w.cell >= self.cells.len() {
                out.push(b);
                continue;
            }
            let along = w.velocity.scale(1.0 / u);
            let (density, g, depth, tension) = (self.density, self.g, self.depth, self.tension);
            let c = &mut self.cells[w.cell];
            let h = significant_height(c.waves, density, g);
            let drag = log_law_drag(w.height.max(1e-3), h);
            let stress = w.air_density * drag * u * u;
            let impulse = stress * w.area * dt;
            b.momentum = along.scale(impulse);
            // The column the stress lands on, and the energy its current took.
            let column = density * depth.max(1e-9) * c.area;
            let before = 0.5 * column * c.velocity.norm2();
            let v = c.velocity + b.momentum.scale(1.0 / column);
            c.velocity = v - c.up.scale(v.dot(c.up));
            b.to_currents = 0.5 * column * c.velocity.norm2() - before;
            let (speed, _) = phase_speed(h, depth, g, tension, density);
            if speed < u {
                // Never past the height at which it would outrun the wind: a
                // step that would take it there takes it there and no further.
                let most = energy_of_height(outrun_height(u, depth, g, tension, density), density, g);
                let gained = (stress * speed * w.area * dt).min(((most - c.waves) * c.area).max(0.0));
                let total = c.waves * c.area + gained;
                // The sea heads the way the energy in it came from.
                let heading = c.heading.scale(c.waves * c.area) + along.scale(gained);
                c.heading = if heading.norm() > 0.0 { heading.unit() } else { Vec3::ZERO };
                c.waves = total / c.area;
                b.to_waves = gained;
            }
            let cap = energy_of_height(MCCOWAN * depth, density, g);
            if c.waves > cap {
                b.broken = (c.waves - cap) * c.area;
                b.to_waves -= b.broken;
                c.waves = cap;
            }
            out.push(b);
        }
        out
    }

    /// Carry the sea across the ocean for `dt`, at its group speed, the way it
    /// is heading: each edge passes the energy upwind of it to the cell
    /// downwind, once, so nothing is made or lost in transit. A sea keeps
    /// going after the wind that raised it drops, which is what swell is.
    pub fn propagate(&mut self, dt: f64) {
        if !(dt > 0.0) {
            return;
        }
        let nc = self.cells.len();
        let group: Vec<f64> = (0..nc)
            .map(|k| {
                let h = significant_height(self.cells[k].waves, self.density, self.g);
                if h > 0.0 { phase_speed(h, self.depth, self.g, self.tension, self.density).1 } else { 0.0 }
            })
            .collect();
        let mut energy = vec![0.0; nc];
        let mut heading = vec![Vec3::ZERO; nc];
        for (k, c) in self.cells.iter().enumerate() {
            energy[k] = c.waves * c.area;
            heading[k] = c.heading.scale(energy[k]);
        }
        for e in &self.edges {
            let (a, b) = (e.a as usize, e.b as usize);
            for (from, to, sign) in [(a, b, 1.0), (b, a, -1.0)] {
                let c = &self.cells[from];
                let out = c.heading.dot(e.toward.scale(sign));
                if !(out > 0.0) || !(c.waves > 0.0) {
                    continue;
                }
                // Never more than the cell holds in a step, whatever the edge.
                let flux = (c.waves * group[from] * out * e.length * dt).min(c.waves * c.area);
                energy[from] -= flux;
                energy[to] += flux;
                heading[to] += c.heading.scale(flux);
            }
        }
        for (k, c) in self.cells.iter_mut().enumerate() {
            c.waves = (energy[k] / c.area).max(0.0);
            let h = heading[k] - c.up.scale(heading[k].dot(c.up));
            c.heading = if c.waves > 0.0 && h.norm() > 0.0 { h.unit() } else { Vec3::ZERO };
        }
    }

    /// The step the sea can be carried at: a group speed never crosses more
    /// than half a cell in one.
    pub fn wave_step(&self) -> f64 {
        let side = self.cells.iter().map(|c| c.area.sqrt()).fold(f64::INFINITY, f64::min);
        let fastest = (0..self.cells.len())
            .map(|k| {
                let h = self.wave_height(k);
                if h > 0.0 { phase_speed(h, self.depth, self.g, self.tension, self.density).1 } else { 0.0 }
            })
            .fold(0.0f64, f64::max);
        if fastest > 0.0 { 0.5 * side / fastest } else { f64::INFINITY }
    }
}
