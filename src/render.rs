//! A headless renderer, so that "the engine works" can be looked at rather
//! than only asserted.
//!
//! This is a *diagnostic* renderer, not a graphics pipeline. It draws members
//! as shaded tubes with a depth buffer, writes a PNG with no dependencies, and
//! runs anywhere the tests run. What it is for is answering questions the test
//! suite cannot: does the tree look like a tree, did the storm take the limbs
//! you would expect, is the damage on the side the wind came from.
//!
//! The real-time path described in `docs/GPU.md` is a different program. This
//! one deliberately trades every scrap of speed for having no dependencies and
//! no window.
//!
//! `docs/VIEWING.md` is the plan for pointing this at the test suite, so a
//! scenario can be watched unfolding rather than only asserted about. It also
//! records what is missing for that, measured — the first being that a node
//! with no topology draws nothing at all, because the loop below is bounded by
//! `topo.base.len()`.

use crate::math::{v3, Vec3};
use crate::state::Body;
use crate::topology::Topology;

/// An RGB image with a depth buffer.
pub struct Canvas {
    pub width: usize,
    pub height: usize,
    pub rgb: Vec<u8>,
    depth: Vec<f32>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Canvas {
        Canvas {
            width,
            height,
            rgb: vec![0; width * height * 3],
            depth: vec![f32::INFINITY; width * height],
        }
    }

    /// Fill with a vertical gradient — a cheap sky that also makes the
    /// silhouette readable without needing a ground plane.
    pub fn sky(&mut self, top: [u8; 3], bottom: [u8; 3]) {
        for y in 0..self.height {
            let t = y as f32 / (self.height.max(2) - 1) as f32;
            let c = [
                (top[0] as f32 + (bottom[0] as f32 - top[0] as f32) * t) as u8,
                (top[1] as f32 + (bottom[1] as f32 - top[1] as f32) * t) as u8,
                (top[2] as f32 + (bottom[2] as f32 - top[2] as f32) * t) as u8,
            ];
            for x in 0..self.width {
                let i = (y * self.width + x) * 3;
                self.rgb[i] = c[0];
                self.rgb[i + 1] = c[1];
                self.rgb[i + 2] = c[2];
            }
        }
    }

    /// Fraction of the canvas that is not still sky.
    ///
    /// Measured off the depth buffer rather than by comparing colours: a body
    /// drawn in exactly the sky's colour is still a body, and a fogged one very
    /// nearly is. Anything written has a finite depth; the sky never is.
    pub fn coverage(&self) -> f64 {
        if self.depth.is_empty() {
            return 0.0;
        }
        let drawn = self.depth.iter().filter(|d| d.is_finite()).count();
        drawn as f64 / self.depth.len() as f64
    }

    #[inline]
    fn put(&mut self, x: i64, y: i64, z: f32, c: [f32; 3]) {
        if x < 0 || y < 0 || x >= self.width as i64 || y >= self.height as i64 {
            return;
        }
        let i = y as usize * self.width + x as usize;
        if z >= self.depth[i] {
            return;
        }
        self.depth[i] = z;
        let o = i * 3;
        for k in 0..3 {
            self.rgb[o + k] = (c[k].clamp(0.0, 1.0) * 255.0 + 0.5) as u8;
        }
    }
}

/// Where the camera is and what it can see.
#[derive(Debug, Clone, Copy)]
pub struct Camera {
    pub eye: Vec3,
    pub target: Vec3,
    pub up: Vec3,
    /// Vertical field of view, radians.
    pub fov: f64,
}

impl Camera {
    /// Motion a structure: pull back far enough to see all of it, looking
    /// slightly downward from a corner so the depth reads.
    pub fn framing(centre: Vec3, radius: f64, azimuth: f64, elevation: f64) -> Camera {
        let fov: f64 = 0.6;
        let dist = radius / (fov * 0.5).tan() * 1.15;
        let dir = v3(
            azimuth.cos() * elevation.cos(),
            azimuth.sin() * elevation.cos(),
            elevation.sin(),
        );
        Camera {
            eye: centre + dir.scale(dist),
            target: centre,
            up: v3(0.0, 0.0, 1.0),
            fov,
        }
    }

    fn basis(&self) -> (Vec3, Vec3, Vec3) {
        let forward = (self.target - self.eye).unit();
        let right = forward.cross(self.up).unit();
        let up = right.cross(forward).unit();
        (right, up, forward)
    }
}

/// A point projected into screen space.
#[derive(Clone, Copy)]
struct Projected {
    x: f64,
    y: f64,
    /// Distance along the view direction. Also the depth key.
    z: f64,
    /// Pixels per metre at this depth, for sizing member radii.
    ppm: f64,
}

/// How to colour what is drawn.
#[derive(Debug, Clone, Copy)]
pub struct Style {
    pub sky_top: [u8; 3],
    pub sky_bottom: [u8; 3],
    /// Base colour of intact structure.
    pub member: [f32; 3],
    /// Colour of detached or destroyed material.
    pub broken: [f32; 3],
    /// Loose unstructured matter — litter, rubble, debris.
    pub litter: [f32; 3],
    pub light: Vec3,
    /// Distance over which the colour fades towards the sky, metres. Depth cue.
    pub fog: f64,
    /// Draw parts with no structural role.
    pub show_litter: bool,
}

impl Default for Style {
    fn default() -> Self {
        Style {
            sky_top: [26, 32, 44],
            sky_bottom: [58, 66, 82],
            member: [0.62, 0.47, 0.32],
            broken: [0.78, 0.29, 0.22],
            litter: [0.38, 0.40, 0.36],
            light: v3(-0.5, -0.7, 0.6),
            fog: 60.0,
            show_litter: false,
        }
    }
}

impl Style {
    /// Warm daylight on a pale sky.
    pub fn daylight() -> Style {
        Style {
            sky_top: [176, 196, 222],
            sky_bottom: [232, 236, 232],
            member: [0.40, 0.30, 0.20],
            broken: [0.72, 0.22, 0.16],
            litter: [0.55, 0.53, 0.45],
            ..Default::default()
        }
    }

    /// Everything scorched.
    pub fn burned() -> Style {
        Style {
            sky_top: [58, 30, 20],
            sky_bottom: [128, 66, 34],
            member: [0.16, 0.13, 0.12],
            broken: [0.85, 0.42, 0.12],
            litter: [0.22, 0.18, 0.16],
            ..Default::default()
        }
    }
}

/// Draw a structure: its members as shaded tubes, optionally its loose matter.
///
/// `intact` marks which parts are still attached; anything else is drawn in the
/// broken colour, which is what makes a damage result legible at a glance.
///
/// The `Some` case of [`draw`], kept under its own name and signature because
/// `examples/damage.rs` and the renderer's own test are written against it.
pub fn draw_structure(
    canvas: &mut Canvas,
    camera: &Camera,
    bodies: &[Body],
    topo: &Topology,
    intact: &[bool],
    style: &Style,
) {
    draw(canvas, camera, bodies, Some(topo), intact, &Paint::Role, style);
}

/// How to colour what is drawn.
///
/// `docs/VIEWING.md`: **a diagnostic image should be a measurement, not a
/// legend.** `Role` is the original behaviour and is the one variant that is a
/// label the caller supplies; every other variant is read off the body being
/// drawn, and the range it was read against is returned by [`draw`] so a colour
/// can be turned back into a number.
///
/// `Measured` is the escape hatch, and it is deliberately not a closure: a
/// quantity like a phase fraction cannot be read off a `Body` at all, because a
/// body carries a substance id and the registry that resolves it lives in
/// `World`. So the module that *can* measure it does, and hands the numbers
/// over. That keeps this file a rasteriser — it knows `math`, `state` and
/// `topology`, and a rasteriser that could reach a `World` would eventually
/// solve something.
#[derive(Debug, Clone, PartialEq)]
pub enum Paint {
    /// Member, broken, litter — what the caller says each part is.
    Role,
    /// Kelvin.
    Temperature { range: Option<(f64, f64)> },
    /// Metres per second, in the frame the bodies are given in.
    Speed { range: Option<(f64, f64)> },
    /// A quantity measured elsewhere, one value per body.
    Measured { values: Vec<f64>, range: Option<(f64, f64)> },
}

impl Paint {
    /// The value this paint reads off one body, or `None` for a label.
    fn value(&self, i: usize, b: &Body) -> Option<f64> {
        match self {
            Paint::Role => None,
            Paint::Temperature { .. } => Some(b.temperature),
            Paint::Speed { .. } => Some(b.vel.norm()),
            Paint::Measured { values, .. } => values.get(i).copied(),
        }
    }

    fn stated_range(&self) -> Option<(f64, f64)> {
        match self {
            Paint::Role => None,
            Paint::Temperature { range } | Paint::Speed { range } => *range,
            Paint::Measured { range, .. } => *range,
        }
    }
}

/// What a frame turned out to be, so it can be read back as numbers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Drawn {
    /// Fraction of the canvas the subject covers. The renderer's own test
    /// asserts on this; a film that is all sky is a camera pointed nowhere.
    pub covered: f64,
    /// The range the colour ramp was read against, where there was one.
    pub range: Option<(f64, f64)>,
    pub parts: usize,
}

/// Blue through white to orange: cold low, hot high, and legible in print and
/// to the common forms of colour blindness because the lightness runs with the
/// value rather than only the hue.
fn ramp(t: f64) -> [f32; 3] {
    let t = t.clamp(0.0, 1.0) as f32;
    if t < 0.5 {
        let u = t * 2.0;
        [0.15 + 0.8 * u, 0.35 + 0.6 * u, 0.75 + 0.25 * u]
    } else {
        let u = (t - 0.5) * 2.0;
        [0.95, 0.95 - 0.55 * u, 1.0 - 0.9 * u]
    }
}

/// Draw whatever a node is holding: members where there are members, and a disc
/// at its own radius for everything else.
///
/// `docs/VIEWING.md`'s first missing piece, measured there: `draw_structure`
/// iterated `bodies.len().min(topo.base.len())`, so a node whose contents are
/// loose produced an empty sky — `show_litter` included, because the litter
/// branch was inside that same loop. Everything Phase 1 built is in that class:
/// a ball in a box, a node's hydro parcels, a fragment in flight, and a ball of
/// matter that has not become a planet yet.
pub fn draw(
    canvas: &mut Canvas,
    camera: &Camera,
    bodies: &[Body],
    topo: Option<&Topology>,
    intact: &[bool],
    paint: &Paint,
    style: &Style,
) -> Drawn {
    canvas.sky(style.sky_top, style.sky_bottom);
    let (right, up, forward) = camera.basis();
    let half_h = (camera.fov * 0.5).tan();
    let aspect = canvas.width as f64 / canvas.height as f64;
    let light = style.light.unit();

    let (cw, ch) = (canvas.width as f64, canvas.height as f64);
    let project = |p: Vec3| -> Option<Projected> {
        let rel = p - camera.eye;
        let z = rel.dot(forward);
        if z <= 1e-6 {
            return None;
        }
        let sx = rel.dot(right) / (z * half_h * aspect);
        let sy = rel.dot(up) / (z * half_h);
        Some(Projected {
            x: (sx * 0.5 + 0.5) * cw,
            y: (0.5 - sy * 0.5) * ch,
            z,
            ppm: ch * 0.5 / (z * half_h),
        })
    };

    // The range the ramp is read against: whatever the caller stated, else the
    // frame's own extremes, which is what makes an unattended film legible
    // without anybody having guessed the scale in advance.
    let range = match paint.stated_range() {
        Some(r) => Some(r),
        None => {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for (i, b) in bodies.iter().enumerate() {
                if let Some(v) = paint.value(i, b) {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }
            if lo.is_finite() && hi > lo { Some((lo, hi)) } else { None }
        }
    };

    // Painter's order is handled by the depth buffer, but drawing far members
    // first still reduces the number of overwritten pixels.
    let mut order: Vec<usize> = (0..bodies.len()).collect();
    order.sort_by(|a, b| {
        let za = (bodies[*a].pos - camera.eye).dot(forward);
        let zb = (bodies[*b].pos - camera.eye).dot(forward);
        zb.partial_cmp(&za).unwrap_or(std::cmp::Ordering::Equal)
    });

    // What a part is coloured: the measurement where there is one, and the
    // caller's label where there is not.
    let colour_of = |i: usize, b: &Body, structural: bool| -> [f32; 3] {
        match (paint.value(i, b), range) {
            (Some(v), Some((lo, hi))) if hi > lo => ramp((v - lo) / (hi - lo)),
            _ if !structural => style.litter,
            _ if intact.get(i).copied().unwrap_or(true) => style.member,
            _ => style.broken,
        }
    };

    let mut parts = 0usize;
    for &i in &order {
        let b = match bodies.get(i) {
            Some(b) => b,
            None => continue,
        };
        let structural = topo.is_some_and(|t| {
            i < t.joints.len() && i < t.base.len() && (t.tip[i] - t.base[i]).norm2() > 0.0
        });
        if !structural {
            // A loose body is a disc at its own radius. Litter inside a
            // structure is still optional — it is what `show_litter` was for —
            // but a node with no topology at all is *all* loose, and refusing
            // to draw it is the defect this branch exists to fix.
            if topo.is_some() && !style.show_litter {
                continue;
            }
            if let Some(p) = project(b.pos) {
                let r = (b.radius * p.ppm).max(0.6);
                disc(canvas, p, r, colour_of(i, b, false), style.fog);
                parts += 1;
            }
            continue;
        }
        let t = match topo {
            Some(t) => t,
            None => continue,
        };
        let (a, bb) = (t.base[i], t.tip[i]);
        let (pa, pb) = match (project(a), project(bb)) {
            (Some(x), Some(y)) => (x, y),
            _ => continue,
        };
        let radius_px = (t.joints[i].radius * pa.ppm).max(0.55);
        tube(canvas, pa, pb, radius_px, (bb - a).unit(), light, colour_of(i, b, true), style.fog);
        parts += 1;
    }

    Drawn { covered: canvas.coverage(), range, parts }
}

/// A shaded cylinder between two projected points.
///
/// Rasterised by distance-to-segment over the bounding box rather than by
/// stepping along the axis and drawing a perpendicular run of pixels at each
/// step. The stepping approach is the obvious one and it leaves diagonal
/// members visibly hatched, because consecutive perpendicular runs overlap on
/// near-vertical lines and separate on diagonal ones.
fn tube(
    canvas: &mut Canvas,
    a: Projected,
    b: Projected,
    radius: f64,
    axis: Vec3,
    light: Vec3,
    colour: [f32; 3],
    fog: f64,
) {
    let r = radius.max(0.55);
    let (x0, y0, x1, y1) = (a.x, a.y, b.x, b.y);
    let lo_x = (x0.min(x1) - r - 1.0).floor().max(0.0) as i64;
    let hi_x = (x0.max(x1) + r + 1.0).ceil().min(canvas.width as f64) as i64;
    let lo_y = (y0.min(y1) - r - 1.0).floor().max(0.0) as i64;
    let hi_y = (y0.max(y1) + r + 1.0).ceil().min(canvas.height as f64) as i64;
    if hi_x <= lo_x || hi_y <= lo_y {
        return;
    }
    let (dx, dy) = (x1 - x0, y1 - y0);
    let len2 = (dx * dx + dy * dy).max(1e-12);

    // A cylinder lit from a direction is brightest where its surface normal
    // faces the light, and that normal sweeps across the tube — which is what
    // gives a drawn branch its roundness instead of looking like a flat stick.
    let axis_light = 1.0 - axis.dot(light).abs();

    for py in lo_y..hi_y {
        for px in lo_x..hi_x {
            let (fx, fy) = (px as f64 + 0.5, py as f64 + 0.5);
            // Parameter of the closest point on the segment.
            let t = (((fx - x0) * dx + (fy - y0) * dy) / len2).clamp(0.0, 1.0);
            let (cx, cy) = (x0 + dx * t, y0 + dy * t);
            let dist = ((fx - cx).powi(2) + (fy - cy).powi(2)).sqrt();
            if dist > r {
                continue;
            }
            let across = dist / r;
            let bulge = (1.0 - across * across).max(0.0).sqrt();
            let z = a.z + (b.z - a.z) * t;
            let shade = (0.28 + 0.72 * bulge * axis_light) as f32;
            let depth = (z - bulge * r / a.ppm.max(1e-9)) as f32;
            let f = (1.0 - (z / fog).min(1.0) * 0.55) as f32;
            let c = [
                colour[0] * shade * f + 0.10 * (1.0 - f),
                colour[1] * shade * f + 0.12 * (1.0 - f),
                colour[2] * shade * f + 0.14 * (1.0 - f),
            ];
            canvas.put(px, py, depth, c);
        }
    }
}

fn disc(canvas: &mut Canvas, p: Projected, radius: f64, colour: [f32; 3], fog: f64) {
    let r = radius.max(0.5);
    let span = r.ceil() as i64;
    for dy in -span..=span {
        for dx in -span..=span {
            let d2 = (dx * dx + dy * dy) as f64 / (r * r);
            if d2 > 1.0 {
                continue;
            }
            let shade = (0.4 + 0.6 * (1.0 - d2).sqrt()) as f32;
            let f = (1.0 - (p.z / fog).min(1.0) * 0.55) as f32;
            let c = [
                colour[0] * shade * f + 0.10 * (1.0 - f),
                colour[1] * shade * f + 0.12 * (1.0 - f),
                colour[2] * shade * f + 0.14 * (1.0 - f),
            ];
            canvas.put(p.x as i64 + dx, p.y as i64 + dy, p.z as f32, c);
        }
    }
}

// ---------------------------------------------------------------------------
// PNG output, with no dependencies
// ---------------------------------------------------------------------------

/// Write the canvas as a PNG.
///
/// The zlib stream uses stored (uncompressed) deflate blocks. That makes the
/// files larger than they need to be and the encoder about forty lines instead
/// of a dependency — the right trade for a diagnostic tool that has to run in
/// any environment the tests do.
pub fn write_png(canvas: &Canvas, path: &str) -> std::io::Result<()> {
    use std::io::Write;
    let mut out = Vec::with_capacity(canvas.rgb.len() + 4096);
    out.extend_from_slice(&[137, 80, 78, 71, 13, 10, 26, 10]);

    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&(canvas.width as u32).to_be_bytes());
    ihdr.extend_from_slice(&(canvas.height as u32).to_be_bytes());
    ihdr.extend_from_slice(&[8, 2, 0, 0, 0]); // 8-bit, truecolour RGB
    chunk(&mut out, b"IHDR", &ihdr);

    // Raw scanlines, each prefixed with filter type 0.
    let mut raw = Vec::with_capacity((canvas.width * 3 + 1) * canvas.height);
    for y in 0..canvas.height {
        raw.push(0);
        let o = y * canvas.width * 3;
        raw.extend_from_slice(&canvas.rgb[o..o + canvas.width * 3]);
    }

    let mut z = vec![0x78, 0x01];
    let mut i = 0;
    while i < raw.len() {
        let n = (raw.len() - i).min(65535);
        let last = if i + n >= raw.len() { 1u8 } else { 0u8 };
        z.push(last);
        z.extend_from_slice(&(n as u16).to_le_bytes());
        z.extend_from_slice(&(!(n as u16)).to_le_bytes());
        z.extend_from_slice(&raw[i..i + n]);
        i += n;
    }
    z.extend_from_slice(&adler32(&raw).to_be_bytes());
    chunk(&mut out, b"IDAT", &z);
    chunk(&mut out, b"IEND", &[]);

    let mut f = std::fs::File::create(path)?;
    f.write_all(&out)
}

fn chunk(out: &mut Vec<u8>, tag: &[u8; 4], data: &[u8]) {
    out.extend_from_slice(&(data.len() as u32).to_be_bytes());
    let start = out.len();
    out.extend_from_slice(tag);
    out.extend_from_slice(data);
    let crc = crc32(&out[start..]);
    out.extend_from_slice(&crc.to_be_bytes());
}

fn crc32(data: &[u8]) -> u32 {
    let mut table = [0u32; 256];
    for (i, t) in table.iter_mut().enumerate() {
        let mut c = i as u32;
        for _ in 0..8 {
            c = if c & 1 != 0 { 0xEDB8_8320 ^ (c >> 1) } else { c >> 1 };
        }
        *t = c;
    }
    let mut c = 0xFFFF_FFFFu32;
    for &b in data {
        c = table[((c ^ b as u32) & 0xFF) as usize] ^ (c >> 8);
    }
    c ^ 0xFFFF_FFFF
}

fn adler32(data: &[u8]) -> u32 {
    let (mut a, mut b) = (1u32, 0u32);
    for &byte in data {
        a = (a + byte as u32) % 65521;
        b = (b + a) % 65521;
    }
    (b << 16) | a
}

// ---------------------------------------------------------------------------
// a sequence, with one camera
// ---------------------------------------------------------------------------

/// A camera, a style, a paint and a size, made once.
///
/// `docs/VIEWING.md`: *one camera for every frame in a set.* `examples/damage.rs`
/// paid for that lesson — auto-framing each shot rescales the subject and hides
/// exactly what the comparison exists to show, and a tree that had lost a third
/// of its height looked identical to one that had not. This is the lesson
/// turned into a type, so the camera is fixed by construction rather than by
/// remembering.
#[derive(Debug, Clone)]
pub struct Shot {
    pub camera: Camera,
    pub style: Style,
    pub paint: Paint,
    pub width: usize,
    pub height: usize,
}

impl Shot {
    /// Frame something of this size, from a corner so the depth reads.
    ///
    /// The fog is set from the framing radius rather than left at its 60 m
    /// default, which is `VIEWING.md`'s third measured gap: at six metres
    /// everything was unfogged and at a kilometre everything was fog. Two and a
    /// half radii puts the far side of the subject at about half fog whatever
    /// the subject is, from a nucleus to a galaxy.
    pub fn framing(centre: Vec3, radius: f64, azimuth: f64, elevation: f64) -> Shot {
        let camera = Camera::framing(centre, radius, azimuth, elevation);
        Shot {
            camera,
            style: Style { fog: (2.5 * radius).max(1e-30), ..Style::daylight() },
            paint: Paint::Role,
            width: 480,
            height: 360,
            }
    }

    pub fn painted(mut self, paint: Paint) -> Shot {
        self.paint = paint;
        self
    }

    pub fn sized(mut self, width: usize, height: usize) -> Shot {
        self.width = width;
        self.height = height;
        self
    }

    pub fn styled(mut self, style: Style) -> Shot {
        let fog = self.style.fog;
        self.style = Style { fog, ..style };
        self
    }

    pub fn canvas(&self) -> Canvas {
        Canvas::new(self.width, self.height)
    }
}

/// A numbered sequence of PNGs in a directory, or nothing at all.
///
/// **Off unless asked.** `PHYS_FILM=<dir>` or [`Film::open`] returns `None` and
/// every call on it is a no-op, so `cargo test` stays silent, fast and
/// identical. With the variable set, every instrumented test drops a filmstrip
/// in that directory.
///
/// Not an assertion. `VIEWING.md` is blunt about why: an image assertion is the
/// archetype of a test that cannot fail, and this repo's rule is that a test is
/// verified against the defect it catches. These are diagnostics; the
/// assertions stay numeric.
pub struct Film {
    dir: std::path::PathBuf,
    name: String,
    frame: usize,
}

impl Film {
    /// Open a filmstrip called `name`, if `PHYS_FILM` names a directory.
    pub fn open(name: &str) -> Option<Film> {
        let dir = std::env::var_os("PHYS_FILM")?;
        let dir = std::path::PathBuf::from(dir);
        std::fs::create_dir_all(&dir).ok()?;
        Some(Film { dir, name: name.to_string(), frame: 0 })
    }

    /// Whether anything is being recorded, so a caller can skip the gathering.
    pub fn recording() -> bool {
        std::env::var_os("PHYS_FILM").is_some()
    }

    pub fn frames(&self) -> usize {
        self.frame
    }

    /// Write the next frame. The number is the film's, not the world's.
    pub fn shoot(&mut self, canvas: &Canvas) {
        let path = self.dir.join(format!("{}-{:04}.png", self.name, self.frame));
        if write_png(canvas, &path.to_string_lossy()).is_ok() {
            self.frame += 1;
        }
    }
}
