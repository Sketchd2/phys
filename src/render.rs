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
    /// Frame a structure: pull back far enough to see all of it, looking
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
pub fn draw_structure(
    canvas: &mut Canvas,
    camera: &Camera,
    bodies: &[Body],
    topo: &Topology,
    intact: &[bool],
    style: &Style,
) {
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

    // Painter's order is handled by the depth buffer, but drawing far members
    // first still reduces the number of overwritten pixels.
    let n = bodies.len().min(topo.base.len());
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|a, b| {
        let za = (bodies[*a].pos - camera.eye).dot(forward);
        let zb = (bodies[*b].pos - camera.eye).dot(forward);
        zb.partial_cmp(&za).unwrap_or(std::cmp::Ordering::Equal)
    });

    for &i in &order {
        let structural = i < topo.bonds.len() && (topo.tip[i] - topo.base[i]).norm2() > 0.0;
        if !structural {
            if !style.show_litter {
                continue;
            }
            if let Some(p) = project(bodies[i].pos) {
                let r = (bodies[i].radius * p.ppm).max(0.6);
                disc(canvas, p, r, style.litter, style.fog);
            }
            continue;
        }
        let (a, b) = (topo.base[i], topo.tip[i]);
        let (pa, pb) = match (project(a), project(b)) {
            (Some(x), Some(y)) => (x, y),
            _ => continue,
        };
        let radius_px = (topo.bonds[i].radius * pa.ppm).max(0.55);
        let colour = if intact.get(i).copied().unwrap_or(true) {
            style.member
        } else {
            style.broken
        };
        tube(canvas, pa, pb, radius_px, (b - a).unit(), light, colour, style.fog);
    }
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
