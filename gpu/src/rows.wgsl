// The non-local correlation's double sum, one row per invocation, in single
// precision: phys::electrons::vdw::CpuRows on the GPU (PLAY.md E7). Each
// dispatch adds the columns [j0, j1) to every row's running sums, with Kahan
// compensation carried across dispatches in `acc`.

struct Params {
    n: u32,
    j0: u32,
    j1: u32,
    m: u32,
    d_max: f32,
    scale: f32,
    asym_from: f32,
    asym_c: f32,
    with_c: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
// x, y, z, w n
@group(0) @binding(1) var<storage, read> pts: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> qs: array<f32>;
@group(0) @binding(3) var<storage, read> table: array<f32>;
// Per row: sums a, b, cx, cy, cz, then their compensations.
@group(0) @binding(4) var<storage, read_write> acc: array<f32>;

const TILE: u32 = 256u;
var<workgroup> tile_p: array<vec4<f32>, 256>;
var<workgroup> tile_q: array<f32, 256>;

fn v(p: i32, q: i32) -> f32 {
    return table[u32(p) * params.m + u32(q)];
}

// The grid continued one node past either edge by linear extrapolation, as
// KernelTable::bicubic does.
fn at(a: i32, b: i32) -> f32 {
    let m = i32(params.m);
    var a0 = a; var a1 = a; var fa = false;
    if (a < 0) { a0 = 0; a1 = 1; fa = true; } else if (a >= m) { a0 = m - 1; a1 = m - 2; fa = true; }
    var b0 = b; var b1 = b; var fb = false;
    if (b < 0) { b0 = 0; b1 = 1; fb = true; } else if (b >= m) { b0 = m - 1; b1 = m - 2; fb = true; }
    let base = v(a0, b0);
    var da = 0.0;
    if (fa) { da = base - v(a1, b0); }
    var db = 0.0;
    if (fb) { db = base - v(a0, b1); }
    return base + da + db;
}

fn cr(t: f32) -> vec4<f32> {
    let t2 = t * t;
    let t3 = t2 * t;
    return vec4<f32>(-0.5 * t3 + t2 - 0.5 * t, 1.5 * t3 - 2.5 * t2 + 1.0, -1.5 * t3 + 2.0 * t2 + 0.5 * t, 0.5 * t3 - 0.5 * t2);
}

fn dcr(t: f32) -> vec4<f32> {
    let t2 = t * t;
    return vec4<f32>(-1.5 * t2 + 2.0 * t - 0.5, 4.5 * t2 - 5.0 * t, -4.5 * t2 + 4.0 * t + 0.5, 1.5 * t2 - t);
}

// phi(d1, d2) and its two partial derivatives: KernelTable::phi_and_slopes.
fn phi3(d1: f32, d2: f32) -> vec3<f32> {
    if (min(d1, d2) >= params.asym_from) {
        let s = d1 * d1 + d2 * d2;
        let f = -params.asym_c / (d1 * d1 * d2 * d2 * s);
        return vec3<f32>(f, f * (-2.0 / d1 - 2.0 * d1 / s), f * (-2.0 / d2 - 2.0 * d2 / s));
    }
    let c1 = min(d1, params.d_max);
    let c2 = min(d2, params.d_max);
    let top = f32(params.m - 1u);
    let x = clamp(c1 / (1.0 + c1) * params.scale, 0.0, top);
    let y = clamp(c2 / (1.0 + c2) * params.scale, 0.0, top);
    let last = f32(params.m - 2u);
    let i = i32(min(floor(x), last));
    let j = i32(min(floor(y), last));
    let tx = x - f32(i);
    let ty = y - f32(j);
    let wx = cr(tx);
    let wy = cr(ty);
    let dwx = dcr(tx);
    let dwy = dcr(ty);
    var s = 0.0;
    var sx = 0.0;
    var sy = 0.0;
    for (var p = 0; p < 4; p = p + 1) {
        for (var q = 0; q < 4; q = q + 1) {
            let val = at(i - 1 + p, j - 1 + q);
            s = s + wx[p] * wy[q] * val;
            sx = sx + dwx[p] * wy[q] * val;
            sy = sy + wx[p] * dwy[q] * val;
        }
    }
    var gx = 0.0;
    if (d1 < params.d_max) { gx = params.scale / ((1.0 + c1) * (1.0 + c1)); }
    var gy = 0.0;
    if (d2 < params.d_max) { gy = params.scale / ((1.0 + c2) * (1.0 + c2)); }
    return vec3<f32>(s, sx * gx, sy * gy);
}

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let i = gid.x;
    let live = i < params.n;
    var me = vec4<f32>(0.0);
    var qi = 0.0;
    if (live) { me = pts[i]; qi = qs[i]; }
    var s: array<f32, 5>;
    var c: array<f32, 5>;
    for (var k = 0u; k < 5u; k = k + 1u) {
        if (live) { s[k] = acc[i * 10u + k]; c[k] = acc[i * 10u + 5u + k]; } else { s[k] = 0.0; c[k] = 0.0; }
    }
    var start = params.j0;
    loop {
        if (start >= params.j1) { break; }
        let j = start + lid.x;
        if (j < params.j1) { tile_p[lid.x] = pts[j]; tile_q[lid.x] = qs[j]; }
        workgroupBarrier();
        let count = min(TILE, params.j1 - start);
        if (live) {
            for (var t = 0u; t < count; t = t + 1u) {
                let o = tile_p[t];
                let d = me.xyz - o.xyz;
                let r = sqrt(dot(d, d));
                let f = phi3(qi * r, tile_q[t] * r);
                var add: array<f32, 5>;
                add[0] = o.w * f.x;
                add[1] = o.w * r * f.y;
                add[2] = 0.0; add[3] = 0.0; add[4] = 0.0;
                if (params.with_c != 0u && r > 0.0) {
                    let radial = o.w * (qi * f.y + tile_q[t] * f.z) / r;
                    add[2] = radial * d.x; add[3] = radial * d.y; add[4] = radial * d.z;
                }
                for (var k = 0u; k < 5u; k = k + 1u) {
                    // Kahan.
                    let yk = add[k] - c[k];
                    let tk = s[k] + yk;
                    c[k] = (tk - s[k]) - yk;
                    s[k] = tk;
                }
            }
        }
        workgroupBarrier();
        start = start + TILE;
    }
    if (live) {
        for (var k = 0u; k < 5u; k = k + 1u) { acc[i * 10u + k] = s[k]; acc[i * 10u + 5u + k] = c[k]; }
    }
}
