// The reciprocal-space sums of `phys::ewald` on the card, in double precision.
//
// One workgroup of 64 threads a wave vector `k` adds over the sites (the structure
// factors, the arm sums), and one a site adds over the wave vectors (the forces and
// the field): each thread takes every 64th term in index order and the 64 partial
// sums are added by a fixed tree, so the result is the same bits every time, and
// differs from the reference's single running sum by round-off. (A thread a wave
// vector looping over 648 sites left the card idle: it took as long as the CPU.)
// The reference's scalar sums over `k` (the energy, the virial's loop) are left to the
// CPU, which gets the per-`k` arrays back, so they are the reference's own loops. The phases are the
// reference's too: the CPU builds each site's table of e^(i m theta) along each axis (cos and
// sin, which WGSL has no f64 version of) and this builds
// `e^(i k.r) = (tx * ty) * tz` from them, as `Ewald::phases` does.
//
// Table layout, charges and dipoles alike: axis 0's block (`ns * n0` complex,
// site-major), then axis 1's (`ns * n1`), then axis 2's; `n_a = nmax_a + 1`.

struct P {
    nc: u32,
    nd: u32,
    nk: u32,
    n0: u32,
    n1: u32,
    n2: u32,
    pad0: u32,
    pad1: u32,
};

@group(0) @binding(0) var<uniform> p: P;
@group(0) @binding(1) var<storage, read_write> tab_c: array<f64>;
@group(0) @binding(2) var<storage, read_write> tab_d: array<f64>;
@group(0) @binding(3) var<storage, read_write> kidx: array<i32>;   // nk*3
@group(0) @binding(4) var<storage, read_write> kv: array<f64>;     // nk*3
@group(0) @binding(5) var<storage, read_write> fk: array<f64>;     // nk
@group(0) @binding(6) var<storage, read_write> q: array<f64>;      // nc
@group(0) @binding(7) var<storage, read_write> carm: array<f64>;   // nc*3
@group(0) @binding(8) var<storage, read_write> mu: array<f64>;     // nd*3
@group(0) @binding(9) var<storage, read_write> darm: array<f64>;   // nd*3
@group(0) @binding(10) var<storage, read_write> sq: array<f64>;    // nk*2
@group(0) @binding(11) var<storage, read_write> sd: array<f64>;    // nk*2
@group(0) @binding(12) var<storage, read_write> stot: array<f64>;  // nk*2
@group(0) @binding(13) var<storage, read_write> tcv: array<f64>;   // nk*2
@group(0) @binding(14) var<storage, read_write> tdv: array<f64>;   // nk*2
@group(0) @binding(15) var<storage, read_write> o_d: array<f64>;   // nd*3: the field or the response
@group(0) @binding(16) var<storage, read_write> o_gch: array<f64>; // nc*3
@group(0) @binding(17) var<storage, read_write> o_gdp: array<f64>; // nd*3
@group(0) @binding(18) var<storage, read_write> sc: array<f64>;    // [8 pi / V, a dipole's own coupling]

fn cmul(a: vec2<f64>, b: vec2<f64>) -> vec2<f64> {
    return vec2<f64>(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

// e^(i m theta_axis) for a site, the negative m by conjugate.
fn tab_get(is_charge: bool, site: u32, axis: u32, m: i32) -> vec2<f64> {
    let ns = select(p.nd, p.nc, is_charge);
    var n = p.n0;
    var off = 0u;
    if (axis == 1u) {
        n = p.n1;
        off = ns * p.n0 * 2u;
    } else if (axis == 2u) {
        n = p.n2;
        off = ns * (p.n0 + p.n1) * 2u;
    }
    let mm = u32(abs(m));
    let base = off + (site * n + mm) * 2u;
    var re: f64;
    var im: f64;
    if (is_charge) {
        re = tab_c[base];
        im = tab_c[base + 1u];
    } else {
        re = tab_d[base];
        im = tab_d[base + 1u];
    }
    if (m < 0) {
        im = -im;
    }
    return vec2<f64>(re, im);
}

// e^(i k.r) for a site and wave vector j: (tx * ty) * tz.
fn phase(is_charge: bool, site: u32, j: u32) -> vec2<f64> {
    let tx = tab_get(is_charge, site, 0u, kidx[3u * j]);
    let ty = tab_get(is_charge, site, 1u, kidx[3u * j + 1u]);
    let tz = tab_get(is_charge, site, 2u, kidx[3u * j + 2u]);
    return cmul(cmul(tx, ty), tz);
}

fn kdot(j: u32, x: f64, y: f64, z: f64) -> f64 {
    return kv[3u * j] * x + kv[3u * j + 1u] * y + kv[3u * j + 2u] * z;
}

var<workgroup> r0: array<f64, 64>;
var<workgroup> r1: array<f64, 64>;
var<workgroup> r2: array<f64, 64>;

// Add the 64 partial sums in r0, r1, r2 by a fixed tree; the totals end in [0].
fn reduce3(t: u32) {
    for (var s = 32u; s > 0u; s = s >> 1u) {
        workgroupBarrier();
        if (t < s) {
            r0[t] = r0[t] + r0[t + s];
            r1[t] = r1[t] + r1[t + s];
            r2[t] = r2[t] + r2[t + s];
        }
    }
    workgroupBarrier();
}

// S_q(k) = sum q e^(ik.r): one workgroup a wave vector.
@compute @workgroup_size(64)
fn k_sq(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let j = wg.x;
    var sr = 0.0lf;
    var si = 0.0lf;
    for (var i = t; i < p.nc; i = i + 64u) {
        let ph = phase(true, i, j);
        let qi = q[i];
        sr = sr + qi * ph.x;
        si = si + qi * ph.y;
    }
    r0[t] = sr;
    r1[t] = si;
    r2[t] = 0.0lf;
    reduce3(t);
    if (t == 0u) {
        sq[2u * j] = r0[0];
        sq[2u * j + 1u] = r1[0];
    }
}

// D(k) = sum i (k.mu) e^(ik.r).
@compute @workgroup_size(64)
fn k_dsum(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let j = wg.x;
    var sr = 0.0lf;
    var si = 0.0lf;
    for (var a = t; a < p.nd; a = a + 64u) {
        let ph = phase(false, a, j);
        let km = kdot(j, mu[3u * a], mu[3u * a + 1u], mu[3u * a + 2u]);
        sr = sr + (-km * ph.y);
        si = si + km * ph.x;
    }
    r0[t] = sr;
    r1[t] = si;
    r2[t] = 0.0lf;
    reduce3(t);
    if (t == 0u) {
        sd[2u * j] = r0[0];
        sd[2u * j + 1u] = r1[0];
    }
}

// The total structure factor, S_q + D (one thread a wave vector).
@compute @workgroup_size(64)
fn k_stot(@builtin(global_invocation_id) g: vec3<u32>) {
    let j = g.x;
    if (j >= p.nk) {
        return;
    }
    if (p.nd > 0u) {
        stot[2u * j] = sq[2u * j] + sd[2u * j];
        stot[2u * j + 1u] = sq[2u * j + 1u] + sd[2u * j + 1u];
    } else {
        stot[2u * j] = sq[2u * j];
        stot[2u * j + 1u] = sq[2u * j + 1u];
    }
}

fn sval(which: u32, j: u32) -> vec2<f64> {
    if (which == 0u) {
        return vec2<f64>(sq[2u * j], sq[2u * j + 1u]);
    }
    if (which == 1u) {
        return vec2<f64>(sd[2u * j], sd[2u * j + 1u]);
    }
    return vec2<f64>(stot[2u * j], stot[2u * j + 1u]);
}

// Ewald::grad_dipole for one site: pref * sum_j fk Re[S* i e^(ik.r)] k; the sum is
// spread over the workgroup, and every thread returns the total.
fn grad_dip(is_charge: bool, site: u32, which: u32, t: u32) -> vec3<f64> {
    var gx = 0.0lf;
    var gy = 0.0lf;
    var gz = 0.0lf;
    for (var j = t; j < p.nk; j = j + 64u) {
        let s = sval(which, j);
        let ph = phase(is_charge, site, j);
        let re = s.x * (-ph.y) + s.y * ph.x;
        let w = fk[j] * re;
        gx = gx + w * kv[3u * j];
        gy = gy + w * kv[3u * j + 1u];
        gz = gz + w * kv[3u * j + 2u];
    }
    r0[t] = gx;
    r1[t] = gy;
    r2[t] = gz;
    reduce3(t);
    let pref = sc[0];
    return vec3<f64>(r0[0] * pref, r1[0] * pref, r2[0] * pref);
}

// The long-range field of the permanent charges at the dipoles: minus grad_dipole(S_q).
@compute @workgroup_size(64)
fn k_field(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let a = wg.x;
    let v = grad_dip(false, a, 0u, t);
    if (t == 0u) {
        o_d[3u * a] = -v.x;
        o_d[3u * a + 1u] = -v.y;
        o_d[3u * a + 2u] = -v.z;
    }
}

// The dipole-dipole coupling less each dipole's own: grad_dipole(D) - mu * own.
@compute @workgroup_size(64)
fn k_resp(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let a = wg.x;
    let v = grad_dip(false, a, 1u, t);
    if (t == 0u) {
        let own = sc[1];
        o_d[3u * a] = v.x - mu[3u * a] * own;
        o_d[3u * a + 1u] = v.y - mu[3u * a + 1u] * own;
        o_d[3u * a + 2u] = v.z - mu[3u * a + 2u] * own;
    }
}

// dE/dr of the charges: grad_dipole(S_q + D) scaled by each charge.
@compute @workgroup_size(64)
fn k_gch(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let i = wg.x;
    let v = grad_dip(true, i, 2u, t);
    if (t == 0u) {
        let qi = q[i];
        o_gch[3u * i] = v.x * qi;
        o_gch[3u * i + 1u] = v.y * qi;
        o_gch[3u * i + 2u] = v.z * qi;
    }
}

// dE/dr of the dipoles: pref * sum_j fk Re[-S* (k.mu) e^(ik.r)] k.
@compute @workgroup_size(64)
fn k_gdpos(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let a = wg.x;
    let m = vec3<f64>(mu[3u * a], mu[3u * a + 1u], mu[3u * a + 2u]);
    var gx = 0.0lf;
    var gy = 0.0lf;
    var gz = 0.0lf;
    for (var j = t; j < p.nk; j = j + 64u) {
        let km = kdot(j, m.x, m.y, m.z);
        let s = vec2<f64>(stot[2u * j], stot[2u * j + 1u]);
        let ph = phase(false, a, j);
        let re = -(s.x * ph.x + s.y * ph.y);
        let w = fk[j] * re * km;
        gx = gx + w * kv[3u * j];
        gy = gy + w * kv[3u * j + 1u];
        gz = gz + w * kv[3u * j + 2u];
    }
    r0[t] = gx;
    r1[t] = gy;
    r2[t] = gz;
    reduce3(t);
    if (t == 0u) {
        let pref = sc[0];
        o_gdp[3u * a] = r0[0] * pref;
        o_gdp[3u * a + 1u] = r1[0] * pref;
        o_gdp[3u * a + 2u] = r2[0] * pref;
    }
}

// sum_i q_i (k.a_i) e^(ik.r_i).
@compute @workgroup_size(64)
fn k_arm_c(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let j = wg.x;
    var sr = 0.0lf;
    var si = 0.0lf;
    for (var i = t; i < p.nc; i = i + 64u) {
        let ph = phase(true, i, j);
        let c = q[i] * kdot(j, carm[3u * i], carm[3u * i + 1u], carm[3u * i + 2u]);
        sr = sr + c * ph.x;
        si = si + c * ph.y;
    }
    r0[t] = sr;
    r1[t] = si;
    r2[t] = 0.0lf;
    reduce3(t);
    if (t == 0u) {
        tcv[2u * j] = r0[0];
        tcv[2u * j + 1u] = r1[0];
    }
}

// sum_a (k.mu_a)(k.b_a) e^(ik.r_a).
@compute @workgroup_size(64)
fn k_arm_d(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let j = wg.x;
    var sr = 0.0lf;
    var si = 0.0lf;
    for (var a = t; a < p.nd; a = a + 64u) {
        let ph = phase(false, a, j);
        let c = kdot(j, mu[3u * a], mu[3u * a + 1u], mu[3u * a + 2u]) * kdot(j, darm[3u * a], darm[3u * a + 1u], darm[3u * a + 2u]);
        sr = sr + c * ph.x;
        si = si + c * ph.y;
    }
    r0[t] = sr;
    r1[t] = si;
    r2[t] = 0.0lf;
    reduce3(t);
    if (t == 0u) {
        tdv[2u * j] = r0[0];
        tdv[2u * j + 1u] = r1[0];
    }
}
