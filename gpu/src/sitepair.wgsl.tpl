// The site-site pair loop of `phys::flex::Flex::forces` on the card, in double
// precision: for each link (a pair of molecules inside the cut-off) the 25
// site-site terms of `SiteSite::site_pair_with` - screened Coulomb with each
// site's own charge and width, the clouds' overlap repulsion, Born-Mayer, damped
// dispersion - the forces on both molecules' sites, the energy's derivatives with
// respect to every charge and width, and the switch's force on the centres.
//
// GENERATED from `sitepair.wgsl.tpl` by `gen_induction.py`; edit the template.
//
// A link `k` joins molecules `i = lij[2k]` and `j = lij[2k+1]` with the
// centres' nearest-image separation `d` and the switch `s` and its slope `ds`.
// A site's place relative to its molecule's centre is `rel`; the separation of
// site `a` of i from site `b` of j is `d + rel_i[a] - rel_j[b]`. Per-link
// results are gathered per molecule in a fixed order (no atomics).

struct P {
    n: u32,
    ns: u32,
    nl: u32,
    ntypes: u32,
};

@group(0) @binding(0) var<uniform> p: P;
@group(0) @binding(1) var<storage, read_write> types: array<u32>;     // ns
@group(0) @binding(2) var<storage, read_write> rel: array<f64>;       // n*ns*3
@group(0) @binding(3) var<storage, read_write> qv: array<f64>;        // n*ns
@group(0) @binding(4) var<storage, read_write> sv: array<f64>;        // n*ns
@group(0) @binding(5) var<storage, read_write> pair: array<f64>;      // ntypes^2*4: A B C6 C8
@group(0) @binding(6) var<storage, read_write> ovk: array<f64>;       // ntypes^2: overlap strength
@group(0) @binding(7) var<storage, read_write> bdamp: array<f64>;     // ntypes^2: dispersion damping (<= 0: use B)
@group(0) @binding(8) var<storage, read_write> sbase: array<f64>;     // ntypes: base width, negative if none
@group(0) @binding(9) var<storage, read_write> lk: array<f64>;        // nl*6: d(3) s ds pad
@group(0) @binding(10) var<storage, read_write> lij: array<u32>;      // nl*2
@group(0) @binding(11) var<storage, read_write> inc: array<u32>;      // [n+1 starts][items k*2+side]
@group(0) @binding(12) var<storage, read_write> fa: array<f64>;       // nl*ns*3: force on i's sites
@group(0) @binding(13) var<storage, read_write> fb: array<f64>;       // nl*ns*3: force on j's sites
@group(0) @binding(14) var<storage, read_write> gqa: array<f64>;      // nl*ns
@group(0) @binding(15) var<storage, read_write> gsa: array<f64>;
@group(0) @binding(16) var<storage, read_write> gqb: array<f64>;
@group(0) @binding(17) var<storage, read_write> gsb: array<f64>;
@group(0) @binding(18) var<storage, read_write> lscal: array<f64>;    // nl*5: s*u, switch force(3), virial
@group(0) @binding(19) var<storage, read_write> out_fsite: array<f64>; // n*ns*3
@group(0) @binding(20) var<storage, read_write> out_fc: array<f64>;    // n*3
@group(0) @binding(21) var<storage, read_write> out_gq: array<f64>;    // n*ns
@group(0) @binding(22) var<storage, read_write> out_gs: array<f64>;    // n*ns

@COMMON@

// The reference's Tang-Toennies damping f_n(x) and its derivative.
fn tt(n: u32, x: f64) -> vec2<f64> {
    var sum = 1.0lf;
    var term = 1.0lf;
    for (var k = 1u; k <= n; k = k + 1u) {
        term = term * x / f64(k);
        sum = sum + term;
    }
    let e = gexp(-x);
    return vec2<f64>(1.0lf - e * sum, e * term);
}

struct Pair {
    u: f64,
    du: f64,
    g0: f64,
    g1: f64,
    g2: f64,
    g3: f64,
};

// `SiteSite::site_pair_with`.
fn site_pair_with(ta: u32, tb: u32, r: f64, qa: f64, qb: f64, sa: f64, sb: f64) -> Pair {
    let nt = p.ntypes;
    let ix = ta * nt + tb;
    let a = pair[4u * ix];
    let b = pair[4u * ix + 1u];
    let c6 = pair[4u * ix + 2u];
    let c8 = pair[4u * ix + 3u];
    let qq = qa * qb;
    let s2 = sa * sa + sb * sb;
    var coul = 0.0lf;
    var dcoul = 0.0lf;
    var g0 = 0.0lf;
    var g1 = 0.0lf;
    var g2 = 0.0lf;
    var g3 = 0.0lf;
    if (s2 > 0.0lf && qq != 0.0lf) {
        let sg = sqrt(s2);
        let x = r / sg;
        let g = gerf(x);
        let e = gauss(x);
        let c = 2.0lf / (SQRT_PI * sg) * e;
        let dsig = -qq * 2.0lf / SQRT_PI * e / (sg * sg * sg);
        g0 = qb * g / r;
        g1 = qa * g / r;
        g2 = dsig * sa;
        g3 = dsig * sb;
        coul = qq * g / r;
        dcoul = qq * (c / r - g / (r * r));
    } else {
        g0 = qb / r;
        g1 = qa / r;
        coul = qq / r;
        dcoul = -qq / (r * r);
    }
    var ov = 0.0lf;
    var dov = 0.0lf;
    let k = ovk[ix];
    let x0 = sbase[ta];
    let y0 = sbase[tb];
    if (k > 0.0lf && x0 >= 0.0lf && y0 >= 0.0lf && x0 * x0 + y0 * y0 > 0.0lf) {
        let t2 = x0 * x0 + y0 * y0;
        ov = k * gexp(-r * r / t2);
        dov = -2.0lf * r / t2 * ov;
    }
    let rep = a * gexp(-b * r);
    var bd = b;
    if (bdamp[ix] > 0.0lf) {
        bd = bdamp[ix];
    }
    let t6 = tt(6u, bd * r);
    let t8 = tt(8u, bd * r);
    let r2 = r * r;
    let r3 = r2 * r;
    let r6 = r3 * r3;
    let r8 = r6 * r2;
    var out: Pair;
    out.u = coul + rep + ov - t6.x * c6 / r6 - t8.x * c8 / r8;
    out.du = dcoul + dov - b * rep - (bd * t6.y * c6 / r6 - 6.0lf * t6.x * c6 / (r6 * r)) - (bd * t8.y * c8 / r8 - 8.0lf * t8.x * c8 / (r8 * r));
    out.g0 = g0;
    out.g1 = g1;
    out.g2 = g2;
    out.g3 = g3;
    return out;
}

@compute @workgroup_size(64)
fn k_pair_link(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= p.nl) {
        return;
    }
    let ns = p.ns;
    let i = lij[2u * k];
    let j = lij[2u * k + 1u];
    let d = vec3<f64>(lk[6u * k], lk[6u * k + 1u], lk[6u * k + 2u]);
    let s = lk[6u * k + 3u];
    let ds = lk[6u * k + 4u];
    let dist = sqrt(dot(d, d));
    var u = 0.0lf;
    var net = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    // Zero the per-site results, then accumulate over the site pairs.
    for (var a = 0u; a < ns; a = a + 1u) {
        for (var c = 0u; c < 3u; c = c + 1u) {
            fa[(k * ns + a) * 3u + c] = 0.0lf;
            fb[(k * ns + a) * 3u + c] = 0.0lf;
        }
        gqa[k * ns + a] = 0.0lf;
        gsa[k * ns + a] = 0.0lf;
        gqb[k * ns + a] = 0.0lf;
        gsb[k * ns + a] = 0.0lf;
    }
    for (var a = 0u; a < ns; a = a + 1u) {
        let ia = (i * ns + a) * 3u;
        let pa = vec3<f64>(rel[ia], rel[ia + 1u], rel[ia + 2u]);
        for (var b = 0u; b < ns; b = b + 1u) {
            let ib = (j * ns + b) * 3u;
            let pb = vec3<f64>(rel[ib], rel[ib + 1u], rel[ib + 2u]);
            let r = d + pa - pb;
            let rr = sqrt(dot(r, r));
            let t = site_pair_with(types[a], types[b], rr, qv[i * ns + a], qv[j * ns + b], sv[i * ns + a], sv[j * ns + b]);
            u = u + t.u;
            let f = r * (-t.du / rr * s);
            net = net + f;
            for (var c = 0u; c < 3u; c = c + 1u) {
                var fc = f.x;
                if (c == 1u) {
                    fc = f.y;
                } else if (c == 2u) {
                    fc = f.z;
                }
                fa[(k * ns + a) * 3u + c] = fa[(k * ns + a) * 3u + c] + fc;
                fb[(k * ns + b) * 3u + c] = fb[(k * ns + b) * 3u + c] - fc;
            }
            gqa[k * ns + a] = gqa[k * ns + a] + s * t.g0;
            gsa[k * ns + a] = gsa[k * ns + a] + s * t.g2;
            gqb[k * ns + b] = gqb[k * ns + b] + s * t.g1;
            gsb[k * ns + b] = gsb[k * ns + b] + s * t.g3;
        }
    }
    var sw = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    if (dist > 0.0lf) {
        sw = d * (-u * ds / dist);
    }
    lscal[5u * k] = s * u;
    lscal[5u * k + 1u] = sw.x;
    lscal[5u * k + 2u] = sw.y;
    lscal[5u * k + 3u] = sw.z;
    lscal[5u * k + 4u] = dot(d, net + sw);
}

var<workgroup> red3: array<vec3<f64>, 32>;

fn reduce3(v: vec3<f64>, li: u32) -> vec3<f64> {
    red3[li] = v;
    workgroupBarrier();
    for (var s = 16u; s > 0u; s = s / 2u) {
        if (li < s) {
            red3[li] = red3[li] + red3[li + s];
        }
        workgroupBarrier();
    }
    let total = red3[0];
    workgroupBarrier();
    return total;
}

// Each molecule's site: the forces, charge and width derivatives of every link it is in.
@compute @workgroup_size(32)
fn k_pair_gather(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let ns = p.ns;
    let m = wg.x / ns;
    let a = wg.x % ns;
    let li = lid.x;
    var f = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    var g = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    var c = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    let start = inc[m];
    let count = inc[m + 1u] - start;
    for (var t = li; t < count; t = t + 32u) {
        let item = inc[p.n + 1u + start + t];
        let k = item / 2u;
        let side = item % 2u;
        let o = (k * ns + a) * 3u;
        if (side == 0u) {
            f = f + vec3<f64>(fa[o], fa[o + 1u], fa[o + 2u]);
            g = g + vec3<f64>(gqa[k * ns + a], gsa[k * ns + a], 0.0lf);
            c = c + vec3<f64>(lscal[5u * k + 1u], lscal[5u * k + 2u], lscal[5u * k + 3u]);
        } else {
            f = f + vec3<f64>(fb[o], fb[o + 1u], fb[o + 2u]);
            g = g + vec3<f64>(gqb[k * ns + a], gsb[k * ns + a], 0.0lf);
            c = c - vec3<f64>(lscal[5u * k + 1u], lscal[5u * k + 2u], lscal[5u * k + 3u]);
        }
    }
    let ft = reduce3(f, li);
    let gt = reduce3(g, li);
    let ct = reduce3(c, li);
    if (li == 0u) {
        let o = (m * ns + a) * 3u;
        out_fsite[o] = ft.x;
        out_fsite[o + 1u] = ft.y;
        out_fsite[o + 2u] = ft.z;
        out_gq[m * ns + a] = gt.x;
        out_gs[m * ns + a] = gt.y;
        if (a == 0u) {
            out_fc[3u * m] = ct.x;
            out_fc[3u * m + 1u] = ct.y;
            out_fc[3u * m + 2u] = ct.z;
        }
    }
}
