// A whole Ewald force evaluation's long-range work on the card, in double precision:
// the induced-dipole solve of `phys::induction::solve_ewald` (the short-range tensors are
// the smeared tensor less the one at the Ewald width `sp`, the molecule's own charge-dipole
// terms taken off, the permanent charges' reciprocal field added to the right-hand side and the
// reciprocal dipole-dipole response added to every matrix-vector product) AND the reciprocal
// sums of `phys::ewald`, so that the whole conjugate-gradient loop, reciprocal response
// included, runs on the card in one submission and only a few small arrays come back.
// GENERATED from `ewald_step.wgsl.tpl` (itself derived from `induction.wgsl.tpl`, which is the
// validated flexible-water solve and is untouched) by `gen_induction.py`; edit the template.
//
// Group 0 is the induction's, group 1 the reciprocal sums' (`ewald.wgsl`'s kernels with their
// names changed). The reciprocal dipoles are the induction's polarisable sites, same order.
//
// The induced-dipole solve of `phys::induction::solve`, in double precision, on
// the card. Every kernel is a gather: a thread or workgroup owns one entry, one
// site or one link and adds what it reads in a fixed order, so the answer is
// the same every run (no atomics, which have no f64 anyway). WGSL has no f64
// `exp`, so `gexp` below builds one; `gerf` is the CPU's own series.
//
// GENERATED: this file is `induction.wgsl.tpl` with its tables filled in by
// `gen_induction.py` (the series' reciprocals and the exponential's
// coefficients are computed there, not typed). Edit the template.
//
// Layout: molecules are identical, `np` polarisable sites and `nch` charges
// each. Site `s = m * np + a`; charge `c = m * nch + ci`. A link `k` joins
// molecules `i = lij[2k]` and `j = lij[2k+1]` with the shift added to j's
// positions. Its entries:
//   dd[k][a][b]   (a on i, b on j)                  index  k*np*np + a*np + b
//   cd1[k][a][c]  (dipole a on i, charge c on j)    index  k*2*np*nch + a*nch + c
//   cd2[k][b][c]  (dipole b on j, charge c on i)    index  k*2*np*nch + np*nch + b*nch + c
// The molecule's own pairs of dipoles are `intra[m][pair]`, a < b, in order.
//
// The kernels that gather over a site's links (e0, the product, the forces)
// give the site a workgroup of 32, each thread taking links t, t+32, ... and
// the partial sums added by a fixed tree.

struct P {
    nm: u32,
    np: u32,
    nch: u32,
    nl: u32,
    n: u32,
    nc: u32,
    npairs: u32,
    max_iterations: u32,
};

@group(0) @binding(0) var<uniform> p: P;
@group(0) @binding(1) var<storage, read_write> pol_pos: array<f64>;   // n*3
@group(0) @binding(2) var<storage, read_write> pol_sig: array<f64>;   // n
@group(0) @binding(3) var<storage, read_write> inv_alpha: array<f64>; // n
@group(0) @binding(4) var<storage, read_write> ch_pos: array<f64>;    // nc*3
@group(0) @binding(5) var<storage, read_write> ch_q: array<f64>;      // nc
@group(0) @binding(6) var<storage, read_write> ch_sig: array<f64>;    // nc
@group(0) @binding(7) var<storage, read_write> lf: array<f64>;        // nl*8: shift(3) w dw d(3)
@group(0) @binding(8) var<storage, read_write> lij: array<u32>;       // nl*2
@group(0) @binding(9) var<storage, read_write> inc: array<u32>;       // [nm+1 starts][items k*2+side]
@group(0) @binding(10) var<storage, read_write> dd_t: array<f64>;     // entries*4: ta tbp da dbp
@group(0) @binding(11) var<storage, read_write> dd_r: array<f64>;     // entries*3
@group(0) @binding(12) var<storage, read_write> cd_h: array<f64>;     // entries*2: h dh
@group(0) @binding(13) var<storage, read_write> cd_r: array<f64>;     // entries*3
@group(0) @binding(14) var<storage, read_write> in_t: array<f64>;     // nm*npairs*4
@group(0) @binding(15) var<storage, read_write> in_r: array<f64>;     // nm*npairs*3
@group(0) @binding(16) var<storage, read_write> e0: array<f64>;       // n*3
@group(0) @binding(17) var<storage, read_write> mu: array<f64>;       // n*3
@group(0) @binding(18) var<storage, read_write> rv: array<f64>;       // n*3
@group(0) @binding(19) var<storage, read_write> pv: array<f64>;       // n*3
@group(0) @binding(20) var<storage, read_write> apv: array<f64>;      // n*3
@group(0) @binding(21) var<storage, read_write> zv: array<f64>;       // n*3
@group(0) @binding(22) var<storage, read_write> sc: array<f64>;       // 0 tol 1 bnorm 2 rz 3 done 4 iterations 5 rn
@group(0) @binding(23) var<storage, read_write> dd_f: array<f64>;     // entries*3
@group(0) @binding(24) var<storage, read_write> cd_f: array<f64>;     // entries*3
@group(0) @binding(25) var<storage, read_write> fpol: array<f64>;     // n*3
@group(0) @binding(26) var<storage, read_write> fch: array<f64>;      // nc*3
@group(0) @binding(27) var<storage, read_write> lforce: array<f64>;   // nl*3
@group(0) @binding(28) var<storage, read_write> lvir: array<f64>;     // nl
@group(0) @binding(29) var<storage, read_write> icd_h: array<f64>;    // nm*np*nch*2: h dh (the molecule's own charge-dipole terms, negated)
@group(0) @binding(30) var<storage, read_write> icd_r: array<f64>;    // nm*np*nch*3
@group(0) @binding(31) var<storage, read_write> e0x: array<f64>;      // n*3: the permanent charges' reciprocal field at the dipoles

// ---- the reciprocal sums (group 1) -------------------------------------------

struct RP {
    nc: u32,
    nk: u32,
    n0: u32,
    n1: u32,
    n2: u32,
    pad0: u32,
    pad1: u32,
    pad2: u32,
};

@group(1) @binding(0) var<uniform> rp: RP;
@group(1) @binding(1) var<storage, read_write> tab_c: array<f64>;
@group(1) @binding(2) var<storage, read_write> tab_d: array<f64>;
@group(1) @binding(3) var<storage, read_write> kidx: array<i32>;   // nk*3
@group(1) @binding(4) var<storage, read_write> kv: array<f64>;     // nk*3
@group(1) @binding(5) var<storage, read_write> fk: array<f64>;     // nk
@group(1) @binding(6) var<storage, read_write> rq: array<f64>;     // nc (the merged reciprocal charges)
@group(1) @binding(7) var<storage, read_write> carm: array<f64>;   // nc*3
@group(1) @binding(8) var<storage, read_write> darm: array<f64>;   // n*3
@group(1) @binding(9) var<storage, read_write> sq: array<f64>;     // nk*2
@group(1) @binding(10) var<storage, read_write> sd: array<f64>;    // nk*2
@group(1) @binding(11) var<storage, read_write> stot: array<f64>;  // nk*2
@group(1) @binding(12) var<storage, read_write> tcv: array<f64>;   // nk*2
@group(1) @binding(13) var<storage, read_write> tdv: array<f64>;   // nk*2
@group(1) @binding(14) var<storage, read_write> o_gch: array<f64>; // nc*3
@group(1) @binding(15) var<storage, read_write> o_gdp: array<f64>; // n*3

@COMMON@

// (ta, tbp, da, dbp) of the reference's `dd_terms`.
fn dd_terms(r: f64, s: f64) -> vec4<f64> {
    let x = r / s;
    let g = gerf(x);
    let c = 2.0lf / (SQRT_PI * s) * gauss(x);
    let ir = 1.0lf / r;
    let ir2 = ir * ir;
    let ir3 = ir2 * ir;
    let is2 = 1.0lf / (s * s);
    let a = g * ir3 - c * ir2;
    let bp = 3.0lf * g * ir3 * ir2 - c * (3.0lf * ir2 * ir2 + 2.0lf * is2 * ir2);
    let da = 3.0lf * c * ir3 - 3.0lf * g * ir2 * ir2 + 2.0lf * c * is2 * ir;
    let dbp = 15.0lf * c * ir3 * ir2 - 15.0lf * g * ir3 * ir3 + 10.0lf * c * is2 * ir3 + 4.0lf * c * is2 * is2 * ir;
    return vec4<f64>(a, bp, da, dbp);
}

// The smeared tensor less the one at the Ewald width (sc[6]).
fn dd_sr(r: f64, s: f64) -> vec4<f64> {
    return dd_terms(r, s) - dd_terms(r, sc[6]);
}

// (h, dh) of `charge_terms`.
fn charge_terms(r: f64, s_a: f64) -> vec2<f64> {
    let x = r / s_a;
    let ex = gauss(x);
    let gg = gerf(x) - 2.0lf / SQRT_PI * x * ex;
    let dgg = 4.0lf * r * r / (SQRT_PI * s_a * s_a * s_a) * ex;
    let ir = 1.0lf / r;
    let ir3 = ir * ir * ir;
    return vec2<f64>(gg * ir3, dgg * ir3 - 3.0lf * gg * ir3 * ir);
}

// The same for a charge and a dipole.
fn cd_sr(r: f64, s_a: f64) -> vec2<f64> {
    return charge_terms(r, s_a) - charge_terms(r, sc[6]);
}

fn pos_of(i: u32) -> vec3<f64> { return vec3<f64>(pol_pos[3u * i], pol_pos[3u * i + 1u], pol_pos[3u * i + 2u]); }
fn chpos_of(i: u32) -> vec3<f64> { return vec3<f64>(ch_pos[3u * i], ch_pos[3u * i + 1u], ch_pos[3u * i + 2u]); }
fn mu_of(i: u32) -> vec3<f64> { return vec3<f64>(mu[3u * i], mu[3u * i + 1u], mu[3u * i + 2u]); }
fn pv_of(i: u32) -> vec3<f64> { return vec3<f64>(pv[3u * i], pv[3u * i + 1u], pv[3u * i + 2u]); }
fn ddr_of(e: u32) -> vec3<f64> { return vec3<f64>(dd_r[3u * e], dd_r[3u * e + 1u], dd_r[3u * e + 2u]); }
fn cdr_of(e: u32) -> vec3<f64> { return vec3<f64>(cd_r[3u * e], cd_r[3u * e + 1u], cd_r[3u * e + 2u]); }
fn icdr_of(e: u32) -> vec3<f64> { return vec3<f64>(icd_r[3u * e], icd_r[3u * e + 1u], icd_r[3u * e + 2u]); }
fn inr_of(e: u32) -> vec3<f64> { return vec3<f64>(in_r[3u * e], in_r[3u * e + 1u], in_r[3u * e + 2u]); }
fn ddf_of(e: u32) -> vec3<f64> { return vec3<f64>(dd_f[3u * e], dd_f[3u * e + 1u], dd_f[3u * e + 2u]); }
fn cdf_of(e: u32) -> vec3<f64> { return vec3<f64>(cd_f[3u * e], cd_f[3u * e + 1u], cd_f[3u * e + 2u]); }

fn norm3(v: vec3<f64>) -> f64 {
    return sqrt(dot(v, v));
}

fn tensor_apply(a: f64, bp: f64, r: vec3<f64>, v: vec3<f64>) -> vec3<f64> {
    return v * a - r * (bp * dot(r, v));
}

fn pair_index(a: u32, b: u32, np: u32) -> u32 {
    return a * (2u * np - a - 1u) / 2u + (b - a - 1u);
}

// ---- term construction ------------------------------------------------------

@compute @workgroup_size(64)
fn k_terms_dd(@builtin(global_invocation_id) gid: vec3<u32>) {
    let e = gid.x;
    let np = p.np;
    if (e >= p.nl * np * np) {
        return;
    }
    let k = e / (np * np);
    let a = (e / np) % np;
    let b = e % np;
    let i = lij[2u * k];
    let j = lij[2u * k + 1u];
    let sh = vec3<f64>(lf[8u * k], lf[8u * k + 1u], lf[8u * k + 2u]);
    let sa = i * np + a;
    let sb = j * np + b;
    let r = pos_of(sa) - (pos_of(sb) + sh);
    let s = sqrt(2.0lf * (pol_sig[sa] * pol_sig[sa] + pol_sig[sb] * pol_sig[sb]));
    let t = dd_sr(norm3(r), s);
    dd_t[4u * e] = t.x;
    dd_t[4u * e + 1u] = t.y;
    dd_t[4u * e + 2u] = t.z;
    dd_t[4u * e + 3u] = t.w;
    dd_r[3u * e] = r.x;
    dd_r[3u * e + 1u] = r.y;
    dd_r[3u * e + 2u] = r.z;
}

@compute @workgroup_size(64)
fn k_terms_cd(@builtin(global_invocation_id) gid: vec3<u32>) {
    let e = gid.x;
    let np = p.np;
    let nch = p.nch;
    let per = np * nch;
    if (e >= p.nl * 2u * per) {
        return;
    }
    let k = e / (2u * per);
    let rest = e % (2u * per);
    let second = rest >= per;
    let q = rest % per;
    let a = q / nch;
    let c = q % nch;
    let i = lij[2u * k];
    let j = lij[2u * k + 1u];
    let sh = vec3<f64>(lf[8u * k], lf[8u * k + 1u], lf[8u * k + 2u]);
    var r: vec3<f64>;
    var sdip: f64;
    var cs: f64;
    if (!second) {
        // dipole a on i, charge c on j
        let sa = i * np + a;
        let cc = j * nch + c;
        r = pos_of(sa) - (chpos_of(cc) + sh);
        sdip = pol_sig[sa];
        cs = ch_sig[cc];
    } else {
        // dipole b (= a here) on j, charge c on i
        let sb = j * np + a;
        let cc = i * nch + c;
        r = (pos_of(sb) + sh) - chpos_of(cc);
        sdip = pol_sig[sb];
        cs = ch_sig[cc];
    }
    // A dipole and a charge at one point (never between two molecules) are left out, as the reference does.
    var t = vec2<f64>(0.0lf, 0.0lf);
    if (norm3(r) >= 1.0e-8lf) {
        t = cd_sr(norm3(r), sqrt(2.0lf * sdip * sdip + cs * cs));
    }
    cd_h[2u * e] = t.x;
    cd_h[2u * e + 1u] = t.y;
    cd_r[3u * e] = r.x;
    cd_r[3u * e + 1u] = r.y;
    cd_r[3u * e + 2u] = r.z;
}

@compute @workgroup_size(64)
fn k_terms_intra(@builtin(global_invocation_id) gid: vec3<u32>) {
    let e = gid.x;
    if (e >= p.nm * p.npairs) {
        return;
    }
    let m = e / p.npairs;
    let pr = e % p.npairs;
    let np = p.np;
    var a = 0u;
    var b = 1u;
    var idx = 0u;
    for (var x = 0u; x < np; x = x + 1u) {
        for (var y = x + 1u; y < np; y = y + 1u) {
            if (idx == pr) {
                a = x;
                b = y;
            }
            idx = idx + 1u;
        }
    }
    let sa = m * np + a;
    let sb = m * np + b;
    let r = pos_of(sa) - pos_of(sb);
    let s = sqrt(2.0lf * (pol_sig[sa] * pol_sig[sa] + pol_sig[sb] * pol_sig[sb]));
    let t = dd_sr(norm3(r), s);
    in_t[4u * e] = t.x;
    in_t[4u * e + 1u] = t.y;
    in_t[4u * e + 2u] = t.z;
    in_t[4u * e + 3u] = t.w;
    in_r[3u * e] = r.x;
    in_r[3u * e + 1u] = r.y;
    in_r[3u * e + 2u] = r.z;
}

// The molecule's own charges act on its own dipoles in the reciprocal sum and not in the liquid:
// minus the tensor at the Ewald width. One entry for each (dipole a, charge c) of a molecule.
@compute @workgroup_size(64)
fn k_terms_icd(@builtin(global_invocation_id) gid: vec3<u32>) {
    let e = gid.x;
    let np = p.np;
    let nch = p.nch;
    if (e >= p.nm * np * nch) {
        return;
    }
    let m = e / (np * nch);
    let a = (e / nch) % np;
    let c = e % nch;
    let r = pos_of(m * np + a) - chpos_of(m * nch + c);
    var t = vec2<f64>(0.0lf, 0.0lf);
    if (norm3(r) >= 1.0e-8lf) {
        t = -charge_terms(norm3(r), sc[6]);
    }
    icd_h[2u * e] = t.x;
    icd_h[2u * e + 1u] = t.y;
    icd_r[3u * e] = r.x;
    icd_r[3u * e + 1u] = r.y;
    icd_r[3u * e + 2u] = r.z;
}

// ---- reductions -------------------------------------------------------------

var<workgroup> red3: array<vec3<f64>, 32>;
var<workgroup> part: array<f64, 256>;

// The sum of 32 threads' vectors by a fixed tree.
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

fn reduce(v: f64, li: u32) -> f64 {
    part[li] = v;
    workgroupBarrier();
    for (var s = 128u; s > 0u; s = s / 2u) {
        if (li < s) {
            part[li] = part[li] + part[li + s];
        }
        workgroupBarrier();
    }
    let total = part[0];
    workgroupBarrier();
    return total;
}

// ---- the field of the other molecules' charges on each dipole ---------------

@compute @workgroup_size(32)
fn k_e0(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let s = wg.x;
    let li = lid.x;
    let np = p.np;
    let nch = p.nch;
    let m = s / np;
    let a = s % np;
    var acc = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    let start = inc[m];
    let count = inc[m + 1u] - start;
    for (var t = li; t < count; t = t + 32u) {
        let item = inc[p.nm + 1u + start + t];
        let k = item / 2u;
        let side = item % 2u;
        let w = lf[8u * k + 3u];
        let owner = lij[2u * k + 1u - side];
        for (var c = 0u; c < nch; c = c + 1u) {
            let e = k * 2u * np * nch + side * np * nch + a * nch + c;
            acc = acc + cdr_of(e) * (w * ch_q[owner * nch + c] * cd_h[2u * e]);
        }
    }
    var tot = reduce3(acc, li);
    if (li == 0u) {
        tot = tot + vec3<f64>(e0x[3u * s], e0x[3u * s + 1u], e0x[3u * s + 2u]);
        for (var c = 0u; c < nch; c = c + 1u) {
            let ee = s * nch + c;
            tot = tot + icdr_of(ee) * (ch_q[m * nch + c] * icd_h[2u * ee]);
        }
        e0[3u * s] = tot.x;
        e0[3u * s + 1u] = tot.y;
        e0[3u * s + 2u] = tot.z;
    }
}

// ---- the matrix-vector product ----------------------------------------------

fn xv(src: u32, i: u32) -> vec3<f64> {
    if (src == 0u) {
        return mu_of(i);
    }
    return pv_of(i);
}

fn matvec_site(s: u32, src: u32, li: u32) -> vec3<f64> {
    let np = p.np;
    let m = s / np;
    let a = s % np;
    var y = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    if (li == 0u) {
        y = xv(src, s) * inv_alpha[s];
        for (var b = 0u; b < np; b = b + 1u) {
            if (b == a) {
                continue;
            }
            var lo = a;
            var hi = b;
            if (b < a) {
                lo = b;
                hi = a;
            }
            let e = m * p.npairs + pair_index(lo, hi, np);
            y = y + tensor_apply(in_t[4u * e], in_t[4u * e + 1u], inr_of(e), xv(src, m * np + b));
        }
    }
    let start = inc[m];
    let count = inc[m + 1u] - start;
    for (var t = li; t < count; t = t + 32u) {
        let item = inc[p.nm + 1u + start + t];
        let k = item / 2u;
        let side = item % 2u;
        let w = lf[8u * k + 3u];
        if (side == 0u) {
            let j = lij[2u * k + 1u];
            for (var b = 0u; b < np; b = b + 1u) {
                let e = k * np * np + a * np + b;
                y = y + tensor_apply(w * dd_t[4u * e], w * dd_t[4u * e + 1u], ddr_of(e), xv(src, j * np + b));
            }
        } else {
            let i = lij[2u * k];
            for (var a2 = 0u; a2 < np; a2 = a2 + 1u) {
                let e = k * np * np + a2 * np + a;
                y = y + tensor_apply(w * dd_t[4u * e], w * dd_t[4u * e + 1u], ddr_of(e), xv(src, i * np + a2));
            }
        }
    }
    return reduce3(y, li);
}

@compute @workgroup_size(32)
fn k_matvec_mu(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let y = matvec_site(wg.x, 0u, lid.x);
    if (lid.x == 0u) {
        apv[3u * wg.x] = y.x;
        apv[3u * wg.x + 1u] = y.y;
        apv[3u * wg.x + 2u] = y.z;
    }
}

@compute @workgroup_size(32)
fn k_matvec_p(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let y = matvec_site(wg.x, 1u, lid.x);
    if (lid.x == 0u && sc[3] == 0.0lf) {
        apv[3u * wg.x] = y.x;
        apv[3u * wg.x + 1u] = y.y;
        apv[3u * wg.x + 2u] = y.z;
    }
}

// ---- conjugate gradients: one workgroup, fixed-order reductions -------------

// r = e0 - A mu (apv holds A mu), z = alpha r, p = z, rz, |e0|; the first
// convergence test, with the iteration count it implies.
@compute @workgroup_size(256)
fn k_cg_init(@builtin(local_invocation_id) lid: vec3<u32>) {
    let li = lid.x;
    var rz = 0.0lf;
    var b2 = 0.0lf;
    var r2 = 0.0lf;
    for (var x = li; x < 3u * p.n; x = x + 256u) {
        let r = e0[x] - apv[x];
        let z = r / inv_alpha[x / 3u];
        rv[x] = r;
        zv[x] = z;
        pv[x] = z;
        rz = rz + r * z;
        r2 = r2 + r * r;
        b2 = b2 + e0[x] * e0[x];
    }
    let rz_t = reduce(rz, li);
    let b2_t = reduce(b2, li);
    let r2_t = reduce(r2, li);
    if (li == 0u) {
        let bnorm = max(sqrt(b2_t), 1.0e-300lf);
        sc[1] = bnorm;
        sc[2] = rz_t;
        sc[5] = sqrt(r2_t);
        if (sqrt(r2_t) <= sc[0] * bnorm) {
            sc[3] = 1.0lf;
            sc[4] = 0.0lf;
        } else {
            sc[3] = 0.0lf;
            sc[4] = 1.0lf;
        }
    }
}

// One conjugate-gradient iteration after the product, and the convergence test
// that follows it (the reference's test at the head of the next).
@compute @workgroup_size(256)
fn k_cg_b(@builtin(local_invocation_id) lid: vec3<u32>) {
    let li = lid.x;
    let live = sc[3] == 0.0lf;
    var pap = 0.0lf;
    for (var x = li; x < 3u * p.n; x = x + 256u) {
        pap = pap + pv[x] * apv[x];
    }
    let pap_t = reduce(pap, li);
    let go = live && pap_t > 0.0lf;
    let rz = sc[2];
    let alpha = rz / pap_t;
    var rz_new = 0.0lf;
    var r2 = 0.0lf;
    for (var x = li; x < 3u * p.n; x = x + 256u) {
        if (go) {
            mu[x] = mu[x] + alpha * pv[x];
            rv[x] = rv[x] - alpha * apv[x];
            zv[x] = rv[x] / inv_alpha[x / 3u];
        }
        rz_new = rz_new + rv[x] * zv[x];
        r2 = r2 + rv[x] * rv[x];
    }
    let rz_t = reduce(rz_new, li);
    let r2_t = reduce(r2, li);
    let beta = rz_t / rz;
    for (var x = li; x < 3u * p.n; x = x + 256u) {
        if (go) {
            pv[x] = zv[x] + beta * pv[x];
        }
    }
    if (li == 0u && live) {
        if (go) {
            sc[2] = rz_t;
            let rn = sqrt(r2_t);
            sc[5] = rn;
            if (rn <= sc[0] * sc[1]) {
                sc[3] = 1.0lf;
            } else if (sc[4] >= f64(p.max_iterations)) {
                sc[3] = 1.0lf;
            } else {
                sc[4] = sc[4] + 1.0lf;
            }
        } else {
            sc[3] = 1.0lf;
        }
    }
}

// ---- forces -----------------------------------------------------------------

@compute @workgroup_size(64)
fn k_link_force(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= p.nl) {
        return;
    }
    let np = p.np;
    let nch = p.nch;
    let i = lij[2u * k];
    let j = lij[2u * k + 1u];
    let w = lf[8u * k + 3u];
    let dw = lf[8u * k + 4u];
    var u = 0.0lf;
    var fi = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    for (var a = 0u; a < np; a = a + 1u) {
        for (var b = 0u; b < np; b = b + 1u) {
            let e = k * np * np + a * np + b;
            let ma = mu_of(i * np + a);
            let mb = mu_of(j * np + b);
            let r = ddr_of(e);
            let rr = norm3(r);
            let ta = dd_t[4u * e];
            let tbp = dd_t[4u * e + 1u];
            let da = dd_t[4u * e + 2u];
            let dbp = dd_t[4u * e + 3u];
            let mamb = dot(ma, mb);
            let mar = dot(ma, r);
            let mbr = dot(mb, r);
            u = u + mamb * ta - tbp * mar * mbr;
            let grad = r * (da / rr * mamb) - r * (dbp / rr * mar * mbr) - ma * (tbp * mbr) - mb * (tbp * mar);
            let f = grad * (-w);
            dd_f[3u * e] = f.x;
            dd_f[3u * e + 1u] = f.y;
            dd_f[3u * e + 2u] = f.z;
            fi = fi + f;
        }
    }
    for (var s = 0u; s < 2u; s = s + 1u) {
        for (var a = 0u; a < np; a = a + 1u) {
            for (var c = 0u; c < nch; c = c + 1u) {
                let e = k * 2u * np * nch + s * np * nch + a * nch + c;
                var m: vec3<f64>;
                var q: f64;
                if (s == 0u) {
                    m = mu_of(i * np + a);
                    q = ch_q[j * nch + c];
                } else {
                    m = mu_of(j * np + a);
                    q = ch_q[i * nch + c];
                }
                let r = cdr_of(e);
                let rr = norm3(r);
                let h = cd_h[2u * e];
                let dh = cd_h[2u * e + 1u];
                let mr = dot(m, r);
                u = u - q * h * mr;
                let grad = r * (-q * dh / rr * mr) - m * (q * h);
                let f = grad * (-w);
                cd_f[3u * e] = f.x;
                cd_f[3u * e + 1u] = f.y;
                cd_f[3u * e + 2u] = f.z;
                if (s == 0u) {
                    fi = fi + f;
                } else {
                    fi = fi - f;
                }
            }
        }
    }
    let d = vec3<f64>(lf[8u * k + 5u], lf[8u * k + 6u], lf[8u * k + 7u]);
    let dist = norm3(d);
    var fw = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    if (dist > 0.0lf) {
        fw = d * (-u * dw / dist);
    }
    lforce[3u * k] = fw.x;
    lforce[3u * k + 1u] = fw.y;
    lforce[3u * k + 2u] = fw.z;
    lvir[k] = dot(d, fi + fw);
}

// Force on each dipole site: the pair terms of every link it is in, and its molecule's own.
@compute @workgroup_size(32)
fn k_gather_pol(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let s = wg.x;
    let li = lid.x;
    let np = p.np;
    let nch = p.nch;
    let m = s / np;
    let a = s % np;
    var f = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    if (li == 0u) {
        // The molecule's own: fpol[a] -= grad, fpol[b] += grad for a < b.
        for (var b = 0u; b < np; b = b + 1u) {
            if (b == a) {
                continue;
            }
            var lo = a;
            var hi = b;
            if (b < a) {
                lo = b;
                hi = a;
            }
            let e = m * p.npairs + pair_index(lo, hi, np);
            let ma = mu_of(m * np + lo);
            let mb = mu_of(m * np + hi);
            let r = inr_of(e);
            let rr = norm3(r);
            let da = in_t[4u * e + 2u];
            let tbp = in_t[4u * e + 1u];
            let dbp = in_t[4u * e + 3u];
            let mamb = dot(ma, mb);
            let mar = dot(ma, r);
            let mbr = dot(mb, r);
            let grad = r * (da / rr * mamb) - r * (dbp / rr * mar * mbr) - ma * (tbp * mbr) - mb * (tbp * mar);
            if (a < b) {
                f = f - grad;
            } else {
                f = f + grad;
            }
        }
        for (var c = 0u; c < nch; c = c + 1u) {
            let ee = s * nch + c;
            let r = icdr_of(ee);
            let rr = norm3(r);
            if (rr >= 1.0e-8lf) {
                let md = mu_of(s);
                let q = ch_q[m * nch + c];
                let grad = r * (-q * icd_h[2u * ee + 1u] / rr * dot(md, r)) - md * (q * icd_h[2u * ee]);
                f = f - grad;
            }
        }
    }
    let start = inc[m];
    let count = inc[m + 1u] - start;
    for (var t = li; t < count; t = t + 32u) {
        let item = inc[p.nm + 1u + start + t];
        let k = item / 2u;
        let side = item % 2u;
        if (side == 0u) {
            for (var b = 0u; b < np; b = b + 1u) {
                f = f + ddf_of(k * np * np + a * np + b);
            }
            for (var c = 0u; c < nch; c = c + 1u) {
                f = f + cdf_of(k * 2u * np * nch + a * nch + c);
            }
        } else {
            for (var a2 = 0u; a2 < np; a2 = a2 + 1u) {
                f = f - ddf_of(k * np * np + a2 * np + a);
            }
            for (var c = 0u; c < nch; c = c + 1u) {
                f = f + cdf_of(k * 2u * np * nch + np * nch + a * nch + c);
            }
        }
    }
    let tot = reduce3(f, li);
    if (li == 0u) {
        fpol[3u * s] = tot.x;
        fpol[3u * s + 1u] = tot.y;
        fpol[3u * s + 2u] = tot.z;
    }
}

// Force on each charge: minus what the dipoles that feel it receive.
@compute @workgroup_size(32)
fn k_gather_charge(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let cc = wg.x;
    let li = lid.x;
    let np = p.np;
    let nch = p.nch;
    let m = cc / nch;
    let c = cc % nch;
    var f = vec3<f64>(0.0lf, 0.0lf, 0.0lf);
    if (li == 0u) {
        // The molecule's own dipoles and this charge: what they receive, this charge gives back.
        for (var a = 0u; a < np; a = a + 1u) {
            let ee = (m * np + a) * nch + c;
            let r = icdr_of(ee);
            let rr = norm3(r);
            if (rr >= 1.0e-8lf) {
                let md = mu_of(m * np + a);
                let q = ch_q[cc];
                let grad = r * (-q * icd_h[2u * ee + 1u] / rr * dot(md, r)) - md * (q * icd_h[2u * ee]);
                f = f + grad;
            }
        }
    }
    let start = inc[m];
    let count = inc[m + 1u] - start;
    for (var t = li; t < count; t = t + 32u) {
        let item = inc[p.nm + 1u + start + t];
        let k = item / 2u;
        let side = item % 2u;
        for (var a = 0u; a < np; a = a + 1u) {
            // The first molecule's charges meet the second's dipoles (cd2), and the reverse.
            f = f - cdf_of(k * 2u * np * nch + (1u - side) * np * nch + a * nch + c);
        }
    }
    let tot = reduce3(f, li);
    if (li == 0u) {
        fch[3u * cc] = tot.x;
        fch[3u * cc + 1u] = tot.y;
        fch[3u * cc + 2u] = tot.z;
    }
}

// ---- the reciprocal sums ----------------------------------------------------
//
// One workgroup of 64 threads a wave vector (the structure factors, the arm sums) or a
// site (the forces and the field); each thread takes every 64th term and the 64 partial sums
// are added by a fixed tree. The phases are products of the CPU's per-axis tables of
// e^(i m theta), as `Ewald::phases` forms them. See `ewald.wgsl`, the standalone version.
// The dipole vector the sums act on is chosen by `src`: 0 the converged or starting dipoles,
// 1 the conjugate-gradient direction.

var<workgroup> r0: array<f64, 64>;
var<workgroup> r1: array<f64, 64>;
var<workgroup> r2: array<f64, 64>;

fn rreduce3(t: u32) {
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

fn cmul(a: vec2<f64>, b: vec2<f64>) -> vec2<f64> {
    return vec2<f64>(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

fn tab_get(is_charge: bool, site: u32, axis: u32, m: i32) -> vec2<f64> {
    let ns = select(p.n, rp.nc, is_charge);
    var n = rp.n0;
    var off = 0u;
    if (axis == 1u) {
        n = rp.n1;
        off = ns * rp.n0 * 2u;
    } else if (axis == 2u) {
        n = rp.n2;
        off = ns * (rp.n0 + rp.n1) * 2u;
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

fn phase(is_charge: bool, site: u32, j: u32) -> vec2<f64> {
    let tx = tab_get(is_charge, site, 0u, kidx[3u * j]);
    let ty = tab_get(is_charge, site, 1u, kidx[3u * j + 1u]);
    let tz = tab_get(is_charge, site, 2u, kidx[3u * j + 2u]);
    return cmul(cmul(tx, ty), tz);
}

fn kdot(j: u32, x: f64, y: f64, z: f64) -> f64 {
    return kv[3u * j] * x + kv[3u * j + 1u] * y + kv[3u * j + 2u] * z;
}

fn rvec(src: u32, a: u32) -> vec3<f64> {
    if (src == 0u) {
        return mu_of(a);
    }
    return pv_of(a);
}

@compute @workgroup_size(64)
fn r_sq(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let j = wg.x;
    var sr = 0.0lf;
    var si = 0.0lf;
    for (var i = t; i < rp.nc; i = i + 64u) {
        let ph = phase(true, i, j);
        let qi = rq[i];
        sr = sr + qi * ph.x;
        si = si + qi * ph.y;
    }
    r0[t] = sr;
    r1[t] = si;
    r2[t] = 0.0lf;
    rreduce3(t);
    if (t == 0u) {
        sq[2u * j] = r0[0];
        sq[2u * j + 1u] = r1[0];
    }
}

fn dsum(j: u32, t: u32, src: u32) {
    var sr = 0.0lf;
    var si = 0.0lf;
    for (var a = t; a < p.n; a = a + 64u) {
        let ph = phase(false, a, j);
        let m = rvec(src, a);
        let km = kdot(j, m.x, m.y, m.z);
        sr = sr + (-km * ph.y);
        si = si + km * ph.x;
    }
    r0[t] = sr;
    r1[t] = si;
    r2[t] = 0.0lf;
    rreduce3(t);
    if (t == 0u) {
        sd[2u * j] = r0[0];
        sd[2u * j + 1u] = r1[0];
    }
}

@compute @workgroup_size(64)
fn r_dsum_mu(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    dsum(wg.x, t, 0u);
}

@compute @workgroup_size(64)
fn r_dsum_p(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    dsum(wg.x, t, 1u);
}

@compute @workgroup_size(64)
fn r_stot(@builtin(global_invocation_id) g: vec3<u32>) {
    let j = g.x;
    if (j >= rp.nk) {
        return;
    }
    stot[2u * j] = sq[2u * j] + sd[2u * j];
    stot[2u * j + 1u] = sq[2u * j + 1u] + sd[2u * j + 1u];
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

fn grad_dip(is_charge: bool, site: u32, which: u32, t: u32) -> vec3<f64> {
    var gx = 0.0lf;
    var gy = 0.0lf;
    var gz = 0.0lf;
    for (var j = t; j < rp.nk; j = j + 64u) {
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
    rreduce3(t);
    let pref = sc[7];
    return vec3<f64>(r0[0] * pref, r1[0] * pref, r2[0] * pref);
}

// The permanent charges' long-range field at the dipoles: minus grad_dipole(S_q).
@compute @workgroup_size(64)
fn r_field(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let a = wg.x;
    let v = grad_dip(false, a, 0u, t);
    if (t == 0u) {
        e0x[3u * a] = -v.x;
        e0x[3u * a + 1u] = -v.y;
        e0x[3u * a + 2u] = -v.z;
    }
}

// Add the long-range dipole-dipole coupling less each dipole's own to the product A x.
fn add_resp(a: u32, t: u32, src: u32) {
    let v = grad_dip(false, a, 1u, t);
    if (t == 0u) {
        let own = sc[8];
        let x = rvec(src, a);
        apv[3u * a] = apv[3u * a] + (v.x - x.x * own);
        apv[3u * a + 1u] = apv[3u * a + 1u] + (v.y - x.y * own);
        apv[3u * a + 2u] = apv[3u * a + 2u] + (v.z - x.z * own);
    }
}

@compute @workgroup_size(64)
fn r_resp_mu(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    add_resp(wg.x, t, 0u);
}

// Only while the iteration is live, as the product itself is.
@compute @workgroup_size(64)
fn r_resp_p(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    if (sc[3] == 0.0lf) {
        add_resp(wg.x, t, 1u);
    }
}

@compute @workgroup_size(64)
fn r_gch(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let i = wg.x;
    let v = grad_dip(true, i, 2u, t);
    if (t == 0u) {
        let qi = rq[i];
        o_gch[3u * i] = v.x * qi;
        o_gch[3u * i + 1u] = v.y * qi;
        o_gch[3u * i + 2u] = v.z * qi;
    }
}

@compute @workgroup_size(64)
fn r_gdpos(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let a = wg.x;
    let m = mu_of(a);
    var gx = 0.0lf;
    var gy = 0.0lf;
    var gz = 0.0lf;
    for (var j = t; j < rp.nk; j = j + 64u) {
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
    rreduce3(t);
    if (t == 0u) {
        let pref = sc[7];
        o_gdp[3u * a] = r0[0] * pref;
        o_gdp[3u * a + 1u] = r1[0] * pref;
        o_gdp[3u * a + 2u] = r2[0] * pref;
    }
}

@compute @workgroup_size(64)
fn r_arm_c(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let j = wg.x;
    var sr = 0.0lf;
    var si = 0.0lf;
    for (var i = t; i < rp.nc; i = i + 64u) {
        let ph = phase(true, i, j);
        let c = rq[i] * kdot(j, carm[3u * i], carm[3u * i + 1u], carm[3u * i + 2u]);
        sr = sr + c * ph.x;
        si = si + c * ph.y;
    }
    r0[t] = sr;
    r1[t] = si;
    r2[t] = 0.0lf;
    rreduce3(t);
    if (t == 0u) {
        tcv[2u * j] = r0[0];
        tcv[2u * j + 1u] = r1[0];
    }
}

@compute @workgroup_size(64)
fn r_arm_d(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) t: u32) {
    let j = wg.x;
    var sr = 0.0lf;
    var si = 0.0lf;
    for (var a = t; a < p.n; a = a + 64u) {
        let ph = phase(false, a, j);
        let m = mu_of(a);
        let c = kdot(j, m.x, m.y, m.z) * kdot(j, darm[3u * a], darm[3u * a + 1u], darm[3u * a + 2u]);
        sr = sr + c * ph.x;
        si = si + c * ph.y;
    }
    r0[t] = sr;
    r1[t] = si;
    r2[t] = 0.0lf;
    rreduce3(t);
    if (t == 0u) {
        tdv[2u * j] = r0[0];
        tdv[2u * j + 1u] = r1[0];
    }
}
