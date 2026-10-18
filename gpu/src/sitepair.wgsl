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

// Shared by the GPU kernels that need the reference's special functions in
// double precision: WGSL has no f64 exp, so `gexp` is built here; `gerf` is the
// reference's own series. Included by gen_induction.py in each template at
// `@COMMON@`; edit it here.

const SQRT_PI = 1.7724538509055160273lf;
const LOG2E = 1.4426950408889634074lf;
const LN2_HI = 0.693147180369123816490lf;
const LN2_LO = 1.90821492927058770002e-10lf;
const INV_ODD = array<f64, 130>(1.0lf, 0.33333333333333331lf, 0.20000000000000001lf, 0.14285714285714285lf, 0.1111111111111111lf, 0.090909090909090912lf, 0.076923076923076927lf, 0.066666666666666666lf, 0.058823529411764705lf, 0.052631578947368418lf, 0.047619047619047616lf, 0.043478260869565216lf, 0.040000000000000001lf, 0.037037037037037035lf, 0.034482758620689655lf, 0.032258064516129031lf, 0.030303030303030304lf, 0.028571428571428571lf, 0.027027027027027029lf, 0.02564102564102564lf, 0.024390243902439025lf, 0.023255813953488372lf, 0.022222222222222223lf, 0.021276595744680851lf, 0.020408163265306121lf, 0.019607843137254902lf, 0.018867924528301886lf, 0.018181818181818181lf, 0.017543859649122806lf, 0.016949152542372881lf, 0.016393442622950821lf, 0.015873015873015872lf, 0.015384615384615385lf, 0.014925373134328358lf, 0.014492753623188406lf, 0.014084507042253521lf, 0.013698630136986301lf, 0.013333333333333334lf, 0.012987012987012988lf, 0.012658227848101266lf, 0.012345679012345678lf, 0.012048192771084338lf, 0.011764705882352941lf, 0.011494252873563218lf, 0.011235955056179775lf, 0.01098901098901099lf, 0.010752688172043012lf, 0.010526315789473684lf, 0.010309278350515464lf, 0.010101010101010102lf, 0.0099009900990099011lf, 0.0097087378640776691lf, 0.0095238095238095247lf, 0.0093457943925233638lf, 0.0091743119266055051lf, 0.0090090090090090089lf, 0.0088495575221238937lf, 0.0086956521739130436lf, 0.0085470085470085479lf, 0.0084033613445378148lf, 0.0082644628099173556lf, 0.008130081300813009lf, 0.0080000000000000002lf, 0.007874015748031496lf, 0.0077519379844961239lf, 0.0076335877862595417lf, 0.0075187969924812026lf, 0.0074074074074074077lf, 0.0072992700729927005lf, 0.0071942446043165471lf, 0.0070921985815602835lf, 0.006993006993006993lf, 0.0068965517241379309lf, 0.0068027210884353739lf, 0.0067114093959731542lf, 0.0066225165562913907lf, 0.0065359477124183009lf, 0.0064516129032258064lf, 0.006369426751592357lf, 0.0062893081761006293lf, 0.0062111801242236021lf, 0.0061349693251533744lf, 0.0060606060606060606lf, 0.0059880239520958087lf, 0.0059171597633136093lf, 0.0058479532163742687lf, 0.0057803468208092483lf, 0.0057142857142857143lf, 0.0056497175141242938lf, 0.0055865921787709499lf, 0.0055248618784530384lf, 0.0054644808743169399lf, 0.0054054054054054057lf, 0.0053475935828877002lf, 0.0052910052910052907lf, 0.005235602094240838lf, 0.0051813471502590676lf, 0.0051282051282051282lf, 0.005076142131979695lf, 0.0050251256281407036lf, 0.0049751243781094526lf, 0.0049261083743842365lf, 0.0048780487804878049lf, 0.004830917874396135lf, 0.0047846889952153108lf, 0.0047393364928909956lf, 0.0046948356807511738lf, 0.0046511627906976744lf, 0.004608294930875576lf, 0.0045662100456621002lf, 0.0045248868778280547lf, 0.0044843049327354259lf, 0.0044444444444444444lf, 0.0044052863436123352lf, 0.0043668122270742356lf, 0.004329004329004329lf, 0.0042918454935622317lf, 0.0042553191489361703lf, 0.0042194092827004216lf, 0.0041841004184100415lf, 0.0041493775933609959lf, 0.00411522633744856lf, 0.0040816326530612249lf, 0.004048582995951417lf, 0.0040160642570281121lf, 0.0039840637450199202lf, 0.003952569169960474lf, 0.0039215686274509803lf, 0.0038910505836575876lf, 0.0038610038610038611lf);

// exp(x) for x <= 0 to full double precision: x = n ln2 + r, |r| <= ln2/2, a
// degree-14 polynomial for e^r, and 2^n by exact repeated squaring.
fn gexp(x: f64) -> f64 {
    if (x < -700.0lf) {
        return 0.0lf;
    }
    let n = floor(x * LOG2E + 0.5lf);
    let r = (x - n * LN2_HI) - n * LN2_LO;
    var poly = 1.1470745597729725e-11lf;
    poly = poly * r + 1.6059043836821613e-10lf;
    poly = poly * r + 2.08767569878681e-09lf;
    poly = poly * r + 2.505210838544172e-08lf;
    poly = poly * r + 2.7557319223985888e-07lf;
    poly = poly * r + 2.7557319223985893e-06lf;
    poly = poly * r + 2.4801587301587302e-05lf;
    poly = poly * r + 0.00019841269841269841lf;
    poly = poly * r + 0.0013888888888888889lf;
    poly = poly * r + 0.0083333333333333332lf;
    poly = poly * r + 0.041666666666666664lf;
    poly = poly * r + 0.16666666666666666lf;
    poly = poly * r + 0.5lf;
    poly = poly * r + 1.0lf;
    poly = poly * r + 1.0lf;
    var m = i32(n);
    var base = 2.0lf;
    if (m < 0) {
        base = 0.5lf;
        m = -m;
    }
    var scale = 1.0lf;
    for (var b = 0; b < 11; b = b + 1) {
        if ((m & 1) != 0) {
            scale = scale * base;
        }
        base = base * base;
        m = m >> 1u;
    }
    return poly * scale;
}

// exp(-x^2), taken as 0 beyond x^2 = 40 (below 5e-18 of one, under the
// reference's own rounding at the places it is used).
fn gauss(x: f64) -> f64 {
    let x2 = x * x;
    if (x2 > 40.0lf) {
        return 0.0lf;
    }
    return gexp(-x2);
}

// The reference's erf: the series for |x| < 5.5, the sign beyond.
fn gerf(x: f64) -> f64 {
    let ax = abs(x);
    if (ax == 0.0lf) {
        return 0.0lf;
    }
    if (ax >= 5.5lf) {
        if (x < 0.0lf) {
            return -1.0lf;
        }
        return 1.0lf;
    }
    let x2 = ax * ax;
    let t2 = 2.0lf * x2;
    var term = ax;
    var sum = ax;
    for (var n = 1u; n < 129u; n = n + 1u) {
        term = term * t2 * INV_ODD[n];
        sum = sum + term;
        if (term <= 1.0e-17lf * sum) {
            break;
        }
    }
    var v = 2.0lf / SQRT_PI * gexp(-x2) * sum;
    if (x < 0.0lf) {
        v = -v;
    }
    return v;
}



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
    // The Ewald remainder: only erf(r / s) - erf(r / s') of the Coulomb term is here (`EwaldLaw`).
    let es = pair[4u * nt * nt];
    if (es > 0.0lf && qq != 0.0lf) {
        let xe = r / es;
        let ge = gerf(xe);
        let ce = 2.0lf / (SQRT_PI * es) * gauss(xe);
        coul = coul - qq * ge / r;
        dcoul = dcoul - qq * (ce / r - ge / (r * r));
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
