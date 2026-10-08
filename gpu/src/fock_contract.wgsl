// x[i][m][k] += sum over n of W(m, n)[k] c[n][i], for the pairs (m, n) of one
// segment of the resident whitened integrals. W is stored for m >= n only, in
// pair order p = m (m + 1) / 2 + n, `nk` doubles a pair, a segment holding the
// pairs p0 .. p1; the sum is over all n and takes the symmetric partner as it
// goes. One workgroup of 64 threads covers 64 of the k and one row m.

struct P {
    n: u32,
    nk: u32,
    count: u32,
    p0: u32,
    p1: u32,
    m0: u32,
    rows: u32,
    pad: u32,
};

@group(0) @binding(0) var<uniform> p: P;
@group(0) @binding(1) var<storage, read> w: array<f64>;
@group(0) @binding(2) var<storage, read> c: array<f64>;
@group(0) @binding(3) var<storage, read_write> x: array<f64>;

const MAXC: u32 = 32u;

fn pidx(a: u32, b: u32) -> u32 {
    if (a >= b) {
        return a * (a + 1u) / 2u + b;
    }
    return b * (b + 1u) / 2u + a;
}

@compute @workgroup_size(64)
fn main(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) li: vec3<u32>) {
    let k = wg.x * 64u + li.x;
    let m = p.m0 + wg.y;
    if (k >= p.nk || m >= p.n) {
        return;
    }
    var acc: array<f64, 32>;
    for (var i = 0u; i < MAXC; i = i + 1u) {
        acc[i] = 0.0lf;
    }
    for (var nn = 0u; nn < p.n; nn = nn + 1u) {
        let pi = pidx(m, nn);
        if (pi < p.p0 || pi >= p.p1) {
            continue;
        }
        let wv = w[(pi - p.p0) * p.nk + k];
        for (var i = 0u; i < p.count; i = i + 1u) {
            acc[i] = acc[i] + wv * c[nn * p.count + i];
        }
    }
    for (var i = 0u; i < p.count; i = i + 1u) {
        let idx = (i * p.n + m) * p.nk + k;
        x[idx] = x[idx] + acc[i];
    }
}
