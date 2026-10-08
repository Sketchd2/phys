// dt[k] += sum over the pairs p of one segment of dp[p] W[p][k]: the density
// fitted onto the whitened directions (the Coulomb matrix's right-hand side).
// `dp[p]` is the density element for the pair, doubled off the diagonal. One
// thread a direction, the pairs walked in order, so the reads along k are
// contiguous across a workgroup.

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
@group(0) @binding(2) var<storage, read> dp: array<f64>;
@group(0) @binding(3) var<storage, read_write> dt: array<f64>;

@compute @workgroup_size(64)
fn main(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) li: vec3<u32>) {
    let k = wg.x * 64u + li.x;
    if (k >= p.nk) {
        return;
    }
    var acc = 0.0lf;
    for (var pi = p.p0; pi < p.p1; pi = pi + 1u) {
        acc = acc + dp[pi] * w[(pi - p.p0) * p.nk + k];
    }
    dt[k] = dt[k] + acc;
}
