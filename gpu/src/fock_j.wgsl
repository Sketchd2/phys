// jp[p] = sum over k of W[p][k] dt[k] for the pairs of one segment: the Coulomb
// matrix once the density has been fitted. One workgroup of 64 threads a pair,
// the sum over k shared out and reduced in workgroup memory.

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
@group(0) @binding(2) var<storage, read> dt: array<f64>;
@group(0) @binding(3) var<storage, read_write> jp: array<f64>;

var<workgroup> partial: array<f64, 64>;

@compute @workgroup_size(64)
fn main(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) li: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let pi = p.p0 + wg.y * nwg.x + wg.x;
    if (pi >= p.p1) {
        return;
    }
    var acc = 0.0lf;
    for (var k = li.x; k < p.nk; k = k + 64u) {
        acc = acc + w[(pi - p.p0) * p.nk + k] * dt[k];
    }
    partial[li.x] = acc;
    workgroupBarrier();
    for (var s = 32u; s > 0u; s = s / 2u) {
        if (li.x < s) {
            partial[li.x] = partial[li.x] + partial[li.x + s];
        }
        workgroupBarrier();
    }
    if (li.x == 0u) {
        jp[pi] = partial[0];
    }
}
