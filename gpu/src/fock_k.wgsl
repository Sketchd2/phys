// kmat[a][b] += sum over orbitals i and directions k of x[i][a][k] x[i][b][k],
// for the lower triangle (a >= b) of 64 x 64 tiles, in double precision: the
// exchange matrix from the half-transformed integrals. The tile rows
// m0 .. m0 + rows only, so that a dispatch stays short. Each orbital's block
// is a product of the same shape as `gemm.wgsl`'s (sixteen k at a time through
// workgroup memory, a 4 x 4 block a thread).

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
@group(0) @binding(1) var<storage, read> x: array<f64>;
@group(0) @binding(2) var<storage, read_write> kmat: array<f64>;

const TILE: u32 = 64u;
const KT: u32 = 16u;
const PITCH: u32 = 17u;

var<workgroup> ta: array<f64, 1088>;
var<workgroup> tb: array<f64, 1088>;

@compute @workgroup_size(16, 16)
fn main(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) li: vec3<u32>) {
    let tx = li.x;
    let ty = li.y;
    let tid = ty * 16u + tx;
    let row_tile = p.m0 + wg.y;
    let col_tile = wg.x;
    if (col_tile > row_tile) {
        return;
    }
    let row_base = row_tile * TILE;
    let col_base = col_tile * TILE;
    var acc: array<f64, 16>;
    for (var q = 0u; q < 16u; q = q + 1u) {
        acc[q] = 0.0lf;
    }
    for (var i = 0u; i < p.count; i = i + 1u) {
        let base = i * p.n * p.nk;
        for (var k0 = 0u; k0 < p.nk; k0 = k0 + KT) {
            for (var l = 0u; l < 4u; l = l + 1u) {
                let idx = tid + l * 256u;
                let r = idx / KT;
                let kk = idx % KT;
                let gk = k0 + kk;
                let ga = row_base + r;
                var va = 0.0lf;
                if (ga < p.n && gk < p.nk) {
                    va = x[base + ga * p.nk + gk];
                }
                ta[r * PITCH + kk] = va;
                let gb = col_base + r;
                var vb = 0.0lf;
                if (gb < p.n && gk < p.nk) {
                    vb = x[base + gb * p.nk + gk];
                }
                tb[r * PITCH + kk] = vb;
            }
            workgroupBarrier();
            for (var kk = 0u; kk < KT; kk = kk + 1u) {
                var ra: array<f64, 4>;
                var rb: array<f64, 4>;
                for (var i2 = 0u; i2 < 4u; i2 = i2 + 1u) {
                    ra[i2] = ta[(ty * 4u + i2) * PITCH + kk];
                    rb[i2] = tb[(tx * 4u + i2) * PITCH + kk];
                }
                for (var i2 = 0u; i2 < 4u; i2 = i2 + 1u) {
                    for (var j2 = 0u; j2 < 4u; j2 = j2 + 1u) {
                        acc[i2 * 4u + j2] = acc[i2 * 4u + j2] + ra[i2] * rb[j2];
                    }
                }
            }
            workgroupBarrier();
        }
    }
    for (var i2 = 0u; i2 < 4u; i2 = i2 + 1u) {
        let ga = row_base + ty * 4u + i2;
        for (var j2 = 0u; j2 < 4u; j2 = j2 + 1u) {
            let gb = col_base + tx * 4u + j2;
            if (ga < p.n && gb < p.n) {
                kmat[ga * p.n + gb] = kmat[ga * p.n + gb] + acc[i2 * 4u + j2];
            }
        }
    }
}
