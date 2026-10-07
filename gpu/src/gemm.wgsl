// out = a b^T in single precision: `a` is m x kd and `b` is n x kd, rows
// contiguous (so both are read along k), and `out` holds rows row0 .. row0 + rows
// of the m x n result. A 64 x 64 tile of the result per workgroup of 16 x 16
// threads, each thread a 4 x 4 block, the operands staged through workgroup
// memory sixteen k at a time (rows padded to 17 so that a thread's four rows
// fall in different banks).

struct Params {
    m: u32,
    n: u32,
    kd: u32,
    row0: u32,
    rows: u32,
    pad0: u32,
    pad1: u32,
    pad2: u32,
};

@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read> a: array<f32>;
@group(0) @binding(2) var<storage, read> b: array<f32>;
@group(0) @binding(3) var<storage, read_write> c: array<f32>;

const TILE: u32 = 64u;
const KT: u32 = 16u;
const PITCH: u32 = 17u;

var<workgroup> ta: array<f32, 1088>;
var<workgroup> tb: array<f32, 1088>;

@compute @workgroup_size(16, 16)
fn main(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_id) li: vec3<u32>) {
    let tx = li.x;
    let ty = li.y;
    let tid = ty * 16u + tx;
    let row_base = p.row0 + wg.y * TILE;
    let col_base = wg.x * TILE;
    var acc: array<f32, 16>;
    for (var q = 0u; q < 16u; q = q + 1u) {
        acc[q] = 0.0;
    }
    for (var k0 = 0u; k0 < p.kd; k0 = k0 + KT) {
        for (var l = 0u; l < 4u; l = l + 1u) {
            let idx = tid + l * 256u;
            let r = idx / KT;
            let kk = idx % KT;
            let gk = k0 + kk;
            let gr = row_base + r;
            var va = 0.0;
            if (gr < p.m && gk < p.kd) {
                va = a[gr * p.kd + gk];
            }
            ta[r * PITCH + kk] = va;
            let gc = col_base + r;
            var vb = 0.0;
            if (gc < p.n && gk < p.kd) {
                vb = b[gc * p.kd + gk];
            }
            tb[r * PITCH + kk] = vb;
        }
        workgroupBarrier();
        for (var kk = 0u; kk < KT; kk = kk + 1u) {
            var ra: array<f32, 4>;
            var rb: array<f32, 4>;
            for (var i = 0u; i < 4u; i = i + 1u) {
                ra[i] = ta[(ty * 4u + i) * PITCH + kk];
                rb[i] = tb[(tx * 4u + i) * PITCH + kk];
            }
            for (var i = 0u; i < 4u; i = i + 1u) {
                for (var j = 0u; j < 4u; j = j + 1u) {
                    acc[i * 4u + j] = acc[i * 4u + j] + ra[i] * rb[j];
                }
            }
        }
        workgroupBarrier();
    }
    let last = min(p.m, p.row0 + p.rows);
    for (var i = 0u; i < 4u; i = i + 1u) {
        let gr = row_base + ty * 4u + i;
        for (var j = 0u; j < 4u; j = j + 1u) {
            let gc = col_base + tx * 4u + j;
            if (gr < last && gc < p.n) {
                c[(gr - p.row0) * p.n + gc] = acc[i * 4u + j];
            }
        }
    }
}
