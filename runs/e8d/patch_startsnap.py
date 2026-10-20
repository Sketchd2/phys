p = r'E:\Projects\phys\src\bin\bulk.rs'
s = open(p, encoding='utf-8').read()

# 1. the rotation helper, before main
helper = '''/// The rotation taking the principal-frame vectors `(b1, b2)` to `(a1, a2)`, as a quaternion
/// (a three-atom molecule's orientation from its atoms; `PHYS_START_SNAP`).
fn rotation_between(b1: Vec3, b2: Vec3, a1: Vec3, a2: Vec3) -> Quat {
    let frame = |u1: Vec3, u2: Vec3| -> [Vec3; 3] {
        let e1 = u1.unit();
        let e3 = u1.cross(u2).unit();
        let e2 = e3.cross(e1);
        [e1, e2, e3]
    };
    let (a, b) = (frame(a1, a2), frame(b1, b2));
    let mut r = [[0.0f64; 3]; 3];
    for k in 0..3 {
        let (ak, bk) = ([a[k].x, a[k].y, a[k].z], [b[k].x, b[k].y, b[k].z]);
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] += ak[i] * bk[j];
            }
        }
    }
    let tr = r[0][0] + r[1][1] + r[2][2];
    let (w, x, y, z);
    if tr > 0.0 {
        let s = (tr + 1.0).sqrt() * 2.0;
        w = 0.25 * s;
        x = (r[2][1] - r[1][2]) / s;
        y = (r[0][2] - r[2][0]) / s;
        z = (r[1][0] - r[0][1]) / s;
    } else if r[0][0] > r[1][1] && r[0][0] > r[2][2] {
        let s = (1.0 + r[0][0] - r[1][1] - r[2][2]).sqrt() * 2.0;
        w = (r[2][1] - r[1][2]) / s;
        x = 0.25 * s;
        y = (r[0][1] + r[1][0]) / s;
        z = (r[0][2] + r[2][0]) / s;
    } else if r[1][1] > r[2][2] {
        let s = (1.0 + r[1][1] - r[0][0] - r[2][2]).sqrt() * 2.0;
        w = (r[0][2] - r[2][0]) / s;
        x = (r[0][1] + r[1][0]) / s;
        y = 0.25 * s;
        z = (r[1][2] + r[2][1]) / s;
    } else {
        let s = (1.0 + r[2][2] - r[0][0] - r[1][1]).sqrt() * 2.0;
        w = (r[1][0] - r[0][1]) / s;
        x = (r[0][2] + r[2][0]) / s;
        y = (r[1][2] + r[2][1]) / s;
        z = 0.25 * s;
    }
    Quat { w, v: Vec3 { x, y, z } }.unit()
}

fn main() {'''
assert "fn main() {" in s
s = s.replace("fn main() {", helper, 1)

# 2. after the placement loop
old = """    if unplaced > 0 {
        println!("  {unplaced} molecules could not be turned clear of their neighbours");
    }
"""
new = old + """    // `PHYS_START_SNAP=<snapshot>`: start from an earlier run's configuration (a three-atom
    // molecule's centre and orientation fitted to its atoms), the centres scaled to this box -
    // a test of whether a structure is the equilibrium one or the one a quench from a random
    // lattice selects. Unset, nothing changes.
    if let Ok(path) = std::env::var("PHYS_START_SNAP") {
        assert_eq!(atoms, 3, "PHYS_START_SNAP is for three-atom molecules");
        let text = std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("no {path}"));
        let mut old_side = 0.0;
        let mut mols: Vec<Vec<Vec3>> = Vec::new();
        for line in text.lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            match w.first().copied() {
                Some("box") => old_side = w[1].parse().unwrap(),
                Some("mol") => {
                    let v: Vec<f64> = w[2..].iter().map(|x| x.parse().unwrap()).collect();
                    mols.push(v.chunks(3).map(|c| Vec3 { x: c[0], y: c[1], z: c[2] }).collect());
                }
                _ => {}
            }
        }
        assert_eq!(mols.len(), n, "the snapshot holds {} molecules, the run {n}", mols.len());
        let sites = l.kind.sites.clone();
        let scale = l.cell[0] / old_side;
        for (i, m) in mols.iter().enumerate() {
            let q = rotation_between(sites[1] - sites[0], sites[2] - sites[0], m[1] - m[0], m[2] - m[0]);
            l.orientation[i] = q;
            l.com[i] = (m[0] - q.rotate(sites[0])).scale(scale);
        }
        println!("  started from {path}, centres scaled by {scale:.4}");
    }
"""
assert old in s
s = s.replace(old, new, 1)
open(p, 'w', encoding='utf-8').write(s)
print('ok')
