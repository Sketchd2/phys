p = r"E:\Projects\phys\src\flex.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
new = open(r"E:\Projects\phys\runs\e8c\scratch\forces_new.rs.txt").read()
a = s.index("    pub fn forces(&self, law: &dyn SiteLaw")
b = s.index("        let mut intra = 0.0;", a)
s = s[:a] + new + s[b:]
helper = """/// The gradient of the H-O-H angle with respect to the three atoms.
pub fn angle_gradient(pos: &[Vec3; 3]) -> [Vec3; 3] {
    let a = pos[1] - pos[0];
    let b = pos[2] - pos[0];
    let (r1, r2) = (a.norm(), b.norm());
    let (u, v) = (a.scale(1.0 / r1), b.scale(1.0 / r2));
    let c = u.dot(v).clamp(-1.0, 1.0);
    let sn = (1.0 - c * c).sqrt();
    let dth1 = (v - u.scale(c)).scale(-1.0 / (r1 * sn));
    let dth2 = (u - v.scale(c)).scale(-1.0 / (r2 * sn));
    [-(dth1 + dth2), dth1, dth2]
}

"""
marker = "/// The force `f` on the virtual bisector site"
assert marker in s
s = s.replace(marker, helper + marker, 1)
open(p, "w", newline="").write(s.replace("\n", nl))
print("spliced")
