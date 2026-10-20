p = r"E:\Projects\phys\src\pairs.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")

# 1. the field
s = s.replace("    pub flex: Option<FlexDraw>,\n}\n", """    pub flex: Option<FlexDraw>,
    /// Draw pairs at close contact (`--wall LO HI`, angstrom): the nearest two atoms of
    /// the molecules between LO and HI apart, which is the repulsive wall that the
    /// uniform and biased draws never reach (see [`Monomer::wall_pair`]).
    pub wall: Option<(f64, f64)>,
}
""", 1)
s = s.replace("bias: None, flex: None }", "bias: None, flex: None, wall: None }", 1)

# 2. the draw
marker = "    /// Pair `k` of the random draw: two molecules' atom positions"
assert marker in s
wall = '''    /// Pair `k` drawn at close contact: orientations as in [`Monomer::random_pair`], the
    /// centres 2.4 to 5.6 bohr apart, kept if the nearest atoms of the two molecules are
    /// between the `--wall` limits. The uniform draw rejects anything nearer than
    /// 0.64 of the van der Waals contact and the biased one anything the law finds
    /// repulsive, so neither puts an energy on the wall; a law fitted to them alone
    /// has no wall of its own (a scan found -577 and -649 kcal/mol at O-O 1.59 A).
    /// The pair is a function of `k` alone.
    pub fn wall_pair(&self, k: usize) -> (Vec<Vec3>, Vec<Vec3>, f64) {
        let (lo, hi) = self.wall.expect("--wall limits");
        let (lo, hi) = (lo / 0.529177210903, hi / 0.529177210903);
        let n = self.z.len();
        let mut s = Stream::at(0x7761_6c6c_5f70 ^ self.name.len() as u64, k as u128, 0, Purpose::Positions);
        loop {
            let qa = random_rotation(&mut s);
            let qb = random_rotation(&mut s);
            let dir = s.direction();
            let sep = s.range(2.4, 5.6);
            let (body_a, body_b) = (self.drawn_body(&mut s), self.drawn_body(&mut s));
            let a: Vec<Vec3> = body_a.iter().map(|p| qa.rotate(*p)).collect();
            let b: Vec<Vec3> = body_b.iter().map(|p| qb.rotate(*p) + dir.scale(sep)).collect();
            let near = (0..n).flat_map(|i| (0..n).map(move |j| (i, j))).map(|(i, j)| (a[i] - b[j]).norm()).fold(f64::INFINITY, f64::min);
            if near >= lo && near <= hi {
                return (a, b, sep);
            }
        }
    }

'''
s = s.replace(marker, wall + marker, 1)

# 3. the driver
s = s.replace("""        let (a, b, sep) = match &mono.bias {
            Some(b) => mono.biased_pair(k, b),
            None => mono.random_pair(k),""", """        let (a, b, sep) = match (&mono.bias, mono.wall) {
            (_, Some(_)) => mono.wall_pair(k),
            (Some(b), None) => mono.biased_pair(k, b),
            (None, None) => mono.random_pair(k),""", 1)
open(p, "w", newline="").write(s.replace("\n", nl))
print("pairs.rs patched")
