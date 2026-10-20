p = r'E:\Projects\phys\src\ewald.rs'
s = open(p, encoding='utf-8').read()


def replace_between(s, start_marker, end_marker, new):
    a = s.index(start_marker)
    b = s.index(end_marker)
    return s[:a] + new + s[b:]


s = s.replace("const PI: f64 = std::f64::consts::PI;", """const PI: f64 = std::f64::consts::PI;

/// How many threads the reciprocal loops use. Each thread owns a disjoint slice
/// of the output and adds in a fixed order, so the numbers do not depend on it.
fn threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).min(12)
}""", 1)

s = replace_between(s, "    pub fn phases(&self, pos: &[Vec3]) -> Vec<C> {", "    /// `S_q(k) = sum q e^(ik.r)`.", """    pub fn phases(&self, pos: &[Vec3]) -> Vec<C> {
        let nk = self.idx.len();
        let mut out = vec![C::default(); pos.len() * nk];
        if pos.is_empty() {
            return out;
        }
        let n = [self.nmax[0] as usize, self.nmax[1] as usize, self.nmax[2] as usize];
        let per = pos.len().div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in out.chunks_mut(per * nk).enumerate() {
                sc.spawn(move || {
                    for (si, row) in chunk.chunks_mut(nk).enumerate() {
                        let p = pos[ci * per + si];
                        // e^(i n theta) for n = 0..=nmax on each axis; negative n by conjugate.
                        let table = |x: f64, a: usize| -> Vec<C> {
                            let th = 2.0 * PI * x / self.cell[a];
                            let step = C { re: th.cos(), im: th.sin() };
                            let mut t = Vec::with_capacity(n[a] + 1);
                            let mut cur = C { re: 1.0, im: 0.0 };
                            for m in 0..=n[a] {
                                // Recompute exactly now and then so the error does not creep.
                                if m % 8 == 0 {
                                    let a2 = th * m as f64;
                                    cur = C { re: a2.cos(), im: a2.sin() };
                                }
                                t.push(cur);
                                cur = cur.mul(step);
                            }
                            t
                        };
                        let (tx, ty, tz) = (table(p.x, 0), table(p.y, 1), table(p.z, 2));
                        let get = |t: &Vec<C>, m: i32| if m >= 0 { t[m as usize] } else { t[(-m) as usize].conj() };
                        for (j, m) in self.idx.iter().enumerate() {
                            row[j] = get(&tx, m[0]).mul(get(&ty, m[1])).mul(get(&tz, m[2]));
                        }
                    }
                });
            }
        });
        out
    }

""")

s = replace_between(s, "    pub fn charge_sum(&self, ph: &[C], q: &[f64]) -> Vec<C> {", "    /// `D(k) = sum i (k.mu) e^(ik.r)`.", """    pub fn charge_sum(&self, ph: &[C], q: &[f64]) -> Vec<C> {
        let nk = self.idx.len();
        let mut s = vec![C::default(); nk];
        let per = nk.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in s.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    let k0 = ci * per;
                    for (i, qi) in q.iter().enumerate() {
                        let row = &ph[i * nk + k0..i * nk + k0 + chunk.len()];
                        for (sk, p) in chunk.iter_mut().zip(row) {
                            sk.re += qi * p.re;
                            sk.im += qi * p.im;
                        }
                    }
                });
            }
        });
        s
    }

""")

s = replace_between(s, "    pub fn dipole_sum(&self, ph: &[C], mu: &[Vec3]) -> Vec<C> {", "    /// `(2 pi / V) sum_(all k) f |S|^2", """    pub fn dipole_sum(&self, ph: &[C], mu: &[Vec3]) -> Vec<C> {
        let nk = self.idx.len();
        let mut s = vec![C::default(); nk];
        let per = nk.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in s.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    let k0 = ci * per;
                    for (a, m) in mu.iter().enumerate() {
                        let row = &ph[a * nk + k0..a * nk + k0 + chunk.len()];
                        for (j, (sk, p)) in chunk.iter_mut().zip(row).enumerate() {
                            let k = self.kv[k0 + j];
                            let km = k[0] * m.x + k[1] * m.y + k[2] * m.z;
                            // i (k.mu) (p.re + i p.im) = (k.mu)(-p.im + i p.re)
                            sk.re += -km * p.im;
                            sk.im += km * p.re;
                        }
                    }
                });
            }
        });
        s
    }

""")

s = replace_between(s, "    pub fn grad_dipole(&self, s: &[C], ph: &[C], n: usize) -> Vec<Vec3> {", "    /// `dE/dr_i` of the charges", """    pub fn grad_dipole(&self, s: &[C], ph: &[C], n: usize) -> Vec<Vec3> {
        let nk = self.idx.len();
        let pref = 8.0 * PI / self.volume;
        let mut out = vec![Vec3::ZERO; n];
        if n == 0 {
            return out;
        }
        let per = n.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in out.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    for (ai, o) in chunk.iter_mut().enumerate() {
                        let a = ci * per + ai;
                        let row = &ph[a * nk..(a + 1) * nk];
                        let (mut gx, mut gy, mut gz) = (0.0, 0.0, 0.0);
                        for j in 0..nk {
                            // S* i e^(ikr): S* = (sr - i si), i e = i(pr + i pi) = (-pi + i pr)
                            let (sr, si) = (s[j].re, s[j].im);
                            let (pr, pi) = (row[j].re, row[j].im);
                            let re = sr * (-pi) + si * pr; // Re[(sr - i si)(-pi + i pr)]
                            let w = self.fk[j] * re;
                            gx += w * self.kv[j][0];
                            gy += w * self.kv[j][1];
                            gz += w * self.kv[j][2];
                        }
                        *o = Vec3 { x: gx * pref, y: gy * pref, z: gz * pref };
                    }
                });
            }
        });
        out
    }

""")

s = replace_between(s, "    pub fn grad_dipole_pos(&self, s: &[C], ph: &[C], mu: &[Vec3]) -> Vec<Vec3> {", "    /// The long-range field of the permanent charges", """    pub fn grad_dipole_pos(&self, s: &[C], ph: &[C], mu: &[Vec3]) -> Vec<Vec3> {
        let nk = self.idx.len();
        let pref = 8.0 * PI / self.volume;
        let mut out = vec![Vec3::ZERO; mu.len()];
        if mu.is_empty() {
            return out;
        }
        let per = mu.len().div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in out.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    for (ai, o) in chunk.iter_mut().enumerate() {
                        let a = ci * per + ai;
                        let m = mu[a];
                        let row = &ph[a * nk..(a + 1) * nk];
                        let (mut gx, mut gy, mut gz) = (0.0, 0.0, 0.0);
                        for j in 0..nk {
                            let k = self.kv[j];
                            let km = k[0] * m.x + k[1] * m.y + k[2] * m.z;
                            // Re[ -(sr - i si)(pr + i pi) ] = -(sr pr + si pi)
                            let re = -(s[j].re * row[j].re + s[j].im * row[j].im);
                            let w = self.fk[j] * re * km;
                            gx += w * k[0];
                            gy += w * k[1];
                            gz += w * k[2];
                        }
                        *o = Vec3 { x: gx * pref, y: gy * pref, z: gz * pref };
                    }
                });
            }
        });
        out
    }

""")
open(p, 'w', encoding='utf-8').write(s)
print('patched')
