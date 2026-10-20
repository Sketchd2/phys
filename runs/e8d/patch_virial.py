p = r'E:\Projects\phys\src\ewald.rs'
s = open(p, encoding='utf-8').read()

marker = "    /// The charges' self energy in the reciprocal sum, to subtract"
assert marker in s
new = '''    /// `sum_i w_i (k.a_i) e^(ik.r_i)` for every `k` (a site's weight and its arm from its
    /// molecule's centre): what the derivative of a site's phase under rigid scaling needs.
    fn arm_sum(&self, ph: &[C], w: &[f64], arm: &[Vec3]) -> Vec<C> {
        let nk = self.idx.len();
        let mut s = vec![C::default(); nk];
        let per = nk.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in s.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    let k0 = ci * per;
                    for (i, (wi, a)) in w.iter().zip(arm).enumerate() {
                        let row = &ph[i * nk + k0..i * nk + k0 + chunk.len()];
                        for (j, (sk, p)) in chunk.iter_mut().zip(row).enumerate() {
                            let k = self.kv[k0 + j];
                            let c = wi * (k[0] * a.x + k[1] * a.y + k[2] * a.z);
                            sk.re += c * p.re;
                            sk.im += c * p.im;
                        }
                    }
                });
            }
        });
        s
    }

    /// `sum_a (k.mu_a)(k.b_a) e^(ik.r_a)` for every `k`.
    fn dipole_arm_sum(&self, ph: &[C], mu: &[Vec3], arm: &[Vec3]) -> Vec<C> {
        let nk = self.idx.len();
        let mut s = vec![C::default(); nk];
        let per = nk.div_ceil(threads());
        std::thread::scope(|sc| {
            for (ci, chunk) in s.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    let k0 = ci * per;
                    for (a, (m, b)) in mu.iter().zip(arm).enumerate() {
                        let row = &ph[a * nk + k0..a * nk + k0 + chunk.len()];
                        for (j, (sk, p)) in chunk.iter_mut().zip(row).enumerate() {
                            let k = self.kv[k0 + j];
                            let c = (k[0] * m.x + k[1] * m.y + k[2] * m.z) * (k[0] * b.x + k[1] * b.y + k[2] * b.z);
                            sk.re += c * p.re;
                            sk.im += c * p.im;
                        }
                    }
                });
            }
        });
        s
    }

    /// The reciprocal part of the molecular virial, `-dE_rec/d ln(lambda)` under rigid scaling,
    /// analytically, from the structure factor already in hand. With `k -> k / lambda`,
    /// `V -> V lambda^3` and a site's phase `k.R + (k / lambda).a`:
    ///
    /// ```text
    /// dE/dlnL = -3E + (4 pi / V) sum f (k^2 s^2 / 2 + 2) |S|^2 + (8 pi / V) sum f Re[ S* dS/dlnL ],
    /// dS/dlnL = -i sum q (k.a) e^(ik.r) - D + sum (k.mu)(k.b) e^(ik.r),
    /// ```
    ///
    /// `D` the dipoles' own sum. Agrees with [`Ewald::virial`] (the central difference) to
    /// rounding; it costs one pass over the sites where that costs two full sets of phases.
    pub fn virial_analytic(&self, s: &[C], ph_charge: &[C], q: &[f64], carm: &[Vec3], ph_dip: &[C], mu: &[Vec3], darm: &[Vec3], d: &[C]) -> f64 {
        let tc = self.arm_sum(ph_charge, q, carm);
        let td = if mu.is_empty() { vec![C::default(); self.idx.len()] } else { self.dipole_arm_sum(ph_dip, mu, darm) };
        let mut e = 0.0;
        let mut damp = 0.0;
        let mut cross = 0.0;
        for j in 0..self.idx.len() {
            let k = self.kv[j];
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            let sj = s[j];
            let s2 = sj.re * sj.re + sj.im * sj.im;
            e += self.fk[j] * s2;
            damp += self.fk[j] * (0.5 * k2 * self.s * self.s + 2.0) * s2;
            // dS/dlnL = -i tc - d + td ; -i (a + ib) = b - i a
            let ds = C { re: tc[j].im - if d.is_empty() { 0.0 } else { d[j].re } + td[j].re, im: -tc[j].re - if d.is_empty() { 0.0 } else { d[j].im } + td[j].im };
            // Re[S* dS] = sr dsr + si dsi
            cross += self.fk[j] * (sj.re * ds.re + sj.im * ds.im);
        }
        let pref = 4.0 * PI / self.volume;
        let de = -3.0 * pref * e + pref * damp + 2.0 * pref * cross;
        -de
    }

'''
s = s.replace(marker, new + marker, 1)
open(p, 'w', encoding='utf-8').write(s)

# liquid.rs: use it, with the numeric one kept behind an environment switch as the oracle
p = r'E:\Projects\phys\src\liquid.rs'
t = open(p, encoding='utf-8').read()
old = "            virial_rec = ew.virial(&ccom, &carm, &cq, &dcom, &darm, &mu);"
new = """            virial_rec = if std::env::var_os("PHYS_EWALD_NUMERIC_VIRIAL").is_some() { ew.virial(&ccom, &carm, &cq, &dcom, &darm, &mu) } else { ew.virial_analytic(&stot, &ph_c, &cq, &carm, &ph_d, &mu, &darm, &d) };"""
assert old in t
t = t.replace(old, new, 1)
old = "            virial_rec = ew.virial(&ccom, &carm, &cq, &[], &[], &[]);"
new = """            virial_rec = if std::env::var_os("PHYS_EWALD_NUMERIC_VIRIAL").is_some() { ew.virial(&ccom, &carm, &cq, &[], &[], &[]) } else { ew.virial_analytic(&sq, &ph_c, &cq, &carm, &[], &[], &[], &[]) };"""
assert old in t
t = t.replace(old, new, 1)
# the induction solve's tolerance, 1e-10 unless asked
old = "let out = crate::induction::solve_ewald(&cl, &links, ew.s, &e0_lr, &response, warm, 1e-10);"
new = "let out = crate::induction::solve_ewald(&cl, &links, ew.s, &e0_lr, &response, warm, ind_tolerance());"
assert old in t
t = t.replace(old, new, 1)
old = "/// A law with its Coulomb term cut down to the short-range remainder"
new = """/// The relative tolerance of the induced-dipole solve in an Ewald run: 1e-10 unless
/// `PHYS_IND_TOL` says otherwise. Forces at fixed dipoles are linear in the dipoles, so
/// 1e-7 is ample for dynamics and saves iterations; the finite-difference checks want the
/// default.
fn ind_tolerance() -> f64 {
    static TOL: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *TOL.get_or_init(|| std::env::var("PHYS_IND_TOL").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-10))
}

/// A law with its Coulomb term cut down to the short-range remainder"""
assert old in t
t = t.replace(old, new, 1)
open(p, 'w', encoding='utf-8').write(t)
print('ok')
