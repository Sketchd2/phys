p = r'E:\Projects\phys\src\induction.rs'
s = open(p, encoding='utf-8').read()

# 1. erf_fast, next to erf
marker = "/// The factor on every smeared dipole's width"
assert marker in s
fast = '''/// `erf(x)` by a seventh-order Taylor expansion about the nearest of a table of nodes (every
/// 1/64, to 6), the table built once from [`erf`] and the Hermite recurrence for its derivatives.
/// Accurate to ~1e-16 like the series, several times faster at the large arguments a liquid's
/// pair links have (the series needs ~100 terms at `x` = 5). Used by the Ewald paths only, so no
/// result that predates it changes.
pub fn erf_fast(x: f64) -> f64 {
    const STEP: f64 = 64.0;
    const NODES: usize = 6 * 64 + 2;
    static TABLE: std::sync::OnceLock<Vec<[f64; 8]>> = std::sync::OnceLock::new();
    let table = TABLE.get_or_init(|| {
        let two_over_sqrt_pi = 2.0 / PI.sqrt();
        (0..NODES)
            .map(|i| {
                let x0 = i as f64 / STEP;
                let g = (-x0 * x0).exp();
                // a_n = d^n/dx^n e^(-x^2) = (-1)^n H_n(x0) g, by H_(n+1) = 2 x H_n - 2 n H_(n-1)
                let mut h = [0.0f64; 7];
                h[0] = 1.0;
                h[1] = 2.0 * x0;
                for n in 1..6 {
                    h[n + 1] = 2.0 * x0 * h[n] - 2.0 * n as f64 * h[n - 1];
                }
                let mut c = [0.0f64; 8];
                c[0] = erf(x0);
                let mut fact = 1.0;
                for n in 0..7 {
                    fact *= (n + 1) as f64;
                    let a_n = if n % 2 == 0 { h[n] } else { -h[n] } * g;
                    c[n + 1] = two_over_sqrt_pi * a_n / fact;
                }
                c
            })
            .collect()
    });
    let ax = x.abs();
    if ax >= 5.99 {
        return x.signum();
    }
    let i = (ax * STEP + 0.5) as usize;
    let dx = ax - i as f64 / STEP;
    let c = &table[i];
    let v = c[0] + dx * (c[1] + dx * (c[2] + dx * (c[3] + dx * (c[4] + dx * (c[5] + dx * (c[6] + dx * c[7]))))));
    if x < 0.0 {
        -v
    } else {
        v
    }
}

'''
s = s.replace(marker, fast + marker, 1)

# 2. fast copies of the three small tensor functions
marker2 = "fn tensor_apply(a: f64, bp: f64, r: Vec3, v: Vec3) -> Vec3 {"
assert marker2 in s
copies = '''/// [`dd_terms`] with [`erf_fast`] (the Ewald paths).
fn dd_terms_fast(r: f64, s: f64) -> (f64, f64, f64, f64) {
    let x = r / s;
    let g = erf_fast(x);
    let c = 2.0 / (PI.sqrt() * s) * (-x * x).exp();
    let r2 = r * r;
    let r3 = r2 * r;
    let a = g / r3 - c / r2;
    let bp = 3.0 * g / (r3 * r2) - c * (3.0 / (r2 * r2) + 2.0 / (s * s * r2));
    let da = 3.0 * c / r3 - 3.0 * g / (r2 * r2) + 2.0 * c / (s * s * r);
    let dbp = 15.0 * c / (r3 * r2) - 15.0 * g / (r3 * r3) + 10.0 * c / (s * s * r3) + 4.0 * c / (s * s * s * s * r);
    (a, bp, da, dbp)
}

/// [`charge_terms`] with [`erf_fast`] (the Ewald paths).
fn charge_terms_fast(r: f64, s_a: f64) -> (f64, f64) {
    let x = r / s_a;
    let ex = (-x * x).exp();
    let gg = erf_fast(x) - 2.0 / PI.sqrt() * x * ex;
    let dgg = 4.0 * r * r / (PI.sqrt() * s_a * s_a * s_a) * ex;
    let r3 = r * r * r;
    (gg / r3, dgg / r3 - 3.0 * gg / (r3 * r))
}

'''
s = s.replace(marker2, copies + marker2, 1)

# 3. use them inside solve_ewald only
a = s.index("pub fn solve_ewald(")
head, tail = s[:a], s[a:]
tail = tail.replace("let (a, b, c, d) = dd_terms(r, s);\n        let (a2, b2, c2, d2) = dd_terms(r, sp);", "let (a, b, c, d) = dd_terms_fast(r, s);\n        let (a2, b2, c2, d2) = dd_terms_fast(r, sp);")
tail = tail.replace("let (h, dh) = charge_terms(r, s);\n        let (h2, dh2) = charge_terms(r, sp);", "let (h, dh) = charge_terms_fast(r, s);\n        let (h2, dh2) = charge_terms_fast(r, sp);")
tail = tail.replace("let (h, dh) = charge_terms(r.norm(), sp);", "let (h, dh) = charge_terms_fast(r.norm(), sp);")
s = head + tail
open(p, 'w', encoding='utf-8').write(s)
print('ok', tail.count('_fast('))
