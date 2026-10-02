//! A spherical atom, solved radially.
//!
//! A free atom with its open shell shared equally among the degenerate levels
//! is spherical, and then every orbital is a radial function times a spherical
//! harmonic. The Kohn-Sham problem falls apart into one small matrix per
//! angular momentum, on radial functions `r^l exp(-a r^2)`; the Coulomb
//! potential is Gauss's law on shells; exchange-correlation sees `rho(r)` and
//! `rho'(r)`. It is the same physics as [`super::scf`] in three dimensions —
//! `tests/electrons.rs` holds the two to each other — at a cost small enough
//! that a basis can be derived by trying hundreds of them.
//!
//! # The Coulomb potential to round-off
//!
//! Gauss's law needs the charge inside each radius, a cumulative integral, and
//! a running sum on the quadrature points is only second-order accurate: for
//! neon's 66 hartree of Coulomb energy that is an error of order 1e-5, too
//! coarse to derive a basis against. In the radial map's angle variable `t`
//! (`x = cos t`), the integrand vanishes at both ends — as `(pi - t)^5` at the
//! nucleus and exponentially far out — so it is a smooth odd function, its sine
//! series converges spectrally, and the series integrates term by term.

use super::functional::{evaluate, Functional};
use super::linalg::{generalised, orthogonaliser, Matrix};

/// The answer for one atom.
#[derive(Debug, Clone)]
pub struct Atom {
    pub energy: f64,
    pub converged: bool,
    pub iterations: usize,
    /// `(l, spin, level, occupation per orbital)`.
    pub levels: Vec<(usize, usize, f64, f64)>,
    /// The converged effective potential, per spin, on the radial grid:
    /// `(Hartree + exchange-correlation, gradient term)`.
    potential: [(Vec<f64>, Vec<f64>); 2],
    /// Occupied orbitals: `(l, spin, level, electrons, radial values)`.
    orbitals: Vec<(usize, usize, f64, f64, Vec<f64>)>,
    /// The same orbitals' coefficients over the radial functions solved in:
    /// `(l, spin, electrons, coefficients)`.
    pub coefficients: Vec<(usize, usize, f64, Vec<f64>)>,
    z: f64,
}

const PI: f64 = std::f64::consts::PI;

/// The radial grid: points `r`, angle variables `t`, and weights for
/// `integral f(r) r^2 dr`.
struct Radial {
    r: Vec<f64>,
    t: Vec<f64>,
    w: Vec<f64>,
    /// `dr/dt`, for the spectral integrals.
    drdt: Vec<f64>,
}

fn radial_grid(n: usize, scale: f64) -> Radial {
    let mut g = Radial { r: Vec::new(), t: Vec::new(), w: Vec::new(), drdt: Vec::new() };
    let h = PI / (n + 1) as f64;
    for i in 1..=n {
        let t = i as f64 * h;
        let x = t.cos();
        let r = scale * (1.0 + x) / (1.0 - x);
        // dr/dt = dr/dx * dx/dt = 2 scale / (1 - x)^2 * (-sin t); integrate in t
        // over (0, pi) with r decreasing, so |dr/dt|.
        let drdt = 2.0 * scale / ((1.0 - x) * (1.0 - x)) * t.sin();
        g.r.push(r);
        g.t.push(t);
        g.drdt.push(drdt);
        g.w.push(h * drdt * r * r);
    }
    g
}

/// `G(t_i) = integral_0^{t_i} g(t) dt` for `g` sampled at `t_i = i pi/(n+1)`
/// and vanishing at both ends, through its sine series. `table` holds
/// `sin(k t_i)` and `(1 - cos(k t_i)) / k`, which depend only on `n`.
fn cumulative(g: &[f64], table: &Series) -> Vec<f64> {
    let n = g.len();
    let mut out = vec![0.0; n];
    for k in 0..n {
        let row = &table.sin[k * n..k * n + n];
        let b: f64 = g.iter().zip(row).map(|(a, b)| a * b).sum::<f64>() * 2.0 / (n + 1) as f64;
        let cum = &table.cum[k * n..k * n + n];
        for (o, c) in out.iter_mut().zip(cum) {
            *o += b * c;
        }
    }
    out
}

/// The sine-series tables for `n` points.
struct Series {
    sin: Vec<f64>,
    cum: Vec<f64>,
}

impl Series {
    fn new(n: usize) -> Series {
        let h = PI / (n + 1) as f64;
        let mut sin = vec![0.0; n * n];
        let mut cum = vec![0.0; n * n];
        for k in 0..n {
            let kf = (k + 1) as f64;
            for i in 0..n {
                let a = kf * (i + 1) as f64 * h;
                sin[k * n + i] = a.sin();
                cum[k * n + i] = (1.0 - a.cos()) / kf;
            }
        }
        Series { sin, cum }
    }
}

/// The Hartree potential of a spherical density, `4 pi [Q(<r)/r + integral_r^inf rho r' dr']`.
fn hartree(g: &Radial, table: &Series, rho: &[f64]) -> Vec<f64> {
    let n = g.r.len();
    // t runs from r = infinity (t -> 0) to r = 0 (t -> pi).
    let q: Vec<f64> = (0..n).map(|i| 4.0 * PI * rho[i] * g.r[i] * g.r[i] * g.drdt[i]).collect();
    let p: Vec<f64> = (0..n).map(|i| 4.0 * PI * rho[i] * g.r[i] * g.drdt[i]).collect();
    let cq = cumulative(&q, table); // charge outside r_i
    let cp = cumulative(&p, table); // integral of rho r' from r_i to infinity
    // Charge inside r_i is the total less what lies outside. The total is the
    // integral over the whole of (0, pi), where the plain sum is spectrally
    // accurate for a smooth function vanishing at both ends.
    let total: f64 = q.iter().sum::<f64>() * PI / (n + 1) as f64;
    (0..n).map(|i| (total - cq[i]) / g.r[i] + cp[i]).collect()
}

/// Normalisation of `r^l exp(-a r^2)` under `integral (.)^2 r^2 dr = 1`.
fn radial_norm(a: f64, l: usize) -> f64 {
    let mut df = 1.0;
    let mut k = 2 * l as i64 + 1;
    while k > 1 {
        df *= k as f64;
        k -= 2;
    }
    let p = 2.0 * a;
    let integral = df / (2f64.powi(l as i32 + 2) * p.powi(l as i32 + 1)) * (PI / p).sqrt();
    1.0 / integral.sqrt()
}

/// Radial points used. Validated against the three-dimensional solver.
pub const RADIAL_POINTS: usize = 600;

/// Solve the atom of charge `z` with `alpha` and `beta` electrons, in the
/// radial functions with exponents `exponents[l]`.
pub fn solve(z: f64, alpha: f64, beta: f64, exponents: &[Vec<f64>], f: Functional, max_iter: usize, tol: f64) -> Atom {
    let functions: Vec<Vec<Vec<(f64, f64)>>> = exponents.iter().map(|ex| ex.iter().map(|&a| vec![(a, 1.0)]).collect()).collect();
    solve_functions(z, alpha, beta, &functions, f, max_iter, tol)
}

/// As [`solve`], with each radial function of angular momentum `l` any
/// combination `sum c_k chi_k` of normalised primitives, given as
/// `functions[l][i] = [(exponent, c), ...]` — a contracted basis.
pub fn solve_functions(z: f64, alpha: f64, beta: f64, functions: &[Vec<Vec<(f64, f64)>>], f: Functional, max_iter: usize, tol: f64) -> Atom {
    let g = radial_grid(RADIAL_POINTS, 1.0);
    let n = g.r.len();
    let table = Series::new(n);
    let lmax = functions.len();
    let mut chi: Vec<Vec<Vec<f64>>> = Vec::new();
    let mut dchi: Vec<Vec<Vec<f64>>> = Vec::new();
    for (l, funcs) in functions.iter().enumerate() {
        let li = l as i32;
        chi.push(funcs.iter().map(|terms| {
            g.r.iter().map(|&x| terms.iter().map(|&(a, c)| c * radial_norm(a, l) * x.powi(li) * (-a * x * x).exp()).sum()).collect()
        }).collect());
        dchi.push(funcs.iter().map(|terms| {
            g.r.iter().map(|&x| {
                let lead = if l > 0 { l as f64 * x.powi(li - 1) } else { 0.0 };
                terms.iter().map(|&(a, c)| c * radial_norm(a, l) * (-a * x * x).exp() * (lead - 2.0 * a * x.powi(li + 1))).sum()
            }).collect()
        }).collect());
    }
    // Kinetic + nuclear per l, and the orthogonaliser of each overlap.
    let mut h_l = Vec::new();
    let mut x_l = Vec::new();
    for l in 0..lmax {
        let m = chi[l].len();
        let mut s = Matrix::zeros(m);
        let mut h = Matrix::zeros(m);
        let cent = (l * (l + 1)) as f64;
        for i in 0..m {
            for j in 0..=i {
                let (mut sv, mut hv) = (0.0, 0.0);
                for p in 0..n {
                    let ab = chi[l][i][p] * chi[l][j][p];
                    sv += g.w[p] * ab;
                    hv += g.w[p] * (0.5 * (dchi[l][i][p] * dchi[l][j][p] + cent / (g.r[p] * g.r[p]) * ab) - z / g.r[p] * ab);
                }
                s.set(i, j, sv);
                s.set(j, i, sv);
                h.set(i, j, hv);
                h.set(j, i, hv);
            }
        }
        x_l.push(orthogonaliser(&s, 1e-10));
        h_l.push(h);
    }
    let counts = [alpha, beta];
    let mut rho = [vec![0.0; n], vec![0.0; n]];
    let mut drho = [vec![0.0; n], vec![0.0; n]];
    let mut levels = Vec::new();
    let mut orbitals: Vec<(usize, usize, f64, f64, Vec<f64>)> = Vec::new();
    let mut coefficients: Vec<(usize, usize, f64, Vec<f64>)> = Vec::new();
    let mut potential = [(vec![0.0; n], vec![0.0; n]), (vec![0.0; n], vec![0.0; n])];
    let mut energy = 0.0;
    let mut last = f64::INFINITY;
    let mut converged = false;
    let mut iterations = 0;
    for it in 0..max_iter {
        iterations = it + 1;
        let total: Vec<f64> = (0..n).map(|p| rho[0][p] + rho[1][p]).collect();
        let vh = hartree(&g, &table, &total);
        let mut vx = [vec![0.0; n], vec![0.0; n]];
        let mut wx = [vec![0.0; n], vec![0.0; n]];
        for p in 0..n {
            let (ga, gb) = (drho[0][p], drho[1][p]);
            let (_, de) = evaluate(f, rho[0][p], rho[1][p], ga * ga, ga * gb, gb * gb);
            vx[0][p] = de[0];
            vx[1][p] = de[1];
            wx[0][p] = 2.0 * de[2] * ga + de[3] * gb;
            wx[1][p] = 2.0 * de[4] * gb + de[3] * ga;
        }
        // Levels of every l and spin.
        let mut all: Vec<(usize, usize, f64, Vec<f64>)> = Vec::new();
        for spin in 0..2 {
            for l in 0..lmax {
                let m = chi[l].len();
                let mut fm = h_l[l].clone();
                for i in 0..m {
                    for j in 0..=i {
                        let mut v = 0.0;
                        for p in 0..n {
                            let ab = chi[l][i][p] * chi[l][j][p];
                            let dab = dchi[l][i][p] * chi[l][j][p] + chi[l][i][p] * dchi[l][j][p];
                            v += g.w[p] * ((vh[p] + vx[spin][p]) * ab + wx[spin][p] * dab);
                        }
                        fm.add(i, j, v);
                        if i != j {
                            fm.add(j, i, v);
                        }
                    }
                }
                let (x, mm) = &x_l[l];
                let (e, c) = generalised(&fm, x, *mm);
                for k in 0..*mm {
                    all.push((l, spin, e[k], (0..m).map(|i| c[i * mm + k]).collect()));
                }
            }
        }
        // Occupy: each level holds 2l+1 orbitals, shared equally when not full.
        let mut new_rho = [vec![0.0; n], vec![0.0; n]];
        let mut new_drho = [vec![0.0; n], vec![0.0; n]];
        let mut one_e = 0.0;
        levels.clear();
        orbitals.clear();
        coefficients.clear();
        for spin in 0..2 {
            potential[spin] = ((0..n).map(|p| vh[p] + vx[spin][p]).collect(), wx[spin].clone());
        }
        for spin in 0..2 {
            let mut idx: Vec<usize> = (0..all.len()).filter(|&k| all[k].1 == spin).collect();
            idx.sort_by(|&a, &b| all[a].2.total_cmp(&all[b].2));
            let mut left = counts[spin];
            for &k in &idx {
                let (l, _, e, ref c) = all[k];
                let deg = (2 * l + 1) as f64;
                let take = left.clamp(0.0, deg);
                left -= take;
                levels.push((l, spin, e, take / deg));
                if take <= 0.0 {
                    continue;
                }
                let m = c.len();
                orbitals.push((l, spin, e, take, (0..n).map(|p| (0..m).map(|i| c[i] * chi[l][i][p]).sum()).collect()));
                coefficients.push((l, spin, take, c.clone()));
                let mut hc = 0.0;
                for i in 0..m {
                    for j in 0..m {
                        hc += c[i] * c[j] * h_l[l].get(i, j);
                    }
                }
                one_e += take * hc;
                for p in 0..n {
                    let (mut rv, mut dv) = (0.0, 0.0);
                    for i in 0..m {
                        rv += c[i] * chi[l][i][p];
                        dv += c[i] * dchi[l][i][p];
                    }
                    new_rho[spin][p] += take * rv * rv / (4.0 * PI);
                    new_drho[spin][p] += take * 2.0 * rv * dv / (4.0 * PI);
                }
            }
        }
        // The energy of the new orbitals' density, every term evaluated on it.
        let tot_new: Vec<f64> = (0..n).map(|p| new_rho[0][p] + new_rho[1][p]).collect();
        let vh_new = hartree(&g, &table, &tot_new);
        let mut eh = 0.0;
        let mut exc = 0.0;
        for p in 0..n {
            eh += 0.5 * 4.0 * PI * g.w[p] * tot_new[p] * vh_new[p];
            let (ga, gb) = (new_drho[0][p], new_drho[1][p]);
            exc += 4.0 * PI * g.w[p] * evaluate(f, new_rho[0][p], new_rho[1][p], ga * ga, ga * gb, gb * gb).0;
        }
        energy = one_e + eh + exc;
        let change = (energy - last).abs();
        last = energy;
        // Damped density mixing; the first step takes the new density whole.
        let mix = if it == 0 { 1.0 } else { 0.4 };
        for s in 0..2 {
            for p in 0..n {
                rho[s][p] = mix * new_rho[s][p] + (1.0 - mix) * rho[s][p];
                drho[s][p] = mix * new_drho[s][p] + (1.0 - mix) * drho[s][p];
            }
        }
        if it > 3 && change < tol {
            converged = true;
            break;
        }
    }
    Atom { energy, converged, iterations, levels, potential, orbitals, coefficients, z }
}

/// How strongly the atom's occupied orbitals of angular momentum
/// `channel - order` respond to a weak field of multipole `order` (1 a uniform
/// field, 2 a field gradient), through functions of angular momentum `channel`
/// with these exponents: the Hylleraas second-order energy, summed over those
/// orbitals and weighted by their electrons.
///
/// **Variational in the functions**: a larger set can only make it more
/// negative, so the polarisation functions a molecule needs can be derived the
/// way the occupied ones are — extended while the response still grows — from
/// the free atom alone. Uncoupled (the potential is held at the ground
/// state's), which is enough to tell which functions are needed. Arbitrary
/// units: the angular factor of the dipole operator is common to every term.
pub fn response(atom: &Atom, channel: usize, order: usize, exponents: &[f64]) -> f64 {
    let functions: Vec<Vec<(f64, f64)>> = exponents.iter().map(|&a| vec![(a, 1.0)]).collect();
    response_functions(atom, channel, order, &functions).0
}

/// As [`response`], over contracted functions `[(exponent, c)]`, also giving
/// the first-order response orbital as coefficients over those functions,
/// averaged over the source orbitals by their electrons.
pub fn response_functions(atom: &Atom, channel: usize, order: usize, functions: &[Vec<(f64, f64)>]) -> (f64, Vec<f64>) {
    let g = radial_grid(RADIAL_POINTS, 1.0);
    let n = g.r.len();
    let l = channel;
    let m = functions.len();
    let li = l as i32;
    let chi: Vec<Vec<f64>> = functions.iter().map(|terms| {
        g.r.iter().map(|&x| terms.iter().map(|&(a, c)| c * radial_norm(a, l) * x.powi(li) * (-a * x * x).exp()).sum()).collect()
    }).collect();
    let dchi: Vec<Vec<f64>> = functions.iter().map(|terms| {
        g.r.iter().map(|&x| {
            let lead = if l > 0 { l as f64 * x.powi(li - 1) } else { 0.0 };
            terms.iter().map(|&(a, c)| c * radial_norm(a, l) * (-a * x * x).exp() * (lead - 2.0 * a * x.powi(li + 1))).sum()
        }).collect()
    }).collect();
    let cent = (l * (l + 1)) as f64;
    let mut s = Matrix::zeros(m);
    for i in 0..m {
        for j in 0..=i {
            let v: f64 = (0..n).map(|p| g.w[p] * chi[i][p] * chi[j][p]).sum();
            s.set(i, j, v);
            s.set(j, i, v);
        }
    }
    let (x, mm) = orthogonaliser(&s, 1e-10);
    let mut total = 0.0;
    let mut orbital = vec![0.0; m];
    let mut weight = 0.0;
    for spin in 0..2 {
        let (v, w) = &atom.potential[spin];
        let mut f = Matrix::zeros(m);
        for i in 0..m {
            for j in 0..=i {
                let mut e = 0.0;
                for p in 0..n {
                    let ab = chi[i][p] * chi[j][p];
                    let dab = dchi[i][p] * chi[j][p] + chi[i][p] * dchi[j][p];
                    e += g.w[p] * (0.5 * (dchi[i][p] * dchi[j][p] + cent / (g.r[p] * g.r[p]) * ab) - atom.z / g.r[p] * ab + v[p] * ab + w[p] * dab);
                }
                f.set(i, j, e);
                f.set(j, i, e);
            }
        }
        let (levels, c) = generalised(&f, &x, mm);
        for (lo, sp, eps, electrons, radial) in &atom.orbitals {
            if *sp != spin || lo + order != l {
                continue;
            }
            // <chi_k | r^order | orbital>
            let b: Vec<f64> = (0..m).map(|k| (0..n).map(|p| g.w[p] * chi[k][p] * g.r[p].powi(order as i32) * radial[p]).sum()).collect();
            let mut psi = vec![0.0; m];
            for j in 0..mm {
                let proj: f64 = (0..m).map(|k| c[k * mm + j] * b[k]).sum();
                let gap = levels[j] - eps;
                if gap > 0.0 {
                    total -= electrons * proj * proj / gap;
                    for k in 0..m {
                        psi[k] -= c[k * mm + j] * proj / gap;
                    }
                }
            }
            // Accumulate with a consistent sign, weighted by electrons.
            let sign = if orbital.iter().zip(&psi).map(|(a, b)| a * b).sum::<f64>() < 0.0 { -1.0 } else { 1.0 };
            for k in 0..m {
                orbital[k] += sign * electrons * psi[k];
            }
            weight += electrons;
        }
    }
    if weight > 0.0 {
        for o in &mut orbital {
            *o /= weight;
        }
    }
    (total, orbital)
}
