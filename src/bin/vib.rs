//! What the molecule's own surface says about its motion — `docs/PLAY.md` E8c.
//!
//! ```sh
//! cargo run --release --bin phys-vib -- pes-water.txt
//! ```
//!
//! Fits the surface `phys-pes` computed (a quartic polynomial in the two bond
//! displacements and the angle's, symmetric in the hydrogens by construction
//! because the data is mirrored), finds its minimum, and from it the three
//! vibrational frequencies (the mass-weighted Hessian of the surface as a
//! function of the atoms' positions, its three largest eigenvalues) to set
//! beside the measured ones, the cubic constants along the normal modes, and the
//! shape of the molecule averaged over its ground state: each normal
//! coordinate's mean moved by the cubic terms, `<Q_k> = -sum_j f_kjj / (2 w_j)
//! / (2 w_k^2)`, and its spread `1 / sqrt(2 w_k)`, sampled and the bond lengths
//! and angle taken of every sample, so that the curvature of the coordinates is in
//! the average. That averaged shape is the rigid molecule of E8c's first
//! treatment.

use phys::electrons::linalg::{eigh, Matrix};

const BOHR_PER_ANGSTROM: f64 = 1.0 / 0.529177210903;
const AMU: f64 = 1822.888486;
const HARTREE_TO_CM: f64 = 219474.6313;

fn terms() -> Vec<(usize, usize, usize)> {
    let mut t = Vec::new();
    for a in 0..=4 {
        for b in 0..=4 - a {
            for c in 0..=4 - a - b {
                t.push((a, b, c));
            }
        }
    }
    t
}

fn basis(t: &[(usize, usize, usize)], x: f64, y: f64, z: f64) -> Vec<f64> {
    // Scaled so that the normal equations stay well conditioned.
    let (x, y, z) = (x / 0.1, y / 0.1, z / 0.2);
    t.iter().map(|&(a, b, c)| x.powi(a as i32) * y.powi(b as i32) * z.powi(c as i32)).collect()
}

fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for c in 0..n {
        let piv = (c..n).max_by(|&i, &j| a[i][c].abs().partial_cmp(&a[j][c].abs()).unwrap()).unwrap();
        a.swap(c, piv);
        b.swap(c, piv);
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            for k in c..n {
                a[r][k] -= f * a[c][k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for c in (0..n).rev() {
        x[c] = (b[c] - (c + 1..n).map(|k| a[c][k] * x[k]).sum::<f64>()) / a[c][c];
    }
    x
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let file = args.first().cloned().unwrap_or_else(|| "pes-water.txt".into());
    let mut data: Vec<(f64, f64, f64, f64)> = Vec::new();
    for l in std::fs::read_to_string(&file).unwrap_or_else(|_| panic!("no {file}")).lines().filter(|l| !l.starts_with('#')) {
        let w: Vec<f64> = l.split_whitespace().filter_map(|x| x.parse().ok()).collect();
        if w.len() >= 5 {
            data.push((w[0], w[1], w[2], w[3] + w[4]));
            if (w[0] - w[1]).abs() > 1e-12 {
                data.push((w[1], w[0], w[2], w[3] + w[4]));
            }
        }
    }
    let t = terms();
    let (r0, th0) = (0.9689 * BOHR_PER_ANGSTROM, 104.21f64.to_radians());
    let e_ref = data.iter().map(|d| d.3).fold(f64::INFINITY, f64::min);
    let n = t.len();
    let mut ata = vec![vec![0.0; n]; n];
    let mut atb = vec![0.0; n];
    for &(r1, r2, th, e) in &data {
        let f = basis(&t, r1 - r0, r2 - r0, th - th0);
        for i in 0..n {
            atb[i] += f[i] * (e - e_ref);
            for j in 0..n {
                ata[i][j] += f[i] * f[j];
            }
        }
    }
    let coef = solve(ata, atb);
    let surface = |r1: f64, r2: f64, th: f64| -> f64 { basis(&t, r1 - r0, r2 - r0, th - th0).iter().zip(&coef).map(|(f, c)| f * c).sum::<f64>() + e_ref };
    let rms = (data.iter().map(|&(r1, r2, th, e)| (surface(r1, r2, th) - e).powi(2)).sum::<f64>() / data.len() as f64).sqrt();
    println!("{} points (mirrored), {} terms: fit rms {:.2e} hartree ({:.2} cm^-1)", data.len(), n, rms, rms * HARTREE_TO_CM);
    // The minimum, by Newton on the fitted surface.
    let mut q = [r0, r0, th0];
    for _ in 0..50 {
        let h = 1e-4;
        let mut g = [0.0; 3];
        let mut hess = [[0.0; 3]; 3];
        for i in 0..3 {
            let mut p = q;
            p[i] += h;
            let mut m = q;
            m[i] -= h;
            g[i] = (surface(p[0], p[1], p[2]) - surface(m[0], m[1], m[2])) / (2.0 * h);
            for j in 0..3 {
                let at = |si: f64, sj: f64| {
                    let mut x = q;
                    x[i] += si * h;
                    x[j] += sj * h;
                    surface(x[0], x[1], x[2])
                };
                hess[i][j] = (at(1.0, 1.0) - at(1.0, -1.0) - at(-1.0, 1.0) + at(-1.0, -1.0)) / (4.0 * h * h);
            }
        }
        let step = solve(hess.iter().map(|r| r.to_vec()).collect(), g.to_vec());
        for i in 0..3 {
            q[i] -= step[i];
        }
        if step.iter().map(|s| s.abs()).fold(0.0, f64::max) < 1e-12 {
            break;
        }
    }
    println!("minimum of the surface: r_OH {:.5} A, theta {:.3} degrees, E {:.8} hartree", q[0] / BOHR_PER_ANGSTROM, q[2].to_degrees(), surface(q[0], q[1], q[2]));
    // Cartesian energy and the mass-weighted Hessian.
    let masses = [15.9949146 * AMU, 1.00782503 * AMU, 1.00782503 * AMU];
    let geometry = |x: &[f64; 9]| -> (f64, f64, f64) {
        let d1 = [x[3] - x[0], x[4] - x[1], x[5] - x[2]];
        let d2 = [x[6] - x[0], x[7] - x[1], x[8] - x[2]];
        let (r1, r2) = ((d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt(), (d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]).sqrt());
        let c = (d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2]) / (r1 * r2);
        (r1, r2, c.clamp(-1.0, 1.0).acos())
    };
    let energy = |x: &[f64; 9]| -> f64 {
        let (r1, r2, th) = geometry(x);
        surface(r1, r2, th)
    };
    let xe: [f64; 9] = [0.0, 0.0, 0.0, q[0], 0.0, 0.0, q[1] * q[2].cos(), q[1] * q[2].sin(), 0.0];
    let h = 1e-3;
    let mut hm = Matrix::zeros(9);
    for i in 0..9 {
        for j in 0..9 {
            let at = |si: f64, sj: f64| {
                let mut x = xe;
                x[i] += si * h;
                x[j] += sj * h;
                energy(&x)
            };
            let v = (at(1.0, 1.0) - at(1.0, -1.0) - at(-1.0, 1.0) + at(-1.0, -1.0)) / (4.0 * h * h);
            hm.set(i, j, v / (masses[i / 3] * masses[j / 3]).sqrt());
        }
    }
    let (vals, vecs) = eigh(&hm);
    let modes: Vec<usize> = (6..9).collect();
    let omega: Vec<f64> = modes.iter().map(|&k| vals[k].max(0.0).sqrt()).collect();
    println!("harmonic frequencies from the surface: {:.0} {:.0} {:.0} cm^-1 (the measured fundamentals are 1595, 3657, 3756; the measured harmonics 1649, 3832, 3943); the smallest of the six others {:.1} cm^-1", omega[0] * HARTREE_TO_CM, omega[1] * HARTREE_TO_CM, omega[2] * HARTREE_TO_CM, vals[5].abs().sqrt() * HARTREE_TO_CM);
    // Potential along the normal coordinates Q (mass-weighted, atomic units).
    let displacement = |qv: &[f64; 3]| -> [f64; 9] {
        let mut x = xe;
        for (m, &k) in modes.iter().enumerate() {
            for i in 0..9 {
                x[i] += vecs.get(i, k) * qv[m] / masses[i / 3].sqrt();
            }
        }
        x
    };
    let v = |qv: &[f64; 3]| energy(&displacement(qv));
    let hq = 0.15;
    // f_kjj by differencing the second derivative along j across a step in k.
    let second = |center: [f64; 3], j: usize| -> f64 {
        let mut p = center;
        p[j] += hq;
        let mut m = center;
        m[j] -= hq;
        (v(&p) + v(&m) - 2.0 * v(&center)) / (hq * hq)
    };
    let mut f3 = [[0.0f64; 3]; 3]; // f3[k][j] = f_kjj
    for k in 0..3 {
        for j in 0..3 {
            let mut p = [0.0; 3];
            p[k] += hq;
            let mut m = [0.0; 3];
            m[k] -= hq;
            f3[k][j] = (second(p, j) - second(m, j)) / (2.0 * hq);
        }
    }
    let shift: Vec<f64> = (0..3).map(|k| -(0..3).map(|j| f3[k][j] / (2.0 * omega[j])).sum::<f64>() / (2.0 * omega[k] * omega[k])).collect();
    let zpe: f64 = omega.iter().map(|w| 0.5 * w).sum();
    println!("zero-point energy {:.0} cm^-1; mean shifts of the normal coordinates {:.4} {:.4} {:.4}", zpe * HARTREE_TO_CM, shift[0], shift[1], shift[2]);
    // Sample the ground state's Gaussian, moved by the shifts.
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut uniform = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let samples = 200000;
    let (mut sr1, mut sr2, mut sth, mut sinv) = (0.0, 0.0, 0.0, 0.0);
    for _ in 0..samples {
        let mut qv = [0.0; 3];
        for k in 0..3 {
            let (u1, u2) = (uniform().max(1e-300), uniform());
            let g = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
            qv[k] = shift[k] + g / (2.0 * omega[k]).sqrt();
        }
        let (r1, r2, th) = geometry(&displacement(&qv));
        sr1 += r1;
        sr2 += r2;
        sth += th;
        sinv += 0.5 * (1.0 / r1 + 1.0 / r2);
    }
    let s = samples as f64;
    let r_avg = 0.5 * (sr1 + sr2) / s;
    println!(
        "averaged over the ground state: <r_OH> {:.5} A, <theta> {:.3} degrees (the minimum {:.5} A, {:.3}); <1/r> gives {:.5} A",
        r_avg / BOHR_PER_ANGSTROM,
        sth / s * 180.0 / std::f64::consts::PI,
        q[0] / BOHR_PER_ANGSTROM,
        q[2].to_degrees(),
        s / sinv / BOHR_PER_ANGSTROM
    );
    println!("(the engine's rigid molecule so far: r_OH 0.9689 A, theta 104.21 degrees)");
}
