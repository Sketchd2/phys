//! Ball-and-stick water: flexible molecules on springs, every atom pair summed
//! by brute force, stepped to a steady state (`docs/PLAY.md` E8d, the owner's
//! idea for a simpler method).
//!
//! ```sh
//! phys-sticks params.txt [temperature K] [molecules] [g/cm3] [ps] [seed]
//! ```
//!
//! `params.txt` is written by `runs/e8d/sticks/fit_sticks.py`: the springs are the
//! Hessian of the engine's own monomer surface (HF + RI-MP2 over `r1, r2, theta`),
//! the atom-atom law is point charges, an exponential repulsion and a damped C6
//! per atom-pair type, fitted to the engine's MP2 pair energies, and nothing else:
//! no clouds, no induced dipoles, no extra sites. Charges are summed with Ewald
//! (`ewald.rs`), the same reciprocal machinery the rigid liquid uses, because a
//! cut-off was found to cost the pressure 1-2 kbar (E8d). Atoms are not wrapped
//! (a molecule stays whole); intermolecular distances use the nearest image. The
//! pressure is the atomic virial: (2 K + sum r.F) / 3V over every force, the
//! springs and the reciprocal sum included.

use phys::ewald::Ewald;
use phys::liquid::{switch, tang_toennies, BAR_PER_HARTREE_PER_BOHR3, K_B};
use phys::math::Vec3;
use phys::rng::{Purpose, Stream};
use std::io::Write;
use std::time::Instant;

const PS: f64 = 2.418884326585747e-5; // one atomic unit of time in ps
const KCAL: f64 = 627.509474;
const AMU: f64 = 1822.888486209; // electron masses per amu
const PI: f64 = std::f64::consts::PI;

struct Params {
    qh: f64,
    a: [f64; 3],
    b: [f64; 3],
    c6: [f64; 3],
    r0: f64,
    th0: f64,
    h: [[f64; 3]; 3],
}

fn read_params(path: &str) -> Params {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| panic!("no {path}"));
    let get = |key: &str| -> Vec<f64> { text.lines().find(|l| l.starts_with(key) && l.split_whitespace().next() == Some(key)).unwrap_or_else(|| panic!("no {key} in {path}")).split_whitespace().skip(1).map(|x| x.parse().unwrap()).collect() };
    let h = get("hessian");
    Params { qh: get("qh")[0], a: get("A").try_into().unwrap(), b: get("B").try_into().unwrap(), c6: get("C6").try_into().unwrap(), r0: get("r0")[0], th0: get("theta0")[0], h: [[h[0], h[1], h[2]], [h[3], h[4], h[5]], [h[6], h[7], h[8]]] }
}

struct World {
    n: usize, // molecules
    cell: f64,
    pos: Vec<Vec3>, // 3 per molecule: O, H1, H2
    mom: Vec<Vec3>,
    mass: [f64; 3],
}

struct Out {
    force: Vec<Vec3>,
    coulomb: f64,
    repulsion: f64,
    dispersion: f64,
    spring: f64,
    virial: f64,
}

fn threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).min(12)
}

fn erfc(x: f64) -> f64 {
    1.0 - phys::induction::erf(x)
}

fn forces(w: &World, p: &Params, ew: &Ewald, r_cut: f64, r_on: f64) -> Out {
    let na = 3 * w.n;
    let s = ew.s;
    let q = |i: usize| if i % 3 == 0 { -2.0 * p.qh } else { p.qh };
    let ty = |i: usize| usize::from(i % 3 != 0); // 0 = O, 1 = H
    let image = |mut d: Vec3| {
        d.x -= w.cell * (d.x / w.cell).round();
        d.y -= w.cell * (d.y / w.cell).round();
        d.z -= w.cell * (d.z / w.cell).round();
        d
    };
    // Intermolecular atom pairs, brute force, interleaved over threads.
    let nt = threads();
    let parts: Vec<(Vec<Vec3>, f64, f64, f64, f64)> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..nt)
            .map(|t| {
                sc.spawn(move || {
                    let mut f = vec![Vec3::ZERO; na];
                    let (mut ec, mut er, mut ed, mut wv) = (0.0, 0.0, 0.0, 0.0);
                    let mut i = t;
                    while i < na {
                        for j in i + 1..na {
                            if i / 3 == j / 3 {
                                continue;
                            }
                            let d = image(w.pos[i] - w.pos[j]);
                            let r2 = d.norm2();
                            if r2 >= r_cut * r_cut {
                                continue;
                            }
                            let r = r2.sqrt();
                            let kind = ty(i) + ty(j);
                            // Coulomb, real-space part.
                            let qq = q(i) * q(j);
                            let x = r / s;
                            let ex = (-x * x).exp();
                            let ec_ij = qq * erfc(x) / r;
                            let dec = qq * (-2.0 / (PI.sqrt() * s) * ex / r - erfc(x) / (r * r));
                            // Repulsion and damped dispersion, switched.
                            let bk = p.b[kind];
                            let rep = p.a[kind] * (-bk * r).exp();
                            let (f6, df6) = tang_toennies(6, bk * r);
                            let r6 = r2 * r2 * r2;
                            let disp = f6 * p.c6[kind] / r6;
                            let (sw, dsw) = switch(r, r_on, r_cut);
                            let drep = -bk * rep;
                            let ddisp = p.c6[kind] * (bk * df6 / r6 - 6.0 * f6 / (r6 * r));
                            let du = dec + sw * (drep - ddisp) + dsw * (rep - disp);
                            ec += ec_ij;
                            er += sw * rep;
                            ed += sw * disp;
                            // Force on i is -dU/dr along d/r; virial r.F = -r dU/dr.
                            let fi = d.scale(-du / r);
                            f[i] += fi;
                            f[j] -= fi;
                            wv += -r * du;
                        }
                        i += nt;
                    }
                    (f, ec, er, ed, wv)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut out = Out { force: vec![Vec3::ZERO; na], coulomb: 0.0, repulsion: 0.0, dispersion: 0.0, spring: 0.0, virial: 0.0 };
    for (f, ec, er, ed, wv) in parts {
        for (a, b) in out.force.iter_mut().zip(f) {
            *a += b;
        }
        out.coulomb += ec;
        out.repulsion += er;
        out.dispersion += ed;
        out.virial += wv;
    }
    // Inside each molecule: the long-range Coulomb the reciprocal sum includes is taken off
    // (the molecule's own charges do not act on each other), and the springs act.
    for m in 0..w.n {
        let (io, i1, i2) = (3 * m, 3 * m + 1, 3 * m + 2);
        for (a, b) in [(io, i1), (io, i2), (i1, i2)] {
            let d = w.pos[a] - w.pos[b];
            let r = d.norm();
            let qq = q(a) * q(b);
            let x = r / s;
            let e = -qq * phys::induction::erf(x) / r;
            let de = -qq * (2.0 / (PI.sqrt() * s) * (-x * x).exp() / r - phys::induction::erf(x) / (r * r));
            out.coulomb += e;
            let fa = d.scale(-de / r);
            out.force[a] += fa;
            out.force[b] -= fa;
            out.virial += -r * de;
        }
        // Springs: E = 1/2 delta^T H delta, delta = (r1 - r0, r2 - r0, theta - theta0).
        let (d1, d2) = (w.pos[i1] - w.pos[io], w.pos[i2] - w.pos[io]);
        let (r1, r2) = (d1.norm(), d2.norm());
        let (u1, u2) = (d1.scale(1.0 / r1), d2.scale(1.0 / r2));
        let cth = u1.dot(u2).clamp(-1.0, 1.0);
        let th = cth.acos();
        let sth = th.sin().max(1e-6);
        let dl = [r1 - p.r0, r2 - p.r0, th - p.th0];
        let mut g = [0.0; 3]; // dE/d(delta)
        let mut e = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                g[a] += p.h[a][b] * dl[b];
                e += 0.5 * dl[a] * p.h[a][b] * dl[b];
            }
        }
        out.spring += e;
        // gradients of the internal coordinates
        let dth_d1 = (u1.scale(cth) - u2).scale(1.0 / (r1 * sth));
        let dth_d2 = (u2.scale(cth) - u1).scale(1.0 / (r2 * sth));
        let f1 = (u1.scale(g[0]) + dth_d1.scale(g[2])).scale(-1.0); // force on H1
        let f2 = (u2.scale(g[1]) + dth_d2.scale(g[2])).scale(-1.0); // force on H2
        out.force[i1] += f1;
        out.force[i2] += f2;
        out.force[io] -= f1 + f2;
        out.virial += d1.dot(f1) + d2.dot(f2);
    }
    // The reciprocal sum over every charge.
    let cpos: Vec<Vec3> = w.pos.clone();
    let cq: Vec<f64> = (0..na).map(q).collect();
    let ph = ew.phases(&cpos);
    let sq = ew.charge_sum(&ph, &cq);
    out.coulomb += ew.energy(&sq) - ew.charge_self(&cq);
    let g = ew.grad_charge_pos(&sq, &ph, &cq);
    for (f, gi) in out.force.iter_mut().zip(g) {
        *f -= gi;
    }
    let zero = vec![Vec3::ZERO; na];
    out.virial += ew.virial(&cpos, &zero, &cq, &[], &[], &[]);
    out
}

fn kinetic(w: &World) -> f64 {
    w.mom.iter().enumerate().map(|(i, p)| p.norm2() / (2.0 * w.mass[i % 3])).sum()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let params = read_params(args.first().expect("a params file"));
    let temperature: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(298.0);
    let n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(216);
    let density: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let total_ps: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(20.0);
    let seed: u64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(1);
    let sp: f64 = std::env::var("PHYS_EWALD_S").ok().and_then(|v| v.parse().ok()).unwrap_or(4.2);
    let acc: f64 = std::env::var("PHYS_EWALD_ACC").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-5);
    let mass = [15.999 * AMU, 1.008 * AMU, 1.008 * AMU];
    let molecule_g = (15.999 + 2.0 * 1.008) * 1.66053906660e-24;
    let per_bohr3 = density / molecule_g * 0.529177210903e-8f64.powi(3);
    let cell = (n as f64 / per_bohr3).cbrt();
    let r_cut = (0.49 * cell).min(16.0);
    let r_on = r_cut - 2.0;
    let ew = Ewald::new(sp, [cell; 3], acc);
    println!("sticks: {n} flexible waters, {density} g/cm3, {temperature} K, box {cell:.2} bohr, cut-off {r_on}-{r_cut}, Ewald s' {sp} ({} k), q_H {:.4}", ew.nk(), params.qh);

    // A lattice of molecules in random orientations, none within 2 A of another.
    let mut rng = Stream::at(0x57_1c_4, seed as u128, 0, Purpose::Positions);
    let side = (n as f64).cbrt().ceil() as usize;
    let spacing = cell / side as f64;
    let mut w = World { n, cell, pos: Vec::new(), mom: vec![Vec3::ZERO; 3 * n], mass };
    let body = [Vec3::ZERO, Vec3 { x: params.r0, y: 0.0, z: 0.0 }, Vec3 { x: params.r0 * params.th0.cos(), y: params.r0 * params.th0.sin(), z: 0.0 }];
    let tau = std::f64::consts::TAU;
    let min_d = 2.0 / 0.529177210903;
    let mut placed: Vec<Vec3> = Vec::new();
    for m in 0..n {
        let centre = Vec3 { x: ((m % side) as f64 + 0.5) * spacing, y: (((m / side) % side) as f64 + 0.5) * spacing, z: ((m / (side * side)) as f64 + 0.5) * spacing };
        let mut best: Option<Vec<Vec3>> = None;
        for _ in 0..2000 {
            let (u1, u2, u3) = (rng.uniform(), rng.uniform(), rng.uniform());
            let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
            let q = phys::math::Quat { w: a * (tau * u2).sin(), v: Vec3 { x: a * (tau * u2).cos(), y: b * (tau * u3).sin(), z: b * (tau * u3).cos() } };
            let atoms: Vec<Vec3> = body.iter().map(|p| centre + q.rotate(*p - body[0].scale(0.0))).collect();
            let ok = placed.iter().all(|pj| {
                let mut d = atoms.iter().map(|pi| {
                    let mut v = *pi - *pj;
                    v.x -= cell * (v.x / cell).round();
                    v.y -= cell * (v.y / cell).round();
                    v.z -= cell * (v.z / cell).round();
                    v.norm()
                });
                d.all(|x| x > min_d)
            });
            if ok {
                best = Some(atoms);
                break;
            }
        }
        let atoms = best.unwrap_or_else(|| body.iter().map(|p| centre + *p).collect());
        placed.extend(atoms.iter().cloned());
        w.pos.extend(atoms);
    }
    // Maxwell-Boltzmann momenta, centre of mass removed.
    let kt = K_B * temperature;
    let mut stream = Stream::at(0xba7_5, seed as u128, 0, Purpose::ThermalNoise);
    for i in 0..3 * n {
        w.mom[i] = stream.normal3().scale((w.mass[i % 3] * kt).sqrt());
    }
    let total_p = w.mom.iter().fold(Vec3::ZERO, |a, p| a + *p);
    let total_m: f64 = (0..3 * n).map(|i| w.mass[i % 3]).sum();
    for i in 0..3 * n {
        let m = w.mass[i % 3];
        w.mom[i] = w.mom[i] - total_p.scale(m / total_m);
    }

    // The tail of the dispersion beyond the cut-off, from the law's own C6 (isotropic, g = 1).
    let rho_o = n as f64 / (cell * cell * cell);
    let rho = [rho_o, 2.0 * rho_o];
    let (mut integral, mut r, h) = (0.0, r_on, 0.001);
    while r < 2000.0 {
        let (sw, _) = switch(r + h / 2.0, r_on, r_cut);
        integral += (1.0 - sw) / (r + h / 2.0).powi(4) * h;
        r += h;
    }
    let mut p_tail = 0.0;
    let mut u_tail_n = 0.0;
    for (s, rs) in rho.iter().enumerate() {
        for (t, rt) in rho.iter().enumerate() {
            p_tail += -2.0 * PI * rs * rt * params.c6[s + t] * integral;
            u_tail_n += -2.0 * PI * rs * rt * params.c6[s + t] * integral;
        }
    }
    let p_tail_bar = p_tail * BAR_PER_HARTREE_PER_BOHR3;
    let u_tail_per_mol = u_tail_n / rho_o * KCAL; // P = rho U/N for the isotropic tail

    let stages: [(&str, f64, usize, f64); 4] = [("easing", 4.0, 3000, 0.02), ("settling", 10.0, 5000, 0.1), ("warming", 20.0, 10000, 0.5), ("production", 20.0, (total_ps / PS / 20.0) as usize, 1.0)];
    let mut out = forces(&w, &params, &ew, r_cut, r_on);
    let t0 = Instant::now();
    let mut samples: Vec<[f64; 7]> = Vec::new();
    let volume = cell * cell * cell;
    let sample = |w: &World, o: &Out| -> [f64; 7] {
        let k = kinetic(w);
        let dof = 3.0 * 3.0 * w.n as f64 - 3.0;
        let tk = 2.0 * k / (dof * K_B);
        let p_bar = (2.0 * k + o.virial) / (3.0 * volume) * BAR_PER_HARTREE_PER_BOHR3;
        let inter = (o.coulomb + o.repulsion - o.dispersion) / w.n as f64 * KCAL;
        [inter, o.spring / w.n as f64 * KCAL, tk, p_bar, o.coulomb / w.n as f64 * KCAL, o.repulsion / w.n as f64 * KCAL, o.dispersion / w.n as f64 * KCAL]
    };
    let mut hist: Vec<f64> = vec![0.0; (0.5 * cell / (0.1 / 0.529177210903)) as usize];
    let mut frames = 0usize;
    let bin = 0.1 / 0.529177210903;
    let file_name = format!("sticks-{temperature}.txt");
    let mut file = std::fs::File::create(&file_name).expect("the output file");
    writeln!(file, "# ps U_inter_per_molecule spring_per_molecule T pressure_bar E_coul E_rep E_disp (kcal/mol a molecule)").ok();
    for (label, dt, steps, tau_ps) in stages {
        let friction = 1.0 / (tau_ps / PS);
        let c = (-friction * dt).exp();
        let every = ((0.1 / PS) / dt) as usize;
        let settle = if label == "production" { steps / 5 } else { usize::MAX };
        for step in 0..steps {
            // BAOAB
            for i in 0..3 * n {
                w.mom[i] += out.force[i].scale(0.5 * dt);
                w.pos[i] += w.mom[i].scale(0.5 * dt / w.mass[i % 3]);
            }
            for i in 0..3 * n {
                let xi = stream.normal3();
                w.mom[i] = w.mom[i].scale(c) + xi.scale(((1.0 - c * c) * w.mass[i % 3] * kt).sqrt());
                w.pos[i] += w.mom[i].scale(0.5 * dt / w.mass[i % 3]);
            }
            out = forces(&w, &params, &ew, r_cut, r_on);
            for i in 0..3 * n {
                w.mom[i] += out.force[i].scale(0.5 * dt);
            }
            if !out.virial.is_finite() {
                println!("  the run blew up in {label}, step {step}");
                return;
            }
            if label == "production" && step >= settle && step % every == 0 {
                let s = sample(&w, &out);
                writeln!(file, "{:.3} {:.4} {:.4} {:.1} {:.1} {:.4} {:.4} {:.4}", step as f64 * dt * PS, s[0], s[1], s[2], s[3], s[4], s[5], s[6]).ok();
                file.flush().ok();
                samples.push(s);
                // g(OO) of the oxygens.
                let os: Vec<Vec3> = (0..n).map(|m| w.pos[3 * m]).collect();
                for a in 0..n {
                    for b in a + 1..n {
                        let mut d = os[b] - os[a];
                        d.x -= cell * (d.x / cell).round();
                        d.y -= cell * (d.y / cell).round();
                        d.z -= cell * (d.z / cell).round();
                        let dist = d.norm();
                        if dist < hist.len() as f64 * bin {
                            hist[(dist / bin) as usize] += 2.0;
                        }
                    }
                }
                frames += 1;
            }
            if step % (steps / 6).max(1) == 0 || step + 1 == steps {
                let s = sample(&w, &out);
                println!("  {label} {:.0}% after {:.0} s: inter U/N {:.2}, spring {:.2} kcal/mol, T {:.0} K, P {:.0} bar", 100.0 * step as f64 / steps as f64, t0.elapsed().as_secs_f64(), s[0], s[1], s[2], s[3]);
                std::io::stdout().flush().ok();
            }
        }
    }
    // Averages with block errors.
    let nb = 5usize;
    let per = samples.len() / nb;
    let mean = |k: usize| samples.iter().map(|s| s[k]).sum::<f64>() / samples.len() as f64;
    let spread = |k: usize| {
        let b: Vec<f64> = (0..nb).map(|c| samples[c * per..(c + 1) * per].iter().map(|s| s[k]).sum::<f64>() / per as f64).collect();
        let m = b.iter().sum::<f64>() / nb as f64;
        (b.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (nb - 1) as f64 / nb as f64).sqrt()
    };
    let rt = K_B * temperature * KCAL;
    println!("sticks at {temperature} K, {n} molecules, {density} g/cm3 ({} samples over {:.0} ps):", samples.len(), total_ps * 0.8);
    println!("  intermolecular U/N {:.3} +- {:.3} kcal/mol (with the dispersion tail {:.3}); heat of vaporisation -U/N + RT = {:.3} ({:.3} with the tail)", mean(0), spread(0), mean(0) + u_tail_per_mol, -mean(0) + rt, -(mean(0) + u_tail_per_mol) + rt);
    println!("  spring energy {:.3} kcal/mol a molecule (equipartition: 3 x 1/2 RT = {:.3})", mean(1), 1.5 * rt);
    println!("  pressure {:.0} +- {:.0} bar (block); with the dispersion tail {:.0} bar [tail {:.0}]", mean(3), spread(3), mean(3) + p_tail_bar, p_tail_bar);
    println!("  temperature {:.1} K; energy parts a molecule: Coulomb {:.2}, repulsion {:.2}, dispersion {:.2}", mean(2), mean(4), mean(5), mean(6));
    // g(OO): first peak, first minimum, coordination.
    let norm = |k: usize| {
        let r_in = k as f64 * bin;
        let r_out = (k + 1) as f64 * bin;
        let shell = 4.0 / 3.0 * PI * (r_out.powi(3) - r_in.powi(3));
        hist[k] / (frames as f64 * n as f64 * shell * rho_o)
    };
    let g: Vec<f64> = (0..hist.len()).map(norm).collect();
    let (mut peak, mut peak_r) = (0.0, 0.0);
    for (k, v) in g.iter().enumerate() {
        let r = (k as f64 + 0.5) * bin * 0.529177210903;
        if r > 2.3 && r < 3.5 && *v > peak {
            peak = *v;
            peak_r = r;
        }
    }
    let mut kmin = 0;
    for k in 0..g.len() {
        let r = (k as f64 + 0.5) * bin * 0.529177210903;
        if r > peak_r + 0.4 && r < 4.6 && (kmin == 0 || g[k] < g[kmin]) {
            kmin = k;
        }
    }
    let mut cn = 0.0;
    for k in 0..=kmin {
        let r_in = k as f64 * bin;
        let r_out = (k + 1) as f64 * bin;
        cn += g[k] * rho_o * 4.0 / 3.0 * PI * (r_out.powi(3) - r_in.powi(3));
    }
    println!("  g(OO): first peak {:.2} at {:.2} A, first minimum {:.2} at {:.2} A, {:.2} neighbours inside it", peak, peak_r, g[kmin], (kmin as f64 + 0.5) * bin * 0.529177210903, cn);
    // The last frame, for scripts/tetra.py.
    let mut snap = std::fs::File::create(format!("sticks-{temperature}.snap")).expect("snapshot");
    writeln!(snap, "box {} {} {}", cell, cell, cell).ok();
    for m in 0..n {
        let (a, b, c) = (w.pos[3 * m], w.pos[3 * m + 1], w.pos[3 * m + 2]);
        writeln!(snap, "mol {m} {} {} {} {} {} {} {} {} {}", a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z).ok();
    }
    println!("  wall time {:.0} s", t0.elapsed().as_secs_f64());
}
