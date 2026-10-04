//! Interaction energies of sampled pairs of one molecule — `PLAY.md` E8,
//! route A's data.
//!
//! ```sh
//! cargo run --release --bin phys-pairs -- water 300 [element|grown]
//! ```
//!
//! The molecule is taken at the shape its basis growth relaxed it to
//! (`grow-<name>.state`), its atoms typed by symmetry class. Each pair is a
//! seeded draw — its index decides it, so any one can be redone — of a
//! relative orientation and a separation of centres, refused if two atoms
//! sit closer than 0.75 of their van der Waals contact. Its interaction is
//! counterpoise corrected: the pair, and each partner in the pair's basis with
//! the other's nuclei and electrons removed. One field is solved, with PBE
//! exchange and the non-local correlation in it on the coarse grid, and two
//! energies are taken from it on the fine grid: the unfitted form's, and the
//! published vdW-DF1's (revPBE exchange) by swapping the semilocal energy on
//! the same density, as Klimes et al. did. Runs on the CPU: the GPU's
//! precision is tested after E8, so E8's data does not rest on it.
//!
//! Each pair appends one line to `pairs-<name>.txt` and a run carries on from
//! the last line written.

use crate::electrons::functional::Functional;
use crate::electrons::grow::{equivalent_atoms, extras, Resume};
use crate::electrons::molecule::Molecule;
use crate::electrons::scf::{exchange_correlation, solve, Batches};
use crate::electrons::vdw::{energy_on_finer_grid, NonlocalSpec, Z_AB_DF1};
use crate::math::{Quat, Vec3};
use crate::rng::{Purpose, Stream};
use std::io::Write;
use std::time::Instant;

/// Separations of centres sampled, bohr: from inside the repulsive wall
/// (refused where atoms collide) to where only the long-range tail is left.
const SEPARATION: (f64, f64) = (4.0, 15.0);

/// Run the driver for `args` (`name count [element|grown]`), writing
/// `pairs-<name><tag>.txt`: the binary passes no tag; a binary that computes the
/// final non-local energy another way passes one, so its pairs sit beside the
/// reference ones to be compared.
pub fn run(args: &[String], tag: &str) {
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let count: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(300);
    let grown = args.get(2).map(|s| s == "grown").unwrap_or(false);
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).unwrap_or_else(|_| panic!("no grow-{name}.state"))).expect("a readable state");
    let z: Vec<u32> = match name.as_str() {
        "water" => vec![8, 1, 1],
        "methane" => vec![6, 1, 1, 1, 1],
        "ammonia" => vec![7, 1, 1, 1],
        "methanol" => vec![6, 8, 1, 1, 1, 1],
        other => panic!("no atoms known for {other}"),
    };
    let mono = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let types = equivalent_atoms(&mono);
    let extra = extras(&mono, &state.ladders, &state.chosen);
    // Centred on the mass, so a separation is between centres.
    let masses: Vec<f64> = z.iter().map(|&zz| crate::chem::elements::Element(zz as u8).mass_kg().expect("a mass")).collect();
    let total: f64 = masses.iter().sum();
    let mut com = [0.0; 3];
    for (p, m) in mono.positions.iter().zip(&masses) {
        for k in 0..3 {
            com[k] += p[k] * m / total;
        }
    }
    let body: Vec<Vec3> = mono.positions.iter().map(|p| Vec3 { x: p[0] - com[0], y: p[1] - com[1], z: p[2] - com[2] }).collect();
    let contact: Vec<f64> = z.iter().map(|&zz| crate::chem::elements::Element(zz as u8).vdw_radius().unwrap_or(1.5e-10) / 0.529177210903e-10).collect();
    let out = format!("pairs-{name}{tag}.txt");
    let done = std::fs::read_to_string(&out).map(|t| t.lines().filter(|l| !l.starts_with('#')).count()).unwrap_or(0);
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    if done == 0 {
        writeln!(file, "# index separation_bohr E_int_pbe_x E_int_revpbe_x (hartree) | then x y z type for each atom of the first molecule and of the second (bohr)").ok();
    }
    println!("{name}: {} atoms, types {types:?}, {} basis; {done} pairs done, {count} wanted", z.len(), if grown { "grown" } else { "per-element" });
    let n = z.len();
    for k in done..count {
        let t = Instant::now();
        // The draw for pair k: its own stream, so it does not depend on what
        // was refused before it.
        let mut s = Stream::at(0x7061_6972_7300 ^ name.len() as u64, k as u128, 0, Purpose::Positions);
        let (place_b, rot_a, rot_b, sep) = loop {
            let qa = random_rotation(&mut s);
            let qb = random_rotation(&mut s);
            let dir = s.direction();
            let sep = s.range(SEPARATION.0, SEPARATION.1);
            let a: Vec<Vec3> = body.iter().map(|p| qa.rotate(*p)).collect();
            let b: Vec<Vec3> = body.iter().map(|p| qb.rotate(*p) + dir.scale(sep)).collect();
            let clash = (0..n).any(|i| (0..n).any(|j| (a[i] - b[j]).norm() < 0.75 * (contact[i] + contact[j])));
            if !clash {
                break (dir.scale(sep), qa, qb, sep);
            }
        };
        let a: Vec<Vec3> = body.iter().map(|p| rot_a.rotate(*p)).collect();
        let b: Vec<Vec3> = body.iter().map(|p| rot_b.rotate(*p) + place_b).collect();
        let pair = Molecule { z: z.iter().chain(&z).cloned().collect(), positions: a.iter().chain(&b).map(|v| [v.x, v.y, v.z]).collect(), charge: 0, unpaired: 0 };
        let base = if grown {
            let ex: Vec<Vec<(usize, f64)>> = extra.iter().chain(&extra).cloned().collect();
            pair.problem_with(Functional::Pbe, Some(&ex))
        } else {
            pair.problem(Functional::Pbe)
        };
        let electrons: f64 = z.iter().map(|&zz| zz as f64).sum();
        let first: Vec<usize> = (0..n).collect();
        let second: Vec<usize> = (n..2 * n).collect();
        let all: Vec<usize> = (0..2 * n).collect();
        let mut e_pbe = [0.0; 3];
        let mut e_rev = [0.0; 3];
        // Where the pair's time goes: the three fields, the fine-grid
        // non-local energies, and the exchange swap for the second partner.
        let (mut t_solve, mut t_fine, mut t_swap) = (0.0, 0.0, 0.0);
        for (slot, (keep, ne)) in [(all, 2.0 * electrons), (first, electrons), (second, electrons)].into_iter().enumerate() {
            let mut p = base.with_ghosts(&keep, ne);
            p.functional = Functional::PbeXLdaC;
            p.nonlocal = Some(NonlocalSpec::in_the_field(Z_AB_DF1));
            let ts = Instant::now();
            let sol = solve(&p, 200, 1e-10);
            assert!(sol.converged, "pair {k}: a field did not converge");
            t_solve += ts.elapsed().as_secs_f64();
            let ts = Instant::now();
            let fine = energy_on_finer_grid(&p, &sol, 50, 12);
            t_fine += ts.elapsed().as_secs_f64();
            let ts = Instant::now();
            let atoms: Vec<([f64; 3], f64)> = p.nuclei.iter().zip(&p.sizes).map(|((_, q), r)| (*q, *r)).collect();
            let grid = crate::electrons::grid::molecular_pruned(&atoms, p.radial, p.theta, p.prune);
            let batches = Batches::new(&p.basis, &grid);
            let (x_pbe, _, _) = exchange_correlation(&p.basis, &batches, Functional::PbeXLdaC, &sol.density_alpha, &sol.density_beta);
            let (x_rev, _, _) = exchange_correlation(&p.basis, &batches, Functional::RevPbeXLdaC, &sol.density_alpha, &sol.density_beta);
            e_pbe[slot] = fine;
            e_rev[slot] = fine - x_pbe + x_rev;
            t_swap += ts.elapsed().as_secs_f64();
        }
        let int_pbe = e_pbe[0] - e_pbe[1] - e_pbe[2];
        let int_rev = e_rev[0] - e_rev[1] - e_rev[2];
        let mut line = format!("{k} {sep:.6} {int_pbe:.10e} {int_rev:.10e} |");
        for (v, ty) in a.iter().zip(&types).chain(b.iter().zip(&types)) {
            line += &format!(" {:.8} {:.8} {:.8} {ty}", v.x, v.y, v.z);
        }
        writeln!(file, "{line}").expect("the line written");
        file.flush().ok();
        println!("pair {k}: R {sep:.2} bohr, E_int {:.4} / {:.4} kcal/mol (PBE x / revPBE x), {:.0} s ({t_solve:.0} s fields, {t_fine:.0} s fine non-local, {t_swap:.0} s exchange swap; {} functions)", int_pbe * 627.509474, int_rev * 627.509474, t.elapsed().as_secs_f64(), base.basis.size);
        std::io::stdout().flush().ok();
    }
}

/// A rotation drawn uniformly (Shoemake's method).
fn random_rotation(s: &mut Stream) -> Quat {
    let (u1, u2, u3) = (s.uniform(), s.uniform(), s.uniform());
    let tau = std::f64::consts::TAU;
    let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
    Quat { w: a * (tau * u2).sin(), v: Vec3 { x: a * (tau * u2).cos(), y: b * (tau * u3).sin(), z: b * (tau * u3).cos() } }
}
