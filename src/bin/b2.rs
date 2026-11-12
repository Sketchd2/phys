//! `phys-b2`: a law's second virial coefficient, by direct integration over
//! the pair — the vapour side of a boiling point (`docs/SOLUTIONS.md` 5.1's
//! pressure rule: a vapour is treated by its second virial coefficient from
//! the same law, not as an ideal gas).
//!
//! ```sh
//! phys-b2 methane law.txt 111 150 200 [--orientations 4000] [--seed 1]
//! ```
//!
//! `B2(T) = -1/2 integral <exp(-u/kT) - 1> 4 pi R^2 dR`, `R` between the two
//! molecules' centres of mass and the average over both orientations. The
//! orientations are drawn once from a seeded stream and used at every radius
//! and temperature, so the curve is smooth in both and a run replays exactly;
//! the same draw split in two gives the error printed. The pair's whole law
//! is used, induction included and nothing switched off at a cut-off, out to
//! 40 bohr; beyond that the dispersion's `-C6/R^6` is integrated in closed
//! form. `--cut ON CUT` (bohr) instead multiplies the pair's energy by the
//! liquid's own switch on the centres' separation and adds no tail: the second
//! virial coefficient of the truncated law a `phys-coexist` slab runs, for
//! comparing the two like for like. Prints `T  B2 (bohr^3)  B2 (cm^3/mol)  +- (cm^3/mol)`.

use phys::electrons::grow::{equivalent_atoms, Resume};
use phys::electrons::molecule::Molecule;
use phys::liquid::{extra_sites_from_text, pair_energy, Kind, Liquid, PairEnergy, Polarisable, SiteLaw, SiteSite, K_B};
use phys::math::{Quat, Vec3};

const CM3_PER_MOL: f64 = 0.529177210903e-8 * 0.529177210903e-8 * 0.529177210903e-8 * 6.02214076e23;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args[0].clone();
    let law_text = std::fs::read_to_string(&args[1]).unwrap_or_else(|_| panic!("no {}", args[1]));
    let flag = |f: &str| args.iter().position(|a| a == f).and_then(|i| args.get(i + 1)).cloned();
    let orientations: usize = flag("--orientations").and_then(|v| v.parse().ok()).unwrap_or(4000);
    let seed: u64 = flag("--seed").and_then(|v| v.parse().ok()).unwrap_or(1);
    let cut: Option<(f64, f64)> = args.iter().position(|a| a == "--cut").map(|i| (args[i + 1].parse().expect("r_on"), args[i + 2].parse().expect("r_cut")));
    let temperatures: Vec<f64> = args[2..].iter().take_while(|a| !a.starts_with("--")).map(|a| a.parse().expect("a temperature")).collect();
    let law = SiteSite::from_text(&law_text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&law_text, law.charge.len());
    let polarised = Polarisable { law: &law, alpha: alpha.clone() };
    let law_ref: &dyn SiteLaw = if alpha.iter().any(|a| *a > 0.0) { &polarised } else { &law };
    let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).expect("a growth state")).expect("readable");
    let z = phys::pairs::atoms_of(&name);
    let mono = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
    let kind = Kind::of_molecule_sites(&z, &state.positions, &equivalent_atoms(&mono), SiteSite::bisector_from_text(&law_text), &extra_sites_from_text(&law_text));
    // The molecule's dispersion coefficient: every site pair's.
    let types = law.charge.len();
    let c6: f64 = kind.types.iter().flat_map(|a| kind.types.iter().map(move |b| (a, b))).map(|(a, b)| law.pair[a * types + b][2]).sum();
    let mut s = Liquid::noise(seed);
    let mut turn = || {
        let (u1, u2, u3) = (s.uniform(), s.uniform(), s.uniform());
        let tau = std::f64::consts::TAU;
        let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
        Quat { w: a * (tau * u2).sin(), v: Vec3 { x: a * (tau * u2).cos(), y: b * (tau * u3).sin(), z: b * (tau * u3).cos() } }
    };
    let draws: Vec<(Quat, Quat)> = (0..orientations).map(|_| (turn(), turn())).collect();
    let (r_max, h) = (40.0f64, 0.1f64);
    let steps = (r_max / h) as usize;
    // sums[half][temperature]
    let mut sums = vec![vec![0.0f64; temperatures.len()]; 2];
    for i in 0..steps {
        let r = (i as f64 + 0.5) * h;
        let shift = Vec3 { x: r, y: 0.0, z: 0.0 };
        for (d, (qa, qb)) in draws.iter().enumerate() {
            let a: Vec<(Vec3, usize)> = kind.sites.iter().zip(&kind.types).map(|(p, t)| (qa.rotate(*p), *t)).collect();
            let b: Vec<(Vec3, usize)> = kind.sites.iter().zip(&kind.types).map(|(p, t)| (qb.rotate(*p) + shift, *t)).collect();
            let u = match cut {
                Some((on, off)) if r >= off => continue,
                Some((on, off)) => phys::liquid::switch(r, on, off).0 * pair_energy(law_ref, &PairEnergy { a, b, energy: 0.0 }),
                None => pair_energy(law_ref, &PairEnergy { a, b, energy: 0.0 }),
            };
            for (k, t) in temperatures.iter().enumerate() {
                // A law with a hole would send this to infinity; it is capped
                // only so that the hole shows as an absurd number, not a NaN.
                let f = (-u / (K_B * t)).min(700.0).exp() - 1.0;
                sums[d % 2][k] += f * 4.0 * std::f64::consts::PI * r * r * h;
            }
        }
    }
    for (k, t) in temperatures.iter().enumerate() {
        let per = orientations as f64 / 2.0;
        let tail = if cut.is_some() { 0.0 } else { -2.0 * std::f64::consts::PI * c6 / (K_B * t) / (3.0 * r_max.powi(3)) };
        let halves = [-0.5 * sums[0][k] / per + tail, -0.5 * sums[1][k] / per + tail];
        let b2 = 0.5 * (halves[0] + halves[1]);
        println!("{t} {b2:.6e} {:.3} {:.3}", b2 * CM3_PER_MOL, 0.5 * (halves[0] - halves[1]).abs() * CM3_PER_MOL);
    }
}
