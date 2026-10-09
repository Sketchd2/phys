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
//! sit closer than 0.75 (from pair 250, 0.64) of their van der Waals contact. Its interaction is
//! counterpoise corrected: the pair, and each partner in the pair's basis with
//! the other's nuclei and electrons removed. One field is solved, with PBE
//! exchange and the non-local correlation in it on the coarse grid, and two
//! energies are taken from it on the fine grid: the unfitted form's, and the
//! published vdW-DF1's (revPBE exchange) by swapping the semilocal energy on
//! the same density, as Klimes et al. did. Runs on the CPU: the GPU's
//! precision is tested after E8, so E8's data does not rest on it.
//!
//! Each pair appends one line to `pairs-<name>.txt` and a run skips the indices
//! the file already holds, whatever order they were written in (`queue` hands
//! them to many machines).

use crate::electrons::functional::Functional;
use crate::electrons::grow::{equivalent_atoms, extras, Resume};
use crate::electrons::molecule::Molecule;
use crate::electrons::scf::{exchange_correlation, solve, Batches};
use crate::electrons::vdw::{energy_on_finer_grid, NonlocalSpec, Z_AB_DF1};
use crate::math::{Quat, Vec3};
use crate::rng::{Purpose, Stream};
use std::io::Write;
use std::time::{Duration, Instant};

/// Separations of centres sampled, bohr: from inside the repulsive wall
/// (refused where atoms collide) to where only the long-range tail is left.
const SEPARATION: (f64, f64) = (4.0, 15.0);

/// The first pairs, drawn over [`SEPARATION`]: of 124, only 27 fell inside 7
/// bohr (where a liquid's neighbours are and the law's error is, 0.5-1.5
/// kcal/mol there against 0.04-0.08 beyond 10) and 54 beyond 10, where the
/// energies are nearly nothing. Pairs from this index on are drawn over the
/// shorter range, and the earlier ones stay exactly as they were drawn; every
/// pair's geometry is on its line either way. The boundary is 126 and not the
/// 124 the analysis was made on: the run in progress had two pairs in flight
/// when the binary changed, and the data is the authority, so the constant
/// is the one that reproduces every line of `pairs-water-gpu.txt` exactly.
const UNIFORM_PAIRS: usize = 126;

/// How close two atoms may come, as a fraction of the sum of their van der
/// Waals radii. 0.75 for the first pairs, which for O...H is 2.04 A and so
/// refused the hydrogen-bonded geometries themselves (a real bond is 1.8-2.0
/// A, the dimer's minimum 1.95); from `CLOSER_FROM` on it is 0.64 (1.74 A for
/// O...H), so the region a liquid's strongest bonds sit in is sampled, with the
/// repulsive wall just inside it.
const CLOSER_FROM: usize = 250;

fn closest_contact(k: usize) -> f64 {
    if k < CLOSER_FROM { 0.75 } else { 0.64 }
}

/// The separations pair `k` is drawn over.
fn separation_for(k: usize) -> (f64, f64) {
    if k < UNIFORM_PAIRS { SEPARATION } else { (4.2, 11.0) }
}

/// What the electronic structure needs of one pair beyond its geometry.
struct Setup<'a> {
    z: &'a [u32],
    extra: &'a [Vec<(usize, f64)>],
    grown: bool,
}

/// One pair's counterpoise interaction energies and where their time went.
struct Interaction {
    pbe: f64,
    rev: f64,
    t_solve: f64,
    t_fine: f64,
    t_swap: f64,
    iterations: Vec<usize>,
    functions: usize,
}

/// The counterpoise interaction energy of molecules `a` and `b` (atom
/// positions, bohr): the pair, and each partner in the pair's basis with the
/// other's nuclei and electrons removed; one field each with PBE exchange and
/// the non-local correlation in it, and two energies from it on the fine grid,
/// the unfitted form's and published vdW-DF1's (revPBE exchange swapped in on
/// the same density).
fn interact(setup: &Setup, k: usize, a: &[Vec3], b: &[Vec3]) -> Interaction {
    let (z, extra, grown) = (setup.z, setup.extra, setup.grown);
    let n = z.len();
    let pair = Molecule { z: z.iter().chain(z.iter()).cloned().collect(), positions: a.iter().chain(b.iter()).map(|v| [v.x, v.y, v.z]).collect(), charge: 0, unpaired: 0 };
    let base = if grown {
        let ex: Vec<Vec<(usize, f64)>> = extra.iter().chain(extra.iter()).cloned().collect();
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
    let mut iterations = Vec::new();
    for (slot, (keep, ne)) in [(all, 2.0 * electrons), (first, electrons), (second, electrons)].into_iter().enumerate() {
        let mut p = base.with_ghosts(&keep, ne);
        p.functional = Functional::PbeXLdaC;
        p.nonlocal = Some(NonlocalSpec::in_the_field(Z_AB_DF1));
        let ts = Instant::now();
        let sol = solve(&p, 200, 1e-10);
        assert!(sol.converged, "pair {k}: a field did not converge");
        t_solve += ts.elapsed().as_secs_f64();
        iterations.push(sol.iterations);
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
    Interaction { pbe: int_pbe, rev: int_rev, t_solve, t_fine, t_swap, iterations, functions: base.basis.size }
}

/// One molecule as the pair work needs it: its atoms, the symmetry class of
/// each, the basis extras its growth chose, and its shape centred on the mass
/// (so a separation is between centres) with the van der Waals contact of each
/// atom. Loaded from `grow-<name>.state`, which every machine doing this work
/// must hold the same copy of.
pub struct Monomer {
    pub name: String,
    z: Vec<u32>,
    types: Vec<usize>,
    extra: Vec<Vec<(usize, f64)>>,
    body: Vec<Vec3>,
    contact: Vec<f64>,
    masses: Vec<f64>,
    grown: bool,
    /// Compute pairs at Hartree-Fock + RI-MP2 instead of with the density
    /// functional (`--mp2`): a line then carries the total interaction (HF
    /// plus correlation) where the PBE-exchange energy goes and the
    /// Hartree-Fock part where the revPBE one does.
    pub mp2: bool,
    /// Draw pairs where a law says the liquid goes (`--bias law.txt [kT]`),
    /// not uniformly: see [`Monomer::biased_pair`].
    pub bias: Option<Bias>,
}

/// A law and a temperature to draw pairs by.
pub struct Bias {
    law: crate::liquid::SiteSite,
    alpha: Vec<f64>,
    bisector: Option<f64>,
    /// kcal/mol.
    kt: f64,
}

impl Bias {
    pub fn from_file(path: &str, kt: f64) -> Bias {
        let text = std::fs::read_to_string(path).unwrap_or_else(|_| panic!("no {path}"));
        let law = crate::liquid::SiteSite::from_text(&text).expect("a readable law");
        let alpha = crate::liquid::SiteSite::alpha_from_text(&text, law.charge.len());
        let bisector = crate::liquid::SiteSite::bisector_from_text(&text).map(|b| b.1);
        Bias { law, alpha, bisector, kt }
    }
}

/// The atoms of a molecule `phys-grow` knows by name.
fn atoms_of(name: &str) -> Vec<u32> {
    match name {
        "water" => vec![8, 1, 1],
        "methane" => vec![6, 1, 1, 1, 1],
        "ammonia" => vec![7, 1, 1, 1],
        "methanol" => vec![6, 8, 1, 1, 1, 1],
        other => panic!("no atoms known for {other}"),
    }
}

impl Monomer {
    pub fn load(name: &str, grown: bool) -> Monomer {
        let state = Resume::from_text(&std::fs::read_to_string(format!("grow-{name}.state")).unwrap_or_else(|_| panic!("no grow-{name}.state"))).expect("a readable state");
        let z = atoms_of(name);
        let mono = Molecule { z: z.clone(), positions: state.positions.clone(), charge: 0, unpaired: 0 };
        let types = equivalent_atoms(&mono);
        let extra = extras(&mono, &state.ladders, &state.chosen);
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
        Monomer { name: name.to_string(), z, types, extra, body, contact, masses, grown, mp2: false, bias: None }
    }

    /// Pair `k` of the random draw: two molecules' atom positions and the
    /// separation of their centres. The draw is its own stream, so it does not
    /// depend on what was refused before it nor on which machine makes it.
    pub fn random_pair(&self, k: usize) -> (Vec<Vec3>, Vec<Vec3>, f64) {
        let n = self.z.len();
        let mut s = Stream::at(0x7061_6972_7300 ^ self.name.len() as u64, k as u128, 0, Purpose::Positions);
        loop {
            let qa = random_rotation(&mut s);
            let qb = random_rotation(&mut s);
            let dir = s.direction();
            let (lo, hi) = separation_for(k);
            let sep = s.range(lo, hi);
            let a: Vec<Vec3> = self.body.iter().map(|p| qa.rotate(*p)).collect();
            let b: Vec<Vec3> = self.body.iter().map(|p| qb.rotate(*p) + dir.scale(sep)).collect();
            let near = closest_contact(k);
            let clash = (0..n).any(|i| (0..n).any(|j| (a[i] - b[j]).norm() < near * (self.contact[i] + self.contact[j])));
            if !clash {
                return (a, b, sep);
            }
        }
    }

    /// Pair `k` drawn where a law puts a liquid's neighbours: positions and
    /// orientations as in [`Monomer::random_pair`], over the separations a
    /// first shell and a second occupy (2.2 to 4.2 A), accepted with
    /// probability `exp(-(E - E0) / kT)` of the law's energy `E` (induction
    /// included), `E0` being 6 kcal/mol below zero so that every pair at least
    /// that bound is kept. The uniform draw put 10 of 219 pairs inside 2.8 A
    /// of oxygen separation, where the hydrogen bonds are and the fitted law
    /// was 0.7-0.8 kcal/mol too attractive; this puts the pairs where the
    /// law is asked to be right. The pair is a function of `k` alone.
    pub fn biased_pair(&self, k: usize, bias: &Bias) -> (Vec<Vec3>, Vec<Vec3>, f64) {
        let n = self.z.len();
        let mut s = Stream::at(0x6269_6173_5f70 ^ self.name.len() as u64, k as u128, 0, Purpose::Positions);
        let atoms = |m: &[Vec3]| -> Vec<(Vec3, usize)> {
            let mut v: Vec<(Vec3, usize)> = m.iter().zip(&self.types).map(|(p, t)| (*p, *t)).collect();
            if let Some(d) = bias.bisector {
                let site = bias.law.charge.len() - 1;
                v = crate::liquid::with_bisector_site(&v, site, d);
            }
            v
        };
        loop {
            let qa = random_rotation(&mut s);
            let qb = random_rotation(&mut s);
            let dir = s.direction();
            let sep = s.range(4.2, 8.0);
            let a: Vec<Vec3> = self.body.iter().map(|p| qa.rotate(*p)).collect();
            let b: Vec<Vec3> = self.body.iter().map(|p| qb.rotate(*p) + dir.scale(sep)).collect();
            let near = 0.64;
            if (0..n).any(|i| (0..n).any(|j| (a[i] - b[j]).norm() < near * (self.contact[i] + self.contact[j]))) {
                continue;
            }
            let pe = crate::liquid::PairEnergy { a: atoms(&a), b: atoms(&b), energy: 0.0 };
            let e = crate::liquid::pair_energy(&crate::liquid::Polarisable { law: &bias.law, alpha: bias.alpha.clone() }, &pe) * 627.509474;
            let p = (-(e + 6.0) / bias.kt).exp().min(1.0);
            if s.uniform() < p {
                return (a, b, sep);
            }
        }
    }

    /// The interaction of the pair `a`, `b`, and the line that records it.
    fn compute(&self, k: usize, a: &[Vec3], b: &[Vec3], sep: f64) -> (Interaction, String) {
        if self.mp2 {
            let r = self.cluster_interaction(&[a.to_vec(), b.to_vec()], false, None, None, true);
            let total = r.pbe + r.rev;
            let mut line = format!("{k} {sep:.6} {total:.10e} {:.10e} |", r.pbe);
            for (v, ty) in a.iter().zip(&self.types).chain(b.iter().zip(&self.types)) {
                line += &format!(" {:.8} {:.8} {:.8} {ty}", v.x, v.y, v.z);
            }
            return (Interaction { pbe: total, rev: r.pbe, t_solve: r.t_solve, t_fine: r.t_fine, t_swap: r.t_swap, iterations: r.iterations, functions: r.functions }, line);
        }
        let setup = Setup { z: &self.z, extra: &self.extra, grown: self.grown };
        let r = interact(&setup, k, a, b);
        let mut line = format!("{k} {sep:.6} {:.10e} {:.10e} |", r.pbe, r.rev);
        for (v, ty) in a.iter().zip(&self.types).chain(b.iter().zip(&self.types)) {
            line += &format!(" {:.8} {:.8} {:.8} {ty}", v.x, v.y, v.z);
        }
        (r, line)
    }
}

/// The header of a pairs file of the random draw.
pub const RANDOM_HEADER: &str = "# index separation_bohr E_int_pbe_x E_int_revpbe_x (hartree) | then x y z type for each atom of the first molecule and of the second (bohr)";

/// The header of a pairs file taken from a snapshot.
pub fn snapshot_header(snapshot: &str) -> String {
    format!("# pairs of molecules taken from {snapshot}: index centroid_separation_bohr E_int_pbe_x E_int_revpbe_x (hartree) | then x y z type for each atom of the first molecule and of the second (bohr)")
}

/// The indices a pairs file already holds. Lines may have been written out of
/// order (a queue hands them to many machines), so a run skips by index and
/// does not count.
pub fn done_indices(path: &str) -> std::collections::HashSet<usize> {
    std::fs::read_to_string(path).map(|t| t.lines().filter(|l| !l.starts_with('#')).filter_map(|l| l.split_whitespace().next()?.parse().ok()).collect()).unwrap_or_default()
}

/// Run the driver for `args` (`name count [element|grown]`), writing
/// `pairs-<name><tag>.txt`: the binary passes no tag; a binary that computes the
/// final non-local energy another way passes one, so its pairs sit beside the
/// reference ones to be compared.
pub fn run(args: &[String], tag: &str) {
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let count: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(300);
    let grown = args.get(2).map(|s| s == "grown").unwrap_or(false);
    let mut mono = Monomer::load(&name, grown);
    mono.mp2 = args.iter().any(|a| a == "--mp2");
    if let Some(i) = args.iter().position(|a| a == "--bias") {
        let kt: f64 = args.get(i + 2).and_then(|v| v.parse().ok()).unwrap_or(2.0);
        mono.bias = Some(Bias::from_file(args.get(i + 1).expect("--bias needs a law file"), kt));
    }
    let out = format!("pairs-{name}{tag}{}{}.txt", if mono.bias.is_some() { "-bias" } else { "" }, if mono.mp2 { "-mp2" } else { "" });
    let done = done_indices(&out);
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    if done.is_empty() {
        if mono.mp2 {
            writeln!(file, "# MP2 (Hartree-Fock + RI-MP2, frozen cores, counterpoise): index separation_bohr E_int_total E_int_hartree_fock (hartree) | then x y z type for each atom of the first molecule and of the second (bohr)").ok();
        } else {
            writeln!(file, "{RANDOM_HEADER}").ok();
        }
    }
    println!("{name}: {} atoms, types {:?}, {} basis; {} pairs done, {count} wanted", mono.z.len(), mono.types, if grown { "grown" } else { "per-element" }, done.len());
    for k in 0..count {
        if done.contains(&k) {
            continue;
        }
        let t = Instant::now();
        let (a, b, sep) = match &mono.bias {
            Some(b) => mono.biased_pair(k, b),
            None => mono.random_pair(k),
        };
        let (r, line) = mono.compute(k, &a, &b, sep);
        writeln!(file, "{line}").expect("the line written");
        file.flush().ok();
        println!("pair {k}: R {sep:.2} bohr, E_int {:.4} / {:.4} kcal/mol (PBE x / revPBE x), {:.0} s ({:.0} s fields, {:.0} s fine non-local, {:.0} s exchange swap; {} functions; field iterations {:?})", r.pbe * 627.509474, r.rev * 627.509474, t.elapsed().as_secs_f64(), r.t_solve, r.t_fine, r.t_swap, r.functions, r.iterations);
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

/// The pairs of molecules a simulated liquid's snapshot holds, in the seeded
/// order they are taken in.
pub struct Snapshot {
    pub path: String,
    /// A hash of the file's text, so a machine holding a different snapshot
    /// under the same name is refused rather than computing other pairs.
    pub fingerprint: u64,
    mols: Vec<Vec<Vec3>>,
    cell: [f64; 3],
    candidates: Vec<(usize, usize, Vec3, f64)>,
}

/// FNV-1a over bytes: a fingerprint, not a security measure.
pub fn fnv(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in bytes {
        h = (h ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

impl Snapshot {
    /// Candidates are every pair of molecules whose centroids lie within 10.5
    /// bohr by the nearest image, put in a seeded order.
    pub fn load(path: &str, atoms: usize) -> Snapshot {
        let text = std::fs::read_to_string(path).unwrap_or_else(|_| panic!("no {path}"));
        let mut cell = [0.0; 3];
        let mut mols: Vec<Vec<Vec3>> = Vec::new();
        for line in text.lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            match w.first().copied() {
                Some("box") => {
                    for k in 0..3 {
                        cell[k] = w[1 + k].parse().expect("a box edge");
                    }
                }
                Some("mol") => {
                    let v: Vec<f64> = w[2..].iter().map(|x| x.parse().expect("a coordinate")).collect();
                    assert_eq!(v.len(), 3 * atoms, "a molecule line with {} numbers", v.len());
                    mols.push(v.chunks(3).map(|c| Vec3 { x: c[0], y: c[1], z: c[2] }).collect());
                }
                _ => {}
            }
        }
        let centroid = |m: &Vec<Vec3>| m.iter().fold(Vec3::ZERO, |acc, p| acc + *p).scale(1.0 / atoms as f64);
        let cents: Vec<Vec3> = mols.iter().map(centroid).collect();
        let image = |mut d: Vec3| {
            d.x -= cell[0] * (d.x / cell[0]).round();
            d.y -= cell[1] * (d.y / cell[1]).round();
            d.z -= cell[2] * (d.z / cell[2]).round();
            d
        };
        let mut candidates: Vec<(usize, usize, Vec3, f64)> = Vec::new();
        for i in 0..mols.len() {
            for j in i + 1..mols.len() {
                let d = image(cents[j] - cents[i]);
                if d.norm() < 10.5 {
                    // The shift that carries j's centroid to its nearest image of i's.
                    candidates.push((i, j, d - (cents[j] - cents[i]), d.norm()));
                }
            }
        }
        let mut order = Stream::at(0x6c69_7175_6964, 0, 0, Purpose::Positions);
        for i in (1..candidates.len()).rev() {
            let j = (order.uniform() * (i + 1) as f64) as usize;
            candidates.swap(i, j.min(i));
        }
        Snapshot { path: path.to_string(), fingerprint: fnv(text.as_bytes()), mols, cell, candidates }
    }

    pub fn len(&self) -> usize {
        self.candidates.len()
    }

    /// Pair `k` of the order: the molecules' indices, their positions (the
    /// second moved to its nearest image of the first) and the centroid
    /// separation.
    pub fn pair(&self, k: usize) -> (usize, usize, Vec<Vec3>, Vec<Vec3>, f64) {
        let (i, j, shift, sep) = self.candidates[k];
        let a = self.mols[i].clone();
        let b: Vec<Vec3> = self.mols[j].iter().map(|p| *p + shift).collect();
        (i, j, a, b, sep)
    }
}

/// The pairs of a simulated liquid: `args` is `name snapshot count`, the
/// snapshot a file `phys-bulk` wrote (`box x y z`, then one `mol` line per
/// molecule with its atoms' positions). Candidates are every pair of molecules
/// whose centroids lie within 10.5 bohr by the nearest image; they are put in
/// a seeded order and taken from the front, so a run resumes, and each is
/// computed exactly as a random pair is and appended to
/// `pairs-<name><tag>.txt` with its centroid separation as its `R`.
///
/// What this is for: the law is fitted to pairs, and a liquid visits
/// arrangements the random draw does not (the first law fitted to random pairs
/// overbound the liquid by 5 kcal/mol a molecule); the pairs a liquid actually
/// has are what its law must be right for.
pub fn run_snapshot(args: &[String], tag: &str) {
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let snapshot = args.get(1).cloned().expect("a snapshot file");
    let count: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(40);
    let grown = args.get(3).map(|s| s == "grown").unwrap_or(false);
    let mut mono = Monomer::load(&name, grown);
    mono.mp2 = args.iter().any(|a| a == "--mp2");
    let snap = Snapshot::load(&snapshot, mono.z.len());
    let out = format!("pairs-{name}{tag}{}.txt", if mono.mp2 { "-mp2" } else { "" });
    let done = done_indices(&out);
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    if done.is_empty() {
        if mono.mp2 {
            writeln!(file, "# MP2 (Hartree-Fock + RI-MP2, frozen cores, counterpoise), pairs of molecules taken from {snapshot}: index centroid_separation_bohr E_int_total E_int_hartree_fock (hartree) | then x y z type for each atom of the first molecule and of the second (bohr)").ok();
        } else {
            writeln!(file, "{}", snapshot_header(&snapshot)).ok();
        }
    }
    println!("{name}: {} candidate pairs in {snapshot} within 10.5 bohr; {} done, {count} wanted", snap.len(), done.len());
    for k in 0..count.min(snap.len()) {
        if done.contains(&k) {
            continue;
        }
        let t = Instant::now();
        let (i, j, a, b, sep) = snap.pair(k);
        let (r, line) = mono.compute(k, &a, &b, sep);
        writeln!(file, "{line}").expect("the line written");
        file.flush().ok();
        println!("pair {k} (molecules {i}, {j}): R {sep:.2} bohr, E_int {:.4} / {:.4} kcal/mol (PBE x / revPBE x), {:.0} s", r.pbe * 627.509474, r.rev * 627.509474, t.elapsed().as_secs_f64());
        std::io::stdout().flush().ok();
    }
}

/// Turn a plan into the tasks still to do. One line each, `#` for comments:
///
/// ```text
/// pair <name> <from> <to> <element|grown> <class> <out>
/// snap <name> <snapshot> <from> <to> <element|grown> <class> <out>
/// ```
///
/// `pair` is the random draw's indices `from..to`; `snap` is the same range of
/// a snapshot's pairs, in its seeded order. Either may end with `weight=W`
/// (the work in one task relative to a pair, default 1) and `mem=GB` (memory a
/// task needs, which otherwise is learned from the first to finish): what the
/// queue sizes machines against. `class` is `cpu` or `gpu` (see
/// `queue`) and `out` the file the lines go to. An index already in `out` is
/// left out, so a plan can be run again after a stop, and two plans that name
/// one file must not overlap in index. Pairs of the random draw are one-to-one
/// with their index only inside a file: a `cpu` and a `gpu` file hold the same
/// draws of the same indices, so a plan that wants *new* pairs from a `cpu`
/// machine starts above the ones the `gpu` file has.
pub fn plan(text: &str) -> Result<Vec<crate::queue::Task>, String> {
    let mut tasks = Vec::new();
    let mut seen_out: std::collections::HashMap<String, (usize, usize)> = std::collections::HashMap::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut w: Vec<&str> = line.split_whitespace().collect();
        let at = |what: &str| format!("plan line {}: {what}", n + 1);
        let mut weight = 1.0;
        let mut need_mb = None;
        while let Some(&last) = w.last() {
            if let Some(v) = last.strip_prefix("weight=") {
                weight = v.parse::<f64>().ok().filter(|x| *x > 0.0).ok_or_else(|| at("weight= wants a positive number"))?;
            } else if let Some(v) = last.strip_prefix("mem=") {
                need_mb = Some((v.parse::<f64>().ok().filter(|x| *x > 0.0).ok_or_else(|| at("mem= wants gigabytes"))? * 1024.0) as u64);
            } else {
                break;
            }
            w.pop();
        }
        let num = |s: &str| s.parse::<usize>().map_err(|_| at(&format!("{s:?} is not a number")));
        let check_basis = |b: &str| if b == "element" || b == "grown" { Ok(()) } else { Err(at("the basis is `element` or `grown`")) };
        let check_class = |c: &str| if c == "cpu" || c == "gpu" { Ok(()) } else { Err(at("the class is `cpu` or `gpu`")) };
        match w.first().copied() {
            Some("pair") if w.len() == 7 => {
                let (name, from, to, basis, class, out) = (w[1], num(w[2])?, num(w[3])?, w[4], w[5], w[6]);
                check_basis(basis)?;
                check_class(class)?;
                if !["water", "methane", "ammonia", "methanol"].contains(&name) {
                    return Err(at(&format!("no molecule {name:?}")));
                }
                if let Some(&(f0, t0)) = seen_out.get(out) {
                    if from < t0 && f0 < to {
                        return Err(at(&format!("{out} is already given indices {f0}..{t0}")));
                    }
                }
                seen_out.insert(out.to_string(), (from, to));
                let done = done_indices(out);
                for k in from..to {
                    if !done.contains(&k) {
                        tasks.push(crate::queue::Task { out: out.to_string(), index: k, class: class.to_string(), spec: format!("pair {name} {basis} {k}"), header: RANDOM_HEADER.to_string(), kind: format!("{name} {basis}"), weight, need_mb });
                    }
                }
            }
            Some("snap") if w.len() == 8 => {
                let (name, path, from, to, basis, class, out) = (w[1], w[2], num(w[3])?, num(w[4])?, w[5], w[6], w[7]);
                check_basis(basis)?;
                check_class(class)?;
                if !["water", "methane", "ammonia", "methanol"].contains(&name) {
                    return Err(at(&format!("no molecule {name:?}")));
                }
                let snap = Snapshot::load(path, atoms_of(name).len());
                if let Some(&(f0, t0)) = seen_out.get(out) {
                    if from < t0 && f0 < to {
                        return Err(at(&format!("{out} is already given indices {f0}..{t0}")));
                    }
                }
                seen_out.insert(out.to_string(), (from, to));
                let done = done_indices(out);
                for k in from..to.min(snap.len()) {
                    if !done.contains(&k) {
                        tasks.push(crate::queue::Task { out: out.to_string(), index: k, class: class.to_string(), spec: format!("snap {name} {basis} {path} {:x} {k}", snap.fingerprint), header: snapshot_header(path), kind: format!("{name} {basis}"), weight, need_mb });
                    }
                }
            }
            _ => return Err(at("expected `pair <name> <from> <to> <basis> <class> <out>` or `snap <name> <snapshot> <from> <to> <basis> <class> <out>`")),
        }
    }
    Ok(tasks)
}

/// Does a queue task: loads each molecule and snapshot once, then computes the
/// pair a spec names. Everything it needs is in the spec and in files the
/// machine holds; it refuses a snapshot that is not the server's.
#[derive(Default)]
pub struct Executor {
    monomers: std::collections::HashMap<(String, bool), Monomer>,
    snapshots: std::collections::HashMap<String, Snapshot>,
}

impl Executor {
    /// `spec` is `pair <name> <basis> <k>` or
    /// `snap <name> <basis> <path> <fingerprint> <k>`; the result is the line
    /// that goes in the pairs file.
    pub fn run(&mut self, spec: &str) -> Result<String, String> {
        let w: Vec<&str> = spec.split_whitespace().collect();
        let bad = || format!("cannot read the task {spec:?}");
        let (name, grown) = match (w.first().copied(), w.get(1), w.get(2)) {
            (Some("pair"), Some(n), Some(b)) | (Some("snap"), Some(n), Some(b)) => (n.to_string(), *b == "grown"),
            _ => return Err(bad()),
        };
        if !std::path::Path::new(&format!("grow-{name}.state")).exists() {
            return Err(format!("no grow-{name}.state on this machine"));
        }
        let mono = self.monomers.entry((name.clone(), grown)).or_insert_with(|| Monomer::load(&name, grown));
        let t = Instant::now();
        match w[0] {
            "pair" => {
                let k: usize = w.get(3).and_then(|s| s.parse().ok()).ok_or_else(bad)?;
                let (a, b, sep) = mono.random_pair(k);
                let (r, line) = mono.compute(k, &a, &b, sep);
                println!("pair {k}: R {sep:.2} bohr, E_int {:.4} / {:.4} kcal/mol, {:.0} s ({} functions)", r.pbe * 627.509474, r.rev * 627.509474, t.elapsed().as_secs_f64(), r.functions);
                Ok(line)
            }
            _ => {
                let (path, finger, k) = match (w.get(3), w.get(4), w.get(5).and_then(|s| s.parse::<usize>().ok())) {
                    (Some(p), Some(f), Some(k)) => (p.to_string(), u64::from_str_radix(f, 16).map_err(|_| bad())?, k),
                    _ => return Err(bad()),
                };
                if !std::path::Path::new(&path).exists() {
                    return Err(format!("no {path} on this machine"));
                }
                let atoms = mono.z.len();
                let snap = self.snapshots.entry(path.clone()).or_insert_with(|| Snapshot::load(&path, atoms));
                if snap.fingerprint != finger {
                    return Err(format!("{path} here is not the server's file (fingerprint {:x}, wanted {finger:x})", snap.fingerprint));
                }
                if k >= snap.len() {
                    return Err(format!("{path} has {} pairs, not {}", snap.len(), k + 1));
                }
                let (i, j, a, b, sep) = snap.pair(k);
                let (r, line) = mono.compute(k, &a, &b, sep);
                println!("pair {k} (molecules {i}, {j}): R {sep:.2} bohr, E_int {:.4} / {:.4} kcal/mol, {:.0} s", r.pbe * 627.509474, r.rev * 627.509474, t.elapsed().as_secs_f64());
                Ok(line)
            }
        }
    }
}

/// The worker's whole `main`: `args` is `host:port [--name N] [--jobs N]
/// [--patience MINUTES] [--mem GB]` (`--mem` offers less than is free, or says
/// how much where it cannot be measured), and the token is `$PHYS_QUEUE_TOKEN` (default
/// `open`). The class is the binary's: `cpu` for `phys-worker`, `gpu` for the
/// GPU one, which has installed its engine before calling this.
pub fn work_main(args: &[String], class: &str) {
    let server = args.first().cloned().unwrap_or_else(|| {
        eprintln!("usage: worker host:port [--name N] [--jobs N] [--patience MINUTES] [--mem GB]");
        std::process::exit(2);
    });
    let flag = |f: &str| args.iter().position(|a| a == f).and_then(|i| args.get(i + 1)).cloned();
    let name = flag("--name").or_else(|| std::env::var("COMPUTERNAME").ok()).or_else(|| std::env::var("HOSTNAME").ok()).unwrap_or_else(|| "unnamed".into());
    let cfg = crate::queue::WorkerConfig {
        server,
        token: std::env::var("PHYS_QUEUE_TOKEN").unwrap_or_else(|_| "open".into()),
        name: name.replace(char::is_whitespace, "-"),
        class: class.to_string(),
        mem_mb: flag("--mem").and_then(|s| s.parse::<f64>().ok()).map(|gb| (gb * 1024.0) as u64),
        max_jobs: flag("--jobs").and_then(|s| s.parse().ok()),
        patience: Duration::from_secs(60 * flag("--patience").and_then(|s| s.parse().ok()).unwrap_or(30)),
    };
    remove_stale_spill();
    println!("{} ({}) working for {}", cfg.name, cfg.class, cfg.server);
    let mut exec = Executor::default();
    let ended = crate::queue::work(&cfg, |spec| exec.run(spec));
    println!("worker stopped: {ended:?}");
    // A wrapper loop restarts a worker that stopped on its job limit and
    // leaves one that was told there is nothing left.
    std::process::exit(match ended {
        crate::queue::Ended::Limit => 10,
        crate::queue::Ended::NothingLeft => 0,
        crate::queue::Ended::ServerGone => 3,
        crate::queue::Ended::TooSmall => 4,
    });
}

/// One line of a pairs file, read back.
struct Recorded {
    index: usize,
    sep: f64,
    pbe: f64,
    rev: f64,
    atoms: Vec<([f64; 3], usize)>,
}

fn parse_recorded(line: &str, atoms: usize) -> Result<Recorded, String> {
    let (head, tail) = line.split_once('|').ok_or("no `|`")?;
    let h: Vec<&str> = head.split_whitespace().collect();
    if h.len() != 4 {
        return Err(format!("{} numbers before the `|`, not 4", h.len()));
    }
    let num = |s: &str| s.parse::<f64>().ok().filter(|x| x.is_finite()).ok_or_else(|| format!("{s:?} is not a finite number"));
    let t: Vec<&str> = tail.split_whitespace().collect();
    if t.len() != 4 * atoms {
        return Err(format!("{} values after the `|`, not {}", t.len(), 4 * atoms));
    }
    let mut at = Vec::new();
    for c in t.chunks(4) {
        at.push(([num(c[0])?, num(c[1])?, num(c[2])?], c[3].parse::<usize>().map_err(|_| format!("{:?} is not a type", c[3]))?));
    }
    Ok(Recorded { index: h[0].parse::<usize>().map_err(|_| "the index is not a number".to_string())?, sep: num(h[1])?, pbe: num(h[2])?, rev: num(h[3])?, atoms: at })
}

impl Monomer {
    /// What a pairs file's lines must satisfy whatever computed them: every
    /// line whole, no index twice, the molecules rigid (the monomer's own
    /// internal distances to 1e-6 bohr, which the 8 decimals written allow),
    /// each atom's type the monomer's, and the separation written the distance
    /// between the molecules' centres of mass. Returns the lines read and what
    /// is wrong, each with its index.
    fn validate(&self, text: &str) -> (Vec<Recorded>, Vec<String>) {
        let n = self.z.len();
        let mut rows = Vec::new();
        let mut bad = Vec::new();
        let mut seen = std::collections::HashSet::new();
        let total: f64 = self.masses.iter().sum();
        let com = |m: &[([f64; 3], usize)]| {
            let mut c = Vec3::ZERO;
            for (a, mass) in m.iter().zip(&self.masses) {
                c = c + Vec3 { x: a.0[0], y: a.0[1], z: a.0[2] }.scale(mass / total);
            }
            c
        };
        for (ln, line) in text.lines().enumerate().filter(|(_, l)| !l.starts_with('#') && !l.trim().is_empty()) {
            let r = match parse_recorded(line, 2 * n) {
                Ok(r) => r,
                Err(e) => {
                    bad.push(format!("line {}: {e}", ln + 1));
                    continue;
                }
            };
            let k = r.index;
            if !seen.insert(k) {
                bad.push(format!("pair {k}: the index appears twice"));
            }
            for (m, name) in [(&r.atoms[..n], "first"), (&r.atoms[n..], "second")] {
                for (i, a) in m.iter().enumerate() {
                    if a.1 != self.types[i] {
                        bad.push(format!("pair {k}: the {name} molecule's atom {i} has type {}, not {}", a.1, self.types[i]));
                    }
                    for j in 0..i {
                        let d = (Vec3 { x: a.0[0], y: a.0[1], z: a.0[2] } - Vec3 { x: m[j].0[0], y: m[j].0[1], z: m[j].0[2] }).norm();
                        let want = (self.body[i] - self.body[j]).norm();
                        if (d - want).abs() > 1e-6 {
                            bad.push(format!("pair {k}: the {name} molecule is not rigid: atoms {j}-{i} are {d:.8} apart, not {want:.8}"));
                        }
                    }
                }
            }
            // Only a random pair's separation is the centres' distance; a
            // snapshot pair's is its centroids'.
            let sep_com = (com(&r.atoms[..n]) - com(&r.atoms[n..])).norm();
            let centroid = |m: &[([f64; 3], usize)]| m.iter().fold(Vec3::ZERO, |c, a| c + Vec3 { x: a.0[0], y: a.0[1], z: a.0[2] }).scale(1.0 / n as f64);
            let sep_centroid = (centroid(&r.atoms[..n]) - centroid(&r.atoms[n..])).norm();
            if (sep_com - r.sep).abs() > 2e-6 && (sep_centroid - r.sep).abs() > 2e-6 {
                bad.push(format!("pair {k}: separation {} on the line, {sep_com:.6} between centres of mass, {sep_centroid:.6} between centroids", r.sep));
            }
            rows.push(r);
        }
        (rows, bad)
    }
}

/// Read a pairs file back and check it against itself and, for a sample of
/// its pairs, against the calculation: `args` is `file name [--sample N]
/// [--seed S] [--indices 124,125] [--tol KCAL] [--grown]`.
///
/// What it is for: the file is the dataset, and it was written over days by
/// processes that were restarted, rebuilt and moved to other machines. The
/// first check is structural (see [`Monomer::validate`]) and needs no solve.
/// The second takes each chosen pair's geometry *as recorded* and solves it
/// again, to say whether the energies on the line are what this binary gets
/// now. Run it with the binary class the file was made with: `phys-recheck` for
/// a CPU file, `phys-recheck-gpu` for a GPU one. The recorded geometry has 8
/// decimals, so a repeat differs from the original by that rounding as well as
/// by anything that changed. Differences are reported as they are; `--tol`
/// (default 0.01 kcal/mol) only decides the exit status.
pub fn recheck(args: &[String]) -> i32 {
    let file = args.first().cloned().expect("a pairs file");
    let name = args.get(1).cloned().expect("a molecule name");
    let flag = |f: &str| args.iter().position(|a| a == f).and_then(|i| args.get(i + 1)).cloned();
    let grown = args.iter().any(|a| a == "--grown");
    let tol: f64 = flag("--tol").and_then(|s| s.parse().ok()).unwrap_or(0.01);
    let mono = Monomer::load(&name, grown);
    let text = std::fs::read_to_string(&file).unwrap_or_else(|_| panic!("no {file}"));
    let (rows, bad) = mono.validate(&text);
    let mut indices: Vec<usize> = rows.iter().map(|r| r.index).collect();
    indices.sort();
    let gaps: Vec<usize> = (0..indices.last().map(|&m| m + 1).unwrap_or(0)).filter(|k| indices.binary_search(k).is_err()).collect();
    println!("{file}: {} pairs read, {} problems{}", rows.len(), bad.len(), if gaps.is_empty() { String::new() } else { format!(", indices missing: {gaps:?}") });
    for b in &bad {
        println!("  PROBLEM {b}");
    }
    let chosen: Vec<usize> = if let Some(list) = flag("--indices") {
        list.split(',').filter_map(|s| s.trim().parse().ok()).collect()
    } else {
        let want: usize = flag("--sample").and_then(|s| s.parse().ok()).unwrap_or(0);
        let seed: u128 = flag("--seed").and_then(|s| s.parse().ok()).unwrap_or(1);
        let mut order = Stream::at(0x6368_6563_6b00, seed, 0, Purpose::Positions);
        let mut pool = indices.clone();
        for i in (1..pool.len()).rev() {
            let j = (order.uniform() * (i + 1) as f64) as usize;
            pool.swap(i, j.min(i));
        }
        pool.truncate(want);
        pool.sort();
        pool
    };
    let n = mono.z.len();
    let (mut worst, mut sum, mut sq, mut count) = (0.0f64, 0.0, 0.0, 0usize);
    let mut over = 0;
    for k in chosen {
        let r = match rows.iter().find(|r| r.index == k) {
            Some(r) => r,
            None => {
                println!("  pair {k}: not in the file");
                over += 1;
                continue;
            }
        };
        let v = |a: &([f64; 3], usize)| Vec3 { x: a.0[0], y: a.0[1], z: a.0[2] };
        let a: Vec<Vec3> = r.atoms[..n].iter().map(v).collect();
        let b: Vec<Vec3> = r.atoms[n..].iter().map(v).collect();
        let t = Instant::now();
        let setup = Setup { z: &mono.z, extra: &mono.extra, grown: mono.grown };
        let again = interact(&setup, k, &a, &b);
        let (dp, dr) = ((again.pbe - r.pbe) * 627.509474, (again.rev - r.rev) * 627.509474);
        for d in [dp, dr] {
            worst = worst.max(d.abs());
            sum += d;
            sq += d * d;
            count += 1;
        }
        if dp.abs() > tol || dr.abs() > tol {
            over += 1;
        }
        println!("  pair {k}: R {:.2} bohr, recorded {:.5} / {:.5}, now {:.5} / {:.5} kcal/mol, difference {dp:+.5} / {dr:+.5}, {:.0} s", r.sep, r.pbe * 627.509474, r.rev * 627.509474, again.pbe * 627.509474, again.rev * 627.509474, t.elapsed().as_secs_f64());
        std::io::stdout().flush().ok();
    }
    if count > 0 {
        println!("{} pairs solved again: largest difference {worst:.5} kcal/mol, mean {:+.5}, rms {:.5}; {over} beyond {tol} kcal/mol", count / 2, sum / count as f64, (sq / count as f64).sqrt());
    }
    if bad.is_empty() && over == 0 { 0 } else { 1 }
}

impl Snapshot {
    /// A cluster of `count` molecules about molecule `centre`: it and the
    /// `count - 1` nearest by centroid (nearest image), each moved to its
    /// image about the centre's so the cluster is in one piece. Returns each
    /// molecule's index in the snapshot and its atoms' positions.
    pub fn cluster(&self, centre: usize, count: usize) -> Vec<(usize, Vec<Vec3>)> {
        let atoms = self.mols[centre].len();
        let centroid = |m: &Vec<Vec3>| m.iter().fold(Vec3::ZERO, |acc, p| acc + *p).scale(1.0 / atoms as f64);
        let c0 = centroid(&self.mols[centre]);
        let image = |mut d: Vec3| {
            d.x -= self.cell[0] * (d.x / self.cell[0]).round();
            d.y -= self.cell[1] * (d.y / self.cell[1]).round();
            d.z -= self.cell[2] * (d.z / self.cell[2]).round();
            d
        };
        let mut others: Vec<(f64, usize, Vec3)> = (0..self.mols.len())
            .filter(|&j| j != centre)
            .map(|j| {
                let raw = centroid(&self.mols[j]) - c0;
                let d = image(raw);
                (d.norm(), j, d - raw)
            })
            .collect();
        others.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut out = vec![(centre, self.mols[centre].clone())];
        for (_, j, shift) in others.into_iter().take(count - 1) {
            out.push((j, self.mols[j].iter().map(|p| *p + shift).collect()));
        }
        out
    }
}

/// What a cluster costs, before anything is solved: the functions in its
/// basis, the auxiliary functions, and the most its three-centre table can
/// hold (every pair of functions against every auxiliary function, in bytes;
/// screening takes some of it away, and a cluster of separate molecules has a
/// lot to take).
pub struct ClusterSize {
    pub functions: usize,
    pub auxiliary: usize,
    pub table_upper_bytes: f64,
}

impl Monomer {
    fn cluster_problem(&self, mols: &[Vec<Vec3>]) -> crate::electrons::scf::Problem {
        let m = mols.len();
        let z: Vec<u32> = (0..m).flat_map(|_| self.z.iter().cloned()).collect();
        let positions: Vec<[f64; 3]> = mols.iter().flatten().map(|v| [v.x, v.y, v.z]).collect();
        let cluster = Molecule { z, positions, charge: 0, unpaired: 0 };
        if self.grown {
            let ex: Vec<Vec<(usize, f64)>> = (0..m).flat_map(|_| self.extra.iter().cloned()).collect();
            cluster.problem_with(Functional::Pbe, Some(&ex))
        } else {
            cluster.problem(Functional::Pbe)
        }
    }

    pub fn cluster_size(&self, mols: &[Vec<Vec3>]) -> ClusterSize {
        let p = self.cluster_problem(mols);
        let n = p.basis.size;
        let aux = p.auxiliary.as_ref().map(|b| b.size).unwrap_or(0);
        ClusterSize { functions: n, auxiliary: aux, table_upper_bytes: (n * (n + 1) / 2) as f64 * aux as f64 * 8.0 }
    }

    /// The counterpoise interaction energy of a cluster of any number of
    /// molecules, as [`interact`] is for two: the cluster, and each molecule
    /// in the whole cluster's basis with the others' nuclei and electrons
    /// removed; one field each, the same energies from it on the fine grid.
    /// The cluster's interaction is the cluster's energy minus the molecules',
    /// subtracted in order, which for two molecules is `interact`'s own sum.
    fn cluster_interaction(&self, mols: &[Vec<Vec3>], first_only: bool, only: Option<usize>, mut log: Option<&mut FieldLog>, mp2_mode: bool) -> Interaction {
        let n = self.z.len();
        let m = mols.len();
        let base = self.cluster_problem(mols);
        let electrons: f64 = self.z.iter().map(|&zz| zz as f64).sum();
        let mut slots: Vec<(Vec<usize>, f64)> = vec![((0..m * n).collect(), m as f64 * electrons)];
        for i in 0..m {
            slots.push(((i * n..(i + 1) * n).collect(), electrons));
        }
        let mut e_pbe = Vec::new();
        let mut e_rev = Vec::new();
        let (mut t_solve, mut t_fine, mut t_swap) = (0.0, 0.0, 0.0);
        let mut iterations = Vec::new();
        for (slot, (keep, ne)) in slots.into_iter().enumerate() {
            // One field asked for: the others are some other machine's.
            if only.is_some_and(|o| o != slot) {
                continue;
            }
            // A field already in the log is not solved again.
            if let Some((pe, re)) = log.as_ref().and_then(|l| l.done.get(&slot).copied()) {
                e_pbe.push(pe);
                e_rev.push(re);
                iterations.push(0);
                println!("  field {} of {} taken from {}", slot + 1, m + 1, log.as_ref().map(|l| l.path.as_str()).unwrap_or(""));
                if only.is_some() {
                    break;
                }
                continue;
            }
            let field_start = Instant::now();
            if mp2_mode {
                // Hartree-Fock and MP2 for this field: the two numbers stored
                // are the Hartree-Fock energy and the MP2 correlation energy.
                let p = base.with_ghosts(&keep, ne);
                let ts = Instant::now();
                let hf = crate::electrons::hf::hartree_fock(&p, 200, 1e-9);
                assert!(hf.converged, "a Hartree-Fock field of the cluster did not converge");
                t_solve += ts.elapsed().as_secs_f64();
                iterations.push(hf.iterations);
                let ts = Instant::now();
                let z_all: Vec<u32> = (0..m).flat_map(|_| self.z.iter().cloned()).collect();
                let corr = crate::electrons::hf::mp2(&p, &hf, crate::electrons::hf::frozen_core(&z_all, &keep));
                t_fine += ts.elapsed().as_secs_f64();
                e_pbe.push(hf.energy);
                e_rev.push(corr.correlation);
                if let Some(l) = log.as_mut() {
                    l.record(slot, hf.energy, corr.correlation, hf.iterations, field_start.elapsed().as_secs_f64());
                }
                println!("  field {} of {} done: {:.0} s so far ({:.0} s Hartree-Fock, {} iterations; {:.0} s MP2)", slot + 1, m + 1, t_solve + t_fine, t_solve, hf.iterations, t_fine);
                std::io::stdout().flush().ok();
                if first_only || only.is_some() {
                    break;
                }
                continue;
            }
            let mut p = base.with_ghosts(&keep, ne);
            p.functional = Functional::PbeXLdaC;
            p.nonlocal = Some(NonlocalSpec::in_the_field(Z_AB_DF1));
            let ts = Instant::now();
            let sol = solve(&p, 200, 1e-10);
            assert!(sol.converged, "a field of the cluster did not converge");
            t_solve += ts.elapsed().as_secs_f64();
            iterations.push(sol.iterations);
            let ts = Instant::now();
            let fine = energy_on_finer_grid(&p, &sol, 50, 12);
            t_fine += ts.elapsed().as_secs_f64();
            let ts = Instant::now();
            let atoms: Vec<([f64; 3], f64)> = p.nuclei.iter().zip(&p.sizes).map(|((_, q), r)| (*q, *r)).collect();
            let grid = crate::electrons::grid::molecular_pruned(&atoms, p.radial, p.theta, p.prune);
            let batches = Batches::new(&p.basis, &grid);
            let (x_pbe, _, _) = exchange_correlation(&p.basis, &batches, Functional::PbeXLdaC, &sol.density_alpha, &sol.density_beta);
            let (x_rev, _, _) = exchange_correlation(&p.basis, &batches, Functional::RevPbeXLdaC, &sol.density_alpha, &sol.density_beta);
            e_pbe.push(fine);
            e_rev.push(fine - x_pbe + x_rev);
            t_swap += ts.elapsed().as_secs_f64();
            if let Some(l) = log.as_mut() {
                l.record(slot, fine, fine - x_pbe + x_rev, *iterations.last().unwrap(), field_start.elapsed().as_secs_f64());
            }
            println!("  field {} of {} done: {:.0} s so far ({:.0} s solving, {} iterations; {:.0} s on the fine grid; {:.0} s swapping exchange)", slot + 1, m + 1, t_solve + t_fine + t_swap, t_solve, iterations.last().unwrap(), t_fine, t_swap);
            std::io::stdout().flush().ok();
            if first_only || only.is_some() {
                break;
            }
        }
        let sub = |e: &[f64]| e[1..].iter().fold(e[0], |acc, x| acc - x);
        Interaction { pbe: sub(&e_pbe), rev: sub(&e_rev), t_solve, t_fine, t_swap, iterations, functions: base.basis.size }
    }
}

/// The electrostatic interaction of two molecules' own charge distributions —
/// the first-order energy of nuclei and electrons as they are in the
/// molecules, before either feels the other — for pairs already computed
/// (`PLAY.md` E8b). Each pair's two monomers are solved in the pair's basis
/// (Hartree-Fock, the ghosts of the counterpoise correction), and the energy is
/// `Z_a Z_b / R` across, the nuclei of each in the other's electrons, and the
/// electrons of each in the other's, with the electron-electron term through
/// the fit as everything else is. A law's electrostatics can then be compared
/// with the density's own, geometry by geometry, apart from its repulsion,
/// dispersion and induction.
///
/// `phys-es-gpu water pairs-water-gpu-bias-mp2.txt [--max-oo 4.6 (A)] [--every n]`
/// appends `k separation E_es` (bohr, hartree) to `es-<name>.txt`, resuming.
pub fn es_main(args: &[String]) {
    use crate::electrons::integrals::one_electron;
    use crate::electrons::scf::Fitted;
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let file = args.get(1).cloned().expect("a pairs file");
    let flag = |f: &str| args.iter().position(|a| a == f).and_then(|i| args.get(i + 1)).cloned();
    let max_oo: f64 = flag("--max-oo").and_then(|v| v.parse().ok()).unwrap_or(4.6);
    let every: usize = flag("--every").and_then(|v| v.parse().ok()).unwrap_or(1);
    let mono = Monomer::load(&name, false);
    let n = mono.z.len();
    let out = format!("es-{name}.txt");
    let done = done_indices(&out);
    let mut outfile = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    let text = std::fs::read_to_string(&file).unwrap_or_else(|_| panic!("no {file}"));
    let mut taken = 0usize;
    for line in text.lines().filter(|l| !l.starts_with('#') && l.contains('|')) {
        let (head, tail) = line.split_once('|').expect("a |");
        let h: Vec<&str> = head.split_whitespace().collect();
        let k: usize = h[0].parse().expect("an index");
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<Vec3> = t.chunks(4).map(|c| Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }).collect();
        let (a, b) = (atoms[..n].to_vec(), atoms[n..].to_vec());
        let oo = (a[0] - b[0]).norm() * 0.529177210903;
        if oo > max_oo || done.contains(&k) {
            continue;
        }
        taken += 1;
        if taken % every != 0 {
            continue;
        }
        let start = Instant::now();
        let mols = vec![a.clone(), b.clone()];
        let base = mono.cluster_problem(&mols);
        let electrons: f64 = mono.z.iter().map(|&zz| zz as f64).sum();
        let nb = base.basis.size;
        let mut dens = Vec::new();
        for i in 0..2 {
            let keep: Vec<usize> = (i * n..(i + 1) * n).collect();
            let p = base.with_ghosts(&keep, electrons);
            let hf = crate::electrons::hf::hartree_fock(&p, 200, 1e-9);
            assert!(hf.converged, "a monomer did not converge");
            let mut d = crate::electrons::linalg::Matrix::zeros(nb);
            for x in 0..nb {
                for y in 0..nb {
                    let s: f64 = (0..hf.occupied).map(|o| hf.c[x * hf.m + o] * hf.c[y * hf.m + o]).sum();
                    d.set(x, y, 2.0 * s);
                }
            }
            dens.push(d);
        }
        let nuc = |m: &[Vec3]| -> Vec<(f64, [f64; 3])> { m.iter().zip(&mono.z).map(|(p, zz)| (*zz as f64, [p.x, p.y, p.z])).collect() };
        let (na, nbn) = (nuc(&a), nuc(&b));
        let mut e = 0.0;
        for (za, pa) in &na {
            for (zb, pb) in &nbn {
                e += za * zb / ((pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2) + (pa[2] - pb[2]).powi(2)).sqrt();
            }
        }
        let (s, _, v_b) = one_electron(&base.basis, &nbn);
        let (_, _, v_a) = one_electron(&base.basis, &na);
        for x in 0..nb * nb {
            e += dens[0].a[x] * v_b.a[x] + dens[1].a[x] * v_a.a[x];
        }
        let aux = base.auxiliary.as_ref().expect("an auxiliary basis");
        let fit = Fitted::new(&base.basis, aux, &s);
        let (j, _) = fit.coulomb_and_energy(&dens[0]);
        for x in 0..nb * nb {
            e += dens[1].a[x] * j.a[x];
        }
        let sep = (a[0] - b[0]).norm();
        writeln!(outfile, "{k} {sep:.6} {e:.10e}").expect("written");
        outfile.flush().ok();
        println!("pair {k}: O-O {oo:.2} A, electrostatic {:.3} kcal/mol, {:.0} s", e * 627.509474, start.elapsed().as_secs_f64());
        std::io::stdout().flush().ok();
    }
}

/// The slope of a pair's MP2 interaction energy along the line between the
/// molecules' centres of mass, orientations fixed — the force a liquid's
/// pressure is made of (`PLAY.md` E8b). A liquid's virial is the sum over pairs
/// of the centres' separation times the force on the centres, which for rigid
/// molecules is minus this derivative; an energy fit constrains the wall's
/// height and not its slope, and a pressure of 1956 bar short of zero is a
/// slope that is a few per cent out. Each pair is computed with the second
/// molecule moved `delta` bohr (0.1 by default) farther and nearer along the
/// line, and the two total counterpoise energies written as
/// `k delta E_plus E_minus` (hartree), for the fit to difference.
///
/// `phys-es-gpu water pairs-file --deriv [--max-oo 4.6 (A)] [--every n] [--delta d]`
/// appends to `deriv-<name>.txt`, resuming.
pub fn deriv_main(args: &[String]) {
    let name = args.first().cloned().unwrap_or_else(|| "water".into());
    let file = args.get(1).cloned().expect("a pairs file");
    let flag = |f: &str| args.iter().position(|a| a == f).and_then(|i| args.get(i + 1)).cloned();
    let max_oo: f64 = flag("--max-oo").and_then(|v| v.parse().ok()).unwrap_or(4.6);
    let every: usize = flag("--every").and_then(|v| v.parse().ok()).unwrap_or(1);
    let delta: f64 = flag("--delta").and_then(|v| v.parse().ok()).unwrap_or(0.1);
    let mut mono = Monomer::load(&name, false);
    mono.mp2 = true;
    let n = mono.z.len();
    let out = format!("deriv-{name}.txt");
    let done = done_indices(&out);
    let mut outfile = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    let text = std::fs::read_to_string(&file).unwrap_or_else(|_| panic!("no {file}"));
    let total: f64 = mono.masses.iter().sum();
    let com = |m: &[Vec3]| -> Vec3 { m.iter().zip(&mono.masses).fold(Vec3::ZERO, |acc, (p, w)| acc + p.scale(*w / total)) };
    let mut taken = 0usize;
    for line in text.lines().filter(|l| !l.starts_with('#') && l.contains('|')) {
        let (head, tail) = line.split_once('|').expect("a |");
        let h: Vec<&str> = head.split_whitespace().collect();
        let k: usize = h[0].parse().expect("an index");
        let t: Vec<&str> = tail.split_whitespace().collect();
        let atoms: Vec<Vec3> = t.chunks(4).map(|c| Vec3 { x: c[0].parse().unwrap(), y: c[1].parse().unwrap(), z: c[2].parse().unwrap() }).collect();
        let (a, b) = (atoms[..n].to_vec(), atoms[n..].to_vec());
        let oo = (a[0] - b[0]).norm() * 0.529177210903;
        if oo > max_oo || done.contains(&k) {
            continue;
        }
        taken += 1;
        if taken % every != 0 {
            continue;
        }
        let start = Instant::now();
        let axis = (com(&b) - com(&a)).unit();
        let mut e = [0.0f64; 2];
        for (i, sgn) in [1.0, -1.0].iter().enumerate() {
            let shifted: Vec<Vec3> = b.iter().map(|p| *p + axis.scale(sgn * delta)).collect();
            let r = mono.cluster_interaction(&[a.clone(), shifted], false, None, None, true);
            e[i] = r.pbe + r.rev;
        }
        writeln!(outfile, "{k} {delta} {:.10e} {:.10e}", e[0], e[1]).expect("written");
        outfile.flush().ok();
        println!("pair {k}: O-O {oo:.2} A, slope {:.3} kcal/mol/bohr, {:.0} s", (e[0] - e[1]) / (2.0 * delta) * 627.509474, start.elapsed().as_secs_f64());
        std::io::stdout().flush().ok();
    }
}

/// The fields of a cluster calculation that are finished, kept in the cluster's
/// own file as they come: one line `field <slot> <PBE-exchange energy>
/// <revPBE-exchange energy> <iterations> <seconds>`, slot 0 the cluster and slot
/// `i + 1` molecule `i` in the cluster's basis. A field takes an hour or more
/// at six molecules, so what is interrupted loses at most the field it was in:
/// a run that finds a slot here does not solve it again.
pub struct FieldLog {
    pub path: String,
    pub done: std::collections::HashMap<usize, (f64, f64)>,
}

impl FieldLog {
    pub fn load(path: &str) -> FieldLog {
        let mut done = std::collections::HashMap::new();
        for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            if w.len() >= 4 && w[0] == "field" {
                if let (Ok(slot), Ok(p), Ok(r)) = (w[1].parse::<usize>(), w[2].parse::<f64>(), w[3].parse::<f64>()) {
                    done.insert(slot, (p, r));
                }
            }
        }
        FieldLog { path: path.to_string(), done }
    }

    pub fn record(&mut self, slot: usize, pbe: f64, rev: f64, iterations: usize, seconds: f64) {
        // Rust prints a float in the shortest form that reads back exactly.
        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(&self.path).expect("the cluster file");
        writeln!(f, "field {slot} {pbe:e} {rev:e} {iterations} {seconds:.0}").expect("the field written");
        f.sync_all().ok();
        self.done.insert(slot, (pbe, rev));
    }
}

/// Is there a process with this id? `None` where that cannot be told.
fn process_alive(pid: u32) -> Option<bool> {
    #[cfg(target_os = "linux")]
    {
        Some(std::path::Path::new(&format!("/proc/{pid}")).exists())
    }
    #[cfg(windows)]
    {
        let out = std::process::Command::new("tasklist").args(["/FI", &format!("PID eq {pid}"), "/NH"]).output().ok()?;
        Some(String::from_utf8_lossy(&out.stdout).split_whitespace().any(|w| w == pid.to_string()))
    }
    #[cfg(not(any(target_os = "linux", windows)))]
    {
        let _ = pid;
        None
    }
}

/// Delete the spill files of processes that are no longer there. A table is
/// removed when the process that wrote it drops it, but a process that was
/// killed leaves 30 GB behind, and a run being restarted is the one that has
/// to find the room. Only `phys-fit-<pid>-<n>.bin` files in the spill
/// directory (`PHYS_SPILL_DIR`) whose process is known not to exist are
/// touched; where that cannot be told nothing is.
pub fn remove_stale_spill() {
    let Some(to) = crate::electrons::scf::SpillTo::from_environment() else { return };
    remove_stale_spill_in(&to.dir, process_alive);
}

/// [`remove_stale_spill`] in `dir`, with `alive` saying whether a process id is
/// in use (`None` for cannot tell).
pub fn remove_stale_spill_in(dir: &std::path::Path, alive: impl Fn(u32) -> Option<bool>) {
    let Ok(entries) = std::fs::read_dir(dir) else { return };
    for e in entries.flatten() {
        let name = e.file_name().to_string_lossy().to_string();
        let Some(rest) = name.strip_prefix("phys-fit-").and_then(|r| r.strip_suffix(".bin")) else { continue };
        let Some(pid) = rest.split('-').next().and_then(|p| p.parse::<u32>().ok()) else { continue };
        if pid != std::process::id() && alive(pid) == Some(false) {
            let size = e.metadata().map(|m| m.len()).unwrap_or(0);
            if std::fs::remove_file(e.path()).is_ok() {
                println!("removed {name}, left by process {pid} ({:.1} GB)", size as f64 / 1e9);
            }
        }
    }
}

/// The key a cluster file's result line is filed under: `pair i j`, `field k`
/// or `cluster`.
fn cluster_line_key(line: &str) -> Option<String> {
    let w: Vec<&str> = line.split_whitespace().collect();
    match w.first().copied() {
        Some("pair") if w.len() >= 3 => Some(format!("pair {} {}", w[1], w[2])),
        Some("field") if w.len() >= 2 => Some(format!("field {}", w[1])),
        Some("cluster") => Some("cluster".to_string()),
        _ => None,
    }
}

/// Add to the cluster file `main` every pair, field and cluster result that
/// `others` hold and it does not. Each of `others` must have been begun on the
/// same cluster (its `# mol` lines are `geometry`); a result that `main`
/// already has is not replaced. Returns how many lines were added. This is how
/// the pieces of a cluster calculation done on several machines come together.
pub fn merge_cluster_files(main: &str, others: &[String], geometry: &[String]) -> usize {
    let mut have: std::collections::HashSet<String> = std::fs::read_to_string(main).unwrap_or_default().lines().filter_map(cluster_line_key).collect();
    let mut added = Vec::new();
    for other in others {
        let text = std::fs::read_to_string(other).unwrap_or_else(|_| panic!("no {other}"));
        let kept: Vec<&str> = text.lines().filter(|l| l.starts_with("# mol ")).collect();
        assert!(kept.len() == geometry.len() && kept.iter().zip(geometry).all(|(a, b)| a == b), "{other} was begun on a different cluster from {main}");
        for line in text.lines() {
            if let Some(key) = cluster_line_key(line) {
                if have.insert(key) {
                    added.push(line.to_string());
                }
            }
        }
    }
    if !added.is_empty() {
        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(main).expect("the cluster file");
        for l in &added {
            writeln!(f, "{l}").expect("a line merged");
        }
        f.sync_all().ok();
    }
    added.len()
}

/// `phys-cluster name snapshot [--centre I] [--size N] [--law FILE] [--dry]
/// [--grown]`: a cluster of `N` molecules (default 6) about molecule `I` of a
/// snapshot, and the question of whether a law that sums pairs gets its
/// energy right.
///
/// It computes the cluster's counterpoise interaction energy, and the same
/// for each of its `N(N-1)/2` pairs, then reports the cluster's energy against
/// the **sum of the pairs'** (the difference is everything that is not
/// additive: three-body polarisation and the rest) and, given a law, against
/// that law's sum over the same pairs. Each piece is appended to
/// `cluster-<name>-<N>-c<I>.txt` as it is finished and a run carries on from
/// what is there. `--dry` prints what the cluster would cost and solves nothing;
/// `--first-field` solves only the cluster's own field and says how long it took.
/// `--method mp2` does it all with Hartree-Fock and RI-MP2 (`electrons::hf`,
/// the engine's own table, no functional) into `cluster-mp2-...` files, the two
/// numbers on each line being the Hartree-Fock energy and the correlation energy.
///
/// **Across machines.** Everything in a cluster after the table is built is
/// independent: its fields (slot 0 the cluster, slot `i + 1` molecule `i` in the
/// cluster's basis) and its pairs. `--field K` does one field and `--pair I,J`
/// one pair, each writing its line to `--out FILE` (default the cluster's own),
/// and `--merge a.txt b.txt` adds what other machines wrote to this one's
/// file before the usual run reports from it. A file is merged only into one
/// begun on the same cluster.
///
/// **Interrupted?** Run the same command again. Every pair and every field of
/// the cluster is on its line in the file the moment it is finished, so a run
/// carries on from there and loses at most the one field it was in (the tables
/// are rebuilt, about five minutes at six molecules); a file begun on another
/// cluster is refused. The spill of a killed run is deleted at the start.
pub fn cluster_main(args: &[String]) {
    let name = args.first().cloned().expect("a molecule name");
    let snapshot = args.get(1).cloned().expect("a snapshot file");
    let flag = |f: &str| args.iter().position(|a| a == f).and_then(|i| args.get(i + 1)).cloned();
    let centre: usize = flag("--centre").and_then(|s| s.parse().ok()).unwrap_or(0);
    let size: usize = flag("--size").and_then(|s| s.parse().ok()).unwrap_or(6);
    let grown = args.iter().any(|a| a == "--grown");
    let dry = args.iter().any(|a| a == "--dry");
    let mp2_mode = flag("--method").map(|m| m == "mp2").unwrap_or(false);
    let mono = Monomer::load(&name, grown);
    let snap = Snapshot::load(&snapshot, mono.z.len());
    // `--check-pair K`: pair K of the snapshot's order through the cluster's
    // own code, to be set beside the line `run_snapshot` wrote for it.
    if let Some(k) = flag("--check-pair").and_then(|s| s.parse::<usize>().ok()) {
        let (i, j, a, b, sep) = snap.pair(k);
        let r = mono.cluster_interaction(&[a.clone(), b.clone()], false, None, None, false);
        let _ = sep;
        println!("pair {k} (molecules {i}, {j}) through the cluster code: {:.10e} / {:.10e} hartree", r.pbe, r.rev);
        let file = flag("--against").unwrap_or_else(|| format!("pairs-{name}-liq-gpu.txt"));
        let line = std::fs::read_to_string(&file).ok().and_then(|t| t.lines().find(|l| l.split_whitespace().next() == Some(&k.to_string())).map(|l| l.to_string()));
        match line.as_deref().map(|l| l.split_whitespace().collect::<Vec<_>>()) {
            Some(w) if w.len() > 4 => println!("  recorded in {file}: {} / {} hartree; differences {:+.2e} / {:+.2e}", w[2], w[3], r.pbe - w[2].parse::<f64>().unwrap_or(f64::NAN), r.rev - w[3].parse::<f64>().unwrap_or(f64::NAN)),
            _ => println!("  no line {k} in {file}"),
        }
        return;
    }
    let mols = snap.cluster(centre, size);
    let ids: Vec<usize> = mols.iter().map(|(i, _)| *i).collect();
    let geom: Vec<Vec<Vec3>> = mols.into_iter().map(|(_, m)| m).collect();
    let cs = mono.cluster_size(&geom);
    let two = mono.cluster_size(&geom[..2]);
    println!("{name} cluster of {size} about molecule {centre} of {snapshot}: molecules {ids:?}");
    println!("  the cluster: {} functions, {} auxiliary, a table of at most {:.1} GB; a pair: {} functions, at most {:.2} GB", cs.functions, cs.auxiliary, cs.table_upper_bytes / 1e9, two.functions, two.table_upper_bytes / 1e9);
    if dry {
        return;
    }
    remove_stale_spill();
    // `--first-field`: only the cluster's own field, to time it.
    if args.iter().any(|a| a == "--first-field") {
        let t = Instant::now();
        mono.cluster_interaction(&geom, true, None, None, mp2_mode);
        println!("the first field of the cluster of {size}: {:.0} s", t.elapsed().as_secs_f64());
        return;
    }
    let out = flag("--out").unwrap_or_else(|| format!("cluster-{}{name}-{size}-c{centre}.txt", if mp2_mode { "mp2-" } else { "" }));
    let have = std::fs::read_to_string(&out).unwrap_or_default();
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    let geometry: Vec<String> = ids
        .iter()
        .zip(&geom)
        .map(|(i, m)| {
            let mut line = format!("# mol {i}");
            for p in m {
                line += &format!(" {:.8} {:.8} {:.8}", p.x, p.y, p.z);
            }
            line
        })
        .collect();
    if have.is_empty() {
        let what = if mp2_mode { "Hartree-Fock energy then MP2 correlation energy (RI, frozen cores)" } else { "PBE-exchange form then revPBE-exchange" };
        writeln!(file, "# cluster of {size} about molecule {centre} of {snapshot}; molecules {ids:?}; energies in hartree, {what}").ok();
        for line in &geometry {
            writeln!(file, "{line}").ok();
        }
    } else {
        // Energies kept from an earlier run belong to the cluster that run
        // had; a different snapshot would give a different one under the same
        // name.
        let kept: Vec<&str> = have.lines().filter(|l| l.starts_with("# mol ")).collect();
        assert!(kept.len() == geometry.len() && kept.iter().zip(&geometry).all(|(a, b)| a == b), "{out} was begun on a different cluster from this snapshot gives; move it aside or use another --centre");
    }
    // `--merge a.txt b.txt ...`: results other machines wrote for this cluster.
    if let Some(i) = args.iter().position(|a| a == "--merge") {
        let others: Vec<String> = args[i + 1..].iter().take_while(|a| !a.starts_with("--")).cloned().collect();
        let n = merge_cluster_files(&out, &others, &geometry);
        println!("  {n} results added to {out} from {others:?}");
    }
    let have = std::fs::read_to_string(&out).unwrap_or_default();
    let mut fields = FieldLog::load(&out);
    let read = |key: &str| -> Option<(f64, f64)> {
        have.lines().find(|l| l.starts_with(key)).and_then(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            Some((w[w.len() - 2].parse().ok()?, w[w.len() - 1].parse().ok()?))
        })
    };
    let law = flag("--law").map(|f| {
        let text = std::fs::read_to_string(&f).unwrap_or_else(|_| panic!("no {f}"));
        (crate::liquid::SiteSite::from_text(&text).expect("a readable law"), crate::liquid::SiteSite::bisector_from_text(&text))
    });
    let law_energy = |a: &[Vec3], b: &[Vec3]| -> Option<f64> {
        let (law, bis) = law.as_ref()?;
        let sites = |m: &[Vec3]| -> Vec<(Vec3, usize)> {
            let atoms: Vec<(Vec3, usize)> = m.iter().cloned().zip(mono.types.iter().cloned()).collect();
            match bis {
                Some((t, d)) => crate::liquid::with_bisector_site(&atoms, *t, *d),
                None => atoms,
            }
        };
        Some(crate::liquid::pair_energy(law, &crate::liquid::PairEnergy { a: sites(a), b: sites(b), energy: 0.0 }))
    };
    // `--field K`: one field, for a machine that is doing only that.
    if let Some(k) = flag("--field").and_then(|s| s.parse::<usize>().ok()) {
        assert!(k <= size, "a cluster of {size} has fields 0 to {size}");
        if fields.done.contains_key(&k) {
            println!("  field {k} is already in {out}");
        } else {
            let t = Instant::now();
            mono.cluster_interaction(&geom, false, Some(k), Some(&mut fields), mp2_mode);
            println!("  field {k} of the cluster of {size}: {:.0} s, written to {out}", t.elapsed().as_secs_f64());
        }
        return;
    }
    // `--pair I,J`: one pair of the cluster, likewise.
    if let Some(ij) = flag("--pair") {
        let v: Vec<usize> = ij.split(',').filter_map(|s| s.trim().parse().ok()).collect();
        assert!(v.len() == 2 && v[0] < v[1] && v[1] < size, "--pair wants I,J with I < J < {size}");
        let (i, j) = (v[0], v[1]);
        if read(&format!("pair {i} {j} ")).is_some() {
            println!("  pair {i} {j} is already in {out}");
        } else {
            let t = Instant::now();
            let r = mono.cluster_interaction(&[geom[i].clone(), geom[j].clone()], false, None, None, mp2_mode);
            writeln!(file, "pair {i} {j} {} {} {:.10e} {:.10e}", ids[i], ids[j], r.pbe, r.rev).ok();
            file.flush().ok();
            println!("  pair {i} {j}: {:.4} / {:.4} kcal/mol, {:.0} s, written to {out}", r.pbe * 627.509474, r.rev * 627.509474, t.elapsed().as_secs_f64());
        }
        return;
    }
    let mut pair_sum = (0.0, 0.0);
    let mut law_sum = 0.0;
    for i in 0..size {
        for j in i + 1..size {
            let key = format!("pair {i} {j} ");
            let (p, r) = match read(&key) {
                Some(v) => v,
                None => {
                    let t = Instant::now();
                    let r = mono.cluster_interaction(&[geom[i].clone(), geom[j].clone()], false, None, None, mp2_mode);
                    let line = format!("pair {i} {j} {} {} {:.10e} {:.10e}", ids[i], ids[j], r.pbe, r.rev);
                    writeln!(file, "{line}").ok();
                    file.flush().ok();
                    println!("  pair {i} {j}: {:.4} / {:.4} kcal/mol, {:.0} s", r.pbe * 627.509474, r.rev * 627.509474, t.elapsed().as_secs_f64());
                    (r.pbe, r.rev)
                }
            };
            pair_sum.0 += p;
            pair_sum.1 += r;
            if let Some(e) = law_energy(&geom[i], &geom[j]) {
                law_sum += e;
            }
        }
    }
    let head = |v: (f64, f64)| if mp2_mode { v.0 + v.1 } else { v.0 };
    println!("  the pairs' sum ({}): {:.4} / {:.4} kcal/mol{}; the law's sum: {}", if mp2_mode { "HF / MP2 correlation" } else { "DFT" }, pair_sum.0 * 627.509474, pair_sum.1 * 627.509474, if mp2_mode { format!(" = {:.4}", head(pair_sum) * 627.509474) } else { String::new() }, if law.is_some() { format!("{:.4} kcal/mol", law_sum * 627.509474) } else { "no law given".into() });
    let whole = match read("cluster ") {
        Some(v) => v,
        None => {
            let t = Instant::now();
            let r = mono.cluster_interaction(&geom, false, None, Some(&mut fields), mp2_mode);
            writeln!(file, "cluster {size} {:.10e} {:.10e}", r.pbe, r.rev).ok();
            file.flush().ok();
            println!("  the cluster: {:.4} / {:.4} kcal/mol, {:.0} s ({} functions)", r.pbe * 627.509474, r.rev * 627.509474, t.elapsed().as_secs_f64(), r.functions);
            (r.pbe, r.rev)
        }
    };
    let k = 627.509474;
    println!("{size} molecules: cluster {:.4} kcal/mol ({}), the {} pairs' sum {:.4}, so what is not additive is {:+.4} kcal/mol ({:+.2} a molecule)", head(whole) * k, if mp2_mode { "HF + MP2" } else { "PBE x" }, if mp2_mode { "MP2" } else { "DFT" }, head(pair_sum) * k, (head(whole) - head(pair_sum)) * k, (head(whole) - head(pair_sum)) * k / size as f64);
    if mp2_mode {
        println!("  of which Hartree-Fock {:+.4} and MP2 correlation {:+.4} kcal/mol", (whole.0 - pair_sum.0) * k, (whole.1 - pair_sum.1) * k);
    }
    if law.is_some() && !mp2_mode {
        println!("  the law's pair sum is {:.4}: {:+.4} kcal/mol from the cluster ({:+.2} a molecule), of which {:+.4} is the law's pairs and {:+.4} the non-additivity", law_sum * k, (law_sum - whole.0) * k, (law_sum - whole.0) * k / size as f64, (law_sum - pair_sum.0) * k, (pair_sum.0 - whole.0) * k);
    }
}
