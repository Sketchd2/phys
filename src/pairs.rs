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
        Monomer { name: name.to_string(), z, types, extra, body, contact, masses, grown }
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

    /// The interaction of the pair `a`, `b`, and the line that records it.
    fn compute(&self, k: usize, a: &[Vec3], b: &[Vec3], sep: f64) -> (Interaction, String) {
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
    let mono = Monomer::load(&name, grown);
    let out = format!("pairs-{name}{tag}.txt");
    let done = done_indices(&out);
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    if done.is_empty() {
        writeln!(file, "{RANDOM_HEADER}").ok();
    }
    println!("{name}: {} atoms, types {:?}, {} basis; {} pairs done, {count} wanted", mono.z.len(), mono.types, if grown { "grown" } else { "per-element" }, done.len());
    for k in 0..count {
        if done.contains(&k) {
            continue;
        }
        let t = Instant::now();
        let (a, b, sep) = mono.random_pair(k);
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
        Snapshot { path: path.to_string(), fingerprint: fnv(text.as_bytes()), mols, candidates }
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
    let mono = Monomer::load(&name, grown);
    let snap = Snapshot::load(&snapshot, mono.z.len());
    let out = format!("pairs-{name}{tag}.txt");
    let done = done_indices(&out);
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(&out).expect("the output file");
    if done.is_empty() {
        writeln!(file, "{}", snapshot_header(&snapshot)).ok();
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
