//! A basis grown for the molecule.
//!
//! The functions a molecule's electrons are expanded in are chosen for *that
//! molecule* (`PLAY.md` E5b, the owner's decision), not per element. A set
//! derived on the free atom's energy (E4) is both larger than a bond needs and,
//! on water, not enough for its length: 214 functions left the O-H bond 0.10%
//! long, because a lone atom's energy cares about the cusp at its nucleus and
//! not about the space between two atoms.
//!
//! So the basis starts from each atom's own contracted orbitals — what the free
//! atom is, and its contracted response to a field — and grows where the
//! molecule needs it:
//!
//! 1. Every candidate function (a primitive at some atom, angular momentum and
//!    exponent) gets a *cheap estimate* of what adding it would do, from the
//!    current density alone: the second-order energy lowering of mixing it into
//!    the occupied orbitals, read off one Fock matrix built over the current
//!    functions and every candidate together. Nothing is solved again.
//! 2. What matters is the shape, so each candidate's estimate is taken with
//!    each internal coordinate (a bond length, a bond angle) stretched and
//!    compressed; the difference is the extra pull the candidate puts on that
//!    coordinate, and over the coordinate's own stiffness, how far it would move
//!    it.
//! 3. The candidates predicted to move the shape most are added, and the
//!    molecule is relaxed again — the *full check*, which says what they
//!    actually did.
//! 4. It stops when one more step moves every length by less than its
//!    tolerance and every angle by less than its tolerance, twice running, and
//!    what the remaining candidates are predicted to do is inside it too: the
//!    tolerance is on the next step, not on any reference (the owner's rule),
//!    with care for a result that swings rather than settles.
//!
//! **Where candidates come from is derived too.** Each atom has a ladder of
//! exponents per angular momentum, seeded from its element's own derived set
//! (E4) and one angular momentum above it; a pick at either end of a ladder
//! extends that end, and a pick on the highest ladder opens the next one —
//! the same extend-while-it-helps rule E4 used, applied to the molecule.

use super::basis::{Basis, Shell};
use super::functional::Functional;
use super::integrals::{eri_three_contracted, one_electron, pairs, Pair};
use super::linalg::{eigh, generalised, orthogonaliser, Matrix};
use super::molecule::{element_basis, relax_in, Molecule};
use super::scf::{exchange_correlation, occupy, parallel_interleaved, solve, Batches, Problem, Solution};

/// A ladder of exponents one atom may draw from at one angular momentum:
/// `first * ratio^k` for `k` in `lo..=hi`.
#[derive(Debug, Clone, PartialEq)]
pub struct Ladder {
    pub atom: usize,
    pub l: usize,
    pub first: f64,
    pub ratio: f64,
    pub lo: i32,
    pub hi: i32,
}

impl Ladder {
    pub fn exponent(&self, k: i32) -> f64 {
        self.first * self.ratio.powi(k)
    }
}

/// One candidate: rung `k` of ladder `ladder`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Rung {
    pub ladder: usize,
    pub k: i32,
}

/// The ladders for a molecule, seeded from each element's derived set: per
/// angular momentum that element has, its own exponents' span and spacing; and
/// one angular momentum above its highest, on the highest one's span.
pub fn ladders(mol: &Molecule, f: Functional) -> Vec<Ladder> {
    let mut out = Vec::new();
    for (atom, &z) in mol.z.iter().enumerate() {
        let b = element_basis(z, f);
        let mut by_l: Vec<Vec<f64>> = Vec::new();
        let mut push = |l: usize, e: f64| {
            if by_l.len() <= l {
                by_l.resize(l + 1, Vec::new());
            }
            by_l[l].push(e);
        };
        for (l, exps) in &b.shells {
            exps.iter().for_each(|e| push(*l, *e));
        }
        for (l, exps) in b.free.iter().chain(&b.polarisation) {
            exps.iter().for_each(|e| push(*l, *e));
        }
        for (l, terms) in &b.contracted {
            terms.iter().for_each(|t| push(*l, t.0));
        }
        let mut last: Option<(f64, f64, i32)> = None;
        for (l, exps) in by_l.iter_mut().enumerate() {
            exps.sort_by(|a, b| a.total_cmp(b));
            exps.dedup_by(|a, b| (*a / *b - 1.0).abs() < 1e-9);
            if exps.is_empty() {
                continue;
            }
            let (lo, hi) = (exps[0], exps[exps.len() - 1]);
            let (ratio, count) = if exps.len() >= 2 {
                let n = exps.len() as i32 - 1;
                ((hi / lo).powf(1.0 / n as f64), n)
            } else {
                (last.map(|x| x.1).unwrap_or(2.0), 0)
            };
            out.push(Ladder { atom, l, first: lo, ratio, lo: 0, hi: count });
            last = Some((lo, ratio, count));
        }
        if let Some((first, ratio, count)) = last {
            out.push(Ladder { atom, l: by_l.len(), first, ratio, lo: 0, hi: count });
        }
    }
    out
}

/// The extra functions per atom that a set of rungs adds, `(l, exponent)`.
pub fn extras(mol: &Molecule, ladders: &[Ladder], chosen: &[Rung]) -> Vec<Vec<(usize, f64)>> {
    let mut out = vec![Vec::new(); mol.z.len()];
    for r in chosen {
        let lad = &ladders[r.ladder];
        out[lad.atom].push((lad.l, lad.exponent(r.k)));
    }
    out
}

/// The Fock matrices over the current basis and a set of candidate shells
/// together, built from the current density: what the cheap estimate reads.
struct UnionFock {
    s: Matrix,
    f: [Matrix; 2],
    n_cur: usize,
    /// Each candidate's function indices in the union.
    ranges: Vec<std::ops::Range<usize>>,
}

/// Build it: one-electron terms; the Coulomb potential of the density the
/// solve already fitted, between every pair (which needs no new fit — it is
/// the fitted density's own potential); exchange-correlation over the union
/// with the density embedded and the candidates empty.
fn union_fock(problem: &Problem, solution: &Solution, cands: &[Shell]) -> UnionFock {
    let n_cur = problem.basis.size;
    let mut shells = problem.basis.shells.clone();
    shells.extend(cands.iter().cloned());
    let union = Basis::new(shells);
    let n = union.size;
    let ranges: Vec<std::ops::Range<usize>> = (problem.basis.shells.len()..union.shells.len()).map(|i| union.offsets[i]..union.offsets[i] + union.shells[i].size()).collect();
    let (s, t, v) = one_electron(&union, &problem.nuclei);
    let aux = problem.auxiliary.as_ref().expect("a grown basis is solved with fitted Coulomb");
    let c = &solution.fitted;
    let unit = Shell::unit();
    let kets: Vec<Vec<Pair>> = aux.shells.iter().map(|p| pairs(p, &unit, 0)).collect();
    let sh = &union.shells;
    let shell_pairs: Vec<(usize, usize)> = (0..sh.len()).flat_map(|a| (0..=a).map(move |b| (a, b))).filter(|&(a, b)| {
        let r2: f64 = (0..3).map(|k| (sh[a].centre[k] - sh[b].centre[k]).powi(2)).sum();
        sh[a].exponents.iter().any(|x| sh[b].exponents.iter().any(|y| (-(x * y) / (x + y) * r2).exp() > 1e-14))
    }).collect();
    let job = |idx: &mut dyn Iterator<Item = usize>| {
        let mut out: Vec<(usize, usize, Vec<f64>)> = Vec::new();
        for (a, b) in idx.map(|i| shell_pairs[i]) {
            let bra = pairs(&sh[a], &sh[b], 0);
            let mut block = vec![0.0; sh[a].size() * sh[b].size()];
            for (ip, p) in aux.shells.iter().enumerate() {
                let cp = &c[aux.offsets[ip]..aux.offsets[ip] + p.size()];
                let add = eri_three_contracted(&bra, &kets[ip], sh[a].l, sh[b].l, p.l, cp);
                for (x, y) in block.iter_mut().zip(add) {
                    *x += y;
                }
            }
            out.push((a, b, block));
        }
        out
    };
    let mut j = Matrix::zeros(n);
    for part in parallel_interleaved(shell_pairs.len(), &job) {
        for (a, b, block) in part {
            let nb = sh[b].size();
            for i in 0..sh[a].size() {
                for k in 0..nb {
                    let (r, col) = (union.offsets[a] + i, union.offsets[b] + k);
                    j.set(r, col, block[i * nb + k]);
                    j.set(col, r, block[i * nb + k]);
                }
            }
        }
    }
    let embed = |d: &Matrix| {
        let mut e = Matrix::zeros(n);
        for i in 0..n_cur {
            e.a[i * n..i * n + n_cur].copy_from_slice(&d.a[i * n_cur..i * n_cur + n_cur]);
        }
        e
    };
    let atoms: Vec<([f64; 3], f64)> = problem.nuclei.iter().zip(&problem.sizes).map(|((_, p), r)| (*p, *r)).collect();
    let grid = super::grid::molecular_pruned(&atoms, problem.radial, problem.theta, problem.prune);
    let batches = Batches::new(&union, &grid);
    let (_, va, vb) = exchange_correlation(&union, &batches, problem.functional, &embed(&solution.density_alpha), &embed(&solution.density_beta));
    let mut fa = t.clone();
    for k in 0..n * n {
        fa.a[k] += v.a[k] + j.a[k];
    }
    let mut fb = fa.clone();
    for k in 0..n * n {
        fa.a[k] += va.a[k];
        fb.a[k] += vb.a[k];
    }
    UnionFock { s, f: [fa, fb], n_cur, ranges }
}

/// Below this fraction of its own norm left once the current functions are
/// taken out, a candidate is the current space over again, and its estimate is
/// round-off amplified by normalising what is left.
const DEPENDENT: f64 = 1e-2;

/// The cheap estimate: for each candidate shell, the second-order lowering of
/// the energy if it were mixed into the occupied orbitals —
/// `-sum_spin sum_(x, i) g_xi^2 / (e_x - e_i)`, with `x` running over the
/// candidate's functions orthogonalised to the current space and diagonal in
/// the Fock matrix, `g` their coupling to occupied orbital `i`. Hartree, each
/// negative or zero.
///
/// The denominator is held at least at the current gap between the highest
/// occupied and lowest empty level: a candidate whose own level falls near or
/// below the occupied ones is one the perturbation cannot describe, and its
/// estimate would otherwise diverge. The gap is the molecule's own.
pub fn estimate(problem: &Problem, solution: &Solution, cands: &[Shell]) -> Vec<f64> {
    let u = union_fock(problem, solution, cands);
    let n = u.s.n;
    let nc = u.n_cur;
    // The current block's overlap, and its pseudo-inverse for projecting out
    // the current space.
    let mut s_cur = Matrix::zeros(nc);
    for i in 0..nc {
        s_cur.a[i * nc..i * nc + nc].copy_from_slice(&u.s.a[i * n..i * n + nc]);
    }
    let (sv, svec) = eigh(&s_cur);
    let top = sv.last().cloned().unwrap_or(1.0);
    let s_inv = {
        let mut m = Matrix::zeros(nc);
        for k in 0..nc {
            if sv[k] <= 1e-10 * top {
                continue;
            }
            for i in 0..nc {
                let a = svec.get(i, k) / sv[k];
                for j in 0..nc {
                    m.a[i * nc + j] += a * svec.get(j, k);
                }
            }
        }
        m
    };
    let (x, m) = orthogonaliser(&s_cur, 1e-8);
    let counts = [problem.alpha, problem.beta];
    // Per spin: occupied orbitals over the union, their levels, the gap, and
    // F times them.
    let mut spins = Vec::new();
    for s in 0..2 {
        if counts[s] <= 0.0 {
            continue;
        }
        let mut f_cur = Matrix::zeros(nc);
        for i in 0..nc {
            f_cur.a[i * nc..i * nc + nc].copy_from_slice(&u.f[s].a[i * n..i * n + nc]);
        }
        let (levels, c) = generalised(&f_cur, &x, m);
        let occ = occupy(&levels, counts[s]);
        let homo = levels.iter().zip(&occ).filter(|(_, o)| **o > 0.0).map(|(e, _)| *e).fold(f64::NEG_INFINITY, f64::max);
        let lumo = levels.iter().zip(&occ).filter(|(_, o)| **o < 1.0).map(|(e, _)| *e).fold(f64::INFINITY, f64::min);
        let gap = if lumo.is_finite() { (lumo - homo).max(0.0) } else { 0.0 };
        let occupied: Vec<usize> = (0..m).filter(|&k| occ[k] > 0.0).collect();
        // F (union) times each occupied orbital embedded in the union.
        let fc: Vec<Vec<f64>> = occupied.iter().map(|&k| {
            (0..n).map(|r| (0..nc).map(|q| u.f[s].a[r * n + q] * c[q * m + k]).sum()).collect()
        }).collect();
        spins.push((s, occupied.iter().map(|&k| (levels[k], occ[k])).collect::<Vec<(f64, f64)>>(), fc, gap));
    }
    u.ranges.iter().map(|range| {
        let sz = range.len();
        // Each candidate function with the current space taken out:
        // y = e_j - S_cur^-1 S_(cur, j).
        let proj: Vec<Vec<f64>> = range.clone().map(|j| (0..nc).map(|q| (0..nc).map(|r| s_inv.a[q * nc + r] * u.s.a[r * n + j]).sum()).collect()).collect();
        let mut mm = Matrix::zeros(sz);
        for (a, ja) in range.clone().enumerate() {
            for (b, jb) in range.clone().enumerate() {
                let corr: f64 = (0..nc).map(|q| u.s.a[ja * n + q] * proj[b][q]).sum();
                mm.set(a, b, u.s.a[ja * n + jb] - corr);
            }
        }
        let (w, wv) = eigh(&mm);
        let kept: Vec<usize> = (0..sz).filter(|&k| w[k] > DEPENDENT).collect();
        if kept.is_empty() {
            return 0.0;
        }
        // Orthonormal remainders, as vectors over the union.
        let ys: Vec<Vec<f64>> = kept.iter().map(|&k| {
            let mut y = vec![0.0; n];
            let norm = 1.0 / w[k].sqrt();
            for (a, ja) in range.clone().enumerate() {
                let coef = wv.get(a, k) * norm;
                y[ja] += coef;
                for q in 0..nc {
                    y[q] -= coef * proj[a][q];
                }
            }
            y
        }).collect();
        let mut total = 0.0;
        for (s, occs, fc, gap) in &spins {
            let f = &u.f[*s];
            let nk = ys.len();
            // F within the remainders, diagonalised.
            let fy: Vec<Vec<f64>> = ys.iter().map(|y| (0..n).map(|r| (0..n).map(|q| f.a[r * n + q] * y[q]).sum()).collect()).collect();
            let mut fyy = Matrix::zeros(nk);
            for a in 0..nk {
                for b in 0..nk {
                    fyy.set(a, b, ys[a].iter().zip(&fy[b]).map(|(p, q)| p * q).sum());
                }
            }
            let (ex, exv) = eigh(&fyy);
            for xi in 0..nk {
                for (i, (ei, oi)) in occs.iter().enumerate() {
                    let g: f64 = (0..nk).map(|a| exv.get(a, xi) * ys[a].iter().zip(&fc[i]).map(|(p, q)| p * q).sum::<f64>()).sum();
                    let denom = (ex[xi] - ei).max(*gap).max(1e-12);
                    total -= oi * g * g / denom;
                }
            }
        }
        total
    }).collect()
}

/// An internal coordinate: a bond's length, or the angle at `b` between `a`
/// and `c`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Coordinate {
    Bond(usize, usize),
    Angle(usize, usize, usize),
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

impl Coordinate {
    /// Its value: bohr for a bond, radians for an angle.
    pub fn value(&self, p: &[[f64; 3]]) -> f64 {
        match *self {
            Coordinate::Bond(a, b) => norm(sub(p[a], p[b])),
            Coordinate::Angle(a, b, c) => {
                let (u, v) = (sub(p[a], p[b]), sub(p[c], p[b]));
                ((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (norm(u) * norm(v))).clamp(-1.0, 1.0).acos()
            }
        }
    }

    /// How far it may move in one step and still count as settled: the
    /// owner's 0.1% of a length, 0.1 degree of an angle.
    pub fn tolerance(&self, p: &[[f64; 3]]) -> f64 {
        match self {
            Coordinate::Bond(..) => LENGTH_TOLERANCE * self.value(p),
            Coordinate::Angle(..) => ANGLE_TOLERANCE_DEGREES.to_radians(),
        }
    }

    /// The finite-difference step its stiffness and pull are read with: five
    /// times its tolerance — 0.5% of a bond, as the prototype used — small
    /// enough to stay where the energy is quadratic, and large enough that the
    /// energy it changes (about 2e-5 hartree on water's stretch) stands far
    /// above the solve's convergence of 1e-10.
    fn step(&self, p: &[[f64; 3]]) -> f64 {
        5.0 * self.tolerance(p)
    }
}

/// Positions with coordinate `which` moved by `d` and, to first order, no
/// other: along `B^+ e_which`, where `B` is every coordinate's derivative with
/// respect to every atom's position and `B^+` its pseudo-inverse. Moving a
/// bond's two ends apart along the bond moves the atoms they share with other
/// bonds, and so every angle and length around them: the pull read off such a
/// step is the pull on all of them at once. Through the pseudo-inverse the
/// step is the smallest motion that changes this coordinate alone — or, in a
/// ring, where the coordinates are not independent, the smallest that changes
/// it with the rest as little as they allow — and it never translates or
/// turns the molecule as a whole.
pub fn displaced(coords: &[Coordinate], p: &[[f64; 3]], which: usize, d: f64) -> Vec<[f64; 3]> {
    let n = 3 * p.len();
    let nq = coords.len();
    // B by central differences: the coordinates are smooth and this is exact
    // to the step's square.
    let h = 1e-6;
    let mut b = vec![0.0; nq * n];
    for x in 0..n {
        let mut plus = p.to_vec();
        let mut minus = p.to_vec();
        plus[x / 3][x % 3] += h;
        minus[x / 3][x % 3] -= h;
        for (qi, c) in coords.iter().enumerate() {
            b[qi * n + x] = (c.value(&plus) - c.value(&minus)) / (2.0 * h);
        }
    }
    // dx = B^T (B B^T)^+ e_which.
    let mut g = Matrix::zeros(nq);
    for i in 0..nq {
        for j in 0..nq {
            g.a[i * nq + j] = (0..n).map(|x| b[i * n + x] * b[j * n + x]).sum();
        }
    }
    let (vals, vecs) = eigh(&g);
    let top = vals.iter().cloned().fold(0.0, f64::max);
    let mut y = vec![0.0; nq];
    for k in 0..nq {
        if vals[k] <= 1e-10 * top {
            continue;
        }
        let w = vecs.get(which, k) / vals[k];
        for i in 0..nq {
            y[i] += w * vecs.get(i, k);
        }
    }
    let mut q = p.to_vec();
    for x in 0..n {
        let dx: f64 = (0..nq).map(|i| b[i * n + x] * y[i]).sum();
        q[x / 3][x % 3] += d * dx;
    }
    q
}

/// The owner's tolerance on a bond length, relative: 0.1%.
pub const LENGTH_TOLERANCE: f64 = 1e-3;
/// The owner's tolerance on a bond angle: 0.1 degree.
pub const ANGLE_TOLERANCE_DEGREES: f64 = 0.1;

/// Bond lengths and the angles between bonds sharing an atom.
pub fn coordinates(bonds: &[(usize, usize)]) -> Vec<Coordinate> {
    let mut out: Vec<Coordinate> = bonds.iter().map(|&(a, b)| Coordinate::Bond(a, b)).collect();
    let n = bonds.iter().map(|&(a, b)| a.max(b) + 1).max().unwrap_or(0);
    for centre in 0..n {
        let nb: Vec<usize> = bonds.iter().filter_map(|&(a, b)| if a == centre { Some(b) } else if b == centre { Some(a) } else { None }).collect();
        for i in 0..nb.len() {
            for j in i + 1..nb.len() {
                out.push(Coordinate::Angle(nb[i], centre, nb[j]));
            }
        }
    }
    out
}

/// One round of growth, kept for the comparison of estimate and check.
#[derive(Debug, Clone)]
pub struct Round {
    pub functions: usize,
    pub energy: f64,
    /// Each coordinate's value after relaxing in this round's basis.
    pub values: Vec<f64>,
    /// What was added after this round, and the move of each coordinate the
    /// estimate predicted for it.
    pub added: Vec<Rung>,
    pub predicted: Vec<f64>,
    /// The largest of the remaining candidates' predicted moves, each over its
    /// coordinate's tolerance, summed: what is left to do, as the estimate
    /// sees it.
    pub left: f64,
}

/// What growing produced.
#[derive(Debug, Clone)]
pub struct Growth {
    pub molecule: Molecule,
    pub ladders: Vec<Ladder>,
    pub chosen: Vec<Rung>,
    pub rounds: Vec<Round>,
    pub coordinates: Vec<Coordinate>,
    pub converged: bool,
}

/// The candidate shells a round considers: every rung of every ladder not yet
/// chosen.
fn candidates(mol: &Molecule, ladders: &[Ladder], chosen: &[Rung]) -> (Vec<Rung>, Vec<Shell>) {
    let mut rungs = Vec::new();
    let mut shells = Vec::new();
    for (i, lad) in ladders.iter().enumerate() {
        for k in lad.lo..=lad.hi {
            let r = Rung { ladder: i, k };
            if chosen.contains(&r) {
                continue;
            }
            rungs.push(r);
            shells.push(Shell::primitive(mol.positions[lad.atom], lad.l, lad.exponent(k)));
        }
    }
    (rungs, shells)
}

/// The candidates a round would consider, exposed for the tests.
pub fn candidates_for_test(mol: &Molecule, ladders: &[Ladder], chosen: &[Rung]) -> (Vec<Rung>, Vec<Shell>) {
    candidates(mol, ladders, chosen)
}

/// Grow a basis for `start`, whose bonds are `bonds`, until the shape stops
/// moving. `max_rounds` bounds it; `picks` is how many distinct candidates
/// (each with its symmetry partners) are added a round.
pub fn grow(start: &Molecule, bonds: &[(usize, usize)], f: Functional, max_rounds: usize, picks: usize) -> Growth {
    grow_reporting(start, bonds, f, max_rounds, picks, &mut |_, _| {})
}

/// As [`grow`], handing each round to `report` as it finishes — a run on a
/// large molecule takes hours, and what it has found should not wait for the
/// end.
pub fn grow_reporting(start: &Molecule, bonds: &[(usize, usize)], f: Functional, max_rounds: usize, picks: usize, report: &mut dyn FnMut(&Round, &[Coordinate])) -> Growth {
    let coords = coordinates(bonds);
    let mut lads = ladders(start, f);
    let mut chosen: Vec<Rung> = Vec::new();
    let mut rounds: Vec<Round> = Vec::new();
    let mut mol = start.clone();
    let mut converged = false;
    for _ in 0..max_rounds {
        let extra = extras(&mol, &lads, &chosen);
        // The full check: relaxed in this basis, stopping well inside the
        // smallest tolerance.
        let tol_min = coords.iter().filter(|c| matches!(c, Coordinate::Bond(..))).map(|c| c.tolerance(&mol.positions)).fold(f64::INFINITY, f64::min);
        let (relaxed, _, _) = relax_in(&mol, f, Some(&extra), 60, (RELAXED_FORCE, tol_min / 20.0));
        mol = relaxed;
        let p0 = mol.problem_with(f, Some(&extra));
        let s0 = solve(&p0, 200, 1e-10);
        let values: Vec<f64> = coords.iter().map(|c| c.value(&mol.positions)).collect();
        let (rungs, shells) = candidates(&mol, &lads, &chosen);
        // Each candidate's predicted move of each coordinate.
        let mut pred = vec![vec![0.0; coords.len()]; rungs.len()];
        for (ci, c) in coords.iter().enumerate() {
            let h = c.step(&mol.positions);
            let mut e = [0.0; 2];
            let mut dq = [0.0; 2];
            let mut gain = [Vec::new(), Vec::new()];
            for (si, sgn) in [1.0, -1.0].iter().enumerate() {
                let moved = Molecule { positions: displaced(&coords, &mol.positions, ci, sgn * h), ..mol.clone() };
                // What the coordinate actually moved: all of `h` when the
                // coordinates are independent, less in a ring.
                dq[si] = c.value(&moved.positions) - values[ci];
                let mut p = moved.problem_with(f, Some(&extra));
                p.guess = Some((s0.density_alpha.clone(), s0.density_beta.clone()));
                let s = solve(&p, 200, 1e-10);
                e[si] = s.energy;
                let moved_shells: Vec<Shell> = shells.iter().zip(&rungs).map(|(sh, r)| Shell { centre: moved.positions[lads[r.ladder].atom], ..sh.clone() }).collect();
                gain[si] = estimate(&p, &s, &moved_shells);
            }
            // Stiffness and pull from the two unequal steps it actually took.
            let (hp, hm) = (dq[0], dq[1]);
            let k = 2.0 * ((e[0] - s0.energy) / hp - (e[1] - s0.energy) / hm) / (hp - hm);
            for (ri, row) in pred.iter_mut().enumerate() {
                let pull = (gain[0][ri] - gain[1][ri]) / (hp - hm);
                row[ci] = if k > 0.0 { -pull / k } else { 0.0 };
            }
        }
        let tols: Vec<f64> = coords.iter().map(|c| c.tolerance(&mol.positions)).collect();
        let score = |row: &Vec<f64>| row.iter().zip(&tols).map(|(p, t)| p.abs() / t).fold(0.0, f64::max);
        let mut order: Vec<usize> = (0..rungs.len()).collect();
        order.sort_by(|&a, &b| score(&pred[b]).total_cmp(&score(&pred[a])).then(rungs[a].cmp(&rungs[b])));
        let left: f64 = order.iter().map(|&i| score(&pred[i])).sum();
        let functions = p0.basis.size;
        // Settled: two steps running inside every tolerance, the swing over
        // the last three inside it too, and what is predicted to be left
        // inside it.
        let settled_steps = rounds.len() >= 2 && {
            let n = rounds.len();
            let step_ok = |a: &Vec<f64>, b: &Vec<f64>| a.iter().zip(b).zip(&tols).all(|((x, y), t)| (x - y).abs() < *t);
            step_ok(&values, &rounds[n - 1].values) && step_ok(&rounds[n - 1].values, &rounds[n - 2].values) && (0..coords.len()).all(|i| {
                let w = [values[i], rounds[n - 1].values[i], rounds[n - 2].values[i]];
                w.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - w.iter().cloned().fold(f64::INFINITY, f64::min) < tols[i]
            })
        };
        if settled_steps && left < 1.0 {
            rounds.push(Round { functions, energy: s0.energy, values, added: Vec::new(), predicted: vec![0.0; coords.len()], left });
            report(rounds.last().expect("just pushed"), &coords);
            converged = true;
            break;
        }
        // The best `picks` distinct candidates, each with any symmetry partner
        // (same element, angular momentum and exponent, and the same predicted
        // effect to a part in a thousand).
        let key = |i: usize| {
            let lad = &lads[rungs[i].ladder];
            (mol.z[lad.atom], lad.l, lad.exponent(rungs[i].k), score(&pred[i]))
        };
        let same = |a: (u32, usize, f64, f64), b: (u32, usize, f64, f64)| a.0 == b.0 && a.1 == b.1 && (a.2 / b.2 - 1.0).abs() < 1e-9 && (a.3 - b.3).abs() <= 1e-3 * a.3.max(b.3);
        let mut groups: Vec<(u32, usize, f64, f64)> = Vec::new();
        let mut added = Vec::new();
        for &i in &order {
            let k = key(i);
            if k.3 <= 0.0 {
                break;
            }
            if !groups.iter().any(|g| same(*g, k)) {
                if groups.len() >= picks {
                    break;
                }
                groups.push(k);
            }
            added.push(rungs[i]);
        }
        let predicted: Vec<f64> = (0..coords.len()).map(|c| added.iter().map(|r| pred[rungs.iter().position(|x| x == r).expect("added from the list")][c]).sum()).collect();
        rounds.push(Round { functions, energy: s0.energy, values, added: added.clone(), predicted, left });
        report(rounds.last().expect("just pushed"), &coords);
        if added.is_empty() {
            break;
        }
        // Extend any ladder picked at an end, and open the next angular
        // momentum above any atom's highest that was picked.
        for r in &added {
            let lad = lads[r.ladder].clone();
            if r.k == lad.lo {
                lads[r.ladder].lo -= 1;
            }
            if r.k == lad.hi {
                lads[r.ladder].hi += 1;
            }
            let top = lads.iter().filter(|x| x.atom == lad.atom).map(|x| x.l).max().unwrap_or(0);
            if lad.l == top {
                lads.push(Ladder { l: top + 1, ..lad.clone() });
            }
        }
        chosen.extend(added);
    }
    Growth { molecule: mol, ladders: lads, chosen, rounds, coordinates: coords, converged }
}

/// The force at which a relaxation inside growth counts as done, hartree per
/// bohr. At a stiffness of about half a hartree per square bohr (water's
/// stretch) it leaves a length 2e-5 bohr from its minimum, about a hundredth
/// of the tolerance on water's bond (1.8e-3 bohr).
const RELAXED_FORCE: f64 = 1e-5;
