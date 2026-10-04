//! A molecule's density shared among its atoms — `docs/PLAY.md` Phase 6, E6.
//!
//! # Why
//!
//! E6's site-site potential (the route the owner had built beside route A's
//! fitted one, to be compared on water) needs what each atom of a molecule
//! is like *inside* it: its charge, its size, how strongly it disperses. None
//! of that is an element's property; it is the molecule's density, shared.
//!
//! # How
//!
//! Hirshfeld's stockholder partition: a point's density goes to each atom in
//! proportion to that atom's own free density there,
//! `w_A(r) = rho_A(|r - R_A|) / sum_B rho_B(|r - R_B|)`. The free densities
//! are the engine's own: each element solved alone in its ground-state spin
//! with its derived basis, spherically averaged (an open-shell atom's
//! density is not quite round), and tabulated once against radius.
//!
//! From the weights: each atom's charge, `Z_A - sum w_A n`; its in-molecule
//! volume, `sum w_A n r_A^3`; and its share of the molecule's dispersion.
//! For that last the engine's own non-local functional (E7) is split rather
//! than a published table read: its `C6` between two molecules is a double
//! sum over their points (`vdw::c6`), and weighting each point by the atom it
//! belongs to divides that sum exactly into atom-pair parts, which add back
//! to the molecule's `C6`.

use super::functional::Functional;
use super::molecule::{element_basis, Molecule};
use super::scf::{solve, Batches};
use super::vdw::{density_and_gradient_on, Site};

/// A free atom's spherically averaged density against radius: `ln rho` at
/// radii spaced evenly in `ln r`.
#[derive(Debug, Clone)]
pub struct FreeAtom {
    pub z: u32,
    ln_r0: f64,
    step: f64,
    ln_rho: Vec<f64>,
}

impl FreeAtom {
    /// Solve the element alone, in its ground-state spin, and average its
    /// density over directions at radii from 1e-4 to 25 bohr.
    pub fn new(z: u32, f: Functional) -> FreeAtom {
        let unpaired = element_basis(z, f).unpaired;
        let atom = Molecule { z: vec![z], positions: vec![[0.0; 3]], charge: 0, unpaired };
        let p = atom.problem(f);
        let s = solve(&p, 300, 1e-10);
        let mut d = s.density_alpha.clone();
        for k in 0..d.a.len() {
            d.a[k] += s.density_beta.a[k];
        }
        let (dirs, wts) = super::grid::sphere(12);
        let (r0, r1, n) = (1e-4f64, 25.0f64, 400usize);
        let step = (r1 / r0).ln() / (n - 1) as f64;
        let nb = p.basis.size;
        let mut vals = vec![0.0; nb];
        let ln_rho = (0..n).map(|i| {
            let r = r0 * (step * i as f64).exp();
            let mut avg = 0.0;
            for (u, w) in dirs.iter().zip(&wts) {
                super::values::at(&p.basis, [u[0] * r, u[1] * r, u[2] * r], &mut vals, None);
                let mut rho = 0.0;
                for a in 0..nb {
                    let x: f64 = (0..nb).map(|b| d.a[a * nb + b] * vals[b]).sum();
                    rho += vals[a] * x;
                }
                avg += w * rho;
            }
            (avg / (4.0 * std::f64::consts::PI)).max(1e-300).ln()
        }).collect();
        FreeAtom { z, ln_r0: r0.ln(), step, ln_rho }
    }

    /// `rho(r)`, interpolating `ln rho` linearly in `ln r`; beyond the table,
    /// the last two points' exponential continued.
    pub fn density(&self, r: f64) -> f64 {
        let x = (r.max(1e-12).ln() - self.ln_r0) / self.step;
        let n = self.ln_rho.len();
        let i = (x.floor().max(0.0) as usize).min(n - 2);
        let t = x - i as f64;
        (self.ln_rho[i] * (1.0 - t) + self.ln_rho[i + 1] * t).exp()
    }
}

/// What the partition gives for one molecule.
#[derive(Debug, Clone)]
pub struct Partition {
    /// Each atom's net charge, elementary charges.
    pub charge: Vec<f64>,
    /// Each atom's in-molecule `<r^3>`, bohr^3.
    pub volume: Vec<f64>,
    /// The molecule's dipole from its density and nuclei, e bohr.
    pub dipole: [f64; 3],
    /// For each grid point kept, each atom's share of it: what an atom-pair
    /// `C6` split needs.
    pub shares: Vec<Vec<f64>>,
    /// The points those shares are of, as the non-local term sees them.
    pub sites: Vec<Site>,
}

/// Partition a solved molecule's density over a grid among its atoms, with
/// `free` holding a free atom for every element present. `z_ab` is the
/// non-local functional's constant, for the points' `q0`.
pub fn partition(mol: &Molecule, problem: &super::scf::Problem, solution: &super::scf::Solution, free: &[FreeAtom], radial: usize, theta: usize, z_ab: f64) -> Partition {
    let atoms: Vec<([f64; 3], f64)> = problem.nuclei.iter().zip(&problem.sizes).map(|((_, p), r)| (*p, *r)).collect();
    let grid = super::grid::molecular_pruned(&atoms, radial, theta, false);
    let batches = Batches::new(&problem.basis, &grid);
    let mut d = solution.density_alpha.clone();
    for k in 0..d.a.len() {
        d.a[k] += solution.density_beta.a[k];
    }
    let dens = density_and_gradient_on(&problem.basis, &batches, &d, grid.points.len());
    let na = mol.z.len();
    let free_of: Vec<&FreeAtom> = mol.z.iter().map(|z| free.iter().find(|f| f.z == *z).expect("a free atom for every element")).collect();
    let mut population = vec![0.0; na];
    let mut volume = vec![0.0; na];
    let mut moment = [0.0; 3];
    let mut shares = Vec::new();
    let mut sites = Vec::new();
    for ((p, w), (n, g)) in grid.points.iter().zip(&grid.weights).zip(&dens) {
        if *n <= 0.0 {
            continue;
        }
        let dist: Vec<f64> = mol.positions.iter().map(|c| ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt()).collect();
        let pro: Vec<f64> = dist.iter().zip(&free_of).map(|(r, fa)| fa.density(*r)).collect();
        let total: f64 = pro.iter().sum();
        let share: Vec<f64> = if total > 0.0 { pro.iter().map(|x| x / total).collect() } else { vec![1.0 / na as f64; na] };
        let wn = w * n;
        for a in 0..na {
            population[a] += share[a] * wn;
            volume[a] += share[a] * wn * dist[a].powi(3);
        }
        for k in 0..3 {
            moment[k] -= wn * p[k];
        }
        if wn.abs() >= super::vdw::FLOOR {
            shares.push(share);
            sites.push(Site { r: *p, wn, q: super::vdw::q0(*n, g[0] * g[0] + g[1] * g[1] + g[2] * g[2], z_ab) });
        }
    }
    for (z, c) in mol.z.iter().zip(&mol.positions) {
        for k in 0..3 {
            moment[k] += *z as f64 * c[k];
        }
    }
    let charge = mol.z.iter().zip(&population).map(|(z, pop)| *z as f64 - pop).collect();
    Partition { charge, volume, dipole: moment, shares, sites }
}

/// The non-local functional's `C6` between two molecules, split into atom
/// pairs: `[a][b]` is atom `a` of the first with atom `b` of the second.
/// They add to `vdw::c6` of the two molecules' points.
pub fn atom_pair_c6(first: &Partition, second: &Partition) -> Vec<Vec<f64>> {
    let (na, nb) = (first.charge.len(), second.charge.len());
    let job = |idx: &mut dyn Iterator<Item = usize>| {
        idx.map(|i| {
            let s = &first.sites[i];
            let q2 = s.q * s.q;
            let mut row = vec![0.0; na * nb];
            for (t, share_b) in second.sites.iter().zip(&second.shares) {
                let p2 = t.q * t.q;
                let g = s.wn * t.wn / (q2 * p2 * (q2 + p2));
                for a in 0..na {
                    let ga = g * first.shares[i][a];
                    for b in 0..nb {
                        row[a * nb + b] += ga * share_b[b];
                    }
                }
            }
            (i, row)
        }).collect::<Vec<_>>()
    };
    let mut rows = vec![vec![0.0; na * nb]; first.sites.len()];
    for part in super::scf::parallel_interleaved(first.sites.len(), &job) {
        for (i, row) in part {
            rows[i] = row;
        }
    }
    let mut out = vec![vec![0.0; nb]; na];
    for row in rows {
        for a in 0..na {
            for b in 0..nb {
                out[a][b] += super::vdw::ASYMPTOTE_C * row[a * nb + b];
            }
        }
    }
    out
}
