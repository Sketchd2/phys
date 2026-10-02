//! A molecule: nuclei, electrons, and the shape the electrons give it.
//!
//! Everything a calculation on a molecule needs is assembled here from the
//! elements' own derived bases (`element`), each derived once per process and
//! kept, and the shape is relaxed by following the forces (`gradient`) until
//! they vanish — which is what replaces a conformer stated from bond lengths
//! and VSEPR angles.

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use super::basis::{Basis, Shell};
use super::element::{derive, ElementBasis};
use super::functional::Functional;
use super::gradient::gradient;
use super::linalg::Matrix;
use super::scf::{solve, Problem, Solution};

/// Bohr per metre.
const BOHR_PER_METRE: f64 = 1.0 / 0.529177210903e-10;

/// An element's derived basis, derived the first time any molecule needs it.
pub fn element_basis(z: u32, f: Functional) -> ElementBasis {
    static CACHE: OnceLock<Mutex<HashMap<(u32, u8), ElementBasis>>> = OnceLock::new();
    let key = (z, f as u8);
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(b) = cache.lock().expect("element cache").get(&key) {
        return b.clone();
    }
    let b = derive(z, f, 1e-5);
    cache.lock().expect("element cache").insert(key, b.clone());
    b
}

/// Nuclei and electrons.
#[derive(Debug, Clone, PartialEq)]
pub struct Molecule {
    pub z: Vec<u32>,
    /// Bohr.
    pub positions: Vec<[f64; 3]>,
    pub charge: i32,
    pub unpaired: u32,
}

impl Molecule {
    /// From a chemical arrangement, with the conformer `chem::geometry::embed`
    /// gives as the starting shape.
    pub fn from_arrangement(arr: &crate::chem::arrange::Arrangement) -> Molecule {
        let pos = crate::chem::geometry::embed(arr);
        let z: Vec<u32> = arr.atoms.iter().map(|e| e.z() as u32).collect();
        let electrons: i64 = z.iter().map(|&x| x as i64).sum::<i64>() - arr.charge as i64;
        Molecule {
            z,
            positions: pos.iter().map(|p| [p.x * BOHR_PER_METRE, p.y * BOHR_PER_METRE, p.z * BOHR_PER_METRE]).collect(),
            charge: arr.charge as i32,
            unpaired: (electrons.rem_euclid(2)) as u32,
        }
    }

    pub fn electrons(&self) -> f64 {
        self.z.iter().map(|&x| x as f64).sum::<f64>() - self.charge as f64
    }

    /// The problem to solve for it: derived bases, an auxiliary set capped at
    /// the molecule's highest angular momentum, grid scales from each atom's
    /// covalent radius, and a starting density made of the free atoms' own.
    pub fn problem(&self, f: Functional) -> Problem {
        let bases: Vec<ElementBasis> = self.z.iter().map(|&z| element_basis(z, f)).collect();
        let lmax = bases.iter().map(|b| b.lmax()).max().unwrap_or(0);
        let mut shells: Vec<Shell> = Vec::new();
        let mut aux: Vec<Shell> = Vec::new();
        let mut sizes = Vec::new();
        // (offset of each contracted function, its l, alpha and beta electrons)
        let mut guess_entries: Vec<(usize, usize, f64, f64)> = Vec::new();
        let mut offset = 0;
        for (b, p) in bases.iter().zip(&self.positions) {
            let these = b.shells_at(*p);
            for (k, sh) in these.iter().enumerate() {
                if k < b.contracted.len() {
                    let (oa, ob) = b.occupation[k];
                    guess_entries.push((offset, sh.l, oa, ob));
                }
                offset += sh.size();
            }
            shells.extend(these);
            aux.extend(b.auxiliary_at(*p, lmax));
            let r = crate::chem::elements::Element(b.z as u8).covalent_radius().unwrap_or(1.0e-10) * BOHR_PER_METRE;
            sizes.push(r.max(0.5));
        }
        let basis = Basis::new(shells);
        let n = basis.size;
        let ne = self.electrons();
        let (alpha, beta) = ((ne + self.unpaired as f64) / 2.0, (ne - self.unpaired as f64) / 2.0);
        // The free atoms' densities, spin-averaged for a closed shell so the
        // guess does not start the molecule off polarised. A level is shared
        // over its Cartesian components; for d that is six rather than five,
        // which a guess can afford.
        let mut da = Matrix::zeros(n);
        let mut db = Matrix::zeros(n);
        let neutral: f64 = self.z.iter().map(|&x| x as f64).sum();
        let scale = if neutral > 0.0 { ne / neutral } else { 1.0 };
        for (off, l, oa, ob) in guess_entries {
            let comps = (l + 1) * (l + 2) / 2;
            let (a, b) = if self.unpaired == 0 { ((oa + ob) / 2.0, (oa + ob) / 2.0) } else { (oa, ob) };
            for c in 0..comps {
                da.set(off + c, off + c, a * scale / comps as f64);
                db.set(off + c, off + c, b * scale / comps as f64);
            }
        }
        Problem {
            basis,
            nuclei: self.z.iter().zip(&self.positions).map(|(&z, p)| (z as f64, *p)).collect(),
            sizes,
            alpha,
            beta,
            functional: f,
            radial: 70,
            theta: 16,
            auxiliary: Some(Basis::new(aux)),
            prune: true,
            guess: Some((da, db)),
        }
    }

    /// Energy, its gradient, and the solution behind them.
    pub fn energy_and_gradient(&self, f: Functional) -> (Solution, Vec<[f64; 3]>) {
        let p = self.problem(f);
        let s = solve(&p, 200, 1e-10);
        let g = gradient(&p, &s);
        (s, g)
    }
}

/// One step of a relaxation, kept for diagnosis.
#[derive(Debug, Clone)]
pub struct Step {
    pub energy: f64,
    pub largest_force: f64,
    pub largest_move: f64,
}

/// Largest force, hartree/bohr, at which a shape counts as relaxed.
pub const FORCE_CONVERGED: f64 = 4.5e-4;
/// Largest move, bohr, at which a shape counts as relaxed.
pub const MOVE_CONVERGED: f64 = 1.8e-3;
/// Largest move of any atom in one step, bohr.
const TRUST: f64 = 0.3;

/// Relax a molecule's shape by following its forces: quasi-Newton (BFGS) on
/// the Cartesian coordinates, driven by gradients alone — the self-consistent
/// energy scatters by a few times 1e-7 hartree between nearby geometries, too
/// much to steer a line search by, while the gradient is smooth.
pub fn relax(start: &Molecule, f: Functional, max_steps: usize) -> (Molecule, Vec<Step>, bool) {
    let mut mol = start.clone();
    let n = 3 * mol.z.len();
    let flat = |g: &[[f64; 3]]| g.iter().flat_map(|v| v.iter().cloned()).collect::<Vec<f64>>();
    // Inverse Hessian, starting at a stiffness of about one hartree per bohr^2.
    let mut hinv = vec![0.0; n * n];
    for i in 0..n {
        hinv[i * n + i] = 1.0;
    }
    let (mut sol, g0) = mol.energy_and_gradient(f);
    let mut g = flat(&g0);
    let mut steps = Vec::new();
    let mut converged = false;
    for _ in 0..max_steps {
        let fmax = g.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        // Step: -H^-1 g, held inside the trust radius per atom.
        let mut dx: Vec<f64> = (0..n).map(|i| -(0..n).map(|j| hinv[i * n + j] * g[j]).sum::<f64>()).collect();
        let worst_atom = (0..n / 3).map(|a| (dx[3 * a].powi(2) + dx[3 * a + 1].powi(2) + dx[3 * a + 2].powi(2)).sqrt()).fold(0.0f64, f64::max);
        if worst_atom > TRUST {
            for v in &mut dx {
                *v *= TRUST / worst_atom;
            }
        }
        let dmax = dx.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        steps.push(Step { energy: sol.energy, largest_force: fmax, largest_move: dmax });
        if fmax < FORCE_CONVERGED && dmax < MOVE_CONVERGED {
            converged = true;
            break;
        }
        for (a, p) in mol.positions.iter_mut().enumerate() {
            for k in 0..3 {
                p[k] += dx[3 * a + k];
            }
        }
        let (s2, g2) = mol.energy_and_gradient(f);
        let g2 = flat(&g2);
        // BFGS update of the inverse Hessian, skipped if the curvature is not
        // positive along the step.
        let y: Vec<f64> = (0..n).map(|i| g2[i] - g[i]).collect();
        let sy: f64 = (0..n).map(|i| dx[i] * y[i]).sum();
        if sy > 1e-10 {
            let hy: Vec<f64> = (0..n).map(|i| (0..n).map(|j| hinv[i * n + j] * y[j]).sum()).collect();
            let yhy: f64 = (0..n).map(|i| y[i] * hy[i]).sum();
            for i in 0..n {
                for j in 0..n {
                    hinv[i * n + j] += (sy + yhy) * dx[i] * dx[j] / (sy * sy) - (hy[i] * dx[j] + dx[i] * hy[j]) / sy;
                }
            }
        }
        sol = s2;
        g = g2;
    }
    (mol, steps, converged)
}
