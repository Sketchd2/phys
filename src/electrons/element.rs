//! Each element's basis, derived from its own free atom.
//!
//! **Nothing here is a published basis set.** For each angular momentum the
//! engine lays down an even-tempered sequence `a, a b, a b^2, ...` and lets the
//! variational principle decide everything about it: the sequence is extended
//! at its diffuse end and at its tight end, and its spacing refined, for as
//! long as doing so lowers the free atom's energy by more than a tolerance.
//! The energy can only fall as the space grows, so the stopping point is where
//! the atom has stopped caring.
//!
//! **Which angular momenta are present is derived the same way.** A spherical
//! atom's orbitals each have one `l`, so an `l` with no electrons in it lowers
//! the energy by exactly nothing, and the search stops at the first `l` that
//! does not help. The shell structure of the periodic table is discovered here,
//! not stated.
//!
//! **And the ground state's spin**, by trying each multiplicity the electron
//! count allows and keeping the lowest — Hund's first rule as a consequence.

use super::atom::{response, response_functions, solve, solve_functions, Atom};
use super::functional::Functional;

/// A derived basis for one element: exponents per angular momentum.
#[derive(Debug, Clone, PartialEq)]
pub struct ElementBasis {
    pub z: u32,
    /// `(l, exponents ascending)`.
    pub shells: Vec<(usize, Vec<f64>)>,
    /// Unpaired electrons in the ground state.
    pub unpaired: u32,
    /// The free atom's energy in this basis, hartree.
    pub energy: f64,
    /// Self-consistent-field runs it took to derive.
    pub evaluations: usize,
    /// What a molecule is given: for each occupied angular momentum, the free
    /// atom's own orbitals as contracted functions, `(l, [(exponent, c)])`...
    pub contracted: Vec<(usize, Vec<(f64, f64)>)>,
    /// ...the most diffuse primitives of that `l` left free, `(l, exponents)`...
    pub free: Vec<(usize, Vec<f64>)>,
    /// ...and, above the occupied angular momenta, the polarisation sets'
    /// diffuse primitives left free beside their contracted response orbital
    /// (which is in `contracted`), `(l, exponents)`.
    pub polarisation: Vec<(usize, Vec<f64>)>,
    /// How far the contracted basis puts the atom's ions (+-1/2 e) from the
    /// uncontracted one, hartree.
    pub contraction_error: f64,
}

impl ElementBasis {
    /// Its shells for a molecule, placed at `centre` (bohr): contracted atomic
    /// orbitals, free diffuse primitives, polarisation sets.
    pub fn shells_at(&self, centre: [f64; 3]) -> Vec<super::basis::Shell> {
        use super::basis::Shell;
        let mut out = Vec::new();
        for (l, terms) in &self.contracted {
            out.push(Shell::contracted(centre, *l, terms.iter().map(|t| t.0).collect(), terms.iter().map(|t| t.1).collect()));
        }
        for (l, exps) in self.free.iter().chain(&self.polarisation) {
            for e in exps {
                out.push(Shell::primitive(centre, *l, *e));
            }
        }
        out
    }

    /// Auxiliary functions for fitting this element's densities: for each
    /// angular momentum `L` a product of two of its functions can carry, an
    /// even-tempered set spanning the exponents such products have
    /// (`a_i + a_j`), at ratio `AUXILIARY_RATIO`.
    ///
    /// **Capped at the orbital basis's own highest angular momentum**, not at
    /// twice it as products could carry. Measured on water (derived basis, f
    /// the highest): auxiliary L <= 6, 1745 functions, 1.4e-6 Ha from exact
    /// Coulomb in 155 s of integrals; L <= 4, 1402, 1.2e-6 in 63 s; L <= 3, 937,
    /// 1.8e-7 in 21 s. The high-L functions were nearly dependent and added
    /// round-off rather than fit.
    ///
    /// The cap is the highest angular momentum in the *whole molecule's*
    /// basis, which a single element does not know: capping hydrogen at its own
    /// d left water 5.4e-6 Ha below exact. So this takes it as an argument.
    pub fn auxiliary_at(&self, centre: [f64; 3], molecule_lmax: usize) -> Vec<super::basis::Shell> {
        self.auxiliary_up_to(centre, molecule_lmax)
    }

    /// The highest angular momentum among its functions.
    pub fn lmax(&self) -> usize {
        self.contracted.iter().map(|c| c.0).chain(self.free.iter().chain(&self.polarisation).map(|c| c.0)).max().unwrap_or(0)
    }

    /// As [`auxiliary_at`], with angular momenta above `cap` left out.
    pub fn auxiliary_up_to(&self, centre: [f64; 3], cap: usize) -> Vec<super::basis::Shell> {
        let mut prims: Vec<(usize, f64)> = Vec::new();
        for (l, terms) in &self.contracted {
            prims.extend(terms.iter().map(|t| (*l, t.0)));
        }
        for (l, exps) in self.free.iter().chain(&self.polarisation) {
            prims.extend(exps.iter().map(|e| (*l, *e)));
        }
        let lmax = prims.iter().map(|p| p.0).max().unwrap_or(0);
        let mut out = Vec::new();
        for big_l in 0..=(2 * lmax).min(cap) {
            let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
            for &(la, a) in &prims {
                for &(lb, b) in &prims {
                    if la + lb >= big_l && (la + lb - big_l) % 2 == 0 {
                        lo = lo.min(a + b);
                        hi = hi.max(a + b);
                    }
                }
            }
            if !(hi > 0.0) {
                continue;
            }
            let count = ((hi / lo).ln() / AUXILIARY_RATIO.ln()).ceil() as usize + 1;
            for k in 0..count {
                out.push(super::basis::Shell::primitive(centre, big_l, lo * AUXILIARY_RATIO.powi(k as i32)));
            }
        }
        out
    }

    /// Every primitive uncontracted: the basis the contraction is measured against.
    pub fn primitive_shells_at(&self, centre: [f64; 3]) -> Vec<super::basis::Shell> {
        let mut out = Vec::new();
        for (l, exps) in &self.shells {
            for e in exps {
                out.push(super::basis::Shell::primitive(centre, *l, *e));
            }
        }
        out
    }
}

/// The radial functions of a contracted basis, per `l`, for the atom solver.
fn radial_functions(contracted: &[(usize, Vec<(f64, f64)>)], free: &[(usize, Vec<f64>)], lmax: usize) -> Vec<Vec<Vec<(f64, f64)>>> {
    let mut out: Vec<Vec<Vec<(f64, f64)>>> = vec![Vec::new(); lmax + 1];
    for (l, terms) in contracted {
        out[*l].push(terms.clone());
    }
    for (l, exps) in free {
        for e in exps {
            out[*l].push(vec![(*e, 1.0)]);
        }
    }
    out
}

/// Spacing of the auxiliary (density-fitting) exponents.
pub const AUXILIARY_RATIO: f64 = 2.0;

/// How closely the contracted basis must reproduce the atom's ions, hartree.
pub const CONTRACTION_TOLERANCE: f64 = 1e-4;

/// How closely a contracted polarisation set must reproduce the ions'
/// response to a field, relatively.
pub const POLARISATION_TOLERANCE: f64 = 1e-3;

/// An even-tempered range: the most diffuse exponent, the ratio, the count.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Range {
    first: f64,
    ratio: f64,
    count: usize,
}

impl Range {
    fn exponents(&self) -> Vec<f64> {
        (0..self.count).map(|k| self.first * self.ratio.powi(k as i32)).collect()
    }
}

fn atom_energy(z: u32, ranges: &[(usize, Range)], unpaired: u32, f: Functional) -> Atom {
    // Exponents per l, l = 0.. in order; an l absent from `ranges` gets none.
    let lmax = ranges.iter().map(|(l, _)| *l).max().unwrap_or(0);
    let mut exps: Vec<Vec<f64>> = vec![Vec::new(); lmax + 1];
    for (l, r) in ranges {
        exps[*l] = r.exponents();
    }
    let n = z as f64;
    solve(n, (n + unpaired as f64) / 2.0, (n - unpaired as f64) / 2.0, &exps, f, 500, 1e-10)
}

/// Derive the basis of element `z` to `tolerance` hartree per refinement.
pub fn derive(z: u32, f: Functional, tolerance: f64) -> ElementBasis {
    let mut evaluations = 0;
    let zf = z as f64;
    // A starting range for any l: from a diffuse exponent that suits a valence
    // electron to a tight one on the scale of the 1s orbital, Z^2.
    let start = |l: usize| Range { first: 0.05 / (1 + l) as f64, ratio: 3.0, count: ((10.0 * zf * zf).ln() / 3f64.ln()).ceil() as usize + 3 };
    let mut eval = |ranges: &[(usize, Range)], unpaired: u32| {
        evaluations += 1;
        let s = atom_energy(z, ranges, unpaired, f);
        if s.converged { s.energy } else { f64::INFINITY }
    };
    // Spin: the lowest of the multiplicities the count allows, on a starting
    // basis with s and p (a p shell is needed for an open p shell to show).
    let parity = z % 2;
    let probe: Vec<(usize, Range)> = (0..=1).map(|l| (l, start(l))).collect();
    let mut unpaired = parity;
    let mut best = f64::INFINITY;
    let mut u = parity;
    while u <= z.min(8) {
        let e = eval(&probe, u);
        if e < best - 1e-7 {
            best = e;
            unpaired = u;
        } else if e > best {
            break;
        }
        u += 2;
    }
    // Angular momenta: add l while it helps.
    let mut ranges: Vec<(usize, Range)> = vec![(0, start(0))];
    let mut energy = eval(&ranges, unpaired);
    for l in 1..=4 {
        let mut trial = ranges.clone();
        trial.push((l, start(l)));
        let e = eval(&trial, unpaired);
        if e < energy - tolerance {
            ranges = trial;
            energy = e;
        } else {
            break;
        }
    }
    // Refine each l: extend the ends and tighten the ratio while it helps.
    loop {
        let mut improved = false;
        for k in 0..ranges.len() {
            // Diffuse end.
            loop {
                let mut t = ranges.clone();
                t[k].1.first /= t[k].1.ratio;
                t[k].1.count += 1;
                let e = eval(&t, unpaired);
                if e < energy - tolerance {
                    ranges = t;
                    energy = e;
                    improved = true;
                } else {
                    break;
                }
            }
            // Tight end.
            loop {
                let mut t = ranges.clone();
                t[k].1.count += 1;
                let e = eval(&t, unpaired);
                if e < energy - tolerance {
                    ranges = t;
                    energy = e;
                    improved = true;
                } else {
                    break;
                }
            }
            // Spacing: the same span with a finer ratio.
            let r = ranges[k].1;
            if r.ratio > 1.6 {
                let span = r.ratio.powi(r.count as i32 - 1);
                let ratio = r.ratio.powf(0.85);
                let count = (span.ln() / ratio.ln()).ceil() as usize + 1;
                let mut t = ranges.clone();
                t[k].1 = Range { first: r.first, ratio, count };
                let e = eval(&t, unpaired);
                if e < energy - tolerance {
                    ranges = t;
                    energy = e;
                    improved = true;
                }
            }
        }
        if !improved {
            break;
        }
    }
    // Polarisation. The free atom has no use for functions above its highest
    // occupied angular momentum — they lower its energy by nothing — but a
    // molecule's field is what they are for, so they are derived from the
    // atom's response to one: one step up from a uniform field (dipole), two
    // from a field gradient (quadrupole). Each set is extended while that
    // response still grows by more than a part in a thousand.
    let atom = atom_energy(z, &ranges, unpaired, f);
    evaluations += 1;
    let top = ranges.iter().map(|(l, _)| *l).max().unwrap_or(0);
    for order in 1..=2 {
        let channel = top + order;
        let mut pol = Range { first: 0.1, ratio: 3.0, count: 3 };
        let relative = 1e-3;
        let mut e2 = response(&atom, channel, order, &pol.exponents());
        loop {
            let mut improved = false;
            let mut trials = vec![
                Range { first: pol.first / pol.ratio, ratio: pol.ratio, count: pol.count + 1 },
                Range { first: pol.first, ratio: pol.ratio, count: pol.count + 1 },
            ];
            if pol.ratio > 1.6 {
                let span = pol.ratio.powi(pol.count as i32 - 1);
                let ratio = pol.ratio.powf(0.85);
                trials.push(Range { first: pol.first, ratio, count: (span.ln() / ratio.ln()).ceil() as usize + 1 });
            }
            for t in trials {
                let e = response(&atom, channel, order, &t.exponents());
                if e < e2 - relative * e2.abs() {
                    pol = t;
                    e2 = e;
                    improved = true;
                    break;
                }
            }
            if !improved {
                break;
            }
        }
        ranges.push((channel, pol));
    }
    // Contraction. Each occupied level becomes one function: its alpha and
    // beta orbitals averaged, weighted by their electrons, signs aligned. Then
    // diffuse primitives are freed, one angular momentum at a time, until the
    // atom's ions — what a molecule does to an atom is mostly move charge on
    // or off it — come out of the contracted basis as they do out of the full
    // one, to `CONTRACTION_TOLERANCE`.
    let occupied: Vec<(usize, Range)> = ranges.iter().filter(|(l, _)| *l <= top).cloned().collect();
    let mut contracted: Vec<(usize, Vec<(f64, f64)>)> = Vec::new();
    for (l, r) in &occupied {
        let exps = r.exponents();
        let alpha: Vec<&(usize, usize, f64, Vec<f64>)> = atom.coefficients.iter().filter(|c| c.0 == *l && c.1 == 0).collect();
        let beta: Vec<&(usize, usize, f64, Vec<f64>)> = atom.coefficients.iter().filter(|c| c.0 == *l && c.1 == 1).collect();
        for (i, a) in alpha.iter().enumerate() {
            let mut c = a.3.iter().map(|x| x * a.2).collect::<Vec<f64>>();
            let mut w = a.2;
            if let Some(b) = beta.get(i) {
                let sign = if a.3.iter().zip(&b.3).map(|(x, y)| x * y).sum::<f64>() < 0.0 { -1.0 } else { 1.0 };
                for (ci, bi) in c.iter_mut().zip(&b.3) {
                    *ci += sign * b.2 * bi;
                }
                w += b.2;
            }
            for ci in &mut c {
                *ci /= w;
            }
            contracted.push((*l, exps.iter().cloned().zip(c).collect()));
        }
    }
    let ion = |q: f64| -> (f64, f64) {
        let n = zf - q;
        let (a0, b0) = ((zf + unpaired as f64) / 2.0, (zf - unpaired as f64) / 2.0);
        let mut b = b0 - q / 2.0;
        let mut a = a0 - q / 2.0;
        if b < 0.0 {
            a += b;
            b = 0.0;
        }
        let _ = n;
        (a, b)
    };
    let full: Vec<Vec<Vec<(f64, f64)>>> = radial_functions(&[], &occupied.iter().map(|(l, r)| (*l, r.exponents())).collect::<Vec<_>>(), top);
    let ions: Vec<Atom> = [0.5, -0.5].iter().map(|&q| {
        let (a, b) = ion(q);
        solve_functions(zf, a, b, &full, f, 500, 1e-10)
    }).collect();
    let reference: Vec<f64> = ions.iter().map(|a| a.energy).collect();
    evaluations += 2;
    let mut free_count: Vec<usize> = vec![1; occupied.len()];
    let free_of = |counts: &[usize]| -> Vec<(usize, Vec<f64>)> {
        occupied.iter().zip(counts).map(|((l, r), &k)| (*l, r.exponents()[..k.min(r.count)].to_vec())).collect()
    };
    let error_of = |counts: &[usize], evaluations: &mut usize| -> f64 {
        let funcs = radial_functions(&contracted, &free_of(counts), top);
        [0.5, -0.5].iter().zip(&reference).map(|(&q, r)| {
            *evaluations += 1;
            let (a, b) = ion(q);
            (solve_functions(zf, a, b, &funcs, f, 500, 1e-10).energy - r).abs()
        }).fold(0.0, f64::max)
    };
    let mut err = error_of(&free_count, &mut evaluations);
    while err > CONTRACTION_TOLERANCE {
        let mut best: Option<(usize, f64)> = None;
        for k in 0..free_count.len() {
            if free_count[k] >= occupied[k].1.count {
                continue;
            }
            let mut t = free_count.clone();
            t[k] += 1;
            let e = error_of(&t, &mut evaluations);
            if best.map_or(true, |(_, b)| e < b) {
                best = Some((k, e));
            }
        }
        match best {
            Some((k, e)) => {
                free_count[k] += 1;
                err = e;
            }
            None => break,
        }
    }
    let free = free_of(&free_count);
    // The polarisation sets contract the same way, into the response orbital
    // itself: the first-order change of the occupied orbitals in a field, which
    // is the optimal function of its angular momentum by construction. On the
    // neutral atom that reproduces the full set's response exactly, so whether
    // it is *flexible* enough is judged on the ions again: diffuse primitives
    // are freed until each ion's response agrees with the full set's to
    // `POLARISATION_TOLERANCE`, relatively.
    let mut polarisation: Vec<(usize, Vec<f64>)> = Vec::new();
    for (order, (channel, r)) in ranges.iter().filter(|(l, _)| *l > top).enumerate() {
        let order = order + 1;
        let exps = r.exponents();
        let prims: Vec<Vec<(f64, f64)>> = exps.iter().map(|&a| vec![(a, 1.0)]).collect();
        let (_, orbital) = response_functions(&atom, *channel, order, &prims);
        let shape: Vec<(f64, f64)> = exps.iter().cloned().zip(orbital).collect();
        let full: Vec<f64> = ions.iter().map(|ion| response_functions(ion, *channel, order, &prims).0).collect();
        let mut k = 0;
        loop {
            let mut funcs = vec![shape.clone()];
            funcs.extend(exps[..k].iter().map(|&a| vec![(a, 1.0)]));
            let worst = ions.iter().zip(&full).map(|(ion, e)| {
                let c = response_functions(ion, *channel, order, &funcs).0;
                ((c - e) / e).abs()
            }).fold(0.0, f64::max);
            if worst <= POLARISATION_TOLERANCE || k >= exps.len() {
                break;
            }
            k += 1;
        }
        // If every primitive had to be freed the contracted function is an
        // exact combination of them and adds nothing but a dependency. (This
        // was first written as the cause of water's 87 SCF iterations; it was
        // not — the same run without it took 87 too.)
        if k < exps.len() {
            contracted.push((*channel, shape));
        }
        polarisation.push((*channel, exps[..k].to_vec()));
    }
    ElementBasis {
        z,
        shells: ranges.iter().map(|(l, r)| (*l, r.exponents())).collect(),
        unpaired,
        energy,
        evaluations,
        contracted,
        free,
        polarisation,
        contraction_error: err,
    }
}
