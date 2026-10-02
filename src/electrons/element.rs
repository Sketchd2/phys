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

use super::atom::{response, solve, Atom};
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
}

impl ElementBasis {
    /// Its shells, placed at `centre` (bohr).
    pub fn shells_at(&self, centre: [f64; 3]) -> Vec<super::basis::Shell> {
        let mut out = Vec::new();
        for (l, exps) in &self.shells {
            for e in exps {
                out.push(super::basis::Shell::primitive(centre, *l, *e));
            }
        }
        out
    }
}

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
    ElementBasis {
        z,
        shells: ranges.iter().map(|(l, r)| (*l, r.exponents())).collect(),
        unpaired,
        energy,
        evaluations,
    }
}
