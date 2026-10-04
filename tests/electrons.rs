//! The electronic-structure numerics (`docs/PLAY.md` Phase 6, stage E1),
//! against results that are exact.

use phys::electrons::basis::{Basis, Shell};
use phys::electrons::boys::boys;
use phys::electrons::integrals::{eri_block, one_electron};
use phys::electrons::linalg::{eigh, generalised, orthogonaliser, Matrix};

const PI: f64 = std::f64::consts::PI;

/// `F_n(T)` by direct quadrature, Simpson's rule on a fine grid.
fn boys_by_quadrature(n: usize, t: f64) -> f64 {
    let m = 200_000;
    let h = 1.0 / m as f64;
    let f = |x: f64| x.powi(2 * n as i32) * (-t * x * x).exp();
    let mut s = f(0.0) + f(1.0);
    for k in 1..m {
        s += f(k as f64 * h) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    s * h / 3.0
}

#[test]
fn the_boys_function_is_the_integral_it_names() {
    let mut worst: f64 = 0.0;
    let mut out = [0.0; 13];
    for &t in &[0.0, 1e-9, 0.3, 1.0, 7.5, 24.0, 49.9, 50.1, 120.0] {
        boys(12, t, &mut out);
        for n in 0..=12 {
            let q = boys_by_quadrature(n, t);
            let e = (out[n] - q).abs() / q.max(1e-300);
            worst = worst.max(e);
        }
    }
    println!("  F_0..F_12 over T in [0, 120]: worst relative error {worst:.2e}");
    assert!(worst < 1e-10, "{worst}");
}

/// One normalised Gaussian on a proton: `E(a) = 3a/2 - 2 sqrt(2a/pi)`, least at
/// `a = 8/(9 pi)` where it is `-4/(3 pi)`.
#[test]
fn one_gaussian_on_a_proton() {
    let a = 8.0 / (9.0 * PI);
    let basis = Basis::new(vec![Shell::primitive([0.0; 3], 0, a)]);
    let (s, t, v) = one_electron(&basis, &[(1.0, [0.0; 3])]);
    let e = t.get(0, 0) + v.get(0, 0);
    println!("  S = {:.15}, E = {e:.15} against {:.15}", s.get(0, 0), -4.0 / (3.0 * PI));
    assert!((s.get(0, 0) - 1.0).abs() < 1e-14);
    assert!((t.get(0, 0) - 1.5 * a).abs() < 1e-14);
    assert!((e + 4.0 / (3.0 * PI)).abs() < 1e-14);
}

fn even_tempered(centre: [f64; 3], l: usize, first: f64, ratio: f64, count: usize) -> Vec<Shell> {
    (0..count).map(|k| Shell::primitive(centre, l, first * ratio.powi(k as i32))).collect()
}

fn lowest(basis: &Basis, nuclei: &[(f64, [f64; 3])]) -> (f64, usize) {
    let (s, t, v) = one_electron(basis, nuclei);
    let mut h = t.clone();
    for k in 0..h.a.len() {
        h.a[k] += v.a[k];
    }
    let (x, m) = orthogonaliser(&s, 1e-9);
    let (e, _) = generalised(&h, &x, m);
    (e[0], m)
}

/// The hydrogen atom, whose ground state is exactly -1/2 hartree, in an
/// even-tempered set of s Gaussians.
#[test]
fn the_hydrogen_atom_converges_to_a_half() {
    for count in [8, 14, 20, 26] {
        let basis = Basis::new(even_tempered([0.0; 3], 0, 0.01, 2.5, count));
        let (e, m) = lowest(&basis, &[(1.0, [0.0; 3])]);
        println!("  {count:2} s functions ({m} kept): E = {e:.10}, error {:.2e}", e + 0.5);
        if count == 26 {
            assert!((e + 0.5).abs() < 1e-6, "{e}");
        }
    }
}

/// Two Gaussian charge clouds of unit charge repel exactly as
/// `erf(sqrt(g) R) / R`, `g = 2a 2c / (2a + 2c)`.
#[test]
fn two_gaussian_charges_repel_by_coulombs_law() {
    let (a, c) = (0.8, 1.7);
    for r in [0.0, 0.4, 1.3, 6.0] {
        let sa = Shell::primitive([0.0; 3], 0, a);
        let sc = Shell::primitive([0.0, 0.0, r], 0, c);
        let j = eri_block(&sa, &sa, &sc, &sc)[0];
        let g: f64 = (2.0 * a) * (2.0 * c) / (2.0 * a + 2.0 * c);
        let exact = if r == 0.0 { 2.0 * (g / PI).sqrt() } else { erf(g.sqrt() * r) / r };
        println!("  R = {r}: (aa|cc) = {j:.14} against {exact:.14}");
        assert!((j - exact).abs() < 1e-12);
    }
}

/// erf by its Taylor series (small argument) and continued-fraction tail.
fn erf(x: f64) -> f64 {
    if x < 3.0 {
        let mut term = x;
        let mut sum = x;
        let mut n = 0;
        loop {
            n += 1;
            term *= -x * x / n as f64;
            let add = term / (2 * n + 1) as f64;
            sum += add;
            if add.abs() < 1e-18 {
                break;
            }
        }
        2.0 / PI.sqrt() * sum
    } else {
        // erfc by Lentz continued fraction.
        let mut f = x;
        let mut c = x;
        let mut d = 0.0;
        for k in 1..300 {
            let an = k as f64 / 2.0;
            d = x + an * d;
            d = 1.0 / d;
            c = x + an / c;
            let delta = c * d;
            f *= delta;
            if (delta - 1.0).abs() < 1e-16 {
                break;
            }
        }
        1.0 - (-x * x).exp() / PI.sqrt() / f
    }
}

/// Permutational symmetry of `(ab|cd)` with `p` and `d` shells, which is where a
/// slip in the Hermite recursion would show.
#[test]
fn electron_repulsion_has_its_eightfold_symmetry() {
    let a = Shell::primitive([0.1, -0.2, 0.3], 1, 0.9);
    let b = Shell::primitive([1.0, 0.4, -0.5], 2, 0.6);
    let c = Shell::primitive([-0.7, 0.2, 0.8], 1, 1.3);
    let d = Shell::primitive([0.3, -1.1, 0.0], 0, 0.45);
    let abcd = eri_block(&a, &b, &c, &d);
    let bacd = eri_block(&b, &a, &c, &d);
    let cdab = eri_block(&c, &d, &a, &b);
    let (na, nb, nc, nd) = (a.size(), b.size(), c.size(), d.size());
    let mut worst: f64 = 0.0;
    for i in 0..na {
        for j in 0..nb {
            for k in 0..nc {
                for l in 0..nd {
                    let x = abcd[((i * nb + j) * nc + k) * nd + l];
                    let y = bacd[((j * na + i) * nc + k) * nd + l];
                    let z = cdab[((k * nd + l) * na + i) * nb + j];
                    worst = worst.max((x - y).abs()).max((x - z).abs());
                }
            }
        }
    }
    println!("  (pd|ps) against (dp|ps) and (ps|pd): worst {worst:.2e}");
    assert!(worst < 1e-13);
}

fn rotate(p: [f64; 3], r: &[[f64; 3]; 3]) -> [f64; 3] {
    [
        r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2],
        r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2],
        r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2],
    ]
}

/// A molecule turned in space has the same spectrum. Cartesian shells span
/// rotations of themselves, so every level of the one-electron Hamiltonian in
/// the overlap's metric must be identical before and after.
///
/// Not the overlap's own eigenvalues: the change of basis a rotation makes
/// among Cartesian `d` functions is not orthogonal (`xx` and `xy` normalise
/// differently), so those move — by 0.31 here — while the physics does not. The
/// first version of this test asserted that and failed for exactly that reason.
#[test]
fn turning_the_molecule_changes_nothing() {
    let atoms = [(8.0, [0.0, 0.0, 0.0]), (1.0, [0.0, 1.43, 1.1]), (1.0, [0.0, -1.43, 1.1])];
    let build = |rot: &[[f64; 3]; 3]| {
        let mut shells = Vec::new();
        let mut nuclei = Vec::new();
        for (z, p) in atoms {
            let q = rotate(p, rot);
            nuclei.push((z, q));
            for l in 0..=2 {
                for e in [0.3, 1.1, 4.0] {
                    shells.push(Shell::primitive(q, l, e));
                }
            }
        }
        let basis = Basis::new(shells);
        let (s, t, v) = one_electron(&basis, &nuclei);
        let mut h = t.clone();
        for k in 0..h.a.len() {
            h.a[k] += v.a[k];
        }
        let (x, m) = orthogonaliser(&s, 1e-9);
        (generalised(&h, &x, m).0, h)
    };
    let (c, s) = (0.6f64.cos(), 0.6f64.sin());
    let (c2, s2) = (1.1f64.cos(), 1.1f64.sin());
    let rx = [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]];
    let rz = [[c2, -s2, 0.0], [s2, c2, 0.0], [0.0, 0.0, 1.0]];
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += rz[i][k] * rx[k][j];
            }
        }
    }
    let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let (e0, h0) = build(&id);
    let (e1, _) = build(&r);
    assert_eq!(e0.len(), e1.len());
    let worst = e0.iter().zip(&e1).map(|(a, b)| (a - b).abs() / a.abs().max(1.0)).fold(0.0, f64::max);
    println!("  {} one-electron levels, lowest {:.12}; turned, the worst moved {worst:.2e}", e0.len(), e0[0]);
    assert!(h0.asymmetry() < 1e-12, "the Hamiltonian is symmetric");
    assert!(worst < 1e-9);
}

/// H2+ at 2 bohr: one electron, two protons, and an exact electronic energy of
/// -1.1026342144949 hartree (Bates, Ledsham and Stewart; Peek).
#[test]
fn the_hydrogen_molecular_ion() {
    let centres = [[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]];
    let nuclei: Vec<(f64, [f64; 3])> = centres.iter().map(|c| (1.0, *c)).collect();
    let mut last = 0.0;
    for (ns, np, nd) in [(10, 0, 0), (14, 4, 0), (18, 6, 4)] {
        let mut shells = Vec::new();
        for c in centres {
            shells.extend(even_tempered(c, 0, 0.02, 2.6, ns));
            shells.extend(even_tempered(c, 1, 0.1, 2.6, np));
            shells.extend(even_tempered(c, 2, 0.2, 2.6, nd));
        }
        let basis = Basis::new(shells);
        let (e, m) = lowest(&basis, &nuclei);
        println!("  s{ns} p{np} d{nd} per centre ({} functions, {m} kept): E = {e:.10}, error {:.2e}", basis.size, e + 1.1026342144949);
        last = e;
    }
    assert!((last + 1.1026342144949).abs() < 2e-4, "{last}");
    assert!(last > -1.1026342144949 - 1e-9, "variational: a basis can only lie above the exact answer");
}

#[test]
fn the_eigensolver_solves() {
    // A random symmetric matrix: A v = l v for every pair, and V orthonormal.
    let n = 40;
    let mut m = Matrix::zeros(n);
    let mut x: u64 = 12345;
    for i in 0..n {
        for j in 0..=i {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let v = ((x >> 11) as f64 / (1u64 << 53) as f64) - 0.5;
            m.set(i, j, v);
            m.set(j, i, v);
        }
    }
    let (vals, vecs) = eigh(&m);
    let mut worst: f64 = 0.0;
    for k in 0..n {
        for i in 0..n {
            let mut av = 0.0;
            for j in 0..n {
                av += m.get(i, j) * vecs.get(j, k);
            }
            worst = worst.max((av - vals[k] * vecs.get(i, k)).abs());
        }
    }
    let vtv = vecs.transpose().mul(&vecs);
    let ortho = (0..n * n).map(|k| (vtv.a[k] - if k % (n + 1) == 0 { 1.0 } else { 0.0 }).abs()).fold(0.0, f64::max);
    println!("  40x40: residual {worst:.2e}, orthonormality {ortho:.2e}");
    assert!(worst < 1e-12 && ortho < 1e-12);
    assert!(vals.windows(2).all(|w| w[0] <= w[1]));
}

/// Stage E2: the grid reproduces the analytic overlap of a water molecule's
/// basis, s to d, every pair — a stricter test than integrating an electron
/// count, which is one sum of these. With Becke's size adjustment the same
/// grid reached only 2.5e-5; see `grid::molecular`.
#[test]
fn the_grid_integrates_a_molecules_functions() {
    use phys::electrons::grid::molecular;
    use phys::electrons::values::at;
    let atoms = [(8.0, [0.0, 0.0, 0.0], 1.2), (1.0, [0.0, 1.43, 1.1], 0.6), (1.0, [0.0, -1.43, 1.1], 0.6)];
    let mut shells = Vec::new();
    for (_, p, _) in atoms {
        for l in 0..=2 {
            for e in [0.15, 0.6, 2.5, 10.0] {
                shells.push(Shell::primitive(p, l, e));
            }
        }
    }
    let basis = Basis::new(shells);
    let (s, _, _) = one_electron(&basis, &atoms.map(|(z, p, _)| (z, p)));
    let n = basis.size;
    for (nr, nt) in [(40, 12), (60, 18), (90, 24)] {
        let grid = molecular(&atoms.map(|(_, p, r)| (p, r)), nr, nt);
        let mut g = vec![0.0; n * n];
        let mut v = vec![0.0; n];
        for (p, w) in grid.points.iter().zip(&grid.weights) {
            at(&basis, *p, &mut v, None);
            for i in 0..n {
                let wi = w * v[i];
                for j in 0..n {
                    g[i * n + j] += wi * v[j];
                }
            }
        }
        let worst = (0..n * n).map(|k| (g[k] - s.a[k]).abs()).fold(0.0, f64::max);
        println!("  {nr} radial x {} angular per atom, {} points: worst overlap error {worst:.2e}", 2 * nt * nt, grid.points.len());
        if nr == 90 {
            assert!(worst < 1e-6, "{worst}");
        }
    }
}

/// Gradients of the basis functions agree with finite differences.
#[test]
fn basis_gradients_are_derivatives() {
    use phys::electrons::values::at;
    let basis = Basis::new(vec![
        Shell::primitive([0.1, 0.2, -0.3], 0, 0.7),
        Shell::primitive([0.1, 0.2, -0.3], 1, 1.3),
        Shell::primitive([-0.4, 0.0, 0.5], 2, 0.9),
        Shell::primitive([0.0, 0.3, 0.0], 3, 0.5),
    ]);
    let n = basis.size;
    let p = [0.35, -0.2, 0.6];
    let mut v = vec![0.0; n];
    let mut g = vec![0.0; 3 * n];
    at(&basis, p, &mut v, Some(&mut g));
    let h = 1e-6;
    let mut worst: f64 = 0.0;
    for d in 0..3 {
        let mut pp = p;
        let mut pm = p;
        pp[d] += h;
        pm[d] -= h;
        let mut vp = vec![0.0; n];
        let mut vm = vec![0.0; n];
        at(&basis, pp, &mut vp, None);
        at(&basis, pm, &mut vm, None);
        for i in 0..n {
            worst = worst.max(((vp[i] - vm[i]) / (2.0 * h) - g[3 * i + d]).abs());
        }
    }
    println!("  s, p, d, f: worst gradient error {worst:.2e}");
    assert!(worst < 1e-8);
}

// ---------------------------------------------------------------------------
// Stage E3: Kohn-Sham, against an independent implementation
// ---------------------------------------------------------------------------
//
// The reference numbers are PySCF 2.14 with libxc ("slater,pw_mod" and
// "pbe,pbe"), unrestricted, Cartesian, on *exactly these basis sets*, grid
// level 9 and conv_tol 1e-11. Same basis, same functional: the only things
// that can differ are the grid and the code, so agreement says the integrals,
// the functional, its potential and the self-consistency are all right. How
// far either is from the basis-set limit is stage E4's question, not this one.

use phys::electrons::functional::Functional;
use phys::electrons::scf::{solve, Problem};

fn atom_basis() -> Basis {
    let mut shells = Vec::new();
    shells.extend(even_tempered([0.0; 3], 0, 0.05, 3.0, 12));
    shells.extend(even_tempered([0.0; 3], 1, 0.1, 3.0, 6));
    shells.extend(even_tempered([0.0; 3], 2, 0.3, 3.0, 2));
    Basis::new(shells)
}

#[test]
fn atoms_agree_with_an_independent_implementation() {
    let cases = [
        ("H", 1.0, 1.0, 0.0, Functional::Lda, -0.47870014),
        ("H", 1.0, 1.0, 0.0, Functional::Pbe, -0.49997908),
        ("He", 2.0, 1.0, 1.0, Functional::Lda, -2.83443563),
        ("He", 2.0, 1.0, 1.0, Functional::Pbe, -2.89291612),
        ("Be", 4.0, 2.0, 2.0, Functional::Lda, -14.44617263),
        ("Be", 4.0, 2.0, 2.0, Functional::Pbe, -14.62961013),
        ("N", 7.0, 5.0, 2.0, Functional::Lda, -54.13047335),
        ("N", 7.0, 5.0, 2.0, Functional::Pbe, -54.53188572),
        ("Ne", 10.0, 5.0, 5.0, Functional::Lda, -128.17827154),
        ("Ne", 10.0, 5.0, 5.0, Functional::Pbe, -128.81428838),
    ];
    let mut worst: f64 = 0.0;
    for (name, z, a, b, f, reference) in cases {
        let p = Problem { basis: atom_basis(), nuclei: vec![(z, [0.0; 3])], sizes: vec![1.0], alpha: a, beta: b, functional: f, radial: 100, theta: 8, auxiliary: None, prune: false, guess: None, nonlocal: None };
        let s = solve(&p, 100, 1e-10);
        assert!(s.converged, "{name} {f:?} did not converge");
        let e = (s.energy - reference).abs();
        worst = worst.max(e);
        println!("  {name:2} {f:?}: {:.8} against {reference:.8} ({e:.1e}), {} iterations", s.energy, s.iterations);
    }
    assert!(worst < 2e-8, "{worst}");
}

/// A molecule: the Becke partition, two-centre integrals and the nuclei's own
/// repulsion all enter. Reference -75.81621481 (PySCF, same basis).
#[test]
fn water_agrees_with_an_independent_implementation() {
    let o = [0.0, 0.0, 0.0];
    let h1 = [0.0, 1.43, 1.1];
    let h2 = [0.0, -1.43, 1.1];
    let mut sh = even_tempered(o, 0, 0.1, 3.0, 8);
    sh.extend(even_tempered(o, 1, 0.15, 3.0, 4));
    sh.push(Shell::primitive(o, 2, 0.8));
    for c in [h1, h2] {
        sh.extend(even_tempered(c, 0, 0.06, 3.0, 5));
        sh.push(Shell::primitive(c, 1, 0.7));
    }
    let p = Problem {
        basis: Basis::new(sh),
        nuclei: vec![(8.0, o), (1.0, h1), (1.0, h2)],
        sizes: vec![1.0, 0.6, 0.6],
        alpha: 5.0,
        beta: 5.0,
        functional: Functional::Pbe,
        radial: 70,
        theta: 16,
        auxiliary: None,
        prune: false,
        guess: None,
        nonlocal: None,
    };
    let s = solve(&p, 100, 1e-10);
    println!("  water, PBE: {:.8} against -75.81621481 ({:.1e}), {} iterations", s.energy, (s.energy + 75.81621481).abs(), s.iterations);
    assert!(s.converged);
    assert!((s.energy + 75.81621481).abs() < 1e-6);
}

#[test]
fn degenerate_levels_share_their_electrons() {
    use phys::electrons::scf::occupy;
    // 1s, 2s, three 2p: carbon's alpha electrons, four of them.
    let occ = occupy(&[-10.0, -0.5, -0.2, -0.2, -0.2, 0.3], 4.0);
    println!("  {occ:?}");
    assert_eq!(occ[0], 1.0);
    assert_eq!(occ[1], 1.0);
    for o in &occ[2..5] {
        assert!((o - 2.0 / 3.0).abs() < 1e-15);
    }
    assert_eq!(occ[5], 0.0);
}

// ---------------------------------------------------------------------------
// Stage E4: the spherical atom, and each element's basis derived from it
// ---------------------------------------------------------------------------

/// The radial solver is the three-dimensional one restricted to a spherical
/// atom, so on the same s and p functions they must agree. Measured: H to
/// 4e-11, Ne to 6e-8 at 600 radial points (9e-7 at 300, which is why it is 600).
#[test]
fn the_radial_atom_is_the_three_dimensional_one() {
    use phys::electrons::atom;
    let s: Vec<f64> = (0..12).map(|k| 0.05 * 3f64.powi(k)).collect();
    let p: Vec<f64> = (0..6).map(|k| 0.1 * 3f64.powi(k)).collect();
    let mut worst: f64 = 0.0;
    for (name, z, a, b) in [("H", 1.0, 1.0, 0.0), ("C", 6.0, 4.0, 2.0), ("Ne", 10.0, 5.0, 5.0)] {
        let mut sh: Vec<Shell> = s.iter().map(|e| Shell::primitive([0.0; 3], 0, *e)).collect();
        sh.extend(p.iter().map(|e| Shell::primitive([0.0; 3], 1, *e)));
        let d3 = solve(&Problem { basis: Basis::new(sh), nuclei: vec![(z, [0.0; 3])], sizes: vec![1.0], alpha: a, beta: b, functional: Functional::Pbe, radial: 100, theta: 8, auxiliary: None, prune: false, guess: None, nonlocal: None }, 200, 1e-11);
        let r = atom::solve(z, a, b, &[s.clone(), p.clone()], Functional::Pbe, 500, 1e-11);
        let d = (r.energy - d3.energy).abs();
        worst = worst.max(d);
        println!("  {name:2}: three-dimensional {:.9}, radial {:.9}, {d:.1e}", d3.energy, r.energy);
    }
    assert!(worst < 2e-7, "{worst}");
}

/// Stage E4's done-when for atoms, and E3's clause against the basis limit:
/// each element's own derived basis gives its free-atom PBE energy to 1e-4 Ha
/// of the limit. The limits are PySCF with even-tempered s/p sets of 32 and 26
/// functions at ratio 1.9 — large enough that the next size changes them by
/// 1e-7 (H) to 1.2e-4 (Ar). Measured over H to Ar: worst 1.7e-5 (B).
///
/// Three things the derivation decides by itself are checked as well: the
/// ground state's spin (Hund's first rule, not told), and which angular momenta
/// an element needs (only s up to beryllium, s and p after: the periodic
/// table's shell structure, not told).
#[test]
fn each_element_derives_its_own_basis() {
    use phys::electrons::element::derive;
    // (Z, basis-limit PBE energy, unpaired, angular momenta)
    let cases: [(u32, f64, u32, usize); 4] = [(1, -0.4999904, 1, 1), (4, -14.6299392, 0, 1), (6, -37.7936834, 2, 2), (10, -128.8664256, 0, 2)];
    for (z, limit, unpaired, ls) in cases {
        let b = derive(z, Functional::Pbe, 1e-5);
        let counts: Vec<String> = b.shells.iter().map(|(l, e)| format!("l={l}: {}", e.len())).collect();
        println!("  Z = {z:2}: {:.7} against a limit of {limit:.7} ({:+.1e}); {} unpaired; {}; {} atoms solved", b.energy, b.energy - limit, b.unpaired, counts.join(", "), b.evaluations);
        assert!((b.energy - limit).abs() < 1e-4);
        assert!(b.energy > limit - 2e-5, "a finite basis cannot go far below the limit");
        assert_eq!(b.unpaired, unpaired, "the ground state's spin");
        // The occupied angular momenta, then two polarisation sets derived from
        // the atom's response to a field and to a field gradient.
        assert_eq!(b.shells.len(), ls + 2, "the angular momenta it needs");
        assert_eq!(b.shells.iter().map(|(l, _)| *l).collect::<Vec<_>>(), (0..ls + 2).collect::<Vec<_>>());
    }
}

/// Stage E4's done-when for molecules: a molecule in the engine's derived,
/// contracted basis, with density fitting and the pruned grid, is no worse
/// than a published quadruple-zeta set. H2 at 1.4 bohr, PBE, against PySCF:
/// aug-cc-pVQZ -1.16654088, aug-cc-pV5Z -1.16667410. Measured: -1.16654368
/// with 108 functions (5Z has 160). Water told the same story: -76.38748642
/// with 214 functions, below aug-cc-pVQZ (-76.38610) and 0.42 mHa above
/// aug-cc-pV5Z (-76.38791, 287 functions); it takes about a minute, so it is
/// recorded in `docs/PLAY.md` rather than run here.
#[test]
fn a_molecule_in_the_derived_basis_is_quadruple_zeta() {
    use phys::electrons::element::derive;
    let h = derive(1, Functional::Pbe, 1e-5);
    let (a, b) = ([0.0, 0.0, -0.7], [0.0, 0.0, 0.7]);
    let lmax = h.lmax();
    let mut sh = h.shells_at(a);
    sh.extend(h.shells_at(b));
    let mut ax = h.auxiliary_at(a, lmax);
    ax.extend(h.auxiliary_at(b, lmax));
    let p = Problem { basis: Basis::new(sh), nuclei: vec![(1.0, a), (1.0, b)], sizes: vec![0.6, 0.6], alpha: 1.0, beta: 1.0, functional: Functional::Pbe, radial: 70, theta: 16, auxiliary: Some(Basis::new(ax)), prune: true, guess: None, nonlocal: None };
    let s = solve(&p, 100, 1e-9);
    let (qz, fz) = (-1.16654088, -1.16667410);
    println!("  H2: {:.8} with {} functions; aug-cc-pVQZ {qz}, aug-cc-pV5Z {fz}", s.energy, p.basis.size);
    assert!(s.converged);
    assert!(s.energy < qz + 1e-6, "no worse than quadruple zeta");
    assert!(s.energy > fz - 1e-5, "and not below a larger basis, which would mean an error");
}

// ---------------------------------------------------------------------------
// Stage E5: forces
// ---------------------------------------------------------------------------

fn small_molecule(atoms: &[(f64, [f64; 3])]) -> Problem {
    let mut sh = Vec::new();
    let mut ax = Vec::new();
    let mut sizes = Vec::new();
    for (z, c) in atoms {
        if *z > 1.5 {
            sh.extend(even_tempered(*c, 0, 0.1, 3.0, 8));
            sh.extend(even_tempered(*c, 1, 0.15, 3.0, 4));
            sh.push(Shell::primitive(*c, 2, 0.8));
            sizes.push(1.0);
        } else {
            sh.extend(even_tempered(*c, 0, 0.06, 3.0, 5));
            sh.push(Shell::primitive(*c, 1, 0.7));
            sizes.push(0.6);
        }
        for l in 0..=2 {
            ax.extend(even_tempered(*c, l, 0.12, 3.0, if l == 0 { 9 } else { 7 }));
        }
    }
    let ne: f64 = atoms.iter().map(|a| a.0).sum();
    Problem { basis: Basis::new(sh), nuclei: atoms.to_vec(), sizes, alpha: ne / 2.0, beta: ne / 2.0, functional: Functional::Pbe, radial: 50, theta: 12, auxiliary: Some(Basis::new(ax)), prune: false, guess: None, nonlocal: None }
}

/// H2: the analytic force against finite differences of the whole
/// self-consistent energy. Measured 4e-6, which is the differences' own noise
/// (their x and y components, zero by symmetry, come out 3e-6).
#[test]
fn the_force_on_a_stretched_bond_is_the_slope_of_its_energy() {
    use phys::electrons::gradient::gradient;
    let atoms = vec![(1.0, [0.0, 0.0, -0.8]), (1.0, [0.0, 0.0, 0.8])];
    let p = small_molecule(&atoms);
    let s = solve(&p, 200, 1e-11);
    let g = gradient(&p, &s);
    let h = 1e-3;
    let e = |dz: f64| {
        let mut a = atoms.clone();
        a[1].1[2] += dz;
        solve(&small_molecule(&a), 200, 1e-11).energy
    };
    let fd = (e(h) - e(-h)) / (2.0 * h);
    println!("  H2 stretched to 1.6 bohr: analytic {:+.7}, finite difference {fd:+.7}", g[1][2]);
    assert!((g[1][2] - fd).abs() < 2e-5);
    assert!(g[1][2] > 0.0, "a bond stretched past its length pulls back");
    assert!((g[0][2] + g[1][2]).abs() < 1e-10, "and pulls both ends equally");
}

/// Each term of the force is the derivative of its own energy at a fixed
/// density — which is how the gradient was debugged, and what would catch it
/// going wrong. The fitted Coulomb energy is the robust form
/// (`Fitted::coulomb_and_energy`): written as `1/2 b.c` it scattered by
/// 2.2e-7 Ha between geometries 2e-3 bohr apart, which once forced this check
/// onto a larger step and a looser bound and then failed it on round-off
/// alone. And the forces sum to zero: translating the molecule does nothing,
/// which held only once the grid's own motion was included.
#[test]
fn every_term_of_the_force_is_the_slope_of_its_energy() {
    use phys::electrons::gradient::gradient_parts;
    use phys::electrons::integrals::{one_electron, one_electron_gradient};
    use phys::electrons::linalg::Matrix;
    use phys::electrons::scf::{exchange_correlation, Batches, Fitted};
    let atoms = vec![(8.0, [0.05, -0.03, 0.02]), (1.0, [0.0, 1.43, 1.1]), (1.0, [0.1, -1.43, 1.0])];
    let p = small_molecule(&atoms);
    let s = solve(&p, 200, 1e-11);
    let n = p.basis.size;
    let mut d = s.density_alpha.clone();
    let mut w = s.weighted_alpha.clone();
    for k in 0..n * n {
        d.a[k] += s.density_beta.a[k];
        w.a[k] += s.weighted_beta.a[k];
    }
    let owner: Vec<usize> = p.basis.shells.iter().map(|sh| atoms.iter().position(|(_, c)| (0..3).all(|k| (c[k] - sh.centre[k]).abs() < 1e-10)).unwrap()).collect();
    let zero = Matrix::zeros(n);
    let gh = one_electron_gradient(&p.basis, &p.nuclei, &owner, &d, &zero);
    let gs = one_electron_gradient(&p.basis, &p.nuclei, &owner, &zero, &w);
    let parts = gradient_parts(&p, &s);
    for (name, part) in ["nuclear", "one-electron", "coulomb", "xc"].iter().zip(&parts) {
        let sum: f64 = (0..3).map(|k| part.iter().map(|v| v[k]).sum::<f64>().abs()).sum();
        println!("  {name}: forces sum to {sum:.1e}");
        assert!(sum < 1e-9, "{name}");
    }
    let terms = |q: &Problem| -> [f64; 4] {
        let (sm, t, v) = one_electron(&q.basis, &q.nuclei);
        let mut hm = t.clone();
        for k in 0..n * n {
            hm.a[k] += v.a[k];
        }
        let fit = Fitted::new(&q.basis, q.auxiliary.as_ref().unwrap(), &sm);
        let at: Vec<([f64; 3], f64)> = q.nuclei.iter().zip(&q.sizes).map(|((_, pp), r)| (*pp, *r)).collect();
        let grid = phys::electrons::grid::molecular_pruned(&at, q.radial, q.theta, q.prune);
        let (exc, _, _) = exchange_correlation(&q.basis, &Batches::new(&q.basis, &grid), q.functional, &s.density_alpha, &s.density_beta);
        [hm.dot(&d), -sm.dot(&w), fit.energy(&d), exc]
    };
    let (a, dd) = (1usize, 0usize);
    let fd = |h: f64| {
        let mut pl = atoms.clone();
        pl[a].1[dd] += h;
        let mut mi = atoms.clone();
        mi[a].1[dd] -= h;
        let (tp, tm) = (terms(&small_molecule(&pl)), terms(&small_molecule(&mi)));
        (0..4).map(|k| (tp[k] - tm[k]) / (2.0 * h)).collect::<Vec<f64>>()
    };
    let f4 = fd(1e-4);
    let checks = [("D dH", gh[a][dd], f4[0], 1e-6), ("-W dS", gs[a][dd], f4[1], 1e-6), ("Coulomb", parts[2][a][dd], f4[2], 1e-6), ("XC", parts[3][a][dd], f4[3], 1e-6)];
    for (name, analytic, numeric, tol) in checks {
        println!("  {name}: analytic {analytic:+.8}, finite difference {numeric:+.8}");
        assert!((analytic - numeric).abs() < tol, "{name}");
    }
}

/// The three-centre integrals reorganised for speed are the same numbers as
/// the general four-centre routine with a unit function, and their contracted
/// derivative is the same as contracting the general derivative — over every
/// pairing of angular momenta up to f on the bra and g on the auxiliary side.
#[test]
fn three_centre_integrals_match_the_general_routine() {
    use phys::electrons::integrals::{eri_derivative, eri_three_block, eri_three_gradient_block};
    let unit = Shell::unit();
    let mut worst = 0.0f64;
    let mut worst_g = 0.0f64;
    for la in 0..=3 {
        for lb in 0..=2 {
            for lp in 0..=4 {
                let a = Shell::contracted([0.1, -0.2, 0.3], la, vec![3.0, 0.7], vec![0.6, 0.5]);
                let b = Shell::primitive([-0.5, 0.4, 0.9], lb, 1.3);
                let p = Shell::primitive([0.3, 0.8, -0.4], lp, 0.9);
                let fast = eri_three_block(&a, &b, &p);
                let slow = eri_block(&a, &b, &p, &unit);
                for (x, y) in fast.iter().zip(&slow) {
                    worst = worst.max((x - y).abs());
                }
                // A density and coefficients with no structure.
                let d: Vec<f64> = (0..a.size() * b.size()).map(|i| ((i * 7 + 3) % 11) as f64 / 11.0 - 0.4).collect();
                let c: Vec<f64> = (0..p.size()).map(|i| ((i * 5 + 2) % 7) as f64 / 7.0 - 0.3).collect();
                let (ga, gb) = eri_three_gradient_block(&a, &b, &p, &d, &c);
                let da = eri_derivative(&a, &b, &p, &unit);
                let db = eri_derivative(&b, &a, &p, &unit);
                let (na, nb, np) = (a.size(), b.size(), p.size());
                for dir in 0..3 {
                    let (mut ra, mut rb) = (0.0, 0.0);
                    for i in 0..na {
                        for j in 0..nb {
                            for k in 0..np {
                                ra += d[i * nb + j] * c[k] * da[dir][(i * nb + j) * np + k];
                                rb += d[i * nb + j] * c[k] * db[dir][(j * na + i) * np + k];
                            }
                        }
                    }
                    worst_g = worst_g.max((ga[dir] - ra).abs()).max((gb[dir] - rb).abs());
                }
            }
        }
    }
    println!("three-centre: worst {worst:.1e}, gradient worst {worst_g:.1e}");
    assert!(worst < 1e-12, "three-centre integrals differ by {worst:e}");
    assert!(worst_g < 1e-11, "three-centre gradient differs by {worst_g:e}");
}

/// The blocked product is the plain triple loop to the last bit — on shapes
/// that are not multiples of its block, and whichever vector width the
/// processor chose — because it sums each element in the same order.
#[test]
fn the_blocked_product_is_the_plain_loop_exactly() {
    use phys::electrons::linalg::product_nt;
    for (m, n, kd) in [(1, 1, 1), (3, 5, 7), (4, 4, 4), (9, 13, 31), (17, 6, 256)] {
        let a: Vec<f64> = (0..m * kd).map(|i| ((i * 37 + 11) % 101) as f64 / 7.0 - 6.0).collect();
        let b: Vec<f64> = (0..n * kd).map(|i| ((i * 53 + 5) % 97) as f64 / 3.0 - 15.0).collect();
        let mut fast = vec![0.0; m * n];
        product_nt(&a, &b, m, n, kd, &mut fast);
        for i in 0..m {
            for j in 0..n {
                let mut t = 0.0;
                for p in 0..kd {
                    t += a[i * kd + p] * b[j * kd + p];
                }
                assert_eq!(fast[i * n + j].to_bits(), t.to_bits(), "{m}x{n}x{kd} at ({i},{j})");
            }
        }
    }
}

/// The fitted density's potential between two functions is the three-centre
/// integrals contracted with the coefficients — checked against contracting
/// the general routine's integrals by hand.
#[test]
fn a_fitted_potential_is_the_contracted_three_centre_integrals() {
    use phys::electrons::integrals::eri_three_contracted_block;
    let unit = Shell::unit();
    let mut worst = 0.0f64;
    for la in 0..=3 {
        for lb in 0..=2 {
            for lp in 0..=4 {
                let a = Shell::contracted([0.1, -0.2, 0.3], la, vec![3.0, 0.7], vec![0.6, 0.5]);
                let b = Shell::primitive([-0.5, 0.4, 0.9], lb, 1.3);
                let p = Shell::primitive([0.3, 0.8, -0.4], lp, 0.9);
                let c: Vec<f64> = (0..p.size()).map(|i| ((i * 5 + 2) % 7) as f64 / 7.0 - 0.3).collect();
                let fast = eri_three_contracted_block(&a, &b, &p, &c);
                let slow = eri_block(&a, &b, &p, &unit);
                let (nb, np) = (b.size(), p.size());
                for i in 0..a.size() {
                    for j in 0..nb {
                        let v: f64 = (0..np).map(|k| c[k] * slow[(i * nb + j) * np + k]).sum();
                        worst = worst.max((fast[i * nb + j] - v).abs());
                    }
                }
            }
        }
    }
    println!("fitted potential: worst {worst:.1e}");
    assert!(worst < 1e-12, "{worst:e}");
}

/// The non-local correlation's force is the slope of its own energy at a
/// fixed density, on the grid the field was solved with — basis functions,
/// points and partition all moving with their atoms, and the points moving
/// apart, which the semilocal term does not have. And its forces sum to zero.
#[test]
fn the_nonlocal_force_is_the_slope_of_its_energy() {
    use phys::electrons::gradient::gradient_parts;
    use phys::electrons::scf::Batches;
    use phys::electrons::vdw::{density_and_gradient_on, kernel_table, nonlocal, NonlocalSpec, Z_AB_DF1};
    let atoms = vec![(8.0, [0.05, -0.03, 0.02]), (1.0, [0.0, 1.43, 1.1]), (1.0, [0.1, -1.43, 1.0])];
    let with_nl = |a: &[(f64, [f64; 3])]| {
        let mut p = small_molecule(a);
        p.functional = Functional::PbeXLdaC;
        p.nonlocal = Some(NonlocalSpec::in_the_field(Z_AB_DF1));
        p
    };
    let p = with_nl(&atoms);
    let s = solve(&p, 200, 1e-11);
    let mut d = s.density_alpha.clone();
    for k in 0..d.a.len() {
        d.a[k] += s.density_beta.a[k];
    }
    let nl_force = gradient_parts(&p, &s)[4].clone();
    let sum: f64 = (0..3).map(|k| nl_force.iter().map(|v| v[k]).sum::<f64>().abs()).sum();
    println!("  non-local forces sum to {sum:.1e}");
    assert!(sum < 1e-9, "the non-local forces do not sum to zero: {sum}");
    let energy_at = |a: &[(f64, [f64; 3])]| {
        let q = with_nl(a);
        let spec = q.nonlocal.unwrap();
        let at: Vec<([f64; 3], f64)> = q.nuclei.iter().zip(&q.sizes).map(|((_, pp), r)| (*pp, *r)).collect();
        let grid = phys::electrons::grid::molecular_pruned(&at, spec.radial, spec.theta, false);
        let b = Batches::new(&q.basis, &grid);
        nonlocal(&grid, &density_and_gradient_on(&q.basis, &b, &d, grid.points.len()), spec.z_ab, spec.floor, kernel_table()).energy
    };
    for (a, dd) in [(1usize, 0usize), (0, 2), (2, 1)] {
        let h = 1e-4;
        let mut pl = atoms.clone();
        pl[a].1[dd] += h;
        let mut mi = atoms.clone();
        mi[a].1[dd] -= h;
        let fd = (energy_at(&pl) - energy_at(&mi)) / (2.0 * h);
        println!("  atom {a}, axis {dd}: analytic {:+.9}, finite difference {fd:+.9}", nl_force[a][dd]);
        assert!((nl_force[a][dd] - fd).abs() < 1e-6, "atom {a} axis {dd}: {} against {fd}", nl_force[a][dd]);
    }
}

/// The whole force with non-local correlation in the field, against the slope
/// of the self-consistent energy: the overlap term reaches the non-local
/// potential through the Fock matrix, and this is what checks it does.
#[test]
fn the_force_with_nonlocal_correlation_is_the_slope_of_its_energy() {
    use phys::electrons::gradient::gradient;
    use phys::electrons::vdw::{NonlocalSpec, Z_AB_DF1};
    let atoms = vec![(8.0, [0.05, -0.03, 0.02]), (1.0, [0.0, 1.43, 1.1]), (1.0, [0.1, -1.43, 1.0])];
    let with_nl = |a: &[(f64, [f64; 3])]| {
        let mut p = small_molecule(a);
        p.functional = Functional::PbeXLdaC;
        p.nonlocal = Some(NonlocalSpec::in_the_field(Z_AB_DF1));
        p
    };
    let p = with_nl(&atoms);
    let s = solve(&p, 200, 1e-11);
    let g = gradient(&p, &s);
    let (a, dd, h) = (1usize, 1usize, 1e-3);
    let e = |dz: f64| {
        let mut m = atoms.clone();
        m[a].1[dd] += dz;
        solve(&with_nl(&m), 200, 1e-11).energy
    };
    let fd = (e(h) - e(-h)) / (2.0 * h);
    println!("  water with non-local correlation: analytic {:+.7}, finite difference {fd:+.7}", g[a][dd]);
    assert!((g[a][dd] - fd).abs() < 2e-5, "{} against {fd}", g[a][dd]);
}

/// The fitted Coulomb table stores what fits in its memory budget and builds
/// the rest again when it is visited, in the same order: the coefficients and
/// the Coulomb matrix come out bit-identical whether all of it, half of it or
/// none of it is stored, and however small the slices it is rebuilt in.
#[test]
fn the_coulomb_fit_is_the_same_however_much_of_it_is_stored() {
    use phys::electrons::scf::Fitted;
    let atoms = vec![(8.0, [0.05, -0.03, 0.02]), (1.0, [0.0, 1.43, 1.1]), (1.0, [0.1, -1.43, 1.0])];
    let p = small_molecule(&atoms);
    let s = solve(&p, 200, 1e-11);
    let mut d = s.density_alpha.clone();
    for k in 0..d.a.len() {
        d.a[k] += s.density_beta.a[k];
    }
    let (sm, _, _) = phys::electrons::integrals::one_electron(&p.basis, &p.nuclei);
    let aux = p.auxiliary.as_ref().unwrap();
    let all = Fitted::new_within(&p.basis, aux, &sm, usize::MAX);
    let (c_all, (j_all, e_all)) = (all.coefficients(&d), all.coulomb_and_energy(&d));
    for (budget, slice) in [(0usize, usize::MAX), (0, 1), (200_000, 50_000)] {
        let mut part = Fitted::new_within(&p.basis, aux, &sm, budget);
        part.set_slice_bytes(slice);
        let (c, (j, e)) = (part.coefficients(&d), part.coulomb_and_energy(&d));
        println!("  budget {budget} B, slices of {slice} B: {:.0}% stored; E_J {e:.15} against {e_all:.15}", 100.0 * part.stored_fraction());
        assert!(c.iter().zip(&c_all).all(|(a, b)| a.to_bits() == b.to_bits()), "the coefficients differ");
        assert!(j.a.iter().zip(&j_all.a).all(|(a, b)| a.to_bits() == b.to_bits()), "the Coulomb matrix differs");
        assert_eq!(e.to_bits(), e_all.to_bits());
    }
    assert_eq!(all.stored_fraction(), 1.0);
}

/// A table spilled to disk is read back to the same bits as one held in
/// memory, whatever is left in memory and however it is sliced; its file is
/// deleted when the table is dropped; a table over the cap or a disk that
/// cannot be written to leaves the fit exactly as if there were no disk; and
/// the three solves of a counterpoise correction share one table.
#[test]
fn a_spilled_table_is_the_same_bits_and_is_shared() {
    use phys::electrons::scf::{FitCache, Fitted, SpillTo};
    let atoms = vec![(8.0, [0.05, -0.03, 0.02]), (1.0, [0.0, 1.43, 1.1]), (1.0, [0.1, -1.43, 1.0])];
    let p = small_molecule(&atoms);
    let s = solve(&p, 200, 1e-11);
    let mut d = s.density_alpha.clone();
    for k in 0..d.a.len() {
        d.a[k] += s.density_beta.a[k];
    }
    let (sm, _, _) = phys::electrons::integrals::one_electron(&p.basis, &p.nuclei);
    let aux = p.auxiliary.as_ref().unwrap();
    let all = Fitted::new_within(&p.basis, aux, &sm, usize::MAX);
    let (c_all, (j_all, e_all)) = (all.coefficients(&d), all.coulomb_and_energy(&d));
    let dir = std::env::temp_dir().join(format!("phys-spill-test-{}", std::process::id()));
    let to = SpillTo { dir: dir.clone(), cap_bytes: u64::MAX };
    for (budget, slice) in [(0usize, usize::MAX), (0, 1), (200_000, 50_000)] {
        let mut part = Fitted::new_within_spill(&p.basis, aux, &sm, budget, Some(&to));
        part.set_slice_bytes(slice);
        assert!(part.is_spilled(), "the rest went to disk");
        let (c, (j, e)) = (part.coefficients(&d), part.coulomb_and_energy(&d));
        println!("  budget {budget} B, slices of {slice} B: {:.0}% in memory, the rest on disk; E_J {e:.15}", 100.0 * part.stored_fraction());
        assert!(c.iter().zip(&c_all).all(|(a, b)| a.to_bits() == b.to_bits()), "the coefficients differ");
        assert!(j.a.iter().zip(&j_all.a).all(|(a, b)| a.to_bits() == b.to_bits()), "the Coulomb matrix differs");
        assert_eq!(e.to_bits(), e_all.to_bits());
    }
    let files = || std::fs::read_dir(&dir).map(|r| r.count()).unwrap_or(0);
    assert_eq!(files(), 0, "a dropped table leaves no file");
    // Over the cap: no file, same answer.
    let capped = SpillTo { dir: dir.clone(), cap_bytes: 10 };
    let c10 = Fitted::new_within_spill(&p.basis, aux, &sm, 0, Some(&capped));
    assert!(!c10.is_spilled() && files() == 0, "a table over the cap is not spilled");
    assert!(c10.coefficients(&d).iter().zip(&c_all).all(|(a, b)| a.to_bits() == b.to_bits()));
    // A disk that cannot be written to: a path under a file.
    let blocker = std::env::temp_dir().join(format!("phys-spill-blocker-{}", std::process::id()));
    std::fs::write(&blocker, b"x").unwrap();
    let bad = SpillTo { dir: blocker.join("under"), cap_bytes: u64::MAX };
    let failed = Fitted::new_within_spill(&p.basis, aux, &sm, 0, Some(&bad));
    assert!(!failed.is_spilled(), "a failed write leaves no spill");
    assert!(failed.coefficients(&d).iter().zip(&c_all).all(|(a, b)| a.to_bits() == b.to_bits()));
    let _ = std::fs::remove_file(&blocker);
    let _ = std::fs::remove_dir_all(&dir);
    // Sharing: the same basis again is the same table, another is not.
    let cache = FitCache::new();
    let first = Fitted::new_cached(&cache, &p.basis, aux, &sm);
    let again = Fitted::new_cached(&cache, &p.basis, aux, &sm);
    assert!(first.shares_table_with(&again), "the same basis shares its table");
    assert!(first.coefficients(&d).iter().zip(&c_all).all(|(a, b)| a.to_bits() == b.to_bits()));
    let moved = small_molecule(&[(8.0, [0.05, -0.03, 0.021]), (1.0, [0.0, 1.43, 1.1]), (1.0, [0.1, -1.43, 1.0])]);
    let other = Fitted::new_cached(&cache, &moved.basis, moved.auxiliary.as_ref().unwrap(), &sm);
    assert!(!first.shares_table_with(&other), "a basis that moved does not");
}
