//! Liquids: what the chemistry derives about a molecule, against what is measured.
//!
//! One dataset serves boiling point, heat of vaporisation, density, and how
//! density moves with temperature, because they are one question: how tightly
//! do these molecules hold each other.
use phys::chem::analyse::analyse;
use phys::chem::arrange::{Arrangement, Bond, Order};
use phys::chem::elements::Element;
use Order::*;

pub fn mol(heavy: &[u8], bonds: &[(usize, usize, Order)]) -> Arrangement {
    let mut atoms: Vec<Element> = heavy.iter().map(|z| Element(*z)).collect();
    let mut bs: Vec<Bond> = bonds.iter().map(|(a, b, o)| Bond::new(*a, *b, *o)).collect();
    for i in 0..heavy.len() {
        let used: usize = bonds.iter().filter(|(a, b, _)| *a == i || *b == i).map(|(_, _, o)| o.slots()).sum();
        let v = atoms[i].valence().unwrap();
        for _ in used..v {
            atoms.push(Element(1));
            bs.push(Bond::new(i, atoms.len() - 1, Order::Single));
        }
    }
    Arrangement::molecule(atoms, bs)
}

fn chain(n: usize) -> Vec<(usize, usize, Order)> {
    (0..n - 1).map(|i| (i, i + 1, Single)).collect()
}
fn ring(n: usize, kek: bool) -> Vec<(usize, usize, Order)> {
    (0..n).map(|i| (i, (i + 1) % n, if kek && i % 2 == 0 { Double } else { Single })).collect()
}

/// A molecule and what is known of it.
pub struct Liquid {
    pub name: &'static str,
    pub arr: Arrangement,
    /// Normal boiling point, K.
    pub tb: f64,
    /// Enthalpy of vaporisation at the boiling point, kJ/mol.
    pub dh: f64,
    /// Density, kg/m^3, at `t_rho`.
    pub rho: f64,
    pub t_rho: f64,
}

pub fn liquids() -> Vec<Liquid> {
    let l = |name, arr, tb, dh, rho, t_rho| Liquid { name, arr, tb, dh, rho, t_rho };
    let at20 = 293.15;
    let mut tol = ring(6, true); tol.push((0, 6, Single));
    let mut phenol = ring(6, true); phenol.push((0, 6, Single));
    let mut chlorob = ring(6, true); chlorob.push((0, 6, Single));
    vec![
        l("water", mol(&[8], &[]), 373.15, 40.65, 998.2, at20),
        l("methane", mol(&[6], &[]), 111.7, 8.19, 422.0, 111.7),
        l("ethane", mol(&[6, 6], &chain(2)), 184.6, 14.7, 546.0, 184.6),
        l("ethylene", mol(&[6, 6], &[(0, 1, Double)]), 169.4, 13.5, 566.0, 169.4),
        l("propane", mol(&[6; 3], &chain(3)), 231.1, 18.8, 582.0, 231.1),
        l("isobutane", mol(&[6; 4], &[(0, 1, Single), (1, 2, Single), (1, 3, Single)]), 261.4, 21.3, 557.0, 261.4),
        l("butane", mol(&[6; 4], &chain(4)), 272.7, 22.4, 601.0, 272.7),
        l("pentane", mol(&[6; 5], &chain(5)), 309.2, 25.8, 626.0, at20),
        l("cyclopentane", mol(&[6; 5], &ring(5, false)), 322.4, 27.3, 745.0, at20),
        l("hexane", mol(&[6; 6], &chain(6)), 341.9, 28.9, 655.0, at20),
        l("cyclohexane", mol(&[6; 6], &ring(6, false)), 353.9, 29.9, 779.0, at20),
        l("heptane", mol(&[6; 7], &chain(7)), 371.6, 31.8, 684.0, at20),
        l("octane", mol(&[6; 8], &chain(8)), 398.8, 34.4, 703.0, at20),
        l("benzene", mol(&[6; 6], &ring(6, true)), 353.2, 30.7, 876.0, at20),
        l("toluene", mol(&[6; 7], &tol), 383.8, 33.2, 867.0, at20),
        l("chlorobenzene", mol(&[6, 6, 6, 6, 6, 6, 17], &chlorob), 404.9, 35.2, 1106.0, at20),
        l("ammonia", mol(&[7], &[]), 239.8, 23.3, 682.0, 239.8),
        l("methanol", mol(&[6, 8], &chain(2)), 337.7, 35.2, 792.0, at20),
        l("ethanol", mol(&[6, 6, 8], &chain(3)), 351.4, 38.6, 789.0, at20),
        l("propanol", mol(&[6, 6, 6, 8], &chain(4)), 370.4, 41.4, 803.0, at20),
        l("isopropanol", mol(&[6, 6, 6, 8], &[(0, 1, Single), (1, 2, Single), (1, 3, Single)]), 355.4, 39.9, 786.0, at20),
        l("butanol", mol(&[6, 6, 6, 6, 8], &chain(5)), 390.9, 43.3, 810.0, at20),
        l("phenol", mol(&[6, 6, 6, 6, 6, 6, 8], &phenol), 455.0, 46.0, 1058.0, 323.15),
        l("ethylene glycol", mol(&[8, 6, 6, 8], &chain(4)), 470.4, 56.9, 1110.0, at20),
        l("glycerol", mol(&[8, 6, 6, 6, 8, 8], &[(0, 1, Single), (1, 2, Single), (2, 3, Single), (3, 4, Single), (2, 5, Single)]), 563.2, 61.0, 1261.0, at20),
        l("dimethyl ether", mol(&[6, 8, 6], &chain(3)), 248.3, 21.5, 735.0, 248.3),
        l("diethyl ether", mol(&[6, 6, 8, 6, 6], &chain(5)), 307.7, 26.5, 713.0, at20),
        l("tetrahydrofuran", mol(&[8, 6, 6, 6, 6], &ring(5, false)), 339.1, 29.9, 889.0, at20),
        l("acetaldehyde", mol(&[6, 6, 8], &[(0, 1, Single), (1, 2, Double)]), 293.3, 25.7, 784.0, 293.3),
        l("acetone", mol(&[6, 6, 6, 8], &[(0, 1, Single), (1, 2, Single), (1, 3, Double)]), 329.2, 29.1, 784.0, at20),
        l("butanone", mol(&[6, 6, 6, 6, 8], &[(0, 1, Single), (1, 2, Single), (2, 3, Single), (1, 4, Double)]), 352.8, 31.3, 805.0, at20),
        l("ethyl acetate", mol(&[6, 6, 8, 8, 6, 6], &[(0, 1, Single), (1, 2, Double), (1, 3, Single), (3, 4, Single), (4, 5, Single)]), 350.2, 31.9, 902.0, at20),
        l("formic acid", mol(&[6, 8, 8], &[(0, 1, Double), (0, 2, Single)]), 373.7, 22.7, 1220.0, at20),
        l("acetic acid", mol(&[6, 6, 8, 8], &[(0, 1, Single), (1, 2, Double), (1, 3, Single)]), 391.1, 23.7, 1049.0, at20),
        l("acetonitrile", mol(&[6, 6, 7], &[(0, 1, Single), (1, 2, Triple)]), 354.8, 29.8, 786.0, at20),
        l("methylamine", mol(&[6, 7], &chain(2)), 266.8, 25.6, 699.0, 266.8),
        l("ethylamine", mol(&[6, 6, 7], &chain(3)), 289.7, 27.4, 683.0, 289.7),
        l("aniline", mol(&[6, 6, 6, 6, 6, 6, 7], &phenol), 457.2, 42.4, 1022.0, at20),
        l("pyridine", mol(&[7, 6, 6, 6, 6, 6], &ring(6, true)), 388.4, 35.1, 982.0, at20),
        l("hydrazine", mol(&[7, 7], &chain(2)), 387.0, 41.8, 1021.0, at20),
        l("hydrogen peroxide", mol(&[8, 8], &chain(2)), 423.4, 51.6, 1450.0, at20),
        l("dichloromethane", mol(&[6, 17, 17], &[(0, 1, Single), (0, 2, Single)]), 313.0, 28.1, 1325.0, at20),
        l("chloroform", mol(&[6, 17, 17, 17], &[(0, 1, Single), (0, 2, Single), (0, 3, Single)]), 334.3, 29.2, 1489.0, at20),
        l("carbon tetrachloride", mol(&[6, 17, 17, 17, 17], &[(0, 1, Single), (0, 2, Single), (0, 3, Single), (0, 4, Single)]), 349.9, 29.8, 1594.0, at20),
    ]
}

/// Density along a liquid's own range: (T, kg/m^3).
pub fn series() -> Vec<(&'static str, Vec<(f64, f64)>)> {
    vec![
        ("water", vec![(273.15, 999.8), (293.15, 998.2), (323.15, 988.0), (348.15, 974.9), (373.15, 958.4)]),
        ("ethanol", vec![(273.15, 806.3), (293.15, 789.3), (323.15, 763.5), (343.15, 744.0)]),
        ("methanol", vec![(273.15, 810.0), (293.15, 791.4), (323.15, 763.6)]),
        ("hexane", vec![(273.15, 677.0), (293.15, 659.0), (323.15, 632.0)]),
        ("benzene", vec![(278.7, 894.0), (293.15, 876.5), (323.15, 847.0), (343.15, 825.0)]),
        ("acetone", vec![(273.15, 812.0), (293.15, 789.9), (323.15, 752.0)]),
        ("carbon tetrachloride", vec![(273.15, 1632.7), (293.15, 1594.0), (323.15, 1539.0)]),
    ]
}

pub fn vdw_volume_and_area(a: &Arrangement) -> (f64, f64) {
    let pi = std::f64::consts::PI;
    let (mut v, mut s) = (0.0, 0.0);
    for e in &a.atoms {
        let r = e.vdw_radius().unwrap();
        v += 4.0 / 3.0 * pi * r.powi(3);
        s += 4.0 * pi * r * r;
    }
    for b in &a.bonds {
        let (ea, eb) = (a.atoms[b.a as usize], a.atoms[b.b as usize]);
        let (ri, rj) = (ea.vdw_radius().unwrap(), eb.vdw_radius().unwrap());
        let d = ea.covalent_radius().unwrap() + eb.covalent_radius().unwrap();
        let h = |r: f64, s: f64| (r - (r * r - s * s + d * d) / (2.0 * d)).clamp(0.0, 2.0 * r);
        let (hi, hj) = (h(ri, rj), h(rj, ri));
        let cap = |r: f64, h: f64| pi * h * h * (3.0 * r - h) / 3.0;
        v -= cap(ri, hi) + cap(rj, hj);
        s -= 2.0 * pi * ri * hi + 2.0 * pi * rj * hj;
    }
    (v, s)
}

fn relative(derived: f64, measured: f64) -> f64 {
    (derived / measured - 1.0).abs()
}

/// A boiling point from the arrangement and nothing else.
///
/// The law has four coefficients fitted to these same molecules, so the error
/// here is in-sample; the error that predicts a molecule the fit never saw is
/// 7.3% on average and 26% at worst (leave-one-out, recorded on
/// `analyse::phase_points`). The law it replaced was 35% and 112% over the same
/// set. The bounds sit just outside what it is now.
#[test]
fn a_boiling_point_is_derived_from_the_arrangement() {
    let set = liquids();
    let (mut total, mut worst) = (0.0, (0.0f64, ""));
    for x in &set {
        let p = analyse(&x.arr).unwrap();
        let e = relative(p.boiling_point, x.tb);
        println!("  {:22} measured {:4.0} derived {:4.0} ({:+.0}%)", x.name, x.tb, p.boiling_point, 100.0 * (p.boiling_point / x.tb - 1.0));
        total += e;
        if e > worst.0 {
            worst = (e, x.name);
        }
    }
    let mean = total / set.len() as f64;
    println!("  mean {mean:.3}, worst {:.3} ({})", worst.0, worst.1);
    assert!(mean < 0.10, "mean error {mean:.3}");
    assert!(worst.0 < 0.32, "{} is out by {:.3}", worst.1, worst.0);
}

/// A liquid's density at the temperature it was measured at.
#[test]
fn a_density_is_derived_from_geometry_and_the_boiling_point() {
    let set = liquids();
    let (mut total, mut worst) = (0.0, (0.0f64, ""));
    for x in &set {
        let p = analyse(&x.arr).unwrap();
        let derived = p.density_at(x.t_rho);
        let e = relative(derived, x.rho);
        println!("  {:22} at {:5.1} K measured {:5.0} derived {:5.0} ({:+.0}%)", x.name, x.t_rho, x.rho, derived, 100.0 * (derived / x.rho - 1.0));
        total += e;
        if e > worst.0 {
            worst = (e, x.name);
        }
    }
    let mean = total / set.len() as f64;
    println!("  mean {mean:.3}, worst {:.3} ({})", worst.0, worst.1);
    assert!(mean < 0.10, "mean error {mean:.3}");
    assert!(worst.0 < 0.30, "{} is out by {:.3}", worst.1, worst.0);
}

/// How density moves with temperature, which is the part a single stored
/// number cannot do. Seven liquids, each measured at three to five
/// temperatures; the error is in the *change* from 293 K, so it is not the
/// error of the density itself.
#[test]
fn a_density_follows_the_temperature() {
    let (mut total, mut count, mut worst) = (0.0, 0.0, (0.0f64, "", 0.0));
    for (name, points) in series() {
        let x = liquids().into_iter().find(|l| l.name == name).unwrap();
        let p = analyse(&x.arr).unwrap();
        let at = |t: f64| points.iter().find(|(tt, _)| (tt - t).abs() < 1e-6).map(|(_, r)| *r).unwrap();
        let (r0, d0) = (at(293.15), p.density_at(293.15));
        for (t, r) in &points {
            let measured = r / r0;
            let derived = p.density_at(*t) / d0;
            let e = (derived - measured).abs();
            println!("  {name:22} {t:6.1} K  measured {measured:.4} of its 293 K value, derived {derived:.4}");
            total += e;
            count += 1.0;
            if e > worst.0 {
                worst = (e, name, *t);
            }
        }
    }
    println!("  mean {:.4}, worst {:.4} ({} at {} K)", total / count, worst.0, worst.1, worst.2);
    assert!(total / count < 0.015, "mean error in the change {:.4}", total / count);
    assert!(worst.0 < 0.06, "{} at {} K is out by {:.4}", worst.1, worst.2, worst.0);
}

/// The case that started this: water, which the engine had at 1653 kg/m^3
/// and then, with the overlap corrected, at 847.
#[test]
fn water_is_about_as_dense_as_water() {
    let w = analyse(&mol(&[8], &[])).unwrap();
    println!("  water: {:.0} at 293 K (998), {:.0} at 373 K (958), boils at {:.0} K (373), expansion {:.2e} /K (real 4.3e-4 mean: water is the one liquid that expands slowly)", w.density_at(293.15), w.density_at(373.15), w.boiling_point, w.expansion);
    assert!(relative(w.density_at(293.15), 998.2) < 0.06);
    assert!(relative(w.boiling_point, 373.15) < 0.15);
    assert!(w.density_at(373.15) < w.density_at(273.15), "water expands as it warms");
}

/// A crystal expands by the Gruneisen relation from numbers the engine already
/// derives, so it inherits their error: iron's stiffness is low because its
/// cohesive energy is (`substances::iron`), and its expansion is high by the
/// same factor. The bound is a factor of four either way, which is wide and
/// says so; what it catches is a unit slip or a sign, not a few per cent.
#[test]
fn a_crystal_expands_by_the_gruneisen_relation() {
    use phys::material::substances;
    for (name, p, real) in [
        ("iron", substances::iron(), 3.5e-5),
        ("quartz", substances::silica(), 3.5e-5),
        ("calcite", substances::calcium_carbonate(), 1.4e-5),
    ] {
        println!("  {name}: volumetric expansion {:.2e} /K (real about {real:.1e}, {:.1}x)", p.expansion, p.expansion / real);
        assert!(p.expansion > real / 4.0 && p.expansion < real * 4.0, "{name}: {:.2e} against {real:.1e}", p.expansion);
    }
    let hexane = analyse(&liquids().into_iter().find(|l| l.name == "hexane").unwrap().arr).unwrap();
    println!("  hexane: {:.2e} /K (real 1.4e-3)", hexane.expansion);
    assert!(hexane.expansion > 1.0e-3 && hexane.expansion < 1.8e-3, "a liquid expands about 40x a crystal");
}

/// The consumers of density read it at the matter's own temperature. Without
/// this a sea told it is at 373 K rests exactly as heavy as one at 273 K, and
/// the `density_at` law above is correct and unused.
#[test]
fn the_equation_of_state_follows_the_temperature() {
    use phys::chem::{Mixture, Phase, Registry};
    use phys::eos::{Condensed, Eos};
    let mut reg = Registry::new();
    let h2o = reg.intern(mol(&[8], &[])).unwrap();
    let mut mix = Mixture::new();
    mix.add(h2o, Phase::Liquid, 1.0);
    let rest = |t: f64| match Eos::of_mixture(&mix, &reg, t) {
        Eos::Condensed(c) => c.rest_density,
        Eos::Gas => panic!("liquid water is not a gas at {t} K"),
    };
    let (cold, hot) = (rest(280.0), rest(370.0));
    let props = &reg.get(h2o).unwrap().props;
    println!("  water rests at {cold:.1} kg/m^3 at 280 K and {hot:.1} at 370 K (density_at: {:.1}, {:.1})", props.density_at(280.0), props.density_at(370.0));
    assert!(hot < 0.95 * cold, "{hot} at 370 K against {cold} at 280 K");
    assert!(relative(cold, props.density_at(280.0)) < 1e-12);
    // And the single-substance form agrees with the mixture form.
    assert!(relative(Condensed::liquid(props, 370.0).unwrap().rest_density, hot) < 1e-12);
}
