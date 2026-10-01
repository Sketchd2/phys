use phys::chem::analyse::analyse;
use phys::chem::arrange::{Arrangement, Bond, Order};
use phys::chem::elements::Element;

fn mol(heavy: &[u8], bonds: &[(usize, usize, Order)]) -> Arrangement {
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
use Order::*;

/// Bondi: sum of atom volumes less the cap each bond shares with its neighbour.
fn vdw_volume(a: &Arrangement) -> f64 {
    let pi = std::f64::consts::PI;
    let mut v = 0.0;
    for e in &a.atoms { let r = e.vdw_radius().unwrap(); v += 4.0 / 3.0 * pi * r.powi(3); }
    for b in &a.bonds {
        if matches!(b.order, Hydrogen) { continue; }
        let (ea, eb) = (a.atoms[b.a as usize], a.atoms[b.b as usize]);
        let (ri, rj) = (ea.vdw_radius().unwrap(), eb.vdw_radius().unwrap());
        let d = ea.covalent_radius().unwrap() + eb.covalent_radius().unwrap();
        let cap = |r: f64, s: f64| { let h = r - (r * r - s * s + d * d) / (2.0 * d); pi * h * h * (3.0 * r - h) / 3.0 };
        v -= cap(ri, rj) + cap(rj, ri);
    }
    v
}

/// What a molecule's density is derived to be, against what it is measured to be.
///
/// Twenty-four molecules at 293 K, or at the boiling point for those that are
/// gases at 293 K. Before the overlap correction the error was 62% on average
/// and 134% at worst; the bounds below are set just outside what it is now
/// (10% and 26%) so that a regression shows and an improvement is not blocked.
/// The misses that remain are systematic — the high-boiling, hydrogen-bonded
/// liquids come out 15 to 26% light — and track T/T_b: given a true boiling
/// point a one-term correction brings the mean to 7%, but this module's own
/// boiling points are 0.7 to 2.5 times out, so the correction waits on them.
#[test]
fn a_molecules_density_is_derived_from_its_geometry() {
    let chain = |n: usize| -> Vec<(usize, usize, Order)> { (0..n - 1).map(|i| (i, i + 1, Single)).collect() };
    let ring = |n: usize, kek: bool| -> Vec<(usize, usize, Order)> {
        (0..n).map(|i| (i, (i + 1) % n, if kek && i % 2 == 0 { Double } else { Single })).collect()
    };
    let mut tol = ring(6, true); tol.push((0, 6, Single));
    // name, arrangement, measured density kg/m3, measured boiling point K
    let set: Vec<(&str, Arrangement, f64, f64)> = vec![
        ("water", mol(&[8], &[]), 998.0, 373.0),
        ("methane", mol(&[6], &[]), 422.0, 112.0),
        ("ethane", mol(&[6, 6], &chain(2)), 546.0, 184.0),
        ("ammonia", mol(&[7], &[]), 682.0, 240.0),
        ("methanol", mol(&[6, 8], &chain(2)), 792.0, 338.0),
        ("ethanol", mol(&[6, 6, 8], &chain(3)), 789.0, 351.0),
        ("propanol", mol(&[6, 6, 6, 8], &chain(4)), 803.0, 370.0),
        ("pentane", mol(&[6; 5], &chain(5)), 626.0, 309.0),
        ("hexane", mol(&[6; 6], &chain(6)), 655.0, 342.0),
        ("octane", mol(&[6; 8], &chain(8)), 703.0, 399.0),
        ("cyclohexane", mol(&[6; 6], &ring(6, false)), 779.0, 354.0),
        ("ether", mol(&[6, 6, 8, 6, 6], &chain(5)), 713.0, 308.0),
        ("acetone", mol(&[6, 6, 6, 8], &[(0, 1, Single), (1, 2, Single), (1, 3, Double)]), 784.0, 329.0),
        ("acetic acid", mol(&[6, 6, 8, 8], &[(0, 1, Single), (1, 2, Double), (1, 3, Single)]), 1049.0, 391.0),
        ("formic acid", mol(&[6, 8, 8], &[(0, 1, Double), (0, 2, Single)]), 1220.0, 374.0),
        ("acetonitrile", mol(&[6, 6, 7], &[(0, 1, Single), (1, 2, Triple)]), 786.0, 355.0),
        ("benzene", mol(&[6; 6], &ring(6, true)), 876.0, 353.0),
        ("toluene", mol(&[6; 7], &tol), 867.0, 384.0),
        ("ethylene glycol", mol(&[8, 6, 6, 8], &chain(4)), 1110.0, 470.0),
        ("glycerol", mol(&[8, 6, 6, 6, 8, 8], &[(0, 1, Single), (1, 2, Single), (2, 3, Single), (3, 4, Single), (2, 5, Single)]), 1261.0, 563.0),
        ("urea", mol(&[7, 6, 7, 8], &[(0, 1, Single), (1, 2, Single), (1, 3, Double)]), 1320.0, 405.0),
        ("H2O2", mol(&[8, 8], &chain(2)), 1450.0, 423.0),
        ("CCl4", mol(&[6, 17, 17, 17, 17], &[(0, 1, Single), (0, 2, Single), (0, 3, Single), (0, 4, Single)]), 1594.0, 350.0),
        ("chloroform", mol(&[6, 17, 17, 17], &[(0, 1, Single), (0, 2, Single), (0, 3, Single)]), 1489.0, 334.0),
    ];
    let mut worst = (0.0f64, "");
    let mut total = 0.0;
    let n = set.len() as f64;
    for (name, a, meas, _tb) in &set {
        let p = analyse(a).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let err = (p.density / meas - 1.0).abs();
        println!("  {name:16} measured {meas:5.0} derived {:5.0} ({:+.0}%)", p.density, 100.0 * (p.density / meas - 1.0));
        total += err;
        if err > worst.0 {
            worst = (err, name);
        }
    }
    println!("  mean {:.3}, worst {:.3} ({})", total / n, worst.0, worst.1);
    assert!(total / n < 0.13, "mean error {:.3}", total / n);
    assert!(worst.0 < 0.30, "{} is out by {:.3}", worst.1, worst.0);
}
