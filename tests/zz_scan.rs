use phys::liquid::{extra_sites_from_text, pair_energy, with_bisector_site, with_frame_sites, PairEnergy, Polarisable, SiteSite};
use phys::math::{Quat, Vec3};
use phys::rng::{Purpose, Stream};

fn rotation(s: &mut Stream) -> Quat {
    let (u1, u2, u3) = (s.uniform(), s.uniform(), s.uniform());
    let tau = std::f64::consts::TAU;
    let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
    Quat { w: a * (tau * u2).sin(), v: Vec3 { x: a * (tau * u2).cos(), y: b * (tau * u3).sin(), z: b * (tau * u3).cos() } }
}

struct Law {
    law: SiteSite,
    alpha: Vec<f64>,
    bis: (usize, f64),
    frame: Vec<phys::liquid::ExtraSite>,
}

fn load(file: &str) -> Law {
    let text = std::fs::read_to_string(file).unwrap();
    let law = SiteSite::from_text(&text).unwrap();
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    Law { bis: SiteSite::bisector_from_text(&text).unwrap(), frame: extra_sites_from_text(&text), law, alpha }
}

fn energy(l: &Law, a: &[(Vec3, usize)], b: &[(Vec3, usize)]) -> f64 {
    let place = |m: &[(Vec3, usize)]| with_frame_sites(&with_bisector_site(m, l.bis.0, l.bis.1), &l.frame);
    let pe = PairEnergy { a: place(a), b: place(b), energy: 0.0 };
    pair_energy(&Polarisable { law: &l.law, alpha: l.alpha.clone() }, &pe) * 627.509474
}

#[test]
fn scan() {
    let names: Vec<String> = std::env::var("SCAN_LAWS").unwrap().split(',').map(|x| x.to_string()).collect();
    let laws: Vec<(String, Law)> = names.iter().map(|n| (n.clone(), load(n))).collect();
    let r = 1.8499284924411517;
    let th = 104.078f64.to_radians();
    let body = [Vec3::ZERO, Vec3 { x: r, y: 0.0, z: 0.0 }, Vec3 { x: r * th.cos(), y: r * th.sin(), z: 0.0 }];
    let mut s = Stream::at(99, 0, 0, Purpose::Positions);
    // random dimers at centre separations 3.0 to 6.5 bohr, atoms never nearer than 1.2 bohr
    let mut worst: Vec<(f64, f64, f64, f64)> = vec![(0.0, 0.0, 0.0, 0.0); laws.len()];
    let mut count = vec![0usize; laws.len()];
    let mut bins = vec![[f64::INFINITY; 8]; laws.len()];
    for _ in 0..400000 {
        let qa = rotation(&mut s);
        let qb = rotation(&mut s);
        let dir = s.direction();
        let sep = s.range(3.0, 6.5);
        let a: Vec<(Vec3, usize)> = body.iter().enumerate().map(|(i, p)| (qa.rotate(*p), if i == 0 { 0 } else { 1 })).collect();
        let b: Vec<(Vec3, usize)> = body.iter().enumerate().map(|(i, p)| (qb.rotate(*p) + dir.scale(sep), if i == 0 { 0 } else { 1 })).collect();
        let near = a.iter().flat_map(|x| b.iter().map(move |y| (x.0 - y.0).norm())).fold(f64::INFINITY, f64::min);
        if near < 1.2 {
            continue;
        }
        for (k, (_, l)) in laws.iter().enumerate() {
            let e = energy(l, &a, &b);
            count[k] += 1;
            let bi = (((near * 0.529177) - 0.6) / 0.2).floor().clamp(0.0, 7.0) as usize;
            if e < bins[k][bi] { bins[k][bi] = e; }
            if e < worst[k].0 {
                worst[k] = (e, near * 0.529177, (a[0].0 - b[0].0).norm() * 0.529177, sep * 0.529177);
            }
        }
    }
    for (k, (name, _)) in laws.iter().enumerate() {
        println!("{name}: lowest energy by nearest-atom distance (0.6-0.8, 0.8-1.0, ... 2.0+ A): {:?}", bins[k].iter().map(|x| (x * 10.0).round() / 10.0).collect::<Vec<_>>());
        println!("{name}: most negative dimer energy {:.2} kcal/mol, nearest atoms {:.2} A, O-O {:.2} A, centres {:.2} A ({} configurations)", worst[k].0, worst[k].1, worst[k].2, worst[k].3, count[k]);
    }
}
