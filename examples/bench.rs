//! Measured cost of every hot path, on whatever machine runs it.
//! These are the numbers `docs/PERFORMANCE.md` reasons from.

use phys::prolong::*;
use phys::solvers::*;
use phys::state::*;
use phys::units::*;
use std::time::Instant;

fn bodies(n: usize, profile: Profile, kind: BodyKind, agg: &Aggregate) -> Vec<Body> {
    let spec = ProlongSpec::new(n, profile, MassSpectrum::Equal, kind);
    prolong(agg, spec, 1, 0x1234, 0).0
}

fn time<F: FnMut()>(reps: usize, mut f: F) -> f64 {
    let t = Instant::now();
    for _ in 0..reps {
        f();
    }
    t.elapsed().as_secs_f64() * 1e6 / reps as f64
}

fn main() {
    println!("# measured on this machine, single core, release build\n");

    println!("## materialisation (prolong)");
    let agg = Aggregate::neutral(1e5 * M_SUN, PARSEC, 30.0, Composition::solar());
    for n in [1_000usize, 10_000, 100_000, 500_000] {
        let spec = ProlongSpec::new(n, Profile::Plummer, MassSpectrum::Equal, BodyKind::Star);
        let us = time(2, || {
            std::hint::black_box(prolong(&agg, spec, 1, 0x99, 0));
        });
        println!("  n={n:>9}  {us:>10.0} us   {:>7.3} us/body   {:>8.1} M bodies/s",
            us / n as f64, n as f64 / us);
    }

    println!("\n## gravity (Barnes-Hut, theta=0.5, one leapfrog step)");
    for n in [1_000usize, 10_000, 50_000] {
        let b0 = bodies(n, Profile::Plummer, BodyKind::Star, &agg);
        let p = gravity::GravityParams { theta: 0.5, softening: PARSEC * 0.01, retarded: false, post_newtonian: false, quadrupole: false };
        let mut b = b0.clone();
        let us = time(2, || {
            gravity::step_leapfrog(&mut b, 1e3 * YEAR, p);
        });
        println!("  n={n:>9}  {us:>10.0} us   {:>7.3} us/body   {:>8.2} M bodies/s",
            us / n as f64, n as f64 / us);
    }

    println!("\n## gravity: cost of the extras");
    let b0 = bodies(50_000, Profile::Plummer, BodyKind::Star, &agg);
    for (name, p) in [
        ("monopole only, theta=0.7", gravity::GravityParams { theta: 0.7, softening: PARSEC * 0.01, retarded: false, post_newtonian: false, quadrupole: false }),
        ("theta=0.5", gravity::GravityParams { theta: 0.5, softening: PARSEC * 0.01, retarded: false, post_newtonian: false, quadrupole: false }),
        ("theta=0.5 + retardation", gravity::GravityParams { theta: 0.5, softening: PARSEC * 0.01, retarded: true, post_newtonian: false, quadrupole: false }),
        ("theta=0.5 + retard + 1PN", gravity::GravityParams { theta: 0.5, softening: PARSEC * 0.01, retarded: true, post_newtonian: true, quadrupole: false }),
        ("theta=0.5 + quadrupole", gravity::GravityParams { theta: 0.5, softening: PARSEC * 0.01, retarded: false, post_newtonian: false, quadrupole: true }),
        ("theta=0.3", gravity::GravityParams { theta: 0.3, softening: PARSEC * 0.01, retarded: false, post_newtonian: false, quadrupole: false }),
    ] {
        let mut b = b0.clone();
        let us = time(2, || {
            gravity::step_leapfrog(&mut b, 1e3 * YEAR, p);
        });
        println!("  {name:<28} {us:>9.0} us   {:>6.3} us/body", us / 50_000.0);
    }

    println!("\n## hydrodynamics (SPH, ~50 neighbours)");
    let gas = Aggregate::neutral(1e30, 1e12, 1e4, Composition::solar());
    for n in [1_000usize, 10_000, 50_000] {
        let mut b = bodies(n, Profile::Uniform, BodyKind::GasParcel, &gas);
        let h = 1e12 / (n as f64).cbrt() * 1.2;
        let p = hydro::HydroParams { h, ..Default::default() };
        let us = time(2, || {
            hydro::step(&mut b, 1.0, p);
        });
        println!("  n={n:>9}  {us:>10.0} us   {:>7.3} us/body   {:>8.2} M bodies/s",
            us / n as f64, n as f64 / us);
    }

    println!("\n## molecular dynamics (LJ, cell lists)");
    for n in [1_000usize, 10_000, 100_000] {
        let side = 4e-9 * (n as f64 / 4096.0).cbrt();
        let mol = Aggregate::neutral(n as f64 * 12.0 * AMU, side, 300.0, Composition::pure(Species::Carbon));
        let mut b = bodies(n, Profile::Uniform, BodyKind::Atom, &mol);
        let p = md::MdParams::default();
        let us = time(2, || {
            md::step(&mut b, 1e-15, p, 1, 1, 0, 0);
        });
        println!("  n={n:>9}  {us:>10.0} us   {:>7.3} us/body   {:>8.2} M bodies/s",
            us / n as f64, n as f64 / us);
    }

    println!("\n## restriction (coarsen)");
    for n in [10_000usize, 100_000, 500_000] {
        let b = bodies(n, Profile::Plummer, BodyKind::Star, &agg);
        let us = time(5, || {
            std::hint::black_box(restrict(&b, 0.0));
        });
        println!("  n={n:>9}  {us:>10.0} us   {:>7.3} us/body", us / n as f64);
    }

    println!("\n## structural analysis and dynamics");
    {
        use phys::morph::{Morphology, Program};
        use phys::prolong::prolong_structured;
        use phys::solvers::structure as st;
        use phys::math::v3;

        for (label, planned) in [("tree (determinate)", false), ("tower (braced)", true)] {
            let mass = if planned { 3.0e6 } else { 900.0 };
            let mut m = if planned {
                let mut m = Morphology::planned(Program::Tower, mass, 11, 0x77);
                m.progress = 1.0;
                m
            } else {
                let mut m = Morphology::new(Program::Tree, 0xACE, 0x1234, 0);
                m.age = 40.0 * YEAR;
                m
            };
            m.built = mass;
            let prog = if planned { Program::Tower } else { Program::Tree };
            let agg = Aggregate::neutral(mass, m.extent(), 291.0, prog.substrate());
            for n in [500usize, 2000, 8000] {
                let (b, topo, _) = prolong_structured(&agg, &m, n, 7, 0x1234, 0);
                let mut field = st::LoadField::new(b.len(), 291.0);
                field.apply(&st::weather::wind(25.0, v3(1.0, 0.0, 0.0)), &b, &topo);
                field.apply(&st::weather::gravity(), &b, &topo);
                let members = b.len();

                let us = time(3, || {
                    std::hint::black_box(st::analyse_with(&b, &topo, &field));
                });
                print!("  {label:<20} n={members:>6}  static {us:>9.0} us");

                match st::dynamic_structure(&b, &topo) {
                    Some(mut ds) => {
                        ds.advance(&field, 0.01);
                        let du = time(5, || {
                            ds.advance(&field, 0.01);
                        });
                        println!("   dynamic {du:>9.0} us/substep   {:.0} substeps in a 50 ms frame",
                            50_000.0 / du.max(1e-9));
                    }
                    None => println!("   dynamic     n/a"),
                }
            }
        }
    }

    println!("\n## memory");
    println!("  Body           {:>4} bytes", std::mem::size_of::<Body>());
    println!("  Aggregate      {:>4} bytes", std::mem::size_of::<Aggregate>());
    println!("  Node           {:>4} bytes", std::mem::size_of::<phys::tree::Node>());
    println!("  Snapshot       {:>4} bytes", std::mem::size_of::<phys::causal::Snapshot>());
    let per_gb = 1e9 / std::mem::size_of::<Body>() as f64;
    println!("  bodies per GB  {:>4.1} M", per_gb / 1e6);
    println!("  6 GB card, 60% for bodies: {:.1} M bodies resident", 6.0 * 0.6 * per_gb / 1e6);
}
