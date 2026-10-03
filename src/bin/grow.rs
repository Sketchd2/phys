//! Grow a basis for one molecule and relax it, reporting every round —
//! `PLAY.md` E5b. A run on anything larger than water takes hours, so it
//! writes where it has got to after every round and carries on from there
//! when started again.
//!
//! ```sh
//! cargo run --release --bin phys-grow -- water
//! cargo run --release --bin phys-grow -- "6,6,8;0-1,1-2" ethanol
//! ```
//!
//! A molecule is its heavy atoms by atomic number and the bonds between them
//! (`-` single, `=` double, `#` triple); hydrogens fill each heavy atom's
//! valence. The second argument names the run and its state file,
//! `grow-<name>.state` in the working directory.

use phys::chem::arrange::{Arrangement, Bond, Order};
use phys::chem::elements::Element;
use phys::electrons::functional::Functional;
use phys::electrons::grow::{grow_from, Coordinate, Resume};
use phys::electrons::molecule::Molecule;
use std::io::Write;
use std::time::Instant;

fn arrangement(heavy: &[u8], bonds: &[(usize, usize, Order)]) -> Arrangement {
    let mut atoms: Vec<Element> = heavy.iter().map(|z| Element(*z)).collect();
    let mut bs: Vec<Bond> = bonds.iter().map(|(a, b, o)| Bond::new(*a, *b, *o)).collect();
    for i in 0..heavy.len() {
        let used: usize = bonds.iter().filter(|(a, b, _)| *a == i || *b == i).map(|(_, _, o)| o.slots()).sum();
        let v = atoms[i].valence().expect("an element with a valence");
        for _ in used..v {
            atoms.push(Element(1));
            bs.push(Bond::new(i, atoms.len() - 1, Order::Single));
        }
    }
    Arrangement::molecule(atoms, bs)
}

fn parse(spec: &str) -> Option<(Vec<u8>, Vec<(usize, usize, Order)>)> {
    let (atoms, bonds) = spec.split_once(';').unwrap_or((spec, ""));
    let heavy = atoms.split(',').map(|z| z.trim().parse().ok()).collect::<Option<Vec<u8>>>()?;
    let mut out = Vec::new();
    for b in bonds.split(',').map(str::trim).filter(|b| !b.is_empty()) {
        let (sep, order) = [('-', Order::Single), ('=', Order::Double), ('#', Order::Triple)].into_iter().find(|(c, _)| b.contains(*c))?;
        let (x, y) = b.split_once(sep)?;
        out.push((x.trim().parse().ok()?, y.trim().parse().ok()?, order));
    }
    Some((heavy, out))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let spec = args.first().map(String::as_str).unwrap_or("water");
    let spec_text = match spec {
        "water" => "8",
        "methane" => "6",
        "ammonia" => "7",
        "methanol" => "6,8;0-1",
        "ethanol" => "6,6,8;0-1,1-2",
        other => other,
    };
    let name = args.get(1).cloned().unwrap_or_else(|| spec.replace([',', ';', '-', '=', '#'], "_"));
    let (heavy, bonds) = parse(spec_text).unwrap_or_else(|| panic!("cannot read molecule {spec_text:?}"));
    let arr = arrangement(&heavy, &bonds);
    let mol = Molecule::from_arrangement(&arr);
    let all_bonds: Vec<(usize, usize)> = arr.bonds.iter().map(|b| (b.a as usize, b.b as usize)).collect();
    let state = format!("grow-{name}.state");
    let resume = std::fs::read_to_string(&state).ok().and_then(|t| Resume::from_text(&t));
    let a = 1.0 / 0.529177210903;
    let mut round = resume.as_ref().map(|r| r.history.len()).unwrap_or(0);
    println!("{name}: {} atoms, {} bonds; {}", mol.z.len(), all_bonds.len(), if resume.is_some() { format!("carrying on from round {round}") } else { "starting".into() });
    let t = Instant::now();
    let g = grow_from(&mol, &all_bonds, Functional::Pbe, 40, 3, resume, &mut |snap| {
        let r = snap.round;
        let fmt = |c: &Coordinate, v: f64, sign: bool| match c {
            Coordinate::Bond(..) => if sign { format!("{:+.5}", v / a) } else { format!("{:.5}", v / a) },
            _ => if sign { format!("{:+.3}", v.to_degrees()) } else { format!("{:.3}", v.to_degrees()) },
        };
        let vals: Vec<String> = snap.coordinates.iter().zip(&r.values).map(|(c, v)| fmt(c, *v, false)).collect();
        let pred: Vec<String> = snap.coordinates.iter().zip(&r.predicted).map(|(c, v)| fmt(c, *v, true)).collect();
        println!("round {round}: {} functions, E {:.8}, {:.0} s", r.functions, r.energy, t.elapsed().as_secs_f64());
        println!("  shape (A, deg): {}", vals.join(" "));
        println!("  added {}, predicted to move: {}; left {:.2}", r.added.len(), pred.join(" "), r.left);
        std::io::stdout().flush().ok();
        let resume = Resume { ladders: snap.ladders.to_vec(), chosen: snap.chosen.to_vec(), positions: snap.positions.to_vec(), history: snap.history.to_vec() };
        let tmp = format!("{state}.tmp");
        if std::fs::write(&tmp, resume.to_text()).and_then(|_| std::fs::rename(&tmp, &state)).is_err() {
            eprintln!("could not write {state}");
        }
        round += 1;
    });
    println!("{name}: {} after {} rounds, {:.0} s", if g.converged { "settled" } else { "not settled" }, round, t.elapsed().as_secs_f64());
}
