//! Cohesion: structures that hold together, carry load, and break for reasons.

use phys::engine::{galaxy, World};
use phys::math::v3;
use phys::morph::*;
use phys::prolong::*;
use phys::solvers::structure::*;
use phys::state::*;
use phys::topology::*;
use phys::units::*;

fn tree(mass: f64) -> (Aggregate, Morphology) {
    let mut m = Morphology::new(Program::Tree, 0xACE, 0x1234, 0);
    m.built = mass;
    m.age = 40.0 * YEAR;
    let mut agg = Aggregate::neutral(mass, m.extent(), 291.0, Program::Tree.substrate());
    agg.chemical_energy = m.stored_energy();
    (agg, m)
}

fn load(agg: &Aggregate, m: &Morphology, budget: usize) -> (Vec<Body>, Topology) {
    let (b, t, _) = prolong_structured(agg, m, budget, 7, 0x1234, 0);
    (b, t)
}

/// Every part except the anchors is held by something, and the support graph
/// is acyclic — otherwise the load accumulation would never terminate.
#[test]
fn the_support_graph_is_a_well_formed_tree() {
    for program in [Program::Tree, Program::Coral, Program::Tower, Program::Wall] {
        let mut m = Morphology::new(program, 1, 2, 0);
        m.built = 5000.0;
        m.design_mass = 5000.0;
        m.progress = 1.0;
        let agg = Aggregate::neutral(5000.0, m.extent(), 290.0, program.substrate());
        let (bodies, topo, report) = prolong_structured(&agg, &m, 3000, 7, 0x2, 0);
        assert!(!topo.is_empty(), "{program:?} produced no joints");

        let n = report.structural_parts;
        let mut anchors = 0;
        for i in 0..n {
            let p = topo.support[i];
            if p == NO_SUPPORT {
                anchors += 1;
                continue;
            }
            // Parents are always emitted before their children, which is what
            // makes the one-pass load accumulation correct.
            assert!(
                (p as usize) < i,
                "{program:?}: part {i} is supported by {p}, which comes after it"
            );
        }
        assert!(anchors >= 1, "{program:?} is anchored to nothing");
        // Joints sit where the parts actually are.
        for i in 0..n {
            let d = (topo.bonds[i].at - bodies[i].pos).norm();
            let len = (topo.tip[i] - topo.base[i]).norm();
            assert!(
                d <= len * 0.51 + 1e-9,
                "{program:?}: joint {i} is {d:.3} m from its part, which is only {len:.3} m long"
            );
        }
        println!("{program:<8?} {n:>5} parts, {anchors} anchored, topology {} bytes",
            topo.bytes());
    }
}

/// A healthy tree stands up under its own weight with a real safety margin —
/// and if it did not, every failure result below would be meaningless.
#[test]
fn a_tree_stands_up() {
    let (agg, m) = tree(900.0);
    let (bodies, topo) = load(&agg, &m, 4000);
    let temps = vec![291.0; bodies.len()];
    let loads = analyse(&bodies, &topo, G_EARTH, None, &temps);
    let peak = loads.iter().fold(0.0f64, |a, l| a.max(l.utilisation));
    println!(
        "{:.1} m tree, trunk radius {:.3} m, peak self-weight utilisation {:.3} (safety factor {:.1})",
        m.tree_height(),
        topo.bonds[0].radius,
        peak,
        1.0 / peak
    );
    assert!(peak < 1.0, "the tree collapsed under its own weight");
    assert!(peak > 0.02, "safety factor {:.0} is implausibly large", 1.0 / peak);
    // The trunk carries everything.
    assert!(
        (loads[0].carried - 900.0).abs() / 900.0 < 1e-9,
        "the trunk carries {:.1} kg of a 900 kg tree",
        loads[0].carried
    );
}

/// The member proportions have to agree with the mass at the material's
/// density, or every stress in the model is wrong by the cube of the error.
#[test]
fn member_geometry_matches_its_mass() {
    let (agg, m) = tree(900.0);
    let (bodies, topo, report) = prolong_structured(&agg, &m, 4000, 7, 0x1234, 0);
    let mut volume = 0.0;
    for i in 0..report.structural_parts {
        let len = (topo.tip[i] - topo.base[i]).norm();
        volume += std::f64::consts::PI * topo.bonds[i].radius.powi(2) * len;
    }
    let expected = m.built / Program::Tree.density();
    let err = (volume - expected).abs() / expected;
    println!(
        "members enclose {volume:.3} m^3, mass implies {expected:.3} m^3 (error {:.1}%)",
        err * 100.0
    );
    assert!(err < 1e-6, "member volume disagrees with mass by {:.1}%", err * 100.0);
    let _ = bodies;
    // And a real tree's trunk is not half a metre thick at 13 m tall.
    assert!(
        (0.05..0.30).contains(&topo.bonds[0].radius),
        "trunk radius {:.3} m",
        topo.bonds[0].radius
    );
}

/// Wind: survivable gale, damaging storm, destructive hurricane.
#[test]
fn wind_damage_scales_with_speed() {
    let (agg, m) = tree(900.0);
    let mut last = 0.0;
    let mut results = Vec::new();
    for speed in [15.0, 25.0, 40.0, 60.0] {
        let (bodies, mut topo) = load(&agg, &m, 4000);
        let (extra, temps, _) = apply_insult(
            &bodies,
            &mut topo,
            Insult::Wind { speed, direction: v3(1.0, 0.0, 0.0) },
            291.0,
        );
        let loads = analyse(&bodies, &topo, G_EARTH, Some(&extra), &temps);
        let r = apply_failures(&bodies, &mut topo, &loads);
        println!(
            "  {speed:>4.0} m/s: utilisation {:.2}, {} joints broke, {:.0} kg down",
            r.peak_utilisation, r.broken_sites.len(), r.detached_mass
        );
        assert!(r.peak_utilisation > last, "damage did not increase with wind speed");
        last = r.peak_utilisation;
        results.push(r.detached_mass);
    }
    assert_eq!(results[0], 0.0, "a 15 m/s breeze should do nothing");
    assert_eq!(results[1], 0.0, "a 25 m/s gale should do nothing");
    assert!(results[2] > 0.0, "a 40 m/s storm should break something");
    assert!(results[3] >= results[2], "a hurricane should be at least as bad");
}

/// Wet snow is what brings limbs down; dry powder is not. A model that cannot
/// tell them apart is not modelling snow, it is modelling weight.
#[test]
fn only_wet_snow_breaks_branches() {
    let (agg, m) = tree(900.0);
    let survives = |depth: f64, density: f64| {
        let (bodies, mut topo) = load(&agg, &m, 4000);
        let (extra, temps, _) = apply_insult(
            &bodies,
            &mut topo,
            Insult::Snow { depth, density, crown_area: m.capture_area() },
            271.0,
        );
        let loads = analyse(&bodies, &topo, G_EARTH, Some(&extra), &temps);
        let r = apply_failures(&bodies, &mut topo, &loads);
        (r.peak_utilisation, r.detached_mass)
    };
    let (u_powder, m_powder) = survives(0.60, 100.0);
    let (u_settled, m_settled) = survives(0.30, 200.0);
    let (u_wet, m_wet) = survives(0.10, 400.0);
    println!("  600 mm powder: util {u_powder:.2}, {m_powder:.0} kg down");
    println!("  300 mm settled: util {u_settled:.2}, {m_settled:.0} kg down");
    println!("  100 mm wet:    util {u_wet:.2}, {m_wet:.0} kg down");
    assert_eq!(m_powder, 0.0, "600 mm of dry powder should not break a tree");
    assert_eq!(m_settled, 0.0, "300 mm of settled snow should not break a tree");
    assert!(m_wet > 0.0, "100 mm of wet snow should bring limbs down");
    assert!(u_wet > u_powder * 3.0, "wetness barely mattered");
}

/// Lightning follows the support chain to ground and destroys what it cannot
/// pass through — the thin members, which have the highest resistance per
/// kilogram.
#[test]
fn lightning_destroys_along_its_path() {
    let (agg, m) = tree(900.0);
    let mut previous = 0;
    for joules in [1e7, 1e8, 1e9] {
        let (bodies, mut topo) = load(&agg, &m, 4000);
        let entry = (bodies.len() / 2) as u32;
        let (extra, temps, insult) = apply_insult(
            &bodies,
            &mut topo,
            Insult::Lightning { joules, entry },
            291.0,
        );
        let loads = analyse(&bodies, &topo, G_EARTH, Some(&extra), &temps);
        let r = apply_failures(&bodies, &mut topo, &loads);
        println!(
            "  {joules:>8.0e} J: {} members destroyed on the channel, {:.1} kg down",
            insult.broken_sites.len(),
            r.detached_mass
        );
        assert!(
            insult.broken_sites.len() >= previous,
            "more energy destroyed less"
        );
        previous = insult.broken_sites.len();
        assert_eq!(insult.energy_delivered, joules);
    }
    assert!(previous > 0, "a gigajoule strike did nothing at all");
}

/// A ground fire consumes the fine fuel and scorches the trunk, because thin
/// members have small thermal mass and reach the flame temperature while thick
/// ones do not.
#[test]
fn fire_consumes_fine_fuel_first() {
    let (agg, m) = tree(900.0);
    let burn = |temperature: f64, duration: f64, height: f64| {
        let (bodies, mut topo) = load(&agg, &m, 4000);
        let (_, temps, insult) = apply_insult(
            &bodies,
            &mut topo,
            Insult::Fire { temperature, duration, height },
            291.0,
        );
        let hottest = temps.iter().cloned().fold(0.0f64, f64::max);
        let trunk = temps[0];
        (insult.consumed_mass, hottest, trunk)
    };
    let (light, _, trunk_light) = burn(700.0, 60.0, 3.0);
    let (severe, hottest, trunk_severe) = burn(1100.0, 600.0, 25.0);
    println!("  light ground fire: {light:.1} kg consumed, trunk reached {trunk_light:.0} K");
    println!("  crown fire:        {severe:.1} kg consumed, trunk {trunk_severe:.0} K, hottest {hottest:.0} K");
    assert_eq!(light, 0.0, "a brief low ground fire should not consume a mature tree");
    assert!(severe > 0.0, "a sustained crown fire should consume the tree");
    assert!(
        trunk_severe < hottest,
        "the trunk heated as fast as the twigs, which is not how thermal mass works"
    );
}

/// Damage is permanent, survives the structure being discarded, and conserves
/// what it must.
#[test]
fn damage_persists_and_conserves() {
    let mut w = World::new(galaxy(0x3333, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let node = w.tree.promote(root, 7, phys::engine::default_spec(Tier::Stellar));
    {
        let n = &mut w.tree.nodes[node.get()];
        n.agg = Aggregate::neutral(900.0, 6.0, 291.0, Program::Tree.substrate());
        n.spec.count = 3000;
    }
    w.plant(node, Program::Tree, Environment::default());
    let mass0 = w.tree.nodes[node.get()].agg.mass;
    let baryon0 = w.tree.nodes[node.get()].agg.baryon_number;
    let built0 = w.tree.nodes[node.get()].morphology.as_ref().unwrap().built;
    let entropy0 = w.tree.nodes[node.get()].agg.total_entropy();

    let out = w.damage(
        node,
        Insult::Snow { depth: 0.25, density: 450.0, crown_area: 30.0 },
    );
    println!(
        "  wet snow: {} joints, {:.0} kg down, peak utilisation {:.2}",
        out.broken_joints, out.detached_mass, out.peak_utilisation
    );
    assert!(out.broken_joints > 0, "the storm did nothing");
    assert!(out.detached_mass > 0.0);

    let n = &w.tree.nodes[node.get()];
    // Mass and baryon number are untouched: the limb is on the ground, not gone.
    assert_eq!(n.agg.mass, mass0, "node mass changed when a limb fell");
    assert!((n.agg.baryon_number - baryon0).abs() / baryon0 < 1e-12);
    assert!(n.morphology.as_ref().unwrap().built < built0, "structure did not lose mass");
    assert!(n.agg.total_entropy() >= entropy0, "total entropy fell");
    assert_eq!(w.rejected_transactions, 0);

    // And it is still deterministic, and still damaged, after regeneration.
    let a = w.tree.refine(node).to_vec();
    w.tree.coarsen(node);
    let b = w.tree.refine(node).to_vec();
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(x.pos, y.pos, "the damaged tree regenerated differently");
    }
    let events = w.tree.nodes[node.get()].morphology.as_ref().unwrap().events.len();
    assert!(events > 0, "the damage left no record");
    println!("  {events} breaks recorded and replayed through a coarsen/refine cycle");
}

/// Burning releases the free energy the structure was holding, and the atoms
/// stay put.
#[test]
fn fire_releases_stored_energy_without_losing_mass() {
    let mut w = World::new(galaxy(0x4444, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let node = w.tree.promote(root, 9, phys::engine::default_spec(Tier::Stellar));
    {
        let n = &mut w.tree.nodes[node.get()];
        n.agg = Aggregate::neutral(900.0, 6.0, 291.0, Program::Tree.substrate());
        n.spec.count = 3000;
    }
    w.plant(node, Program::Tree, Environment::default());
    let mass0 = w.tree.nodes[node.get()].agg.mass;
    let chem0 = w.tree.nodes[node.get()].agg.chemical_energy;
    let internal0 = w.tree.nodes[node.get()].agg.internal_energy;

    let out = w.damage(
        node,
        Insult::Fire { temperature: 1100.0, duration: 600.0, height: 30.0 },
    );
    let n = &w.tree.nodes[node.get()];
    println!(
        "  {:.0} kg consumed, {:.3e} J of chemical energy released as heat",
        out.consumed_mass, out.energy_released
    );
    assert!(out.consumed_mass > 0.0, "nothing burned");
    assert_eq!(n.agg.mass, mass0, "combustion lost mass; the atoms have to go somewhere");
    assert!(n.agg.chemical_energy < chem0, "burning released no stored energy");
    assert!(n.agg.internal_energy > internal0, "the fire produced no heat");
    assert_eq!(w.rejected_transactions, 0, "a combustion transaction failed to balance");
}

/// Topology costs little enough to be worth having on every structure.
#[test]
fn topology_is_cheap() {
    let (agg, m) = tree(900.0);
    let (bodies, topo, _) = prolong_structured(&agg, &m, 8000, 7, 0x1234, 0);
    let geometry = bodies.len() * std::mem::size_of::<Body>();
    let cohesion = topo.bytes();
    println!(
        "  {} parts: {:.2} MB of geometry, {:.2} MB of joints ({:.0}% overhead)",
        bodies.len(),
        geometry as f64 / 1e6,
        cohesion as f64 / 1e6,
        100.0 * cohesion as f64 / geometry as f64
    );
    assert!(
        cohesion < geometry,
        "the joints cost more than the parts they join"
    );
}

/// The renderer produces a real image of the real geometry.
#[test]
fn renderer_draws_the_structure() {
    use phys::render::*;
    let (agg, m) = tree(900.0);
    let (bodies, topo) = load(&agg, &m, 3000);
    let cam = Camera::framing(v3(0.0, 0.0, 0.0), m.extent() * 1.2, 0.6, 0.1);
    let intact = vec![true; bodies.len()];
    let mut canvas = Canvas::new(320, 260);
    draw_structure(&mut canvas, &cam, &bodies, &topo, &intact, &Style::daylight());

    // Something was actually drawn: the tree's colours are not the sky's.
    let sky = Style::daylight().sky_top;
    let drawn = canvas
        .rgb
        .chunks(3)
        .filter(|p| {
            let d = (p[0] as i32 - sky[0] as i32).abs()
                + (p[1] as i32 - sky[1] as i32).abs()
                + (p[2] as i32 - sky[2] as i32).abs();
            d > 90
        })
        .count();
    let coverage = drawn as f64 / (320.0 * 260.0);
    println!("  tree covers {:.1}% of the frame", coverage * 100.0);
    assert!(coverage > 0.02, "the renderer drew almost nothing");
    assert!(coverage < 0.9, "the renderer filled the frame");

    let path = std::env::temp_dir().join("phys_render_test.png");
    write_png(&canvas, path.to_str().unwrap()).expect("png");
    let bytes = std::fs::read(&path).expect("read back");
    // A valid PNG: signature, then IHDR, and ending in IEND.
    assert_eq!(&bytes[..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
    assert_eq!(&bytes[12..16], b"IHDR");
    assert_eq!(&bytes[bytes.len() - 8..bytes.len() - 4], b"IEND");
    println!("  wrote a valid {} byte PNG", bytes.len());
    let _ = std::fs::remove_file(&path);
}
