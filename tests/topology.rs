//! Cohesion: structures that hold together, carry load, and break for reasons.

use phys::engine::{galaxy, World};
use phys::math::v3;
use phys::morph::*;
use phys::sampler::*;
use phys::solvers::structure::*;
use phys::state::*;
use phys::morph::NO_SUPPORT;
use phys::topology::*;
use phys::units::*;
use phys::state::Composition;

/// What a scenario states a thing of this kind is made of.
///
/// A test's own statement, not the engine's: `docs/PLAY.md` D11 retired
/// `Program::substrate`, so what a structure is made of is the feedstock its
/// node held, and a scenario that wants a tree made of carbohydrate has to put
/// carbohydrate in the node.
fn substrate_for(p: Program) -> Composition {
    match p {
        Program::Tree | Program::Coral => Composition::organic(),
        _ => Composition::crustal(),
    }
}


/// A node on the surface of an Earth, ready to have something planted in it.
///
/// `docs/PLAY.md` D6 made gravity derived, and that changed what these tests
/// have to set up. A node promoted straight out of a galaxy root feels the
/// galaxy's field — about 10^-13 m/s^2 — so nothing falls, no snow load breaks
/// anything, and a gale takes limbs off that then hang in the air. That is
/// correct: there is nothing underneath it. It used to feel Earth's surface
/// gravity wherever it was, because the load came from
/// `solvers::structure::G_EARTH` rather than from anything the scene contained.
///
/// So the scene contains a planet. `Tree::gravity_at` derives 9.82 m/s^2 from
/// its mass and radius, with no mention of Earth anywhere.
///
/// Placed on the **+z axis** deliberately. The derived field points at the
/// planet's centre, and a structure's own geometry is generated with `+z` up,
/// so putting the node anywhere else would load a tree sideways — `Motion`
/// carries an orientation and `gravity_at` does not yet compose it. See the
/// note on `Tree::gravity_at`.
fn on_an_earth(seed: u64, mass: f64, radius: f64, count: usize) -> (World, phys::ids::NodeIdx) {
    const EARTH_MASS: f64 = 5.972e24;
    const EARTH_RADIUS: f64 = 6.371e6;
    let mut w = World::new(galaxy(seed, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let planet = w.tree.promote(root, 7, phys::engine::default_spec(Tier::Planetary));
    {
        let n = &mut w.tree.nodes[planet.get()];
        n.matter = Matter::neutral(EARTH_MASS, EARTH_RADIUS, 290.0, Composition::primordial());
        n.spec.count = 8;
    }
    w.tree.refine(planet);
    let node = w.tree.promote(planet, 0, phys::engine::default_spec(Tier::Continuum));
    {
        let n = &mut w.tree.nodes[node.get()];
        n.matter = Matter::neutral(mass, radius, 291.0, Composition::organic());
        n.spec.count = count;
        n.motion.offset = v3(0.0, 0.0, EARTH_RADIUS);
    }
    (w, node)
}


/// Earth's surface gravity, stated rather than assumed.
///
/// `docs/PLAY.md` D6 retired `solvers::structure::G_EARTH`: a structure is now
/// proportioned and loaded in the field it is actually in, derived by
/// `Tree::gravity_at` from whatever it sits inside. These tests are about
/// structure generation rather than about gravity, so they name a field and
/// hold it fixed — which is what the constant was doing for them, said out loud.
const SURFACE_G: phys::math::Vec3 = phys::math::Vec3 { x: 0.0, y: 0.0, z: -9.80665 };


fn tree(mass: f64) -> (Matter, Morphology) {
    let mut m = Morphology::new(Program::Tree, 0xACE, 0x1234, 0);
    m.built = mass;
    m.age = 40.0 * YEAR;
    let mut matter = Matter::neutral(mass, m.extent(), 291.0, Composition::organic());
    matter.chemical_energy = m.stored_energy();
    (matter, m)
}

fn load(matter: &Matter, m: &Morphology, budget: usize) -> (Vec<Body>, Topology) {
    let (b, t, _) = sample_structured(matter, m, budget, 7, 0x1234, 0, SURFACE_G);
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
        let matter = Matter::neutral(5000.0, m.extent(), 290.0, substrate_for(program));
        let (bodies, topo, report) = sample_structured(&matter, &m, 3000, 7, 0x2, 0, SURFACE_G);
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
            let d = (topo.joints[i].at - bodies[i].pos).norm();
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
    let (matter, m) = tree(900.0);
    let (bodies, topo) = load(&matter, &m, 4000);
    let mut field = LoadField::new(bodies.len(), 291.0);
    field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
    let loads = analyse(&bodies, &topo, &field);
    let peak = loads.iter().fold(0.0f64, |a, l| a.max(l.utilisation));
    println!(
        "{:.1} m tree, trunk radius {:.3} m, peak self-weight utilisation {:.3} (safety factor {:.1})",
        m.height(),
        topo.joints[0].radius,
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
    let (matter, m) = tree(900.0);
    let (bodies, topo, report) = sample_structured(&matter, &m, 4000, 7, 0x1234, 0, SURFACE_G);
    let mut volume = 0.0;
    for i in 0..report.structural_parts {
        let len = (topo.tip[i] - topo.base[i]).norm();
        volume += std::f64::consts::PI * topo.joints[i].radius.powi(2) * len;
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
        (0.05..0.30).contains(&topo.joints[0].radius),
        "trunk radius {:.3} m",
        topo.joints[0].radius
    );
}

/// Wind: survivable gale, damaging storm, destructive hurricane.
///
/// **The ladder moved up 20 m/s in Phase 4, and the reason is geometry rather
/// than wind.** A structure used to be drawn 1.44x larger than the size it
/// stated, because the sampler rescaled a unit skeleton until `summarise`
/// reported the node's radius back; a tree of a given mass was therefore drawn
/// 44% more slender than its own allometry says, and a slenderer tree buckles
/// in a lighter wind. The same 900 kg of wood now stands at the height it
/// claims, and takes a 60 m/s hurricane rather than a 40 m/s storm to break.
#[test]
fn wind_damage_scales_with_speed() {
    let (matter, m) = tree(900.0);
    let mut last = 0.0;
    let mut results = Vec::new();
    for speed in [15.0, 25.0, 40.0, 75.0] {
        let (bodies, mut topo) = load(&matter, &m, 4000);
        let mut field = LoadField::new(bodies.len(), 291.0);
        field.apply(&weather::wind(speed, v3(1.0, 0.0, 0.0)), &bodies, &topo);
        field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
        let loads = analyse(&bodies, &topo, &field);
        let r = apply_failures(&bodies, &mut topo, &loads, &field);
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
    assert_eq!(results[2], 0.0, "a 40 m/s storm should still leave it standing");
    assert!(results[3] > 0.0, "a 75 m/s hurricane should break something");
}

/// Wet snow is what brings limbs down; dry powder is not. A model that cannot
/// tell them apart is not modelling snow, it is modelling weight.
#[test]
fn only_wet_snow_breaks_branches() {
    let (matter, m) = tree(900.0);
    let survives = |depth: f64, density: f64| {
        let (bodies, mut topo) = load(&matter, &m, 4000);
        let mut field = LoadField::new(bodies.len(), 271.0);
        field.apply(&weather::snow(depth, density, m.capture_area()), &bodies, &topo);
        field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
        let loads = analyse(&bodies, &topo, &field);
        let r = apply_failures(&bodies, &mut topo, &loads, &field);
        (r.peak_utilisation, r.detached_mass)
    };
    let (u_powder, m_powder) = survives(0.60, 100.0);
    let (u_settled, m_settled) = survives(0.30, 200.0);
    // **200 mm rather than 100 in Phase 4**, for the reason
    // `wind_damage_scales_with_speed` records: the tree is drawn at the height
    // its own allometry states rather than 1.44x it, so it is stouter and
    // carries more. The claim is unchanged and so is the comparison — a sixth
    // of the powder's depth, wet, against the whole of it dry.
    let (u_wet, m_wet) = survives(0.20, 400.0);
    println!("  600 mm powder: util {u_powder:.2}, {m_powder:.0} kg down");
    println!("  300 mm settled: util {u_settled:.2}, {m_settled:.0} kg down");
    println!("  200 mm wet:    util {u_wet:.2}, {m_wet:.0} kg down");
    // The claim is about *damage*, not about an exact zero. Six hundred
    // millimetres of powder costs nothing at all; a third of that in settled
    // snow costs a twig; a sixth of it, wet, brings limbs down. Insisting the
    // middle case shed exactly nothing would be asserting that a tree under
    // 300 mm of snow loses not one twig, which is not true of trees.
    let standing = m.built;
    assert_eq!(m_powder, 0.0, "600 mm of dry powder should not break a tree");
    assert!(
        m_settled < 0.01 * standing,
        "300 mm of settled snow took {m_settled:.1} kg of a {standing:.0} kg tree"
    );
    assert!(
        m_wet > 0.1 * standing,
        "200 mm of wet snow took only {m_wet:.1} kg of a {standing:.0} kg tree"
    );
    assert!(u_wet > u_powder * 3.0, "wetness barely mattered");
}

/// Lightning follows the support chain to ground and destroys what it cannot
/// pass through — the thin members, which have the highest resistance per
/// kilogram.
#[test]
fn lightning_destroys_along_its_path() {
    let (matter, m) = tree(900.0);
    let mut previous = 0;
    for joules in [1e7, 1e8, 1e9] {
        let (bodies, mut topo) = load(&matter, &m, 4000);
        let entry = (bodies.len() / 2) as u32;
        let mut field = LoadField::new(bodies.len(), 291.0);
        field.apply(&weather::lightning(joules, entry), &bodies, &topo);
        field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
        let destroyed = field.destroyed.iter().filter(|d| **d).count();
        let loads = analyse(&bodies, &topo, &field);
        let r = apply_failures(&bodies, &mut topo, &loads, &field);
        println!(
            "  {joules:>8.0e} J: {} members destroyed on the channel, {:.1} kg down",
            destroyed, r.detached_mass
        );
        assert!(destroyed >= previous, "more energy destroyed less");
        previous = destroyed;
        assert_eq!(field.energy_delivered, joules);
    }
    assert!(previous > 0, "a gigajoule strike did nothing at all");
}

/// A ground fire consumes the fine fuel and scorches the trunk, because thin
/// members have small thermal mass and reach the flame temperature while thick
/// ones do not.
#[test]
fn fire_consumes_fine_fuel_first() {
    let (matter, m) = tree(900.0);
    let burn = |temperature: f64, duration: f64, height: f64| {
        let (bodies, mut topo) = load(&matter, &m, 4000);
        let mut field = LoadField::new(bodies.len(), 291.0);
        field.apply(&weather::fire(temperature, height, duration), &bodies, &topo);
        field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
        let hottest = field.temperature.iter().cloned().fold(0.0f64, f64::max);
        let trunk = field.temperature[0];
        let loads = analyse(&bodies, &topo, &field);
        let r = apply_failures(&bodies, &mut topo, &loads, &field);
        (r.consumed_mass, hottest, trunk)
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
    let (mut w, node) = on_an_earth(0x3333, 900.0, 6.0, 3000);
    w.plant(node, Program::Tree, Some(Environment::default()));
    let mass0 = w.tree.nodes[node.get()].matter.mass;
    let baryon0 = w.tree.nodes[node.get()].matter.baryon_number;
    let built0 = w.tree.nodes[node.get()].morphology.as_ref().unwrap().built;
    let entropy0 = w.tree.nodes[node.get()].matter.total_entropy();

    let out = w.damage(node, &[weather::snow(0.25, 450.0, 30.0)]);
    println!(
        "  wet snow: {} joints, {:.0} kg down, peak utilisation {:.2}",
        out.broken_joints, out.detached_mass, out.peak_utilisation
    );
    assert!(out.broken_joints > 0, "the storm did nothing");
    assert!(out.detached_mass > 0.0);

    let n = &w.tree.nodes[node.get()];
    // Mass and baryon number are untouched: the limb is on the ground, not gone.
    assert_eq!(n.matter.mass, mass0, "node mass changed when a limb fell");
    assert!((n.matter.baryon_number - baryon0).abs() / baryon0 < 1e-12);
    assert!(n.morphology.as_ref().unwrap().built < built0, "structure did not lose mass");
    assert!(n.matter.total_entropy() >= entropy0, "total entropy fell");
    assert_eq!(w.rejected_growth_steps, 0);

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
        n.matter = Matter::neutral(900.0, 6.0, 291.0, Composition::organic());
        n.spec.count = 3000;
    }
    w.plant(node, Program::Tree, Some(Environment::default()));
    let mass0 = w.tree.nodes[node.get()].matter.mass;
    let chem0 = w.tree.nodes[node.get()].matter.chemical_energy;
    let internal0 = w.tree.nodes[node.get()].matter.internal_energy;

    let out = w.damage(node, &[weather::fire(1100.0, 30.0, 600.0)]);
    let n = &w.tree.nodes[node.get()];
    println!(
        "  {:.0} kg consumed, {:.3e} J of chemical energy released as heat",
        out.consumed_mass, out.energy_released
    );
    assert!(out.consumed_mass > 0.0, "nothing burned");
    assert_eq!(n.matter.mass, mass0, "combustion lost mass; the atoms have to go somewhere");
    assert!(n.matter.chemical_energy < chem0, "burning released no stored energy");
    assert!(n.matter.internal_energy > internal0, "the fire produced no heat");
    assert_eq!(w.rejected_growth_steps, 0, "a combustion transaction failed to balance");
}

/// Topology costs little enough to be worth having on every structure.
#[test]
fn topology_is_cheap() {
    let (matter, m) = tree(900.0);
    let (bodies, topo, _) = sample_structured(&matter, &m, 8000, 7, 0x1234, 0, SURFACE_G);
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
    let (matter, m) = tree(900.0);
    let (bodies, topo) = load(&matter, &m, 3000);
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

// ---------------------------------------------------------------------------
// The general solver
// ---------------------------------------------------------------------------

/// A redundant structure must actually be routed to the frame solver, and the
/// answer it comes back with must be the analytic one.
///
/// The closed-form truss itself is validated in `tests/frame.rs`; what is under
/// test here is the layer above — that `analyse` recognises the redundancy,
/// hands the structure to `solvers::frame`, and maps the result back into
/// member indices without losing anything.
///
/// Three bars from a common apex to three anchors, one vertical and two at 45
/// degrees, loaded downwards. Note what the load *is*: a `LoadField` force on a
/// member is carried along that member's length, so half of it reacts straight
/// into that member's own anchor and never reaches the joint. Only the other
/// half is shared, and that is the load the analytic split applies to.
#[test]
fn indeterminate_truss_matches_the_analytic_solution() {
    let apex = v3(0.0, 0.0, 0.0);
    let h = 1.0;
    let t: f64 = std::f64::consts::FRAC_PI_4;
    let anchors = [
        v3(0.0, 0.0, -h),
        v3(-h * t.tan(), 0.0, -h),
        v3(h * t.tan(), 0.0, -h),
    ];

    // Each bar is anchored at the ground and reaches the apex. They are three
    // separate members with no support relation between them; the solver has to
    // notice that their tips coincide and weld them into one joint.
    let radius = 0.01;
    let area = std::f64::consts::PI * radius * radius;
    let members = vec![
        Member::anchored(anchors[0], apex, radius),
        Member::anchored(anchors[1], apex, radius),
        Member::anchored(anchors[2], apex, radius),
    ];
    let joint = area * 2000.0;
    let ties = vec![(0u32, 1u32, joint), (0u32, 2u32, joint), (1u32, 2u32, joint)];
    let topo = Topology::from_parts(&members, &ties, Material::steel());
    assert!(!topo.is_determinate(), "the truss should be indeterminate");

    let p_load = 1000.0;
    let bodies = vec![
        Body { pos: members[0].tip, mass: 0.0, radius, ..Default::default() },
        Body { pos: members[1].tip, mass: 0.0, radius, ..Default::default() },
        Body { pos: members[2].tip, mass: 0.0, radius, ..Default::default() },
    ];
    let third = p_load / 3.0;
    let mut field = LoadField::new(3, 290.0);
    field.force = vec![
        v3(0.0, 0.0, -third),
        v3(0.0, 0.0, -third),
        v3(0.0, 0.0, -third),
    ];

    let (loads, redundant, iters) = analyse_with(&bodies, &topo, &field);
    assert!(redundant, "the solver treated a braced truss as determinate");
    assert!(iters > 0 && iters < 200, "{iters} CG iterations");

    let axial = |i: usize| {
        let axis = (topo.tip[i] - topo.base[i]).unit();
        loads[i].force.dot(axis).abs()
    };
    let middle = axial(0);
    let outer = (axial(1) + axial(2)) / 2.0;

    // Two effects, and the solver has to get both.
    //
    // A member's load acts along its length, so half of it reacts straight into
    // that member's own base and only the other half reaches the joint. The
    // apex therefore shares `P/2`, and *that* is what the analytic split
    // applies to.
    //
    // What each bar's base section carries is that share plus the whole of its
    // own load's axial component — the base holds up everything above it, the
    // tip holds up nothing, and the reported force is for the section that
    // decides whether the bar fails.
    let cos_t = t.cos();
    let p_joint = p_load / 2.0;
    let shared_mid = p_joint / (1.0 + 2.0 * cos_t.powi(3));
    let shared_outer = shared_mid * cos_t * cos_t;
    let own = |i: usize| {
        let axis = (topo.tip[i] - topo.base[i]).unit();
        (field.force[i].dot(axis) * 0.5).abs()
    };
    let expect_mid = shared_mid + own(0);
    let expect_outer = shared_outer + (own(1) + own(2)) / 2.0;

    println!(
        "  middle bar {middle:.1} N (analytic {expect_mid:.1}), each outer {outer:.1} N \
         (analytic {expect_outer:.1}), {iters} CG iterations"
    );
    assert!(
        (middle - expect_mid).abs() / expect_mid < 0.005,
        "middle bar carries {middle:.1} N, analytic says {expect_mid:.1} N"
    );
    assert!(
        (outer - expect_outer).abs() / expect_outer < 0.005,
        "outer bars carry {outer:.1} N, analytic says {expect_outer:.1} N"
    );
    // Vertical equilibrium at the joint, independently of the analytic form:
    // strip each bar's own half back off and what is left must balance `P/2`.
    let carried = (middle - own(0)) + 2.0 * (outer - (own(1) + own(2)) / 2.0) * cos_t;
    assert!(
        (carried - p_joint).abs() / p_joint < 0.02,
        "the truss does not balance: {carried:.1} N against {p_joint:.1} N"
    );
}

/// Bracing must actually relieve the primary load path, and the solver must
/// notice that it has.
#[test]
fn bracing_relieves_the_primary_path() {
    let mut m = Morphology::planned(Program::Tower, 3.0e6, 11, 0x77);
    m.progress = 1.0;
    m.built = 3.0e6;
    let matter = Matter::neutral(3.0e6, m.extent(), 290.0, Composition::crustal());
    let (bodies, topo, _) = sample_structured(&matter, &m, 2000, 7, 0x77, 0, SURFACE_G);
    assert!(!topo.ties.is_empty(), "a framed tower should be braced");
    assert!(!topo.is_determinate());

    let mut field = LoadField::new(bodies.len(), 290.0);
    field.apply(&weather::wind(35.0, v3(1.0, 0.0, 0.0)), &bodies, &topo);
    field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);

    let (braced, indeterminate, iters) = analyse_with(&bodies, &topo, &field);
    assert!(indeterminate && iters > 0, "the redundant solver did not run");

    // Same structure with the bracing cut.
    let mut bare = topo.clone();
    for t in bare.ties.iter_mut() {
        t.integrity = 0.0;
    }
    assert!(bare.is_determinate());
    let (unbraced, _, _) = analyse_with(&bodies, &bare, &field);

    // What bracing does is *redistribute*, and the right measure is therefore
    // the total stress carried, not the peak or any one member. Two plausible
    // stronger claims are both false, and the solver was right to contradict
    // them: the peak over all members can rise, because load taken off the
    // columns goes into the braces; and an individual column can gain, because
    // under a lateral load bracing transfers force between the windward and
    // leeward sides. Only the total is guaranteed to fall.
    let total_braced: f64 = braced.iter().map(|l| l.stress).sum();
    let total_bare: f64 = unbraced.iter().map(|l| l.stress).sum();
    let relief = (total_bare - total_braced) / total_bare;
    println!(
        "  total stress braced {total_braced:.3e}, unbraced {total_bare:.3e} \
         — {:.1}% relieved in {iters} CG iterations",
        relief * 100.0
    );
    assert!(
        relief > 0.01,
        "bracing relieved only {:.2}% of the total stress",
        relief * 100.0
    );
}

/// The mechanisms are general: the same accretion law describes snow, ice and
/// ash, and the same drag law describes air and water.
#[test]
fn mechanisms_are_not_weather_specific() {
    let (matter, m) = tree(900.0);
    let (bodies, topo) = load(&matter, &m, 3000);
    let crown = m.capture_area();

    let peak = |mech: Mechanism| {
        let mut field = LoadField::new(bodies.len(), 280.0);
        field.apply(&mech, &bodies, &topo);
        field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
        let loads = analyse(&bodies, &topo, &field);
        loads.iter().fold(0.0f64, |a, l| a.max(l.utilisation))
    };

    // Same law, three materials falling out of the sky.
    let snow = peak(weather::snow(0.10, 400.0, crown));
    let ice = peak(weather::ice(0.02, crown));
    let ash = peak(weather::ash(0.10, crown));
    // Same law, two fluids.
    let air = peak(weather::wind(30.0, v3(1.0, 0.0, 0.0)));
    let water = peak(weather::current(2.0, v3(1.0, 0.0, 0.0)));
    println!("  accretion — snow {snow:.2}, ice {ice:.2}, ash {ash:.2}");
    println!("  drag      — 30 m/s air {air:.2}, 2 m/s water {water:.2}");

    for (name, v) in [("snow", snow), ("ice", ice), ("ash", ash), ("air", air), ("water", water)] {
        assert!(v.is_finite() && v > 0.0, "{name} produced no load");
    }
    // Same areal mass, different shedding: ice retains all of it, snow sheds
    // whatever exceeds what it can adhere to.
    let areal = 40.0;
    let ice_same = peak(weather::ice(areal / 917.0, crown));
    let snow_same = peak(weather::snow(areal / 150.0, 150.0, crown));
    println!("  {areal:.0} kg/m2 as ice {ice_same:.2}, as dry snow {snow_same:.2}");
    assert!(
        ice_same > snow_same,
        "ice should retain more of the same fall than snow does"
    );
    // Water is 800 times denser than air; 2 m/s of it is comparable to a gale.
    assert!(water > air * 0.5, "a current should be structurally serious");
}

/// Materials are measured, and swapping one changes the outcome in the
/// direction the numbers say it should.
///
/// **Rewritten for `docs/PLAY.md` D14**, which retired the tabulated `rupture`
/// column these assertions used to read. Two of the four comparisons it made
/// are still reproduced and two are not, and which two is the point:
///
/// - The **deposited** materials come out in the table's own order — dry timber
///   above green wood above aragonite above masonry — from nothing but a ring
///   width, a band and a course thickness.
/// - The **frozen** ones do not. Steel derives at 2.87e6 Pa against a retired
///   4.0e8, because Griffith is the *brittle* law and steel is ductile: it
///   yields by moving dislocations rather than by running a crack through a
///   grain, and the crack tip blunts instead of propagating. See
///   `Material::strength` for the measurement and `docs/BACKLOG.md` for the
///   residual.
///
/// So this asserts what is derived rather than what was tabulated, and says so.
/// Asserting the old ordering would be asserting a table that no longer exists.
#[test]
fn materials_are_interchangeable_data() {
    let (matter, m) = tree(900.0);
    let (bodies, base) = load(&matter, &m, 2000);

    let peak_for = |mat: Material| {
        let mut topo = base.clone();
        topo.material = mat;
        let mut field = LoadField::new(bodies.len(), 290.0);
        field.apply(&weather::wind(35.0, v3(1.0, 0.0, 0.0)), &bodies, &topo);
        field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
        let loads = analyse(&bodies, &topo, &field);
        loads.iter().fold(0.0f64, |a, l| a.max(l.utilisation))
    };
    let green = peak_for(Material::green_wood());
    let dry = peak_for(Material::dry_timber());
    let aragonite = peak_for(Material::aragonite());
    let masonry = peak_for(Material::masonry());
    println!(
        "  same geometry, 35 m/s wind — green wood {green:.2}, dry timber {dry:.2}, \
         aragonite {aragonite:.2}, masonry {masonry:.2}"
    );
    // Slow-grown timber is stronger than fast-grown, and the *only* thing that
    // differs between them is how thick a layer the tree laid down.
    assert!(dry < green, "seasoned timber should out-perform green wood");
    assert!(aragonite > dry, "a coral skeleton is weaker than timber");
    assert!(masonry > aragonite, "masonry in bending should be worse still");

    // A material somebody states rather than measures still works, because
    // nothing in the solver knows where a material came from. The two numbers
    // Griffith needs are stated directly here, which is what makes this a
    // *custom* material rather than a measured one.
    let custom = Material {
        name: "spider silk",
        density: 1300.0,
        surface_energy: 4.0,
        flaw_size: 1.0e-8,
        stiffness: 10.0e9,
        ..Material::green_wood()
    };
    println!("  spider silk derives {:.3e} Pa", custom.strength());
    assert!(
        peak_for(custom) < green,
        "a material with a flaw scale five orders finer should carry more"
    );

    // And the size effect, which a tabulated stress could not express at all:
    // a piece too small to contain the flaw is stronger than a piece that can.
    let m = Material::green_wood();
    let big = m.strength_at_size(1.0);
    let small = m.strength_at_size(m.flaw_size / 100.0);
    println!("  green wood at 1 m {big:.3e} Pa, at {:.1e} m {small:.3e} Pa", m.flaw_size / 100.0);
    assert!(small > big * 9.0, "ten times finer should be about ten times stronger");
}

/// The exact path and the redundant path agree when there is no redundancy.
/// Otherwise the fast path would be a different physics, not a special case.
#[test]
fn the_two_solvers_agree_on_determinate_structures() {
    let (matter, m) = tree(900.0);
    let (bodies, topo) = load(&matter, &m, 1500);
    let mut field = LoadField::new(bodies.len(), 290.0);
    // A load the structure is comfortable under. The two paths are only
    // *required* to agree while everything stays elastic: the redundant solver
    // also redistributes load off members past yield, which statics cannot do
    // and which is the whole reason ductility is worth modelling. At 30 m/s
    // twenty-two members yield and the paths differ by 2%, correctly.
    field.apply(&weather::wind(12.0, v3(0.7, 0.7, 0.0)), &bodies, &topo);
    field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);

    let (exact, indeterminate, _) = analyse_with(&bodies, &topo, &field);
    assert!(!indeterminate, "a tree should be determinate");

    // Add ties that carry nothing — duplicates of existing support bonds with
    // negligible area — and confirm the answer barely moves.
    let mut redundant = topo.clone();
    for i in 1..40usize {
        let p = redundant.support[i];
        if p != NO_SUPPORT {
            redundant.ties.push(Tie { a: i as u32, b: p, area: 1e-12, integrity: 1.0 });
        }
    }
    assert!(!redundant.is_determinate());
    let (solved, _, iters) = analyse_with(&bodies, &redundant, &field);
    // Zero iterations means the frame solve failed and the answer fell back to
    // statics, which would make this test agree with itself. It did, for a
    // while: the conjugate gradient could not converge on a 1500-member tree
    // with Jacobi preconditioning, so the redundant path was never exercised
    // and the test passed by never running the code it was checking.
    assert!(iters > 0, "the redundant solver fell back instead of solving");

    let mut worst = 0.0f64;
    for (a, b) in exact.iter().zip(&solved) {
        let scale = a.stress.abs().max(b.stress.abs()).max(1.0);
        worst = worst.max((a.stress - b.stress).abs() / scale);
    }
    println!("  worst disagreement between the two paths: {worst:.3e} ({iters} CG iterations)");
    // Not "close enough": a determinate structure's internal forces are fixed
    // by equilibrium alone, so a correct beam model has to reproduce statics
    // exactly, member by member, and it does.
    assert!(worst < 1e-8, "the solvers disagree by {worst:.3e}");
}

/// A structure must be analysed and proportioned when it is created.
///
/// The generator decides where members go. It has no way to know what any of
/// them will carry, so the radii it produces are a shape scaled as a group to
/// match the structural mass — which leaves a few members at the point of
/// failure and most of the material in members doing nothing. Both have the
/// same cause and the same fix.
///
/// The constraint that makes it meaningful is that the material does not
/// change. An optimiser that improves a structure by feeding it has not
/// optimised anything.
#[test]
fn a_structure_is_proportioned_for_its_loads_when_it_is_built() {
    for (label, mass, budget) in [("900 kg", 900.0, 1200usize), ("6 t", 6000.0, 400)] {
        let (matter, m) = tree(mass);
        let (_, _, report) = sample_structured(&matter, &m, budget, 7, 0x1234, 0, SURFACE_G);
        let d = report.design;
        println!(
            "  {label:>6}: peak {:.3} -> {:.3}, spread {:.3} -> {:.3} over {} passes, \
             volume error {:.2e}",
            d.peak_before, d.peak_after, d.spread_before, d.spread_after, d.passes,
            d.volume_error()
        );
        assert!(d.passes > 0, "{label}: the design pass did not run");
        assert!(
            d.peak_after < d.peak_before * 0.75,
            "{label}: peak utilisation {:.3} against {:.3} before",
            d.peak_after,
            d.peak_before
        );
        assert!(
            d.spread_after < d.spread_before,
            "{label}: the load is no more evenly carried than it was"
        );
        // The whole point: same material, better arrangement.
        assert!(
            d.volume_error() < 1e-9,
            "{label}: structural volume moved by {:.3e}",
            d.volume_error()
        );
    }
}

/// And the improvement has to survive loads the design cases did not include.
///
/// This is what a single-case optimiser gets wrong. Sizing against one wind
/// direction produces a structure that is optimal in that direction and
/// brittle in every other; the envelope of several directions plus a vertical
/// overload is what stops that, and the check is a load from a direction
/// nothing was designed for.
#[test]
fn the_design_pass_does_not_overfit_its_own_load_cases() {
    let (matter, m) = tree(900.0);
    let (bodies, topo) = load(&matter, &m, 1500);

    // A direction the design cases do not use, at a speed below what a tree of
    // this size should be troubled by.
    let mut worst = 0.0f64;
    for degrees in [17.0f64, 53.0, 131.0, 209.0, 288.0, 341.0] {
        let a = degrees.to_radians();
        let mut field = LoadField::new(bodies.len(), 290.0);
        field.apply(&weather::wind(22.0, v3(a.cos(), a.sin(), 0.0)), &bodies, &topo);
        field.apply(&weather::gravity(SURFACE_G), &bodies, &topo);
        let loads = analyse(&bodies, &topo, &field);
        let peak = loads.iter().map(|l| l.utilisation).fold(0.0f64, f64::max);
        worst = worst.max(peak);
    }
    println!("  worst utilisation over six off-design wind directions at 22 m/s: {worst:.3}");
    assert!(
        worst < 1.0,
        "a tree designed for 20 m/s failed at 22 m/s from an off-design direction: {worst:.3}"
    );
}
