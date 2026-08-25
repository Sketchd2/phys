//! Render a tree, then load it to failure four different ways.
//!
//! Nothing here scripts what damage looks like. Each insult produces forces,
//! temperatures or deposited energy; the same stress calculation decides what
//! survives; and the renderer draws whatever is left.

use phys::engine::{default_spec, galaxy, World};
use phys::math::v3;
use phys::morph::*;
use phys::render::*;
use phys::solvers::structure::*;
use phys::state::*;
use phys::units::*;

fn plant(seed: u64, mass: f64) -> (World, phys::ids::NodeIdx) {
    let mut w = World::new(galaxy(seed, 1e9), 20.0);
    let root = w.tree.root;
    w.tree.refine(root);
    let node = w.tree.promote(root, 7, default_spec(Tier::Stellar));
    {
        let n = &mut w.tree.nodes[node.get()];
        // A patch of ground: soil, water and air, of which the tree will build
        // itself a few tonnes. The reservoir is what bounds how big it gets.
        n.agg = Aggregate::neutral(mass, 6.0, 291.0, Program::Tree.substrate());
        n.spec.count = 9000;
    }
    w.plant(node, Program::Tree, Environment::default());
    (w, node)
}

/// Render whatever the node currently holds.
fn shoot(
    w: &mut World,
    node: phys::ids::NodeIdx,
    path: &str,
    style: Style,
    caption: &str,
    fixed: Option<Camera>,
) -> Camera {
    let bodies = w.tree.refine(node).to_vec();
    let topo = match w.tree.nodes[node.get()].topology.clone() {
        Some(t) => t,
        None => return fixed.unwrap_or(Camera::framing(v3(0.0, 0.0, 0.0), 1.0, 0.6, 0.12)),
    };
    // Frame on the structure's own extent.
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    let mut radius: f64 = 0.0;
    for (i, b) in bodies.iter().enumerate() {
        if i < topo.bonds.len() && (topo.tip[i] - topo.base[i]).norm2() > 0.0 {
            lo = lo.min(b.pos.z);
            hi = hi.max(b.pos.z);
            radius = radius.max(b.pos.norm());
        }
    }
    // One camera for every frame in the set. Auto-framing each shot rescales
    // the subject and hides exactly what the comparison is meant to show — a
    // tree that lost a third of its height looked identical.
    let cam = fixed.unwrap_or_else(|| {
        let centre = v3(0.0, 0.0, (lo + hi) * 0.5);
        Camera::framing(centre, radius.max(1.0) * 1.05, 0.6, 0.12)
    });
    let intact = vec![true; bodies.len()];
    let mut canvas = Canvas::new(760, 620);
    draw_structure(&mut canvas, &cam, &bodies, &topo, &intact, &style);
    write_png(&canvas, path).expect("write png");
    let m = w.tree.nodes[node.get()].morphology.as_ref().unwrap();
    println!(
        "  {:<24} {:>5} members  {:>8.1} kg standing  {:>6.1} m  -> {}",
        caption,
        topo.bonds.len(),
        m.built,
        m.tree_height(),
        path
    );
    cam
}

fn main() {
    println!("Growing a tree for 90 years, then loading it to failure.\n");

    // --- healthy ---------------------------------------------------------
    let (mut w, node) = plant(0xA11CE, 60_000.0);
    for _ in 0..(90 * 12) {
        w.grow_node(node, YEAR / 12.0);
    }
    let cam = shoot(&mut w, node, "render_healthy.png", Style::daylight(), "healthy, 90 years", None);

    let grown = w.tree.nodes[node.get()].morphology.as_ref().unwrap().built;
    let crown = w.tree.nodes[node.get()].morphology.as_ref().unwrap().capture_area();

    // --- wet snow --------------------------------------------------------
    let (mut w2, n2) = plant(0xA11CE, 60_000.0);
    for _ in 0..(90 * 12) {
        w2.grow_node(n2, YEAR / 12.0);
    }
    let snow = w2.damage(n2, &[weather::snow(0.35, 480.0, crown)]);
    shoot(&mut w2, n2, "render_snow.png", Style::daylight(), "after 350 mm wet snow", Some(cam));
    println!(
        "      {} joints failed, {:.0} kg down, peak utilisation {:.2}",
        snow.broken_joints, snow.detached_mass, snow.peak_utilisation
    );

    // --- lightning -------------------------------------------------------
    let (mut w3, n3) = plant(0xA11CE, 60_000.0);
    for _ in 0..(90 * 12) {
        w3.grow_node(n3, YEAR / 12.0);
    }
    // Enter at a structural member well up in the crown.
    w3.tree.refine(n3);
    let structural = w3.tree.nodes[n3.get()]
        .topology
        .as_ref()
        .map(|t| t.bonds.iter().filter(|b| b.radius > 0.0).count())
        .unwrap_or(1);
    let strike = w3.damage(n3, &[weather::lightning(2.0e9, (structural / 3) as u32)]);
    shoot(&mut w3, n3, "render_lightning.png", Style::default(), "after a 2 GJ strike", Some(cam));
    println!(
        "      {} joints destroyed on the channel, {:.0} kg down, {:.2e} J delivered",
        strike.broken_joints, strike.detached_mass, strike.energy_delivered
    );

    // --- fire ------------------------------------------------------------
    let (mut w4, n4) = plant(0xA11CE, 60_000.0);
    for _ in 0..(90 * 12) {
        w4.grow_node(n4, YEAR / 12.0);
    }
    let fire = w4.damage(n4, &[weather::fire(1050.0, 30.0, 420.0)]);
    shoot(&mut w4, n4, "render_fire.png", Style::burned(), "after a crown fire", Some(cam));
    println!(
        "      {:.0} kg consumed, {:.3e} J of stored energy released as heat",
        fire.consumed_mass, fire.energy_released
    );

    println!("\n  Grown mass before damage: {grown:.0} kg, crown {crown:.0} m^2");
    println!("  Every outcome above came from the stress calculation, not a script.");
}
