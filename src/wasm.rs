//! A C ABI over the engine, so a browser can drive the real thing.
//!
//! Not a port and not a recording. The `.wasm` this compiles to *is* the engine
//! — the same growth model, the same conservation projection, the same
//! structural solver validated against the analytic truss. A viewer built on it
//! is watching physics happen, not replaying frames somebody exported.
//!
//! Deliberately a raw C ABI rather than a binding generator: nothing here needs
//! a dependency, and the interface is a dozen functions and a few slices of
//! linear memory. The host reads geometry straight out of the module's memory
//! with a typed array.
//!
//! Everything is `f32` on the wire — positions in a node's local frame span a
//! few orders of magnitude, which `f32` holds comfortably, and halving the
//! bytes crossing the boundary matters more than digits nobody will look at.

#![cfg(target_arch = "wasm32")]

use crate::engine::{default_spec, galaxy, World};
use crate::ids::NodeIdx;
use crate::math::{v3, Vec3};
use crate::morph::{Environment, Program};
use crate::solvers::structure::{self as st, LoadField, Mechanism};
use crate::state::Aggregate;
use crate::units::{Tier, YEAR};

/// Everything the host needs to keep between calls.
pub struct Session {
    world: World,
    node: NodeIdx,
    /// Interleaved member geometry: base xyz, tip xyz, radius, stress ratio.
    geometry: Vec<f32>,
    /// Scalar readouts, see `readout_*` indices below.
    readouts: [f32; 16],
    /// Standing conditions, re-applied every step.
    wind: f64,
    wind_dir: Vec3,
    snow_depth: f64,
    snow_density: f64,
    fire_temperature: f64,
    fire_height: f64,
    budget: usize,
    dirty: bool,
    /// Joints broken since the world was created. The per-call count is useless
    /// to watch: a standing load breaks what it can on the first frame and
    /// nothing after, so the live number reads zero while the structure is
    /// visibly coming apart.
    total_broken: u32,
}

static mut SESSION: Option<Session> = None;

fn session() -> &'static mut Session {
    // Single-threaded by construction: WebAssembly modules instantiated by the
    // viewer have one thread, and every entry point below runs to completion
    // before the host can call again.
    unsafe {
        let s = &raw mut SESSION;
        (*s).as_mut().expect("create() first")
    }
}

/// Create a world with one structure in it, ready to grow.
#[unsafe(no_mangle)]
pub extern "C" fn create(seed: u32, program: u32, reservoir_kg: f32, budget: u32) {
    let mut world = World::new(galaxy(seed as u64 | 0x5EED_0000, 1e9), 20.0);
    let root = world.tree.root;
    world.tree.refine(root);
    let node = world.tree.promote(root, 7, default_spec(Tier::Stellar));
    let prog = match program {
        1 => Program::Coral,
        2 => Program::Tower,
        3 => Program::Wall,
        _ => Program::Tree,
    };
    {
        let n = &mut world.tree.nodes[node.get()];
        n.agg = Aggregate::neutral(reservoir_kg as f64, 6.0, 291.0, prog.substrate());
        n.spec.count = budget as usize;
    }
    // The environment override is what the node actually grows in, so a
    // planned structure's labour rate has to go *in* it. Setting
    // `world.labour_rate` alone was silently ignored, and the tower never left
    // its foundations.
    let env = Environment {
        labour: if prog.is_planned() { 1.0 / (90.0 * 86400.0) } else { 0.0 },
        ..Environment::default()
    };
    world.plant(node, prog, env);
    if prog.is_planned() {
        world.labour_rate = env.labour;
        if let Some(m) = world.tree.nodes[node.get()].morphology.as_mut() {
            m.design_mass = reservoir_kg as f64 * 0.5;
        }
    }
    unsafe {
        let s = &raw mut SESSION;
        *s = Some(Session {
            world,
            node,
            geometry: Vec::new(),
            readouts: [0.0; 16],
            wind: 0.0,
            wind_dir: v3(1.0, 0.0, 0.0),
            snow_depth: 0.0,
            snow_density: 300.0,
            fire_temperature: 0.0,
            fire_height: 0.0,
            budget: budget as usize,
            dirty: true,
            total_broken: 0,
        });
    }
}

/// Advance the world by `years` of simulated time, then apply standing
/// conditions and resolve any failures.
#[unsafe(no_mangle)]
pub extern "C" fn step(years: f32) {
    let s = session();
    let dt = years as f64 * YEAR;
    if dt > 0.0 {
        // Growth runs on the aggregate — the structure is not materialised for
        // this, however large it is.
        s.world.grow_node(s.node, dt);
        s.dirty = true;
    }

    let mut mechanisms: Vec<Mechanism> = Vec::new();
    if s.wind > 0.0 {
        mechanisms.push(st::weather::wind(s.wind, s.wind_dir));
    }
    if s.snow_depth > 0.0 {
        let crown = s
            .world
            .tree
            .nodes[s.node.get()]
            .morphology
            .as_ref()
            .map(|m| m.capture_area())
            .unwrap_or(0.0);
        mechanisms.push(st::weather::snow(s.snow_depth, s.snow_density, crown));
    }
    if s.fire_temperature > 0.0 {
        mechanisms.push(st::weather::fire(s.fire_temperature, s.fire_height, dt.max(1.0).min(60.0)));
    }
    if !mechanisms.is_empty() {
        let before = s.world.tree.nodes[s.node.get()]
            .morphology
            .as_ref()
            .map(|m| m.built)
            .unwrap_or(0.0);
        let out = s.world.damage(s.node, &mechanisms);
        let after = s.world.tree.nodes[s.node.get()]
            .morphology
            .as_ref()
            .map(|m| m.built)
            .unwrap_or(0.0);
        s.total_broken += out.broken_joints as u32;
        s.readouts[5] = out.peak_utilisation as f32;
        s.readouts[7] += (before - after).max(0.0) as f32;
        s.readouts[11] = if out.indeterminate { 1.0 } else { 0.0 };
        s.readouts[12] = out.solver_iterations as f32;
        s.dirty = true;
    }
    refresh(s);
}

/// Fire a discharge into the structure at a fraction of the way up.
#[unsafe(no_mangle)]
pub extern "C" fn strike(joules: f32, height_fraction: f32) {
    let s = session();
    s.world.tree.refine(s.node);
    let structural = s.world.tree.nodes[s.node.get()]
        .topology
        .as_ref()
        .map(|t| t.bonds.iter().filter(|b| b.radius > 0.0).count())
        .unwrap_or(1);
    let entry = ((structural as f32 * height_fraction.clamp(0.0, 0.99)) as usize).min(structural.saturating_sub(1));
    let before = s.world.tree.nodes[s.node.get()]
        .morphology.as_ref().map(|m| m.built).unwrap_or(0.0);
    let out = s
        .world
        .damage(s.node, &[st::weather::lightning(joules as f64, entry as u32)]);
    let after = s.world.tree.nodes[s.node.get()]
        .morphology.as_ref().map(|m| m.built).unwrap_or(0.0);
    s.total_broken += out.broken_joints as u32;
    s.readouts[7] += (before - after).max(0.0) as f32;
    s.dirty = true;
    refresh(s);
}

/// Set the standing wind.
#[unsafe(no_mangle)]
pub extern "C" fn set_wind(speed: f32, angle: f32) {
    let s = session();
    s.wind = speed as f64;
    s.wind_dir = v3(angle.cos() as f64, angle.sin() as f64, 0.0);
}

/// Set standing snow. Density distinguishes powder from wet snow, and it is the
/// distinction that decides whether anything breaks.
#[unsafe(no_mangle)]
pub extern "C" fn set_snow(depth: f32, density: f32) {
    let s = session();
    s.snow_depth = depth as f64;
    s.snow_density = density.max(50.0) as f64;
}

/// Set a fire front. Zero temperature puts it out.
#[unsafe(no_mangle)]
pub extern "C" fn set_fire(temperature: f32, flame_height: f32) {
    let s = session();
    s.fire_temperature = temperature as f64;
    s.fire_height = flame_height as f64;
}

/// Recompute geometry and readouts if anything changed.
fn refresh(s: &mut Session) {
    if !s.dirty {
        return;
    }
    s.dirty = false;
    let node = s.node;
    let bodies = s.world.tree.refine(node).to_vec();
    let topo = match s.world.tree.nodes[node.get()].topology.clone() {
        Some(t) => t,
        None => return,
    };

    // Current stress, so the viewer can colour by how hard each member is
    // working rather than only showing what has already failed.
    let ambient = s.world.tree.nodes[node.get()].agg.temperature;
    let mut field = LoadField::new(bodies.len(), ambient);
    if s.wind > 0.0 {
        field.apply(&st::weather::wind(s.wind, s.wind_dir), &bodies, &topo);
    }
    if s.snow_depth > 0.0 {
        let crown = s.world.tree.nodes[node.get()]
            .morphology
            .as_ref()
            .map(|m| m.capture_area())
            .unwrap_or(0.0);
        field.apply(
            &st::weather::snow(s.snow_depth, s.snow_density, crown),
            &bodies,
            &topo,
        );
    }
    field.apply(&st::weather::gravity(), &bodies, &topo);
    let (loads, indeterminate, iters) = st::analyse_with(&bodies, &topo, &field);

    s.geometry.clear();
    let mut peak = 0.0f32;
    for i in 0..bodies.len().min(topo.bonds.len()) {
        if topo.bonds[i].radius <= 0.0 {
            continue;
        }
        let u = loads.get(i).map(|l| l.utilisation).unwrap_or(0.0) as f32;
        peak = peak.max(u);
        s.geometry.extend_from_slice(&[
            topo.base[i].x as f32,
            topo.base[i].y as f32,
            topo.base[i].z as f32,
            topo.tip[i].x as f32,
            topo.tip[i].y as f32,
            topo.tip[i].z as f32,
            topo.bonds[i].radius as f32,
            u,
        ]);
    }

    let n = &s.world.tree.nodes[node.get()];
    let m = n.morphology.as_ref();
    s.readouts[0] = m.map(|m| m.built).unwrap_or(0.0) as f32;
    s.readouts[1] = m.map(|m| m.height()).unwrap_or(0.0) as f32;
    s.readouts[2] = m.map(|m| m.age / YEAR).unwrap_or(0.0) as f32;
    s.readouts[3] = (s.geometry.len() / 8) as f32;
    s.readouts[4] = m.map(|m| m.state_bytes()).unwrap_or(0) as f32;
    s.readouts[5] = peak;
    s.readouts[8] = n.agg.chemical_energy as f32;
    s.readouts[9] = n.agg.entropy_exported as f32;
    s.readouts[10] = m.map(|m| m.progress).unwrap_or(0.0) as f32;
    s.readouts[11] = if indeterminate { 1.0 } else { 0.0 };
    s.readouts[12] = iters as f32;
    s.readouts[13] = (n.agg.mass - m.map(|m| m.built).unwrap_or(0.0)) as f32;
    s.readouts[6] = s.total_broken as f32;
    s.readouts[14] = m.map(|m| m.events.len()).unwrap_or(0) as f32;
    s.readouts[15] = n.agg.temperature as f32;
}

/// Pointer to the interleaved member geometry in linear memory.
#[unsafe(no_mangle)]
pub extern "C" fn geometry_ptr() -> *const f32 {
    session().geometry.as_ptr()
}

/// Number of members currently in the geometry buffer.
#[unsafe(no_mangle)]
pub extern "C" fn geometry_len() -> u32 {
    (session().geometry.len() / 8) as u32
}

/// Pointer to the readout array.
#[unsafe(no_mangle)]
pub extern "C" fn readouts_ptr() -> *const f32 {
    session().readouts.as_ptr()
}

/// Materialisation budget, so the host can trade detail for frame time.
#[unsafe(no_mangle)]
pub extern "C" fn set_budget(budget: u32) {
    let s = session();
    s.budget = budget as usize;
    s.world.tree.nodes[s.node.get()].spec.count = s.budget;
    s.world.tree.nodes[s.node.get()].bodies.clear();
    s.world.tree.nodes[s.node.get()].topology = None;
    s.dirty = true;
    refresh(s);
}

/// Discard the structure's detail entirely — the engine's normal behaviour when
/// nobody is looking. The next frame regenerates it from the developmental
/// state, which is the claim worth being able to watch.
#[unsafe(no_mangle)]
pub extern "C" fn coarsen() -> u32 {
    let s = session();
    let before = s.world.tree.nodes[s.node.get()].bodies.len() as u32;
    s.world.tree.coarsen(s.node);
    s.dirty = true;
    before
}

/// Bytes of detail currently materialised.
#[unsafe(no_mangle)]
pub extern "C" fn detail_bytes() -> u32 {
    session().world.tree.detail_bytes() as u32
}
