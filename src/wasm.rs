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
    readouts: [f32; 24],
    /// Standing conditions, re-applied every step.
    wind: f64,
    wind_dir: Vec3,
    snow_depth: f64,
    snow_density: f64,
    fire_temperature: f64,
    fire_height: f64,
    budget: usize,
    dirty: bool,
    /// Whether the structure is being integrated through real time.
    dynamic: bool,
    /// Wall-clock seconds of dynamics elapsed, for the gust model.
    shake_time: f64,
    /// Current turbulent fluctuation about the mean wind, m/s.
    gust: f64,
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
            readouts: [0.0; 24],
            wind: 0.0,
            wind_dir: v3(1.0, 0.0, 0.0),
            snow_depth: 0.0,
            snow_density: 300.0,
            fire_temperature: 0.0,
            fire_height: 0.0,
            budget: budget as usize,
            dirty: true,
            dynamic: false,
            shake_time: 0.0,
            gust: 0.0,
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


// ---------------------------------------------------------------------------
// Real time
// ---------------------------------------------------------------------------

/// Turn dynamics on or off.
///
/// With it off the viewer shows the structure where a quasi-static analysis
/// puts it, which is where it would sit under a load held forever. With it on
/// the structure has mass, and what is drawn is where it actually is at this
/// instant — leaning into a gust, swinging back out of one, ringing after a
/// limb goes.
#[unsafe(no_mangle)]
pub extern "C" fn set_dynamic(on: u32) {
    let s = session();
    s.dynamic = on != 0;
    if !s.dynamic {
        s.world.settle();
        s.dirty = true;
        refresh(s);
    }
}

/// Advance the structure by `seconds` of real time.
///
/// Separate from [`step`], which advances *growth* by years. The two run on
/// wildly different clocks and it would be a mistake to couple them: a tree
/// sways with a period of seconds and grows over decades, and a viewer wants
/// to watch one while fast-forwarding the other.
#[unsafe(no_mangle)]
pub extern "C" fn shake(seconds: f32) {
    let s = session();
    if !s.dynamic || seconds <= 0.0 {
        return;
    }
    let dt = (seconds as f64).min(0.25);
    s.shake_time += dt;

    // Turbulence as an Ornstein-Uhlenbeck process: the standard first-order
    // model of a gust, and three lines. The fluctuation decays towards zero
    // with a time constant of a few seconds and is kicked by noise whose size
    // is set by the mean speed, which is what gives the low-frequency energy a
    // real boundary layer has and a sine wave does not.
    if s.wind > 0.0 {
        let tau = 4.0;
        let intensity = 0.18 * s.wind;
        let mut stream = crate::rng::Stream::at(
            s.world.tree.world_seed,
            0,
            0,
            crate::rng::Purpose::ThermalNoise,
        )
        .split((s.shake_time * 1000.0) as u64);
        let decay = (-dt / tau).exp();
        s.gust = s.gust * decay
            + intensity * (1.0 - decay * decay).sqrt() * stream.normal();
    } else {
        s.gust = 0.0;
    }

    let mut mechanisms: Vec<Mechanism> = Vec::new();
    let speed = (s.wind + s.gust).max(0.0);
    if speed > 0.0 {
        mechanisms.push(st::weather::wind(speed, s.wind_dir));
    }
    if s.snow_depth > 0.0 {
        let crown = s.world.tree.nodes[s.node.get()]
            .morphology
            .as_ref()
            .map(|m| m.capture_area())
            .unwrap_or(0.0);
        mechanisms.push(st::weather::snow(s.snow_depth, s.snow_density, crown));
    }

    let out = s.world.shake(s.node, &mechanisms, dt);
    s.total_broken += out.broken_joints as u32;
    s.readouts[12] = out.iterations as f32;
    s.readouts[16] = speed as f32;
    s.readouts[17] = out.displacement as f32;
    s.readouts[18] = (out.kinetic + out.strain) as f32;
    s.readouts[19] = out.displacement_ratio as f32;
    s.readouts[6] = s.total_broken as f32;
    if out.broken_joints > 0 {
        s.dirty = true;
        refresh(s);
        return;
    }
    write_deformed(s);
}

/// Overwrite the geometry buffer with where the structure actually is.
///
/// Only the endpoints move; the utilisation in the eighth slot is left at
/// whatever the last analysis put there, so the colouring stays meaningful
/// between analyses rather than flickering with every substep.
fn write_deformed(s: &mut Session) {
    let Some(ds) = s.world.shaken(s.node) else {
        return;
    };
    let members = ds.deformed_members();
    let topo = match s.world.tree.nodes[s.node.get()].topology.as_ref() {
        Some(t) => t,
        None => return,
    };
    let mut slot = 0usize;
    for i in 0..members.len().min(topo.bonds.len()) {
        if topo.bonds[i].radius <= 0.0 {
            continue;
        }
        let base = slot * 8;
        slot += 1;
        if base + 5 >= s.geometry.len() {
            break;
        }
        let (b, t) = members[i];
        s.geometry[base] = b.x as f32;
        s.geometry[base + 1] = b.y as f32;
        s.geometry[base + 2] = b.z as f32;
        s.geometry[base + 3] = t.x as f32;
        s.geometry[base + 4] = t.y as f32;
        s.geometry[base + 5] = t.z as f32;
    }
}

// ---------------------------------------------------------------------------
// A stand of trees
// ---------------------------------------------------------------------------

/// A field of independently grown trees under one wind.
///
/// The single-structure session above answers "what does this thing do". This
/// answers a different question, and the one that actually shows whether the
/// engine's claims hold: given twenty structures that were never designed,
/// grown from twenty different seeds to twenty different sizes, does a *shared*
/// gust produce twenty different responses — and are the differences the ones
/// the physics would give?
///
/// Nothing here is per-tree tuning. Each tree gets a seed and a reservoir; its
/// height, its taper, the period it sways at and the wind that takes its limbs
/// all follow from that.
pub struct Forest {
    world: World,
    trees: Vec<Stand>,
    /// Interleaved member geometry: base xyz, tip xyz, radius, utilisation.
    geometry: Vec<f32>,
    readouts: [f32; 24],
    /// Mean wind speed, m/s.
    wind: f64,
    /// Turbulence intensity as a fraction of the mean.
    turbulence: f64,
    wind_dir: Vec3,
    /// Current gust, m/s, one per tree — a gust front is not everywhere at once.
    gust: Vec<f64>,
    elapsed: f64,
    budget: usize,
    dynamic: bool,
    total_broken: u32,
    shed: f64,
    /// Floats of standing geometry, before any debris is appended.
    standing_floats: usize,
}

struct Stand {
    node: NodeIdx,
    /// Where the tree stands, metres.
    at: Vec3,
    /// Members currently drawn for it.
    drawn: usize,
    height: f32,
    sway: f32,
    utilisation: f32,
}

static mut FOREST: Option<Forest> = None;

fn forest() -> &'static mut Forest {
    unsafe {
        let f = &raw mut FOREST;
        (*f).as_mut().expect("create_forest() first")
    }
}

/// Plant `count` trees, scattered over a square `extent` metres across.
///
/// Positions and reservoirs come from the world's own deterministic stream, so
/// the same seed gives the same field — including which tree is the big one on
/// the windward edge that goes first.
#[unsafe(no_mangle)]
pub extern "C" fn create_forest(seed: u32, count: u32, extent: f32, budget: u32) {
    use crate::rng::{Purpose, Stream};

    let world_seed = seed as u64 | 0x5EED_0000;
    let mut world = World::new(galaxy(world_seed, 1e9), 20.0);
    let root = world.tree.root;
    world.tree.refine(root);

    let count = count.clamp(1, 64) as usize;
    let per_tree = (budget as usize / count).max(24);
    let mut trees = Vec::with_capacity(count);
    let mut stream = Stream::at(world_seed, 0, 0, Purpose::Structure);

    for i in 0..count {
        let node = world.tree.promote(root, i, default_spec(Tier::Stellar));
        if node.is_none() {
            continue;
        }
        // Poisson-ish scatter: a jittered grid, so trees do not pile up and do
        // not look planted either.
        let side = (count as f64).sqrt().ceil();
        let cell = extent as f64 / side;
        let (gx, gy) = ((i as f64 % side), (i as f64 / side).floor());
        let jx = (stream.uniform() - 0.5) * cell * 0.8;
        let jy = (stream.uniform() - 0.5) * cell * 0.8;
        let at = v3(
            (gx + 0.5) * cell - extent as f64 * 0.5 + jx,
            (gy + 0.5) * cell - extent as f64 * 0.5 + jy,
            0.0,
        );
        // A real stand is not uniform: a few big trees, many small ones.
        let reservoir = 8000.0 * (0.35 + 1.9 * stream.uniform().powi(2));
        {
            let n = &mut world.tree.nodes[node.get()];
            n.agg = Aggregate::neutral(reservoir, 6.0, 291.0, Program::Tree.substrate());
            n.spec.count = per_tree;
        }
        world.plant(node, Program::Tree, Environment::default());
        trees.push(Stand {
            node,
            at,
            drawn: 0,
            height: 0.0,
            sway: 0.0,
            utilisation: 0.0,
        });
    }

    unsafe {
        let f = &raw mut FOREST;
        *f = Some(Forest {
            world,
            trees,
            geometry: Vec::new(),
            readouts: [0.0; 24],
            wind: 0.0,
            turbulence: 0.25,
            wind_dir: v3(1.0, 0.0, 0.0),
            gust: vec![0.0; count],
            elapsed: 0.0,
            budget: per_tree,
            dynamic: true,
            total_broken: 0,
            shed: 0.0,
            standing_floats: 0,
        });
    }
    refresh_forest(forest());
}

/// Advance every tree's growth by `years`.
#[unsafe(no_mangle)]
pub extern "C" fn grow_forest(years: f32) {
    let f = forest();
    let dt = years as f64 * YEAR;
    if dt <= 0.0 {
        return;
    }
    let nodes: Vec<NodeIdx> = f.trees.iter().map(|t| t.node).collect();
    for node in nodes {
        f.world.grow_node(node, dt);
    }
    f.world.settle();
    refresh_forest(f);
}

/// Set the wind: a mean speed, a turbulence intensity, and a direction.
#[unsafe(no_mangle)]
pub extern "C" fn set_forest_wind(speed: f32, turbulence: f32, angle: f32) {
    let f = forest();
    f.wind = speed.max(0.0) as f64;
    f.turbulence = turbulence.clamp(0.0, 1.0) as f64;
    f.wind_dir = v3(angle.cos() as f64, angle.sin() as f64, 0.0);
}

/// Turn the dynamics on or off.
#[unsafe(no_mangle)]
pub extern "C" fn set_forest_dynamic(on: u32) {
    let f = forest();
    f.dynamic = on != 0;
    if !f.dynamic {
        f.world.settle();
        refresh_forest(f);
    }
}

/// Advance the whole stand by `seconds` of real time under the current wind.
#[unsafe(no_mangle)]
pub extern "C" fn shake_forest(seconds: f32) {
    use crate::rng::{Purpose, Stream};
    let f = forest();
    if seconds <= 0.0 {
        return;
    }
    let dt = (seconds as f64).min(0.2);
    f.elapsed += dt;

    // Turbulence as an Ornstein-Uhlenbeck process per tree, correlated by a
    // shared gust front travelling downwind. A single fluctuation applied to
    // every tree at once would make the whole stand lean together, which is
    // exactly what a real gust does not look like.
    let tau = 3.5;
    let decay = (-dt / tau).exp();
    let front = self_gust_position(f);
    for (i, tree) in f.trees.iter().enumerate() {
        let mut stream = Stream::at(f.world.tree.world_seed, i as u128, 0, Purpose::ThermalNoise)
            .split((f.elapsed * 1000.0) as u64);
        let intensity = f.turbulence * f.wind;
        let along = tree.at.dot(f.wind_dir);
        // Trees the gust front has already passed are inside it.
        let exposure = 0.55 + 0.45 * ((front - along) * 0.08).sin();
        let target = intensity * exposure;
        f.gust[i] = f.gust[i] * decay
            + target * (1.0 - decay * decay).sqrt() * stream.normal();
    }

    if !f.dynamic || f.wind <= 0.0 {
        if !f.dynamic {
            return;
        }
    }

    let mut worst_sway = 0.0f32;
    let mut worst_util = 0.0f32;
    let mut broken = 0usize;
    for i in 0..f.trees.len() {
        let node = f.trees[i].node;
        let speed = (f.wind + f.gust[i]).max(0.0);
        let mut mechanisms: Vec<Mechanism> = Vec::new();
        if speed > 0.0 {
            mechanisms.push(st::weather::wind(speed, f.wind_dir));
        }
        let out = f.world.shake(node, &mechanisms, dt);
        broken += out.broken_joints;
        f.shed += out.detached_mass as f32 as f64;
        f.trees[i].sway = out.displacement as f32;
        worst_sway = worst_sway.max(out.displacement as f32);
        worst_util = worst_util.max(f.trees[i].utilisation);
    }
    // Whatever came away is falling now, and what it lands on has to answer
    // for it. This is where a stand becomes more than twenty independent
    // trees: a limb off one of them arrives on the branches below at a speed
    // the ordinary stress calculation has an opinion about.
    let fall = f.world.drop_fragments(dt);
    f.total_broken += (broken + fall.secondary_breaks) as u32;
    f.shed += fall.secondary_mass;
    f.readouts[20] = fall.still_falling as f32;
    f.readouts[21] = fall.struck_members as f32;
    f.readouts[22] += fall.secondary_breaks as f32;
    if broken > 0 || fall.secondary_breaks > 0 {
        refresh_forest(f);
    } else {
        write_forest_geometry(f);
    }
    write_debris(f);
    f.readouts[16] = (f.wind + f.gust.iter().sum::<f64>() / f.gust.len().max(1) as f64) as f32;
    f.readouts[17] = worst_sway;
    f.readouts[5] = worst_util;
    f.readouts[6] = f.total_broken as f32;
    f.readouts[7] = f.shed as f32;
}

/// Where the gust front has travelled to, in metres along the wind.
fn self_gust_position(f: &Forest) -> f64 {
    // A front moves downwind at roughly the mean speed.
    f.elapsed * f.wind.max(1.0)
}

/// Rebuild the geometry and the per-tree readouts from scratch.
fn refresh_forest(f: &mut Forest) {
    f.geometry.clear();
    let mut worst_util = 0.0f32;
    let mut standing = 0.0f64;
    for i in 0..f.trees.len() {
        let node = f.trees[i].node;
        let at = f.trees[i].at;
        let bodies = f.world.tree.refine(node).to_vec();
        let topo = match f.world.tree.nodes[node.get()].topology.clone() {
            Some(t) => t,
            None => {
                f.trees[i].drawn = 0;
                continue;
            }
        };
        let ambient = f.world.tree.nodes[node.get()].agg.temperature;
        let mut field = LoadField::new(bodies.len(), ambient);
        let speed = (f.wind + f.gust.get(i).copied().unwrap_or(0.0)).max(0.0);
        if speed > 0.0 {
            field.apply(&st::weather::wind(speed, f.wind_dir), &bodies, &topo);
        }
        field.apply(&st::weather::gravity(), &bodies, &topo);
        let loads = st::analyse(&bodies, &topo, &field);

        let start = f.geometry.len();
        let mut peak = 0.0f32;
        for m in 0..bodies.len().min(topo.bonds.len()) {
            if topo.bonds[m].radius <= 0.0 {
                continue;
            }
            let u = loads.get(m).map(|l| l.utilisation).unwrap_or(0.0) as f32;
            peak = peak.max(u);
            let (b, t) = (topo.base[m] + at, topo.tip[m] + at);
            f.geometry.extend_from_slice(&[
                b.x as f32,
                b.y as f32,
                b.z as f32,
                t.x as f32,
                t.y as f32,
                t.z as f32,
                topo.bonds[m].radius as f32,
                u,
            ]);
        }
        f.trees[i].drawn = (f.geometry.len() - start) / 8;
        f.trees[i].utilisation = peak;
        worst_util = worst_util.max(peak);
        let m = f.world.tree.nodes[node.get()].morphology.as_ref();
        f.trees[i].height = m.map(|m| m.height()).unwrap_or(0.0) as f32;
        standing += m.map(|m| m.built).unwrap_or(0.0);
    }
    f.standing_floats = f.geometry.len();
    f.readouts[0] = standing as f32;
    f.readouts[1] = f
        .trees
        .iter()
        .map(|t| t.height)
        .fold(0.0f32, f32::max);
    f.readouts[2] = f.trees
        .first()
        .and_then(|t| f.world.tree.nodes[t.node.get()].morphology.as_ref())
        .map(|m| (m.age / YEAR) as f32)
        .unwrap_or(0.0);
    f.readouts[3] = (f.geometry.len() / 8) as f32;
    f.readouts[5] = worst_util;
    f.readouts[18] = f.trees.len() as f32;
    f.readouts[19] = f.budget as f32;
}

/// Append the falling debris to the geometry buffer.
///
/// Debris is drawn at full utilisation, which is not a stress reading — it has
/// none, it is not carrying anything — but it is the honest colour for it: the
/// thing an observer needs to see is that this wood is no longer part of a
/// structure.
fn write_debris(f: &mut Forest) {
    f.geometry.truncate(f.standing_floats);
    for (node, frag) in f.world.falling() {
        let Some(at) = f.trees.iter().find(|t| t.node == *node).map(|t| t.at) else {
            continue;
        };
        for (a, b) in frag.dynamics.deformed_members() {
            if a == b {
                continue;
            }
            let (a, b) = (a + at, b + at);
            f.geometry.extend_from_slice(&[
                a.x as f32,
                a.y as f32,
                a.z as f32,
                b.x as f32,
                b.y as f32,
                b.z as f32,
                0.03,
                1.0,
            ]);
        }
    }
}

/// Overwrite the geometry with where every tree actually is.
fn write_forest_geometry(f: &mut Forest) {
    let mut slot = 0usize;
    for i in 0..f.trees.len() {
        let node = f.trees[i].node;
        let at = f.trees[i].at;
        let drawn = f.trees[i].drawn;
        let Some(ds) = f.world.shaken(node) else {
            slot += drawn;
            continue;
        };
        let members = ds.deformed_members();
        let topo = match f.world.tree.nodes[node.get()].topology.as_ref() {
            Some(t) => t,
            None => {
                slot += drawn;
                continue;
            }
        };
        let mut written = 0usize;
        for m in 0..members.len().min(topo.bonds.len()) {
            if topo.bonds[m].radius <= 0.0 {
                continue;
            }
            let base = (slot + written) * 8;
            written += 1;
            if base + 5 >= f.geometry.len() {
                break;
            }
            let (b, t) = members[m];
            let (b, t) = (b + at, t + at);
            f.geometry[base] = b.x as f32;
            f.geometry[base + 1] = b.y as f32;
            f.geometry[base + 2] = b.z as f32;
            f.geometry[base + 3] = t.x as f32;
            f.geometry[base + 4] = t.y as f32;
            f.geometry[base + 5] = t.z as f32;
        }
        slot += drawn;
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn forest_geometry_ptr() -> *const f32 {
    forest().geometry.as_ptr()
}

#[unsafe(no_mangle)]
pub extern "C" fn forest_geometry_len() -> u32 {
    (forest().geometry.len() / 8) as u32
}

#[unsafe(no_mangle)]
pub extern "C" fn forest_readouts_ptr() -> *const f32 {
    forest().readouts.as_ptr()
}

/// Per-tree state, for a table: height, sway, utilisation, members.
#[unsafe(no_mangle)]
pub extern "C" fn forest_tree_stat(index: u32, which: u32) -> f32 {
    let f = forest();
    let Some(t) = f.trees.get(index as usize) else {
        return 0.0;
    };
    match which {
        0 => t.height,
        1 => t.sway,
        2 => t.utilisation,
        3 => t.drawn as f32,
        4 => t.at.x as f32,
        5 => t.at.y as f32,
        _ => 0.0,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn forest_count() -> u32 {
    forest().trees.len() as u32
}

// ---------------------------------------------------------------------------
// The explorer: any scenario, any scale
// ---------------------------------------------------------------------------

/// A world being looked at, at whatever scale the observer has descended to.
///
/// # Why one viewer works at every scale
///
/// Because nothing here knows what scale it is at. Positions are emitted in
/// units of the node's own radius, so the numbers a renderer sees are of order
/// one whether the node is fifteen kiloparsecs across or four femtometres. That
/// is not a convenience for the viewer: it is the engine's own argument about
/// precision, which is that there is no global coordinate system and a body's
/// position only ever means something relative to the node holding it.
///
/// The solver, the timestep and the refinement policy are all chosen from the
/// tier, which the tree already knows. So descending is one call, and the frame
/// after it is drawn by exactly the same code as the frame before.
pub struct Explorer {
    world: World,
    /// Root first, current node last. Ascending is a pop.
    path: Vec<NodeIdx>,
    /// Interleaved: x, y, z (in units of the node radius), radius, mass,
    /// temperature, speed, kind.
    points: Vec<f32>,
    readouts: [f32; 32],
    scenario: usize,
    /// The seed the scenario was built from, so the viewer can offer another
    /// world without the caller having to remember what it asked for.
    seed: u32,
    budget: usize,
    /// Simulated seconds elapsed at the current node.
    elapsed: f64,
    dirty: bool,
}

/// Floats per body in the point buffer.
const POINT_STRIDE: usize = 8;

static mut EXPLORER: Option<Explorer> = None;

fn explorer() -> &'static mut Explorer {
    unsafe {
        let e = &raw mut EXPLORER;
        (*e).as_mut().expect("scene_create() first")
    }
}

/// How many scenarios are on the shelf.
#[unsafe(no_mangle)]
pub extern "C" fn scene_count() -> u32 {
    crate::scenario::ALL.len() as u32
}

/// Load one. Everything below it is generated on demand from here.
#[unsafe(no_mangle)]
pub extern "C" fn scene_create(index: u32, seed: u32, budget: u32) {
    let s = crate::scenario::by_index(index as usize);
    let mut world = World::new(s.build(seed as u64 | 0x5EED_0000), 20.0);
    let root = world.tree.root;
    if s.tier < Tier::Atomic {
        world.tree.nodes[root.get()].spec.count = (budget as usize).clamp(64, 40_000);
    }
    unsafe {
        let e = &raw mut EXPLORER;
        *e = Some(Explorer {
            world,
            path: vec![root],
            points: Vec::new(),
            readouts: [0.0; 32],
            scenario: index as usize,
            seed,
            budget: budget as usize,
            elapsed: 0.0,
            dirty: true,
        });
    }
    refresh_scene(explorer());
}

/// Descend into one of the bodies on show.
///
/// Promoting is what makes a body into a node: the statistical stand-in becomes
/// a thing with its own contents, generated from its own address. Nothing about
/// the body was stored — the child is built from the aggregate the body
/// represents, at the refinement policy its new tier calls for.
#[unsafe(no_mangle)]
pub extern "C" fn scene_descend(index: u32) -> u32 {
    let e = explorer();
    let Some(&here) = e.path.last() else {
        return 0;
    };
    let tier = e.world.tree.nodes[here.get()].tier;
    if tier == Tier::Nuclear {
        return 0;
    }
    e.world.tree.refine(here);
    if index as usize >= e.world.tree.nodes[here.get()].bodies.len() {
        return 0;
    }
    let spec = crate::engine::budgeted_spec(tier.finer(), e.budget);
    let child = e.world.tree.promote(here, index as usize, spec);
    if child.is_none() {
        return 0;
    }
    e.world.tree.nodes[child.get()].residency = crate::tree::Residency::Observed;
    e.path.push(child);
    e.elapsed = 0.0;
    e.dirty = true;
    refresh_scene(e);
    1
}

/// Descend into whichever body is largest, which is usually the one worth
/// looking at and is always a defined choice.
#[unsafe(no_mangle)]
pub extern "C" fn scene_descend_largest() -> u32 {
    let e = explorer();
    let Some(&here) = e.path.last() else {
        return 0;
    };
    e.world.tree.refine(here);
    let best = {
        let n = &e.world.tree.nodes[here.get()];
        let mut bi = 0usize;
        let mut bm = f64::NEG_INFINITY;
        for (i, b) in n.bodies.iter().enumerate() {
            if b.mass > bm {
                bm = b.mass;
                bi = i;
            }
        }
        bi
    };
    scene_descend(best as u32)
}

/// Back out to the node above.
#[unsafe(no_mangle)]
pub extern "C" fn scene_ascend() -> u32 {
    let e = explorer();
    if e.path.len() <= 1 {
        return 0;
    }
    let leaving = e.path.pop().unwrap();
    // Coarsening on the way out is the point of the whole architecture: the
    // detail is thrown away, and going back in regenerates it.
    e.world.tree.coarsen(leaving);
    e.elapsed = 0.0;
    e.dirty = true;
    refresh_scene(e);
    1
}

/// Advance the node being watched by `seconds` of *its own* time.
///
/// Each tier runs at its own rate — a galaxy steps in megayears and a nucleus in
/// zeptoseconds — so the argument is a number of the node's own timesteps
/// rather than a wall-clock interval, which would be meaningless across
/// forty-five orders of magnitude.
#[unsafe(no_mangle)]
pub extern "C" fn scene_step(steps: f32) {
    let e = explorer();
    let Some(&here) = e.path.last() else {
        return;
    };
    if steps <= 0.0 {
        return;
    }
    let dt = e.world.node_dt(here);
    if !(dt > 0.0) || !dt.is_finite() {
        return;
    }
    // Whole steps, one at a time. Handing the solver eight times its own
    // timestep is not the same as eight steps: every stability limit in the
    // engine is a limit on *dt*, and the SPH solver answers a 8 ms step at the
    // continuum tier by heating its parcels to two thirds of light speed.
    let whole = (steps as f64).clamp(0.0, 8.0);
    let count = whole.floor() as u32;
    let mut drift = 0.0f64;
    let mut momentum = 0.0f64;
    for _ in 0..count.max(1) {
        let report = e.world.advance_node(here, dt);
        drift = drift.max(report.drift());
        momentum = momentum.max(report.momentum_drift());
        e.elapsed += dt;
    }
    e.readouts[20] = drift as f32;
    e.readouts[21] = momentum as f32;
    e.dirty = true;
    refresh_scene(e);
}

/// Throw away the detail and rebuild it, which is the claim worth watching.
#[unsafe(no_mangle)]
pub extern "C" fn scene_regenerate() -> u32 {
    let e = explorer();
    let Some(&here) = e.path.last() else {
        return 0;
    };
    let before = e.world.tree.nodes[here.get()].bodies.len() as u32;
    e.world.tree.coarsen(here);
    e.dirty = true;
    refresh_scene(e);
    before
}

/// Rebuild the same scenario from a different seed — another galaxy of the
/// same kind, another cloud, another nucleus.
#[unsafe(no_mangle)]
pub extern "C" fn scene_reseed(step: u32) {
    let e = explorer();
    let (scenario, budget) = (e.scenario as u32, e.budget as u32);
    let seed = e.seed.wrapping_add(step.max(1));
    scene_create(scenario, seed, budget);
}

/// Set how many bodies a node materialises into.
#[unsafe(no_mangle)]
pub extern "C" fn scene_set_budget(budget: u32) {
    let e = explorer();
    e.budget = (budget as usize).clamp(64, 40_000);
    let Some(&here) = e.path.last() else {
        return;
    };
    e.world.tree.nodes[here.get()].spec.count = e.budget;
    e.world.tree.coarsen(here);
    e.dirty = true;
    refresh_scene(e);
}

fn refresh_scene(e: &mut Explorer) {
    if !e.dirty {
        return;
    }
    e.dirty = false;
    let Some(&here) = e.path.last() else {
        return;
    };
    let bodies = e.world.tree.refine(here).to_vec();
    let node_radius = e.world.tree.nodes[here.get()].agg.radius.max(1e-300);
    let inv = 1.0 / node_radius;

    e.points.clear();
    let mut fastest = 0.0f64;
    let mut hottest = 0.0f64;
    let mut heaviest = 0.0f64;
    for b in &bodies {
        let speed = b.vel.norm();
        fastest = fastest.max(speed);
        hottest = hottest.max(b.temperature);
        heaviest = heaviest.max(b.mass);
        e.points.extend_from_slice(&[
            (b.pos.x * inv) as f32,
            (b.pos.y * inv) as f32,
            (b.pos.z * inv) as f32,
            (b.radius * inv) as f32,
            b.mass as f32,
            b.temperature as f32,
            speed as f32,
            b.kind as u32 as f32,
        ]);
    }

    let n = &e.world.tree.nodes[here.get()];
    e.readouts[0] = n.agg.mass as f32;
    e.readouts[1] = node_radius as f32;
    e.readouts[2] = n.agg.temperature as f32;
    e.readouts[3] = bodies.len() as f32;
    e.readouts[4] = n.tier as u32 as f32;
    e.readouts[5] = e.path.len() as f32 - 1.0;
    e.readouts[6] = heaviest as f32;
    e.readouts[7] = fastest as f32;
    e.readouts[8] = hottest as f32;
    e.readouts[9] = e.world.node_dt(here) as f32;
    e.readouts[10] = e.elapsed as f32;
    e.readouts[11] = e.world.tree.detail_bytes() as f32;
    e.readouts[12] = n.agg.internal_energy as f32;
    e.readouts[13] = n.agg.binding_energy as f32;
    e.readouts[14] = n.agg.luminosity as f32;
    e.readouts[15] = e.scenario as f32;
    e.readouts[16] = e.budget as f32;
    e.readouts[17] = crate::solvers::for_tier(n.tier) as u32 as f32;
    e.readouts[18] = n.agg.charge as f32;
    e.readouts[19] = e.world.tree.nodes.iter().filter(|n| n.alive && n.is_materialised()).count() as f32;
}

/// One quantity about the current node, in full precision.
///
/// A galaxy weighs 10^39 kg and a nucleon 10^-27, which single precision cannot
/// hold at both ends: the readout array is `f32` for the bulk of the numbers a
/// renderer needs, and a galaxy's mass overflows it to infinity and prints as
/// zero. Anything spanning the whole ladder comes back through here instead.
#[unsafe(no_mangle)]
pub extern "C" fn scene_value(which: u32) -> f64 {
    let e = explorer();
    let Some(&here) = e.path.last() else {
        return 0.0;
    };
    let n = &e.world.tree.nodes[here.get()];
    match which {
        0 => n.agg.mass,
        1 => n.agg.radius,
        2 => n.agg.temperature,
        3 => e.world.node_dt(here),
        4 => e.elapsed,
        5 => n.agg.internal_energy,
        6 => n.agg.binding_energy,
        7 => n.agg.luminosity,
        8 => e.points
            .chunks_exact(POINT_STRIDE)
            .map(|p| p[6] as f64)
            .fold(0.0f64, f64::max),
        9 => e.world.tree.detail_bytes() as f64,
        _ => 0.0,
    }
}

/// Radius in metres of the node `up` levels above the current one, in full
/// precision. See [`scene_value`].
#[unsafe(no_mangle)]
pub extern "C" fn scene_trail_metres(up: u32) -> f64 {
    let e = explorer();
    let len = e.path.len();
    if (up as usize) >= len {
        return 0.0;
    }
    let idx = e.path[len - 1 - up as usize];
    e.world.tree.nodes[idx.get()].agg.radius
}

#[unsafe(no_mangle)]
pub extern "C" fn scene_points_ptr() -> *const f32 {
    explorer().points.as_ptr()
}

#[unsafe(no_mangle)]
pub extern "C" fn scene_points_len() -> u32 {
    (explorer().points.len() / POINT_STRIDE) as u32
}

#[unsafe(no_mangle)]
pub extern "C" fn scene_readouts_ptr() -> *const f32 {
    explorer().readouts.as_ptr()
}

/// Depth of the current node below the scenario's root.
#[unsafe(no_mangle)]
pub extern "C" fn scene_depth() -> u32 {
    explorer().path.len() as u32 - 1
}

/// Tier of the node `up` levels above the current one, as an index. Returns
/// 255 past the top, which is how the viewer knows where the trail ends.
#[unsafe(no_mangle)]
pub extern "C" fn scene_trail_tier(up: u32) -> u32 {
    let e = explorer();
    let len = e.path.len();
    if (up as usize) >= len {
        return 255;
    }
    let idx = e.path[len - 1 - up as usize];
    e.world.tree.nodes[idx.get()].tier as u32
}

/// Radius in metres of the node `up` levels above the current one.
#[unsafe(no_mangle)]
pub extern "C" fn scene_trail_scale(up: u32) -> f32 {
    let e = explorer();
    let len = e.path.len();
    if (up as usize) >= len {
        return 0.0;
    }
    let idx = e.path[len - 1 - up as usize];
    e.world.tree.nodes[idx.get()].agg.radius as f32
}
