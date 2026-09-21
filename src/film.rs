//! Gathering a scene out of the tree, so a test can be watched rather than
//! only asserted about.
//!
//! `docs/VIEWING.md`'s fourth piece. `render.rs` knows `math`, `state` and
//! `topology` and should keep knowing nothing else: it is a rasteriser, and a
//! rasteriser that can reach a `World` will eventually solve something. So the
//! gathering — take a node's contents, walk its promoted children, place each
//! one through `Tree::separation`, flatten to one body list — lives here
//! instead, beside it, where knowing about `World` is allowed.
//!
//! It lives in the library rather than in `tests/` because every integration
//! test is its own crate and a helper in one is invisible to the rest.
//!
//! # The instrument does not disturb what it measures
//!
//! `VIEWING.md` sketched this as "refine the node, take its topology, walk its
//! promoted children". It does not refine, and the difference matters in this
//! engine more than it would in most: materialising bumps nothing but it does
//! set `last_disturbed`, allocate bodies against the frame's byte budget, and
//! change which nodes the scheduler then finds materialised. A diagnostic that
//! changes the LOD of the thing it is drawing is measuring the picture rather
//! than the world.
//!
//! So this draws what is *there*. A node holding bodies draws its bodies; a
//! node holding none draws as one disc at its own radius, which is the honest
//! picture of a node nobody has resolved and is exactly what makes a ball of
//! matter becoming a planet watchable from the first frame.

use crate::ids::NodeIdx;
use crate::math::Vec3;
use crate::render::{Canvas, Drawn, Paint, Shot};
use crate::state::Body;
use crate::topology::Topology;
use crate::World;

/// One node's contents, flattened into the frame of the node that was asked
/// about.
pub struct Scene {
    pub bodies: Vec<Body>,
    /// Present only when the whole scene came from one node that has one. A
    /// composed scene has several frames' worth of members and no single index
    /// space to put them in, so members are drawn for the frame node alone.
    pub topology: Option<Topology>,
    pub intact: Vec<bool>,
    /// How far the contents reach from the frame node's origin, metres — what
    /// a camera should frame.
    pub reach: f64,
    /// Worst accumulated round-off in any placement, metres.
    ///
    /// `Tree::separation` returns a `Bounded` and this is its `err`, carried
    /// rather than discarded: it is the honest statement of when a composed
    /// scene is below its own precision. The trap it names is real and has
    /// happened — two children 50 m apart at 2x10^20 m, drawn on top of each
    /// other.
    pub err: f64,
}

impl Scene {
    pub fn is_empty(&self) -> bool {
        self.bodies.is_empty()
    }
}

/// Everything `node` holds, in `node`'s frame.
///
/// Promoted children are placed through `Tree::separation`, which walks to the
/// lowest common ancestor and back, so a child of a child lands in the right
/// place with the right error bound. Orientation is composed too, by
/// `Tree::axes_from` — a child's contents are in the child's own axes, and
/// before Phase 4 nothing in the tree composed a rotation across more than one
/// level.
pub fn of_node(world: &World, node: NodeIdx) -> Scene {
    let mut scene = Scene {
        bodies: Vec::new(),
        topology: None,
        intact: Vec::new(),
        reach: 0.0,
        err: 0.0,
    };
    if node.is_none() || node.get() >= world.tree.nodes.len() || !world.tree.nodes[node.get()].alive
    {
        return scene;
    }

    // The frame node's own contents, in its own frame already.
    {
        let n = &world.tree.nodes[node.get()];
        scene.bodies.extend_from_slice(&n.bodies);
        scene.topology = n.topology.clone();
        if scene.bodies.is_empty() {
            // Nothing resolved: the node is one body of its own stated size.
            scene.bodies.push(Body {
                pos: Vec3::ZERO,
                radius: n.matter.radius,
                mass: n.matter.mass,
                temperature: n.matter.temperature,
                ..Default::default()
            });
            scene.topology = None;
        }
        scene.intact = n.structural_mask().unwrap_or_else(|| vec![true; scene.bodies.len()]);
    }

    // Everything promoted below, wherever it is. A promoted child's own body is
    // already in the list as a stand-in, so the child replaces it rather than
    // being added beside it — otherwise a resolved thing is drawn twice, once
    // where it is and once as the smear it left behind.
    let mut placed: Vec<(usize, Vec<Body>)> = Vec::new();
    let mut descend: Vec<NodeIdx> = world.tree.nodes[node.get()].children.clone();
    let slots: Vec<usize> = (0..descend.len()).collect();
    for (slot, child) in slots.into_iter().zip(descend.drain(..)) {
        if child.is_none() || !world.tree.nodes[child.get()].alive {
            continue;
        }
        let sub = of_node(world, child);
        let sep = world.tree.separation(node, Vec3::ZERO, child, Vec3::ZERO);
        let q = world.tree.axes_from(node, child);
        scene.err = scene.err.max(sep.err);
        let moved: Vec<Body> = sub
            .bodies
            .iter()
            .map(|b| Body { pos: sep.value + q.rotate(b.pos), ..*b })
            .collect();
        placed.push((slot, moved));
    }
    // Applied back to front so the earlier splices do not move the later ones.
    placed.sort_by_key(|(slot, _)| std::cmp::Reverse(*slot));
    for (slot, bodies) in placed {
        if slot < scene.bodies.len() {
            scene.bodies.splice(slot..slot + 1, bodies.iter().copied());
            let mask = vec![true; bodies.len()];
            if slot < scene.intact.len() {
                scene.intact.splice(slot..slot + 1, mask);
            }
            // A composed scene has members from more than one index space, and
            // there is no honest way to keep one topology across them.
            scene.topology = None;
        } else {
            scene.bodies.extend(bodies.iter().copied());
            scene.intact.extend(std::iter::repeat_n(true, bodies.len()));
        }
    }

    scene.reach = scene
        .bodies
        .iter()
        .map(|b| b.pos.norm() + b.radius)
        .fold(0.0f64, f64::max);
    scene
}

/// Draw one node onto a canvas, and say what the frame turned out to be.
pub fn shoot(world: &World, node: NodeIdx, shot: &Shot) -> (Canvas, Drawn) {
    let scene = of_node(world, node);
    let mut canvas = shot.canvas();
    let drawn = crate::render::draw(
        &mut canvas,
        &shot.camera,
        &scene.bodies,
        scene.topology.as_ref(),
        &scene.intact,
        &shot.paint,
        &shot.style,
    );
    (canvas, drawn)
}

/// A shot that frames whatever the node is currently holding.
///
/// Made **once**, at the start of a sequence, and then reused — see
/// [`crate::render::Shot`] for why that is a type rather than a convention.
pub fn framing(world: &World, node: NodeIdx, azimuth: f64, elevation: f64) -> Shot {
    let scene = of_node(world, node);
    let radius = scene
        .reach
        .max(world.tree.nodes[node.get()].matter.radius)
        .max(1e-30);
    Shot::framing(Vec3::ZERO, radius, azimuth, elevation)
}

/// Each body's liquid-and-gas fraction, measured rather than labelled.
///
/// The one quantity `Paint` cannot read off a `Body`: a body carries a
/// substance id and the registry that resolves it lives in `World`. Zero is
/// wholly solid, one is wholly fluid, and a body nobody has speciated takes its
/// node's own mixture, which is the same answer at a coarser resolution.
pub fn fluidity(world: &World, node: NodeIdx, bodies: &[Body]) -> Vec<f64> {
    use crate::chem::Phase;
    let mix = world.mixture_of(node);
    let node_fluid = mix.in_phase(Phase::Liquid) + mix.in_phase(Phase::Gas);
    bodies
        .iter()
        .map(|b| {
            if b.substance == crate::chem::SubstanceId::UNSPECIATED {
                return node_fluid;
            }
            let mut fluid = 0.0;
            let mut total = 0.0;
            for e in mix.entries() {
                if e.substance != b.substance {
                    continue;
                }
                total += e.fraction;
                if e.phase != Phase::Solid {
                    fluid += e.fraction;
                }
            }
            if total > 0.0 {
                fluid / total
            } else {
                node_fluid
            }
        })
        .collect()
}

/// A film of one node, taken at whatever cadence the caller steps the world.
///
/// The camera is made once from the first frame and never moved. Nothing
/// happens at all unless `PHYS_FILM` names a directory.
pub struct NodeFilm {
    film: Option<crate::render::Film>,
    shot: Shot,
    node: NodeIdx,
}

impl NodeFilm {
    /// Open a film of `node`, framed at `radius` metres.
    pub fn open(name: &str, node: NodeIdx, radius: f64, paint: Paint) -> NodeFilm {
        NodeFilm {
            film: crate::render::Film::open(name),
            shot: Shot::framing(Vec3::ZERO, radius.max(1e-30), 0.6, 0.35).painted(paint),
            node,
        }
    }

    /// Whether this film is actually recording.
    pub fn recording(&self) -> bool {
        self.film.is_some()
    }

    pub fn frames(&self) -> usize {
        self.film.as_ref().map(|f| f.frames()).unwrap_or(0)
    }

    /// Take a frame. Free when nothing is recording.
    pub fn take(&mut self, world: &World) -> Option<Drawn> {
        let film = self.film.as_mut()?;
        let (canvas, drawn) = shoot(world, self.node, &self.shot);
        film.shoot(&canvas);
        Some(drawn)
    }
}
