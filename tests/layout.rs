//! A layout is derived, not selected.
//!
//! `docs/PLAY.md` D11's complaint about `morph::Program` is that it is a
//! species table: a variant per kind of thing and a column per property. The
//! owner's call on Phase 4 went past removing the columns — *"Program should be
//! a function that takes in the current and historical data, and decides on how
//! it is laid out. The program should be derived using the axioms, not
//! engineered."*
//!
//! So this is the engine deriving one, under its own power, for the case the
//! owner named: molten rock cooling. Nothing here plants anything, names a
//! program or states a shape. A ball of melt is put somewhere dark and left to
//! radiate, and what it is laid out as when it has frozen is a consequence of
//! the physics it went through.

use phys::chem::{Mixture, Phase};
use phys::engine::World;
use phys::ids::NodeIdx;
use phys::recipe::Recipe;
use phys::sampler::{MassSpectrum, Profile, SampleSpec};
use phys::state::{BodyKind, Composition, Matter};
use phys::tree::Tree;
use phys::units::*;

/// A ball of silicate of this radius, at this temperature, alone in the dark.
///
/// The scenario states four things — a size, a temperature, a mass consistent
/// with the two, and what it is made of — and none of them is a layout. Whether
/// the silica is stated liquid or solid is the *only* difference between the
/// two cases below, and it is a statement about the matter, not about the
/// shape.
fn a_ball_of_silicate(seed: u64, radius: f64, temperature: f64, phase: Phase) -> (World, NodeIdx) {
    let spec = SampleSpec::new(64, Profile::Uniform, MassSpectrum::Equal, BodyKind::Grain);
    // Mass from the substance's own density, so the node's bulk density is the
    // density of what it is made of and nothing is being smuggled in.
    let probe = phys::chem::analyse(&phys::material::substances::silica_arrangement())
        .expect("silica analyses");
    let mass = probe.density * 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    let matter = Matter::neutral(mass, radius, temperature, Composition::primordial());
    let mut w = World::new(Tree::new(seed, matter, Tier::Continuum, spec), 20.0);
    w.pace_fixed(60.0);
    let rock = w.tree.root;
    let silica = w
        .substances
        .intern(phys::material::substances::silica_arrangement())
        .expect("silica interns");
    let mut mix = Mixture::new();
    mix.add(silica, phase, 1.0);
    w.set_mixture(rock, mix);
    (w, rock)
}

/// Silica's melting point, derived like everything else about it.
fn silica_melts_at() -> f64 {
    phys::chem::analyse(&phys::material::substances::silica_arrangement())
        .expect("silica analyses")
        .melting_point
}

/// Run until the node has a layout, or give up. Returns the frame it got one.
fn run_until_laid_out(w: &mut World, idx: NodeIdx, frames: usize) -> Option<usize> {
    for f in 0..frames {
        w.step_frame(20_000.0);
        if w.tree.nodes[idx.get()].morphology.is_some() {
            return Some(f);
        }
    }
    None
}

/// Molten rock, left to cool, lays itself out as grains — and says how big.
#[test]
fn the_engine_derives_a_program_for_molten_rock_cooling() {
    let melts = silica_melts_at();
    let (mut w, rock) = a_ball_of_silicate(0x1A7A, 1.0, 2200.0, Phase::Liquid);
    assert!(
        w.tree.nodes[rock.get()].morphology.is_none(),
        "it starts as a ball of melt and nothing else"
    );
    assert!(2200.0 > melts, "the melt has to start molten: silica melts at {melts:.0} K");

    // Nobody is watching it. That is deliberate: a layout that appeared only
    // for an observed node would be a fact about the observer.
    assert!(w.observers.is_empty());

    let frame = run_until_laid_out(&mut w, rock, 400).expect("the melt never laid itself out");
    let n = &w.tree.nodes[rock.get()];
    let Some(Recipe::Granular(g)) = n.morphology.as_ref().and_then(|m| m.recipe.as_ref()) else {
        panic!("a freezing melt should lay down grains, not {:?}",
            n.morphology.as_ref().and_then(|m| m.recipe.as_ref()).map(|r| r.habit()));
    };
    println!(
        "  a ball of melt at 2200 K froze through {melts:.0} K at frame {frame}, t = {:.0} s, \
         and laid down {:.3e} m grains — {:.2e} of them across {:.3} m",
        w.time,
        g.grain,
        g.grains_across(),
        g.side
    );

    assert!(
        n.matter.temperature < melts,
        "it should be below its melting point: {:.1} K",
        n.matter.temperature
    );
    assert!(
        n.matter.mixture.in_phase(Phase::Solid) > 0.5,
        "and mostly solid: {:.3}",
        n.matter.mixture.in_phase(Phase::Solid)
    );
    assert!(g.grain > 0.0 && g.grain.is_finite(), "the grain is a length: {:.3e}", g.grain);
    assert!(
        g.grain < g.side,
        "and a grain is smaller than the thing it is a grain of: {:.3e} against {:.3}",
        g.grain,
        g.side
    );

    // The rule draws bodies, and they are the grains at whatever resolution is
    // asked for rather than one body per grain.
    w.tree.refine(rock);
    let drawn = w.tree.nodes[rock.get()].bodies.len();
    println!("  and it draws {drawn} cells of it at this node's budget");
    assert!(drawn > 0, "the derived rule draws nothing");
}

/// The gate is the freezing, not the being solid.
///
/// Without this the first test proves nothing: a derivation that fired on "is
/// this node mostly solid and below its melting point" would hand a granular
/// layout to every rock in the world, including one that was never molten. Run
/// against that version, this fails on the first frame.
#[test]
fn rock_that_was_never_molten_gets_no_layout() {
    let melts = silica_melts_at();
    let (mut w, rock) = a_ball_of_silicate(0x1A7B, 1.0, 290.0, Phase::Solid);
    assert!(290.0 < melts, "it has to start cold");
    let laid = run_until_laid_out(&mut w, rock, 400);
    let n = &w.tree.nodes[rock.get()];
    println!(
        "  a cold rock of the same size, same mixture, 400 frames: {:.1} K, solid {:.3}, layout {:?}",
        n.matter.temperature,
        n.matter.mixture.in_phase(Phase::Solid),
        n.morphology.as_ref().and_then(|m| m.recipe.as_ref()).map(|r| r.habit())
    );
    assert!(
        n.matter.mixture.in_phase(Phase::Solid) > 0.5,
        "it is solid throughout, which is the whole point of the control"
    );
    assert_eq!(laid, None, "a rock that was always cold has no formation to derive a layout from");
}

/// And the grain is a consequence of the cooling, not a number.
///
/// A small ball has more surface per kilogram than a big one, so it cools
/// faster through the same melting point, so nucleation wins more of the race
/// against growth and the grains come out finer. Nothing states the
/// relationship; it falls out of the same nucleation march `Formation::Frozen`
/// already ran for Phase 4's flaw law.
#[test]
fn a_faster_freeze_lays_down_a_finer_grain() {
    let mut measured: Vec<(f64, f64)> = Vec::new();
    for (seed, radius) in [(0x1A80u64, 0.25f64), (0x1A81, 1.0), (0x1A82, 4.0)] {
        let (mut w, rock) = a_ball_of_silicate(seed, radius, 2200.0, Phase::Liquid);
        run_until_laid_out(&mut w, rock, 2000).expect("a melt that never froze");
        let Some(Recipe::Granular(g)) =
            w.tree.nodes[rock.get()].morphology.as_ref().and_then(|m| m.recipe.as_ref())
        else {
            panic!("no granular layout at {radius} m");
        };
        println!("  a {radius:>5} m ball of melt freezes into {:.4e} m grains", g.grain);
        measured.push((radius, g.grain));
    }
    for pair in measured.windows(2) {
        let (r0, g0) = pair[0];
        let (r1, g1) = pair[1];
        assert!(
            g1 > g0,
            "the {r1} m ball cooled more slowly than the {r0} m one and should be coarser: \
             {g1:.4e} against {g0:.4e}"
        );
    }
}
