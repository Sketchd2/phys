//! A world runs at one second per second.
//!
//! `docs/PLAY.md` D1: that is what a shared world *is*, and it retires the
//! engine's most quotable property. Pacing the clock to whatever is being
//! watched — a galaxy advancing millennia per frame, a carbon atom advancing
//! femtoseconds — is a fine answer to a single observer's question and the
//! wrong answer to a shared one, because a player inspecting a rifle bolt must
//! not slow down the war. `pace_to` survives as a tool for single-player
//! exploration and offline study.
//!
//! Worth saying what actually moved, because the field itself did not: a
//! freshly built `World` already had `pace == 1.0`. What changed is that
//! `World::new` no longer ends with an implicit `pace_to(root)`, and that
//! `time_throttle` — the *second* mechanism that slows a world's clock under
//! load — now applies only in `Follow`. Half a fix would have been a world that
//! claims one second per second and still slows down when it is busy.

use phys::engine::{default_spec, galaxy, PaceMode, World};
use phys::units::{Tier, YEAR};

#[test]
fn a_fresh_world_runs_at_one_second_per_second() {
    let mut w = World::new(galaxy(0x15EC, 1e9), 20.0);
    assert_eq!(w.pace_mode, PaceMode::Fixed, "a world is a fixed-pace world");
    assert_eq!(w.pace, 1.0);

    let before = w.time;
    w.step_frame(50_000.0);
    let span = w.time - before;
    assert!(
        (span - 1.0).abs() < 1e-9,
        "one frame of a galaxy advanced {span} s of world time, not 1 s"
    );
}

/// Materialising what is being watched does not drag the clock.
///
/// This is the property D1 retires, stated so that reinstating the old default
/// is caught rather than merely noticed. Under the old behaviour, refining a
/// galaxy into its stars shortened its characteristic time and the pace came
/// down with it — by two orders of magnitude, which `README.md` used to cite as
/// the engine's headline.
#[test]
fn refining_what_is_watched_does_not_change_the_clock() {
    let mut w = World::new(galaxy(0x15ED, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 256;
    let root = w.tree.root;

    let before = w.time;
    w.step_frame(50_000.0);
    let coarse = w.time - before;

    w.tree.refine(root);
    let tier = w.tree.nodes[root.get()].tier;
    assert!(!w.tree.promote(root, 0, default_spec(tier.finer())).is_none());

    let before = w.time;
    w.step_frame(50_000.0);
    let fine = w.time - before;

    assert!(
        (fine - coarse).abs() < 1e-9,
        "a frame covered {coarse} s before the galaxy was materialised and \
         {fine} s after"
    );
}

/// Overload makes a fixed-pace world *staler*, not slower.
///
/// The clock-slowing mechanism is `time_throttle`, and leaving it applying
/// everywhere would have been the misleading half of this change: the pace
/// would read 1.0 while the frame advanced a twentieth of a second. Measured
/// before it was scoped to `Follow` — a world drilled to the nucleus on a
/// 50 ms budget advanced 0.05 s per frame.
#[test]
fn a_starved_world_keeps_its_clock_and_falls_behind_instead() {
    let mut w = World::new(galaxy(0x15EF, 1e9), 20.0);
    let root = w.tree.root;
    let path = w.drill_to(root, Tier::Nuclear.max_radius(), &default_spec);
    assert!(path.len() >= 6, "expected a deep ladder, got {}", path.len());

    for _ in 0..10 {
        let before = w.time;
        // A budget far too small for the work, which is the point.
        w.step_frame(200.0);
        let span = w.time - before;
        assert!(
            (span - 1.0).abs() < 1e-9,
            "a starved frame advanced {span} s instead of holding the clock at 1 s"
        );
    }
    assert_eq!(w.time_throttle, 1.0, "the throttle must not move in Fixed");
    // What gives way instead. The engine is behind, and says so rather than
    // hiding it in the clock.
    assert!(
        w.stats.worst_lateness > 0.0,
        "a starved world should be visibly late, and reports {}",
        w.stats.worst_lateness
    );
}

/// `pace_to` still works, because it is still wanted — just asked for. And the
/// throttle goes with it: a galactic pace is exactly the case the throttle's
/// own reasoning is about, since nothing's trajectory can be integrated across
/// a frame of millennia.
#[test]
fn pacing_to_a_subject_is_available_on_request() {
    let mut w = World::new(galaxy(0x15EE, 1e9), 20.0);
    let root = w.tree.root;
    w.pace_to(root);
    assert_eq!(w.pace_mode, PaceMode::Follow);
    assert!(
        w.pace > 1.0e3 * YEAR,
        "paced to a galaxy a frame should cover an astronomical span, and \
         covers {} s",
        w.pace
    );

    let before = w.time;
    w.step_frame(50_000.0);
    assert!(
        w.time - before > 1.0e3 * YEAR,
        "the clock did not follow the subject it was paced to"
    );

    // And a world can be handed back its own clock.
    w.pace_fixed(1.0);
    assert_eq!(w.pace_mode, PaceMode::Fixed);
    let before = w.time;
    w.step_frame(50_000.0);
    assert!((w.time - before - 1.0).abs() < 1e-9);
}
