//! The frame-rate guarantee, and how it degrades.

use phys::budget::*;
use phys::engine::{default_spec, galaxy, World};
use phys::ids::NodeIdx;
use phys::math::v3;
use phys::observe::Observer;
use phys::units::*;

fn task(node: u32, cost: f64, value: f64, bytes: i64) -> Task {
    Task {
        node: NodeIdx(node),
        kind: TaskKind::Step,
        cost_us: cost,
        salience: value,
        urgency: 1.0,
        error: 1.0,
        novelty: 0.0,
        bytes,
    }
}

/// The planner takes the most valuable work per microsecond first.
#[test]
fn planner_prefers_value_density() {
    let b = FrameBudget::ups(20.0);
    let tasks = vec![
        task(0, 1000.0, 1.0, 0),   // density 0.001
        task(1, 10.0, 5.0, 0),     // density 0.5   <- best
        task(2, 100.0, 10.0, 0),   // density 0.1
    ];
    let plan = b.plan(tasks, 0);
    assert_eq!(plan.accepted[0].node, NodeIdx(1));
    assert_eq!(plan.accepted[1].node, NodeIdx(2));
    assert_eq!(plan.accepted[2].node, NodeIdx(0));
}

/// Planning is deterministic, including for ties — otherwise replay diverges.
#[test]
fn planning_is_deterministic() {
    let b = FrameBudget::ups(20.0);
    let make = || (0..200).map(|i| task(i, 100.0, 1.0, 0)).collect::<Vec<_>>();
    let a: Vec<u32> = b.plan(make(), 0).accepted.iter().map(|t| t.node.0).collect();
    let c: Vec<u32> = b.plan(make(), 0).accepted.iter().map(|t| t.node.0).collect();
    assert_eq!(a, c, "planner must be order-stable");
}

/// Work beyond the budget is deferred, not run. The frame time is the
/// invariant.
#[test]
fn excess_work_is_deferred() {
    let b = FrameBudget::ups(20.0);
    let budget = b.sim_budget_us();
    let tasks: Vec<Task> = (0..100).map(|i| task(i, budget / 10.0, 1.0, 0)).collect();
    let plan = b.plan(tasks, 0);
    assert!(plan.planned_us <= budget * 1.001, "planned {:.0} us over budget {budget:.0}", plan.planned_us);
    assert_eq!(plan.accepted.len(), 10);
    assert_eq!(plan.deferred, 90);
    assert!(plan.unmet_value > 0.0, "unserved demand must be reported");
}

/// Coarsening frees resources, so it is never deferred for cost.
#[test]
fn freeing_work_is_always_accepted() {
    let b = FrameBudget::ups(20.0);
    let mut tasks: Vec<Task> = (0..100).map(|i| task(i, b.sim_budget_us(), 1.0, 1 << 20)).collect();
    tasks.push(Task {
        kind: TaskKind::Coarsen,
        ..task(999, 1.0, 0.001, -(1 << 20))
    });
    let plan = b.plan(tasks, 0);
    assert!(
        plan.accepted.iter().any(|t| t.node == NodeIdx(999)),
        "a task that frees memory must not be starved"
    );
}

/// The memory cap is respected independently of the time budget: on a 6 GB
/// card, memory runs out first.
#[test]
fn memory_cap_is_enforced() {
    let mut b = FrameBudget::ups(20.0);
    b.memory_cap = 10 << 20; // 10 MB
    let tasks: Vec<Task> = (0..100).map(|i| task(i, 1.0, 1.0, 1 << 20)).collect();
    let plan = b.plan(tasks, 0);
    assert!(plan.planned_bytes <= 10 << 20, "exceeded memory cap");
    assert!(plan.deferred >= 90);
}

/// The cost model must converge on a machine whose speed it did not know.
#[test]
fn calibration_converges() {
    let mut b = FrameBudget::ups(20.0);
    // Pretend everything actually costs 7x the estimate.
    for _ in 0..40 {
        let planned = 10_000.0;
        b.observe_frame(planned, planned * 7.0 / b.calibration());
    }
    let c = b.calibration();
    println!("calibration converged to {c:.3} (true factor 7)");
    assert!((c - 7.0).abs() / 7.0 < 0.25, "calibration {c:.3} did not converge to 7");
}

/// End to end: whatever a user asks for, the frame budget holds.
#[test]
fn frames_stay_within_budget() {
    let mut w = World::new(galaxy(0x1234, 1e9), 20.0);
    let root = w.tree.root;
    w.add_observer(Observer {
        anchor: root,
        offset: v3(8.0 * KPC, 0.0, 0.0),
        look: v3(-1.0, 0.0, 0.0),
        field: 3.2,
        angular_resolution: 1e-9, // an absurdly demanding observer
        horizon: 1e4 * YEAR,
        priority: 1000.0,
        ..Default::default()
    });
    w.gate = phys::causal::CausalGate::new(1e4 * YEAR);
    let path = w.drill(root, Tier::Molecular, &default_spec);
    for &n in &path {
        w.tree.pin(n);
    }
    let target_us = 50_000.0;
    let mut worst: f64 = 0.0;
    let mut debt_seen = false;
    // The first frames pay materialisation costs the model has not calibrated
    // for yet; measure the steady state.
    for i in 0..25 {
        let plan = w.step_frame(target_us);
        if i >= 5 {
            worst = worst.max(w.stats.last_frame_us);
        }
        debt_seen |= plan.deferred > 0;
    }
    println!(
        "steady-state worst frame {:.1} ms (target {:.0} ms), calibration {:.2}",
        worst / 1e3,
        target_us / 1e3,
        w.budget.calibration()
    );
    assert!(debt_seen, "an over-demanding observer should produce detail debt");
    assert!(
        worst < target_us * 1.6,
        "frame overran: {:.1} ms vs {:.0} ms target",
        worst / 1e3,
        target_us / 1e3
    );
    // And it must still be doing useful work, not just refusing everything.
    assert!(w.stats.bodies_stepped > 0, "engine did no work at all");
}
