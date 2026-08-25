//! Structural dynamics against closed-form results.
//!
//! A dynamic solver has two things to get right that a static one does not: the
//! rate at which the structure moves, and the fact that a load arriving
//! suddenly is worse than the same load standing still. Both have exact
//! answers, so both are checked here rather than eyeballed in the viewer.

use phys::math::{v3, Vec3};
use phys::solvers::dynamics::Dynamics;
use phys::solvers::frame::{Dof, Frame};
use phys::topology::Material;

/// A cantilever of `n` equal elements, built in at node 0, running along +x.
fn cantilever(material: Material, length: f64, radius: f64, n: usize) -> Dynamics {
    let mut frame = Frame::new(material);
    let mut prev = frame.add_node(v3(0.0, 0.0, 0.0), true);
    for i in 1..=n {
        let x = length * i as f64 / n as f64;
        let node = frame.add_node(v3(x, 0.0, 0.0), false);
        frame.add_beam(prev, node, radius);
        prev = node;
    }
    Dynamics::new(frame)
}

fn undamped(d: &mut Dynamics) {
    d.mass_damping = 0.0;
    d.stiff_damping = 0.0;
}

/// Analytic fundamental period of a uniform cantilever.
///
/// `omega_1 = (1.875104 / L)^2 sqrt(EI / rho A)`, the first root of
/// `cos(kL) cosh(kL) = -1`.
fn analytic_period(m: Material, length: f64, radius: f64) -> f64 {
    let area = std::f64::consts::PI * radius * radius;
    let inertia = std::f64::consts::PI * radius.powi(4) / 4.0;
    let k = 1.8751040687119611 / length;
    let omega = k * k * (m.stiffness * inertia / (m.density * area)).sqrt();
    2.0 * std::f64::consts::PI / omega
}

/// The Rayleigh quotient on a static deflected shape must reproduce the
/// cantilever's fundamental period.
///
/// This is the cheap estimate the scheduler uses to decide whether a timestep
/// still resolves the motion, so it has to be right to better than the margin
/// it is asked to defend.
#[test]
fn the_natural_period_matches_beam_theory() {
    let (l, r, n) = (10.0, 0.15, 16);
    let m = Material::STEEL;
    let mut d = cantilever(m, l, r, n);
    undamped(&mut d);

    // Trial shape: the static deflection under a tip load, which for a
    // cantilever is within a fraction of a percent of the true first mode.
    let mut load = vec![Dof::default(); d.frame.nodes.len()];
    load[n].t = v3(0.0, 0.0, -1000.0);
    let s = d.frame.solve(&load);
    assert!(s.converged);
    let shape: Vec<Dof> = s
        .translation
        .iter()
        .zip(&s.rotation)
        .map(|(t, r)| Dof { t: *t, r: *r })
        .collect();

    let period = d.dominant_period(&shape);
    let expect = analytic_period(m, l, r);
    println!("  Rayleigh period {period:.6} s, analytic {expect:.6} s");
    assert!(
        (period - expect).abs() / expect < 0.02,
        "period {period:.6} s against analytic {expect:.6} s"
    );
}

/// Released from a deflected shape, the structure must actually oscillate, at
/// the period beam theory says.
///
/// Backward Euler is unconditionally stable and dissipative: it will not blow
/// up, but it lengthens the period by `O((omega h)^2)` and bleeds amplitude. A
/// timestep at a fortieth of the period keeps both under a percent, and the
/// test asserts the decay is monotone so a solver that *gains* energy fails
/// here rather than in the viewer.
#[test]
fn a_released_cantilever_oscillates_at_its_natural_period() {
    let (l, r, n) = (10.0, 0.15, 12);
    let m = Material::STEEL;
    let expect = analytic_period(m, l, r);
    let mut d = cantilever(m, l, r, n);
    undamped(&mut d);

    // Deflect statically under a tip load, then let go.
    let mut load = vec![Dof::default(); d.frame.nodes.len()];
    load[n].t = v3(0.0, 0.0, -2.0e4);
    let s = d.frame.solve(&load);
    assert!(s.converged);
    for i in 0..d.frame.nodes.len() {
        d.displacement[i] = Dof { t: s.translation[i], r: s.rotation[i] };
    }
    let start = d.displacement[n].t.z;
    assert!(start < 0.0, "the tip should start deflected downwards");

    let h = expect / 40.0;
    let free = vec![Dof::default(); d.frame.nodes.len()];
    let mut crossings = Vec::new();
    let mut previous = start;
    let mut peak: f64 = 0.0;
    for step in 1..200 {
        let rep = d.step(&free, h);
        assert!(rep.converged, "step {step} did not converge");
        assert!(rep.broken.is_empty(), "nothing should break under free decay");
        let tip = d.displacement[n].t.z;
        peak = peak.max(tip.abs());
        if (previous < 0.0) != (tip < 0.0) {
            // Linear interpolation to the crossing time.
            let frac = previous / (previous - tip);
            crossings.push((step as f64 - 1.0 + frac) * h);
        }
        previous = tip;
    }

    assert!(crossings.len() >= 4, "only {} zero crossings", crossings.len());
    // Successive crossings are half a period apart.
    let measured = 2.0 * (crossings[3] - crossings[1]) / 2.0;
    println!(
        "  measured period {measured:.6} s, analytic {expect:.6} s, peak {peak:.4} m \
         from a start of {:.4} m",
        start.abs()
    );
    assert!(
        (measured - expect).abs() / expect < 0.05,
        "measured {measured:.6} s against analytic {expect:.6} s"
    );
    // Free decay: no step may leave the structure with more energy than a
    // dissipative integrator started it with.
    assert!(
        peak <= start.abs() * 1.001,
        "amplitude grew from {:.4} to {peak:.4} m",
        start.abs()
    );
}

/// A load that arrives suddenly does about twice the work of the same load
/// standing still. This is the whole reason dynamics is worth solving.
///
/// For an undamped system a step load gives a dynamic load factor of exactly 2:
/// the structure overshoots the static position by as much as it deflected to
/// reach it. A quasi-static analysis of a gust therefore under-reads the stress
/// in it by a factor of two, which is the difference between a member at 60%
/// utilisation and one that has already failed.
#[test]
fn a_suddenly_applied_load_doubles_the_deflection() {
    let (l, r, n) = (10.0, 0.12, 12);
    let m = Material::STEEL;
    let mut d = cantilever(m, l, r, n);
    undamped(&mut d);

    let mut load = vec![Dof::default(); d.frame.nodes.len()];
    load[n].t = v3(0.0, 0.0, -3.0e3);
    let s = d.frame.solve(&load);
    assert!(s.converged);
    let statik = s.translation[n].z.abs();

    let h = analytic_period(m, l, r) / 200.0;
    let mut worst: f64 = 0.0;
    for _ in 0..400 {
        let rep = d.step(&load, h);
        assert!(rep.converged);
        worst = worst.max(d.displacement[n].t.z.abs());
    }
    let factor = worst / statik;
    println!("  static {statik:.5} m, dynamic peak {worst:.5} m, factor {factor:.3}");
    assert!(
        (factor - 2.0).abs() < 0.1,
        "dynamic load factor {factor:.3}, theory says 2"
    );
}

/// Damping must remove energy, and the report must say how much.
///
/// Both Rayleigh terms are checked separately: mass-proportional damping is
/// drag on the structure as a whole and stiffness-proportional damping is
/// internal friction, and a sign error in either is invisible when they are
/// only ever used together.
#[test]
fn damping_removes_energy_and_the_report_accounts_for_it() {
    let (l, r, n) = (8.0, 0.1, 10);
    let m = Material::STEEL;
    let period = analytic_period(m, l, r);

    for (label, mass_damping, stiff_damping) in
        [("mass", 1.5, 0.0), ("stiffness", 0.0, 2.0e-3), ("none", 0.0, 0.0)]
    {
        let mut d = cantilever(m, l, r, n);
        d.mass_damping = mass_damping;
        d.stiff_damping = stiff_damping;

        let mut load = vec![Dof::default(); d.frame.nodes.len()];
        load[n].t = v3(0.0, 0.0, -1.0e4);
        let s = d.frame.solve(&load);
        for i in 0..d.frame.nodes.len() {
            d.displacement[i] = Dof { t: s.translation[i], r: s.rotation[i] };
        }
        let start = d.strain_energy();

        let h = period / 60.0;
        let free = vec![Dof::default(); d.frame.nodes.len()];
        let mut dissipated = 0.0;
        for _ in 0..240 {
            let rep = d.step(&free, h);
            assert!(rep.converged);
            dissipated += rep.dissipated;
        }
        let left = d.kinetic_energy() + d.strain_energy();
        let fraction = left / start;
        println!(
            "  {label:>9} damping: {:.1}% of {start:.3e} J left after 4 cycles, \
             {dissipated:.3e} J accounted",
            fraction * 100.0
        );
        // The ledger must close: what is gone is what was reported gone.
        assert!(
            (start - left - dissipated).abs() <= 1e-6 * start.abs().max(left.abs()),
            "{label}: {start:.6e} - {left:.6e} != {dissipated:.6e}"
        );
        if mass_damping > 0.0 || stiff_damping > 0.0 {
            assert!(fraction < 0.5, "{label} damping left {:.1}%", fraction * 100.0);
            assert!(dissipated > 0.0, "{label} damping added energy");
        } else {
            // Undamped, only the integrator's own dissipation, which over four
            // cycles at this timestep is a few percent.
            assert!(fraction > 0.85, "undamped left only {:.1}%", fraction * 100.0);
        }
    }
}

/// A structure that loses a member must respond to having lost it.
///
/// The static solver could report that a member had failed. It could not show
/// what the rest of the structure then did, which is the part an observer
/// actually sees: the load the broken member was carrying arrives at its
/// neighbours all at once, and they move.
#[test]
fn breaking_a_member_moves_the_rest() {
    let (l, r, n) = (6.0, 0.02, 8);
    let m = Material::DRY_TIMBER;
    let mut d = cantilever(m, l, r, n);

    // A distributed load, ramped in over a couple of seconds so the failure is
    // the structure running out of strength rather than a shock wave from a
    // point force landing on one light node. Gravity is the shape of it; the
    // multiplier is what grows.
    let gravity = v3(0.0, 0.0, -9.80665);
    let lumped: Vec<f64> = d.frame.lumped.iter().map(|x| x.t.z).collect();
    let h = 5.0e-3;
    let ramp = 2.0;

    let mut broke_at = None;
    let mut speed_before = 0.0;
    let mut load = vec![Dof::default(); d.frame.nodes.len()];
    for step in 0..600 {
        let scale = 1.0 + 220.0 * (step as f64 * h / ramp).min(1.0);
        for (i, w) in lumped.iter().enumerate() {
            load[i].t = gravity.scale(w * scale);
        }
        let before = d.velocity[n].t.norm();
        let rep = d.step(&load, h);
        assert!(rep.converged, "step {step} did not converge");
        if !rep.broken.is_empty() && broke_at.is_none() {
            broke_at = Some((step, rep.broken.clone(), rep.released));
            speed_before = before;
        }
    }
    let (step, broken, released) = broke_at.expect("the cantilever should have failed");
    println!(
        "  failed at step {step} ({:.3} s), elements {broken:?}, released {released:.3} J, \
         tip speed {speed_before:.3} -> {:.3} m/s",
        step as f64 * h,
        d.velocity[n].t.norm()
    );
    // The root is where a cantilever's bending stress peaks, so that is what
    // must go first.
    assert!(broken.contains(&0), "the root should fail first, not {broken:?}");
    assert!(released > 0.0, "a member snapped holding no energy");
    assert!(
        d.velocity[n].t.norm() > speed_before * 2.0,
        "the tip should accelerate once the root lets go: {speed_before:.3} -> {:.3} m/s",
        d.velocity[n].t.norm()
    );
    assert!(
        d.displacement[n].t.z < -1.0,
        "the tip should fall, not hold at {:.4} m",
        d.displacement[n].t.z
    );
}

/// Rigid-body motion must cost nothing.
///
/// Translating an unsupported structure bodily produces no strain, so a
/// dynamic solver that charges it stiffness is producing forces out of nothing
/// — the classic way an assembled operator goes wrong without any test noticing.
#[test]
fn a_free_structure_drifts_without_straining() {
    let mut frame = Frame::new(Material::STEEL);
    let a = frame.add_node(v3(0.0, 0.0, 0.0), false);
    let b = frame.add_node(v3(2.0, 0.0, 0.0), false);
    let c = frame.add_node(v3(2.0, 2.0, 0.0), false);
    frame.add_beam(a, b, 0.05);
    frame.add_beam(b, c, 0.05);
    let mut d = Dynamics::new(frame);
    undamped(&mut d);

    let drift = v3(1.0, -0.5, 0.25);
    for v in d.velocity.iter_mut() {
        v.t = drift;
    }
    let free = vec![Dof::default(); d.frame.nodes.len()];
    let ke = d.kinetic_energy();
    for _ in 0..50 {
        let rep = d.step(&free, 0.01);
        assert!(rep.converged);
    }
    println!(
        "  strain after drifting {:.3e} J against {ke:.3e} J of kinetic energy",
        d.strain_energy()
    );
    assert!(
        d.strain_energy() <= 1e-9 * ke,
        "rigid translation strained the structure by {:.3e} J",
        d.strain_energy()
    );
    for i in 0..d.frame.nodes.len() {
        let v = d.velocity[i].t;
        assert!(
            (v - drift).norm() < 1e-9 * drift.norm(),
            "node {i} drifted at {v:?} instead of {drift:?}"
        );
    }
    let _ = Vec3::ZERO;
}

// ---------------------------------------------------------------------------
// Real structures, from the generators
// ---------------------------------------------------------------------------

use phys::morph::{Morphology, Program};
use phys::prolong::prolong_structured;
use phys::solvers::structure::*;
use phys::state::Aggregate;

/// The two solvers must agree.
///
/// Hold a load steady for long enough and the dynamics has to settle to exactly
/// the deflection the static analysis predicts. If it does not, one of them is
/// wrong, and the whole point of running the dynamic step on the static
/// operator was that they cannot be wrong separately.
#[test]
fn a_held_load_settles_to_the_static_answer() {
    let mut m = Morphology::planned(Program::Tower, 3.0e6, 11, 0x77);
    m.progress = 1.0;
    m.built = 3.0e6;
    let agg = Aggregate::neutral(3.0e6, m.extent(), 290.0, Program::Tower.substrate());
    let (bodies, topo, _) = prolong_structured(&agg, &m, 600, 7, 0x77, 0);

    let mut field = LoadField::new(bodies.len(), 290.0);
    field.apply(&weather::wind(22.0, v3(1.0, 0.0, 0.0)), &bodies, &topo);
    field.apply(&weather::gravity(), &bodies, &topo);

    // Static answer.
    let built = build_frame(&topo, bodies.len());
    let mut ds = dynamic_structure(&bodies, &topo).expect("a tower has members");
    let load = ds.nodal_loads(&field);
    let statik = built.frame.solve_with(&load, false);
    assert!(statik.converged, "the static solve did not converge");
    let apex = statik
        .translation
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.norm().total_cmp(&b.1.norm()))
        .map(|(i, _)| i)
        .unwrap();
    let target = statik.translation[apex];

    // Dynamic: heavily damped, so it settles rather than ringing forever.
    ds.dynamics.mass_damping = 6.0;
    ds.dynamics.stiff_damping = 1.0e-2;
    let mut worst_break = Vec::new();
    for _ in 0..400 {
        let rep = ds.advance(&field, 0.02);
        assert!(rep.converged, "the dynamic step did not converge");
        worst_break.extend(rep.broken);
    }
    assert!(worst_break.is_empty(), "the tower should survive a 22 m/s wind");

    let settled = ds.dynamics.displacement[apex].t;
    let err = (settled - target).norm() / target.norm().max(1e-12);
    println!(
        "  settled {:.5} m against static {:.5} m ({:.3}% apart)",
        settled.norm(),
        target.norm(),
        err * 100.0
    );
    assert!(
        err < 0.02,
        "settled at {settled:?} against a static answer of {target:?}"
    );
}

/// A tree in a gust must sway, and must still be swaying a moment later.
///
/// This is the behaviour the static path could never produce, and the reason
/// the open problem was worth closing: the structure has a period of its own,
/// and what an observer sees is that period, not the load.
#[test]
fn a_tree_sways_and_rings_down() {
    let mut m = Morphology::new(Program::Tree, 0xACE, 0x1234, 0);
    m.built = 900.0;
    m.age = 40.0 * phys::units::YEAR;
    let mut agg = Aggregate::neutral(900.0, m.extent(), 291.0, Program::Tree.substrate());
    agg.chemical_energy = m.stored_energy();
    let (bodies, topo, _) = prolong_structured(&agg, &m, 400, 7, 0x1234, 0);

    let mut ds = dynamic_structure(&bodies, &topo).expect("a tree has members");
    let tip = (0..ds.dynamics.frame.nodes.len())
        .max_by(|&a, &b| {
            ds.dynamics.frame.nodes[a].z.total_cmp(&ds.dynamics.frame.nodes[b].z)
        })
        .unwrap();

    let gust = {
        let mut f = LoadField::new(bodies.len(), 291.0);
        f.apply(&weather::wind(18.0, v3(1.0, 0.0, 0.0)), &bodies, &topo);
        f
    };
    let calm = LoadField::new(bodies.len(), 291.0);

    let h = 0.02;
    let mut peak: f64 = 0.0;
    for _ in 0..50 {
        let rep = ds.advance(&gust, h);
        assert!(rep.converged);
        peak = peak.max(ds.dynamics.displacement[tip].t.x);
    }
    assert!(peak > 1e-4, "the tree did not move in an 18 m/s gust: {peak:.3e} m");

    // Gust drops. The crown must swing back through its rest position rather
    // than stopping where it was pushed to.
    let mut crossed = 0;
    let mut previous = ds.dynamics.displacement[tip].t.x;
    let mut swing_back: f64 = 0.0;
    for _ in 0..200 {
        let rep = ds.advance(&calm, h);
        assert!(rep.converged);
        let x = ds.dynamics.displacement[tip].t.x;
        if (previous > 0.0) != (x > 0.0) {
            crossed += 1;
        }
        swing_back = swing_back.min(x);
        previous = x;
    }
    println!(
        "  crown pushed {peak:.4} m, swung back to {swing_back:.4} m, crossed rest \
         {crossed} times, residual {:.4} m",
        ds.dynamics.displacement[tip].t.x.abs()
    );
    assert!(crossed >= 1, "the crown never came back through its rest position");
    assert!(swing_back < 0.0, "the crown never overshot the other way");
    // And it must eventually stop: material damping is not optional.
    assert!(
        ds.dynamics.displacement[tip].t.x.abs() < peak,
        "the tree is still swinging as far as it was pushed"
    );
}
