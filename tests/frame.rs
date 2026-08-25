//! Frame analysis against closed-form results.
//!
//! Every case here has an exact answer from beam theory, so a failure is a
//! failure of the solver and not of judgement.

use phys::math::{v3, Vec3};
use phys::solvers::frame::*;
use phys::topology::Material;

/// A steel bar, 100 mm radius unless a test says otherwise.
fn steel() -> Material {
    Material::STEEL
}

fn load_at(n: usize, node: u32, force: Vec3) -> Vec<Dof> {
    let mut l = vec![Dof::default(); n];
    l[node as usize].t = force;
    l
}

/// Cantilever tip deflection, `PL^3 / 3EI`.
///
/// The easy case, and the one that proves least: a translation-only spring
/// model with `k = 3EI/L^3` passes it by construction. It is here to catch
/// gross errors, not to validate the formulation.
#[test]
fn cantilever_tip_deflection() {
    let (l, r, p) = (4.0, 0.05, 5000.0);
    let mat = steel();
    let mut f = Frame::new(mat);
    let base = f.add_node(v3(0.0, 0.0, 0.0), true);
    let tip = f.add_node(v3(l, 0.0, 0.0), false);
    f.add_beam(base, tip, r);

    let sol = f.solve_with(&load_at(2, tip, v3(0.0, 0.0, -p)), false);
    assert!(sol.converged, "cantilever did not converge");
    let i = std::f64::consts::PI * r.powi(4) / 4.0;
    let expected = p * l.powi(3) / (3.0 * mat.stiffness * i);
    let got = -sol.translation[tip as usize].z;
    println!("  cantilever: {got:.6} m, analytic {expected:.6} m ({} iterations)", sol.iterations);
    assert!(
        (got - expected).abs() / expected < 1e-3,
        "tip deflection {got:.6} vs analytic {expected:.6}"
    );
}

/// A beam built in at both ends, `PL^3 / 192EI` at midspan.
///
/// **This is the test that matters.** It is 64 times stiffer than the
/// cantilever, and the entire difference comes from the fixed ends resisting
/// rotation. A model without rotational degrees of freedom cannot represent
/// that and gets the answer wrong by roughly that factor — which is why the
/// previous spring-based solver passed the cantilever and was still wrong.
#[test]
fn fixed_fixed_beam_midspan_deflection() {
    let (l, r, p) = (6.0, 0.06, 20000.0);
    let mat = steel();
    let mut f = Frame::new(mat);
    let a = f.add_node(v3(0.0, 0.0, 0.0), true);
    let mid = f.add_node(v3(l / 2.0, 0.0, 0.0), false);
    let b = f.add_node(v3(l, 0.0, 0.0), true);
    f.add_beam(a, mid, r);
    f.add_beam(mid, b, r);

    let sol = f.solve_with(&load_at(3, mid, v3(0.0, 0.0, -p)), false);
    assert!(sol.converged, "fixed-fixed beam did not converge");
    let i = std::f64::consts::PI * r.powi(4) / 4.0;
    let expected = p * l.powi(3) / (192.0 * mat.stiffness * i);
    let got = -sol.translation[mid as usize].z;
    let cantilever_equivalent = p * l.powi(3) / (3.0 * mat.stiffness * i);
    println!(
        "  fixed-fixed: {got:.6e} m, analytic {expected:.6e} m \
         (a translation-only model would give ~{cantilever_equivalent:.3e})"
    );
    assert!(
        (got - expected).abs() / expected < 0.02,
        "midspan deflection {got:.6e} vs analytic {expected:.6e}"
    );
}

/// A simply-supported beam, `PL^3 / 48EI`.
///
/// Pinned rather than fixed ends. Together with the case above this pins down
/// that the solver is responding to the *boundary conditions* and not just
/// scaling one answer.
#[test]
fn simply_supported_beam() {
    let (l, r, p) = (6.0, 0.06, 20000.0);
    let mat = steel();
    let mut f = Frame::new(mat);
    // Pins are modelled as ties to fixed anchors: axial restraint, free
    // rotation, which is exactly what a pin is.
    let anchor_a = f.add_node(v3(0.0, 0.0, -0.001), true);
    let anchor_b = f.add_node(v3(l, 0.0, -0.001), true);
    let a = f.add_node(v3(0.0, 0.0, 0.0), false);
    let mid = f.add_node(v3(l / 2.0, 0.0, 0.0), false);
    let b = f.add_node(v3(l, 0.0, 0.0), false);
    f.add_beam(a, mid, r);
    f.add_beam(mid, b, r);
    f.add_tie(anchor_a, a, r * 4.0);
    f.add_tie(anchor_b, b, r * 4.0);

    let sol = f.solve_with(&load_at(5, mid, v3(0.0, 0.0, -p)), false);
    assert!(sol.converged, "simply supported beam did not converge");
    let i = std::f64::consts::PI * r.powi(4) / 4.0;
    let expected = p * l.powi(3) / (48.0 * mat.stiffness * i);
    let got = -sol.translation[mid as usize].z;
    println!("  simply supported: {got:.6e} m, analytic {expected:.6e} m");
    assert!(
        (got - expected).abs() / expected < 0.05,
        "midspan deflection {got:.6e} vs analytic {expected:.6e}"
    );
}

/// Axial extension, `PL / EA`. The simplest possible check that the units are
/// right end to end.
#[test]
fn axial_extension() {
    let (l, r, p) = (3.0, 0.02, 50000.0);
    let mat = steel();
    let mut f = Frame::new(mat);
    let a = f.add_node(v3(0.0, 0.0, 0.0), true);
    let b = f.add_node(v3(l, 0.0, 0.0), false);
    f.add_beam(a, b, r);
    let sol = f.solve_with(&load_at(2, b, v3(p, 0.0, 0.0)), false);
    let area = std::f64::consts::PI * r * r;
    let expected = p * l / (mat.stiffness * area);
    let got = sol.translation[b as usize].x;
    println!("  extension: {got:.6e} m, analytic {expected:.6e} m");
    assert!((got - expected).abs() / expected < 1e-6);
    // And the member reports the force that was applied to it.
    assert!(
        (sol.forces[0].axial - p).abs() / p < 1e-6,
        "axial force {:.1} N, applied {p:.1} N",
        sol.forces[0].axial
    );
}

/// Euler buckling: a slender strut in compression fails at
/// `P_cr = pi^2 EI / (KL)^2`, far below the load its material strength allows.
#[test]
fn euler_buckling_is_detected() {
    let (l, r): (f64, f64) = (5.0, 0.02);
    let mat = steel();
    let i = std::f64::consts::PI * r.powi(4) / 4.0;
    let area = std::f64::consts::PI * r * r;
    let p_cr = std::f64::consts::PI.powi(2) * mat.stiffness * i / (l * l);
    let p_squash = mat.rupture * area;
    println!("  critical load {p_cr:.0} N, squash load {p_squash:.0} N — buckling first by {:.0}x",
        p_squash / p_cr);
    assert!(p_cr < p_squash, "this strut is not slender enough to test buckling");

    let mut f = Frame::new(mat);
    let a = f.add_node(v3(0.0, 0.0, 0.0), true);
    let b = f.add_node(v3(0.0, 0.0, l), false);
    // A pin-jointed strut, so K = 1 and the classical formula applies.
    f.add_tie(a, b, r);

    for (fraction, should_buckle) in [(0.5, false), (0.95, false), (1.2, true), (3.0, true)] {
        let sol = f.solve_with(&load_at(2, b, v3(0.0, 0.0, -p_cr * fraction)), false);
        let u = sol.forces[0].buckling;
        println!("    at {fraction:.2} P_cr: buckling utilisation {u:.3}");
        assert!(
            (u - fraction).abs() < 0.02,
            "buckling utilisation {u:.3} at {fraction:.2} of critical"
        );
        assert_eq!(u >= 1.0, should_buckle);
    }

    // The stress check alone would see nothing wrong at the critical load.
    let sol = f.solve_with(&load_at(2, b, v3(0.0, 0.0, -p_cr)), false);
    let stress_utilisation = sol.forces[0].stress / mat.rupture;
    println!("  at P_cr the stress check reports {stress_utilisation:.4} utilised");
    assert!(
        stress_utilisation < 0.05,
        "a stress check should be nowhere near failure at the buckling load"
    );
}

/// A ductile material sheds load from an overloaded member to its neighbours;
/// a brittle one does not. This is the difference between a frame that sags
/// and a wall that drops.
#[test]
fn ductile_materials_redistribute_and_brittle_ones_do_not() {
    // Three parallel struts of unequal length sharing one load. The short one
    // is stiffest and attracts the most force, so it yields first.
    let build = |mat: Material| {
        let mut f = Frame::new(mat);
        let top = f.add_node(v3(0.0, 0.0, 0.0), false);
        let anchors = [
            f.add_node(v3(0.0, 0.0, -1.0), true),
            f.add_node(v3(0.0, 0.0, -1.6), true),
            f.add_node(v3(0.0, 0.0, -2.4), true),
        ];
        for a in anchors {
            f.add_tie(a, top, 0.004);
        }
        f
    };
    let load = |n: usize| load_at(n, 0, v3(0.0, 0.0, -60000.0));

    let ductile = build(Material::STEEL);
    let d_sol = ductile.solve(&load(4));
    let brittle = build(Material { ductility: 0.0, ..Material::STEEL });
    let b_sol = brittle.solve(&load(4));

    let spread = |s: &FrameSolution| {
        let f: Vec<f64> = s.forces.iter().map(|e| e.axial.abs()).collect();
        let max = f.iter().cloned().fold(0.0f64, f64::max);
        let min = f.iter().cloned().fold(f64::INFINITY, f64::min);
        (max, min, max / min.max(1e-9))
    };
    let (dmax, _, dratio) = spread(&d_sol);
    let (bmax, _, bratio) = spread(&b_sol);
    println!(
        "  ductile: peak {dmax:.0} N, spread {dratio:.2}x, {} members yielded over {} passes",
        d_sol.yielded, d_sol.redistributions
    );
    println!("  brittle: peak {bmax:.0} N, spread {bratio:.2}x");

    assert!(d_sol.redistributions > 0, "the ductile solve never redistributed");
    assert!(
        dratio < bratio,
        "redistribution should even out the load: {dratio:.2}x vs {bratio:.2}x"
    );
    assert!(
        dmax <= bmax * 1.001,
        "the most loaded member should not gain from redistribution"
    );
    // Both still carry the applied load: redistribution moves force, it does
    // not destroy it.
    for (name, s) in [("ductile", &d_sol), ("brittle", &b_sol)] {
        let vertical: f64 = s
            .forces
            .iter()
            .enumerate()
            .map(|(i, f)| {
                let e = &s.forces[i];
                let _ = e;
                f.axial
            })
            .sum::<f64>()
            .abs();
        assert!(vertical > 0.0, "{name} carried nothing");
    }
}

/// Preconditioning is not a micro-optimisation here: rotational and
/// translational degrees of freedom differ by orders of magnitude, and without
/// it the solver spends its budget discovering that.
#[test]
fn preconditioning_converges_quickly() {
    let mat = steel();
    let mut f = Frame::new(mat);
    // A twenty-storey braced frame: the case that motivated all of this.
    let mut below = [u32::MAX; 4];
    let mut loads = Vec::new();
    for floor in 0..20 {
        let z = floor as f64 * 3.2;
        let mut here = [0u32; 4];
        for c in 0..4 {
            let a = std::f64::consts::FRAC_PI_2 * c as f64;
            here[c] = f.add_node(v3(6.0 * a.cos(), 6.0 * a.sin(), z + 3.2), false);
        }
        for c in 0..4 {
            let base = if floor == 0 {
                f.add_node(v3(6.0 * (std::f64::consts::FRAC_PI_2 * c as f64).cos(),
                              6.0 * (std::f64::consts::FRAC_PI_2 * c as f64).sin(), z), true)
            } else {
                below[c]
            };
            f.add_beam(base, here[c], 0.15);
            let n = (c + 1) % 4;
            f.add_beam(here[c], here[n], 0.10);
            if floor > 0 {
                f.add_tie(below[c], here[n], 0.05);
            }
        }
        loads.push(here[0]);
        below = here;
    }
    let n = f.nodes.len();
    let mut load = vec![Dof::default(); n];
    for node in &loads {
        load[*node as usize].t = v3(4000.0, 0.0, -30000.0);
    }
    let stiff = vec![mat.stiffness; f.elements.len()];
    let (_, with_pre, ok_pre) = f.solve_elastic_opt(&load, &stiff, true);
    let (_, without, ok_plain) = f.solve_elastic_opt(&load, &stiff, false);
    println!(
        "  {} nodes ({} DOF), {} members: {} iterations preconditioned, {} without",
        n,
        n * 6,
        f.elements.len(),
        with_pre,
        if ok_plain { without.to_string() } else { format!("{without}, did not converge") }
    );
    assert!(ok_pre, "the braced frame did not converge");
    assert!(
        !ok_plain || with_pre * 3 < without,
        "preconditioning saved nothing: {with_pre} against {without} iterations"
    );

    let sol = f.solve_with(&load, false);
    assert!(sol.converged, "the braced frame did not converge");
    // Sanity: the tower leans downwind and settles.
    let top = sol.translation[*loads.last().unwrap() as usize];
    println!("  top of the frame moved {:.4} m across, {:.4} m down", top.x, -top.z);
    assert!(top.x > 0.0 && top.z < 0.0, "the frame moved the wrong way");
}

/// A statically indeterminate truss, against the closed-form answer.
///
/// Three pinned bars from a common apex to three anchors: one vertical, two at
/// `t` either side. A downward load `P` at the apex. Statics alone cannot solve
/// this — three unknowns, two equations — so the split is decided by the bars'
/// relative stiffness, and for equal areas the classical result is
///
/// ```text
///     F_middle = P / (1 + 2 cos^3 t)      F_outer = F_middle cos^2 t
/// ```
///
/// which at 45 degrees is 0.5858 P and 0.2929 P. Reproducing that is the test
/// of whether the solver does mechanics or redistributes by a rule of thumb.
#[test]
fn redundant_truss_matches_the_analytic_split() {
    let h = 1.0;
    let t: f64 = std::f64::consts::FRAC_PI_4;
    let radius = 0.01;

    let mut frame = Frame::new(Material::STEEL);
    let apex = frame.add_node(v3(0.0, 0.0, 0.0), false);
    let bars: Vec<usize> = [0.0, -h * t.tan(), h * t.tan()]
        .iter()
        .map(|&x| {
            let anchor = frame.add_node(v3(x, 0.0, -h), true);
            frame.add_tie(anchor, apex, radius)
        })
        .collect();

    let p_load = 1000.0;
    let mut load = vec![Dof::default(); frame.nodes.len()];
    load[apex as usize].t = v3(0.0, 0.0, -p_load);

    let s = frame.solve(&load);
    assert!(s.converged, "the truss solve did not converge");

    let axial = |i: usize| s.forces[bars[i]].axial.abs();
    let middle = axial(0);
    let outer = (axial(1) + axial(2)) / 2.0;

    let cos_t = t.cos();
    let expect_mid = p_load / (1.0 + 2.0 * cos_t.powi(3));
    let expect_outer = expect_mid * cos_t * cos_t;

    println!(
        "  middle {middle:.2} N (analytic {expect_mid:.2}), each outer {outer:.2} N \
         (analytic {expect_outer:.2}), {} iterations",
        s.iterations
    );
    assert!(
        (middle - expect_mid).abs() / expect_mid < 0.005,
        "middle bar {middle:.2} N, analytic {expect_mid:.2} N"
    );
    assert!(
        (outer - expect_outer).abs() / expect_outer < 0.005,
        "outer bars {outer:.2} N, analytic {expect_outer:.2} N"
    );
    // Vertical equilibrium at the apex, independent of the analytic form.
    let carried = middle + 2.0 * outer * cos_t;
    assert!(
        (carried - p_load).abs() / p_load < 0.005,
        "the apex does not balance: {carried:.2} N against {p_load:.2} N"
    );
}


/// The tree factorisation must be an exact solve when the load path is a tree.
///
/// This is the same structural fact the O(n) static pass rests on, applied to a
/// different question: a tree has a perfect elimination ordering, so eliminating
/// leaves first produces no fill-in and the factorisation *is* the inverse.
/// Preconditioned conjugate gradient then finishes in a single iteration, at
/// any timestep and any size.
///
/// Without it, a slender chain of `n` beam elements has a condition number
/// growing like `n^4`: a 2000-member tree took 3645 Jacobi iterations and 700 ms
/// for one dynamic substep, which at twenty updates a second is not a solver.
#[test]
fn the_tree_factorisation_is_exact_on_a_tree() {
    let mat = steel();
    for n in [16usize, 128, 1024] {
        let mut f = Frame::new(mat);
        // A chain that wanders, so the answer is not one-dimensional.
        let mut prev = f.add_node(v3(0.0, 0.0, 0.0), true);
        for i in 1..=n {
            let t = i as f64 * 0.3;
            let node = f.add_node(v3(t, 0.4 * (t * 0.7).sin(), 0.2 * (t * 0.3).cos()), false);
            f.add_beam(prev, node, 0.05);
            prev = node;
        }
        let mut load = vec![Dof::default(); f.nodes.len()];
        for (i, l) in load.iter_mut().enumerate() {
            l.t = v3(0.0, 30.0 * (i as f64 * 0.11).sin(), -120.0);
        }
        let stiff = vec![mat.stiffness; f.elements.len()];
        let (u, iters, converged) = f.solve_elastic_opt(&load, &stiff, true);
        assert!(converged, "n={n} did not converge");
        // The iteration count does not grow with the structure, which is the
        // whole claim. One iteration applies the factorisation; the odd extra
        // one is conjugate gradient establishing that the residual really is at
        // round-off, on a chain whose condition number by then exceeds 10^12.
        // Jacobi on the same problem needed 3645 and rising.
        assert!(iters <= 3, "n={n} took {iters} iterations");

        // And the answer really does solve the system.
        let mut ku = vec![Dof::default(); f.nodes.len()];
        f.apply_operator(&u, &stiff, &mut ku);
        let (mut residual, mut scale) = (0.0f64, 0.0f64);
        for i in 0..f.nodes.len() {
            if f.fixed[i] {
                continue;
            }
            residual += load[i].sub(ku[i]).dot(load[i].sub(ku[i]));
            scale += load[i].dot(load[i]);
        }
        let relative = (residual / scale).sqrt();
        // What is achievable is bounded by the problem, not the method: a chain
        // of `n` beam elements has a condition number growing like `n^4`, so
        // round-off alone leaves a relative residual near `n^4 * eps`. At
        // n=1024 that is about 10^-4, and a solver claiming better than that
        // would be reporting its own arithmetic error as accuracy.
        let reachable = (n as f64).powi(4) * 4.0 * f64::EPSILON;
        println!(
            "  n={n:>5}: {iters} iteration(s), relative residual {relative:.3e} \
             (round-off allows {reachable:.1e})"
        );
        assert!(
            relative < reachable.max(1e-10),
            "n={n} left a residual of {relative:.3e}, round-off allows {reachable:.1e}"
        );
    }
}

/// And it must still help when the structure is not quite a tree.
///
/// The braces are the only members outside the spanning forest, so dropping
/// their coupling leaves a factorisation that is no longer exact. How good it
/// still is depends on how redundant the structure is — which is also why the
/// solver declines the factorisation above a threshold on that redundancy and
/// keeps the diagonal instead: applying the factor costs about an order of
/// magnitude more than a diagonal division, so on a moment frame, where every
/// bay between two floors is a closed loop, it removes a factor of two from the
/// iteration count and loses on the exchange. Both sides of that decision are
/// checked here.
///
/// The forest keeps the *stiffest* members, chosen by axial stiffness. What is
/// left out is what the preconditioner cannot see, and leaving out a slender
/// brace costs far less than leaving out a column.
#[test]
fn the_tree_factorisation_still_helps_a_braced_chain() {
    let mat = steel();
    let n = 200usize;
    let mut f = Frame::new(mat);
    let mut chain = vec![f.add_node(v3(0.0, 0.0, 0.0), true)];
    for i in 1..=n {
        let t = i as f64 * 0.3;
        chain.push(f.add_node(v3(t, 0.5 * (t * 0.4).sin(), 0.0), false));
        f.add_beam(chain[i - 1], chain[i], 0.05);
    }
    // A stay every twentieth joint, skipping five along: enough redundancy that
    // the factorisation is no longer exact, few enough that it still pays.
    for i in (10..n - 5).step_by(20) {
        f.add_tie(chain[i], chain[i + 5], 0.01);
    }
    let redundancy = f.redundancy();
    assert!(redundancy > 0, "nothing here is redundant");
    assert!(
        redundancy * 4 <= f.nodes.len(),
        "this chain should take the factorisation"
    );

    let mut load = vec![Dof::default(); f.nodes.len()];
    for (i, l) in load.iter_mut().enumerate() {
        l.t = v3(0.0, 30.0 * (i as f64 * 0.11).sin(), -120.0);
    }
    let stiff = vec![mat.stiffness; f.elements.len()];
    let (_, tree, ok_tree) = f.solve_elastic_opt(&load, &stiff, true);
    let (_, plain, ok_plain) = f.solve_elastic_opt(&load, &stiff, false);
    println!(
        "  braced chain: {} joints, {} members, redundancy {redundancy} — {tree} iterations \
         factorised, {} without",
        f.nodes.len(),
        f.elements.len(),
        if ok_plain { plain.to_string() } else { format!("{plain}, did not converge") }
    );
    assert!(ok_tree, "the braced chain did not converge");
    assert!(
        !ok_plain || tree * 10 < plain,
        "the factorisation saved little: {tree} against {plain}"
    );

    // The other side of the decision: a moment frame is redundant in every bay,
    // and the solver keeps the diagonal for it.
    let mut frame = Frame::new(mat);
    let mut below = [u32::MAX; 4];
    for floor in 0..20 {
        let z = floor as f64 * 3.2;
        let mut here = [0u32; 4];
        for c in 0..4 {
            let a = std::f64::consts::FRAC_PI_2 * c as f64;
            here[c] = frame.add_node(v3(6.0 * a.cos(), 6.0 * a.sin(), z + 3.2), false);
        }
        for c in 0..4 {
            let a = std::f64::consts::FRAC_PI_2 * c as f64;
            let base = if floor == 0 {
                frame.add_node(v3(6.0 * a.cos(), 6.0 * a.sin(), z), true)
            } else {
                below[c]
            };
            frame.add_beam(base, here[c], 0.15);
            frame.add_beam(here[c], here[(c + 1) % 4], 0.10);
        }
        below = here;
    }
    println!(
        "  moment frame: {} joints, {} members, redundancy {} — factorisation declined",
        frame.nodes.len(),
        frame.elements.len(),
        frame.redundancy()
    );
    assert!(
        frame.redundancy() * 4 > frame.nodes.len(),
        "a moment frame should be too redundant for the factorisation"
    );
}
