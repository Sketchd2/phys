//! Where the atoms actually sit.
//!
//! # Why an analyser needs geometry
//!
//! Two bond dipoles of equal size either cancel or add, and which one decides
//! whether a substance dissolves in water. Carbon dioxide and water have
//! comparably polar bonds; carbon dioxide is linear so they cancel exactly and
//! it is a gas that barely dissolves, and water is bent at 104.5 degrees so
//! they add and it is the solvent everything else is measured against. No
//! property of the *graph* separates those two cases. Only the shape does.
//!
//! The first attempt tried to avoid building geometry by asking whether the
//! atoms were interchangeable, which said water's two hydrogens were a
//! symmetric pair and gave it a dipole of exactly zero. Interchangeable is not
//! the same as opposed, and there is no way to patch around that: the shape has
//! to be built.
//!
//! # How the shape is derived
//!
//! VSEPR, which is a first-principles argument rather than a table: electron
//! pairs around an atom repel, so they sit as far apart as they can get. The
//! number of pairs is bonds plus lone pairs, and lone pairs come from the
//! element's own valence electron count — measured, per element, like its mass.
//!
//! Two domains put bonds at 180 degrees, three at 120, four at 109.47, five and
//! six at their own arrangements. Lone pairs occupy domains without being
//! bonds, which is why oxygen with two bonds and two lone pairs is bent rather
//! than linear, and they push harder than bonding pairs, which is why the angle
//! closes from 109.47 to about 104.5.
//!
//! Bond lengths come from `analyse::bond_length`. So the conformer is built
//! entirely out of per-element constants and the graph.
//!
//! # What it is not
//!
//! A relaxed structure. This is a plausible conformer built by walking the
//! bond graph outward, with any ring then pulled closed onto the same derived
//! distances (`close_rings`), and for anything with a ring or a long flexible
//! chain it will differ from the real minimum. It is good enough for the two things
//! it is for — a dipole, and a starting configuration for the molecular
//! dynamics solver to relax properly — and `md.rs` is where a structure that
//! has to be right goes next.

use super::arrange::{Arrangement, Order};
use crate::math::Vec3;

/// How far apart electron domains sit, radians, for `n` of them.
///
/// The VSEPR arrangements. Beyond six the geometry stops being a single angle
/// and the tetrahedral value is used as a stand-in, which no arrangement this
/// engine builds reaches.
pub fn domain_angle(domains: usize, lone_pairs: usize) -> f64 {
    let base: f64 = match domains {
        0 | 1 => 180.0,
        2 => 180.0,
        3 => 120.0,
        4 => 109.4712206,
        5 => 102.0,
        _ => 90.0,
    };
    // A lone pair is held closer to the nucleus than a bonding pair and so
    // pushes the bonds together: about two and a half degrees each. That is
    // the difference between methane's 109.5, ammonia's 107 and water's 104.5,
    // and it is the whole reason water is as polar as it is.
    (base - 2.5 * lone_pairs as f64).to_radians()
}

/// Electron domains and lone pairs at one atom.
pub fn domains_at(arr: &Arrangement, atom: usize) -> (usize, usize) {
    let adj = arr.neighbours();
    let element = arr.atoms[atom];
    let bonds = adj[atom].iter().filter(|(_, o)| *o != Order::Hydrogen).count();
    let used: usize = adj[atom]
        .iter()
        .map(|(_, o)| match o {
            Order::Single | Order::Ionic => 1,
            Order::Double => 2,
            Order::Triple => 3,
            Order::Hydrogen => 0,
        })
        .sum();
    let electrons = element.valence_electrons().unwrap_or(4);
    // Every valence electron not in a bond is half of a lone pair.
    let lone = electrons.saturating_sub(used) / 2;
    (bonds + lone, lone)
}

/// A conformer: one position per atom, metres, centred on the first atom.
///
/// Built by walking outward from the most connected atom, placing each new
/// neighbour at its bond length and at the VSEPR angle to the bonds already
/// placed there. Deterministic, so the same arrangement always gives the same
/// shape.
pub fn embed(arr: &Arrangement) -> Vec<Vec3> {
    let n = arr.atoms.len();
    let mut pos = vec![Vec3::ZERO; n];
    if n == 0 {
        return pos;
    }
    let adj = arr.neighbours();

    // Start from the most connected atom, so the first shell is the one the
    // geometry is most constrained by.
    let root = (0..n).max_by_key(|i| adj[*i].len()).unwrap_or(0);
    let mut placed = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    placed[root] = true;
    let mut queue = std::collections::VecDeque::from([root]);

    while let Some(i) = queue.pop_front() {
        let (domains, lone) = domains_at(arr, i);
        let angle = domain_angle(domains, lone);

        // Directions already spoken for at this atom: back toward whatever
        // placed it, plus anything else already positioned.
        let mut taken: Vec<Vec3> = adj[i]
            .iter()
            .filter(|(j, _)| placed[*j])
            .map(|(j, _)| (pos[*j] - pos[i]).unit())
            .filter(|d| d.is_finite() && d.norm() > 0.5)
            .collect();

        let pending: Vec<(usize, Order)> =
            adj[i].iter().copied().filter(|(j, _)| !placed[*j]).collect();
        let mut left = pending.len();

        // Which way round the bond to the parent the new bonds turn. Left to a
        // fixed helper axis, every atom of a chain turned the same way and the
        // chain curled back into itself — pentane put two carbons 0.62 A
        // apart, octane two atoms 0.17 A apart — so the new bonds are set
        // anti to the bond the parent came from: the zig-zag of a stretched
        // chain, which is also as far apart as their electron pairs can get.
        let reference = (parent[i] != usize::MAX)
            .then(|| {
                let p = parent[i];
                adj[p].iter().find(|(g, _)| *g != i && placed[*g]).map(|(g, _)| pos[p] - pos[*g])
            })
            .flatten();

        for (j, order) in pending {
            // Each new direction avoids every direction already chosen —
            // including the ones chosen a moment ago in this same loop. Not
            // accumulating them was what put both of carbon dioxide's oxygens
            // at the same point, giving a linear molecule a dipole.
            let dir = direction_avoiding(&taken, angle, left, reference);
            let length =
                super::analyse::bond_length(arr.atoms[i], arr.atoms[j], order).unwrap_or(1.5e-10);
            pos[j] = pos[i] + dir.scale(length);
            placed[j] = true;
            parent[j] = i;
            taken.push(dir);
            left -= 1;
            queue.push_back(j);
        }
    }

    // Anything the walk could not reach — a lattice's formula unit is not
    // required to be connected — is laid out on a line so it at least has
    // distinct positions.
    let mut spare = 1.0;
    for i in 0..n {
        if !placed[i] {
            pos[i] = Vec3 { x: spare * 3e-10, y: 0.0, z: 0.0 };
            spare += 1.0;
        }
    }
    if has_ring(arr) {
        close_rings(arr, &mut pos);
    }
    pos
}

/// Whether the bond graph has a cycle: more edges than a forest of the same
/// atoms could hold.
fn has_ring(arr: &Arrangement) -> bool {
    let n = arr.atoms.len();
    let mut root: Vec<usize> = (0..n).collect();
    fn find(root: &mut [usize], mut i: usize) -> usize {
        while root[i] != i {
            root[i] = root[root[i]];
            i = root[i];
        }
        i
    }
    for (i, adj) in arr.neighbours().iter().enumerate() {
        for &(j, _) in adj {
            if j < i {
                continue;
            }
            let (a, b) = (find(&mut root, i), find(&mut root, j));
            if a == b {
                return true;
            }
            root[a] = b;
        }
    }
    false
}

/// Pull the walk's conformer onto every bond it was meant to have.
///
/// The walk places atoms along a spanning tree of the bond graph, so a bond
/// that closes a ring is never looked at: benzene's last bond came out 5.1 A
/// long, which gave it a dipole of 0.52 D, and an electronic-structure
/// relaxation started from that shape has no bond there to follow. So when
/// there is a ring, the conformer is relaxed by least squares onto the
/// distances the walk already derives — each bond at `bond_length`, and each
/// pair of atoms sharing a neighbour at the distance the VSEPR angle there puts
/// them — every one of them, ring-closing bonds included. Nothing new is
/// stated: the targets are the walk's own. For a six-membered ring at 120
/// degrees those distances admit only the flat hexagon, which is benzene.
///
/// Least squares on the relative error, by BFGS.
/// Deterministic, so the same arrangement always gives the same shape.
fn close_rings(arr: &Arrangement, pos: &mut [Vec3]) {
    let adj = arr.neighbours();
    let mut targets: Vec<(usize, usize, f64)> = Vec::new();
    for b in &arr.bonds {
        let (i, j) = (b.a as usize, b.b as usize);
        if i == j || i >= pos.len() || j >= pos.len() {
            continue;
        }
        let d = super::analyse::bond_length(arr.atoms[i], arr.atoms[j], b.order).unwrap_or(1.5e-10);
        targets.push((i, j, d));
    }
    for (c, around) in adj.iter().enumerate() {
        let (domains, lone) = domains_at(arr, c);
        let angle = domain_angle(domains, lone);
        for (x, &(j, oj)) in around.iter().enumerate() {
            for &(k, ok) in &around[x + 1..] {
                if j == k {
                    continue;
                }
                let a = super::analyse::bond_length(arr.atoms[c], arr.atoms[j], oj).unwrap_or(1.5e-10);
                let b = super::analyse::bond_length(arr.atoms[c], arr.atoms[k], ok).unwrap_or(1.5e-10);
                let d = (a * a + b * b - 2.0 * a * b * angle.cos()).max(0.0).sqrt();
                targets.push((j, k, d));
            }
        }
    }

    // Work in units of the first target so the arithmetic is of order one.
    let scale = targets.first().map(|t| t.2).unwrap_or(1e-10);
    let x: Vec<Vec3> = pos.iter().map(|p| p.scale(1.0 / scale)).collect();
    let t: Vec<(usize, usize, f64)> = targets.iter().map(|&(i, j, d)| (i, j, d / scale)).collect();

    let cost = |x: &[Vec3], grad: Option<&mut Vec<Vec3>>| -> f64 {
        let mut e = 0.0;
        let mut g = grad;
        if let Some(g) = g.as_deref_mut() {
            g.iter_mut().for_each(|v| *v = Vec3::ZERO);
        }
        for &(i, j, d0) in &t {
            let r = x[j] - x[i];
            let d = r.norm().max(1e-9);
            let u = (d - d0) / d0;
            e += u * u;
            if let Some(g) = g.as_deref_mut() {
                let f = r.scale(2.0 * u / (d0 * d));
                g[j] += f;
                g[i] -= f;
            }
        }
        e
    };

    // BFGS on the flattened coordinates, with a backtracking line search.
    // Steepest descent was tried first and crawls along the soft twist of a
    // ring: benzene was still 8e-4 off after 20,000 steps.
    let m = 3 * x.len();
    let flat = |x: &[Vec3]| -> Vec<f64> { x.iter().flat_map(|v| [v.x, v.y, v.z]).collect() };
    let unflat = |f: &[f64]| -> Vec<Vec3> {
        f.chunks(3).map(|c| Vec3 { x: c[0], y: c[1], z: c[2] }).collect()
    };
    let mut grad = vec![Vec3::ZERO; x.len()];
    let mut e = cost(&x, Some(&mut grad));
    let mut xf = flat(&x);
    let mut gf = flat(&grad);
    let mut h = vec![0.0; m * m];
    for k in 0..m {
        h[k * m + k] = 1.0;
    }
    for _ in 0..10 * m + 100 {
        let g2: f64 = gf.iter().map(|g| g * g).sum();
        if g2 < 1e-26 {
            break;
        }
        let mut dir: Vec<f64> = (0..m).map(|r| -(0..m).map(|c| h[r * m + c] * gf[c]).sum::<f64>()).collect();
        let mut slope: f64 = dir.iter().zip(&gf).map(|(d, g)| d * g).sum();
        if slope >= 0.0 {
            // Not a descent direction: start the curvature over.
            h.iter_mut().for_each(|v| *v = 0.0);
            for k in 0..m {
                h[k * m + k] = 1.0;
            }
            dir = gf.iter().map(|g| -g).collect();
            slope = -g2;
        }
        let mut step = 1.0;
        let mut moved = None;
        while step > 1e-12 {
            let trial: Vec<f64> = xf.iter().zip(&dir).map(|(a, d)| a + step * d).collect();
            let mut tg = vec![Vec3::ZERO; x.len()];
            let et = cost(&unflat(&trial), Some(&mut tg));
            if et <= e + 1e-4 * step * slope {
                moved = Some((trial, et, flat(&tg)));
                break;
            }
            step *= 0.5;
        }
        let Some((nx, ne, ng)) = moved else { break };
        let sv: Vec<f64> = nx.iter().zip(&xf).map(|(a, b)| a - b).collect();
        let yv: Vec<f64> = ng.iter().zip(&gf).map(|(a, b)| a - b).collect();
        let sy: f64 = sv.iter().zip(&yv).map(|(a, b)| a * b).sum();
        if sy > 1e-30 {
            let hy: Vec<f64> = (0..m).map(|r| (0..m).map(|c| h[r * m + c] * yv[c]).sum()).collect();
            let yhy: f64 = yv.iter().zip(&hy).map(|(a, b)| a * b).sum();
            let rho = 1.0 / sy;
            for r in 0..m {
                for c in 0..m {
                    h[r * m + c] += (1.0 + yhy * rho) * rho * sv[r] * sv[c]
                        - rho * (hy[r] * sv[c] + sv[r] * hy[c]);
                }
            }
        }
        xf = nx;
        e = ne;
        gf = ng;
    }
    let x = unflat(&xf);
    for (p, q) in pos.iter_mut().zip(&x) {
        *p = q.scale(scale);
    }
}

/// A unit direction at the domain angle to everything in `taken`.
///
/// `left` is how many directions still have to fit, which decides how far
/// around the available cone this one is rotated. `reference`, where there is
/// one, is the way the first of them leans when only the bond back to the
/// parent is spoken for.
fn direction_avoiding(taken: &[Vec3], angle: f64, left: usize, reference: Option<Vec3>) -> Vec3 {
    match taken.len() {
        // Free choice. Deterministic, so the same arrangement always gives the
        // same shape.
        0 => Vec3 { x: 0.0, y: 0.0, z: 1.0 },
        // One direction spoken for: sit at the domain angle to it, fanning
        // around it if more than one still has to fit.
        1 => {
            let axis = taken[0];
            // The part of the reference square to the axis; a reference along
            // the axis (a linear parent) gives no lean, and the helper does.
            let perp = reference
                .filter(|r| r.norm() > 0.0)
                .map(|r| {
                    let r = r.unit();
                    r - axis.scale(r.dot(axis))
                })
                .filter(|p| p.norm() > 1e-6)
                .unwrap_or_else(|| perpendicular_to(axis));
            spread(axis, perp, 0, angle, left)
        }
        // Two or more: away from their mean, half the domain angle off it, and
        // rotated so the remaining ones share the cone.
        _ => {
            let mut mean = Vec3::ZERO;
            for t in taken {
                mean += *t;
            }
            let axis = if mean.norm() > 1e-12 { -mean.unit() } else { -taken[0] };
            // The reference for the rotation is the plane the taken directions
            // already span, so the new ones come out of it rather than into it.
            let plane = taken[0].cross(taken[1]);
            let perp = if plane.norm() > 1e-12 {
                let p = plane - axis.scale(plane.dot(axis));
                if p.norm() > 1e-12 { p.unit() } else { perpendicular_to(axis) }
            } else {
                perpendicular_to(axis)
            };
            // How far off the negative mean depends on how much is already
            // spoken for. With two placed, the two remaining sit at half the
            // domain angle either side; with three, the last one is exactly
            // opposite their mean, which is what closes a tetrahedron.
            let tilt = if taken.len() >= 3 { 0.0 } else { angle * 0.5 };
            spread(axis, perp, 0, tilt, left)
        }
    }
}

fn perpendicular_to(axis: Vec3) -> Vec3 {
    let helper = if axis.x.abs() < 0.9 {
        Vec3 { x: 1.0, y: 0.0, z: 0.0 }
    } else {
        Vec3 { x: 0.0, y: 1.0, z: 0.0 }
    };
    let p = axis.cross(helper);
    if p.norm() > 1e-12 {
        p.unit()
    } else {
        Vec3 { x: 0.0, y: 1.0, z: 0.0 }
    }
}

/// A unit vector `tilt` radians off `axis`, rotated around it so that `left`
/// such vectors would share the cone evenly.
fn spread(axis: Vec3, perp: Vec3, slot: usize, tilt: f64, left: usize) -> Vec3 {
    let axis = if axis.norm() > 1e-12 { axis.unit() } else { Vec3 { x: 0.0, y: 0.0, z: 1.0 } };
    let perp = {
        let p = perp - axis.scale(perp.dot(axis));
        if p.norm() > 1e-12 { p.unit() } else { perpendicular_to(axis) }
    };
    let third = axis.cross(perp);
    let around = if left > 1 {
        std::f64::consts::TAU * slot as f64 / left as f64
    } else {
        0.0
    };
    let radial = perp.scale(around.cos()) + third.scale(around.sin());
    (axis.scale(tilt.cos()) + radial.scale(tilt.sin())).unit()
}

/// Net dipole moment, coulomb-metres, from the conformer and the bond
/// ionicities.
///
/// Each bond carries charge `ionicity * e` displaced along its own length, from
/// the less electronegative atom toward the more. Summing them *as vectors* is
/// what makes carbon dioxide non-polar and water polar, and the difference
/// between those two answers is most of what "does it dissolve" means.
pub fn dipole(arr: &Arrangement, pos: &[Vec3]) -> f64 {
    let mut total = Vec3::ZERO;
    for b in &arr.bonds {
        let (i, j) = (b.a as usize, b.b as usize);
        if i >= pos.len() || j >= pos.len() {
            continue;
        }
        let (ea, eb) = (arr.atoms[i], arr.atoms[j]);
        let q = super::analyse::ionicity(ea, eb) * crate::units::E_CHARGE;
        let (xa, xb) = (
            ea.electronegativity().unwrap_or(0.0),
            eb.electronegativity().unwrap_or(0.0),
        );
        // Points from positive to negative: toward the more electronegative
        // atom.
        let sep = pos[j] - pos[i];
        let toward = if xb >= xa { sep } else { -sep };
        total += toward.scale(q);
    }
    total.norm()
}
