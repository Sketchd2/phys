//! Relativistic frames, nested coordinates, and honest precision accounting.
//!
//! # Why there is no global coordinate system
//!
//! The engine spans 10^21 m down to 10^-15 m. Representing a femtometre offset
//! at galactic radius needs ~36 significant digits; `f64` has 15.95 and even
//! `f128` would not be enough. Fixed-point at 10^-15 m resolution over 10^21 m
//! needs a 120-bit integer per axis, which is 45 bytes of position per particle
//! — three times the rest of the state, and unusable on a GPU.
//!
//! So the engine never forms a global coordinate. Every position is an offset
//! from its parent node's centre, and any two positions are compared by walking
//! to their lowest common ancestor. The precision you get is then *exactly the
//! precision the question deserves*: two nucleons in the same nucleus share a
//! deep ancestor, so their separation is computed from small offsets and is
//! good to 10^-31 m; a nucleon and a star on the other side of the galaxy share
//! only the root, so their separation is good to ~10^5 m — which is fine,
//! because no physical process couples them more tightly than that.
//!
//! This is not a workaround. It is the same statement as the causal-gating
//! rule: things that are far apart cannot interact sharply, so they do not need
//! to be located sharply relative to one another.

use crate::math::{v3, Quat, Vec3};
use crate::units::C;

/// A vector together with a bound on its accumulated round-off.
///
/// Carried through every cross-frame computation so the engine can *prove* it
/// is not reporting a number finer than its own arithmetic supports. When
/// `err` exceeds the length scale of a physical process, the engine refuses to
/// evaluate that process at that separation and defers to a coarser tier.
#[derive(Debug, Clone, Copy)]
pub struct Located {
    pub value: Vec3,
    /// Absolute error bound in metres (1-sigma-ish; conservative).
    pub err: f64,
}

impl Located {
    pub fn exact(value: Vec3) -> Located {
        Located { value, err: 0.0 }
    }

    pub fn new(value: Vec3, err: f64) -> Located {
        Located { value, err }
    }

    /// Sum two offsets, growing the error bound by the round-off of the
    /// larger magnitude. `f64::EPSILON/2` is the unit round-off.
    pub fn add(self, o: Located) -> Located {
        let v = self.value + o.value;
        let mag = self.value.max_abs().max(o.value.max_abs()).max(v.max_abs());
        Located {
            value: v,
            err: self.err + o.err + mag * (f64::EPSILON * 0.5),
        }
    }

    pub fn sub(self, o: Located) -> Located {
        self.add(Located {
            value: -o.value,
            err: o.err,
        })
    }

    pub fn scale(self, s: f64) -> Located {
        Located {
            value: self.value.scale(s),
            err: self.err * s.abs(),
        }
    }

    /// Is this separation resolved well enough to evaluate a process whose
    /// characteristic length is `scale`?
    pub fn resolves(&self, scale: f64) -> bool {
        self.err < scale * 0.01
    }

    /// Relative precision of the separation itself.
    pub fn relative_error(&self) -> f64 {
        let n = self.value.norm();
        if n > 0.0 {
            self.err / n
        } else {
            f64::INFINITY
        }
    }
}

// ---------------------------------------------------------------------------
// Special relativity
// ---------------------------------------------------------------------------

/// Lorentz factor, clamped so that a numerically superluminal velocity (which
/// can only arise from a user injecting one) saturates instead of producing
/// NaN and poisoning the whole tree.
#[inline]
pub fn gamma(v: Vec3) -> f64 {
    let b2 = v.norm2() / (C * C);
    if b2 >= 1.0 - 1e-18 {
        1e9
    } else {
        1.0 / (1.0 - b2).sqrt()
    }
}

#[inline]
pub fn beta(v: Vec3) -> Vec3 {
    v.scale(1.0 / C)
}

/// Relativistic velocity composition: the velocity of an object measured in
/// frame S, given its velocity `u` in frame S' and S' moving at `v` in S.
///
/// Nested frames are the norm here — a nucleon inside an atom inside a cell
/// inside a planet inside a star inside the galaxy is five boosts deep — so
/// naive addition would accumulate a real error for relativistic species
/// (cosmic rays, jet material, decay products).
pub fn velocity_add(v: Vec3, u: Vec3) -> Vec3 {
    let c2 = C * C;
    let v2 = v.norm2();
    if v2 < 1e-12 {
        // Non-relativistic frame: exact to O(v^2/c^2), and avoids a 0/0.
        return v + u;
    }
    let g = gamma(v);
    let vu = v.dot(u);
    let denom = 1.0 + vu / c2;
    let par = v.scale(1.0 + (g / (g + 1.0)) * vu / c2);
    let perp = u.scale(1.0 / g);
    (par + perp).scale(1.0 / denom)
}

/// Proper-time increment for coordinate-time step `dt` at velocity `v`.
/// Every node carries its own proper time; there is no global "now".
#[inline]
pub fn proper_time_step(dt: f64, v: Vec3) -> f64 {
    dt / gamma(v)
}

/// Gravitational time dilation factor at potential `phi` (negative), to first
/// post-Newtonian order. Used near compact objects; the engine switches to a
/// full metric solver inside `r < 100 r_s` (see `solvers/gravity.rs`).
#[inline]
pub fn gravitational_dilation(phi: f64) -> f64 {
    let x = 1.0 + 2.0 * phi / (C * C);
    if x <= 1e-9 {
        1e-9f64.sqrt()
    } else {
        x.sqrt()
    }
}

/// Doppler factor for light emitted by a source with velocity `v_src`, seen
/// along unit direction `n` (source → observer).
pub fn doppler(v_src: Vec3, n: Vec3) -> f64 {
    let g = gamma(v_src);
    let bn = beta(v_src).dot(n);
    1.0 / (g * (1.0 - bn))
}

// ---------------------------------------------------------------------------
// Retarded time
// ---------------------------------------------------------------------------

/// Solve the retardation condition
///
/// ```text
///     |x_obs - x_src(t_ret)| = c (t_obs - t_ret)
/// ```
///
/// by Newton's method. `src_at` returns the source's position at a given time
/// (interpolated from its history ring, see `causal.rs`).
///
/// The obvious implementation is a fixed-point iteration — set
/// `t_ret = t_obs - d/c`, recompute `d`, repeat. That is a contraction with
/// factor `|v|/c`, which is excellent for a star (10^-3) and useless for
/// anything fast: at `v = 0.9c` it needs ~200 iterations to converge, and a
/// jet, a cosmic ray, or a decay product will hit that. Newton's method
/// converges quadratically at any speed because the derivative
///
/// ```text
///     f'(t) = c - n . v_src
/// ```
///
/// is bounded away from zero by `c(1 - |v|/c)` for any subluminal source. Six
/// iterations covers everything up to 0.999c.
pub fn retarded_time<F: Fn(f64) -> Vec3>(x_obs: Vec3, t_obs: f64, src_at: &F) -> (f64, Vec3, f64) {
    let mut t_ret = t_obs;
    let mut pos = src_at(t_ret);
    let mut dist = (x_obs - pos).norm();
    // First guess from the instantaneous distance.
    t_ret = t_obs - dist / C;
    for _ in 0..24 {
        pos = src_at(t_ret);
        let sep = x_obs - pos;
        dist = sep.norm();
        let f = dist - C * (t_obs - t_ret);
        if f.abs() <= 1e-13 * dist.max(1.0) {
            break;
        }
        // Numerical derivative of the source position, so this works for any
        // history interpolation without needing an analytic velocity.
        let h = (t_obs - t_ret).abs().max(1e-9) * 1e-6;
        let v = (src_at(t_ret + h) - src_at(t_ret - h)).scale(0.5 / h);
        let n = if dist > 0.0 { sep.scale(1.0 / dist) } else { Vec3::ZERO };
        // f'(t) = c - n.v  (the derivative of |x_obs - x_src(t)| is -n.v)
        let fp = C - n.dot(v);
        let step = if fp.abs() > 1e-30 { f / fp } else { 0.0 };
        t_ret -= step;
        if !t_ret.is_finite() {
            t_ret = t_obs - dist / C;
            break;
        }
    }
    pos = src_at(t_ret);
    dist = (x_obs - pos).norm();
    (t_ret, pos, dist)
}

/// The engine's fundamental lookahead: nothing at distance `d` can influence
/// you for at least this long. Every scheduling decision leans on it.
#[inline]
pub fn causal_delay(d: f64) -> f64 {
    d / C
}

/// Aberration: direction a moving observer sees light arriving from, given the
/// direction in the rest frame. Applied to every rendered photon.
pub fn aberrate(n_rest: Vec3, v_obs: Vec3) -> Vec3 {
    let b = beta(v_obs);
    let b2 = b.norm2();
    if b2 < 1e-24 {
        return n_rest;
    }
    let g = gamma(v_obs);
    let bn = b.dot(n_rest);
    let denom = 1.0 - bn;
    let par = b.scale(-1.0 + (g / (g + 1.0)) * bn / 1.0);
    let n = (n_rest.scale(1.0 / g) + par).scale(1.0 / denom);
    n.unit()
}

/// A rest frame: origin offset and velocity relative to the parent frame.
#[derive(Debug, Clone, Copy, Default)]
pub struct Frame {
    pub offset: Vec3,
    pub velocity: Vec3,
    /// Which way the node is pointing.
    ///
    /// A node has always had angular momentum; it had nowhere to record what
    /// that momentum had *done*. Without an orientation there is no way to draw
    /// a rotating planet, and no way to ask where a point on its surface is
    /// between one solve and the next — so rotation had to be either integrated
    /// at the rate you wanted to watch it, or not seen at all.
    pub orientation: Quat,
    /// Angular velocity, rad/s, in the node's own frame.
    pub spin_rate: Vec3,
    /// Proper time elapsed in this frame since the simulation epoch.
    pub proper_time: f64,
}

impl Frame {
    pub fn at_rest(offset: Vec3) -> Frame {
        Frame {
            offset,
            velocity: Vec3::ZERO,
            orientation: Quat::IDENTITY,
            spin_rate: Vec3::ZERO,
            proper_time: 0.0,
        }
    }

    /// Compose with the parent frame to express this frame in the grandparent.
    pub fn compose(self, parent: Frame) -> Frame {
        Frame {
            offset: parent.offset + self.offset,
            velocity: velocity_add(parent.velocity, self.velocity),
            orientation: self.orientation.then(parent.orientation).unit(),
            spin_rate: parent.spin_rate + parent.orientation.rotate(self.spin_rate),
            proper_time: self.proper_time,
        }
    }

    /// Carry the frame forward in closed form.
    ///
    /// This is what makes a slow update cadence affordable. Position under
    /// constant velocity and orientation under constant angular velocity are
    /// both exact solutions, so between two solves a node can be asked where it
    /// is at *any* instant and answer without approximation. The solver's job
    /// is to work out when the velocity and the spin rate change; carrying them
    /// forward is free.
    ///
    /// Which is the difference between a planet that turns and a planet that
    /// jumps: Earth rotates ten degrees in one signal-crossing time and Jupiter
    /// nearly twice around, so integrating orientation at the cadence the
    /// *dynamics* need would alias the rotation entirely.
    pub fn advance(&mut self, dt: f64) {
        self.proper_time += proper_time_step(dt, self.velocity);
        self.offset += self.velocity.scale(dt);
        if self.spin_rate != Vec3::ZERO {
            self.orientation = Quat::from_rate(self.spin_rate, dt).then(self.orientation).unit();
        }
    }

    /// Where a point fixed in this node's body appears, in the parent's frame.
    pub fn body_to_parent(&self, local: Vec3) -> Vec3 {
        self.offset + self.orientation.rotate(local)
    }
}

/// Convert a direction and distance to a unit direction plus light travel time.
pub fn light_delay_along(sep: Vec3) -> (Vec3, f64) {
    let d = sep.norm();
    (sep.unit(), d / C)
}

/// Approximate angular size of an object of radius `r` at distance `d`.
pub fn angular_size(r: f64, d: f64) -> f64 {
    if d <= r {
        std::f64::consts::PI
    } else {
        2.0 * (r / d).asin().abs()
    }
}

pub fn spherical_to_cart(theta: f64, phi: f64) -> Vec3 {
    v3(
        theta.sin() * phi.cos(),
        theta.sin() * phi.sin(),
        theta.cos(),
    )
}
