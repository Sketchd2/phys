//! Generated detail has to be *statistically* right, not merely conservative.
//! A cloud whose parcels conserve energy but follow the wrong velocity
//! distribution is detectable by anyone with a spectrograph.

use phys::prolong::*;
use phys::rng::{Purpose, Stream};
use phys::state::*;
use phys::units::*;

/// Sampled velocities must follow Maxwell-Boltzmann.
///
/// Checked on the moments: for a Maxwellian, `<v^4>/<v^2>^2 = 5/3` in 3D,
/// independent of temperature. That ratio is a sharp test — it fails
/// immediately for a uniform or a top-hat distribution, which both pass a
/// naive "the mean speed looks right" check.
#[test]
fn velocities_are_maxwellian() {
    let agg = Aggregate::neutral(1e30, 1e12, 1e4, Composition::primordial());
    let spec = ProlongSpec::new(50_000, Profile::Uniform, MassSpectrum::Equal, BodyKind::GasParcel);
    let (bodies, _) = prolong(&agg, spec, 4242, 0xAA, 0);
    let vbulk = total_momentum(&bodies).scale(1.0 / agg.mass);
    let n = bodies.len() as f64;
    let m2: f64 = bodies.iter().map(|b| (b.vel - vbulk).norm2()).sum::<f64>() / n;
    let m4: f64 = bodies.iter().map(|b| (b.vel - vbulk).norm2().powi(2)).sum::<f64>() / n;
    let kurtosis = m4 / (m2 * m2);
    println!("<v^4>/<v^2>^2 = {kurtosis:.4} (Maxwell-Boltzmann: 1.6667)");
    assert!((kurtosis - 5.0 / 3.0).abs() < 0.05, "not Maxwellian: {kurtosis:.4}");

    // Isotropy: no axis may be preferred.
    let vx: f64 = bodies.iter().map(|b| (b.vel.x - vbulk.x).powi(2)).sum::<f64>() / n;
    let vy: f64 = bodies.iter().map(|b| (b.vel.y - vbulk.y).powi(2)).sum::<f64>() / n;
    let vz: f64 = bodies.iter().map(|b| (b.vel.z - vbulk.z).powi(2)).sum::<f64>() / n;
    let spread = (vx.max(vy).max(vz) - vx.min(vy).min(vz)) / m2;
    println!("axis anisotropy {spread:.4}");
    assert!(spread < 0.05, "velocity field is anisotropic: {spread:.4}");
}

/// The Kroupa IMF must come out with the right slopes, or the engine makes the
/// wrong stars and therefore the wrong galaxy.
#[test]
fn imf_slope_is_kroupa() {
    let agg = Aggregate::neutral(1e5 * M_SUN, PARSEC, 30.0, Composition::solar());
    let spec = ProlongSpec::new(
        200_000,
        Profile::Plummer,
        MassSpectrum::Kroupa { min_msun: 0.08, max_msun: 60.0 },
        BodyKind::Star,
    );
    let (bodies, _) = prolong(&agg, spec, 11, 0xBB, 0);
    // The masses are rescaled to hit the total exactly, so recover the slope
    // from the *shape* rather than absolute values.
    let mut masses: Vec<f64> = bodies.iter().map(|b| b.mass).collect();
    masses.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m_lo = masses[0];
    // Count in log bins above the 0.5 Msun break, where the slope is -2.3.
    let scale = masses[masses.len() / 2];
    let bins = |lo: f64, hi: f64| masses.iter().filter(|&&m| m >= lo && m < hi).count() as f64;
    let a = bins(4.0 * scale, 8.0 * scale);
    let b = bins(8.0 * scale, 16.0 * scale);
    assert!(a > 20.0 && b > 5.0, "not enough massive stars: {a}, {b}");
    // dN/dm ~ m^-2.3  =>  N per octave ~ m^-1.3
    let slope = (b / a).ln() / 2f64.ln();
    println!("high-mass slope: N per octave ~ m^{slope:.2} (Kroupa: -1.3)");
    assert!((slope + 1.3).abs() < 0.45, "IMF slope {slope:.2}");
    assert!(m_lo > 0.0);
}

/// A Plummer sphere must actually have a Plummer density profile.
#[test]
fn plummer_profile_is_correct() {
    let agg = Aggregate::neutral(1e5 * M_SUN, PARSEC, 30.0, Composition::solar());
    let spec = ProlongSpec::new(100_000, Profile::Plummer, MassSpectrum::Equal, BodyKind::Star);
    let (bodies, report) = prolong(&agg, spec, 13, 0xCC, 0);
    let com = bodies.iter().fold(phys::math::Vec3::ZERO, |a, b| a + b.pos).scale(1.0 / bodies.len() as f64);
    let mut r: Vec<f64> = bodies.iter().map(|b| (b.pos - com).norm()).collect();
    r.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // For a Plummer sphere, M(<r)/M = r^3 / (r^2 + a^2)^(3/2). Solving for the
    // quantiles: r50 = 1.3048 a and r25 = 0.8112 a, so the ratio is 1.6086 —
    // a shape test that is independent of the scale the sampler was asked for.
    let r_half = r[r.len() / 2];
    let r_quarter = r[r.len() / 4];
    let ratio = r_half / r_quarter;
    println!("r50/r25 = {ratio:.4} (Plummer: 1.6086)");
    assert!((ratio - 1.6086).abs() < 0.03, "profile ratio {ratio:.4}");
    assert!(report.realised_radius > 0.0);
}

/// Poisson counting statistics: a detector must show sqrt(N) noise, because
/// that is what a real detector shows and it is how an observer would catch a
/// simulation faking its photon counts.
#[test]
fn counting_noise_is_poisson() {
    let mut s = Stream::at(1, 1, 0, Purpose::PhotonEmission);
    let lambda = 25.0;
    let n = 40000;
    let samples: Vec<f64> = (0..n).map(|_| s.poisson(lambda) as f64).collect();
    let mean = samples.iter().sum::<f64>() / n as f64;
    let var = samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    println!("Poisson({lambda}): mean {mean:.3}, variance {var:.3}");
    assert!((mean - lambda).abs() < 0.2, "mean {mean:.3}");
    assert!((var - lambda).abs() / lambda < 0.06, "variance {var:.3} should equal the mean");
}

/// Uniform draws must fill [0,1) without bias in any bin.
#[test]
fn uniform_draws_are_uniform() {
    let mut s = Stream::at(7, 7, 0, Purpose::Positions);
    let n = 400_000;
    let bins = 64;
    let mut hist = vec![0u32; bins];
    for _ in 0..n {
        let u = s.uniform();
        assert!((0.0..1.0).contains(&u));
        hist[(u * bins as f64) as usize] += 1;
    }
    let expected = n as f64 / bins as f64;
    let chi2: f64 = hist.iter().map(|&c| (c as f64 - expected).powi(2) / expected).sum();
    println!("chi^2 = {chi2:.1} for {} bins (expect ~{})", bins, bins - 1);
    // 63 dof: the 0.999 critical value is ~112.
    assert!(chi2 < 112.0, "uniform draws are biased: chi^2 = {chi2:.1}");
}

/// Isotropic directions must have no preferred axis.
#[test]
fn directions_are_isotropic() {
    let mut s = Stream::at(3, 3, 0, Purpose::Spin);
    let n = 200_000;
    let mut sum = phys::math::Vec3::ZERO;
    let mut zz = 0.0;
    for _ in 0..n {
        let d = s.direction();
        assert!((d.norm() - 1.0).abs() < 1e-9);
        sum += d;
        zz += d.z * d.z;
    }
    let drift = sum.norm() / n as f64;
    println!("mean direction magnitude {drift:.5} (should be ~1/sqrt(n) = {:.5})", 1.0 / (n as f64).sqrt());
    assert!(drift < 5.0 / (n as f64).sqrt(), "directions are biased");
    let mean_zz = zz / n as f64;
    assert!((mean_zz - 1.0 / 3.0).abs() < 0.01, "<z^2> = {mean_zz:.4}, expected 1/3");
}

/// Composition scatter must vary the children while leaving the mass-weighted
/// mean exactly equal to the parent's.
#[test]
fn composition_scatter_preserves_the_mean() {
    let agg = Aggregate::neutral(1e30, 1e10, 1e4, Composition::solar());
    let spec = ProlongSpec {
        composition_scatter: 0.3,
        ..ProlongSpec::new(5000, Profile::Uniform, MassSpectrum::Equal, BodyKind::GasParcel)
    };
    let (bodies, _) = prolong(&agg, spec, 21, 0xDD, 0);
    let total: f64 = bodies.iter().map(|b| b.mass).sum();
    for s in Species::ALL {
        let mean: f64 = bodies.iter().map(|b| b.mass * b.composition.get(s)).sum::<f64>() / total;
        let want = agg.composition.get(s);
        let err = (mean - want).abs() / want.max(1e-12);
        assert!(err < 1e-9, "{}: mean fraction {mean:.9e} vs {want:.9e}", s.name());
    }
    // And they must actually differ from one another.
    let spread: f64 = {
        let vals: Vec<f64> = bodies.iter().map(|b| b.composition.get(Species::Iron)).collect();
        let m = vals.iter().sum::<f64>() / vals.len() as f64;
        (vals.iter().map(|v| (v - m).powi(2)).sum::<f64>() / vals.len() as f64).sqrt() / m
    };
    println!("iron fraction scatter: {:.1}%", spread * 100.0);
    assert!(spread > 0.05, "composition scatter had no effect");
}
