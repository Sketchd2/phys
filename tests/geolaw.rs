//! The angle-following cloud's pieces: its text, the pair function that takes
//! charges and widths as arguments agreeing with the table's, and the
//! correction a fit takes off.

use phys::geolaw::{correction, x_of_law, GeoCloud};
use phys::liquid::{extra_sites_from_text, SiteLaw, SiteSite};
use phys::math::Vec3;

fn law() -> (SiteSite, Vec<f64>, String) {
    let text = std::fs::read_to_string("law-water-mp2-cloud-final.txt").expect("the cloud law");
    let law = SiteSite::from_text(&text).expect("a readable law");
    let alpha = SiteSite::alpha_from_text(&text, law.charge.len());
    (law, alpha, text)
}

#[test]
fn the_text_round_trips() {
    let mut g = GeoCloud::constant([6.0, 1.0], [0.5, -0.3, -1.2, 0.3, 1.6, 5.4], 1.8165);
    g.coef[1][1] = 0.07;
    g.coef[4][2] = -0.2;
    let back = GeoCloud::from_text(&g.to_text()).expect("geo lines");
    assert_eq!(back, g);
    assert!(GeoCloud::from_text("charges 1 2 3").is_none());
}

#[test]
fn the_pair_function_with_arguments_is_the_table_one_and_its_gradients_are_right() {
    let (law, _, _) = law();
    let types = law.charge.len();
    for ta in 0..types {
        for tb in 0..types {
            for r in [2.0, 3.7, 6.0] {
                let (u, du) = law.site_pair(ta, tb, r);
                let (sa, sb) = (law.sigma.get(ta).copied().unwrap_or(0.0), law.sigma.get(tb).copied().unwrap_or(0.0));
                let (v, dv, _) = law.site_pair_with(ta, tb, r, law.charge[ta], law.charge[tb], sa, sb);
                assert!((u - v).abs() <= 1e-13 * u.abs().max(1e-3), "energy {ta} {tb} {r}: {u} against {v}");
                assert!((du - dv).abs() <= 1e-12 * du.abs().max(1e-3), "slope {ta} {tb} {r}: {du} against {dv}");
            }
        }
    }
    // The derivatives with respect to charges and widths, against differences.
    let (ta, tb, r) = (3usize, 4usize, 3.1);
    let (qa, qb, sa, sb) = (-5.4, -0.45, 0.8, 0.3);
    let (_, _, g) = law.site_pair_with(ta, tb, r, qa, qb, sa, sb);
    let at = |qa: f64, qb: f64, sa: f64, sb: f64| law.site_pair_with(ta, tb, r, qa, qb, sa, sb).0;
    let h = 1e-6;
    let numeric = [
        (at(qa + h, qb, sa, sb) - at(qa - h, qb, sa, sb)) / (2.0 * h),
        (at(qa, qb + h, sa, sb) - at(qa, qb - h, sa, sb)) / (2.0 * h),
        (at(qa, qb, sa + h, sb) - at(qa, qb, sa - h, sb)) / (2.0 * h),
        (at(qa, qb, sa, sb + h) - at(qa, qb, sa, sb - h)) / (2.0 * h),
    ];
    for k in 0..4 {
        assert!((g[k] - numeric[k]).abs() < 1e-7 * (1.0 + numeric[k].abs()), "gradient {k}: {} against {}", g[k], numeric[k]);
    }
}

#[test]
fn the_correction_is_nothing_at_the_reference_angle_and_something_away_from_it() {
    let (law, alpha, text) = law();
    let bis = SiteSite::bisector_from_text(&text).expect("a bisector").1;
    let frame = extra_sites_from_text(&text);
    let r = frame.iter().find(|s| s.ty == 4).map(|s| s.frame[0]).expect("hydrogen clouds");
    let mut g = GeoCloud::constant([law.charge[0], law.charge[1]], x_of_law(&law, bis), 104.08f64.to_radians());
    g.coef[4][1] = 0.2;
    g.coef[0][1] = 0.1;
    let molecule = |deg: f64, at: Vec3| -> Vec<(Vec3, usize)> {
        let th = deg.to_radians();
        vec![(at, 0), (at + Vec3 { x: r, y: 0.0, z: 0.0 }, 1), (at + Vec3 { x: r * th.cos(), y: r * th.sin(), z: 0.0 }, 1)]
    };
    let b0 = Vec3 { x: 1.0, y: 1.0, z: 5.5 };
    let c_ref = correction(&g, &frame, &alpha, &molecule(104.08, Vec3::ZERO), &molecule(104.08, b0));
    assert!(c_ref.abs() < 1e-12, "no correction at the reference angle: {c_ref:e}");
    let c = correction(&g, &frame, &alpha, &molecule(112.0, Vec3::ZERO), &molecule(98.0, b0));
    assert!(c.abs() > 1e-6, "a correction away from it: {c:e}");
}
