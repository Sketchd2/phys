//! The cluster run's pieces that need no solve: how a cluster is cut, what a
//! finished field leaves behind, and what a restart clears away.

use phys::pairs::{remove_stale_spill_in, FieldLog, Snapshot};

fn scratch(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("phys-cluster-test-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A cluster is one piece: the centre first, then the others nearest first,
/// each a different molecule of the snapshot, each moved to its image about
/// the centre without being bent.
#[test]
fn a_cluster_is_cut_in_one_piece_from_whole_molecules() {
    let snap = Snapshot::load("bulk-water-298.snap", 3);
    let cluster = snap.cluster(5, 6);
    assert_eq!(cluster.len(), 6);
    assert_eq!(cluster[0].0, 5);
    let ids: std::collections::HashSet<usize> = cluster.iter().map(|c| c.0).collect();
    assert_eq!(ids.len(), 6, "a molecule appeared twice");
    let centroid = |m: &[phys::math::Vec3]| m.iter().fold(phys::math::Vec3::ZERO, |a, p| a + *p).scale(1.0 / 3.0);
    let c0 = centroid(&cluster[0].1);
    let dist: Vec<f64> = cluster.iter().map(|c| (centroid(&c.1) - c0).norm()).collect();
    assert!(dist.windows(2).all(|w| w[0] <= w[1] + 1e-12), "not nearest first: {dist:?}");
    assert!(dist[5] < 10.5, "the sixth is {} bohr away, which is not a neighbour", dist[5]);
    // Whole molecules: the same internal distances as every other water in the snapshot.
    let oh = (cluster[0].1[0] - cluster[0].1[1]).norm();
    for (_, m) in &cluster {
        assert!(((m[0] - m[1]).norm() - oh).abs() < 1e-6, "a molecule was bent");
    }
    // Nothing is ever cut in two by the box: the sixth's image is near the first, not a box away.
    assert!(cluster.iter().all(|(_, m)| (centroid(m) - c0).norm() < 20.0));
}

/// What a finished field writes, a restart reads back exactly.
#[test]
fn a_finished_field_is_read_back_exactly() {
    let d = scratch("fields");
    let path = d.join("c.txt").to_string_lossy().to_string();
    std::fs::write(&path, "# a header\n# mol 3 1 2 3\npair 0 1 3 4 -1.5e-3 -1.2e-3\n").unwrap();
    let mut log = FieldLog::load(&path);
    assert!(log.done.is_empty());
    let (p, r) = (-5.678_901_234_567_89e2, -5.678_901_234_567_91e2);
    log.record(0, p, r, 13, 4648.4);
    log.record(3, 1.0 / 3.0, -2.0 / 3.0, 9, 10.0);
    let again = FieldLog::load(&path);
    assert_eq!(again.done.len(), 2, "the pair line is not a field");
    assert_eq!(again.done[&0], (p, r), "the energy did not come back bit for bit");
    assert_eq!(again.done[&3], (1.0 / 3.0, -2.0 / 3.0));
}

/// A restart clears the spill of a process that is gone, and only that: not a
/// live process's, not one it cannot tell about, not this process's own, not a
/// file that is not a spill.
#[test]
fn a_restart_clears_a_dead_runs_spill_and_nothing_else() {
    let d = scratch("spill");
    let me = std::process::id();
    let files = [
        (d.join("phys-fit-222-0.bin"), true),
        (d.join("phys-fit-222-1.bin"), true),
        (d.join("phys-fit-111-0.bin"), false),
        (d.join("phys-fit-333-0.bin"), false),
        (d.join(format!("phys-fit-{me}-7.bin")), false),
        (d.join("notes.txt"), false),
        (d.join("phys-fit-abc-0.bin"), false),
    ];
    for (f, _) in &files {
        std::fs::write(f, b"x").unwrap();
    }
    // 111 is running, 222 is not, and nothing can be said of 333.
    remove_stale_spill_in(&d, |pid| match pid {
        111 => Some(true),
        222 => Some(false),
        _ => None,
    });
    for (f, removed) in &files {
        assert_eq!(!f.exists(), *removed, "{} was {}", f.display(), if *removed { "left" } else { "removed" });
    }
}
