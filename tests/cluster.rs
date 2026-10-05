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

/// How much of the table stays in memory: half of what is free unless a
/// fraction in range is named, and a setting out of range is not trusted.
#[test]
fn the_tables_share_of_memory_is_half_unless_told_otherwise() {
    use phys::electrons::scf::table_budget;
    let free = Some(64_000_000_000u64);
    assert_eq!(table_budget(free, None), 32_000_000_000, "the default is half, as before");
    assert_eq!(table_budget(free, Some("0.5")), 32_000_000_000);
    assert_eq!(table_budget(free, Some("0.8")), 51_200_000_000);
    assert_eq!(table_budget(free, Some(" 0.25 ")), 16_000_000_000);
    for bad in ["0.99", "1", "0", "-1", "lots", ""] {
        assert_eq!(table_budget(free, Some(bad)), 32_000_000_000, "{bad:?} was taken");
    }
    assert_eq!(table_budget(None, Some("0.8")), usize::MAX, "no limit where free memory is not known");
    assert_eq!(table_budget(Some(7), None), 3, "an odd count rounds down as before");
}

/// Pieces done on other machines come together in one file: what is missing is
/// added, what is there is not replaced, and a file begun on another cluster
/// is refused.
#[test]
fn the_pieces_of_a_cluster_come_together() {
    use phys::pairs::merge_cluster_files;
    let d = scratch("merge");
    let path = |n: &str| d.join(n).to_string_lossy().to_string();
    let geometry = vec!["# mol 3 1 2 3".to_string(), "# mol 9 4 5 6".to_string()];
    let head = "# cluster\n# mol 3 1 2 3\n# mol 9 4 5 6\n";
    std::fs::write(path("main.txt"), format!("{head}pair 0 1 3 9 -1.0e-3 -2.0e-3\nfield 0 -10.5e0 -10.6e0 12 100\n")).unwrap();
    std::fs::write(path("a.txt"), format!("{head}field 0 -99.0e0 -99.0e0 1 1\nfield 1 -5.25e0 -5.3e0 9 50\n")).unwrap();
    std::fs::write(path("b.txt"), format!("{head}field 1 -77.0e0 -77.0e0 1 1\nfield 2 -5.5e0 -5.6e0 9 60\ncluster 2 -1.0e-3 -2.0e-3\n")).unwrap();
    let added = merge_cluster_files(&path("main.txt"), &[path("a.txt"), path("b.txt")], &geometry);
    assert_eq!(added, 3, "field 1 from a, field 2 and the cluster from b");
    let merged = std::fs::read_to_string(path("main.txt")).unwrap();
    assert!(merged.contains("field 0 -10.5e0"), "a result already there was replaced");
    assert!(!merged.contains("-99.0e0") && !merged.contains("-77.0e0"), "a duplicate was added");
    assert!(merged.contains("field 1 -5.25e0") && merged.contains("field 2 -5.5e0") && merged.contains("cluster 2 "));
    assert_eq!(merge_cluster_files(&path("main.txt"), &[path("a.txt"), path("b.txt")], &geometry), 0, "a second merge adds nothing");
    std::fs::write(path("other.txt"), "# cluster\n# mol 3 1 2 3\n# mol 8 4 5 6\nfield 3 -1.0e0 -1.0e0 1 1\n").unwrap();
    let refused = std::panic::catch_unwind(|| merge_cluster_files(&path_of(&d, "main.txt"), &[path_of(&d, "other.txt")], &["# mol 3 1 2 3".to_string(), "# mol 9 4 5 6".to_string()]));
    assert!(refused.is_err(), "a file from another cluster was merged");
    assert!(!std::fs::read_to_string(path("main.txt")).unwrap().contains("field 3"));
}

fn path_of(d: &std::path::Path, n: &str) -> String {
    d.join(n).to_string_lossy().to_string()
}
