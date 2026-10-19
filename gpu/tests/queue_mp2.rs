//! The queue's MP2 tasks (`mp2` and `unlike` plan lines, `pairs::Executor`): a plan names them, a worker with
//! the same molecule shapes computes them, and a worker holding another shape refuses. The strong check is that
//! `mp2 water 0` reproduces the first line of a stored MP2 pairs file, geometry and energies, since the draw is a
//! function of the index alone. Skipped where the stored data is not.

use phys::pairs::{plan, Executor};

#[test]
fn the_queues_mp2_tasks_reproduce_stored_pairs_and_are_planned() {
    let stored = "../runs/e8d/core/pairs-water-gpu-mp2.txt";
    let (Ok(text), Ok(state)) = (std::fs::read_to_string(stored), std::fs::read_to_string("../runs/e8d/core/grow-water.state")) else {
        eprintln!("no stored pairs; skipped");
        return;
    };
    // A scratch directory holding the shape that data was made at.
    let dir = std::env::temp_dir().join("phys-queue-mp2-test");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("grow-water.state"), &state).unwrap();
    if let Ok(a) = std::fs::read_to_string("../runs/e8e/grow/grow-ammonia.state") {
        std::fs::write(dir.join("grow-ammonia.state"), a).unwrap();
    }
    std::env::set_current_dir(&dir).unwrap();
    // The plan.
    let out = dir.join("plan-out.txt");
    let _ = std::fs::remove_file(&out);
    let p = plan(&format!("mp2 water 0 3 any {}\nmp2 water 3 5 gpu {} wall=1.0,2.3\n", out.display(), dir.join("plan-wall.txt").display())).expect("a plan");
    assert_eq!(p.len(), 5);
    assert!(p[0].spec.starts_with("mp2 water ") && p[0].spec.ends_with(" 0"), "{}", p[0].spec);
    assert!(p[3].spec.ends_with(" 3 1 2.3"), "{}", p[3].spec);
    assert_eq!(p[0].kind, "mp2 water");
    assert!(plan("mp2 water 0 3 gpu out.txt wall=oops").is_err());
    assert!(plan("unlike water nothing 0 3 gpu out.txt").is_err());
    if dir.join("grow-ammonia.state").exists() {
        let u = plan(&format!("unlike water ammonia 0 2 gpu {}\n", dir.join("plan-unlike.txt").display())).expect("an unlike plan");
        assert_eq!(u.len(), 2);
        assert!(u[0].spec.starts_with("unlike water ammonia "), "{}", u[0].spec);
    }
    // A worker holding a different shape refuses.
    let mut exec = Executor::default();
    let wrong = p[0].spec.replacen(" mp2 water ", " mp2 water ", 1);
    let spec_bad_shape = {
        let mut w: Vec<&str> = wrong.split_whitespace().collect();
        w[2] = "deadbeef";
        w.join(" ")
    };
    assert!(exec.run(&spec_bad_shape).unwrap_err().contains("not the server's shape"));
    // And the first stored pair is reproduced.
    if phys_gpu::install_fock().is_err() {
        eprintln!("no GPU Fock; the CPU does it");
    }
    let line = exec.run(&p[0].spec).expect("pair 0");
    let stored_line = text.lines().find(|l| !l.starts_with('#')).expect("a stored pair");
    let num = |l: &str| -> Vec<f64> { l.replace('|', " ").split_whitespace().map(|x| x.parse::<f64>().unwrap_or(f64::NAN)).collect() };
    let (a, b) = (num(&line), num(stored_line));
    assert_eq!(a.len(), b.len(), "the same number of fields");
    let mut worst_geometry: f64 = 0.0;
    for k in 4..a.len() {
        if !(a[k].is_nan() || b[k].is_nan()) {
            worst_geometry = worst_geometry.max((a[k] - b[k]).abs());
        }
    }
    println!("pair 0: energy {:.10e} against stored {:.10e}; Hartree-Fock {:.10e} against {:.10e}; worst geometry difference {worst_geometry:.1e}", a[2], b[2], a[3], b[3]);
    assert!(worst_geometry < 2e-8, "the drawn geometry is the stored one");
    assert!((a[2] - b[2]).abs() < 5e-8 && (a[3] - b[3]).abs() < 5e-8, "the energies are the stored ones");
}
