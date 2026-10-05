//! The work queue (`src/queue.rs`): tasks reach exactly one result each, a
//! machine that vanishes costs only its lease, and the files are the state.

use phys::queue::{serve, work, Ended, Task, WorkerConfig};
use std::net::TcpListener;
use std::time::Duration;

fn scratch(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("phys-queue-test-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn tasks(out: &str, class: &str, range: std::ops::Range<usize>) -> Vec<Task> {
    range.map(|k| Task { out: out.to_string(), index: k, class: class.to_string(), spec: format!("fake {k}"), header: "# header".into(), kind: "fake".into(), ..Default::default() }).collect()
}

fn config(addr: &str, name: &str, class: &str) -> WorkerConfig {
    WorkerConfig { server: addr.to_string(), token: "t".into(), name: name.into(), class: class.into(), mem_mb: Some(0), max_jobs: None, patience: Duration::from_secs(5) }
}

fn lines(path: &str) -> Vec<String> {
    std::fs::read_to_string(path).unwrap().lines().map(|l| l.to_string()).collect()
}

/// Two machines of one class and one of another share a plan: every task of a
/// class is done once, by a machine of that class, and the file holds one
/// header and one line a task.
#[test]
fn every_task_is_done_once_by_a_machine_of_its_class() {
    let d = scratch("share");
    let (cpu, gpu) = (d.join("cpu.txt").to_string_lossy().to_string(), d.join("gpu.txt").to_string_lossy().to_string());
    let mut all = tasks(&cpu, "cpu", 0..12);
    all.extend(tasks(&gpu, "gpu", 5..9));
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let addr = listener.local_addr().unwrap().to_string();
    let server = std::thread::spawn(move || serve(listener, all, Duration::from_secs(60), "t", None, Duration::from_secs(3)).unwrap());
    let workers: Vec<_> = [("a", "cpu"), ("b", "cpu"), ("c", "gpu")]
        .into_iter()
        .map(|(name, class)| {
            let cfg = config(&addr, name, class);
            std::thread::spawn(move || {
                let end = work(&cfg, |spec| {
                    let k: usize = spec.strip_prefix("fake ").unwrap().parse().unwrap();
                    std::thread::sleep(Duration::from_millis(20));
                    Ok(format!("{k} done-by-{} | x", cfg.class))
                });
                assert_eq!(end, Ended::NothingLeft);
            })
        })
        .collect();
    for w in workers {
        w.join().unwrap();
    }
    assert_eq!(server.join().unwrap(), 16);
    let cpu_lines = lines(&cpu);
    assert_eq!(cpu_lines[0], "# header");
    let mut got: Vec<usize> = cpu_lines[1..].iter().map(|l| l.split_whitespace().next().unwrap().parse().unwrap()).collect();
    got.sort();
    assert_eq!(got, (0..12).collect::<Vec<_>>());
    assert!(cpu_lines[1..].iter().all(|l| l.contains("done-by-cpu")));
    let gpu_lines = lines(&gpu);
    assert_eq!(gpu_lines.len(), 5);
    assert!(gpu_lines[1..].iter().all(|l| l.contains("done-by-gpu")));
}

/// A machine that takes a task and is never heard from loses it when its lease
/// runs out; the machine that finishes it first is the one kept, and the late
/// one is told `DUP` and writes nothing.
#[test]
fn a_lost_lease_is_handed_on_and_a_late_result_is_dropped() {
    use std::io::{Read, Write};
    let d = scratch("lease");
    let out = d.join("o.txt").to_string_lossy().to_string();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let addr = listener.local_addr().unwrap().to_string();
    let t = tasks(&out, "cpu", 0..1);
    let server = std::thread::spawn(move || serve(listener, t, Duration::from_millis(400), "t", None, Duration::from_secs(3)).unwrap());
    let ask = |lines: &[&str]| {
        let mut s = std::net::TcpStream::connect(&addr).unwrap();
        for l in lines {
            writeln!(s, "{l}").unwrap();
        }
        let mut r = String::new();
        s.read_to_string(&mut r).unwrap();
        r
    };
    assert!(ask(&["t GET vanished cpu 0 4"]).starts_with("JOB "));
    // Held, so a second machine is told to wait.
    assert!(ask(&["t GET other cpu 0 4"]).starts_with("WAIT"));
    std::thread::sleep(Duration::from_millis(500));
    assert!(ask(&["t GET other cpu 0 4"]).starts_with("JOB "), "the lease had run out");
    assert_eq!(ask(&[&format!("t PUT {out} 0 other 1.0 0"), "0 first | x"]), "OK\n");
    assert_eq!(ask(&[&format!("t PUT {out} 0 vanished 9.0 0"), "0 second | x"]), "DUP\n");
    assert_eq!(ask(&["wrong GET a cpu 0 4"]), "ERR token\n");
    server.join().unwrap();
    assert_eq!(lines(&out), vec!["# header", "0 first | x"]);
}

/// A line that is not the task's is refused and leaves the task open.
#[test]
fn a_result_for_another_index_is_refused() {
    use std::io::{Read, Write};
    let d = scratch("wrong");
    let out = d.join("o.txt").to_string_lossy().to_string();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let addr = listener.local_addr().unwrap().to_string();
    let t = tasks(&out, "cpu", 3..4);
    let server = std::thread::spawn(move || serve(listener, t, Duration::from_secs(60), "t", None, Duration::from_secs(3)).unwrap());
    let ask = |lines: &[&str]| {
        let mut s = std::net::TcpStream::connect(&addr).unwrap();
        for l in lines {
            writeln!(s, "{l}").unwrap();
        }
        let mut r = String::new();
        s.read_to_string(&mut r).unwrap();
        r
    };
    assert!(ask(&["t GET m cpu 0 4"]).starts_with("JOB "));
    assert!(ask(&[&format!("t PUT {out} 3 m 1.0 0"), "4 wrong | x"]).starts_with("ERR"));
    assert!(ask(&[&format!("t PUT {out} 3 m 1.0 0"), "# comment"]).starts_with("ERR"));
    assert_eq!(ask(&[&format!("t PUT {out} 3 m 1.0 0"), "3 right | x"]), "OK\n");
    server.join().unwrap();
    assert_eq!(lines(&out), vec!["# header", "3 right | x"]);
}

/// A task that fails on three machines is not handed out a fourth time, and
/// the server, with nothing else left, stops rather than waiting for it.
#[test]
fn a_task_that_keeps_failing_is_left_alone() {
    let d = scratch("fail");
    let out = d.join("o.txt").to_string_lossy().to_string();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let addr = listener.local_addr().unwrap().to_string();
    let mut t = tasks(&out, "cpu", 0..2);
    t[0].spec = "bad".into();
    let server = std::thread::spawn(move || serve(listener, t, Duration::from_secs(60), "t", None, Duration::from_secs(3)).unwrap());
    let cfg = config(&addr, "w", "cpu");
    let mut tries = 0;
    let end = work(&cfg, |spec| {
        if spec == "bad" {
            tries += 1;
            if tries == 2 {
                panic!("a field did not converge");
            }
            Err("it did not converge".into())
        } else {
            Ok("1 fine | x".into())
        }
    });
    assert_eq!(end, Ended::NothingLeft);
    assert_eq!(tries, 3, "tried the failing task three times and no more");
    assert_eq!(server.join().unwrap(), 1);
    assert_eq!(lines(&out), vec!["# header", "1 fine | x"]);
}

/// A worker stops at its job limit, and a server that is not there is given
/// up on after its patience.
#[test]
fn a_worker_stops_at_its_limit_and_gives_up_on_a_missing_server() {
    let d = scratch("limit");
    let out = d.join("o.txt").to_string_lossy().to_string();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let addr = listener.local_addr().unwrap().to_string();
    let t = tasks(&out, "cpu", 0..5);
    let server = std::thread::spawn(move || serve(listener, t, Duration::from_secs(60), "t", None, Duration::from_secs(3)).unwrap());
    let mut cfg = config(&addr, "w", "cpu");
    cfg.max_jobs = Some(2);
    assert_eq!(work(&cfg, |s| Ok(format!("{} ok | x", s.strip_prefix("fake ").unwrap()))), Ended::Limit);
    cfg.max_jobs = None;
    assert_eq!(work(&cfg, |s| Ok(format!("{} ok | x", s.strip_prefix("fake ").unwrap()))), Ended::NothingLeft);
    server.join().unwrap();
    let nothing = TcpListener::bind("127.0.0.1:0").unwrap();
    let gone = nothing.local_addr().unwrap().to_string();
    drop(nothing);
    let mut cfg = config(&gone, "w", "cpu");
    cfg.patience = Duration::from_millis(1);
    let before = std::env::current_dir().unwrap();
    std::env::set_current_dir(&d).unwrap();
    assert_eq!(work(&cfg, |_| Ok("0 x | y".into())), Ended::ServerGone);
    std::env::set_current_dir(before).unwrap();
}

/// A plan names the work, leaves out what a file already holds (whatever order
/// it was written in), and refuses a range that overlaps another's.
#[test]
fn a_plan_leaves_out_what_is_done_and_refuses_overlaps() {
    let d = scratch("plan");
    let out = d.join("p.txt").to_string_lossy().to_string();
    std::fs::write(&out, "# header\n3 0 0 0 | x\n1 0 0 0 | x\n").unwrap();
    let t = phys::pairs::plan(&format!("# a plan\n\npair water 0 6 element cpu {out}\n")).unwrap();
    let indices: Vec<usize> = t.iter().map(|t| t.index).collect();
    assert_eq!(indices, vec![0, 2, 4, 5]);
    assert_eq!(t[0].spec, "pair water element 0");
    assert_eq!(t[0].class, "cpu");
    let err = phys::pairs::plan(&format!("pair water 0 6 element cpu {out}\npair water 5 9 element cpu {out}\n")).unwrap_err();
    assert!(err.contains("plan line 2") && err.contains("already given"), "{err}");
    assert!(phys::pairs::plan("pair water 0 6 element tpu x.txt\n").unwrap_err().contains("class"));
    assert!(phys::pairs::plan("pair helium 0 6 element cpu x.txt\n").unwrap_err().contains("helium"));
    assert!(phys::pairs::plan("pair water zero 6 element cpu x.txt\n").unwrap_err().contains("not a number"));
    assert!(phys::pairs::plan("jump water\n").unwrap_err().contains("expected"));
}

fn ask(addr: &str, lines: &[&str]) -> String {
    use std::io::{Read, Write};
    let mut s = std::net::TcpStream::connect(addr).unwrap();
    for l in lines {
        writeln!(s, "{l}").unwrap();
    }
    let mut r = String::new();
    s.read_to_string(&mut r).unwrap();
    r
}

fn sized(out: &str, class: &str, k: usize, weight: f64, need_mb: Option<u64>, kind: &str) -> Task {
    Task { out: out.to_string(), index: k, class: class.to_string(), spec: format!("fake {k}"), header: "# header".into(), kind: kind.into(), weight, need_mb }
}

fn start(tasks: Vec<Task>, lease_ms: u64) -> (String, std::thread::JoinHandle<usize>) {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let addr = listener.local_addr().unwrap().to_string();
    let h = std::thread::spawn(move || serve(listener, tasks, Duration::from_millis(lease_ms), "t", None, Duration::from_secs(3)).unwrap());
    (addr, h)
}

/// A task goes only to a machine with the memory it needs, the machine is told
/// when everything left is too big for it, and a task that has no stated need
/// takes the one the first of its kind was measured to use.
#[test]
fn a_task_goes_only_to_a_machine_it_fits_in() {
    let d = scratch("fit");
    let out = d.join("o.txt").to_string_lossy().to_string();
    let (addr, server) = start(vec![sized(&out, "cpu", 0, 1.0, Some(3000), "big"), sized(&out, "cpu", 1, 1.0, Some(500), "small"), sized(&out, "cpu", 2, 1.0, None, "learned"), sized(&out, "cpu", 3, 1.0, None, "learned")], 60_000);
    // 1000 MB free: the small one, then the first of the unmeasured kind.
    assert!(ask(&addr, &["t GET little cpu 1000 4"]).starts_with("JOB "));
    assert!(ask(&addr, &["t GET little cpu 1000 4"]).starts_with("JOB "), "nothing is known yet of the kind, so it is tried");
    // That one used 2000 MB: the other of its kind no longer fits 1000, and
    // the 3000 MB task never did.
    assert_eq!(ask(&addr, &[&format!("t PUT {out} 2 little 10 2000"), "2 b | x"]), "OK
");
    assert_eq!(ask(&addr, &[&format!("t PUT {out} 1 little 10 400"), "1 a | x"]), "OK
");
    assert_eq!(ask(&addr, &["t GET little cpu 1000 4"]), "NOFIT 3000
");
    // A machine with 4000 MB takes both, the learned one first.
    assert!(ask(&addr, &["t GET large cpu 4000 4"]).starts_with("JOB "));
    assert!(ask(&addr, &["t GET large cpu 4000 4"]).starts_with("JOB "));
    assert_eq!(ask(&addr, &["t GET large cpu 4000 4"]), "WAIT 1\n");
    assert_eq!(ask(&addr, &[&format!("t PUT {out} 0 large 10 100"), "0 c | x"]), "OK\n");
    assert_eq!(ask(&addr, &[&format!("t PUT {out} 3 large 10 100"), "3 d | x"]), "OK\n");
    server.join().unwrap();
}

/// What the machines have shown decides who gets the heavy tasks: the fastest
/// takes the heaviest, the slowest the lightest, and equal weights go in order.
#[test]
fn heavy_tasks_go_to_fast_machines_and_light_ones_to_slow() {
    let d = scratch("speed");
    let out = d.join("o.txt").to_string_lossy().to_string();
    let mut t = vec![sized(&out, "warm", 100, 1.0, None, "warm"), sized(&out, "warm", 101, 1.0, None, "warm")];
    for (k, w) in [1.0, 5.0, 3.0, 3.0].into_iter().enumerate() {
        t.push(sized(&out, "cpu", k, w, None, "work"));
    }
    let (addr, server) = start(t, 60_000);
    // Both machines are measured on the warm-up: one takes 1 s, the other 50.
    assert!(ask(&addr, &["t GET fast warm 0 4"]).starts_with("JOB "));
    assert!(ask(&addr, &["t GET slow warm 0 4"]).starts_with("JOB "));
    assert_eq!(ask(&addr, &[&format!("t PUT {out} 100 fast 1 0"), "100 w | x"]), "OK\n");
    assert_eq!(ask(&addr, &[&format!("t PUT {out} 101 slow 50 0"), "101 w | x"]), "OK\n");
    let fast = ask(&addr, &["t GET fast cpu 0 4"]);
    assert_eq!(fast, format!("JOB {out} 1 fake 1\n"), "the fast machine takes the weight-5 task");
    let slow = ask(&addr, &["t GET slow cpu 0 4"]);
    assert_eq!(slow, format!("JOB {out} 0 fake 0\n"), "the slow machine takes the weight-1 task");
    // A machine nobody has measured is in the middle: of 3 and 3, the first in order.
    let new = ask(&addr, &["t GET new cpu 0 4"]);
    assert_eq!(new, format!("JOB {out} 2 fake 2\n"));
    for (k, host) in [(1, "fast"), (0, "slow"), (2, "new")] {
        assert_eq!(ask(&addr, &[&format!("t PUT {out} {k} {host} 1 0"), &format!("{k} r | x")]), "OK\n");
    }
    assert!(ask(&addr, &["t GET new cpu 0 4"]).starts_with("JOB "));
    ask(&addr, &[&format!("t PUT {out} 3 new 1 0"), "3 r | x"]);
    server.join().unwrap();
}

/// A lease is as long as the machine's own speed says the task needs: a
/// machine measured at 100 s per unit keeps a task past a lease that would have
/// run out for a faster one.
#[test]
fn a_slow_machine_is_given_a_longer_lease() {
    let d = scratch("slow-lease");
    let out = d.join("o.txt").to_string_lossy().to_string();
    let t = vec![sized(&out, "warm", 100, 1.0, None, "warm"), sized(&out, "cpu", 0, 1.0, None, "work"), sized(&out, "cpu", 1, 1.0, None, "work")];
    let (addr, server) = start(t, 300);
    assert!(ask(&addr, &["t GET slow warm 0 4"]).starts_with("JOB "));
    assert_eq!(ask(&addr, &[&format!("t PUT {out} 100 slow 100 0"), "100 w | x"]), "OK\n");
    assert!(ask(&addr, &["t GET slow cpu 0 4"]).starts_with("JOB "));
    assert!(ask(&addr, &["t GET quick cpu 0 4"]).starts_with("JOB "), "an unmeasured machine gets the other task on the plain lease");
    std::thread::sleep(Duration::from_millis(500));
    // The plain lease has run out and the slow machine's has not.
    let reply = ask(&addr, &["t GET third cpu 0 4"]);
    assert_eq!(reply, format!("JOB {out} 1 fake 1\n"), "{reply}");
    assert_eq!(ask(&addr, &["t GET fourth cpu 0 4"]), "WAIT 1\n");
    ask(&addr, &[&format!("t PUT {out} 0 slow 1 0"), "0 r | x"]);
    ask(&addr, &[&format!("t PUT {out} 1 third 1 0"), "1 r | x"]);
    server.join().unwrap();
}

/// A worker that offers too little memory for anything left stops and says so
/// instead of asking again.
#[test]
fn a_worker_too_small_for_what_is_left_stops() {
    let d = scratch("toosmall");
    let out = d.join("o.txt").to_string_lossy().to_string();
    let (addr, server) = start(vec![sized(&out, "cpu", 0, 1.0, Some(5000), "big")], 60_000);
    let mut cfg = config(&addr, "pi", "cpu");
    cfg.mem_mb = Some(2000);
    assert_eq!(work(&cfg, |_| panic!("must not be asked")), Ended::TooSmall);
    cfg.mem_mb = Some(8000);
    assert_eq!(work(&cfg, |_| Ok("0 r | x".into())), Ended::NothingLeft);
    server.join().unwrap();
}

/// `weight=` and `mem=` in a plan size its tasks.
#[test]
fn a_plan_can_state_weight_and_memory() {
    let d = scratch("plan-sized");
    let out = d.join("p.txt").to_string_lossy().to_string();
    let t = phys::pairs::plan(&format!("pair water 0 2 element cpu {out} weight=20 mem=6.5\n")).unwrap();
    assert_eq!(t.len(), 2);
    assert_eq!((t[0].weight, t[0].need_mb), (20.0, Some(6656)));
    assert_eq!(t[0].kind, "water element");
    let plain = phys::pairs::plan(&format!("pair water 0 1 element cpu {out}\n")).unwrap();
    assert_eq!((plain[0].weight, plain[0].need_mb), (1.0, None));
    assert!(phys::pairs::plan(&format!("pair water 0 2 element cpu {out} weight=-1\n")).unwrap_err().contains("weight"));
    assert!(phys::pairs::plan(&format!("pair water 0 2 element cpu {out} mem=lots\n")).unwrap_err().contains("mem"));
}
