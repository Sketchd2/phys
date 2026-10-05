//! A work queue for many machines — `PLAY.md` E8's data is hundreds of
//! independent pair calculations, and the machines that can do them differ in
//! speed by two orders of magnitude and may vanish mid-job.
//!
//! One machine runs [`serve`] (`phys-queue`), holding a list of [`Task`]s; any
//! number run [`work`] (`phys-worker`, `phys-worker-gpu`), each asking for a
//! task, doing it, and handing the result line back. Nothing is shared but
//! that conversation, so a machine needs only the same commit, the same
//! `grow-<name>.state` and snapshot files, and a route to the server.
//!
//! **Leases, not assignment.** A task handed out is leased for a time; if the
//! result has not come back by then (a Pi's battery died, a laptop was shut)
//! the task is handed to the next machine that asks, and the first result to
//! arrive is the one kept. A result for a task already done is dropped
//! (`DUP`). That is the whole fault tolerance: the server never needs to know
//! a machine has gone.
//!
//! **Class.** A worker names a class, and only tasks of that class go to it.
//! The GPU's final non-local energy is single precision and the CPU's double;
//! which of the two a dataset holds is a decision about the data, so a plan
//! says it per task and the queue cannot mix them.
//!
//! **Resuming.** The server's state is the output files themselves: a task is
//! done when its index is in its file, and [`crate::pairs::plan`] leaves out
//! the ones that are. A restarted server loses only the leases in flight.
//!
//! The protocol is one line of text, a reply of text, over TCP, with no
//! dependencies. There is **no authentication beyond a shared token** and no
//! encryption: use it on a network you trust, and do not point it at one you
//! do not.
//!
//! ```text
//! <token> GET <host> <class>             ->  JOB <out> <index> <spec...> | WAIT <secs> | DONE
//! <token> PUT <out> <index> <host> <secs>
//! <result line>                          ->  OK | DUP | ERR <why>
//! <token> FAIL <out> <index> <host> <why> ->  OK
//! <token> STATUS                         ->  a few lines of text
//! ```

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream, ToSocketAddrs};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// One unit of work: the line it will add to `out`, as index `index`.
#[derive(Clone, Debug)]
pub struct Task {
    pub out: String,
    pub index: usize,
    pub class: String,
    /// What a worker needs to do it, in words; see `pairs::Executor`.
    pub spec: String,
    /// Written first if `out` is new.
    pub header: String,
}

enum Status {
    Pending,
    Leased { at: Instant, host: String },
    Done,
}

/// A task that fails this many times is left alone and reported: a failure
/// that repeats on different machines is a fact about the task.
const MAX_FAILURES: usize = 3;

struct Shared {
    tasks: Vec<Task>,
    status: Vec<Status>,
    failures: Vec<usize>,
    index: HashMap<(String, usize), usize>,
    lease: Duration,
    token: String,
    started: Instant,
    completed: usize,
    per_host: BTreeMap<String, (usize, f64)>,
    headered: HashSet<String>,
    log: Option<String>,
    finished: Option<Instant>,
}

impl Shared {
    fn remaining(&self) -> usize {
        self.status.iter().filter(|s| !matches!(s, Status::Done)).count()
    }

    fn leased(&self) -> usize {
        self.status.iter().filter(|s| matches!(s, Status::Leased { .. })).count()
    }

    fn summary(&self) -> String {
        let total = self.tasks.len();
        let left = self.remaining();
        let done = total - left;
        let elapsed = self.started.elapsed().as_secs_f64();
        let eta = if self.completed > 0 && left > 0 { format!(", about {:.0} min to go at this rate", elapsed / self.completed as f64 * left as f64 / 60.0) } else { String::new() };
        let mut s = format!("{done}/{total} done, {} leased, {} failed for good{eta}", self.leased(), self.failures.iter().filter(|&&f| f >= MAX_FAILURES).count());
        for (host, (n, secs)) in &self.per_host {
            s += &format!("\n  {host}: {n} done, {:.0} s each", secs / *n as f64);
        }
        for (st, t) in self.status.iter().zip(&self.tasks) {
            if let Status::Leased { at, host } = st {
                s += &format!("\n  {} #{} held by {host} for {:.0} min", t.out, t.index, at.elapsed().as_secs_f64() / 60.0);
            }
        }
        s
    }

    /// The task a worker of `class` should do next, and mark it leased.
    fn lease_for(&mut self, host: &str, class: &str) -> Reply {
        let now = Instant::now();
        let mut waiting = false;
        let mut stolen = None;
        for i in 0..self.tasks.len() {
            if self.tasks[i].class != class || self.failures[i] >= MAX_FAILURES {
                continue;
            }
            match &self.status[i] {
                Status::Done => {}
                Status::Pending => {
                    self.status[i] = Status::Leased { at: now, host: host.to_string() };
                    return Reply::Job(i);
                }
                Status::Leased { at, .. } => {
                    if now.duration_since(*at) >= self.lease {
                        stolen.get_or_insert(i);
                    } else {
                        waiting = true;
                    }
                }
            }
        }
        if let Some(i) = stolen {
            self.status[i] = Status::Leased { at: now, host: host.to_string() };
            return Reply::Job(i);
        }
        if waiting {
            Reply::Wait
        } else {
            Reply::Done
        }
    }
}

enum Reply {
    Job(usize),
    Wait,
    Done,
}

fn unix_now() -> u64 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0)
}

/// Serve `tasks` on `listener` until every one is done (and `grace` longer, so
/// the last workers are told so), appending each result to its task's file.
/// `log`, if given, gets one line a result: when, who, class, file, index,
/// seconds it took. Returns how many were completed in this run.
pub fn serve(listener: TcpListener, tasks: Vec<Task>, lease: Duration, token: &str, log: Option<String>, grace: Duration) -> std::io::Result<usize> {
    let mut index = HashMap::new();
    for (i, t) in tasks.iter().enumerate() {
        index.insert((t.out.clone(), t.index), i);
    }
    let status = tasks.iter().map(|_| Status::Pending).collect();
    let failures = vec![0; tasks.len()];
    let finished = if tasks.is_empty() { Some(Instant::now()) } else { None };
    let shared = Arc::new(Mutex::new(Shared { tasks, status, failures, index, lease, token: token.to_string(), started: Instant::now(), completed: 0, per_host: BTreeMap::new(), headered: HashSet::new(), log, finished }));
    listener.set_nonblocking(true)?;
    loop {
        match listener.accept() {
            Ok((stream, _)) => {
                let shared = Arc::clone(&shared);
                std::thread::spawn(move || {
                    let _ = handle(&shared, stream);
                });
            }
            Err(e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                std::thread::sleep(Duration::from_millis(50));
            }
            Err(e) => return Err(e),
        }
        let mut s = shared.lock().unwrap();
        if s.finished.is_none() && s.remaining() == 0 {
            s.finished = Some(Instant::now());
            println!("[{}] every task is done", unix_now());
        }
        // A task that has failed for good counts as finished for the purpose of
        // stopping: nothing will ever hand it out again.
        let stuck = s.status.iter().zip(&s.failures).all(|(st, &f)| matches!(st, Status::Done) || f >= MAX_FAILURES);
        if s.finished.is_none() && stuck {
            s.finished = Some(Instant::now());
            println!("[{}] nothing left that can be handed out", unix_now());
        }
        if let Some(t) = s.finished {
            if t.elapsed() > grace {
                return Ok(s.completed);
            }
        }
    }
}

fn handle(shared: &Mutex<Shared>, stream: TcpStream) -> std::io::Result<()> {
    stream.set_read_timeout(Some(Duration::from_secs(20)))?;
    stream.set_write_timeout(Some(Duration::from_secs(20)))?;
    let mut reader = BufReader::new(stream.try_clone()?);
    let mut out = stream;
    let mut line = String::new();
    reader.read_line(&mut line)?;
    let w: Vec<&str> = line.split_whitespace().collect();
    let reply = {
        let mut s = shared.lock().unwrap();
        if w.first().copied() != Some(s.token.as_str()) {
            "ERR token\n".to_string()
        } else {
            match w.get(1).copied() {
                Some("GET") if w.len() >= 4 => match s.lease_for(w[2], w[3]) {
                    Reply::Job(i) => {
                        let t = &s.tasks[i];
                        format!("JOB {} {} {}\n", t.out, t.index, t.spec)
                    }
                    Reply::Wait => "WAIT 60\n".to_string(),
                    Reply::Done => "DONE\n".to_string(),
                },
                Some("STATUS") => format!("{}\n", s.summary()),
                Some("FAIL") if w.len() >= 5 => {
                    if let Some(&i) = w[3].parse::<usize>().ok().and_then(|k| s.index.get(&(w[2].to_string(), k))) {
                        if !matches!(s.status[i], Status::Done) {
                            s.failures[i] += 1;
                            s.status[i] = Status::Pending;
                            println!("[{}] {} failed {} {}: {}", unix_now(), w[4], w[2], w[3], w[5..].join(" "));
                        }
                    }
                    "OK\n".to_string()
                }
                Some("PUT") if w.len() >= 6 => {
                    let mut result = String::new();
                    reader.read_line(&mut result)?;
                    let result = result.trim_end().to_string();
                    put(&mut s, w[2], w[3], w[4], w[5], &result)
                }
                _ => "ERR unknown request\n".to_string(),
            }
        }
    };
    out.write_all(reply.as_bytes())?;
    Ok(())
}

/// Record one result: it must be the line for the task it names, and the first
/// for that task.
fn put(s: &mut Shared, file: &str, index: &str, host: &str, secs: &str, result: &str) -> String {
    let k: usize = match index.parse() {
        Ok(k) => k,
        Err(_) => return "ERR index\n".into(),
    };
    let i = match s.index.get(&(file.to_string(), k)) {
        Some(&i) => i,
        None => return "DUP\n".into(),
    };
    if matches!(s.status[i], Status::Done) {
        return "DUP\n".into();
    }
    if result.split_whitespace().next().and_then(|t| t.parse::<usize>().ok()) != Some(k) || result.starts_with('#') {
        return "ERR the line is not this task's\n".into();
    }
    let path = std::path::Path::new(file);
    let fresh = !s.headered.contains(file) && std::fs::metadata(path).map(|m| m.len() == 0).unwrap_or(true);
    let mut f = match std::fs::OpenOptions::new().create(true).append(true).open(path) {
        Ok(f) => f,
        Err(e) => return format!("ERR {e}\n"),
    };
    if fresh {
        let _ = writeln!(f, "{}", s.tasks[i].header);
    }
    s.headered.insert(file.to_string());
    if writeln!(f, "{result}").is_err() {
        return "ERR write\n".into();
    }
    let _ = f.flush();
    s.status[i] = Status::Done;
    s.completed += 1;
    let seconds: f64 = secs.parse().unwrap_or(0.0);
    let e = s.per_host.entry(host.to_string()).or_insert((0, 0.0));
    e.0 += 1;
    e.1 += seconds;
    if let Some(log) = &s.log {
        if let Ok(mut l) = std::fs::OpenOptions::new().create(true).append(true).open(log) {
            let _ = writeln!(l, "{} {host} {} {file} {k} {seconds:.1}", unix_now(), s.tasks[i].class);
        }
    }
    let first = s.summary().lines().next().unwrap_or("").to_string();
    println!("[{}] {file} #{k} from {host} in {seconds:.0} s: {first}", unix_now());
    "OK\n".into()
}

/// What a worker is told to do and how long to keep trying.
pub struct WorkerConfig {
    /// `host:port` of the server.
    pub server: String,
    pub token: String,
    /// What this machine calls itself in the server's log.
    pub name: String,
    pub class: String,
    /// Stop after this many tasks (a long-lived process slows down; a wrapper
    /// starts a fresh one).
    pub max_jobs: Option<usize>,
    /// How long to keep retrying a server that cannot be reached.
    pub patience: Duration,
}

/// Why a worker stopped.
#[derive(Debug, PartialEq)]
pub enum Ended {
    /// The server said nothing is left for this class.
    NothingLeft,
    /// `max_jobs` were done.
    Limit,
    /// The server could not be reached for `patience`.
    ServerGone,
}

/// Send `lines` to the server and read its whole reply.
fn exchange(server: &str, lines: &[String]) -> std::io::Result<String> {
    let addr = server.to_socket_addrs()?.next().ok_or_else(|| std::io::Error::new(std::io::ErrorKind::NotFound, "no address"))?;
    let mut stream = TcpStream::connect_timeout(&addr, Duration::from_secs(10))?;
    stream.set_read_timeout(Some(Duration::from_secs(30)))?;
    stream.set_write_timeout(Some(Duration::from_secs(30)))?;
    for l in lines {
        stream.write_all(l.as_bytes())?;
        stream.write_all(b"\n")?;
    }
    let mut reply = String::new();
    stream.read_to_string(&mut reply)?;
    Ok(reply)
}

/// [`exchange`], retried every 15 seconds for `patience`.
fn exchange_patiently(cfg: &WorkerConfig, lines: &[String]) -> Option<String> {
    let start = Instant::now();
    loop {
        match exchange(&cfg.server, lines) {
            Ok(r) => return Some(r),
            Err(e) => {
                if start.elapsed() >= cfg.patience {
                    eprintln!("the server {} cannot be reached: {e}", cfg.server);
                    return None;
                }
                eprintln!("the server {} cannot be reached ({e}); trying again", cfg.server);
                std::thread::sleep(Duration::from_secs(15));
            }
        }
    }
}

/// The server's status text, for `phys-queue --status`.
pub fn status(server: &str, token: &str) -> std::io::Result<String> {
    exchange(server, &[format!("{token} STATUS")])
}

/// Ask for tasks and do them with `run` (which turns a spec into the line to
/// record, or says why it cannot) until told there are none, `max_jobs` are
/// done, or the server is gone. A task that panics is reported failed, not
/// fatal: a field that does not converge is a fact about that pair.
pub fn work(cfg: &WorkerConfig, mut run: impl FnMut(&str) -> Result<String, String>) -> Ended {
    let mut jobs = 0;
    loop {
        if cfg.max_jobs.map(|m| jobs >= m).unwrap_or(false) {
            return Ended::Limit;
        }
        let reply = match exchange_patiently(cfg, &[format!("{} GET {} {}", cfg.token, cfg.name, cfg.class)]) {
            Some(r) => r,
            None => return Ended::ServerGone,
        };
        let w: Vec<&str> = reply.split_whitespace().collect();
        match w.first().copied() {
            Some("DONE") => return Ended::NothingLeft,
            Some("WAIT") => {
                let secs: u64 = w.get(1).and_then(|s| s.parse().ok()).unwrap_or(60);
                std::thread::sleep(Duration::from_secs(secs));
            }
            Some("JOB") if w.len() >= 4 => {
                let (file, index, spec) = (w[1], w[2], w[3..].join(" "));
                let t = Instant::now();
                let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| run(&spec))).unwrap_or_else(|_| Err("the calculation panicked".into()));
                let seconds = t.elapsed().as_secs_f64();
                match outcome {
                    Ok(line) => {
                        let put = [format!("{} PUT {file} {index} {} {seconds:.1}", cfg.token, cfg.name), line.clone()];
                        match exchange_patiently(cfg, &put) {
                            Some(r) => println!("{file} #{index}: {seconds:.0} s, server says {}", r.trim()),
                            None => {
                                // Keep what cost an hour: it can be added by hand.
                                let keep = format!("unsent-{}.txt", cfg.name);
                                if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&keep) {
                                    let _ = writeln!(f, "{file} {line}");
                                }
                                eprintln!("the result could not be delivered; kept in {keep}");
                                return Ended::ServerGone;
                            }
                        }
                    }
                    Err(why) => {
                        eprintln!("{file} #{index} failed: {why}");
                        let why = why.replace('\n', " ");
                        let _ = exchange_patiently(cfg, &[format!("{} FAIL {file} {index} {} {why}", cfg.token, cfg.name)]);
                    }
                }
                jobs += 1;
            }
            _ => {
                eprintln!("the server said {:?}; stopping", reply.trim());
                return Ended::ServerGone;
            }
        }
    }
}
