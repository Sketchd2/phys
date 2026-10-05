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
//! **Classifying machines.** A worker says what it is on every request (its
//! class, its cores, the memory it has free) and the server learns what it
//! can do by what it does: seconds per unit of task weight, per machine, and
//! the most memory any task of a kind has been seen to use (a worker reports
//! its peak with each result). From those a task goes only to a machine it
//! fits in; the heaviest tasks go to the fastest machines and the lightest to
//! the slowest, so a slow machine is not holding the last big task while the
//! rest sit idle; and a lease is as long as the machine's own speed says the
//! task needs. None of it is declared by the machine, except what it says it
//! has free, which is measured by the machine's own operating system. What a
//! plan may say is a task's `weight` and, before anything has been measured,
//! its `mem`.
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
//! <token> GET <host> <class> <mem_mb> <cores>  ->  JOB <out> <index> <spec...> | WAIT <secs> | NOFIT <mb> | DONE
//! <token> PUT <out> <index> <host> <secs> <peak_mb>
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
#[derive(Clone, Debug, Default)]
pub struct Task {
    pub out: String,
    pub index: usize,
    pub class: String,
    /// What a worker needs to do it, in words; see `pairs::Executor`.
    pub spec: String,
    /// Written first if `out` is new.
    pub header: String,
    /// What sort of task this is, for what is learned about it: tasks of one
    /// kind are taken to need the same memory.
    pub kind: String,
    /// How much work, relative to a pair (1). Sizes the lease and decides
    /// which machine gets the heavy ones; 0 is read as 1.
    pub weight: f64,
    /// Memory it needs, in MB, if the plan says. Otherwise the most any task
    /// of its kind has been seen to use, and until one has finished, nothing.
    pub need_mb: Option<u64>,
}

impl Task {
    fn weight(&self) -> f64 {
        if self.weight > 0.0 { self.weight } else { 1.0 }
    }
}

enum Status {
    Pending,
    Leased { at: Instant, host: String, lease: Duration },
    Done,
}

/// What the server knows of one machine: what it said it is, and what it has
/// shown by doing. The speed is measured (seconds per unit of task weight) and
/// never declared: a machine's own account of how fast it is does not decide
/// who gets the heavy tasks.
#[derive(Default)]
struct Host {
    class: String,
    cores: usize,
    mem_mb: u64,
    done: usize,
    seconds: f64,
    weight: f64,
}

impl Host {
    /// Seconds per unit of weight so far, if it has finished anything.
    fn rate(&self) -> Option<f64> {
        if self.done > 0 && self.weight > 0.0 { Some(self.seconds / self.weight) } else { None }
    }
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
    hosts: BTreeMap<String, Host>,
    /// The most memory (MB) any finished task of a kind has used.
    peaks: HashMap<String, u64>,
    headered: HashSet<String>,
    log: Option<String>,
    finished: Option<Instant>,
}

impl Shared {
    /// How long a worker with nothing to take is told to wait before asking
    /// again: a fiftieth of the lease, between a second and a minute. It must
    /// be well inside the grace [`serve`] gives at the end, or a worker told
    /// to wait outlives the server and gives up on it.
    fn wait_secs(&self) -> u64 {
        (self.lease.as_secs() / 50).clamp(1, 60)
    }

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
        for (name, h) in &self.hosts {
            let rate = h.rate().map(|r| format!("{r:.0} s per unit of weight")).unwrap_or_else(|| "not yet measured".into());
            s += &format!("\n  {name} ({}): {} cores, {} MB free when it last asked, {} done, {rate}", h.class, h.cores, h.mem_mb, h.done);
        }
        for (kind, mb) in &self.peaks {
            s += &format!("\n  {kind}: most memory seen {mb} MB");
        }
        for (st, t) in self.status.iter().zip(&self.tasks) {
            if let Status::Leased { at, host, .. } = st {
                s += &format!("\n  {} #{} held by {host} for {:.0} min", t.out, t.index, at.elapsed().as_secs_f64() / 60.0);
            }
        }
        s
    }

    /// How much memory (MB) a task needs: what the plan said, else the most a
    /// finished task of its kind used, else nothing known.
    fn need(&self, i: usize) -> Option<u64> {
        self.tasks[i].need_mb.or_else(|| self.peaks.get(&self.tasks[i].kind).copied())
    }

    /// Where `host` stands among the machines that have finished anything, as
    /// 1.0 for the fastest and 0.0 for the slowest; 0.5 if it has not been
    /// measured, 1.0 if it is the only one.
    fn speed_of(&self, host: &str) -> f64 {
        let mine = match self.hosts.get(host).and_then(|h| h.rate()) {
            Some(r) => r,
            None => return 0.5,
        };
        let rates: Vec<f64> = self.hosts.values().filter_map(|h| h.rate()).collect();
        if rates.len() < 2 {
            return 1.0;
        }
        let faster = rates.iter().filter(|&&r| r < mine).count();
        1.0 - faster as f64 / (rates.len() - 1) as f64
    }

    /// The task a worker should do next, and mark it leased. It is of the
    /// worker's class and fits in the memory it has free. Among those, the
    /// heaviest go to the fastest machines and the lightest to the slowest
    /// (equal weights go in order); a task whose lease has run out is taken
    /// only when nothing is waiting for a first try.
    fn lease_for(&mut self, host: &str, class: &str, mem_mb: u64) -> Reply {
        let now = Instant::now();
        let fits = |need: Option<u64>| mem_mb == 0 || need.map(|n| n <= mem_mb).unwrap_or(true);
        let mut pending = Vec::new();
        let mut expired = Vec::new();
        let mut waiting = false;
        let mut too_big: Option<u64> = None;
        for i in 0..self.tasks.len() {
            if self.tasks[i].class != class || self.failures[i] >= MAX_FAILURES {
                continue;
            }
            let need = self.need(i);
            let fit = fits(need);
            match &self.status[i] {
                Status::Done => {}
                Status::Pending if fit => pending.push(i),
                Status::Leased { at, lease, .. } if fit => {
                    if now.duration_since(*at) >= *lease {
                        expired.push(i);
                    } else {
                        waiting = true;
                    }
                }
                _ => too_big = Some(too_big.unwrap_or(0).max(need.unwrap_or(0))),
            }
        }
        let chosen = if !pending.is_empty() {
            let mut by_weight = pending.clone();
            by_weight.sort_by(|&a, &b| self.tasks[b].weight().partial_cmp(&self.tasks[a].weight()).unwrap());
            let at = ((1.0 - self.speed_of(host)) * (by_weight.len() - 1) as f64).round() as usize;
            let target = self.tasks[by_weight[at]].weight();
            pending.into_iter().find(|&i| self.tasks[i].weight() == target)
        } else {
            expired.first().copied()
        };
        match chosen {
            Some(i) => {
                // A lease long enough for this machine to finish this task, at
                // the speed it has shown, with room to spare.
                let mut lease = self.lease;
                if let Some(r) = self.hosts.get(host).and_then(|h| h.rate()) {
                    lease = lease.max(Duration::from_secs_f64(4.0 * r * self.tasks[i].weight()));
                }
                self.status[i] = Status::Leased { at: now, host: host.to_string(), lease };
                Reply::Job(i)
            }
            None if waiting => Reply::Wait,
            None => match too_big {
                Some(n) => Reply::NoFit(n),
                None => Reply::Done,
            },
        }
    }
}

enum Reply {
    Job(usize),
    Wait,
    NoFit(u64),
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
    // A file that cannot be written is found now, not when the first result
    // has cost an hour.
    for out in tasks.iter().map(|t| t.out.as_str()).collect::<HashSet<_>>() {
        std::fs::OpenOptions::new().create(true).append(true).open(out).map_err(|e| std::io::Error::new(e.kind(), format!("cannot write {out}: {e}")))?;
    }
    let mut index = HashMap::new();
    for (i, t) in tasks.iter().enumerate() {
        index.insert((t.out.clone(), t.index), i);
    }
    let status = tasks.iter().map(|_| Status::Pending).collect();
    let failures = vec![0; tasks.len()];
    let finished = if tasks.is_empty() { Some(Instant::now()) } else { None };
    let shared = Arc::new(Mutex::new(Shared { tasks, status, failures, index, lease, token: token.to_string(), started: Instant::now(), completed: 0, hosts: BTreeMap::new(), peaks: HashMap::new(), headered: HashSet::new(), log, finished }));
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
                Some("GET") if w.len() >= 6 => {
                    let mem: u64 = w[4].parse().unwrap_or(0);
                    let h = s.hosts.entry(w[2].to_string()).or_default();
                    h.class = w[3].to_string();
                    h.mem_mb = mem;
                    h.cores = w[5].parse().unwrap_or(0);
                    match s.lease_for(w[2], w[3], mem) {
                        Reply::Job(i) => {
                            let t = &s.tasks[i];
                            format!("JOB {} {} {}\n", t.out, t.index, t.spec)
                        }
                        Reply::Wait => format!("WAIT {}\n", s.wait_secs()),
                        Reply::NoFit(n) => format!("NOFIT {n}\n"),
                        Reply::Done => "DONE\n".to_string(),
                    }
                }
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
                Some("PUT") if w.len() >= 7 => {
                    let mut result = String::new();
                    reader.read_line(&mut result)?;
                    let result = result.trim_end().to_string();
                    put(&mut s, w[2], w[3], w[4], w[5], w[6], &result)
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
fn put(s: &mut Shared, file: &str, index: &str, host: &str, secs: &str, peak: &str, result: &str) -> String {
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
    let weight = s.tasks[i].weight();
    let h = s.hosts.entry(host.to_string()).or_default();
    h.done += 1;
    h.seconds += seconds;
    h.weight += weight;
    let peak: u64 = peak.parse().unwrap_or(0);
    if peak > 0 {
        let kind = s.tasks[i].kind.clone();
        let e = s.peaks.entry(kind).or_insert(0);
        *e = (*e).max(peak);
    }
    if let Some(log) = &s.log {
        if let Ok(mut l) = std::fs::OpenOptions::new().create(true).append(true).open(log) {
            let _ = writeln!(l, "{} {host} {} {file} {k} {seconds:.1}", unix_now(), s.tasks[i].class);
        }
    }
    let first = s.summary().lines().next().unwrap_or("").to_string();
    println!("[{}] {file} #{k} from {host} in {seconds:.0} s: {first}", unix_now());
    "OK\n".into()
}

/// Memory as the operating system reports it: what is free to use now, and
/// the most this process has held. Zero means it cannot be measured here.
mod sys {
    #[cfg(target_os = "linux")]
    fn proc_kb(file: &str, key: &str) -> u64 {
        std::fs::read_to_string(file).ok().and_then(|t| t.lines().find(|l| l.starts_with(key)).and_then(|l| l.split_whitespace().nth(1)?.parse().ok())).unwrap_or(0)
    }

    #[cfg(target_os = "linux")]
    pub fn available_mb() -> u64 {
        proc_kb("/proc/meminfo", "MemAvailable:") / 1024
    }

    #[cfg(target_os = "linux")]
    pub fn peak_mb() -> u64 {
        proc_kb("/proc/self/status", "VmHWM:") / 1024
    }

    #[cfg(windows)]
    #[repr(C)]
    struct MemoryStatusEx {
        length: u32,
        load: u32,
        total_phys: u64,
        avail_phys: u64,
        total_page: u64,
        avail_page: u64,
        total_virtual: u64,
        avail_virtual: u64,
        avail_extended: u64,
    }

    #[cfg(windows)]
    #[repr(C)]
    struct ProcessMemoryCounters {
        cb: u32,
        page_faults: u32,
        peak_working_set: usize,
        working_set: usize,
        quota_peak_paged: usize,
        quota_paged: usize,
        quota_peak_nonpaged: usize,
        quota_nonpaged: usize,
        pagefile: usize,
        peak_pagefile: usize,
    }

    #[cfg(windows)]
    extern "system" {
        fn GlobalMemoryStatusEx(buffer: *mut MemoryStatusEx) -> i32;
        fn GetCurrentProcess() -> isize;
        fn K32GetProcessMemoryInfo(process: isize, counters: *mut ProcessMemoryCounters, cb: u32) -> i32;
    }

    #[cfg(windows)]
    pub fn available_mb() -> u64 {
        let mut m = MemoryStatusEx { length: std::mem::size_of::<MemoryStatusEx>() as u32, load: 0, total_phys: 0, avail_phys: 0, total_page: 0, avail_page: 0, total_virtual: 0, avail_virtual: 0, avail_extended: 0 };
        // SAFETY: `m` is a correctly sized, initialised struct whose length field is set, as the call requires.
        if unsafe { GlobalMemoryStatusEx(&mut m) } != 0 { m.avail_phys / (1024 * 1024) } else { 0 }
    }

    #[cfg(windows)]
    pub fn peak_mb() -> u64 {
        let mut c = ProcessMemoryCounters { cb: std::mem::size_of::<ProcessMemoryCounters>() as u32, page_faults: 0, peak_working_set: 0, working_set: 0, quota_peak_paged: 0, quota_paged: 0, quota_peak_nonpaged: 0, quota_nonpaged: 0, pagefile: 0, peak_pagefile: 0 };
        // SAFETY: the pseudo-handle for this process is always valid, and `c` is a correctly sized struct with `cb` set.
        if unsafe { K32GetProcessMemoryInfo(GetCurrentProcess(), &mut c, c.cb) } != 0 { (c.peak_working_set / (1024 * 1024)) as u64 } else { 0 }
    }

    #[cfg(not(any(target_os = "linux", windows)))]
    pub fn available_mb() -> u64 {
        0
    }

    #[cfg(not(any(target_os = "linux", windows)))]
    pub fn peak_mb() -> u64 {
        0
    }
}

/// What a worker is told to do and how long to keep trying.
pub struct WorkerConfig {
    /// `host:port` of the server.
    pub server: String,
    pub token: String,
    /// What this machine calls itself in the server's log.
    pub name: String,
    pub class: String,
    /// The memory this machine offers, in MB; `None` measures what is free at
    /// each request (0, where that cannot be measured, means no limit).
    pub mem_mb: Option<u64>,
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
    /// Everything that is left needs more memory than this machine has.
    TooSmall,
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
        let mem = cfg.mem_mb.unwrap_or_else(sys::available_mb);
        let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
        let reply = match exchange_patiently(cfg, &[format!("{} GET {} {} {mem} {cores}", cfg.token, cfg.name, cfg.class)]) {
            Some(r) => r,
            None => return Ended::ServerGone,
        };
        let w: Vec<&str> = reply.split_whitespace().collect();
        match w.first().copied() {
            Some("DONE") => return Ended::NothingLeft,
            Some("NOFIT") => {
                eprintln!("what is left needs {} MB and this machine offers {mem}", w.get(1).unwrap_or(&"more"));
                return Ended::TooSmall;
            }
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
                        let put = [format!("{} PUT {file} {index} {} {seconds:.1} {}", cfg.token, cfg.name, sys::peak_mb()), line.clone()];
                        let keep = |why: &str| {
                            // Keep what cost an hour: it can be added by hand.
                            let keep = format!("unsent-{}.txt", cfg.name);
                            if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&keep) {
                                let _ = writeln!(f, "{file} {line}");
                            }
                            eprintln!("the result was not accepted ({why}); kept in {keep}");
                        };
                        match exchange_patiently(cfg, &put) {
                            Some(r) if r.starts_with("OK") || r.starts_with("DUP") => println!("{file} #{index}: {seconds:.0} s, server says {}", r.trim()),
                            Some(r) => keep(r.trim()),
                            None => {
                                keep("the server could not be reached");
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
