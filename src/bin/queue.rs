//! The work queue's server — see `phys::queue` for what it does and why.
//!
//! ```sh
//! phys-queue plan.txt [--listen 0.0.0.0:7878] [--lease MINUTES] [--log queue.log]
//! phys-queue --status host:port
//! ```
//!
//! `plan.txt` is the plan `phys::pairs::plan` reads: one `pair` or `snap`
//! line per range of work, naming the class of machine that may do it and the
//! file its results go to. The lease is how long a machine may hold a task
//! before it is handed to another (default 240 minutes: a late result is not
//! wasted, only duplicated). The token every request carries is
//! `$PHYS_QUEUE_TOKEN` (default `open`); there is no other protection.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let token = std::env::var("PHYS_QUEUE_TOKEN").unwrap_or_else(|_| "open".into());
    let flag = |f: &str| args.iter().position(|a| a == f).and_then(|i| args.get(i + 1)).cloned();
    if let Some(server) = flag("--status") {
        match phys::queue::status(&server, &token) {
            Ok(s) => print!("{s}"),
            Err(e) => {
                eprintln!("cannot reach {server}: {e}");
                std::process::exit(1);
            }
        }
        return;
    }
    let plan_file = args.first().filter(|a| !a.starts_with("--")).cloned().unwrap_or_else(|| {
        eprintln!("usage: phys-queue plan.txt [--listen addr] [--lease minutes] [--log file] | --status host:port");
        std::process::exit(2);
    });
    let text = std::fs::read_to_string(&plan_file).unwrap_or_else(|e| panic!("cannot read {plan_file}: {e}"));
    let tasks = match phys::pairs::plan(&text) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };
    let listen = flag("--listen").unwrap_or_else(|| "0.0.0.0:7878".into());
    let lease = std::time::Duration::from_secs(60 * flag("--lease").and_then(|s| s.parse().ok()).unwrap_or(240));
    let listener = std::net::TcpListener::bind(&listen).unwrap_or_else(|e| panic!("cannot listen on {listen}: {e}"));
    let mut by_class: std::collections::BTreeMap<(String, String), usize> = std::collections::BTreeMap::new();
    for t in &tasks {
        *by_class.entry((t.class.clone(), t.out.clone())).or_default() += 1;
    }
    println!("{} tasks to do, listening on {listen}, lease {} min", tasks.len(), lease.as_secs() / 60);
    for ((class, out), n) in &by_class {
        println!("  {n} for {class} machines into {out}");
    }
    match phys::queue::serve(listener, tasks, lease, &token, flag("--log").or(Some("queue.log".into())), std::time::Duration::from_secs(90)) {
        Ok(n) => println!("finished: {n} results this run"),
        Err(e) => {
            eprintln!("the server stopped: {e}");
            std::process::exit(1);
        }
    }
}
