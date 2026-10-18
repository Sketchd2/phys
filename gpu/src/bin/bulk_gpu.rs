//! `phys-bulk` with the Ewald sum's reciprocal part on the GPU, in double precision
//! (`phys_gpu::ewald`); falls back to the CPU's where no card offers `f64`, or with `--cpu`.
//! The work is `phys::bulk::run`.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let off = args.iter().any(|a| a == "--cpu");
    let args: Vec<String> = args.into_iter().filter(|a| a != "--cpu").collect();
    let engine = if off { None } else { phys_gpu::GpuEwald::new().ok().map(std::sync::Arc::new) };
    match &engine {
        Some(e) => {
            println!("Ewald reciprocal sums on the {}", e.name);
            phys::ewald::set_reciprocal_engine(Some(e.clone()));
        }
        None => println!("Ewald reciprocal sums on the CPU"),
    }
    phys::bulk::run(args);
    if let Some(e) = engine {
        let (secs, calls) = e.timing();
        if calls > 0 {
            println!("  GPU reciprocal sums: {calls} submissions, {:.2} ms each on the card (submit to readback)", secs * 1000.0 / calls as f64);
        }
    }
}
