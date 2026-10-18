//! `phys-bulk` with the Ewald sum's reciprocal part on the GPU, in double precision
//! (`phys_gpu::ewald`); falls back to the CPU's where no card offers `f64`, or with `--cpu`.
//! The work is `phys::bulk::run`.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let off = args.iter().any(|a| a == "--cpu");
    let recip = args.iter().any(|a| a == "--recip");
    let args: Vec<String> = args.into_iter().filter(|a| a != "--cpu" && a != "--recip").collect();
    // The whole step on the card unless asked for the reciprocal sums alone (or for none).
    let step = if off || recip { None } else { phys_gpu::GpuEwaldStep::new().ok().map(std::sync::Arc::new) };
    match &step {
        Some(e) => {
            println!("Ewald step on the {}", e.name);
            phys::ewald::set_step_engine(Some(e.clone()));
        }
        None => {}
    }
    let engine = if off || step.is_some() { None } else { phys_gpu::GpuEwald::new().ok().map(std::sync::Arc::new) };
    match &engine {
        Some(e) => {
            println!("Ewald reciprocal sums on the {}", e.name);
            phys::ewald::set_reciprocal_engine(Some(e.clone()));
        }
        None => println!("Ewald reciprocal sums on the CPU"),
    }
    phys::bulk::run(args);
    if let Some(e) = step {
        let t = e.timing.lock().expect("the timing");
        if t.2 > 0 {
            println!("  GPU Ewald step: {} calls, {:.2} ms each on the card, {:.2} ms each with preparation and reading back", t.2, t.0 * 1000.0 / t.2 as f64, t.1 * 1000.0 / t.2 as f64);
        }
    }
    if let Some(e) = engine {
        let (secs, calls) = e.timing();
        if calls > 0 {
            println!("  GPU reciprocal sums: {calls} submissions, {:.2} ms each on the card (submit to readback)", secs * 1000.0 / calls as f64);
        }
    }
}
