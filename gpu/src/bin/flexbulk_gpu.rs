//! `phys-flexbulk` with the induced-dipole solve on the GPU, in double
//! precision (`phys_gpu::induction`); falls back to the CPU's where no card
//! offers `f64`. The work is `phys::flexbulk::run`.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let off = args.iter().any(|a| a == "--cpu");
    let args: Vec<String> = args.into_iter().filter(|a| a != "--cpu").collect();
    let engine = if off { None } else { phys_gpu::GpuInduction::new().ok().map(std::sync::Arc::new) };
    match &engine {
        Some(e) => {
            println!("induced dipoles on the {}", e.name);
            phys::induction::set_induction_engine(Some(e.clone()));
        }
        None => println!("induced dipoles on the CPU"),
    }
    let pairs = if off || args.iter().any(|a| a == "--cpu-pairs") { None } else { phys_gpu::GpuPairs::new().ok().map(std::sync::Arc::new) };
    match &pairs {
        Some(e) => {
            println!("pair loop on the {}", e.name);
            phys::liquid::set_pair_engine(Some(e.clone()));
        }
        None => println!("pair loop on the CPU"),
    }
    let args: Vec<String> = args.into_iter().filter(|a| a != "--cpu-pairs").collect();
    phys::flexbulk::run(args);
    if let Some(e) = engine {
        let t = e.timing.lock().expect("the timing");
        if t.2 > 0 {
            println!("  GPU solve: {} calls, {:.2} ms each on the card (submit to readback), {:.2} ms each with preparation and reading back", t.2, t.0 * 1000.0 / t.2 as f64, t.1 * 1000.0 / t.2 as f64);
        }
    }
    if let Some(e) = pairs {
        let t = e.timing.lock().expect("the timing");
        if t.2 > 0 {
            println!("  GPU pair loop: {} calls, {:.2} ms each on the card, {:.2} ms each with preparation and reading back", t.2, t.0 * 1000.0 / t.2 as f64, t.1 * 1000.0 / t.2 as f64);
        }
    }
}
