//! A worker for the work queue with the final fine-grid non-local energy on the
//! GPU, in single precision — class `gpu`; see `phys::queue`. The self-
//! consistent fields stay on the CPU.
//!
//! ```sh
//! phys-worker-gpu host:port [--name desk] [--jobs 10] [--patience 30] [--mem GB]
//! ```

fn main() {
    match phys_gpu::install_fine() {
        Ok(name) => println!("final non-local energies on {name}"),
        Err(e) => {
            eprintln!("no GPU: {e}");
            std::process::exit(1);
        }
    }
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::work_main(&args, "gpu");
}
