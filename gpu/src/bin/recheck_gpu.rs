//! `phys-recheck` with the final fine-grid non-local energy on the GPU, in
//! single precision: for a pairs file the GPU made. See `phys::pairs::recheck`.
//!
//! ```sh
//! phys-recheck-gpu pairs-water-gpu.txt water --sample 20 --seed 1
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
    std::process::exit(phys::pairs::recheck(&args));
}
