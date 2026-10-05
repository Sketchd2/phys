//! Pairs taken from a simulated liquid, with the final fine-grid non-local
//! energy on the GPU, in single precision — `PLAY.md` E8. Writes
//! `pairs-<name>-liq-gpu.txt`. See `phys::pairs::run_snapshot`.
//!
//! ```sh
//! cargo run --release -p phys-gpu --bin phys-pairs-liquid-gpu -- water bulk-water-298.snap 40
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
    phys::pairs::run_snapshot(&args, "-liq-gpu");
}
