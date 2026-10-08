//! `phys-pairs` with the final fine-grid non-local energy on the GPU, in
//! single precision — `PLAY.md` E8. Writes `pairs-<name>-gpu.txt` beside the
//! CPU's `pairs-<name>.txt`, whose pairs are the same draws, so the two can be
//! compared pair by pair. The self-consistent fields stay on the CPU.
//!
//! ```sh
//! cargo run --release -p phys-gpu --bin phys-pairs-gpu -- water 10
//! ```

fn main() {
    match phys_gpu::install_fine() {
        Ok(name) => println!("final non-local energies on {name}"),
        Err(e) => {
            eprintln!("no GPU: {e}");
            std::process::exit(1);
        }
    }
    // For `--mp2`: Hartree-Fock and MP2 on the GPU in double precision, as
    // `phys-cluster-gpu` does (`PHYS_NO_GPU_FOCK=1` leaves them on the CPU).
    if std::env::var_os("PHYS_NO_GPU_FOCK").is_none() {
        match phys_gpu::install_fock() {
            Ok(name) => println!("Hartree-Fock passes on {name}"),
            Err(e) => eprintln!("no GPU for Hartree-Fock: {e}"),
        }
    }
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::run(&args, "-gpu");
}
