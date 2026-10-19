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
    // Hartree-Fock and MP2 on the card in double precision for the `mp2` and `unlike` tasks (as
    // `phys-pairs-gpu --mp2` does); `PHYS_NO_GPU_FOCK=1` leaves them on the CPU.
    if std::env::var_os("PHYS_NO_GPU_FOCK").is_none() {
        match phys_gpu::install_fock() {
            Ok(name) => println!("Hartree-Fock passes on {name}"),
            Err(e) => eprintln!("no GPU for Hartree-Fock: {e}"),
        }
    }
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::work_main(&args, "gpu");
}
