//! `phys::pairs::run_unlike` with Hartree-Fock and MP2 on the GPU in double precision
//! (`PHYS_NO_GPU_FOCK=1` leaves them on the CPU): MP2 interaction energies of pairs of two
//! different molecules (`docs/SOLUTIONS.md` S2).
//!
//! ```sh
//! phys-pairs-unlike-gpu water ammonia 40 [--wall 1.0 2.3]
//! ```

fn main() {
    if std::env::var_os("PHYS_NO_GPU_FOCK").is_none() {
        match phys_gpu::install_fock() {
            Ok(name) => println!("Hartree-Fock passes on {name}"),
            Err(e) => eprintln!("no GPU for Hartree-Fock: {e}"),
        }
    }
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::run_unlike(&args);
}
