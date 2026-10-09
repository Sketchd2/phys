//! `phys-es-gpu`: the exact electrostatic interaction of pairs already computed,
//! with Hartree-Fock on the GPU. See `phys::pairs::es_main`.

fn main() {
    if std::env::var_os("PHYS_NO_GPU_FOCK").is_none() {
        match phys_gpu::install_fock() {
            Ok(name) => println!("Hartree-Fock passes on {name}"),
            Err(e) => eprintln!("no GPU for Hartree-Fock: {e}"),
        }
    }
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--deriv") {
        phys::pairs::deriv_main(&args);
    } else {
        phys::pairs::es_main(&args);
    }
}
