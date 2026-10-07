//! `phys-cluster` with the final fine-grid non-local energy on the GPU, in
//! single precision, as the pair data was made: the cluster and its pairs must
//! be computed the same way to be compared. See `phys::pairs::cluster_main`.
//!
//! ```sh
//! phys-cluster-gpu water bulk-water-298.snap --size 6 --law law.txt
//! ```

fn main() {
    match phys_gpu::install_fine() {
        Ok(name) => println!("final non-local energies on {name}"),
        Err(e) => {
            eprintln!("no GPU: {e}");
            std::process::exit(1);
        }
    }
    // For `--method mp2`: the exchange matrix and MP2's products on the GPU, in
    // single precision. `PHYS_NO_GPU_PRODUCTS=1` leaves them on the CPU in
    // double, which is how the precision of the GPU's is checked.
    if std::env::var_os("PHYS_NO_GPU_PRODUCTS").is_none() {
        match phys_gpu::install_product() {
            Ok(name) => println!("Hartree-Fock and MP2 products on {name}"),
            Err(e) => eprintln!("no GPU for products: {e}"),
        }
    }
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::cluster_main(&args);
}
