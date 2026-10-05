//! A cluster of molecules cut from a simulated liquid, against the sum of its
//! pairs — see `phys::pairs::cluster_main`.
//!
//! ```sh
//! phys-cluster water bulk-water-298.snap --size 6 --centre 0 --law law.txt [--dry]
//! ```

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::cluster_main(&args);
}
