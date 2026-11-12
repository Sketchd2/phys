//! `phys-overlap`: the overlap of two molecules' electron densities for pairs
//! already computed, beside their exact electrostatic energy. The work is
//! `phys::pairs::overlap_main`; see there.
//!
//! ```sh
//! cargo run --release --bin phys-overlap -- water pairs-water-mp2.txt [--every n] [--max n] [--out file]
//! ```

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::overlap_main(&args);
}
