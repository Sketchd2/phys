//! Interaction energies of sampled pairs of one molecule — `PLAY.md` E8,
//! route A's data. The work is `phys::pairs::run`; see there.
//!
//! ```sh
//! cargo run --release --bin phys-pairs -- water 300 [element|grown]
//! ```

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::run(&args, "");
}
