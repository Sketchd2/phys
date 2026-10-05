//! A worker for the work queue, on the CPU — see `phys::queue`.
//!
//! ```sh
//! phys-worker host:port [--name pi-017] [--jobs 10] [--patience 30] [--mem GB]
//! ```
//!
//! Needs, in its working directory, the same `grow-<name>.state` (and any
//! snapshot the plan names) as the server. Exits 10 after `--jobs` tasks so a
//! wrapper can start a fresh process (a long-lived one slows down), 0 when
//! the server says nothing is left, 3 when it cannot be reached, 4 when
//! everything left needs more memory than it has free (`--mem GB` offers less).

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    phys::pairs::work_main(&args, "cpu");
}
