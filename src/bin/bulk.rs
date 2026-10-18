//! `phys-bulk`: a liquid of rigid molecules in a periodic box; the work is `phys::bulk::run`.

fn main() {
    phys::bulk::run(std::env::args().skip(1).collect());
}
