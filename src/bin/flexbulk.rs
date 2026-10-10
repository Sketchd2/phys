//! A liquid whose molecules move atom by atom: the work is `phys::flexbulk::run`; see there.

fn main() {
    phys::flexbulk::run(std::env::args().skip(1).collect());
}
