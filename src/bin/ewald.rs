//! Checks and measurements of the Ewald sum on a liquid snapshot; the work is `phys::ewaldcheck::run`.

fn main() {
    phys::ewaldcheck::run(std::env::args().skip(1).collect());
}
