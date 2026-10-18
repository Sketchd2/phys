//! `phys-ewald` with the Ewald step on the GPU: the whole of it (`phys_gpu::ewald_step`), or with
//! `--recip` only the reciprocal sums (`phys_gpu::ewald`); `--cpu` for neither. The work is
//! `phys::ewaldcheck::run`.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let flag = |f: &str| args.iter().any(|a| a == f);
    let (cpu, recip) = (flag("--cpu"), flag("--recip"));
    let args: Vec<String> = args.into_iter().filter(|a| !a.starts_with("--")).collect();
    if cpu {
        println!("Ewald step on the CPU");
    } else if recip {
        match phys_gpu::install_ewald() {
            Ok(n) => println!("Ewald reciprocal sums on the {n}"),
            Err(e) => println!("no GPU ({e}); the CPU's"),
        }
    } else {
        match phys_gpu::install_ewald_step() {
            Ok(n) => println!("Ewald step on the {n}"),
            Err(e) => println!("no GPU ({e}); the CPU's"),
        }
    }
    phys::ewaldcheck::run(args);
}
