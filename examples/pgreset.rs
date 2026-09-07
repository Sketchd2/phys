//! Drop and rebuild a PostgreSQL world at the current schema.
//!
//! `PostgresStore::connect` refuses a database built against a different layout
//! rather than misreading it. There is no migration while the project is
//! pre-alpha, so this is the whole answer: drop the tables and rebuild the
//! world. Destructive by design.
//!
//! ```sh
//! PHYS_PG='host=127.0.0.1 port=5432 user=phys dbname=phys' \
//!   cargo run --release --features postgres --example pgreset
//! ```
#[cfg(not(feature = "postgres"))]
fn main() {
    eprintln!("build with --features postgres");
    std::process::exit(2);
}

#[cfg(feature = "postgres")]
fn main() {
    let Ok(url) = std::env::var("PHYS_PG") else {
        eprintln!("set PHYS_PG to the database to reset");
        std::process::exit(2);
    };
    match phys::store_pg::PostgresStore::reset(&url) {
        Ok(_) => println!("reset to schema version {}", phys::store_pg::SCHEMA_VERSION),
        Err(e) => {
            eprintln!("reset failed: {e}");
            std::process::exit(1);
        }
    }
}
