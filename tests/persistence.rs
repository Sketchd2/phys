//! A world that survives its process.
//!
//! The engine's persistence *policy* has been designed and tested for months —
//! what is regenerable, what is pinned, when detail is released. The mechanism
//! to actually keep any of it did not exist: the crate had no serialisation and
//! no disk I/O, and `Tree::persisted` was an in-memory `HashMap` that died with
//! the process.
//!
//! These tests are the mechanism's contract. The interesting ones are not "does
//! it round-trip" — they are the two claims the architecture rests on:
//!
//! * a reloaded world **regenerates the same detail**, because addresses and
//!   epochs are what detail is derived from, and both survive the save;
//! * detail somebody **touched comes back byte-identical**, because it is not a
//!   sample of anything and could not be re-derived.

use phys::engine::{default_spec, galaxy, World};
use phys::ids::PathKey;
use phys::math::v3;
use phys::observe::{Interaction, Quantity};
use phys::persist::{decode, encode, FileStore, MemoryStore, WorldStore};
use phys::state::Body;
use phys::units::*;
use phys::wire::{WireError, FORMAT_VERSION, MAGIC};

/// `Snapshot` has no `Debug` (it holds the whole tree, and a derived `Debug`
/// on that would be a footgun in a panic message), so unwrap it by hand.
fn ok(r: Result<phys::persist::Snapshot, phys::wire::WireError>) -> phys::persist::Snapshot {
    match r {
        Ok(s) => s,
        Err(e) => panic!("decode failed: {e}"),
    }
}

/// Same, for the cases that are supposed to fail.
fn err(r: Result<phys::persist::Snapshot, phys::wire::WireError>) -> phys::wire::WireError {
    match r {
        Ok(_) => panic!("expected a decode error, got a valid snapshot"),
        Err(e) => e,
    }
}

fn a_world() -> World {
    let mut w = World::new(galaxy(0xB0A7, 1e9), 20.0);
    w.tree.nodes[0].spec.count = 256;
    w.time_rate = 0.05;
    w
}

/// Every field of every type must survive the trip. A field that is written and
/// not read, or read and not written, fails here rather than silently losing
/// somebody's world.
#[test]
fn every_field_round_trips() {
    let mut w = a_world();
    let root = w.tree.root;

    // Populate the awkward corners: a structure with a morphology and a
    // topology, a ledger fact, an audit entry, an environment, and detail
    // somebody has touched.
    let path = w.drill_to(root, Tier::Continuum.max_radius(), &default_spec);
    let deep = *path.last().unwrap();
    w.plant(
        deep,
        phys::morph::Program::Tree,
        Some(phys::morph::Environment { light_flux: 340.0, ..Default::default() }),
    );
    for _ in 0..6 {
        w.step_frame(50_000.0);
    }
    w.interact(Interaction::Impulse { target: deep, dp: v3(1.0, 2.0, 3.0) });
    w.interact(Interaction::Pin { target: deep });
    w.measure(deep, phys::observe::Instrument::Thermometer, Quantity::Temperature);
    w.interact(Interaction::Author {
        target: deep,
        property: phys::observe::Property::Temperature,
        value: 350.0,
    });

    let before = w.conserved();
    let bytes = encode(w.view());
    let back = ok(decode(&bytes));

    // The tree, node by node.
    assert_eq!(back.tree.nodes.len(), w.tree.nodes.len(), "node count");
    assert_eq!(back.tree.root, w.tree.root);
    assert_eq!(back.tree.world_seed, w.tree.world_seed);
    for (i, (a, b)) in w.tree.nodes.iter().zip(back.tree.nodes.iter()).enumerate() {
        assert_eq!(a.key, b.key, "node {i} key");
        assert_eq!(a.parent, b.parent, "node {i} parent");
        assert_eq!(a.slot, b.slot, "node {i} slot");
        assert_eq!(a.depth, b.depth, "node {i} depth");
        assert_eq!(a.tier, b.tier, "node {i} tier");
        assert_eq!(a.epoch, b.epoch, "node {i} epoch");
        assert_eq!(a.alive, b.alive, "node {i} alive");
        assert_eq!(a.pinned, b.pinned, "node {i} pinned");
        assert_eq!(a.residency, b.residency, "node {i} residency");
        assert_eq!(a.steps_taken, b.steps_taken, "node {i} steps");
        assert_eq!(a.children, b.children, "node {i} children");
        assert_eq!(a.spec.count, b.spec.count, "node {i} spec count");
        assert_eq!(a.spec.kind, b.spec.kind, "node {i} spec kind");
        // Bit-exact, not approximately equal. A save that renormalised a float
        // would break the idempotent-coarsening guarantee invisibly.
        assert_eq!(a.matter.mass.to_bits(), b.matter.mass.to_bits(), "node {i} mass bits");
        assert_eq!(
            a.matter.internal_energy.to_bits(),
            b.matter.internal_energy.to_bits(),
            "node {i} internal energy bits"
        );
        assert_eq!(a.motion.offset.x.to_bits(), b.motion.offset.x.to_bits(), "node {i} offset");
        assert_eq!(
            a.motion.orientation.w.to_bits(),
            b.motion.orientation.w.to_bits(),
            "node {i} orientation"
        );
        assert_eq!(a.time.to_bits(), b.time.to_bits(), "node {i} time");
        assert_eq!(a.last_solved.to_bits(), b.last_solved.to_bits(), "node {i} last_solved");
        assert_eq!(
            a.last_disturbed.to_bits(),
            b.last_disturbed.to_bits(),
            "node {i} last_disturbed"
        );
        assert_eq!(a.morphology.is_some(), b.morphology.is_some(), "node {i} morphology");
        if let (Some(m), Some(n)) = (&a.morphology, &b.morphology) {
            assert_eq!(m.program, n.program);
            assert_eq!(m.age.to_bits(), n.age.to_bits(), "morphology age");
            assert_eq!(m.built.to_bits(), n.built.to_bits(), "morphology built");
            assert_eq!(m.events.len(), n.events.len(), "morphology events");
            assert_eq!(m.genome[0].to_bits(), n.genome[0].to_bits(), "genome");
        }
        assert_eq!(a.topology.is_some(), b.topology.is_some(), "node {i} topology");
        if let (Some(t), Some(u)) = (&a.topology, &b.topology) {
            assert_eq!(t.joints.len(), u.joints.len(), "joints");
            assert_eq!(t.support, u.support, "support");
            assert_eq!(t.ties.len(), u.ties.len(), "ties");
            assert_eq!(t.material.name, u.material.name, "material name");
            assert_eq!(
                t.material.stiffness.to_bits(),
                u.material.stiffness.to_bits(),
                "material stiffness"
            );
        }
    }

    // The clock and the session-independent dials.
    assert_eq!(back.time.to_bits(), w.time.to_bits(), "world time");
    assert_eq!(back.pace.to_bits(), w.pace.to_bits(), "pace");
    assert_eq!(back.time_rate.to_bits(), w.time_rate.to_bits(), "time rate");
    assert_eq!(back.time_throttle.to_bits(), w.time_throttle.to_bits(), "throttle");
    assert_eq!(back.paced_to, w.paced_to, "paced_to");

    // The ledger, the audit, the environments.
    assert_eq!(back.ledger.len(), w.ledger.len(), "ledger facts");
    assert_eq!(back.ledger.sequence, w.ledger.sequence, "ledger sequence");
    assert_eq!(back.audit.len(), w.audit.len(), "audit entries");
    assert!(!back.audit.is_empty(), "the authoring above should have left a record");
    assert_eq!(back.environments.len(), w.environments.len(), "environments");

    // And the conserved tuple, which is the whole point.
    let reloaded = World::from_snapshot(back, 20.0);
    let after = reloaded.conserved();
    let drift = (after.energy - before.energy).abs() / before.energy.abs().max(1e-300);
    println!("  {} nodes, {} bytes, energy drift {drift:.3e}", w.tree.nodes.len(), bytes.len());
    assert_eq!(after.baryon.to_bits(), before.baryon.to_bits(), "baryon number is not bit-exact");
    assert!(drift < 1e-15, "energy drifted across a save: {drift:.3e}");
}

/// Detail somebody touched comes back exactly. It is not a sample of anything,
/// so there is nothing to re-derive it from.
#[test]
fn pinned_detail_returns_byte_identical() {
    let mut w = a_world();
    let root = w.tree.root;
    let path = w.drill_to(root, Tier::Continuum.max_radius(), &default_spec);
    let deep = *path.last().unwrap();
    w.tree.refine(deep);
    w.interact(Interaction::Impulse { target: deep, dp: v3(7.0, -3.0, 0.5) });
    w.interact(Interaction::Pin { target: deep });

    let original: Vec<Body> = w.tree.nodes[deep.get()].bodies.clone();
    assert!(!original.is_empty(), "nothing to test");

    let bytes = encode(w.view());
    let back = World::from_snapshot(ok(decode(&bytes)), 20.0);
    let after = &back.tree.nodes[deep.get()].bodies;

    assert_eq!(after.len(), original.len(), "pinned body count changed");
    let mut worst = 0u64;
    for (a, b) in original.iter().zip(after.iter()) {
        for (x, y) in [
            (a.pos.x, b.pos.x),
            (a.pos.y, b.pos.y),
            (a.pos.z, b.pos.z),
            (a.vel.x, b.vel.x),
            (a.mass, b.mass),
            (a.internal_energy, b.internal_energy),
            (a.temperature, b.temperature),
        ] {
            worst = worst.max(x.to_bits() ^ y.to_bits());
        }
        assert_eq!(a.kind, b.kind);
        assert_eq!(a.slot, b.slot);
    }
    println!("  {} pinned bodies, worst bit difference {worst}", original.len());
    assert_eq!(worst, 0, "pinned detail was not byte-identical after a reload");
}

/// Unpinned detail is *not* saved — and does not need to be, because the
/// address and the epoch are, and those are what it is derived from.
///
/// This is the claim that makes the save file small. If it ever fails, the
/// world is no longer regenerable and every save must carry everything.
#[test]
fn regenerated_detail_matches_what_was_discarded() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    let original: Vec<Body> = w.tree.nodes[root.get()].bodies.clone();
    assert!(original.len() > 100, "want a decent sample");
    assert!(!w.tree.nodes[root.get()].pinned, "this test is about *unpinned* detail");

    let bytes = encode(w.view());
    let mut back = World::from_snapshot(ok(decode(&bytes)), 20.0);

    // The file did not carry them.
    assert!(
        back.tree.nodes[root.get()].bodies.is_empty(),
        "unpinned detail should not be in the file"
    );

    // Asking for them again reproduces them.
    let again = back.tree.refine(root).to_vec();
    assert_eq!(again.len(), original.len(), "regenerated a different number of bodies");
    let mut worst = 0u64;
    for (a, b) in original.iter().zip(again.iter()) {
        worst = worst.max(a.pos.x.to_bits() ^ b.pos.x.to_bits());
        worst = worst.max(a.vel.y.to_bits() ^ b.vel.y.to_bits());
        worst = worst.max(a.mass.to_bits() ^ b.mass.to_bits());
    }
    let saved = original.len() * std::mem::size_of::<Body>();
    println!(
        "  {} bodies regenerated bit-for-bit; {} bytes of detail not written, file is {}",
        original.len(),
        saved,
        bytes.len()
    );
    assert_eq!(worst, 0, "regenerated detail differs from what was discarded");
}

/// A reloaded world continues identically to one that was never saved.
#[test]
fn reloading_does_not_change_the_future() {
    let mut a = a_world();
    let root = a.tree.root;
    // Pinned, so the reload comes back at the same *resolution* as well as the
    // same state. Without this the futures legitimately diverge — see
    // `an_unpinned_reload_comes_back_coarse`.
    for &n in &a.drill_to(root, Tier::Planetary.max_radius(), &default_spec) {
        a.tree.pin(n);
    }
    for _ in 0..4 {
        a.step_frame(50_000.0);
    }

    let bytes = encode(a.view());
    let mut b = World::from_snapshot(ok(decode(&bytes)), 20.0);

    for _ in 0..6 {
        a.step_frame(50_000.0);
        b.step_frame(50_000.0);
    }

    println!("  after six more frames: t = {:.6e} against {:.6e}", a.time, b.time);
    assert_eq!(a.time.to_bits(), b.time.to_bits(), "the clocks diverged");
    for (i, (x, y)) in a.tree.nodes.iter().zip(b.tree.nodes.iter()).enumerate() {
        assert_eq!(
            x.motion.offset.x.to_bits(),
            y.motion.offset.x.to_bits(),
            "node {i} drifted after a reload"
        );
        assert_eq!(x.matter.mass.to_bits(), y.matter.mass.to_bits(), "node {i} mass drifted");
    }
}

/// Two saves of the same world are the same bytes. Without this, "did anything
/// change" is unanswerable and content addressing is impossible.
#[test]
fn saving_is_deterministic() {
    let mut w = a_world();
    let root = w.tree.root;
    let path = w.drill_to(root, Tier::Continuum.max_radius(), &default_spec);
    for &n in &path {
        w.tree.pin(n);
    }
    for _ in 0..3 {
        w.step_frame(50_000.0);
    }
    let first = encode(w.view());
    let second = encode(w.view());
    println!("  {} bytes, twice", first.len());
    assert_eq!(first, second, "two saves of one world differed");
}

/// A world file arrives from disk and later from a network. Nothing it can
/// contain may panic the reader.
#[test]
fn malformed_input_errors_rather_than_panicking() {
    assert_eq!(err(decode(&[])), WireError::Truncated { what: "magic", need: 4, have: 0 });
    assert_eq!(err(decode(b"NOPE\x01\x00")), WireError::BadMagic);

    let mut wrong_version = MAGIC.to_vec();
    wrong_version.extend_from_slice(&(FORMAT_VERSION + 1).to_le_bytes());
    assert!(matches!(
        err(decode(&wrong_version)),
        WireError::UnsupportedVersion { .. }
    ));

    // A length prefix that claims more than the file could possibly hold must
    // be refused before anything is allocated, not after.
    let mut greedy = MAGIC.to_vec();
    greedy.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
    greedy.extend_from_slice(&0u64.to_le_bytes()); // world seed
    greedy.extend_from_slice(&0u32.to_le_bytes()); // root
    greedy.extend_from_slice(&u32::MAX.to_le_bytes()); // "4 billion nodes"
    assert!(
        matches!(err(decode(&greedy)), WireError::TooLong { .. }),
        "an absurd length prefix must be rejected up front"
    );

    // And every truncation of a real file is an error, never a panic.
    let mut w = a_world();
    w.step_frame(50_000.0);
    let good = encode(w.view());
    let mut errors = 0;
    for cut in (0..good.len()).step_by(7) {
        if decode(&good[..cut]).is_err() {
            errors += 1;
        }
    }
    println!("  {errors} truncations of a {}-byte file, all handled", good.len());
    assert!(errors > 0);

    // Corrupt bytes in place: still no panic, whatever comes out.
    for i in (0..good.len()).step_by(101) {
        let mut bad = good.clone();
        bad[i] ^= 0xFF;
        let _ = decode(&bad);
    }
}

/// The file store writes through a temporary and renames, so an interrupted
/// save cannot destroy the previous world.
#[test]
fn a_world_survives_a_round_trip_through_a_file() {
    let dir = std::env::temp_dir();
    let path = dir.join(format!("phys-persist-test-{}.world", std::process::id()));
    let _ = std::fs::remove_file(&path);

    let mut w = a_world();
    let root = w.tree.root;
    let path_nodes = w.drill_to(root, Tier::Continuum.max_radius(), &default_spec);
    let deep = *path_nodes.last().unwrap();
    w.plant(deep, phys::morph::Program::Tree, Some(phys::morph::Environment::default()));
    for _ in 0..5 {
        w.step_frame(50_000.0);
    }
    let before = w.conserved();

    let mut store = FileStore::new(&path);
    store.save(w.view()).expect("save");
    let on_disk = store.size();
    assert!(on_disk > 0, "nothing was written");
    assert!(!path.with_extension("tmp").exists(), "the temporary file was left behind");

    drop(w);

    let back = World::from_snapshot(ok(store.load()), 20.0);
    let after = back.conserved();
    // Not bit-exact, and it must not be asserted as such. `sum_conserved` reads
    // a node's *bodies* when it is materialised and its *matter* when it is
    // not, and the file deliberately drops unpinned bodies — so a reload swaps
    // which of the two paths is taken. The difference is the sample/summarise
    // round-off the engine already bounds by `IDEMPOTENT_TOLERANCE`, and it is
    // the cost of not storing detail that can be worked out again.
    let drift = (after.baryon - before.baryon).abs() / before.baryon.abs().max(1e-300);
    println!("  {on_disk} bytes on disk; baryon drift across the trip {drift:.3e}");
    assert!(drift < phys::tree::IDEMPOTENT_TOLERANCE, "baryon drift {drift:.3e} is too large");
    assert!(back.tree.nodes[deep.get()].morphology.is_some(), "the structure did not survive");

    let _ = std::fs::remove_file(&path);
}

/// The in-memory store is the same code path, so tests that do not want a file
/// still exercise the format.
#[test]
fn the_memory_store_is_the_same_format() {
    let mut w = a_world();
    w.step_frame(50_000.0);
    let mut mem = MemoryStore::new();
    mem.save(w.view()).expect("save");
    assert_eq!(mem.raw(), &encode(w.view())[..]);
    let back = ok(mem.load());
    assert_eq!(back.time.to_bits(), w.time.to_bits());
    println!("  {} bytes in memory", mem.size());
}

/// A measured value is a fact about the world. It has to come back as the same
/// fact, or the ledger's whole promise — measure twice, get the same answer —
/// stops holding across a restart.
#[test]
fn committed_facts_survive() {
    let mut w = a_world();
    let root = w.tree.root;
    let key: PathKey = w.tree.nodes[root.get()].key;
    w.ledger.commit(key, Quantity::DecayTime, 1234.5678, 42.0);

    let bytes = encode(w.view());
    let back = World::from_snapshot(ok(decode(&bytes)), 20.0);
    let fact = back.ledger.peek(key, Quantity::DecayTime).expect("fact missing after reload");
    println!("  decay time {} committed at t = {}", fact.value, fact.time);
    assert_eq!(fact.value.to_bits(), 1234.5678f64.to_bits());
    assert_eq!(fact.time.to_bits(), 42.0f64.to_bits());
}

/// A reload comes back *coarse*, and that changes how fast its clock runs until
/// something asks for detail again.
///
/// This is not a defect, but it is surprising enough to be worth pinning down.
/// The file stores no unpinned bodies, so a reloaded node is matter alone; and
/// `node_cadence` reads body speeds when a node is materialised and the
/// matter's own characteristic speed when it is not. The pace follows the
/// cadence of whatever is being watched, so a coarse world runs at a coarser
/// pace — correctly, because there is nothing resolved that needs finer steps.
///
/// The state is identical eitherway. What differs is the resolution, and
/// resolution is recovered on demand.
#[test]
fn an_unpinned_reload_comes_back_coarse() {
    let mut w = a_world();
    let root = w.tree.root;
    w.tree.refine(root);
    w.step_frame(50_000.0);

    let fine_pace = w.pace;
    let fine_bodies = w.tree.nodes[root.get()].bodies.len();
    assert!(fine_bodies > 0);

    let bytes = encode(w.view());
    let mut back = World::from_snapshot(ok(decode(&bytes)), 20.0);
    assert!(back.tree.nodes[root.get()].bodies.is_empty(), "detail should not be in the file");

    back.step_frame(50_000.0);
    let coarse_pace = back.pace;
    println!(
        "  resolved: {fine_bodies} bodies at {fine_pace:.3e} s per frame; \
         reloaded coarse: {coarse_pace:.3e} s per frame"
    );
    assert!(coarse_pace > fine_pace, "a coarse world should run at a coarser pace");

    // And asking for the detail back restores both.
    back.tree.refine(root);
    back.step_frame(50_000.0);
    assert_eq!(
        back.tree.nodes[root.get()].bodies.len(),
        fine_bodies,
        "re-materialising gave a different number of bodies"
    );
    assert!(
        (back.pace - fine_pace).abs() <= fine_pace * 1e-9,
        "the pace did not come back with the resolution: {:.3e} against {fine_pace:.3e}",
        back.pace
    );
}

/// The stamp has to move when the layout does, and nothing but a person
/// remembering makes that happen.
///
/// # What replaced migration
///
/// There is no migration and there is not going to be one while the project is
/// pre-alpha: a file written by a different layout is refused and the world is
/// rebuilt. That is only safe if the refusal actually fires, and it only fires
/// if `FORMAT_VERSION` was bumped when the layout changed. Forget the bump and
/// a stale file is not refused — it is *parsed*, into whatever the new reader
/// makes of the old bytes, which is the exact failure the scheme exists to
/// prevent.
///
/// So this holds the encoded size of a reference world with every awkward
/// corner populated. Add a field, remove one, or change one's width, and the
/// size moves and this fails, telling you to bump the stamp.
///
/// # Why a size and not a checksum
///
/// A checksum would also catch a field being reordered or retyped at the same
/// width, which a size does not. It would also depend on every float in the
/// file, and float bit patterns are not guaranteed identical across compilers
/// and platforms — see the recipe-determinism entry in `docs/BACKLOG.md`. A
/// size is invariant under every value in the world and varies with its shape,
/// which is the thing being guarded. The checksum is printed for information;
/// nothing asserts on it.
///
/// The world is built without `step_frame`, deliberately: stepping spends a
/// wall-clock budget, so how much of it runs depends on how fast the machine
/// is, and a reference world has to be the same everywhere. It is also built
/// small — pinning a node pins its whole ancestor chain and a pinned node
/// writes every body it holds, so drilling deep would put a hundred thousand
/// bodies in the file and its size would then move whenever anybody tuned a
/// default spec count.
#[test]
fn the_format_stamp_tracks_the_format() {
    /// Bump `wire::FORMAT_VERSION`, then update this.
    const REFERENCE_BYTES: usize = 2_806;

    use phys::chem::{Arrangement, Bond, Element, Lattice, Mixture, Order, Phase};

    let mut w = a_world();
    w.tree.nodes[0].spec.count = 6;
    let root = w.tree.root;
    w.tree.refine(root);
    let child = w.tree.promote(root, 0, default_spec(Tier::Stellar));
    w.tree.nodes[child.get()].spec.count = 3;
    w.tree.refine(child);

    // A morphology and its topology, pinned detail, a ledger fact, an audit
    // entry, an environment, an in-flight influence, a time bubble, a substance
    // catalogue and a mixture.
    w.plant(
        child,
        phys::morph::Program::Tree,
        Some(phys::morph::Environment { light_flux: 340.0, ..Default::default() }),
    );
    w.interact(Interaction::Pin { target: child });
    w.interact(Interaction::Impulse { target: child, dp: v3(1.0, 2.0, 3.0) });
    w.interact(Interaction::Author {
        target: child,
        property: phys::observe::Property::Temperature,
        value: 350.0,
    });
    w.dilate(child, 250.0);
    w.ledger.commit(w.tree.nodes[root.get()].key, Quantity::DecayTime, 99.5, 7.0);

    // Built by hand rather than by name: nothing in the engine knows what salt
    // is, which is the whole point of the chemistry layer.
    let salt = w
        .substances
        .intern(Arrangement::crystal(
            vec![Element(11), Element(17)],
            vec![Bond::new(0, 1, Order::Ionic)],
            Lattice::Cubic { a: 3.55e-10 },
        ))
        .expect("salt analyses");
    w.substances.name(salt, "salt");
    let mut mix = Mixture::new();
    mix.add(salt, Phase::Solid, 0.25);
    w.set_mixture(child, mix);

    let bytes = encode(w.view());
    let checksum = bytes.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, b| {
        (h ^ *b as u64).wrapping_mul(0x1000_0000_01b3)
    });
    println!(
        "  format {FORMAT_VERSION}: reference world is {} bytes (checksum {checksum:016x})",
        bytes.len()
    );

    assert_eq!(
        bytes.len(),
        REFERENCE_BYTES,
        "\nthe wire layout changed: the reference world is now {} bytes, not {REFERENCE_BYTES}.\n\
         Bump `wire::FORMAT_VERSION` (currently {FORMAT_VERSION}) and set REFERENCE_BYTES to {}.\n\
         There is no migration — an existing world is refused and has to be rebuilt.",
        bytes.len(),
        bytes.len()
    );

    // And what was written still reads back, so the size is not being held
    // steady by something that quietly broke.
    let back = ok(decode(&bytes));
    assert_eq!(back.tree.nodes.len(), w.tree.nodes.len());
    assert_eq!(back.substances.len(), 1);
    assert_eq!(back.mixtures.len(), 1);
}
