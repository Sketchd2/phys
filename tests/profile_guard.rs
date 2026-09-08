//! The test profile is optimised. These assert that optimising it did not
//! quietly turn off the checks the project relies on.
//!
//! `opt-level = 2` was set on `[profile.test]` because the suite took twenty
//! minutes unoptimised. The danger of that change is silent: `overflow-checks`
//! defaults *off* in an optimised profile, and an integer overflow that used to
//! panic in a test run would start wrapping instead — which is exactly the
//! failure mode that hid the neighbour-grid bug in release for as long as it
//! did. So the profile sets them explicitly, and this asserts they took.

/// Integer overflow must still panic, not wrap.
#[test]
#[should_panic(expected = "attempt to add with overflow")]
fn overflow_checks_are_on() {
    let mut x = i64::MAX;
    // Through a black box so the optimiser cannot fold it at compile time.
    x += std::hint::black_box(1i64);
    println!("{x}");
}

/// `debug_assert!` must still run.
#[test]
fn debug_assertions_are_on() {
    assert!(
        cfg!(debug_assertions),
        "debug_assertions is off in the test profile: every debug_assert! in the \
         engine is now dead code, including the ones guarding node invariants"
    );
}
