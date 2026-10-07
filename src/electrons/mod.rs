//! Electrons in molecules: the electronic-structure route of `docs/PLAY.md`
//! Phase 6, item 2.
//!
//! What a substance does in bulk — how it boils, how dense it is, how it bends
//! light — is set by its electrons, and the engine had no way to ask them. This
//! module is that way: density-functional theory, solved once per substance when
//! it is interned and stored as the substance's shortcut, from the Schrodinger
//! equation, a non-empirical functional, and a basis the engine derives itself.
//!
//! Atomic units throughout (bohr, hartree), converted at the boundary.

pub mod atom;
pub mod basis;
pub mod element;
pub mod boys;
pub mod functional;
pub mod gradient;
pub mod grow;
pub mod hf;
pub mod grid;
pub mod integrals;
pub mod linalg;
pub mod molecule;
pub mod partition;
pub mod scf;
pub mod values;
pub mod vdw;
