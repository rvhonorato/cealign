//! # cealign
//!
//! A Rust implementation of the [Combinatorial Extension (CE) algorithm](https://doi.org/10.1093/protein/11.9.739)
//! for protein structure alignment.
//!
//! ## Quick start
//!
//! ```no_run
//! let (mobile, _) = pdbtbx::open("mobile.pdb").unwrap();
//! let (target, _) = pdbtbx::open("target.pdb").unwrap();
//!
//! let (_aligned_mobile, _reference, rmsd, n_aligned) = cealign::ce::align(mobile, target, false);
//! println!("{:.3} Å over {} residues", rmsd, n_aligned);
//! ```

pub mod ce;
pub mod structure;
pub mod utils;
#[cfg(feature = "plot")]
pub mod visualization;
