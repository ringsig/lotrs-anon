//! # lotrs — Practical Post-Quantum Threshold Ring Signatures from Lattices
//!
//! Rust implementation of the **LoTRS** scheme, test-vector compatible
//! with the Python reference in `../lotrs-py`.
//!
//! See the crate-level `README.md` for a status matrix of which layers
//! are complete.  The arithmetic stack ([`aux_ntt`], [`ring`], [`sample`],
//! [`codec`]) is fully implemented and tested against the Python
//! reference; the scheme layer (`lotrs`) is staged.

#![deny(unsafe_code)]
#![warn(rust_2018_idioms)]

pub mod aux_ntt;
pub mod cdt;
pub mod codec;
pub mod lotrs;
pub mod params;
pub mod ring;
pub mod sample;

pub use codec::{CodecError, LoTRSCodec};
pub use lotrs::{SignTimings, VerifyTimings};
pub use params::{
    LoTRSParams, MaskSamplerKind, BENCH_4OF32, BENCH_PARAMS, PRODUCTION_PARAMS, TEST_PARAMS,
};
pub use ring::{Poly, PolyMat, PolyVec, Ring};
pub use sample::{Tag, Xof};
