// This will be the module for Decaf over Ed448
// This is the newer version of the Decaf strategy, which looks simpler

pub mod constants;
pub mod decaf;

pub use decaf::{CompressedDecaf, DecafPoint};
