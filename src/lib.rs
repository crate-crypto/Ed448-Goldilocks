#![forbid(unsafe_code)]
// XXX: Change this to deny later on
#![warn(unused_attributes, unused_imports, unused_mut, unused_must_use)]
#![allow(non_snake_case)]

// Internal macros. Must come first!
#[macro_use]
pub(crate) mod macros;

// As usual, we will use this file to carefully define the API/ what we expose to the user
pub mod constants;
pub mod curve;
pub mod decaf;
mod field;
pub mod ristretto;

pub use field::Scalar;
