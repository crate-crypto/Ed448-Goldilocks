// radix 28
pub mod u32;
pub use crate::field::u32::Scalar;

// Following Dalek, this should be behind a feature gate and would be the default radix
pub type FieldElement = crate::field::u32::Fq;
