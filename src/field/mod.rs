// radix 28
#[cfg(feature = "u32_backend")]
pub mod u32;

pub use crate::field::u32::Scalar;

#[cfg(feature = "u32_backend")]
pub type FieldElement = crate::field::u32::Fq;
