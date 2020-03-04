pub(crate) mod base;
pub(crate) mod scalar;

// Components not in the field module, can only access the two defined Field structs
// They should not have access to the Limbs etc and the internal structure
pub use base::Fq;
pub use scalar::Scalar;
