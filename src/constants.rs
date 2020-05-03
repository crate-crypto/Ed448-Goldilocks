use crate::curve::twedwards::extended::ExtendedPoint;
use crate::decaf::DecafPoint;
use crate::field::FieldElement;

#[cfg(feature = "u32_backend")]
pub(crate) use crate::field::u32::constants::*;

pub(crate) const DECAF_BASEPOINT: DecafPoint = DecafPoint(TWISTED_EDWARDS_BASE_POINT);
