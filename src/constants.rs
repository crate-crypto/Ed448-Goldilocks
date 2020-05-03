use crate::curve::twedwards::extended::ExtendedPoint;
use crate::decaf::DecafPoint;
use crate::field::FieldElement;

#[cfg(feature = "u32_backend")]
pub(crate) use crate::field::u32::constants::*;

// / Since each backend has a specific set of constants
// / Following Dalek, we will have to use a feature to import the correct constants for the specific backend.
// /
// Refactor this to use the Edwards basepoint internally
pub(crate) const DECAF_BASEPOINT: DecafPoint = DecafPoint(ExtendedPoint {
    X: FieldElement::from_raw_slice([
        268435456, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 134217727,
        268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 134217727,
    ]),
    Y: FieldElement::from_raw_slice([
        266160740, 101161805, 74490312, 12706731, 149232027, 72184820, 68425752, 84169329,
        64300076, 80170041, 105082960, 37781586, 19953866, 222875756, 82854534, 139496929,
    ]),
    Z: FieldElement::from_raw_slice([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    T: FieldElement::from_raw_slice([
        202998021, 238846317, 66379923, 102789507, 54662147, 81652110, 85576069, 171023191,
        104342404, 127188629, 141403663, 236837931, 109226495, 84812757, 24364708, 114517662,
    ]),
});
