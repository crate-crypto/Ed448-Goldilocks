use crate::curve::twedwards::extended::ExtendedPoint;
use crate::decaf::DecafPoint;
use crate::field::FieldElement;

// / Since each backend has a specific set of constants
// / Following Dalek, we will have to use a feature to import the correct constants for the specific backend.
// /
// Refactor this to use the Edwards basepoint internally
pub const DECAF_BASEPOINT: DecafPoint = DecafPoint(ExtendedPoint {
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

/// INVSQRT(a*d-1)
pub const DECAF_FACTOR: FieldElement = FieldElement::from_raw_slice([
    0x05572736, 0x042ef0f4, 0x00ce5296, 0x07bf6aa2, 0x0ed26033, 0x0f4fd6ed, 0x0a839a66, 0x0968c14b,
    0x04a2d780, 0x0b8d54b6, 0x01a7b8a5, 0x06aa0a1f, 0x0d722fa2, 0x0683bf68, 0x0beb24f7, 0x022d962f,
]);

// -4 * Twised_D = -4 * (EDWARDS_D-1)
pub const NEG_FOUR_TIMES_TWISTED_D: FieldElement = FieldElement::from_raw_slice([
    156327, 268435456, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435454,
    268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
]);

// TODO: Differentiate between Twisted edwards d(a=-1) which is -39082 and edwards d(a=1) which is -39081

/// Edwards `d`, equals to -39081
pub const EDWARDS_D: FieldElement = FieldElement::from_raw_slice([
    268396374, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
    268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
]);

/// Neg_Edwards `-d`, equals to 39081
pub const NEG_EDWARDS_D: FieldElement =
    FieldElement::from_raw_slice([39081, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

/// Edwards `d-1`, equals to -39082
pub const D_MINUS_ONE: FieldElement = FieldElement::from_raw_slice([
    268396373, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
    268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
]);

/// Edwards `2(d-1)`, equals to -78164
pub const TWO_D_MINUS_ONE: FieldElement = FieldElement::from_raw_slice([
    268357291, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
    268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
]);

/// Edwards `1-d`, equals to 39082
pub const ONE_MINUS_D: FieldElement =
    FieldElement::from_raw_slice([39082, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

/// Edwards `2(1-d)`, equals to 78164
pub const TWO_ONE_MINUS_D: FieldElement =
    FieldElement::from_raw_slice([78164, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
