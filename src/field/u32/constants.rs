// Move all curve specific constants into this file
// as they will have different representations depending on the backend
use crate::field::u32::Fq;

// -4 * Twised_D = -4 * (EDWARDS_D-1)
pub const NEG_FOUR_TIMES_TWISTED_D: Fq = Fq([
    156327, 268435456, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435454,
    268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
]);

// TODO: Differentiate between Twisted edwards d(a=-1) which is -39082 and edwards d(a=1) which is -39081

/// Edwards `d`, equals to -39081
pub const EDWARDS_D: Fq = Fq([
    268396374, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
    268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
]);

/// Neg_Edwards `-d`, equals to 39081
pub const NEG_EDWARDS_D: Fq = Fq([39081, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

/// Edwards `d-1`, equals to -39082
pub const D_MINUS_ONE: Fq = Fq([
    268396373, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
    268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
]);

/// Edwards `2(d-1)`, equals to -78164
pub const TWO_D_MINUS_ONE: Fq = Fq([
    268357291, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
    268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
]);

/// Edwards `1-d`, equals to 39082
pub const ONE_MINUS_D: Fq = Fq([39082, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

/// Edwards `2(1-d)`, equals to 78164
pub const TWO_ONE_MINUS_D: Fq = Fq([78164, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

/// INVSQRT(a*d-1)
pub const DECAF_FACTOR: Fq = Fq([
    0x05572736, 0x042ef0f4, 0x00ce5296, 0x07bf6aa2, 0x0ed26033, 0x0f4fd6ed, 0x0a839a66, 0x0968c14b,
    0x04a2d780, 0x0b8d54b6, 0x01a7b8a5, 0x06aa0a1f, 0x0d722fa2, 0x0683bf68, 0x0beb24f7, 0x022d962f,
]);
/// 39082 used in the doubling procedure in montgomery ladder
pub const A_PLUS_TWO_OVER_FOUR: Fq = ONE_MINUS_D;
