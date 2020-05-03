use crate::field::FieldElement;

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
