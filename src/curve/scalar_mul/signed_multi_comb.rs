use super::window::wnaf;
use super::window::wnaf::LookupTable;
use crate::curve::twedwards::{extended::ExtendedPoint, extensible::ExtensiblePoint};
use crate::field::Scalar;
use subtle::{Choice, ConditionallyNegatable};

const SCALAR_ADJUSTMENT_FACTOR: Scalar = Scalar([
    0x4a7bb0cf, 0xc873d6d5, 0x23a70aad, 0xe933d8d7, 0x129c96fd, 0xbb124b65, 0x335dc163, 0x00000008,
    0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
]);

// XXX: This is _almost_ unreadable, we will slowly refactor for interoperability
pub fn signed_multi_comb(point: &ExtendedPoint, s: &Scalar) -> ExtendedPoint {
    let mut result = ExtensiblePoint::identity();
    // Recode Scalar
    // XXX: Better to recode to radix_16 as I don't think that this strategy would have a significant speedup?
    let mut scalar_one_x = *s + SCALAR_ADJUSTMENT_FACTOR;
    scalar_one_x = scalar_one_x.halve();

    let lookup = LookupTable::from(point);

    let loop_end = 446 - (445 % wnaf::WINDOW) - 1;
    let loop_start = 0;

    for i in (loop_start..=loop_end).rev().step_by(wnaf::WINDOW) {
        let mut bits = scalar_one_x[i / 32] >> (i % 32);
        if (i % 32 >= 32 - wnaf::WINDOW) && (i / 32 < (14 - 1)) {
            bits ^= scalar_one_x[(i / 32) + 1] << (32 - (i % 32));
        }
        bits &= wnaf::WINDOW_MASK as u32;
        let inv = (bits >> (wnaf::WINDOW - 1)).wrapping_sub(1);
        bits ^= inv;

        let mut neg_P = lookup.select(bits & (wnaf::WINDOW_T_MASK as u32));
        neg_P.conditional_negate(Choice::from((inv & 1) as u8));

        for _ in 0..wnaf::WINDOW {
            result = result.double();
        }
        result = result.add_projective_niels(&neg_P);
    }
    result.to_extended()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::constants::TWISTED_EDWARDS_BASE_POINT;
    use crate::curve::scalar_mul::double_and_add;

    #[test]
    fn test_scalar_mul() {
        // XXX: In the future use known multiples from Sage in bytes form?
        let twisted_point = TWISTED_EDWARDS_BASE_POINT;
        let scalar = Scalar([
            0x6ee372b7, 0xe128ae78, 0x1533427c, 0xad0b7015, 0x307f665e, 0xde8026c1, 0xb64629d1,
            0xab454c66, 0x3fe5bf1a, 0x083f8304, 0x3c003777, 0xdef437f6, 0xee2e1b73, 0x05ca185a,
        ]);

        let got = signed_multi_comb(&twisted_point, &scalar);

        let got2 = double_and_add(&twisted_point, &scalar);
        assert_eq!(got, got2);

        // Lets see if this is conserved over the isogenies
        let edwards_point = twisted_point.to_untwisted();
        let got_untwisted_point = edwards_point.scalar_mul(&scalar);
        let expected_untwisted_point = got.to_untwisted();
        assert_eq!(got_untwisted_point, expected_untwisted_point);
    }

    #[test]
    fn test_simple_scalar_mul_identities() {
        let x = TWISTED_EDWARDS_BASE_POINT;

        // Test that 1 * P = P
        let exp = signed_multi_comb(&x, &Scalar::from(1));
        assert!(x == exp);
        // Test that 2 * (P + P) = 4 * P
        let x_ext = x.to_extensible();
        let expected_two_x = x_ext.add_extensible(&x_ext).double();
        let got = signed_multi_comb(&x, &Scalar::from(4));
        assert!(expected_two_x.to_extended() == got);
    }
}
