use super::window::wnaf;
use super::window::BASE_TABLE;
use crate::curve::twedwards::extended::ExtendedPoint;
use crate::curve::twedwards::extensible::ExtensiblePoint;
use crate::field::Scalar;
use subtle::Choice;

// XXX: This is _almost_ unreadable, we will slowly refactor for interoperability
pub fn signed_multi_comb(point: &ExtendedPoint, s: &Scalar) -> ExtendedPoint {
    let mut result = ExtensiblePoint::identity();
    // Recode Scalar
    // XXX: Better to recode to radix_16 as I don't think that this strategy would have a significant speedup?
    let mut scalar_one_x = *s + BASE_TABLE.scalar_adjustment;
    // XXX: Halve method does not seem to be working correctly, so we just invert 2 for now
    scalar_one_x = scalar_one_x * Scalar::from(2).invert();

    use wnaf::LookupTable;
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

        let mut neg_P = lookup.select(bits);
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
    use crate::curve::scalar_mul::double_and_add;
    use crate::field::FieldElement;
    #[test]
    fn test_scalr_mul() {
        let twisted_point = ExtendedPoint {
            X: FieldElement::from_raw_slice([
                0x0a7f964a, 0x0db033f8, 0x062b9f0b, 0x07bff7d6, 0x05e755a2, 0x013b6f8b, 0x0f080bdc,
                0x0a112ac0, 0x0416988a, 0x03404b2f, 0x00561ea3, 0x01df752c, 0x070e0b1c, 0x0e73a0c4,
                0x078245d5, 0x09a42df0,
            ]),
            Y: FieldElement::from_raw_slice([
                0x0c2e6c3d, 0x0a03c3f2, 0x0fd16e97, 0x0bab4ec6, 0x08ddba78, 0x091638ef, 0x0b0add85,
                0x070c212d, 0x04bcd337, 0x0c828579, 0x0712cfff, 0x09c1534a, 0x0119cafe, 0x08e72ee0,
                0x0f14ff19, 0x0d0c7e25,
            ]),
            Z: FieldElement::from_raw_slice([
                0x0a0d6be1, 0x0bcd9788, 0x00f9ca8a, 0x038cf839, 0x00912da2, 0x0a3c503a, 0x056fe7e0,
                0x03db9a49, 0x0f19d062, 0x052ac631, 0x01cbda35, 0x02967214, 0x0eed2db2, 0x0a948ce0,
                0x05f7a3a7, 0x0fa35bc2,
            ]),
            T: FieldElement::from_raw_slice([
                0x0fc9f32d, 0x0e442931, 0x065e50ff, 0x04be230d, 0x0dc923c2, 0x0000467c, 0x08fc8902,
                0x0e034cfb, 0x0126370c, 0x06ec706d, 0x06ff07ad, 0x0a27cd65, 0x060f214f, 0x0eb7756d,
                0x0b694dc7, 0x015705ad,
            ]),
        };
        let scalar = Scalar([
            0x6ee372b7, 0xe128ae78, 0x1533427c, 0xad0b7015, 0x307f665e, 0xde8026c1, 0xb64629d1,
            0xab454c66, 0x3fe5bf1a, 0x083f8304, 0x3c003777, 0xdef437f6, 0xee2e1b73, 0x05ca185a,
        ]);

        let expected_twisted_point = ExtendedPoint {
            X: FieldElement::from_raw_slice([
                0x08630007, 0x0bd755e6, 0x0f76b928, 0x070d9694, 0x0b952009, 0x0cf85b12, 0x0c3a6e9c,
                0x0e2d860e, 0x02fd2901, 0x09a73726, 0x02aa2d4c, 0x06913ea9, 0x090da66d, 0x06a5c6f1,
                0x04cc7a13, 0x0eb24ed8,
            ]),
            Y: FieldElement::from_raw_slice([
                0x0bb37152, 0x0a3a36b3, 0x0a720c7f, 0x0e29095f, 0x04e76cf4, 0x0cfad965, 0x07439798,
                0x0f4b7ba4, 0x0316ba61, 0x09389566, 0x07f96104, 0x07bdc39c, 0x0f019987, 0x05416850,
                0x0612c6c8, 0x0e231baa,
            ]),
            Z: FieldElement::from_raw_slice([
                0x0179c756, 0x04130eef, 0x07f43255, 0x0cc1534d, 0x03e347fd, 0x0c745e4d, 0x068d7bf5,
                0x020b8465, 0x0356d2f1, 0x069b22fd, 0x0b6cf87f, 0x0edf9761, 0x034f512f, 0x0411b43f,
                0x033f0755, 0x06195e97,
            ]),
            T: FieldElement::from_raw_slice([
                0x0866187a, 0x035622be, 0x0b9e2e78, 0x0cae26c6, 0x041c2c41, 0x07296c68, 0x03343d3e,
                0x062c0927, 0x0cf5d263, 0x08db465d, 0x033382d6, 0x0c5e6eff, 0x0c0ded8d, 0x037837bf,
                0x03780cc6, 0x0e2360df,
            ]),
        };
        let got = signed_multi_comb(&twisted_point, &scalar);
        assert!(expected_twisted_point == got);

        let got2 = double_and_add(&twisted_point, &scalar);
        assert_eq!(got2, got2);

        // Lets see if this is conserved over the isogenies
        let edwards_point = twisted_point.to_untwisted();
        let got_untwisted_point = edwards_point.scalar_mul(&scalar);
        let expected_untwisted_point = expected_twisted_point.to_untwisted();
        assert_eq!(got_untwisted_point, expected_untwisted_point);
    }

    #[test]
    fn test_simple_scalar_mul_identities() {
        let x = ExtendedPoint {
            X: FieldElement::from_raw_slice([
                0x034365c8, 0x06b2a874, 0x0eb875d7, 0x0ae4c7a7, 0x0785df04, 0x09929351, 0x01fe8c3b,
                0x0f2a0e5f, 0x0111d39c, 0x07ab52ba, 0x01df4552, 0x01d87566, 0x0f297be2, 0x027c090f,
                0x0a81b155, 0x0d1a562b,
            ]),
            Y: FieldElement::from_raw_slice([
                0x00da9708, 0x0e7d583e, 0x0dbcc099, 0x0d2dad89, 0x05a49ce4, 0x01cb4ddc, 0x0928d395,
                0x0098d91d, 0x0bff16ce, 0x06f02f9a, 0x0ce27cc1, 0x0dab5783, 0x0b553d94, 0x03251a0c,
                0x064d70fb, 0x07fe3a2f,
            ]),
            Z: FieldElement::from_raw_slice([
                0x0d5237cc, 0x0319d105, 0x02ab2df5, 0x022e9736, 0x0d79742f, 0x00688712, 0x012d3a65,
                0x0ef4925e, 0x0bd0d260, 0x0832b532, 0x05faef27, 0x01ffe567, 0x0161ce73, 0x07bda0f5,
                0x035d04f1, 0x0930f532,
            ]),
            T: FieldElement::from_raw_slice([
                0x01f6cc27, 0x09be7b8a, 0x0226da79, 0x0f6202f1, 0x0e7264dc, 0x0d25aeb1, 0x06c81f07,
                0x03c32cdc, 0x0923c854, 0x0cfc9865, 0x055b2fed, 0x05bdcc90, 0x01a99835, 0x0ea08056,
                0x0abbf763, 0x03826c2f,
            ]),
        };

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
