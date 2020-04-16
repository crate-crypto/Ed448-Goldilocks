use crate::curve::edwards::{extended::ExtendedPoint, extensible::ExtensiblePoint};
use crate::field::base::Fq;
pub struct AffinePoint {
    pub(crate) y: Fq,
    pub(crate) x: Fq,
}

impl AffinePoint {
    pub fn to_extensible(&self) -> ExtensiblePoint {
        ExtensiblePoint {
            X: self.x,
            Y: self.y,
            Z: Fq::one(),
            T1: self.x,
            T2: self.y,
        }
    }
}
// ((y+x)/2, (y-x)/2, dxy)
pub struct AffineNielsPoint {
    pub(crate) y_plus_x: Fq,
    pub(crate) y_minus_x: Fq,
    pub(crate) td: Fq,
}

impl AffineNielsPoint {
    pub fn equals(&self, other: &AffineNielsPoint) -> bool {
        self.y_minus_x.equals(&other.y_minus_x)
            && self.y_plus_x.equals(&other.y_plus_x)
            && self.td.equals(&other.td)
    }

    pub fn to_extended(&self) -> ExtendedPoint {
        ExtendedPoint {
            X: self.y_plus_x - self.y_minus_x,
            Y: self.y_minus_x + self.y_plus_x,
            Z: Fq::one(),
            T: self.y_plus_x * self.y_minus_x,
        }
    }

    pub fn conditional_negate(&mut self, neg: u32) {
        self.y_minus_x.conditional_swap(&mut self.y_plus_x, neg);
        self.td.conditional_negate(neg);
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_to_extensible() {
        let a = AffinePoint {
            x: Fq([
                0x034365c8, 0x06b2a874, 0x0eb875d7, 0x0ae4c7a7, 0x0785df04, 0x09929351, 0x01fe8c3b,
                0x0f2a0e5f, 0x0111d39c, 0x07ab52ba, 0x01df4552, 0x01d87566, 0x0f297be2, 0x027c090f,
                0x0a81b155, 0x0d1a562b,
            ]),
            y: Fq([
                0x00da9708, 0x0e7d583e, 0x0dbcc099, 0x0d2dad89, 0x05a49ce4, 0x01cb4ddc, 0x0928d395,
                0x0098d91d, 0x0bff16ce, 0x06f02f9a, 0x0ce27cc1, 0x0dab5783, 0x0b553d94, 0x03251a0c,
                0x064d70fb, 0x07fe3a2f,
            ]),
        };

        let expected_b = ExtensiblePoint {
            X: Fq([
                0x034365c8, 0x06b2a874, 0x0eb875d7, 0x0ae4c7a7, 0x0785df04, 0x09929351, 0x01fe8c3b,
                0x0f2a0e5f, 0x0111d39c, 0x07ab52ba, 0x01df4552, 0x01d87566, 0x0f297be2, 0x027c090f,
                0x0a81b155, 0x0d1a562b,
            ]),
            Y: Fq([
                0x00da9708, 0x0e7d583e, 0x0dbcc099, 0x0d2dad89, 0x05a49ce4, 0x01cb4ddc, 0x0928d395,
                0x0098d91d, 0x0bff16ce, 0x06f02f9a, 0x0ce27cc1, 0x0dab5783, 0x0b553d94, 0x03251a0c,
                0x064d70fb, 0x07fe3a2f,
            ]),
            Z: Fq::one(),
            T1: Fq([
                0x034365c8, 0x06b2a874, 0x0eb875d7, 0x0ae4c7a7, 0x0785df04, 0x09929351, 0x01fe8c3b,
                0x0f2a0e5f, 0x0111d39c, 0x07ab52ba, 0x01df4552, 0x01d87566, 0x0f297be2, 0x027c090f,
                0x0a81b155, 0x0d1a562b,
            ]),
            T2: Fq([
                0x00da9708, 0x0e7d583e, 0x0dbcc099, 0x0d2dad89, 0x05a49ce4, 0x01cb4ddc, 0x0928d395,
                0x0098d91d, 0x0bff16ce, 0x06f02f9a, 0x0ce27cc1, 0x0dab5783, 0x0b553d94, 0x03251a0c,
                0x064d70fb, 0x07fe3a2f,
            ]),
        };

        let b = a.to_extensible();
        assert!(b.equals(&expected_b));
    }

    use std::num::ParseIntError;
    fn decode_hex(s: &str) -> Result<Vec<u8>, ParseIntError> {
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16))
            .collect()
    }

    fn slice_to_fixed_array(b: &[u8]) -> [u8; 56] {
        let mut a: [u8; 56] = [0; 56];
        a.copy_from_slice(&b);
        a
    }

    fn hex_to_fq(hex: &str) -> Fq {
        let mut bytes = decode_hex(hex).unwrap();
        bytes.reverse();
        Fq::from_bytes(&slice_to_fixed_array(&bytes))
    }
    #[test]
    fn test_conditional_negate() {
        let neg_one = 0xffffffff;

        let y_minus_x = hex_to_fq("4b8a632c1feab72769cd96e7aaa577861871b3613945c802b89377e8b85331ecc0ffb1cb20169bfc9c27274d38b0d01e87a1d5d851770bc8");
        let y_plus_x = hex_to_fq("81a45f02f41053f8d7d2a1f176a340529b33b7ee4d3fa84de384b750b35a54c315bf36c41d023ade226449916e668396589ea2145da09b95");
        let td = hex_to_fq("5f5a2b06a2dbf7136f8dc979fd54d631ca7de50397250a196d3be2a721ab7cbaa92c545d9b15b5319e11b64bc031666049d8637e13838b3b");

        let mut n = AffineNielsPoint {
            y_plus_x,
            y_minus_x,
            td,
        };

        let expected_neg_n = AffineNielsPoint {
            y_plus_x: y_minus_x,
            y_minus_x: y_plus_x,
            td: td.negate(),
        };

        n.conditional_negate(neg_one);

        assert!(expected_neg_n.y_plus_x.equals(&n.y_plus_x));
        assert!(expected_neg_n.y_minus_x.equals(&n.y_minus_x));
        assert!(expected_neg_n.td.equals(&n.td));
    }
}
