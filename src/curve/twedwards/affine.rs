use crate::curve::constants::EDWARDS_D;
use crate::curve::twedwards::{extended::ExtendedPoint, extensible::ExtensiblePoint};
use crate::field::FieldElement;
use subtle::{Choice, ConditionallySelectable};
#[derive(PartialEq, Eq)]
pub struct AffinePoint {
    pub(crate) x: FieldElement,
    pub(crate) y: FieldElement,
}

impl Default for AffinePoint {
    fn default() -> AffinePoint {
        AffinePoint::identity()
    }
}

impl AffinePoint {
    pub(crate) fn identity() -> AffinePoint {
        AffinePoint {
            x: FieldElement::zero(),
            y: FieldElement::one(),
        }
    }
    fn is_on_curve(&self) -> bool {
        let xx = self.x.square();
        let yy = self.y.square();

        xx + yy == FieldElement::one() + (EDWARDS_D * xx * yy)
    }
    pub fn negate(&self) -> AffinePoint {
        AffinePoint {
            x: self.x.negate(),
            y: self.y,
        }
    }
    fn add(&self, other: &AffinePoint) -> AffinePoint {
        let y_numerator = self.y * other.y - self.x * other.x;
        let y_denominator = FieldElement::one() - EDWARDS_D * self.x * other.x * self.y * other.y;

        let x_numerator = self.x * other.y + self.y * other.x;
        let x_denominator = FieldElement::one() + EDWARDS_D * self.x * other.x * self.y * other.y;

        let x = x_numerator * x_denominator.invert();
        let y = y_numerator * y_denominator.invert();
        AffinePoint { x, y }
    }
    pub fn to_extensible(&self) -> ExtensiblePoint {
        ExtensiblePoint {
            X: self.x,
            Y: self.y,
            Z: FieldElement::one(),
            T1: self.x,
            T2: self.y,
        }
    }
    pub fn to_affine_niels(&self) -> AffineNielsPoint {
        AffineNielsPoint {
            y_plus_x: self.y + self.x,
            y_minus_x: self.y - self.x,
            td: self.x * self.y * EDWARDS_D,
        }
    }

    pub fn to_extended(&self) -> ExtendedPoint {
        ExtendedPoint {
            X: self.x,
            Y: self.y,
            Z: FieldElement::one(),
            T: self.x * self.y,
        }
    }
}
// ((y+x)/2, (y-x)/2, dxy)
#[derive(Copy, Clone)]
pub struct AffineNielsPoint {
    pub(crate) y_plus_x: FieldElement,
    pub(crate) y_minus_x: FieldElement,
    pub(crate) td: FieldElement,
}

impl ConditionallySelectable for AffineNielsPoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        AffineNielsPoint {
            y_plus_x: FieldElement::conditional_select(&a.y_plus_x, &b.y_plus_x, choice),
            y_minus_x: FieldElement::conditional_select(&a.y_minus_x, &b.y_minus_x, choice),
            td: FieldElement::conditional_select(&a.td, &b.td, choice),
        }
    }
}

impl AffineNielsPoint {
    pub fn equals(&self, other: &AffineNielsPoint) -> bool {
        (self.y_minus_x == other.y_minus_x)
            && (self.y_plus_x == other.y_plus_x)
            && (self.td == other.td)
    }

    pub fn identity() -> AffineNielsPoint {
        AffineNielsPoint {
            y_plus_x: FieldElement::one(),
            y_minus_x: FieldElement::one(),
            td: FieldElement::zero(),
        }
    }

    pub fn to_extended(&self) -> ExtendedPoint {
        ExtendedPoint {
            X: self.y_plus_x - self.y_minus_x,
            Y: self.y_minus_x + self.y_plus_x,
            Z: FieldElement::one(),
            T: self.y_plus_x * self.y_minus_x,
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use hex::decode as hex_decode;

    #[test]
    fn test_to_extensible() {
        let a = AffinePoint {
            x: FieldElement::from_raw_slice([
                0x034365c8, 0x06b2a874, 0x0eb875d7, 0x0ae4c7a7, 0x0785df04, 0x09929351, 0x01fe8c3b,
                0x0f2a0e5f, 0x0111d39c, 0x07ab52ba, 0x01df4552, 0x01d87566, 0x0f297be2, 0x027c090f,
                0x0a81b155, 0x0d1a562b,
            ]),
            y: FieldElement::from_raw_slice([
                0x00da9708, 0x0e7d583e, 0x0dbcc099, 0x0d2dad89, 0x05a49ce4, 0x01cb4ddc, 0x0928d395,
                0x0098d91d, 0x0bff16ce, 0x06f02f9a, 0x0ce27cc1, 0x0dab5783, 0x0b553d94, 0x03251a0c,
                0x064d70fb, 0x07fe3a2f,
            ]),
        };

        let expected_b = ExtensiblePoint {
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
            Z: FieldElement::one(),
            T1: FieldElement::from_raw_slice([
                0x034365c8, 0x06b2a874, 0x0eb875d7, 0x0ae4c7a7, 0x0785df04, 0x09929351, 0x01fe8c3b,
                0x0f2a0e5f, 0x0111d39c, 0x07ab52ba, 0x01df4552, 0x01d87566, 0x0f297be2, 0x027c090f,
                0x0a81b155, 0x0d1a562b,
            ]),
            T2: FieldElement::from_raw_slice([
                0x00da9708, 0x0e7d583e, 0x0dbcc099, 0x0d2dad89, 0x05a49ce4, 0x01cb4ddc, 0x0928d395,
                0x0098d91d, 0x0bff16ce, 0x06f02f9a, 0x0ce27cc1, 0x0dab5783, 0x0b553d94, 0x03251a0c,
                0x064d70fb, 0x07fe3a2f,
            ]),
        };

        let b = a.to_extensible();
        assert!(b.equals(&expected_b));
    }

    fn slice_to_fixed_array(b: &[u8]) -> [u8; 56] {
        let mut a: [u8; 56] = [0; 56];
        a.copy_from_slice(&b);
        a
    }

    fn hex_to_field(data: &str) -> FieldElement {
        let mut bytes = hex_decode(data).unwrap();
        bytes.reverse();
        FieldElement::from_bytes(&slice_to_fixed_array(&bytes))
    }
    #[test]
    fn test_negation() {
        // Taken from paper
        let x  = hex_to_field("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa955555555555555555555555555555555555555555555555555555555");
        let y  = hex_to_field("ae05e9634ad7048db359d6205086c2b0036ed7a035884dd7b7e36d728ad8c4b80d6565833a2a3098bbbcb2bed1cda06bdaeafbcdea9386ed");
        let a = AffinePoint { x, y };
        assert!(a.is_on_curve());

        let neg_a = a.negate();

        let got = neg_a.add(&a);

        assert!(got == AffinePoint::identity());
    }
}
