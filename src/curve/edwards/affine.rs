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
}
