use crate::curve::edwards::ExtendedPoint;
use crate::field::base::Fq;
use subtle::{Choice, ConditionallyNegatable, ConditionallySelectable};

impl Default for ProjectiveNielsPoint {
    fn default() -> ProjectiveNielsPoint {
        ProjectiveNielsPoint::identity()
    }
}

// Its a variant of Niels, where a Z coordinate is added for unmixed readdition
// ((y+x)/2, (y-x)/2, dxy, Z)
#[derive(Copy, Clone)]
pub struct ProjectiveNielsPoint {
    pub(crate) Y_plus_X: Fq,
    pub(crate) Y_minus_X: Fq,
    pub(crate) Td: Fq,
    pub(crate) Z: Fq,
}

impl ProjectiveNielsPoint {
    pub fn identity() -> ProjectiveNielsPoint {
        ProjectiveNielsPoint {
            Y_plus_X: Fq::one(),
            Y_minus_X: Fq::one(),
            Td: Fq::zero(),
            Z: Fq::one(),
        }
    }
    pub fn to_extended(&self) -> ExtendedPoint {
        let two_y = self.Y_plus_X + self.Y_minus_X;
        let two_x = self.Y_plus_X - self.Y_minus_X;

        assert_ne!(self.Z, Fq::zero());

        let T = two_y * two_x;
        let X = self.Z * two_x;
        let Y = self.Z * two_y;
        let Z = self.Z.square();

        ExtendedPoint { X, Y, Z, T }
    }

    // XXX: Code duplication here as AffineNielsPoint also has this.
    // We can switch it out once ScalarMul get refactored
    pub fn conditional_negate(&mut self, choice: Choice) {
        Fq::conditional_swap(&mut self.Y_minus_X, &mut self.Y_plus_X, choice);
        self.Td.conditional_negate(choice);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use hex::decode as hex_decode;

    fn slice_to_fixed_array(b: &[u8]) -> [u8; 56] {
        let mut a: [u8; 56] = [0; 56];
        a.copy_from_slice(&b);
        a
    }

    fn hex_to_fq(data: &str) -> Fq {
        let mut bytes = hex_decode(data).unwrap();
        bytes.reverse();
        Fq::from_bytes(&slice_to_fixed_array(&bytes))
    }

    #[test]
    fn test_conditional_negate() {
        let Y_minus_X = hex_to_fq("4b8a632c1feab72769cd96e7aaa577861871b3613945c802b89377e8b85331ecc0ffb1cb20169bfc9c27274d38b0d01e87a1d5d851770bc8");
        let Y_plus_X = hex_to_fq("81a45f02f41053f8d7d2a1f176a340529b33b7ee4d3fa84de384b750b35a54c315bf36c41d023ade226449916e668396589ea2145da09b95");
        let Td = hex_to_fq("5f5a2b06a2dbf7136f8dc979fd54d631ca7de50397250a196d3be2a721ab7cbaa92c545d9b15b5319e11b64bc031666049d8637e13838b3b");
        let Z = Fq::one();

        let mut n = ProjectiveNielsPoint {
            Y_plus_X,
            Y_minus_X,
            Td,
            Z,
        };

        let expected_neg_n = ProjectiveNielsPoint {
            Y_plus_X: Y_minus_X,
            Y_minus_X: Y_plus_X,
            Td: Td.negate(),
            Z: Z,
        };

        n.conditional_negate(1.into());

        assert!(expected_neg_n.Y_plus_X == n.Y_plus_X);
        assert!(expected_neg_n.Y_minus_X == n.Y_minus_X);
        assert!(expected_neg_n.Td == n.Td);
    }
}
