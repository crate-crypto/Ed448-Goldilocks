use crate::curve::constants::EDWARDS_D;
use crate::curve::edwards::affine::AffinePoint;
use crate::curve::montgomery::montgomery::MontgomeryPoint; // XXX: need to fix this path
use crate::curve::twedwards::extended::ExtendedPoint as TwistedExtendedPoint;
use crate::field::FieldElement;
use crate::field::Scalar;
use subtle::{Choice, ConditionallyNegatable, ConstantTimeEq};
/// Represent points on the (untwisted) edwards curve using Extended Homogenous Projective Co-ordinates
/// (x, y) -> (X/Z, Y/Z, Z, T)
/// a = 1, d = -39081
/// XXX: Make this more descriptive
/// XXX: WE could probably add the basic double and add methods here for untwisted points, and only use the isogeny when we need to do scalar mul, but we would then have two formulas for doubling etc
/// XXX: THen again we will most likely need doubling to check the isogeny, so maybe add the double and add methods here because if we want to add a+b, then we dont want to force the isogeny which would cost more
/// Should this be renamed to EdwardsPoint so that we are consistent with Dalek crypto?
#[derive(Copy, Clone, Debug)]
pub struct ExtendedPoint {
    pub(crate) X: FieldElement,
    pub(crate) Y: FieldElement,
    pub(crate) Z: FieldElement,
    pub(crate) T: FieldElement,
}

pub struct CompressedEdwardsY(pub [u8; 57]);

impl CompressedEdwardsY {
    pub fn decompress(&self) -> Option<ExtendedPoint> {
        if let Some((sign, b)) = self.0.split_last() {
            let sign = *sign >> 7;

            let mut y_bytes: [u8; 56] = [0; 56];
            y_bytes.copy_from_slice(&b);

            // Recover x
            let y = FieldElement::from_bytes(&y_bytes);
            let yy = y.square();
            let dyy = EDWARDS_D * yy;
            let numerator = FieldElement::one() - yy;
            let denominator = FieldElement::one() - dyy;

            let (mut x, is_res) = FieldElement::sqrt_ratio(&numerator, &denominator);
            if !is_res {
                return None;
            }
            x.strong_reduce();

            // Compute correct sign of x
            let compressed_sign_bit = Choice::from(sign >> 7);
            x.conditional_negate(compressed_sign_bit);

            return Some(ExtendedPoint {
                X: x,
                Y: y,
                Z: FieldElement::one(),
                T: x * y,
            });
        } else {
            return None;
        };
    }
}

impl ConstantTimeEq for ExtendedPoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        let XZ = self.X * other.Z;
        let ZX = self.Z * other.X;

        let YZ = self.Y * other.Z;
        let ZY = self.Z * other.Y;

        (XZ.ct_eq(&ZX)) & (YZ.ct_eq(&ZY))
    }
}

impl PartialEq for ExtendedPoint {
    fn eq(&self, other: &ExtendedPoint) -> bool {
        self.ct_eq(other).into()
    }
}
impl Eq for ExtendedPoint {}

impl Default for ExtendedPoint {
    fn default() -> ExtendedPoint {
        ExtendedPoint::identity()
    }
}

impl ExtendedPoint {
    /// Identity point
    pub fn identity() -> ExtendedPoint {
        ExtendedPoint {
            X: FieldElement::zero(),
            Y: FieldElement::one(),
            Z: FieldElement::one(),
            T: FieldElement::zero(),
        }
    }

    pub fn to_montgomery(&self) -> MontgomeryPoint {
        // u = Y^2 / X^2

        let affine = self.to_affine();

        let xx = affine.x.square();
        let yy = affine.y.square();

        let u = yy * xx.invert();

        MontgomeryPoint(u.to_bytes())
    }

    // Since we do not know whether the point p is in the prime subgroup, we need to mul by floor(s/4)
    // XXX: We also need to add (s mod 4) * P . This strategy does not involve Ristretto/Decaf
    pub fn scalar_mul(&self, scalar: &Scalar) -> ExtendedPoint {
        let adjusted_scalar = *scalar * Scalar::from(4).invert();
        let twisted_point = self.to_twisted();
        use crate::curve::scalar_mul;
        let partial_result =
            scalar_mul::scalar_mul(&twisted_point, &adjusted_scalar).to_untwisted();

        // // Compute s mod 4
        // let s_mod_four = Scalar::from(scalar[0] & 3);

        partial_result
    }

    // Standard compression; store Y and sign of X
    // XXX: This needs more docs and is `compress` the conventional function name? I think to_bytes/encode is?
    pub fn compress(&self) -> CompressedEdwardsY {
        let affine = self.to_affine();

        let affine_x = affine.x;
        let affine_y = affine.y;

        let mut compressed_bytes = [0u8; 57];

        let sign = affine_x.is_negative().unwrap_u8() << 7;

        let y_bytes = affine_y.to_bytes();
        for i in 0..y_bytes.len() {
            compressed_bytes[i] = y_bytes[i]
        }
        *compressed_bytes.last_mut().unwrap() = sign as u8;
        CompressedEdwardsY(compressed_bytes)
    }

    //https://iacr.org/archive/asiacrypt2008/53500329/53500329.pdf (3.1)
    // These formulas are unified, so for now we can use it for doubling. Will refactor later for speed
    pub fn add(&self, other: &ExtendedPoint) -> ExtendedPoint {
        let aXX = self.X * other.X; // aX1X2
        let dTT = EDWARDS_D * self.T * other.T; // dT1T2
        let ZZ = self.Z * other.Z; // Z1Z2
        let YY = self.Y * other.Y;

        let X = {
            let x_1 = (self.X * other.Y) + (self.Y * other.X);
            let x_2 = ZZ - dTT;
            x_1 * x_2
        };
        let Y = {
            let y_1 = YY - aXX;
            let y_2 = ZZ + dTT;
            y_1 * y_2
        };

        let T = {
            let t_1 = YY - aXX;
            let t_2 = (self.X * other.Y) + (self.Y * other.X);
            t_1 * t_2
        };

        let Z = { (ZZ - dTT) * (ZZ + dTT) };

        ExtendedPoint { X, Y, Z, T }
    }

    // XXX: See comment on addition, the formula is unified, so this will do for now
    //https://iacr.org/archive/asiacrypt2008/53500329/53500329.pdf (3.1)
    pub fn double(&self) -> ExtendedPoint {
        self.add(&self)
    }

    pub(crate) fn is_on_curve(&self) -> bool {
        let XY = self.X * self.Y;
        let ZT = self.Z * self.T;

        // Y^2 + X^2 == Z^2 - T^2 * D

        let YY = self.Y.square();
        let XX = self.X.square();
        let ZZ = self.Z.square();
        let TT = self.T.square();
        let lhs = YY + XX;
        let rhs = ZZ + TT * EDWARDS_D;

        (XY == ZT) && (lhs == rhs)
    }

    pub fn to_affine(&self) -> AffinePoint {
        let INV_Z = self.Z.invert();

        let mut x = self.X * INV_Z;
        x.strong_reduce();

        let mut y = self.Y * INV_Z;
        y.strong_reduce();

        AffinePoint { x, y }
    }
    /// Edwards_Isogeny is derived from the doubling formula
    /// XXX: There is a duplicate method in the twisted edwards module to compute the dual isogeny
    /// XXX: Not much point trying to make it generic I think. So what we can do is optimise each respective isogeny method for a=1 or a = -1 (currently, I just made it really slow and simple)
    fn edwards_isogeny(&self, a: FieldElement) -> TwistedExtendedPoint {
        // Convert to affine now, then derive extended version later
        let affine = self.to_affine();
        let x = affine.x;
        let y = affine.y;

        // Compute x
        let xy = x * y;
        let x_numerator = xy + xy;
        let x_denom = y.square() - (a * x.square());
        let new_x = x_numerator * x_denom.invert();

        // Compute y
        let y_numerator = y.square() + (a * x.square());
        let y_denom = (FieldElement::one() + FieldElement::one()) - y.square() - (a * x.square());
        let new_y = y_numerator * y_denom.invert();

        TwistedExtendedPoint {
            X: new_x,
            Y: new_y,
            Z: FieldElement::one(),
            T: new_x * new_y,
        }
    }
    pub fn to_twisted(&self) -> TwistedExtendedPoint {
        self.edwards_isogeny(FieldElement::one())
    }

    pub fn negate(&self) -> ExtendedPoint {
        ExtendedPoint {
            X: self.X.negate(),
            Y: self.Y,
            Z: self.Z,
            T: self.T.negate(),
        }
    }

    pub fn torque(&self) -> ExtendedPoint {
        ExtendedPoint {
            X: self.X.negate(),
            Y: self.Y.negate(),
            Z: self.Z,
            T: self.T,
        }
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

    fn hex_to_field(data: &str) -> FieldElement {
        let mut bytes = hex_decode(data).unwrap();
        bytes.reverse();
        FieldElement::from_bytes(&slice_to_fixed_array(&bytes))
    }
    #[test]
    fn test_isogeny() {
        let x  = hex_to_field("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa955555555555555555555555555555555555555555555555555555555");
        let y  = hex_to_field("ae05e9634ad7048db359d6205086c2b0036ed7a035884dd7b7e36d728ad8c4b80d6565833a2a3098bbbcb2bed1cda06bdaeafbcdea9386ed");
        let a = AffinePoint { x, y }.to_extended();
        let twist_a = a.to_twisted().to_untwisted();
        assert!(twist_a == a.double().double())
    }
    #[test]
    fn test_is_on_curve() {
        let x  = hex_to_field("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa955555555555555555555555555555555555555555555555555555555");
        let y  = hex_to_field("ae05e9634ad7048db359d6205086c2b0036ed7a035884dd7b7e36d728ad8c4b80d6565833a2a3098bbbcb2bed1cda06bdaeafbcdea9386ed");
        let gen = AffinePoint { x, y }.to_extended();
        assert!(gen.is_on_curve());
    }
    #[test]
    fn test_compress_decompress() {
        let x  = hex_to_field("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa955555555555555555555555555555555555555555555555555555555");
        let y  = hex_to_field("ae05e9634ad7048db359d6205086c2b0036ed7a035884dd7b7e36d728ad8c4b80d6565833a2a3098bbbcb2bed1cda06bdaeafbcdea9386ed");
        let gen = AffinePoint { x, y }.to_extended();

        let decompressed_point = gen.compress().decompress();
        assert!(decompressed_point.is_some());

        assert!(gen == decompressed_point.unwrap());

        // XXX: Add a point that should not be on the untwisted curve
    }
}
