use crate::curve::constants::{TWO_D_MINUS_ONE, TWO_ONE_MINUS_D};
use crate::curve::twedwards::{
    affine::AffineNielsPoint, extended::ExtendedPoint, projective::ProjectiveNielsPoint,
};
use crate::field::FieldElement;
use subtle::{Choice, ConstantTimeEq};

// In affine (x,y) is the extensible point (X, Y, Z, T1, T2)
// Where x = X/Z , y = Y/Z , T1 * T2 = T
// XXX: I think we have too many point representations,
// But let's not remove any yet
pub struct ExtensiblePoint {
    pub(crate) X: FieldElement,
    pub(crate) Y: FieldElement,
    pub(crate) Z: FieldElement,
    pub(crate) T1: FieldElement,
    pub(crate) T2: FieldElement,
}

impl ConstantTimeEq for ExtensiblePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        let XZ = self.X * other.Z;
        let ZX = self.Z * other.X;

        let YZ = self.Y * other.Z;
        let ZY = self.Z * other.Y;

        (XZ.ct_eq(&ZX)) & (YZ.ct_eq(&ZY))
    }
}
impl PartialEq for ExtensiblePoint {
    fn eq(&self, other: &ExtensiblePoint) -> bool {
        self.ct_eq(other).into()
    }
}
impl Eq for ExtensiblePoint {}

impl ExtensiblePoint {
    pub fn identity() -> ExtensiblePoint {
        ExtensiblePoint {
            X: FieldElement::zero(),
            Y: FieldElement::one(),
            Z: FieldElement::one(),
            T1: FieldElement::zero(),
            T2: FieldElement::zero(),
        }
    }

    pub fn equals(&self, other: &ExtensiblePoint) -> bool {
        let ZX = self.Z * other.X;
        let XZ = self.X * other.Z;

        let ZY = self.Z * other.Y;
        let YZ = self.Y * other.Z;

        (ZX == XZ) && (ZY == YZ)
    }

    pub fn double(&self) -> ExtensiblePoint {
        let XX = self.X.square();
        let YY = self.Y.square();

        let XX_plus_YY = XX.add_no_reduce(&YY);
        let YY_minus_XX = YY - XX;
        let Y_plus_X = self.Y.add_no_reduce(&self.X);
        let Y_plus_X2 = Y_plus_X.square();

        let T1 = Y_plus_X2 - (XX_plus_YY);

        let ZZ = self.Z.square();

        let ZZ_plus_ZZ = ZZ + (ZZ);

        let ZZ_YY_XX = ZZ_plus_ZZ - (YY_minus_XX);

        let Z = ZZ_YY_XX * YY_minus_XX;
        let X = ZZ_YY_XX * T1;
        let Y = YY_minus_XX * XX_plus_YY;

        ExtensiblePoint {
            X,
            Y,
            Z,
            T1,
            T2: XX_plus_YY,
        }
    }

    pub fn add_extensible(&self, other: &ExtensiblePoint) -> ExtensiblePoint {
        self.add_extended(&other.to_extended())
    }
    pub fn add_extended(&self, other: &ExtendedPoint) -> ExtensiblePoint {
        // Compute T
        let self_T = self.T1 * self.T2;
        let other_T = other.T;

        let mut result = ExtensiblePoint::identity();

        let (mut a, mut b, mut c, mut d) = (
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
        );

        b = self.Y - self.X;
        c = other.Y - other.X;
        d = other.Y.add_no_reduce(&other.X);
        a = c * b;
        b = self.Y.add_no_reduce(&self.X);
        result.Y = d * b;
        b = other_T * self_T;
        result.X = b * TWO_ONE_MINUS_D;
        b = a.add_no_reduce(&result.Y);
        c = result.Y - a;
        a = self.Z * other.Z;
        a = a.add_no_reduce(&a);
        result.Y = a.add_no_reduce(&result.X);
        a = a - result.X;
        result.Z = a * result.Y;
        result.X = result.Y * c;
        result.Y = a * b;
        result.T1 = b;
        result.T2 = c;

        result
    }
    pub fn sub_extended(&self, other: &ExtendedPoint) -> ExtensiblePoint {
        // Compute T
        let self_T = self.T1 * self.T2;
        let other_T = other.T;

        let mut result = ExtensiblePoint::identity();

        let (mut a, mut b, mut c, mut d) = (
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
        );

        b = self.Y - self.X;
        d = other.Y - other.X;
        c = other.Y.add_no_reduce(&other.X);
        a = c * b;
        b = self.Y.add_no_reduce(&self.X);
        result.Y = d * b;
        b = other_T * self_T;
        result.X = b * TWO_ONE_MINUS_D;
        b = a.add_no_reduce(&result.Y);
        c = result.Y - a;
        a = self.Z * other.Z;
        a = a.add_no_reduce(&a);
        result.Y = a - result.X;
        a = a.add_no_reduce(&result.X);
        result.Z = a * result.Y;
        result.X = result.Y * c;
        result.Y = a * b;
        result.T1 = b;
        result.T1 = c;

        result
    }

    pub fn add_affine_niels(&self, other: AffineNielsPoint) -> ExtensiblePoint {
        let mut a = FieldElement::zero();
        let mut b = FieldElement::zero();
        let mut c = FieldElement::zero();

        let mut X = self.X;
        let mut Y = self.Y;
        let mut Z = self.Z;
        let mut T1 = self.T1;
        let mut T2 = self.T2;

        b = Y - X;
        a = other.y_minus_x * b;
        b = X.add_no_reduce(&Y);
        Y = other.y_plus_x * b;
        X = other.td * T1 * T2;
        c = a.add_no_reduce(&Y);
        b = Y - a;
        Y = Z - X;
        a = X.add_no_reduce(&Z);
        Z = a * Y;
        X = Y * b;
        Y = a * c;

        T1 = b;
        T2 = c;

        ExtensiblePoint { X, Y, Z, T1, T2 }
    }

    pub fn add_projective_niels(&mut self, other: &ProjectiveNielsPoint) -> ExtensiblePoint {
        self.Z = self.Z * other.Z;

        self.add_affine_niels(AffineNielsPoint {
            y_plus_x: other.Y_plus_X,
            y_minus_x: other.Y_minus_X,
            td: other.Td,
        })
    }
    pub fn to_extended(&self) -> ExtendedPoint {
        ExtendedPoint {
            X: self.X,
            Y: self.Y,
            Z: self.Z,
            T: self.T1 * self.T2,
        }
    }

    pub fn to_projective_niels(&self) -> ProjectiveNielsPoint {
        ProjectiveNielsPoint {
            Y_plus_X: self.X + self.Y,
            Y_minus_X: self.Y - self.X,
            Z: self.Z + self.Z,
            Td: self.T1 * self.T2 * TWO_D_MINUS_ONE,
        }
    }
}
#[cfg(test)]
mod tests {

    use super::*;
    #[test]
    fn test_basic_double() {
        let iden = ExtensiblePoint::identity();
        assert!(iden.equals(&iden.double()));
    }
}
