use crate::curve::constants::EDWARDS_D;
use crate::curve::edwards::affine::AffineNielsPoint;
use crate::curve::edwards::ExtendedPoint;
use crate::field::base::Fq;
pub struct ProjectivePoint {
    pub(crate) X: Fq,
    pub(crate) Y: Fq,
    pub(crate) Z: Fq,
}

impl Default for ProjectivePoint {
    fn default() -> ProjectivePoint {
        ProjectivePoint::identity()
    }
}

impl ProjectivePoint {
    pub fn identity() -> ProjectivePoint {
        ProjectivePoint {
            X: Fq::zero(),
            Y: Fq::one(),
            Z: Fq::one(),
        }
    }
    pub(crate) fn is_on_curve(&self) -> bool {
        let XX = self.X.square();
        let YY = self.Y.square();
        let ZZ = self.Z.square();
        let ZZZZ = ZZ.square();

        let lhs = (XX + YY) * ZZ;
        let rhs = ZZZZ - (EDWARDS_D * XX * YY);

        lhs.equals(&rhs)
    }
}

impl Default for ProjectiveNielsPoint {
    fn default() -> ProjectiveNielsPoint {
        ProjectiveNielsPoint::identity()
    }
}

// Check this
// Its a variant of Niels, where a Z coordinate is added for unmixed readdition
// XXX: Name it better
// ((y+x)/2, (y-x)/2, dxy, Z)
#[derive(Copy, Clone)]
pub struct ProjectiveNielsPoint {
    pub(crate) y_plus_x: Fq,
    pub(crate) y_minus_x: Fq,
    pub(crate) td: Fq,
    pub(crate) Z: Fq,
}

impl ProjectiveNielsPoint {
    pub fn identity() -> ProjectiveNielsPoint {
        ProjectiveNielsPoint {
            y_plus_x: Fq::one(),
            y_minus_x: Fq::one(),
            td: Fq::zero(),
            Z: Fq::one(),
        }
    }
    pub fn to_extended(&self) -> ExtendedPoint {
        let two_y = self.y_plus_x + self.y_minus_x;
        let two_x = self.y_plus_x - self.y_minus_x;

        assert_ne!(self.Z, Fq::zero());

        let T = two_y * two_x;
        let X = self.Z * two_x;
        let Y = self.Z * two_y;
        let Z = self.Z.square();

        ExtendedPoint { X, Y, Z, T }
    }

    pub fn to_affine_niels(&self) -> AffineNielsPoint {
        AffineNielsPoint {
            y_plus_x: self.y_plus_x,
            y_minus_x: self.y_minus_x,
            td: self.td,
        }
    }
    // XXX: Code duplication here as AffineNielsPoint also has this.
    // We can switch it out once ScalarMul get refactored
    pub fn conditional_negate(&mut self, neg: u32) {
        self.y_minus_x.conditional_swap(&mut self.y_plus_x, neg);
        self.td.conditional_negate(neg);
    }
}
