use crate::curve::constants::EDWARDS_D;
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
            Z: Fq::one()
        }
    }
    pub(crate) fn is_on_curve(&self) -> bool {
        // x^2 + y^2 = 1 + d*x^2*y^2
        // Homogenised: (X^2 + Y^2)Z^4 = Z^4 + dX^2Y^2
        let XX = self.X.square();
        let YY = self.Y.square();
        let ZZ = self.Z.square();
        let ZZZZ = ZZ.square();

        let lhs = (XX + YY) * ZZZZ;
        let rhs = ZZZZ + (EDWARDS_D * XX * YY);
        lhs.equals(&rhs)
    }
}

impl Default for ProjectiveNielsPoint {
    fn default() -> ProjectiveNielsPoint {
        ProjectiveNielsPoint::identity()
    }
}

pub struct ProjectiveNielsPoint {
    pub(crate) Y_plus_X: Fq,
    pub(crate) Y_minus_X: Fq,
    pub(crate) Z: Fq,
    pub(crate) T2d: Fq,
}

impl ProjectiveNielsPoint {
    pub fn identity() -> ProjectiveNielsPoint {
        ProjectiveNielsPoint {
            Y_plus_X: Fq::one(),
            Y_minus_X: Fq::one(),
            Z: Fq::one(),
            T2d: Fq::zero()
        }
    }
}
