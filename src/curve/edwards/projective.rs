use crate::curve::constants::EDWARDS_D;
use crate::field::base::Fq;
pub struct ProjectivePoint {
    pub(crate) X: Fq,
    pub(crate) Y: Fq,
    pub(crate) Z: Fq,
}

impl ProjectivePoint {
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

pub struct ProjectiveNielsPoint {
    pub(crate) Y_plus_X: Fq,
    pub(crate) Y_minus_X: Fq,
    pub(crate) Z: Fq,
    pub(crate) T2d: Fq,
}
