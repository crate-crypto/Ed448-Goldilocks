use crate::curve::constants::EDWARDS_D;
use crate::field::base::Fq;
pub struct ProjectivePoint {
    pub(crate) X: Fq,
    pub(crate) Y: Fq,
    pub(crate) Z: Fq,
}

impl ProjectivePoint {
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
}
