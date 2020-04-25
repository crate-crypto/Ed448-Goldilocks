use crate::curve::edwards::ExtendedPoint;
use crate::field::Fq;
// Affine point on untwisted curve
// XXX: This is only really needed for convenience in extended.rs . Will remove it sooner or later
pub struct AffinePoint {
    pub(crate) x: Fq,
    pub(crate) y: Fq,
}

impl AffinePoint {
    pub fn identity() -> AffinePoint {
        AffinePoint {
            x: Fq::zero(),
            y: Fq::one(),
        }
    }
    pub fn to_extended(&self) -> ExtendedPoint {
        ExtendedPoint {
            X: self.x,
            Y: self.y,
            Z: Fq::one(),
            T: self.x * self.y,
        }
    }
}
