use std::ops::{Add, Mul, Neg, Sub};

use crate::{Scalar, curve::scalar_mul::double_and_add};

use super::DecafPoint;

/// Scalar Mul Operations
impl<'s, 'p> Mul<&'s Scalar> for &'p DecafPoint {
    type Output = DecafPoint;
    fn mul(self, scalar: &'s Scalar) -> DecafPoint {
        // XXX: We can do better than double and add
        DecafPoint(double_and_add(&self.0, &scalar))
    }
}
impl<'p, 's> Mul<&'p DecafPoint> for &'s Scalar {
    type Output = DecafPoint;
    fn mul(self, point: &'p DecafPoint) -> DecafPoint {
        DecafPoint(double_and_add(&point.0,self))
    }
}
impl Mul<DecafPoint> for Scalar {
    type Output = DecafPoint;
    fn mul(self, point: DecafPoint) -> DecafPoint {
        DecafPoint(double_and_add(&point.0,&self))
    }
}
impl Mul<Scalar> for DecafPoint {
    type Output = DecafPoint;
    fn mul(self, scalar : Scalar) -> DecafPoint {
        DecafPoint(double_and_add(&self.0, &scalar))
    }
}

// Point addition 

impl<'a, 'b> Add<&'a DecafPoint> for &'b DecafPoint {
    type Output = DecafPoint;
    fn add(self, other: &'a DecafPoint) -> DecafPoint {
        DecafPoint(self.0.to_extensible().add_extended(&other.0).to_extended())
    }
}
impl Add<DecafPoint> for DecafPoint {
    type Output = DecafPoint;
    fn add(self, other: DecafPoint) -> DecafPoint {
        (&self).add(&other)
    }
}

// Point Subtraction 

impl<'a, 'b> Sub<&'a DecafPoint> for &'b DecafPoint {
    type Output = DecafPoint;
    fn sub(self, other: &'a DecafPoint) -> DecafPoint {
        DecafPoint(self.0.to_extensible().sub_extended(&other.0).to_extended())
    }
}
impl Sub<DecafPoint> for DecafPoint {
    type Output = DecafPoint;
    fn sub(self, other: DecafPoint) -> DecafPoint {
        (&self).sub(&other)
    }
}

// Point Negation 

impl<'b> Neg for &'b DecafPoint {
    type Output = DecafPoint;
    fn neg(self) -> DecafPoint {
        DecafPoint(self.0.negate())
    }
}
impl Neg for DecafPoint {
    type Output = DecafPoint;
    fn neg(self) -> DecafPoint {
        (&self).neg()
    }
}

