use crate::field::base::Fq;

// In affine (x,y) is the extensible point (X, Y, Z, T1, T2)
// Where x = X/Z , y = Y/Z , T1 * T2 = T
// XXX: I think we have too many point representations,
// But let's not remove any yet
pub struct ExtensiblePoint {
    pub(crate) X: Fq,
    pub(crate) Y: Fq,
    pub(crate) Z: Fq,
    pub(crate) T1: Fq,
    pub(crate) T2: Fq,
}

impl ExtensiblePoint {
    pub fn identity() -> ExtensiblePoint {
        ExtensiblePoint {
            X: Fq::zero(),
            Y: Fq::one(),
            Z: Fq::one(),
            T1: Fq::zero(),
            T2: Fq::zero(),
        }
    }

    pub fn equals(&self, other: &ExtensiblePoint) -> bool {
        let ZX = self.Z * other.X;
        let XZ = self.X * other.Z;

        let ZY = self.Z * other.Y;
        let YZ = self.Y * other.Z;

        ZX.equals(&XZ) && ZY.equals(&YZ)
    }

    pub fn double(&self) -> ExtensiblePoint {
        let mut X = self.X;
        let mut Y = self.Y;
        let mut Z = self.Z;
        let mut T1 = self.T1;
        let mut T2 = self.T2;

        let L2 = X.square();
        let mut L0 = Y.square();
        T2 = L2.add_no_reduce(&L0);
        T1 = Y.add_no_reduce(&X);
        let mut L1 = T1.square();
        T1 = L1.sub_no_reduce(&T2);
        T1.bias(3);
        T1.weak_reduce();
        L1 = L0 - L2;
        X = Z.square();
        X.bias(1);
        Z = X.add_no_reduce(&X);
        L0 = Z.sub_no_reduce(&L1);
        L0.weak_reduce();
        Z = L0 * L1;
        X = L0 * T1;
        Y = L1 * T2;

        ExtensiblePoint { X, Y, Z, T1, T2 }
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
