use crate::curve::constants::D_MINUS_ONE;
use crate::curve::edwards::ExtendedPoint as EdwardsExtendedPoint;
use crate::curve::twedwards::affine::AffinePoint;
use crate::curve::twedwards::extensible::ExtensiblePoint;
use crate::field::base::Fq;
use subtle::{Choice, ConstantTimeEq};

#[derive(Copy, Clone, Debug)]
//XXX: ExtendedProjectivePoint
pub struct ExtendedPoint {
    pub(crate) X: Fq,
    pub(crate) Y: Fq,
    pub(crate) Z: Fq,
    pub(crate) T: Fq,
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
            X: Fq::zero(),
            Y: Fq::one(),
            Z: Fq::one(),
            T: Fq::zero(),
        }
    }

    pub(crate) fn double(&self) -> ExtendedPoint {
        self.to_extensible().double().to_extended()
    }

    pub fn to_extensible(&self) -> ExtensiblePoint {
        ExtensiblePoint {
            X: self.X,
            Y: self.Y,
            Z: self.Z,
            T1: self.T,
            T2: Fq::one(),
        }
    }

    /// XXX: We initially took this from the go version, according to this function
    /// it was treating the Extended Coordinates as Extended Jacobian Coordinates.
    /// Would be interesting to check what the cost of Extended Jacobian is   
    pub(crate) fn to_affine(&self) -> AffinePoint {
        // Points to consider:
        // - All points where Z=0, translate to (0,0)
        // - The identity point has z=1, so it is not a problem

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
    fn edwards_isogeny(&self, a: Fq) -> EdwardsExtendedPoint {
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
        let y_denom = Fq::from(2) - y.square() - (a * x.square());
        let new_y = y_numerator * y_denom.invert();

        EdwardsExtendedPoint {
            X: new_x,
            Y: new_y,
            Z: Fq::one(),
            T: new_x * new_y,
        }
    }

    pub fn to_untwisted(&self) -> EdwardsExtendedPoint {
        self.edwards_isogeny(Fq::one().negate())
    }

    pub(crate) fn is_on_curve(&self) -> bool {
        let XY = self.X * self.Y;
        let ZT = self.Z * self.T;

        // Y^2 - X^2 == Z^2 + T^2 * (D - 1) // XXX: Put this in a constants file D-1 = TWISTED_D

        let YY = self.Y.square();
        let XX = self.X.square();
        let ZZ = self.Z.square();
        let TT = self.T.square();
        let lhs = YY - XX;
        let rhs = ZZ + TT * D_MINUS_ONE;

        (XY == ZT) && (lhs == rhs)
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

    fn hex_to_fq(data: &str) -> Fq {
        let mut bytes = hex_decode(data).unwrap();
        bytes.reverse();
        Fq::from_bytes(&slice_to_fixed_array(&bytes))
    }
    #[test]
    fn test_isogeny() {
        let x  = hex_to_fq("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa955555555555555555555555555555555555555555555555555555555");
        let y  = hex_to_fq("ae05e9634ad7048db359d6205086c2b0036ed7a035884dd7b7e36d728ad8c4b80d6565833a2a3098bbbcb2bed1cda06bdaeafbcdea9386ed");
        let a = AffinePoint { x, y }.to_extended();
        let twist_a = a.to_untwisted().to_twisted();
        assert!(twist_a == a.double().double())
    }

    #[test]
    fn test_is_on_curve() {
        let X = Fq([
            0x034365c8, 0x06b2a874, 0x0eb875d7, 0x0ae4c7a7, 0x0785df04, 0x09929351, 0x01fe8c3b,
            0x0f2a0e5f, 0x0111d39c, 0x07ab52ba, 0x01df4552, 0x01d87566, 0x0f297be2, 0x027c090f,
            0x0a81b155, 0x0d1a562b,
        ]);
        let Y = Fq([
            0x00da9708, 0x0e7d583e, 0x0dbcc099, 0x0d2dad89, 0x05a49ce4, 0x01cb4ddc, 0x0928d395,
            0x0098d91d, 0x0bff16ce, 0x06f02f9a, 0x0ce27cc1, 0x0dab5783, 0x0b553d94, 0x03251a0c,
            0x064d70fb, 0x07fe3a2f,
        ]);
        let Z = Fq([
            0x0d5237cc, 0x0319d105, 0x02ab2df5, 0x022e9736, 0x0d79742f, 0x00688712, 0x012d3a65,
            0x0ef4925e, 0x0bd0d260, 0x0832b532, 0x05faef27, 0x01ffe567, 0x0161ce73, 0x07bda0f5,
            0x035d04f1, 0x0930f532,
        ]);
        let T = Fq([
            0x01f6cc27, 0x09be7b8a, 0x0226da79, 0x0f6202f1, 0x0e7264dc, 0x0d25aeb1, 0x06c81f07,
            0x03c32cdc, 0x0923c854, 0x0cfc9865, 0x055b2fed, 0x05bdcc90, 0x01a99835, 0x0ea08056,
            0x0abbf763, 0x03826c2f,
        ]);

        let valid_point = ExtendedPoint { X, Y, Z, T };
        assert!(valid_point.is_on_curve());

        let X = Fq([
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff,
        ]);
        let Y = Fq([
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff,
        ]);
        let Z = Fq([
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff,
        ]);
        let T = Fq([
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff,
        ]);

        let invalid_point = ExtendedPoint { X, Y, Z, T };
        assert!(!invalid_point.is_on_curve());
    }

    #[test]
    fn test_point_add() {
        let X = Fq([
            0x0a7f964a, 0x0db033f8, 0x062b9f0b, 0x07bff7d6, 0x05e755a2, 0x013b6f8b, 0x0f080bdc,
            0x0a112ac0, 0x0416988a, 0x03404b2f, 0x00561ea3, 0x01df752c, 0x070e0b1c, 0x0e73a0c4,
            0x078245d5, 0x09a42df0,
        ]);
        let Y = Fq([
            0x0c2e6c3d, 0x0a03c3f2, 0x0fd16e97, 0x0bab4ec6, 0x08ddba78, 0x091638ef, 0x0b0add85,
            0x070c212d, 0x04bcd337, 0x0c828579, 0x0712cfff, 0x09c1534a, 0x0119cafe, 0x08e72ee0,
            0x0f14ff19, 0x0d0c7e25,
        ]);
        let Z = Fq([
            0x0a0d6be1, 0x0bcd9788, 0x00f9ca8a, 0x038cf839, 0x00912da2, 0x0a3c503a, 0x056fe7e0,
            0x03db9a49, 0x0f19d062, 0x052ac631, 0x01cbda35, 0x02967214, 0x0eed2db2, 0x0a948ce0,
            0x05f7a3a7, 0x0fa35bc2,
        ]);
        let T = Fq([
            0x0fc9f32d, 0x0e442931, 0x065e50ff, 0x04be230d, 0x0dc923c2, 0x0000467c, 0x08fc8902,
            0x0e034cfb, 0x0126370c, 0x06ec706d, 0x06ff07ad, 0x0a27cd65, 0x060f214f, 0x0eb7756d,
            0x0b694dc7, 0x015705ad,
        ]);
        let a = ExtendedPoint { X, Y, Z, T };

        let X = Fq([
            0x075ee82f, 0x0078024b, 0x0a72cc37, 0x07b7b942, 0x01dc37cd, 0x05b2ca96, 0x0fa5deaf,
            0x071020de, 0x09122cbe, 0x01bdbe1d, 0x0eeb69f3, 0x073d88cf, 0x0777b71a, 0x0aa1660a,
            0x0c4476bf, 0x08e2cf30,
        ]);
        let Y = Fq([
            0x0aab8840, 0x0f0932b3, 0x0088be9e, 0x0c4d55d6, 0x01926f40, 0x01c112e0, 0x0884dc6d,
            0x0e66b50b, 0x09120ee4, 0x0750ee39, 0x048c6ce2, 0x00f9fe35, 0x0f74988e, 0x05693a13,
            0x0c1d267c, 0x052d5ba0,
        ]);
        let Z = Fq([
            0x043efd14, 0x07ce59a0, 0x0f9b7154, 0x05663cbd, 0x055ba08f, 0x0525f2b7, 0x0e1d908d,
            0x06d7d26a, 0x0c4cee28, 0x02039ee7, 0x0a733b28, 0x01be5db0, 0x056e9a37, 0x0db1b9b6,
            0x088880cd, 0x03d26863,
        ]);
        let T = Fq([
            0x0b3765ab, 0x0ed3e150, 0x02134041, 0x0ec8f519, 0x0acb91c3, 0x0f916fd5, 0x099a9e35,
            0x0e44da01, 0x0c16e971, 0x0c1b213b, 0x0e824448, 0x0b197385, 0x07988fd4, 0x0ab877a6,
            0x0d658e39, 0x0cf66684,
        ]);
        let b = ExtendedPoint { X, Y, Z, T };

        let X = Fq([
            0x0d6a2637, 0x0cee14c7, 0x0b626f81, 0x02a8151b, 0x01d9b4a2, 0x00c23d77, 0x0859f8bc,
            0x0e02e853, 0x0314bf95, 0x04447761, 0x09fb76bd, 0x0dd35230, 0x03b3f440, 0x017dc316,
            0x09bf7799, 0x054f1cc4,
        ]);
        let Y = Fq([
            0x021dd329, 0x019d887a, 0x0029b32a, 0x06a01e7d, 0x06081030, 0x036080a2, 0x05c8240b,
            0x0c11a3bc, 0x00a2ecfe, 0x045ecf89, 0x08e0d084, 0x06f80067, 0x0b9d1318, 0x0b8bfeb3,
            0x07487524, 0x04e1609a,
        ]);
        let Z = Fq([
            0x067eb923, 0x0d7bbefe, 0x0cf769ff, 0x05725ec8, 0x0c23c0ad, 0x091bcba2, 0x0de48aec,
            0x02c71185, 0x0e607ca2, 0x042ba874, 0x08a195fa, 0x04386d91, 0x079778f4, 0x045ecac6,
            0x02c493ab, 0x050614a7,
        ]);
        let T = Fq([
            0x0f69cb7d, 0x0b9c1d3c, 0x087a352c, 0x06c4b483, 0x025fb591, 0x0bf2bd90, 0x08f3ffa4,
            0x008659f1, 0x062cc1fd, 0x0892fc5b, 0x05f37db7, 0x017d17a8, 0x01b92f20, 0x088ba2bd,
            0x00e3ed5f, 0x0230bb39,
        ]);
        let expected_c = ExtendedPoint { X, Y, Z, T };

        // A + B = B + A = C
        let c = a.to_extensible().add_extended(&b).to_extended();
        assert!(c == expected_c);
        let c = b.to_extensible().add_extended(&a).to_extended();
        assert!(c == expected_c);

        // Adding identity point should not change result
        let c = c.to_extensible().add_extended(&ExtendedPoint::identity());
        let c_commute = b.to_extensible().add_extended(&a);
        assert!(c == c_commute);
    }

    #[test]
    fn test_point_sub() {
        let X = Fq([
            0x0a7f964a, 0x0db033f8, 0x062b9f0b, 0x07bff7d6, 0x05e755a2, 0x013b6f8b, 0x0f080bdc,
            0x0a112ac0, 0x0416988a, 0x03404b2f, 0x00561ea3, 0x01df752c, 0x070e0b1c, 0x0e73a0c4,
            0x078245d5, 0x09a42df0,
        ]);
        let Y = Fq([
            0x0c2e6c3d, 0x0a03c3f2, 0x0fd16e97, 0x0bab4ec6, 0x08ddba78, 0x091638ef, 0x0b0add85,
            0x070c212d, 0x04bcd337, 0x0c828579, 0x0712cfff, 0x09c1534a, 0x0119cafe, 0x08e72ee0,
            0x0f14ff19, 0x0d0c7e25,
        ]);
        let Z = Fq([
            0x0a0d6be1, 0x0bcd9788, 0x00f9ca8a, 0x038cf839, 0x00912da2, 0x0a3c503a, 0x056fe7e0,
            0x03db9a49, 0x0f19d062, 0x052ac631, 0x01cbda35, 0x02967214, 0x0eed2db2, 0x0a948ce0,
            0x05f7a3a7, 0x0fa35bc2,
        ]);
        let T = Fq([
            0x0fc9f32d, 0x0e442931, 0x065e50ff, 0x04be230d, 0x0dc923c2, 0x0000467c, 0x08fc8902,
            0x0e034cfb, 0x0126370c, 0x06ec706d, 0x06ff07ad, 0x0a27cd65, 0x060f214f, 0x0eb7756d,
            0x0b694dc7, 0x015705ad,
        ]);
        let a = ExtendedPoint { X, Y, Z, T };

        let X = Fq([
            0x075ee82f, 0x0078024b, 0x0a72cc37, 0x07b7b942, 0x01dc37cd, 0x05b2ca96, 0x0fa5deaf,
            0x071020de, 0x09122cbe, 0x01bdbe1d, 0x0eeb69f3, 0x073d88cf, 0x0777b71a, 0x0aa1660a,
            0x0c4476bf, 0x08e2cf30,
        ]);
        let Y = Fq([
            0x0aab8840, 0x0f0932b3, 0x0088be9e, 0x0c4d55d6, 0x01926f40, 0x01c112e0, 0x0884dc6d,
            0x0e66b50b, 0x09120ee4, 0x0750ee39, 0x048c6ce2, 0x00f9fe35, 0x0f74988e, 0x05693a13,
            0x0c1d267c, 0x052d5ba0,
        ]);
        let Z = Fq([
            0x043efd14, 0x07ce59a0, 0x0f9b7154, 0x05663cbd, 0x055ba08f, 0x0525f2b7, 0x0e1d908d,
            0x06d7d26a, 0x0c4cee28, 0x02039ee7, 0x0a733b28, 0x01be5db0, 0x056e9a37, 0x0db1b9b6,
            0x088880cd, 0x03d26863,
        ]);
        let T = Fq([
            0x0b3765ab, 0x0ed3e150, 0x02134041, 0x0ec8f519, 0x0acb91c3, 0x0f916fd5, 0x099a9e35,
            0x0e44da01, 0x0c16e971, 0x0c1b213b, 0x0e824448, 0x0b197385, 0x07988fd4, 0x0ab877a6,
            0x0d658e39, 0x0cf66684,
        ]);
        let b = ExtendedPoint { X, Y, Z, T };

        let X = Fq([
            0x0e264012, 0x0c218ff9, 0x06393c0f, 0x0864e62d, 0x05f0e534, 0x0267756d, 0x0ce40403,
            0x0e9e240d, 0x09597584, 0x027844b2, 0x0f8842bf, 0x01b5f03d, 0x05fbfd9a, 0x0e4ed5e3,
            0x07087964, 0x07dc52d0,
        ]);
        let Y = Fq([
            0x05484245, 0x0f3c416d, 0x083a1e46, 0x05e6a9d8, 0x05bfedad, 0x0a9a7379, 0x00b489c3,
            0x0de89d6b, 0x04e7709d, 0x0149bd11, 0x017eb71a, 0x0223de4a, 0x00d9bd0d, 0x093c76a6,
            0x072fe435, 0x0d6fd2c5,
        ]);
        let Z = Fq([
            0x067eb923, 0x0d7bbefe, 0x0cf769ff, 0x05725ec8, 0x0c23c0ad, 0x091bcba2, 0x0de48aec,
            0x02c71185, 0x0e607ca2, 0x042ba874, 0x08a195fa, 0x04386d91, 0x079778f4, 0x045ecac6,
            0x02c493ab, 0x050614a7,
        ]);
        let T = Fq([
            0x0455a73b, 0x0cfbe5f2, 0x0cdb56a2, 0x06477b21, 0x0fda6909, 0x07f6faeb, 0x04ebea8b,
            0x0d1e04b7, 0x00307c2a, 0x0e926e5c, 0x0efdf04c, 0x038841bb, 0x09be04e8, 0x001137e1,
            0x0515b17d, 0x0ea27de2,
        ]);
        let expected_c = ExtendedPoint { X, Y, Z, T };

        let c = a.to_extensible().sub_extended(&b);
        assert!(c.to_extended() == expected_c);
    }

    #[test]
    fn test_negate() {
        let X = Fq([
            0x0a7f964a, 0x0db033f8, 0x062b9f0b, 0x07bff7d6, 0x05e755a2, 0x013b6f8b, 0x0f080bdc,
            0x0a112ac0, 0x0416988a, 0x03404b2f, 0x00561ea3, 0x01df752c, 0x070e0b1c, 0x0e73a0c4,
            0x078245d5, 0x09a42df0,
        ]);
        let Y = Fq([
            0x0c2e6c3d, 0x0a03c3f2, 0x0fd16e97, 0x0bab4ec6, 0x08ddba78, 0x091638ef, 0x0b0add85,
            0x070c212d, 0x04bcd337, 0x0c828579, 0x0712cfff, 0x09c1534a, 0x0119cafe, 0x08e72ee0,
            0x0f14ff19, 0x0d0c7e25,
        ]);
        let Z = Fq([
            0x0a0d6be1, 0x0bcd9788, 0x00f9ca8a, 0x038cf839, 0x00912da2, 0x0a3c503a, 0x056fe7e0,
            0x03db9a49, 0x0f19d062, 0x052ac631, 0x01cbda35, 0x02967214, 0x0eed2db2, 0x0a948ce0,
            0x05f7a3a7, 0x0fa35bc2,
        ]);
        let T = Fq([
            0x0fc9f32d, 0x0e442931, 0x065e50ff, 0x04be230d, 0x0dc923c2, 0x0000467c, 0x08fc8902,
            0x0e034cfb, 0x0126370c, 0x06ec706d, 0x06ff07ad, 0x0a27cd65, 0x060f214f, 0x0eb7756d,
            0x0b694dc7, 0x015705ad,
        ]);
        let a = ExtendedPoint { X, Y, Z, T };

        let X = Fq([
            0x058069b5, 0x024fcc07, 0x09d460f4, 0x08400829, 0x0a18aa5d, 0x0ec49074, 0x00f7f423,
            0x05eed53f, 0x0be96774, 0x0cbfb4d0, 0x0fa9e15c, 0x0e208ad3, 0x08f1f4e3, 0x018c5f3b,
            0x087dba2a, 0x065bd20f,
        ]);
        let Y = Fq([
            0x0c2e6c3d, 0x0a03c3f2, 0x0fd16e97, 0x0bab4ec6, 0x08ddba78, 0x091638ef, 0x0b0add85,
            0x070c212d, 0x04bcd337, 0x0c828579, 0x0712cfff, 0x09c1534a, 0x0119cafe, 0x08e72ee0,
            0x0f14ff19, 0x0d0c7e25,
        ]);
        let Z = Fq([
            0x0a0d6be1, 0x0bcd9788, 0x00f9ca8a, 0x038cf839, 0x00912da2, 0x0a3c503a, 0x056fe7e0,
            0x03db9a49, 0x0f19d062, 0x052ac631, 0x01cbda35, 0x02967214, 0x0eed2db2, 0x0a948ce0,
            0x05f7a3a7, 0x0fa35bc2,
        ]);
        let T = Fq([
            0x00360cd2, 0x01bbd6ce, 0x09a1af00, 0x0b41dcf2, 0x0236dc3d, 0x0fffb983, 0x070376fd,
            0x01fcb304, 0x0ed9c8f2, 0x09138f92, 0x0900f852, 0x05d8329a, 0x09f0deb0, 0x01488a92,
            0x0496b238, 0x0ea8fa52,
        ]);
        let expected_neg_a = ExtendedPoint { X, Y, Z, T };

        let neg_a = a.negate();
        assert!(neg_a == expected_neg_a);

        assert!(a.to_extensible().add_extended(&neg_a) == ExtensiblePoint::identity());
    }

    #[test]
    fn test_double() {
        let X = Fq([
            0x08354b7a, 0x0895b3e8, 0x06ae5175, 0x0644b394, 0x0b7faf9e, 0x0c5237db, 0x013a0c90,
            0x08f5bce0, 0x09a3d79b, 0x00f17559, 0x0de8f041, 0x073e222f, 0x0dc2b7ee, 0x005ac354,
            0x0766db38, 0x065631fe,
        ]);
        let Y = Fq([
            0x00398885, 0x055c9bed, 0x0ae443ca, 0x0fd70ea4, 0x09e2a7d2, 0x04ac2e9d, 0x00678287,
            0x0294768e, 0x0b604cea, 0x07b49317, 0x0dc2a6d9, 0x0e44a6fb, 0x09db3965, 0x049d3bf5,
            0x03e655fe, 0x003a9c02,
        ]);
        let Z = Fq([
            0x0fd57162, 0x0a39f768, 0x03009756, 0x065d735f, 0x0d1da282, 0x0589ecd7, 0x003196b1,
            0x0c001dfe, 0x019f1050, 0x0152e8d2, 0x0c14ff38, 0x00f7a446, 0x028053f6, 0x0f8a91e9,
            0x05a8d694, 0x09d5ae86,
        ]);
        let T = Fq([
            0x04198f2e, 0x0d82440f, 0x0fce100e, 0x0af4829d, 0x0d5c3516, 0x0094a0da, 0x078cdb39,
            0x0e738836, 0x01ec536d, 0x06dfd1e9, 0x0ee16173, 0x0addc8c0, 0x0797fb1d, 0x059741a3,
            0x0a7f9c34, 0x088fe0a6,
        ]);
        let a = ExtendedPoint { X, Y, Z, T };

        let X = Fq([
            0x00d8f04c, 0x03e54689, 0x0eb4db2b, 0x0887ba34, 0x0a5b4ebc, 0x0f6c0261, 0x03bfa803,
            0x0408ff02, 0x03b4ef26, 0x0465c028, 0x0cd47378, 0x064c55b4, 0x08245850, 0x01912682,
            0x0dcbf92c, 0x07a7fa30,
        ]);
        let Y = Fq([
            0x0d94d1a6, 0x0f7306e8, 0x0278b336, 0x04362b7b, 0x0faf02b9, 0x06b01d18, 0x07a597da,
            0x0bd6add0, 0x047afa98, 0x0e64e897, 0x0bbf88e6, 0x01d0a534, 0x04a52b9d, 0x0af374e0,
            0x05091d54, 0x00fcf1a5,
        ]);
        let Z = Fq([
            0x042318ce, 0x04aecdae, 0x0e8f196b, 0x0019d2e3, 0x045d147c, 0x060b153e, 0x0adf2c37,
            0x0419cdd8, 0x06d19046, 0x00d18821, 0x06c7b9c2, 0x0c0ffd68, 0x0b7e4ca2, 0x06da0d56,
            0x0952b40f, 0x03008395,
        ]);
        let T = Fq([
            0x04643593, 0x000e0fdd, 0x013f29f3, 0x0bb8992d, 0x0a30d344, 0x09151eec, 0x0d12bb82,
            0x05c7a054, 0x0103c2c6, 0x08a61fe2, 0x0aced4bf, 0x0f76d481, 0x0db774be, 0x065ef8a8,
            0x0ff47a71, 0x0f49f73e,
        ]);
        let expected_double_a = ExtendedPoint { X, Y, Z, T };

        let double_a = a.double();

        assert!(double_a == expected_double_a);
    }
}
