use crate::curve::constants::{ONE_MINUS_D, TWO_D_MINUS_ONE, TWO_ONE_MINUS_D};
use crate::curve::edwards::affine::{AffineNielsPoint, AffinePoint};
use crate::curve::edwards::projective::{ProjectiveNielsPoint, ProjectivePoint};
use crate::field::{base::Fq, scalar::Scalar};

/// Represent points on the twisted edwards curve using Extended Homogenous Projective Co-ordinates
/// (X : Y : T : Z) this corresponds to the affine (X/Z, Y/Z, T/Z) with Z != 0 and T = XY
/// E_d : (Y/Z)^2 - (X/Z)^2 = 1 + d(X/Z)^2(Y/Z)^2
#[derive(Copy, Clone, Debug)]
pub struct ExtendedPoint {
    pub(crate) X: Fq,
    pub(crate) Y: Fq,
    pub(crate) Z: Fq,
    pub(crate) T: Fq,
}

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
    // XXX: This is _almost_ unreadable, we will slowly refactor for interopability
    pub fn scalar_mul(&self, s: &Scalar) -> ExtendedPoint {
        use crate::window::wnaf;
        use crate::window::BASE_TABLE;

        let mut result = ExtendedPoint {
            X: Fq::zero(),
            Y: Fq::zero(),
            Z: Fq::zero(),
            T: Fq::zero(),
        };

        let mut scalar_one_x = *s + BASE_TABLE.scalar_adjustment;
        // XXX: Halve method does not seem to be working correctly
        scalar_one_x = scalar_one_x * Scalar::from(2).invert();

        let multiples = wnaf::prepare_fixed_window(self);

        let mut first = true;
        let loop_end = 446 - (445 % wnaf::WINDOW) - 1;
        let loop_start = 0;

        for i in (loop_start..=loop_end).rev().step_by(wnaf::WINDOW) {
            let mut bits = scalar_one_x[i / 32] >> (i % 32);
            if (i % 32 >= 32 - wnaf::WINDOW) && (i / 32 < (14 - 1)) {
                bits ^= scalar_one_x[(i / 32) + 1] << (32 - (i % 32));
            }
            bits &= wnaf::WINDOW_MASK as u32;
            let inv = (bits >> (wnaf::WINDOW - 1)).wrapping_sub(1);
            bits ^= inv;

            let mut neg_P = wnaf::lookup(multiples, bits & (wnaf::WINDOW_T_MASK as u32));
            neg_P.conditional_negate(inv);
            if first {
                result = neg_P.to_extended();
                first = false;
            } else {
                for _ in 0..wnaf::WINDOW - 1 {
                    result = result.double_internal(true);
                }
                result = result.double_internal(false);
                result = result.add_projective_niels(&neg_P, false);
            }
        }

        result
    }

    pub fn equals(&self, other: &ExtendedPoint) -> bool {
        // Assumes (0,0,0,0) is not possible
        let XZ = self.X * other.Z;
        let ZX = self.Z * other.X;

        let YZ = self.Y * other.Z;
        let ZY = self.Z * other.Y;

        XZ.equals(&ZX) & YZ.equals(&ZY)
    }
    pub fn to_projective(&self) -> ProjectivePoint {
        ProjectivePoint {
            X: self.X,
            Y: self.Y,
            Z: self.Z,
        }
    }
    pub fn to_affine(&self) -> AffinePoint {
        if self.equals(&ExtendedPoint::identity()) || self.Z.equals(&Fq::zero()) {
            return AffinePoint {
                x: Fq::zero(),
                y: Fq::zero(),
            };
        }
        let r = self.Z.invert();
        let s = r.square();

        let mut x = self.X * s;
        x.strong_reduce();

        let t = self.Y * s;

        let mut y = t * r;
        y.strong_reduce();

        AffinePoint { x, y }
    }
    pub fn to_projective_niels(&self) -> ProjectiveNielsPoint {
        ProjectiveNielsPoint {
            y_plus_x: self.X + self.Y,
            y_minus_x: self.Y - self.X,
            Z: self.Z + self.Z,
            td: self.T * TWO_D_MINUS_ONE,
        }
    }

    pub fn add_affine_niels(&self, other: AffineNielsPoint, before_double: bool) -> ExtendedPoint {
        let mut a = Fq::zero();
        let mut b = Fq::zero();
        let mut c = Fq::zero();

        let mut X = self.X;
        let mut Y = self.Y;
        let mut Z = self.Z;
        let mut T = self.T;

        b = Y - X;
        a = other.y_minus_x * b;
        b = X.add_no_reduce(&Y);
        Y = other.y_plus_x * b;
        X = other.td * T;
        c = a.add_no_reduce(&Y);
        b = Y - a;
        Y = Z - X;
        a = X.add_no_reduce(&Z);
        Z = a * Y;
        X = Y * b;
        Y = a * c;
        if !before_double {
            T = b * c;
        }

        ExtendedPoint { X, Y, Z, T }
    }
    pub fn add_projective_niels(
        &self,
        other: &ProjectiveNielsPoint,
        before_double: bool,
    ) -> ExtendedPoint {
        let mut result = self.clone();
        result.Z = other.Z * result.Z;

        result.add_affine_niels(
            AffineNielsPoint {
                y_plus_x: other.y_plus_x,
                y_minus_x: other.y_minus_x,
                td: other.td,
            },
            before_double,
        )
    }
    pub(crate) fn is_on_curve(&self) -> bool {
        let XY = self.X * self.Y;
        let ZT = self.Z * self.T;

        // y^2 - x^2 == z^2 - t^2 * (1 - D)
        // XX: Check^

        let YY = self.Y.square();
        let XX = self.X.square();
        let ZZ = self.Z.square();
        let TT = self.T.square();
        let lhs = YY - XX;
        let rhs = ZZ - TT * ONE_MINUS_D;

        (XY == ZT) && (lhs == rhs)
    }
    pub fn add(&self, other: &ExtendedPoint) -> ExtendedPoint {
        let mut result = ExtendedPoint {
            X: Fq::zero(),
            Y: Fq::zero(),
            Z: Fq::zero(),
            T: Fq::zero(),
        };

        let (mut a, mut b, mut c, mut d) = (Fq::zero(), Fq::zero(), Fq::zero(), Fq::zero());

        b = self.Y - self.X;
        c = other.Y - other.X;
        d = other.Y.add_no_reduce(&other.X);
        a = c * b;
        b = self.Y.add_no_reduce(&self.X);
        result.Y = d * b;
        b = other.T * self.T;
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
        result.T = b * c;

        result
    }
    pub fn negate(&self) -> ExtendedPoint {
        ExtendedPoint {
            X: self.X.negate(),
            Y: self.Y,
            Z: self.Z,
            T: self.T.negate(),
        }
    }
    pub fn sub(&self, other: &ExtendedPoint) -> ExtendedPoint {
        let mut result = ExtendedPoint {
            X: Fq::zero(),
            Y: Fq::zero(),
            Z: Fq::zero(),
            T: Fq::zero(),
        };

        let (mut a, mut b, mut c, mut d) = (Fq::zero(), Fq::zero(), Fq::zero(), Fq::zero());

        b = self.Y - self.X;
        d = other.Y - other.X;
        c = other.Y.add_no_reduce(&other.X);
        a = c * b;
        b = self.Y.add_no_reduce(&self.X);
        result.Y = d * b;
        b = other.T * self.T;
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
        result.T = b * c;

        result
    }

    pub(crate) fn double_internal(&self, before_double: bool) -> ExtendedPoint {
        let mut result = ExtendedPoint {
            X: self.X,
            Y: self.Y,
            Z: self.Z,
            T: self.T,
        };

        let (mut a, mut b, mut c, mut d) = (Fq::zero(), Fq::zero(), Fq::zero(), Fq::zero());

        c = result.X.square();
        a = result.Y.square();
        d = c.add_no_reduce(&a);
        result.T = result.Y.add_no_reduce(&result.X);
        b = result.T.square();
        b = b.sub_no_reduce(&d);
        b.bias(3);
        b.weak_reduce();
        result.T = a - c;
        result.X = result.Z.square();
        result.Z = result.X.add_no_reduce(&result.X);
        a = result.Z.sub_no_reduce(&result.T);
        a.bias(4);
        a.weak_reduce();
        result.X = a * b;
        result.Z = result.T * a;
        result.Y = result.T * d;
        if !before_double {
            result.T = b * d;
        }

        result
    }

    pub fn double(&self) -> ExtendedPoint {
        self.double_internal(false)
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
        let c = a.add(&b);
        assert!(c.equals(&expected_c));
        let c = b.add(&a);
        assert!(c.equals(&expected_c));
        // Add identity point
        let c = c.add(&ExtendedPoint::identity());
        assert!(c.equals(&c));
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

        let c = a.sub(&b);
        assert!(c.equals(&expected_c));
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
        assert!(neg_a.equals(&expected_neg_a));

        assert!(a.add(&neg_a).equals(&ExtendedPoint::identity()));
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

        let double_a = a.double_internal(false);

        assert!(double_a.equals(&expected_double_a));
    }

    #[test]
    fn test_add_affine_niels() {
        let p = ExtendedPoint::identity();
        let q = AffineNielsPoint {
            y_minus_x: Fq::from(0x068d5b74u32),
            y_plus_x: Fq::from(0x068d5b74u32),
            td: Fq::from(0x068d5b74u32),
        };

        let X = Fq([
            0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff,
            0x0fffffff, 0x0ffffffe, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff,
            0x0fffffff, 0x0fffffff,
        ]);

        let Y = Fq([
            0x0d1ab6e7, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff,
            0x0fffffff, 0x0fffffff,
        ]);

        let Z = Fq([
            0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff, 0x0fffffff,
            0x0fffffff, 0x0fffffff,
        ]);

        // XXX: Check, this test seems to be invariant under the value of T
        let T = Fq::zero();

        let expected_result = ExtendedPoint { X, Y, Z, T };

        let result = p.add_affine_niels(q, true);

        assert!(result.equals(&expected_result));
    }

    #[test]
    fn test_scalr_mul() {
        let p = ExtendedPoint {
            X: Fq([
                0x0a7f964a, 0x0db033f8, 0x062b9f0b, 0x07bff7d6, 0x05e755a2, 0x013b6f8b, 0x0f080bdc,
                0x0a112ac0, 0x0416988a, 0x03404b2f, 0x00561ea3, 0x01df752c, 0x070e0b1c, 0x0e73a0c4,
                0x078245d5, 0x09a42df0,
            ]),
            Y: Fq([
                0x0c2e6c3d, 0x0a03c3f2, 0x0fd16e97, 0x0bab4ec6, 0x08ddba78, 0x091638ef, 0x0b0add85,
                0x070c212d, 0x04bcd337, 0x0c828579, 0x0712cfff, 0x09c1534a, 0x0119cafe, 0x08e72ee0,
                0x0f14ff19, 0x0d0c7e25,
            ]),
            Z: Fq([
                0x0a0d6be1, 0x0bcd9788, 0x00f9ca8a, 0x038cf839, 0x00912da2, 0x0a3c503a, 0x056fe7e0,
                0x03db9a49, 0x0f19d062, 0x052ac631, 0x01cbda35, 0x02967214, 0x0eed2db2, 0x0a948ce0,
                0x05f7a3a7, 0x0fa35bc2,
            ]),
            T: Fq([
                0x0fc9f32d, 0x0e442931, 0x065e50ff, 0x04be230d, 0x0dc923c2, 0x0000467c, 0x08fc8902,
                0x0e034cfb, 0x0126370c, 0x06ec706d, 0x06ff07ad, 0x0a27cd65, 0x060f214f, 0x0eb7756d,
                0x0b694dc7, 0x015705ad,
            ]),
        };
        let scalar = Scalar([
            0x6ee372b7, 0xe128ae78, 0x1533427c, 0xad0b7015, 0x307f665e, 0xde8026c1, 0xb64629d1,
            0xab454c66, 0x3fe5bf1a, 0x083f8304, 0x3c003777, 0xdef437f6, 0xee2e1b73, 0x05ca185a,
        ]);

        let expected_point = ExtendedPoint {
            X: Fq([
                0x08630007, 0x0bd755e6, 0x0f76b928, 0x070d9694, 0x0b952009, 0x0cf85b12, 0x0c3a6e9c,
                0x0e2d860e, 0x02fd2901, 0x09a73726, 0x02aa2d4c, 0x06913ea9, 0x090da66d, 0x06a5c6f1,
                0x04cc7a13, 0x0eb24ed8,
            ]),
            Y: Fq([
                0x0bb37152, 0x0a3a36b3, 0x0a720c7f, 0x0e29095f, 0x04e76cf4, 0x0cfad965, 0x07439798,
                0x0f4b7ba4, 0x0316ba61, 0x09389566, 0x07f96104, 0x07bdc39c, 0x0f019987, 0x05416850,
                0x0612c6c8, 0x0e231baa,
            ]),
            Z: Fq([
                0x0179c756, 0x04130eef, 0x07f43255, 0x0cc1534d, 0x03e347fd, 0x0c745e4d, 0x068d7bf5,
                0x020b8465, 0x0356d2f1, 0x069b22fd, 0x0b6cf87f, 0x0edf9761, 0x034f512f, 0x0411b43f,
                0x033f0755, 0x06195e97,
            ]),
            T: Fq([
                0x0866187a, 0x035622be, 0x0b9e2e78, 0x0cae26c6, 0x041c2c41, 0x07296c68, 0x03343d3e,
                0x062c0927, 0x0cf5d263, 0x08db465d, 0x033382d6, 0x0c5e6eff, 0x0c0ded8d, 0x037837bf,
                0x03780cc6, 0x0e2360df,
            ]),
        };
        let got = p.scalar_mul(&scalar);
        assert!(expected_point.equals(&got));
    }

    #[test]
    fn test_simple_scalar_mul_identities() {
        let x = ExtendedPoint {
            X: Fq([
                0x034365c8, 0x06b2a874, 0x0eb875d7, 0x0ae4c7a7, 0x0785df04, 0x09929351, 0x01fe8c3b,
                0x0f2a0e5f, 0x0111d39c, 0x07ab52ba, 0x01df4552, 0x01d87566, 0x0f297be2, 0x027c090f,
                0x0a81b155, 0x0d1a562b,
            ]),
            Y: Fq([
                0x00da9708, 0x0e7d583e, 0x0dbcc099, 0x0d2dad89, 0x05a49ce4, 0x01cb4ddc, 0x0928d395,
                0x0098d91d, 0x0bff16ce, 0x06f02f9a, 0x0ce27cc1, 0x0dab5783, 0x0b553d94, 0x03251a0c,
                0x064d70fb, 0x07fe3a2f,
            ]),
            Z: Fq([
                0x0d5237cc, 0x0319d105, 0x02ab2df5, 0x022e9736, 0x0d79742f, 0x00688712, 0x012d3a65,
                0x0ef4925e, 0x0bd0d260, 0x0832b532, 0x05faef27, 0x01ffe567, 0x0161ce73, 0x07bda0f5,
                0x035d04f1, 0x0930f532,
            ]),
            T: Fq([
                0x01f6cc27, 0x09be7b8a, 0x0226da79, 0x0f6202f1, 0x0e7264dc, 0x0d25aeb1, 0x06c81f07,
                0x03c32cdc, 0x0923c854, 0x0cfc9865, 0x055b2fed, 0x05bdcc90, 0x01a99835, 0x0ea08056,
                0x0abbf763, 0x03826c2f,
            ]),
        };

        // Test that 1 * P = P
        let exp = x.scalar_mul(&Scalar::from(1));
        assert!(x.equals(&exp));

        // Test that 2 * (P + P) = 4 * P
        let expected_two_x = x.add(&x).double();
        let got = x.scalar_mul(&Scalar::from(4));
        assert!(expected_two_x.equals(&got));
    }
}
