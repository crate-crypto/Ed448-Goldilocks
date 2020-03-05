use std::ops::Add;

// This is the scalar field
// size = 4q = 2^446 - 0x8335dc163bb124b65129c96fde933d8d723a70aadc873d6d54a7bb0d
// We can therefore use 14 saturated 32-bit limbs
// in LE format
/// note: Montgomery reduction and barret are both used: montgomery when multiplying and barret when decoding
/// Probably easier to stick to montgomery I think
#[derive(Debug, PartialEq, Eq)]
pub struct Scalar([u32; 14]);

const MODULUS: Scalar = Scalar([
    0xab5844f3, 0x2378c292, 0x8dc58f55, 0x216cc272, 0xaed63690, 0xc44edb49, 0x7cca23e9, 0xffffffff,
    0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x3fffffff,
]);
impl From<u32> for Scalar {
    fn from(a: u32) -> Scalar {
        Scalar([a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
}

impl Scalar {
    pub fn one() -> Scalar {
        Scalar::from(1)
    }
    pub fn zero() -> Scalar {
        Scalar::from(0)
    }
    pub fn invert(&self) -> Self {
        todo!()
    }

    fn sub_extra(minuend: &Scalar, subtrahend: &Scalar, carry: u32) -> Scalar {
        let mut result = Scalar::zero();
        let mut chain = 0i64;

        for i in 0..14 {
            chain += minuend.0[i] as i64 - subtrahend.0[i] as i64;
            result.0[i] = chain as u32;
            chain >>= 32;
        }

        let borrow = chain + carry as i64;
        chain = 0;

        for i in 0..14 {
            chain += result.0[i] as i64 + (MODULUS.0[i] as i64) & borrow;
            result.0[i] = chain as u32;
            chain >>= 32;
        }

        result
    }
}

impl Add<Scalar> for Scalar {
    type Output = Scalar;
    fn add(self, rhs: Scalar) -> Self::Output {
        let mut chain = 0u64;
        let mut result = Scalar::zero();
        for i in 0..14 {
            chain += self.0[i] as u64 + rhs.0[i] as u64;
            result.0[i] = chain as u32;
            chain >>= 32
        }
        Self::sub_extra(&result, &MODULUS, chain as u32)
    }
}
#[test]
fn test_basic_add() {
    let five = Scalar::from(5);
    let six = Scalar::from(6);

    assert_eq!(five + six, Scalar::from(11))
}
