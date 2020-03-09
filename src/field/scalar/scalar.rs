use super::constants::NUM_LIMBS;
use std::ops::{Add, Mul, Sub};
/// This is the scalar field
/// size = 4q = 2^446 - 0x8335dc163bb124b65129c96fde933d8d723a70aadc873d6d54a7bb0d
/// We can therefore use 14 saturated 32-bit limbs
/// in LE format

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct Scalar([u32; 14]);

const MODULUS: Scalar = Scalar([
    0xab5844f3, 0x2378c292, 0x8dc58f55, 0x216cc272, 0xaed63690, 0xc44edb49, 0x7cca23e9, 0xffffffff,
    0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x3fffffff,
]);
// Montomomery R^2
const R2: Scalar = Scalar([
    0x049b9b60, 0xe3539257, 0xc1b195d9, 0x7af32c4b, 0x88ea1859, 0x0d66de23, 0x5ee4d838, 0xae17cf72,
    0xa3c47c44, 0x1a9cc14b, 0xe4d070af, 0x2052bcb7, 0xf823b729, 0x3402a939,
]);
impl From<u32> for Scalar {
    fn from(a: u32) -> Scalar {
        Scalar([a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
}

// Trait implementations

impl Add<Scalar> for Scalar {
    type Output = Scalar;
    fn add(self, rhs: Scalar) -> Self::Output {
        add(&self, &rhs)
    }
}
impl Mul<Scalar> for Scalar {
    type Output = Scalar;
    fn mul(self, rhs: Scalar) -> Self::Output {
        let unreduced = montgomery_multiply(&self, &rhs);
        montgomery_multiply(&unreduced, &R2)
    }
}
impl Sub<Scalar> for Scalar {
    type Output = Scalar;
    fn sub(self, rhs: Scalar) -> Self::Output {
        sub_extra(&self, &rhs, 0)
    }
}

impl Scalar {
    pub fn one() -> Scalar {
        Scalar::from(1)
    }
    pub fn zero() -> Scalar {
        Scalar::from(0)
    }
    fn square(&self) -> Scalar {
        *self * *self
    }
    pub fn invert(&self) -> Self {
        todo!()
    }
    /// Halves a Scalar
    // XXX: What is expected output on odd Scalars
    pub fn halve(&self) -> Self {
        let mut result = Scalar::zero();

        let mask = -((self.0[0] & 1) as i64);
        let mut chain = 0u64;

        for i in 0..14 {
            chain += (self.0[i] as u64) + ((MODULUS.0[i] as i64 & mask) as u64);
            result.0[i] = chain as u32;
            chain >>= 32
        }

        for i in 0..13 {
            result.0[i] = (result.0[i] >> 1) | (result.0[i] << 31);
        }
        result.0[13] = (result.0[13] >> 1) | ((chain << 31) as u32);
        result
    }
}
/// Computes a + b mod p
pub fn add(a: &Scalar, b: &Scalar) -> Scalar {
    // First add the two Scalars together
    // Since our limbs are saturated, the result of each
    // limb being added can be a 33-bit integer so we propagate the carry bit
    let mut result = Scalar::zero();

    // a + b
    let mut chain = 0u64;
    for i in 0..14 {
        chain += (a.0[i] as u64) + (b.0[i] as u64);
        // Low 32 bits are the results
        result.0[i] = chain as u32;
        // 33rd bit is the carry
        chain >>= 32;
    }

    // Now reduce the results
    sub_extra(&result, &MODULUS, chain as u32)
}

/// Compute a - b mod p
/// Computes a - b and conditionally computes the modulus if the result was negative
fn sub_extra(a: &Scalar, b: &Scalar, carry: u32) -> Scalar {
    let mut result = Scalar::zero();

    // a - b
    let mut chain = 0i64;
    for i in 0..14 {
        chain += a.0[i] as i64 - b.0[i] as i64;
        // Low 32 bits are the results
        result.0[i] = chain as u32;
        // 33rd bit is the borrow
        chain >>= 32
    }

    // if the result of a-b was negative and carry was zero
    // then borrow will be 0xfff..fff and the modulus will be added conditionally to the result
    // If the carry was 1 and a-b was not negative, then the borrow will be 0x00000...001 ( this should not happen)
    // Since the borrow should never be more than 0, the carry should never be more than 1;
    // XXX: Explain why the case of borrow == 1 should never happen
    let borrow = chain + (carry as i64);
    assert!(borrow == -1 || borrow == 0);

    chain = 0i64;
    for i in 0..14 {
        chain += (result.0[i] as i64) + ((MODULUS.0[i] as i64) & borrow);
        // Low 32 bits are the results
        result.0[i] = chain as u32;
        // 33rd bit is the carry
        chain >>= 32;
    }

    result
}

fn montgomery_multiply(x: &Scalar, y: &Scalar) -> Scalar {
    const MONTGOMERY_FACTOR: u32 = 0xae918bc5;

    let mut result = Scalar::zero();
    let mut carry = 0u32;

    // (a * b ) + c
    let mul_add = |a: u32, b: u32, c: u32| -> u64 { ((a as u64) * (b as u64)) + (c as u64) };

    for i in 0..14 {
        let mut chain = 0u64;
        for j in 0..14 {
            chain += mul_add(x.0[i], y.0[j], result.0[j]);
            result.0[j] = chain as u32;
            chain >>= 32;
        }

        let saved = chain as u32;
        let multiplicand = result.0[0].wrapping_mul(MONTGOMERY_FACTOR);
        chain = 0u64;

        for j in 0..14 {
            chain += mul_add(multiplicand, MODULUS.0[j], result.0[j]);
            if j > 0 {
                result.0[j - 1] = chain as u32;
            }
            chain >>= 32;
        }
        chain += (saved as u64) + (carry as u64);
        result.0[NUM_LIMBS - 1] = chain as u32;
        carry = (chain >> 32) as u32;
    }

    sub_extra(&result, &MODULUS, carry)
}

#[test]
fn test_basic_add() {
    let five = Scalar::from(5);
    let six = Scalar::from(6);

    assert_eq!(five + six, Scalar::from(11))
}

#[test]
fn test_basic_sub() {
    let ten = Scalar::from(10);
    let five = Scalar::from(5);
    assert_eq!(ten - five, Scalar::from(5))
}

#[test]
fn test_basic_mul() {
    let ten = Scalar::from(10);
    let five = Scalar::from(5);

    assert_eq!(ten * five, Scalar::from(50))
}

#[test]
fn test_basic_square() {
    let five = Scalar::from(5);
    assert_eq!(five.square(), Scalar::from(25))
}
#[test]
fn test_basic_halving() {
    let eight = Scalar::from(8);
    let four = Scalar::from(4);
    let two = Scalar::from(2);
    assert_eq!(eight.halve(), four);
    assert_eq!(four.halve(), two);
    assert_eq!(two.halve(), Scalar::one());
}
