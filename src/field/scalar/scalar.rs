use std::ops::{Add, Index, IndexMut, Mul, Sub};
/// This is the scalar field
/// size = 4q = 2^446 - 0x8335dc163bb124b65129c96fde933d8d723a70aadc873d6d54a7bb0d
/// We can therefore use 14 saturated 32-bit limbs
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

impl Index<usize> for Scalar {
    type Output = u32;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}
impl IndexMut<usize> for Scalar {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
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
fn word_is_zero(word: u32) -> bool {
    // If the word is zero, then when we minus 1, we should get 0xffffffff
    word.wrapping_sub(1) == 0xffffffff
}
impl Scalar {
    pub fn one() -> Scalar {
        Scalar::from(1)
    }
    pub fn zero() -> Scalar {
        Scalar::from(0)
    }
    fn square(&self) -> Scalar {
        montgomery_multiply(&self, &self)
    }
    pub fn invert(&self) -> Self {
        let mut pre_comp: Vec<Scalar> = vec![Scalar::zero(); 8];
        let mut result = Scalar::zero();

        let scalar_window_bits = 3;
        let last = (1 << scalar_window_bits) - 1;

        // precompute [a^1, a^3,,..]
        pre_comp[0] = montgomery_multiply(self, &R2);

        if last > 0 {
            pre_comp[last] = montgomery_multiply(&pre_comp[0], &pre_comp[0]);
        }

        for i in 1..=last {
            pre_comp[i] = montgomery_multiply(&pre_comp[i - 1], &pre_comp[last])
        }

        // Sliding window
        let mut residue: usize = 0;
        let mut trailing: usize = 0;
        let mut started: usize = 0;

        // XXX: This can definitely be refactored to be readable
        let loop_start = -scalar_window_bits as isize;
        let loop_end = 446 - 1;
        for i in (loop_start..=loop_end).rev() {
            if started != 0 {
                result = result.square()
            }

            let mut w: u32;
            if i >= 0 {
                w = MODULUS[(i / 32) as usize];
            } else {
                w = 0;
            }

            if i >= 0 && i < 32 {
                w -= 2
            }

            residue = (((residue as u32) << 1) | ((w >> ((i as u32) % 32)) & 1)) as usize;
            if residue >> scalar_window_bits != 0 {
                trailing = residue;
                residue = 0
            }

            if trailing > 0 && (trailing & ((1 << scalar_window_bits) - 1)) == 0 {
                if started != 0 {
                    result = montgomery_multiply(
                        &result,
                        &pre_comp[trailing >> (scalar_window_bits + 1)],
                    )
                } else {
                    result = pre_comp[trailing >> (scalar_window_bits + 1)];
                    started = 1
                }
                trailing = 0
            }
            trailing <<= 1
        }

        // de-montgomerize and return result

        montgomery_multiply(&result, &Scalar::one())
    }
    fn equals(&self, rhs: &Scalar) -> bool {
        let mut diff = 0u32;
        for i in 0..14 {
            diff |= self[i] ^ rhs[i]
        }
        word_is_zero(diff)
    }
    /// Halves a Scalar
    // XXX: What is expected output on odd Scalars
    pub fn halve(&self) -> Self {
        let mut result = Scalar::zero();

        let mask = 0u32.wrapping_sub(self[0] & 1);
        let mut chain = 0u64;

        for i in 0..14 {
            chain += (self[i] as u64) + ((MODULUS[i] & mask) as u64);
            result[i] = chain as u32;
            chain >>= 32
        }

        for i in 0..13 {
            result[i] = (result[i] >> 1) | (result[i] << 31);
        }
        result[13] = (result[13] >> 1) | ((chain << 31) as u32);
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
        chain += (a[i] as u64) + (b[i] as u64);
        // Low 32 bits are the results
        result[i] = chain as u32;
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
        chain += a[i] as i64 - b[i] as i64;
        // Low 32 bits are the results
        result[i] = chain as u32;
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
        chain += (result[i] as i64) + ((MODULUS[i] as i64) & borrow);
        // Low 32 bits are the results
        result[i] = chain as u32;
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
            chain += mul_add(x[i], y[j], result[j]);
            result[j] = chain as u32;
            chain >>= 32;
        }

        let saved = chain as u32;
        let multiplicand = result[0].wrapping_mul(MONTGOMERY_FACTOR);
        chain = 0u64;

        for j in 0..14 {
            chain += mul_add(multiplicand, MODULUS[j], result[j]);
            if j > 0 {
                result[j - 1] = chain as u32;
            }
            chain >>= 32;
        }
        chain += (saved as u64) + (carry as u64);
        result[14 - 1] = chain as u32;
        carry = (chain >> 32) as u32;
    }

    sub_extra(&result, &MODULUS, carry)
}
#[cfg(test)]
mod test {
    use super::*;
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
    fn test_mul() {
        let a = Scalar([
            0xffb823a3, 0xc96a3c35, 0x7f8ed27d, 0x087b8fb9, 0x1d9ac30a, 0x74d65764, 0xc0be082e,
            0xa8cb0ae8, 0xa8fa552b, 0x2aae8688, 0x2c3dc273, 0x47cf8cac, 0x3b089f07, 0x1e63e807,
        ]);

        let b = Scalar([
            0xd8bedc42, 0x686eb329, 0xe416b899, 0x17aa6d9b, 0x1e30b38b, 0x188c6b1a, 0xd099595b,
            0xbc343bcb, 0x1adaa0e7, 0x24e8d499, 0x8e59b308, 0x0a92de2d, 0xcae1cb68, 0x16c5450a,
        ]);

        let exp = Scalar([
            0xa18d010a, 0x1f5b3197, 0x994c9c2b, 0x6abd26f5, 0x08a3a0e4, 0x36a14920, 0x74e9335f,
            0x07bcd931, 0xf2d89c1e, 0xb9036ff6, 0x203d424b, 0xfccd61b3, 0x4ca389ed, 0x31e055c1,
        ]);

        assert_eq!(a * b, exp)
    }
    #[test]
    fn test_basic_square() {
        let five = Scalar::from(5);
        assert_eq!(five.square(), Scalar::from(25))
    }

    #[test]
    fn test_sanity_check_index_mut() {
        let mut x = Scalar::one();
        x[0] = 2u32;
        assert_eq!(x, Scalar::from(2))
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

    #[test]
    fn test_equals() {
        let a = Scalar::from(5);
        let b = Scalar::from(5);
        let c = Scalar::from(10);
        assert!(a.equals(&b));
        assert!(!a.equals(&c))
    }

    #[test]
    fn test_basic_inversion() {
        // Inversion of one is one
        let one = Scalar::from(1);
        let expected_one = one.invert();
        assert_eq!(expected_one, one);

        // Test inversion from 2 to 100
        for i in 2..=100 {
            let x = Scalar::from(i);
            let x_inv = x.invert();
            assert_eq!(x_inv * x, Scalar::one())
        }

        // Inversion of zero is zero
        let zero = Scalar::zero();
        let expected_zero = zero.invert();
        assert_eq!(expected_zero, zero)
    }
}
