use super::karatsuba;
use super::limb28::Limb28;
use std::ops::{Add, Index, IndexMut, Mul};
/// Fq is a field element defined over the Goldilocks prime
/// q = 2^448 - 2^224 -1
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Fq(pub(crate) [Limb28; 16]);

impl From<u32> for Fq {
    fn from(a: u32) -> Fq {
        Fq([
            Limb28::from(a),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
            Limb28::zero(),
        ])
    }
}
impl From<u16> for Fq {
    fn from(a: u16) -> Fq {
        Fq::from(a as u32)
    }
}
impl From<u8> for Fq {
    fn from(a: u8) -> Fq {
        Fq::from(a as u32)
    }
}

impl Index<usize> for Fq {
    type Output = Limb28;
    fn index(&self, a: usize) -> &Self::Output {
        &self.0[a]
    }
}
impl IndexMut<usize> for Fq {
    fn index_mut(&mut self, a: usize) -> &mut Self::Output {
        &mut self.0[a]
    }
}
impl Mul<Fq> for Fq {
    type Output = Fq;
    fn mul(self, rhs: Fq) -> Self::Output {
        karatsuba::mul(&self, &rhs)
    }
}

impl Add<Fq> for Fq {
    type Output = Fq;
    fn add(self, rhs: Fq) -> Self::Output {
        let mut inter_res = self.add_no_reduce(&rhs);
        inter_res.weak_reduce();
        inter_res
    }
}
impl Fq {
    pub(crate) fn zero() -> Fq {
        Fq::from(0u8)
    }
    pub(crate) fn one() -> Fq {
        Fq::from(1u8)
    }
    fn bias(&mut self, b: Limb28) {
        const MASK: u32 = (1 << 28) - 1;

        let co1 = Limb28::from_checked_u64(b * MASK);
        let co2 = co1 - b;

        let lo = [co1; 4];
        let hi = [co2, co1, co1, co1];

        self[0] += lo[0];
        self[1] += lo[1];
        self[2] += lo[2];
        self[3] += lo[3];
        self[4] += lo[0];
        self[5] += lo[1];
        self[6] += lo[2];
        self[7] += lo[3];
        self[8] += hi[0];
        self[9] += hi[1];
        self[10] += hi[2];
        self[11] += hi[3];
        self[12] += lo[0];
        self[13] += lo[1];
        self[14] += lo[2];
        self[15] += lo[3];
    }

    fn weak_reduce(&mut self) {
        const MASK: u32 = (1 << 28) - 1;

        let nLimbs = 16;
        let radix = 28;

        let limb16_mod = self[nLimbs - 1] >> radix;
        self[nLimbs / 2] += limb16_mod;

        self[15] = (self[15] & MASK) + (self[14] >> radix);
        self[14] = (self[14] & MASK) + (self[13] >> radix);
        self[13] = (self[13] & MASK) + (self[12] >> radix);
        self[12] = (self[12] & MASK) + (self[11] >> radix);
        self[11] = (self[11] & MASK) + (self[10] >> radix);
        self[10] = (self[10] & MASK) + (self[9] >> radix);
        self[9] = (self[9] & MASK) + (self[8] >> radix);
        self[8] = (self[8] & MASK) + (self[7] >> radix);
        self[7] = (self[7] & MASK) + (self[6] >> radix);
        self[6] = (self[6] & MASK) + (self[5] >> radix);
        self[5] = (self[5] & MASK) + (self[4] >> radix);
        self[4] = (self[4] & MASK) + (self[3] >> radix);
        self[3] = (self[3] & MASK) + (self[2] >> radix);
        self[2] = (self[2] & MASK) + (self[1] >> radix);
        self[1] = (self[1] & MASK) + (self[0] >> radix);

        self[0] = (self[0] & MASK) + limb16_mod;
    }

    fn strong_reduce(&mut self) {
        const MASK: u32 = (1 << 28) - 1;
        let MODULUS = Fq([
            Limb28::from(0xffffff),
            Limb28::from(0xffffffff),
            Limb28::from(0xffffff),
            Limb28::from(0xffffffff),
            Limb28::from(0xffffff),
            Limb28::from(0xffffffff),
            Limb28::from(0xffffff),
            Limb28::from(0xffffffff),
            Limb28::from(0xfffffe),
            Limb28::from(0xffffffff),
            Limb28::from(0xffffff),
            Limb28::from(0xffffffff),
            Limb28::from(0xffffff),
            Limb28::from(0xffffffff),
            Limb28::from(0xffffff),
            Limb28::from(0xffffffff),
        ]);

        // clear high
        self.weak_reduce();

        fn sub_carry(a: Limb28, b: i64) -> i64 {
            let x = a.value() as i64;
            x - b
        }
        fn add_carry(a: Limb28, b: u64) -> u64 {
            let x = a.value() as u64;
            x + b
        }

        // total is less than 2p
        // compute total_value - p.  No need to reduce mod p.

        let mut scarry = 0i64;
        scarry += sub_carry(self[0], 0xfffffff);
        self[0] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[1], 0xfffffff);
        self[1] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[2], 0xfffffff);
        self[2] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[3], 0xfffffff);
        self[3] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[4], 0xfffffff);
        self[4] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[5], 0xfffffff);
        self[5] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[6], 0xfffffff);
        self[6] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[7], 0xfffffff);
        self[7] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[8], 0xffffffe);
        self[8] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[9], 0xfffffff);
        self[9] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[10], 0xfffffff);
        self[10] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[11], 0xfffffff);
        self[11] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[12], 0xfffffff);
        self[12] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[13], 0xfffffff);
        self[13] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[14], 0xfffffff);
        self[14] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        scarry += sub_carry(self[15], 0xfffffff);
        self[15] = Limb28::from_checked_i64(scarry) & MASK;
        scarry >>= 28;

        // uncommon case: it was >= p, so now scarry = 0 and this = x
        // common case: it was < p, so now scarry = -1 and this = x - p + 2^255
        // so let's add back in p.  will carry back off the top for 2^255.
        // it can be asserted:
        assert!(scarry == 0 || scarry + 1 == 0);

        let scarryMask = (scarry as u32) & MASK;
        let mut carry = 0u64;
        let m = scarryMask as u64;

        carry += add_carry(self[0], m);
        self[0] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[1], m);
        self[1] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[2], m);
        self[2] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[3], m);
        self[3] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[4], m);
        self[4] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[5], m);
        self[5] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[6], m);
        self[6] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[7], m);
        self[7] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[8], m & 0xfffffffffffffffe);
        self[8] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[9], m);
        self[9] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[10], m);
        self[10] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[11], m);
        self[11] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[12], m);
        self[12] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[13], m);
        self[13] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[14], m);
        self[14] = Limb28::from_checked_u64(carry) & MASK;
        carry >>= 28;

        carry += add_carry(self[15], m);
        self[15] = Limb28::from_checked_u64(carry) & MASK;

        // assert!(carry < 2 && (carry as u32) + (scarry as u32) == 0);
    }

    fn add_no_reduce(&self, rhs: &Fq) -> Fq {
        let mut result = Fq::zero();

        result[0] = self[0] + rhs[0];
        result[1] = self[1] + rhs[1];
        result[2] = self[2] + rhs[2];
        result[3] = self[3] + rhs[3];
        result[4] = self[4] + rhs[4];
        result[5] = self[5] + rhs[5];
        result[6] = self[6] + rhs[6];
        result[7] = self[7] + rhs[7];
        result[8] = self[8] + rhs[8];
        result[9] = self[9] + rhs[9];
        result[10] = self[10] + rhs[10];
        result[11] = self[11] + rhs[11];
        result[12] = self[12] + rhs[12];
        result[13] = self[13] + rhs[13];
        result[14] = self[14] + rhs[14];
        result[15] = self[15] + rhs[15];

        result
    }
    fn sub_no_reduce(&self, rhs: &Fq) -> Fq {
        let mut result = Fq::zero();
        result[0] = self[0] - rhs[0];
        result[1] = self[1] - rhs[1];
        result[2] = self[2] - rhs[2];
        result[3] = self[3] - rhs[3];
        result[4] = self[4] - rhs[4];
        result[5] = self[5] - rhs[5];
        result[6] = self[6] - rhs[6];
        result[7] = self[7] - rhs[7];
        result[8] = self[8] - rhs[8];
        result[9] = self[9] - rhs[9];
        result[10] = self[10] - rhs[10];
        result[11] = self[11] - rhs[11];
        result[12] = self[12] - rhs[12];
        result[13] = self[13] - rhs[13];
        result[14] = self[14] - rhs[14];
        result[15] = self[15] - rhs[15];
        result
    }
    // temporary equal function, ideally would like to avoid cloning
    fn equals(&self, rhs: &Fq) -> bool {
        let mut r = 0u32;
        let mut x = self.clone();
        x.strong_reduce();
        let mut y = rhs.clone();
        y.strong_reduce();

        r |= x[0] ^ y[0];
        r |= x[1] ^ y[1];
        r |= x[2] ^ y[2];
        r |= x[3] ^ y[3];
        r |= x[4] ^ y[4];
        r |= x[5] ^ y[5];
        r |= x[6] ^ y[6];
        r |= x[7] ^ y[7];
        r |= x[8] ^ y[8];
        r |= x[9] ^ y[9];
        r |= x[10] ^ y[10];
        r |= x[11] ^ y[11];
        r |= x[12] ^ y[12];
        r |= x[13] ^ y[13];
        r |= x[14] ^ y[14];
        r |= x[15] ^ y[15];

        r == 0
    }

    // Currently this does not check if the encoding is canonical (ie if the Field number is reduced)
    // We will parse in chunks of 56 bytes
    // The first 28 bytes will contain the i'th limb
    // The second 28 bytes will contain the (2i+1)'th limb
    pub(crate) fn from_bytes(bytes: &[u8; 56]) -> Fq {
        let load7 = |input: &[u8]| -> u64 {
            (input[0] as u64)
                | ((input[1] as u64) << 8)
                | ((input[2] as u64) << 16)
                | ((input[3] as u64) << 24)
                | ((input[4] as u64) << 32)
                | ((input[5] as u64) << 40)
                | ((input[6] as u64) << 48)
        };

        const MASK: u64 = (1 << 28) - 1;
        let mut res = Fq::zero();
        for i in 0..8 {
            // Load i'th 56 bytes
            let out = load7(&bytes[i + 7..]);
            // Process two 28-bit limbs
            res[2 * i] = Limb28::from_checked_u64(out & MASK);
            res[2 * i + 1] = Limb28::from_checked_u64(out >> 28);
        }

        res
    }

    // We encode the Field element by storing each consecutive into a u64
    pub(crate) fn to_bytes(&self) -> [u8; 56] {
        // First we reduce the element so that to_bytes always produces a canonical encoding
        let mut limbs = self.clone();
        limbs.strong_reduce();

        let mut res = [0u8; 56];

        for i in 0..8 {
            let mut l = (limbs[2 * i].0 as u64) + ((limbs[2 * i + 1].0 as u64) << 28);

            for j in 0..7 {
                res[7 * i + j] = l as u8;
                l >>= 8;
            }
        }
        res
    }
}

#[test]
fn test_add() {
    let a = Fq::from(8u8);
    let b = a + a;
    let mut c = Fq::from(16u8);
    assert!(b.equals(&c));
}
#[test]
fn test_add_bias() {
    let mut a = Fq::from(5u8);
    let mut b = (a + a);
    b.bias(Limb28::from(2u32));
    let mut c = Fq::from(10u8);
    assert!(b.equals(&c));
}

#[test]
fn test_bytes_function() {
    let bytes = [1; 56];
    let a = Fq::from_bytes(&bytes);
    let new_a = Fq::from_bytes(&a.to_bytes());
    assert!(a.equals(&new_a));
}
