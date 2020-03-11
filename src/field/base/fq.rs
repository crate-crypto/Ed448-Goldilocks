use super::karatsuba;
use std::ops::{Add, Index, IndexMut, Mul};
/// Fq represents an element in the field
/// q = 2^448 - 2^224 -1
///
/// Fq is represented using radix 2^28 as 16 u32s. We therefore represent
/// a field element `x` as x_0 * 2^{28 * 0} + x_1 * 2^{28 * 1} + .... + x_15 * 2^{28 * 15}
///
/// XXX: Compute the wiggle room
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Fq(pub(crate) [u32; 16]);

impl From<u32> for Fq {
    fn from(a: u32) -> Fq {
        Fq([a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
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
    type Output = u32;
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
    fn is_zero() -> bool {
        todo!()
    }
    fn is_residue() -> bool {
        todo!()
    }
    pub(crate) fn one() -> Fq {
        Fq::from(1u8)
    }
    /// Bias adds a multiple of `p` to self
    fn bias(&mut self, b: u32) {
        const MASK: u32 = (1 << 28) - 1;

        let co1 = b * MASK;
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

        let limb16_mod = self[16 - 1] >> 28;
        self[16 / 2] += limb16_mod;

        self[15] = (self[15] & MASK) + (self[14] >> 28);
        self[14] = (self[14] & MASK) + (self[13] >> 28);
        self[13] = (self[13] & MASK) + (self[12] >> 28);
        self[12] = (self[12] & MASK) + (self[11] >> 28);
        self[11] = (self[11] & MASK) + (self[10] >> 28);
        self[10] = (self[10] & MASK) + (self[9] >> 28);
        self[9] = (self[9] & MASK) + (self[8] >> 28);
        self[8] = (self[8] & MASK) + (self[7] >> 28);
        self[7] = (self[7] & MASK) + (self[6] >> 28);
        self[6] = (self[6] & MASK) + (self[5] >> 28);
        self[5] = (self[5] & MASK) + (self[4] >> 28);
        self[4] = (self[4] & MASK) + (self[3] >> 28);
        self[3] = (self[3] & MASK) + (self[2] >> 28);
        self[2] = (self[2] & MASK) + (self[1] >> 28);
        self[1] = (self[1] & MASK) + (self[0] >> 28);

        self[0] = (self[0] & MASK) + limb16_mod;
    }

    fn strong_reduce(&mut self) {
        const MASK: u32 = (1 << 28) - 1;

        // After weak reducing, we know can make the guarantee that
        // 0 < self < 2p
        self.weak_reduce();

        fn sub_carry(a: u32, b: i64) -> i64 {
            let x = a as i64;
            x - b
        }
        fn add_carry(a: u32, b: u64) -> u64 {
            let x = a as u64;
            x + b
        }

        // We know that  0 <= self < 2p
        // This next section computes self - p
        // Our range for self would become -p <= self < p

        let mut scarry = 0i64;
        scarry += sub_carry(self[0], 0xfffffff);
        self[0] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[1], 0xfffffff);
        self[1] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[2], 0xfffffff);
        self[2] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[3], 0xfffffff);
        self[3] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[4], 0xfffffff);
        self[4] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[5], 0xfffffff);
        self[5] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[6], 0xfffffff);
        self[6] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[7], 0xfffffff);
        self[7] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[8], 0xffffffe);
        self[8] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[9], 0xfffffff);
        self[9] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[10], 0xfffffff);
        self[10] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[11], 0xfffffff);
        self[11] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[12], 0xfffffff);
        self[12] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[13], 0xfffffff);
        self[13] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[14], 0xfffffff);
        self[14] = (scarry as u32) & MASK;
        scarry >>= 28;
        scarry += sub_carry(self[15], 0xfffffff);
        self[15] = (scarry as u32) & MASK;
        scarry >>= 28;

        // There are two cases to consider; either the value was >= p or it was <less than> p
        // Case 1:
        // If the value was more than p, then we will not have a borrow bit, when we subtracted the modulus
        // In this case, scarry = 0
        // Case 2:
        // The other case, is that our value was less than p.
        // We will therefore have a borrow bit of - 1 signifying that our total is negative.

        // if borrow is not -1 or 0, then we have a problem
        assert!(scarry == 0 || scarry + 1 == 0);

        let scarry_mask = (scarry as u32) & MASK;
        let mut carry = 0u64;
        let m = scarry_mask as u64;

        carry += add_carry(self[0], m);
        self[0] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[1], m);
        self[1] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[2], m);
        self[2] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[3], m);
        self[3] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[4], m);
        self[4] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[5], m);
        self[5] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[6], m);
        self[6] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[7], m);
        self[7] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[8], m & 0xfffffffffffffffe);
        self[8] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[9], m);
        self[9] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[10], m);
        self[10] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[11], m);
        self[11] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[12], m);
        self[12] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[13], m);
        self[13] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[14], m);
        self[14] = (carry as u32) & MASK;
        carry >>= 28;
        carry += add_carry(self[15], m);
        self[15] = (carry as u32) & MASK;

        // This last shift is only needed for the asserts that follow
        carry >>= 28;
        // carry can only be 1 or 0 depending on the previous borrow bit
        assert!(carry < 2);
        // carry + borrow should cancel out
        assert!((carry as i64) + (scarry as i64) == 0);
    }

    // Adds the two field elements together
    // Result is not reduced
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
    // Subtracts the two field elements from each other
    // Result is not reduced
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

    fn equals(&self, rhs: &Fq) -> bool {
        // Subtract self from rhs
        let mut difference = self.sub_no_reduce(rhs);
        difference.bias(2);
        difference.strong_reduce();

        let mut r = 0u32;

        r |= difference[0];
        r |= difference[1];
        r |= difference[2];
        r |= difference[3];
        r |= difference[4];
        r |= difference[5];
        r |= difference[6];
        r |= difference[7];
        r |= difference[8];
        r |= difference[9];
        r |= difference[10];
        r |= difference[11];
        r |= difference[12];
        r |= difference[13];
        r |= difference[14];

        word_is_zero(r)
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
            res[2 * i] = (out & MASK) as u32;
            res[2 * i + 1] = (out >> 28) as u32;
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
            let mut l = (limbs[2 * i] as u64) + ((limbs[2 * i + 1] as u64) << 28);

            for j in 0..7 {
                res[7 * i + j] = l as u8;
                l >>= 8;
            }
        }
        res
    }
}

fn word_is_zero(word: u32) -> bool {
    // If the word is zero, then when we minus 1, we should get 0xffffffff
    word.wrapping_sub(1) == 0xffffffff
}

#[test]
fn test_add() {
    let a = Fq::from(8u8);
    let b = a + a;
    let c = Fq::from(16u8);
    assert!(b.equals(&c));
}
#[test]
fn test_bias() {
    let mut a = Fq::from(5u8);
    a.bias(2);
    let b = Fq::from(5u8);
    assert!(a.equals(&b));
}

#[test]
#[should_panic]
fn test_bias_more_than_headroom() {
    let mut a = Fq::from(5u8);
    a.bias(17);
}

#[test]
fn test_sub() {
    let x = Fq::from(255u8);
    let y = Fq::from(255u8);
    let MODULUS = Fq([
        0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff,
        0xffffffe, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff,
    ]);
    // First subtract without reducing
    let mut z = x.sub_no_reduce(&y);
    // Then bias z by 2
    z.bias(2);
    // Then clear high bits
    z.weak_reduce();

    assert!(z.equals(&MODULUS));
}

#[test]
fn test_bytes_function() {
    let bytes = [1; 56];
    let a = Fq::from_bytes(&bytes);
    let new_a = Fq::from_bytes(&a.to_bytes());
    assert!(a.equals(&new_a));
}
