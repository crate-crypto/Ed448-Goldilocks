use super::karatsuba;
use std::ops::{Add, Index, IndexMut, Mul, Sub};
use subtle::{Choice, ConstantTimeEq};
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
impl Mul<&Fq> for Fq {
    type Output = Fq;
    fn mul(self, rhs: &Fq) -> Self::Output {
        karatsuba::mul(&self, rhs)
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
impl Sub<Fq> for Fq {
    type Output = Fq;
    fn sub(self, rhs: Fq) -> Self::Output {
        let mut inter_res = self.sub_no_reduce(&rhs);
        inter_res.bias(2);
        inter_res.weak_reduce();
        inter_res
    }
}
impl ConstantTimeEq for Fq {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0[0].ct_eq(&other.0[0])
            & self.0[1].ct_eq(&other.0[1])
            & self.0[2].ct_eq(&other.0[2])
            & self.0[3].ct_eq(&other.0[3])
            & self.0[4].ct_eq(&other.0[4])
            & self.0[5].ct_eq(&other.0[5])
            & self.0[6].ct_eq(&other.0[6])
            & self.0[7].ct_eq(&other.0[7])
            & self.0[8].ct_eq(&other.0[8])
            & self.0[9].ct_eq(&other.0[9])
            & self.0[10].ct_eq(&other.0[10])
            & self.0[11].ct_eq(&other.0[11])
            & self.0[12].ct_eq(&other.0[12])
            & self.0[13].ct_eq(&other.0[13])
            & self.0[14].ct_eq(&other.0[14])
            & self.0[15].ct_eq(&other.0[15])
    }
}
impl Fq {
    pub(crate) fn zero() -> Fq {
        Fq::from(0u8)
    }
    fn is_zero(&self) -> Choice {
        self.ct_eq(&Fq::zero())
    }
    fn invert(&self) -> Fq {
        let mut t1 = self.square();
        let (mut t2, _) = t1.inverse_square_root();
        t1 = t2.square();
        t2 = t1 * self;
        t2
    }
    // TODO: Implement karatsuba square
    pub(crate) fn square(&self) -> Fq {
        karatsuba::mul(self, self)
    }
    // Squares self n times
    fn square_n(&self, mut n: u32) -> Fq {
        let mut result = self.square();
        n = n - 1;
        for _ in 0..n {
            result = result.square();
        }

        result
    }
    fn inverse_square_root(&self) -> (Fq, bool) {
        let (mut l0, mut l1, mut l2) = (Fq::zero(), Fq::zero(), Fq::zero());

        l1 = self.square();
        l2 = l1 * self;
        l1 = l2.square();
        l2 = l1 * self;
        l1 = l2.square_n(3);
        l0 = l2 * l1;
        l1 = l0.square_n(3);
        l0 = l2 * l1;
        l2 = l0.square_n(9);
        l1 = l0 * l2;
        l0 = l1 * l1;
        l2 = l0 * self;
        l0 = l2.square_n(18);
        l2 = l1 * l0;
        l0 = l2.square_n(37);
        l1 = l2 * l0;
        l0 = l1.square_n(37);
        l1 = l2 * l0;
        l0 = l1.square_n(111);
        l2 = l1 * l0;
        l0 = l2.square();
        l1 = l0 * self;
        l0 = l1.square_n(223);
        l1 = l2 * l0;
        l2 = l1.square();
        l0 = l2 * self;

        let is_residue = l0.equals(&Fq::one());
        (l1, is_residue)
    }

    pub(crate) fn one() -> Fq {
        Fq::from(1u8)
    }

    pub(crate) fn negate(&self) -> Fq {
        Fq::zero() - *self
    }
    /// Bias adds 'b' multiples of `p` to self
    pub(crate) fn bias(&mut self, b: u32) {
        const MASK: u32 = (1 << 28) - 1;

        let co1 = b * MASK;
        let co2 = co1 - b;

        let lo = [co1; 4];
        let hi = [co2, co1, co1, co1];

        self[0] = self[0].wrapping_add(lo[0]);
        self[1] = self[1].wrapping_add(lo[1]);
        self[2] = self[2].wrapping_add(lo[2]);
        self[3] = self[3].wrapping_add(lo[3]);

        self[4] = self[4].wrapping_add(lo[0]);
        self[5] = self[5].wrapping_add(lo[1]);
        self[6] = self[6].wrapping_add(lo[2]);
        self[7] = self[7].wrapping_add(lo[3]);

        self[8] = self[8].wrapping_add(hi[0]);
        self[9] = self[9].wrapping_add(hi[1]);
        self[10] = self[10].wrapping_add(hi[2]);
        self[11] = self[11].wrapping_add(hi[3]);

        self[12] = self[12].wrapping_add(lo[0]);
        self[13] = self[13].wrapping_add(lo[1]);
        self[14] = self[14].wrapping_add(lo[2]);
        self[15] = self[15].wrapping_add(lo[3]);
    }

    pub(crate) fn weak_reduce(&mut self) {
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
    pub(crate) fn add_no_reduce(&self, rhs: &Fq) -> Fq {
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
    pub(crate) fn sub_no_reduce(&self, rhs: &Fq) -> Fq {
        let mut result = Fq::zero();
        result[0] = self[0].wrapping_sub(rhs[0]);
        result[1] = self[1].wrapping_sub(rhs[1]);
        result[2] = self[2].wrapping_sub(rhs[2]);
        result[3] = self[3].wrapping_sub(rhs[3]);
        result[4] = self[4].wrapping_sub(rhs[4]);
        result[5] = self[5].wrapping_sub(rhs[5]);
        result[6] = self[6].wrapping_sub(rhs[6]);
        result[7] = self[7].wrapping_sub(rhs[7]);
        result[8] = self[8].wrapping_sub(rhs[8]);
        result[9] = self[9].wrapping_sub(rhs[9]);
        result[10] = self[10].wrapping_sub(rhs[10]);
        result[11] = self[11].wrapping_sub(rhs[11]);
        result[12] = self[12].wrapping_sub(rhs[12]);
        result[13] = self[13].wrapping_sub(rhs[13]);
        result[14] = self[14].wrapping_sub(rhs[14]);
        result[15] = self[15].wrapping_sub(rhs[15]);
        result
    }

    pub fn equals(&self, rhs: &Fq) -> bool {
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

    pub(crate) fn mul_double_limb(x: &Fq, w: u64) -> Fq {
        let radix = 28;
        let radix_mask = 0xfffffff as u64;
        let mut result = Fq::zero();
        let whi = (w >> radix) as u32;
        let wlo = (w & (radix_mask)) as u32;
        let mut accum0: u64;
        let mut accum8: u64;

        accum0 = (wlo as u64) * (x[0] as u64);
        accum8 = (wlo as u64) * (x[8] as u64);
        accum0 += (whi as u64) * (x[15] as u64);
        accum8 += (whi as u64) * ((x[15] + x[7]) as u64);
        result[0] = (accum0 & radix_mask) as u32;
        accum0 >>= radix;
        result[8] = (accum8 & (radix_mask)) as u32;
        accum8 >>= radix;
        // 1
        accum0 += (wlo as u64) * (x[1] as u64);
        accum8 += (wlo as u64) * (x[9] as u64);
        accum0 += (whi as u64) * (x[0] as u64);
        accum8 += (whi as u64) * (x[8] as u64);
        result[1] = (accum0 & (radix_mask)) as u32;
        accum0 >>= radix;
        result[9] = (accum8 & (radix_mask)) as u32;
        accum8 >>= radix;
        // 2
        accum0 += (wlo as u64) * (x[2] as u64);
        accum8 += (wlo as u64) * (x[10] as u64);
        accum0 += (whi as u64) * (x[1] as u64);
        accum8 += (whi as u64) * (x[9] as u64);
        result[2] = (accum0 & (radix_mask)) as u32;
        accum0 >>= radix;
        result[10] = (accum8 & (radix_mask)) as u32;
        accum8 >>= radix;
        // 3
        accum0 += (wlo as u64) * (x[3] as u64);
        accum8 += (wlo as u64) * (x[11] as u64);
        accum0 += (whi as u64) * (x[2] as u64);
        accum8 += (whi as u64) * (x[10] as u64);
        result[3] = (accum0 & (radix_mask)) as u32;
        accum0 >>= radix;
        result[11] = (accum8 & (radix_mask)) as u32;
        accum8 >>= radix;
        // 4
        accum0 += (wlo as u64) * (x[4] as u64);
        accum8 += (wlo as u64) * (x[12] as u64);
        accum0 += (whi as u64) * (x[3] as u64);
        accum8 += (whi as u64) * (x[11] as u64);
        result[4] = (accum0 & (radix_mask)) as u32;
        accum0 >>= radix;
        result[12] = (accum8 & (radix_mask)) as u32;
        accum8 >>= radix;
        // 5
        accum0 += (wlo as u64) * (x[5] as u64);
        accum8 += (wlo as u64) * (x[13] as u64);
        accum0 += (whi as u64) * (x[4] as u64);
        accum8 += (whi as u64) * (x[12] as u64);
        result[5] = (accum0 & (radix_mask)) as u32;
        accum0 >>= radix;
        result[13] = (accum8 & (radix_mask)) as u32;
        accum8 >>= radix;
        // 6
        accum0 += (wlo as u64) * (x[6] as u64);
        accum8 += (wlo as u64) * (x[14] as u64);
        accum0 += (whi as u64) * (x[5] as u64);
        accum8 += (whi as u64) * (x[13] as u64);
        result[6] = (accum0 & (radix_mask)) as u32;
        accum0 >>= radix;
        result[14] = (accum8 & (radix_mask)) as u32;
        accum8 >>= radix;
        // 7
        accum0 += (wlo as u64) * (x[7] as u64);
        accum8 += (wlo as u64) * (x[15] as u64);
        accum0 += (whi as u64) * (x[6] as u64);
        accum8 += (whi as u64) * (x[14] as u64);
        result[7] = (accum0 & (radix_mask)) as u32;
        accum0 >>= radix;
        result[15] = (accum8 & (radix_mask)) as u32;
        accum8 >>= radix;

        // finish
        accum0 += accum8 + (result[8] as u64);
        result[8] = (accum0 & (radix_mask)) as u32;
        result[9] += (accum0 >> radix) as u32;

        accum8 += result[0] as u64;
        result[0] = (accum8 & (radix_mask)) as u32;
        result[1] += (accum8 >> radix) as u32;
        result
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

        const MASK: u32 = (1 << 28) - 1;
        let mut res = Fq::zero();
        for i in 0..8 {
            // Load i'th 56 bytes
            let out = load7(&bytes[i * 7..]);
            // Process two 28-bit limbs
            res[2 * i] = (out as u32) & MASK;
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

#[cfg(test)]
mod test {
    use super::*;
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

    #[test]
    fn test_isr() {
        let bytes: [u8; 56] = [
            0x9f, 0x93, 0xed, 0x0a, 0x84, 0xde, 0xf0, 0xc7, 0xa0, 0x4b, 0x3f, 0x03, 0x70, 0xc1,
            0x96, 0x3d, 0xc6, 0x94, 0x2d, 0x93, 0xf3, 0xaa, 0x7e, 0x14, 0x96, 0xfa, 0xec, 0x9c,
            0x70, 0xd0, 0x59, 0x3c, 0x5c, 0x06, 0x5f, 0x24, 0x33, 0xf7, 0xad, 0x26, 0x6a, 0x3a,
            0x45, 0x98, 0x60, 0xf4, 0xaf, 0x4f, 0x1b, 0xff, 0x92, 0x26, 0xea, 0xa0, 0x7e, 0x29,
        ];
        let a = Fq::from_bytes(&bytes);
        let (k, _) = a.inverse_square_root();
        let should_be: [u8; 56] = [
            0x2, 0x4f, 0xc6, 0xd1, 0x17, 0x97, 0x83, 0x70, 0xb1, 0x1d, 0xf0, 0x72, 0xee, 0x70,
            0x7e, 0xbf, 0xd0, 0xdf, 0x77, 0xa6, 0x32, 0x28, 0x6f, 0xc6, 0xab, 0xc4, 0x99, 0xb6,
            0xc2, 0x9d, 0xee, 0x6d, 0xe5, 0x76, 0x90, 0xa0, 0x3, 0x82, 0x26, 0x6, 0x34, 0x4a, 0x2a,
            0xb0, 0x47, 0x42, 0xdf, 0x2f, 0x5, 0xbe, 0x4b, 0xa3, 0x13, 0x7d, 0x2, 0x4,
        ];

        for i in 0..56 {
            assert_eq!(should_be[i], k.to_bytes()[i])
        }
    }

    #[test]
    fn test_square_n() {
        let bytes: [u8; 56] = [
            0x9f, 0x93, 0xed, 0x0a, 0x84, 0xde, 0xf0, 0xc7, 0xa0, 0x4b, 0x3f, 0x03, 0x70, 0xc1,
            0x96, 0x3d, 0xc6, 0x94, 0x2d, 0x93, 0xf3, 0xaa, 0x7e, 0x14, 0x96, 0xfa, 0xec, 0x9c,
            0x70, 0xd0, 0x59, 0x3c, 0x5c, 0x06, 0x5f, 0x24, 0x33, 0xf7, 0xad, 0x26, 0x6a, 0x3a,
            0x45, 0x98, 0x60, 0xf4, 0xaf, 0x4f, 0x1b, 0xff, 0x92, 0x26, 0xea, 0xa0, 0x7e, 0x29,
        ];
        let x = Fq::from_bytes(&bytes);
        let y = x.square_n(200);

        let should_be: [u8; 56] = [
            0x75, 0x3e, 0x43, 0xee, 0xb7, 0xc3, 0x1a, 0x28, 0x50, 0x7f, 0xc4, 0x3, 0xe4, 0xba,
            0x79, 0x3a, 0xeb, 0x3b, 0x7f, 0xb3, 0xfc, 0x62, 0xc1, 0x97, 0x1d, 0xec, 0x4b, 0x91,
            0x30, 0x98, 0xa8, 0xe0, 0xa5, 0xca, 0xed, 0xbb, 0x24, 0x79, 0x16, 0x2, 0xeb, 0xa6,
            0xb1, 0x8c, 0xcc, 0x96, 0x70, 0xe2, 0x78, 0xf1, 0xef, 0xf5, 0x56, 0x4b, 0x31, 0x6e,
        ];

        for i in 0..56 {
            assert_eq!(should_be[i], y.to_bytes()[i])
        }
    }

    #[test]
    fn test_basic_mul_double_limb() {
        let c = Fq::mul_double_limb(&Fq::from(9u8), 100);

        assert_eq!(c, Fq::from(900u16));
    }

    #[test]
    fn test_d_min_one() {
        let d = Fq([
            268396374, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
            268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
        ]);

        let mut d_min_one = Fq::from(2u8) * (d - Fq::one());
        d_min_one.strong_reduce();
        dbg!(d_min_one);
    }

    #[test]
    fn test_is_zero() {
        let a = Fq::from(0u8);
        let b = Fq::from(0u16);
        let c = Fq::from(0u32);
        let d = Fq::from(1u32);
        let e = Fq([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 268435455]);

        assert_eq!(a.is_zero().unwrap_u8(), 1u8);
        assert_eq!(b.is_zero().unwrap_u8(), 1u8);
        assert_eq!(c.is_zero().unwrap_u8(), 1u8);
        assert_eq!(d.is_zero().unwrap_u8(), 0u8);
        assert_eq!(e.is_zero().unwrap_u8(), 0u8);
    }
}
