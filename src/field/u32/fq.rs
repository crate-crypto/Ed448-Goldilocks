use super::karatsuba;
use std::ops::{Add, Index, IndexMut, Mul, Sub};
use subtle::{Choice, ConditionallyNegatable, ConditionallySelectable, ConstantTimeEq};
/// Fq represents an element in the field
/// q = 2^448 - 2^224 -1
///
/// Fq is represented using radix 2^28 as 16 u32s. We therefore represent
/// a field element `x` as x_0 + x_1 * 2^{28 * 1} + .... + x_15 * 2^{28 * 15}

#[derive(Copy, Clone, Debug)]
pub struct Fq(pub(crate) [u32; 16]);

////
/// Trait Implementations
///

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
    fn sub(mut self, rhs: Fq) -> Self::Output {
        self = Fq::bias(&self, 2);
        let mut inter_res = self.sub_no_reduce(&rhs);
        inter_res.weak_reduce();
        inter_res
    }
}

impl ConditionallyNegatable for Fq {
    fn conditional_negate(&mut self, choice: Choice) {
        let self_neg = self.clone().negate();
        self.conditional_assign(&self_neg, choice);
    }
}

impl ConditionallySelectable for Fq {
    fn conditional_select(a: &Fq, b: &Fq, choice: Choice) -> Fq {
        Fq([
            u32::conditional_select(&a.0[0], &b.0[0], choice),
            u32::conditional_select(&a.0[1], &b.0[1], choice),
            u32::conditional_select(&a.0[2], &b.0[2], choice),
            u32::conditional_select(&a.0[3], &b.0[3], choice),
            u32::conditional_select(&a.0[4], &b.0[4], choice),
            u32::conditional_select(&a.0[5], &b.0[5], choice),
            u32::conditional_select(&a.0[6], &b.0[6], choice),
            u32::conditional_select(&a.0[7], &b.0[7], choice),
            u32::conditional_select(&a.0[8], &b.0[8], choice),
            u32::conditional_select(&a.0[9], &b.0[9], choice),
            u32::conditional_select(&a.0[10], &b.0[10], choice),
            u32::conditional_select(&a.0[11], &b.0[11], choice),
            u32::conditional_select(&a.0[12], &b.0[12], choice),
            u32::conditional_select(&a.0[13], &b.0[13], choice),
            u32::conditional_select(&a.0[14], &b.0[14], choice),
            u32::conditional_select(&a.0[15], &b.0[15], choice),
        ])
    }

    fn conditional_assign(&mut self, other: &Fq, choice: Choice) {
        self.0[0].conditional_assign(&other.0[0], choice);
        self.0[1].conditional_assign(&other.0[1], choice);
        self.0[2].conditional_assign(&other.0[2], choice);
        self.0[3].conditional_assign(&other.0[3], choice);
        self.0[4].conditional_assign(&other.0[4], choice);
        self.0[5].conditional_assign(&other.0[5], choice);
        self.0[6].conditional_assign(&other.0[6], choice);
        self.0[7].conditional_assign(&other.0[7], choice);
        self.0[8].conditional_assign(&other.0[8], choice);
        self.0[9].conditional_assign(&other.0[9], choice);
        self.0[10].conditional_assign(&other.0[10], choice);
        self.0[11].conditional_assign(&other.0[11], choice);
        self.0[12].conditional_assign(&other.0[12], choice);
        self.0[13].conditional_assign(&other.0[13], choice);
        self.0[14].conditional_assign(&other.0[14], choice);
        self.0[15].conditional_assign(&other.0[15], choice);
    }

    fn conditional_swap(a: &mut Fq, b: &mut Fq, choice: Choice) {
        u32::conditional_swap(&mut a.0[0], &mut b.0[0], choice);
        u32::conditional_swap(&mut a.0[1], &mut b.0[1], choice);
        u32::conditional_swap(&mut a.0[2], &mut b.0[2], choice);
        u32::conditional_swap(&mut a.0[3], &mut b.0[3], choice);
        u32::conditional_swap(&mut a.0[4], &mut b.0[4], choice);
        u32::conditional_swap(&mut a.0[5], &mut b.0[5], choice);
        u32::conditional_swap(&mut a.0[6], &mut b.0[6], choice);
        u32::conditional_swap(&mut a.0[7], &mut b.0[7], choice);
        u32::conditional_swap(&mut a.0[8], &mut b.0[8], choice);
        u32::conditional_swap(&mut a.0[9], &mut b.0[9], choice);
        u32::conditional_swap(&mut a.0[10], &mut b.0[10], choice);
        u32::conditional_swap(&mut a.0[11], &mut b.0[11], choice);
        u32::conditional_swap(&mut a.0[12], &mut b.0[12], choice);
        u32::conditional_swap(&mut a.0[13], &mut b.0[13], choice);
        u32::conditional_swap(&mut a.0[14], &mut b.0[14], choice);
        u32::conditional_swap(&mut a.0[15], &mut b.0[15], choice);
    }
}

impl ConstantTimeEq for Fq {
    fn ct_eq(&self, other: &Self) -> Choice {
        let mut difference = *self - *other;
        difference.strong_reduce();

        let zero = Fq::zero();

        difference[0].ct_eq(&zero[0])
            & difference[1].ct_eq(&zero[1])
            & difference[2].ct_eq(&zero[2])
            & difference[3].ct_eq(&zero[3])
            & difference[4].ct_eq(&zero[4])
            & difference[5].ct_eq(&zero[5])
            & difference[6].ct_eq(&zero[6])
            & difference[7].ct_eq(&zero[7])
            & difference[8].ct_eq(&zero[8])
            & difference[9].ct_eq(&zero[9])
            & difference[10].ct_eq(&zero[10])
            & difference[11].ct_eq(&zero[11])
            & difference[12].ct_eq(&zero[12])
            & difference[13].ct_eq(&zero[13])
            & difference[14].ct_eq(&zero[14])
            & difference[15].ct_eq(&zero[15])
    }
}

impl PartialEq for Fq {
    fn eq(&self, other: &Fq) -> bool {
        self.ct_eq(&other).into()
    }
}
impl Eq for Fq {}

impl Default for Fq {
    fn default() -> Fq {
        Fq::zero()
    }
}

///
/// Constants
///

impl Fq {
    pub const fn zero() -> Fq {
        Fq([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
    pub const fn one() -> Fq {
        Fq([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
    pub fn minus_one() -> Fq {
        Fq([
            268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
            268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
        ])
    }
}

///
/// Checks
///
impl Fq {
    pub fn is_zero(&self) -> Choice {
        self.ct_eq(&Fq::zero())
    }
    pub fn is_negative(&self) -> Choice {
        let bytes = self.to_bytes();
        (bytes[0] & 1).into()
    }
}

///
/// Serialisation
///
impl Fq {
    /// Helper function for internally constructing a field element
    pub(crate) const fn from_raw_slice(slice: [u32; 16]) -> Fq {
        Fq(slice)
    }

    /// This does not check if the encoding is canonical (ie if the input is reduced)
    /// We parse in chunks of 56 bytes, the first 28 bytes will contain the i'th limb
    /// and the second 28 bytes will contain the (2i+1)'th limb
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
        // Reduce element to a canonical representation.
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

impl Fq {
    /// Inverts a field element
    pub(crate) fn invert(&self) -> Fq {
        let mut t1 = self.square();
        let (mut t2, _) = t1.inverse_square_root();
        t1 = t2.square();
        t2 = t1 * self;
        t2
    }
    /// Sqaures a field element
    pub(crate) fn square(&self) -> Fq {
        karatsuba::square(self)
    }
    /// Squares a field element  `n` times
    fn square_n(&self, mut n: u32) -> Fq {
        let mut result = self.square();

        // Decrease value by 1 since we just did a squaring
        n = n - 1;

        for _ in 0..n {
            result = result.square();
        }

        result
    }

    /// Computes the inverse square root of a field element
    /// Returns the result and a boolean to indicate whether self
    /// was a Quadratic residue
    pub(crate) fn inverse_square_root(&self) -> (Fq, bool) {
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

        let is_residue = l0 == Fq::one();
        (l1, is_residue)
    }

    /// Computes the square root ratio of two elements
    pub(crate) fn sqrt_ratio(u: &Fq, v: &Fq) -> (Fq, bool) {
        let x = *u * v;
        let (inv_sqrt_x, is_res) = x.inverse_square_root();
        (inv_sqrt_x * u, is_res)
    }

    /// Negates a field element
    pub(crate) fn negate(&self) -> Fq {
        Fq::zero() - *self
    }

    /// Bias adds 'b' multiples of `p` to self
    pub(crate) fn bias(a: &Fq, b: u32) -> Fq {
        const MASK: u32 = (1 << 28) - 1;

        let co1 = b * MASK;
        let co2 = co1 - b;

        let lo = [co1; 4];
        let hi = [co2, co1, co1, co1];

        Fq([
            a[0].wrapping_add(lo[0]),
            a[1].wrapping_add(lo[1]),
            a[2].wrapping_add(lo[2]),
            a[3].wrapping_add(lo[3]),
            a[4].wrapping_add(lo[0]),
            a[5].wrapping_add(lo[1]),
            a[6].wrapping_add(lo[2]),
            a[7].wrapping_add(lo[3]),
            a[8].wrapping_add(hi[0]),
            a[9].wrapping_add(hi[1]),
            a[10].wrapping_add(hi[2]),
            a[11].wrapping_add(hi[3]),
            a[12].wrapping_add(lo[0]),
            a[13].wrapping_add(lo[1]),
            a[14].wrapping_add(lo[2]),
            a[15].wrapping_add(lo[3]),
        ])
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

    /// Reduces the field element to a canonical representation
    /// This is used when checking equality between two field elements and
    /// when encoding a field element
    pub(crate) fn strong_reduce(&mut self) {
        const MASK: u32 = (1 << 28) - 1;

        // After weak reducing, we can make the guarantee that
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

        // There are two cases to consider; either the value was >= p or it was <less than p
        // Case 1:
        // If the value was more than p, then the final borrow will be zero. This is scarry.
        // Case 2:
        // If the  value was less than p, the final borrow will be -1.

        // The only two possibilities for the borrow bit is -1 or 0.
        assert!(scarry == 0 || scarry + 1 == 0);

        let scarry_mask = (scarry as u32) & MASK;
        let mut carry = 0u64;
        let m = scarry_mask as u64;

        // Carry propagation

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

    /// Adds two field elements together
    /// Result is not reduced
    pub(crate) fn add_no_reduce(&self, rhs: &Fq) -> Fq {
        Fq([
            self[0] + rhs[0],
            self[1] + rhs[1],
            self[2] + rhs[2],
            self[3] + rhs[3],
            self[4] + rhs[4],
            self[5] + rhs[5],
            self[6] + rhs[6],
            self[7] + rhs[7],
            self[8] + rhs[8],
            self[9] + rhs[9],
            self[10] + rhs[10],
            self[11] + rhs[11],
            self[12] + rhs[12],
            self[13] + rhs[13],
            self[14] + rhs[14],
            self[15] + rhs[15],
        ])
    }

    /// Subtracts the two field elements from each other
    /// Result is not reduced
    pub(crate) fn sub_no_reduce(&self, rhs: &Fq) -> Fq {
        Fq([
            self[0].wrapping_sub(rhs[0]),
            self[1].wrapping_sub(rhs[1]),
            self[2].wrapping_sub(rhs[2]),
            self[3].wrapping_sub(rhs[3]),
            self[4].wrapping_sub(rhs[4]),
            self[5].wrapping_sub(rhs[5]),
            self[6].wrapping_sub(rhs[6]),
            self[7].wrapping_sub(rhs[7]),
            self[8].wrapping_sub(rhs[8]),
            self[9].wrapping_sub(rhs[9]),
            self[10].wrapping_sub(rhs[10]),
            self[11].wrapping_sub(rhs[11]),
            self[12].wrapping_sub(rhs[12]),
            self[13].wrapping_sub(rhs[13]),
            self[14].wrapping_sub(rhs[14]),
            self[15].wrapping_sub(rhs[15]),
        ])
    }
}

pub(crate) fn from_u32(a: u32) -> Fq {
    Fq([a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
}
#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_add() {
        let a = from_u32(8);
        let b = a + a;
        let c = from_u32(16);
        assert!(b == c);
    }
    #[test]
    fn test_bias() {
        let mut a = from_u32(5);
        a = Fq::bias(&a, 2);
        let b = from_u32(5);
        assert!(a == b);
    }

    #[test]
    #[should_panic]
    fn test_bias_more_than_headroom() {
        let mut a = from_u32(5);
        a = Fq::bias(&a, 17);
    }

    #[test]
    fn test_equals() {
        let a = from_u32(10);
        let b = from_u32(20);
        assert!(a != b);

        let c = from_u32(99);
        let d = from_u32(95);

        let ab = a * b;
        let ba = b * a;
        let cd = c * d;
        assert!(ab != cd);
        assert!(ab == ba);
    }

    #[test]
    fn test_sub() {
        let x = from_u32(255);
        let y = from_u32(255);
        let MODULUS = Fq([
            0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff,
            0xffffffe, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff,
        ]);
        // First subtract without reducing
        let mut z = x.sub_no_reduce(&y);
        // Then bias z by 2
        z = Fq::bias(&z, 2);
        // Then clear high bits
        z.weak_reduce();

        assert!(z == MODULUS);
    }

    #[test]
    fn test_bytes_function() {
        let bytes = [1; 56];
        let a = Fq::from_bytes(&bytes);
        let new_a = Fq::from_bytes(&a.to_bytes());
        assert!(a == new_a);
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
    fn test_d_min_one() {
        let d = Fq([
            268396374, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
            268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
        ]);

        let mut d_min_one = (Fq::one() + Fq::one()) * (d - Fq::one());
        d_min_one.strong_reduce();
        dbg!(d_min_one);
    }

    #[test]
    fn test_is_zero() {
        let a = from_u32(0);
        let b = from_u32(0);
        let c = from_u32(0);
        let d = from_u32(1);
        let e = Fq([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 268435455]);

        assert_eq!(a.is_zero().unwrap_u8(), 1u8);
        assert_eq!(b.is_zero().unwrap_u8(), 1u8);
        assert_eq!(c.is_zero().unwrap_u8(), 1u8);
        assert_eq!(d.is_zero().unwrap_u8(), 0u8);
        assert_eq!(e.is_zero().unwrap_u8(), 0u8);
    }
    #[test]
    fn test_conditional_swap() {
        let mut x = from_u32(1908);
        let x_old = from_u32(1908);
        let mut y = from_u32(200);
        let y_old = from_u32(200);

        // x and y should swap value
        Fq::conditional_swap(&mut x, &mut y, Choice::from(1));
        assert!(x == y_old);
        assert!(y == x_old);

        // x and y should stay the same
        Fq::conditional_swap(&mut x, &mut y, Choice::from(0));
        assert!(x == y_old);
        assert!(y == x_old);
    }

    #[test]
    fn test_conditional_negate() {
        let mut a = from_u32(100);
        let a_neg = a.negate();
        a.conditional_negate(Choice::from(1));
        assert!(a == a_neg);
        //
        let mut b = from_u32(200);
        b.conditional_negate(Choice::from(0));
        assert!(b == from_u32(200));
    }

    #[test]
    fn test_sqrt_ratio() {
        let ten = from_u32(10);
        let twenty = from_u32(20);
        let (a, is_res) = Fq::sqrt_ratio(&ten, &twenty);
        assert!(is_res);

        let (inv_ten_sqrt, is_res) = ten.inverse_square_root();
        assert!(is_res);
        let (inv_twenty_sqrt, is_res) = twenty.inverse_square_root();
        assert!(is_res);

        let expected = inv_ten_sqrt.invert() * inv_twenty_sqrt;
        assert_eq!(expected, a)
    }
}
