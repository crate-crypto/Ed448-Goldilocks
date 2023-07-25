use fiat_crypto::p448_solinas_64::*;

use std::ops::{Add, Index, IndexMut, Mul, Sub};
use subtle::{Choice, ConditionallyNegatable, ConditionallySelectable};
/// FieldElement56 represents an element in the field
/// q = 2^448 - 2^224 -1
///
/// FieldElement56 is represented using radix 2^56 as 8 u64s. We therefore represent
/// a field element `x` as x_0 + x_1 * 2^{28 * 1} + .... + x_15 * 2^{28 * 15}

// XXX: Check if the Serialisation procedure in FieldElement56 is consistent with FieldElement28

#[derive(Copy, Clone)]
pub struct FieldElement56(pub(crate) fiat_p448_tight_field_element);

impl std::fmt::Debug for FieldElement56 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("FieldElement56").field(&self.0 .0).finish()
    }
}

////
/// Trait Implementations
///

impl Index<usize> for FieldElement56 {
    type Output = u64;
    fn index(&self, a: usize) -> &Self::Output {
        &self.0[a]
    }
}

impl IndexMut<usize> for FieldElement56 {
    fn index_mut(&mut self, a: usize) -> &mut Self::Output {
        &mut self.0[a]
    }
}
impl Mul<&FieldElement56> for &FieldElement56 {
    type Output = FieldElement56;
    fn mul(self, rhs: &FieldElement56) -> Self::Output {
        let mut result = FieldElement56::zero();
        let mut self_loose = fiat_p448_loose_field_element([0; 8]);
        fiat_p448_relax(&mut self_loose, &self.0);
        let mut rhs_loose = fiat_p448_loose_field_element([0; 8]);
        fiat_p448_relax(&mut rhs_loose, &rhs.0);
        fiat_p448_carry_mul(&mut result.0, &self_loose, &rhs_loose);
        result
    }
}
impl Mul<&FieldElement56> for FieldElement56 {
    type Output = FieldElement56;
    fn mul(self, rhs: &FieldElement56) -> Self::Output {
        &self * rhs
    }
}
impl Mul<FieldElement56> for FieldElement56 {
    type Output = FieldElement56;
    fn mul(self, rhs: FieldElement56) -> Self::Output {
        &self * &rhs
    }
}

impl Add<FieldElement56> for FieldElement56 {
    type Output = FieldElement56;
    fn add(self, rhs: FieldElement56) -> Self::Output {
        let mut result_loose = fiat_p448_loose_field_element([0; 8]);
        fiat_p448_add(&mut result_loose, &self.0, &rhs.0);
        let mut result = FieldElement56::zero();
        fiat_p448_carry(&mut result.0, &result_loose);
        result
    }
}
impl Sub<FieldElement56> for FieldElement56 {
    type Output = FieldElement56;
    fn sub(self, rhs: FieldElement56) -> Self::Output {
        let mut result_loose = fiat_p448_loose_field_element([0; 8]);
        fiat_p448_sub(&mut result_loose, &self.0, &rhs.0);
        let mut result = FieldElement56::zero();
        fiat_p448_carry(&mut result.0, &result_loose);
        result
    }
}

impl ConditionallyNegatable for FieldElement56 {
    fn conditional_negate(&mut self, choice: Choice) {
        let self_neg = self.clone().negate();
        self.conditional_assign(&self_neg, choice);
    }
}

impl ConditionallySelectable for FieldElement56 {
    fn conditional_select(
        a: &FieldElement56,
        b: &FieldElement56,
        choice: Choice,
    ) -> FieldElement56 {
        let mut result = FieldElement56::zero();
        fiat_p448_selectznz(&mut (result.0).0, choice.unwrap_u8(), &(a.0).0, &(b.0).0);
        result
    }
}

impl Default for FieldElement56 {
    fn default() -> FieldElement56 {
        FieldElement56::zero()
    }
}

///
/// Constants
///

impl FieldElement56 {
    pub const fn zero() -> FieldElement56 {
        FieldElement56(fiat_p448_tight_field_element([0; 8]))
    }
    pub const fn one() -> FieldElement56 {
        FieldElement56(fiat_p448_tight_field_element([1, 0, 0, 0, 0, 0, 0, 0]))
    }
    pub fn minus_one() -> FieldElement56 {
        FieldElement56(fiat_p448_tight_field_element([
            144115188075855869,
            144115188075855870,
            144115188075855870,
            144115188075855870,
            144115188075855868,
            144115188075855870,
            144115188075855870,
            144115188075855870,
        ]))
    }
}

///
/// Checks
///
impl FieldElement56 {
    pub fn is_negative(&self) -> Choice {
        let bytes = self.to_bytes();
        (bytes[0] & 1).into()
    }
}

///
/// Serialisation
///
impl FieldElement56 {
    /// Helper function for internally constructing a field element
    pub(crate) const fn from_raw_slice(slice: [u64; 8]) -> FieldElement56 {
        FieldElement56(fiat_p448_tight_field_element(slice))
    }

    /// This does not check if the encoding is canonical (ie if the input is reduced)
    /// We parse in chunks of 56 bytes, the first 28 bytes will contain the i'th limb
    /// and the second 28 bytes will contain the (2i+1)'th limb
    pub(crate) fn from_bytes(bytes: &[u8; 56]) -> FieldElement56 {
        let mut res = FieldElement56::zero();
        fiat_p448_from_bytes(&mut res.0, bytes);
        res
    }

    // We encode the Field element by storing each consecutive into a u64
    pub(crate) fn to_bytes(&self) -> [u8; 56] {
        let mut res = [0u8; 56];
        fiat_p448_to_bytes(&mut res, &self.0);
        res
    }
}

impl FieldElement56 {
    /// Squares a field element
    pub(crate) fn square(&self) -> FieldElement56 {
        let mut self_loose = fiat_p448_loose_field_element([0; 8]);
        fiat_p448_relax(&mut self_loose, &self.0);
        let mut result = FieldElement56::zero();
        fiat_p448_carry_square(&mut result.0, &self_loose);
        result
    }
    /// Negates a field element
    pub(crate) fn negate(&self) -> FieldElement56 {
        let mut result_loose = fiat_p448_loose_field_element([0; 8]);
        fiat_p448_opp(&mut result_loose, &self.0);
        let mut result = FieldElement56::zero();
        fiat_p448_carry(&mut result.0, &result_loose);
        result
    }

    /// Reduces the field element to a canonical representation
    /// This is used when checking equality between two field elements and
    /// when encoding a field element
    pub(crate) fn strong_reduce(&mut self) {
        let mut self_loose = fiat_p448_loose_field_element([0; 8]);
        fiat_p448_relax(&mut self_loose, &self.0);
        fiat_p448_carry(&mut self.0, &self_loose);
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_negate() {
        let x = FieldElement56::zero();
        let y = x.negate();
        assert_eq!(y.to_bytes(), [0u8; 56]);
    }
}
