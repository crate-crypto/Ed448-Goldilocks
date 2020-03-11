use std::ops::{Add, AddAssign, BitAnd, BitXor, Mul, Shr, Sub};

/// Field elements are split into chunks of smaller integers
/// These integers are called limbs
/// Currently, this library uses (unsaturated) 28-bit limbs to represent field elements
/// to reduce the number of carries.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct Limb28(pub(crate) u32);

impl Limb28 {
    pub(crate) fn zero() -> Limb28 {
        Limb28(0)
    }
    pub(crate) fn value(&self) -> u32 {
        self.0
    }
    /// converts a u64 to a u32
    /// Fails if the u64 cannot safely fit into a u32
    pub(crate) fn from_checked_u64(a: u64) -> Limb28 {
        // Check if this number can be safely casted to a u64
        if a > std::u32::MAX as u64 {
            panic!("cannot cast this to an unsaturated 28-bit limb safely")
        }
        Limb28::from(a as u32)
    }
    /// converts a i64 to a u32
    /// Fails if the u64 cannot safely fit into a u32
    pub(crate) fn from_checked_i64(a: i64) -> Limb28 {
        // Check if this number can be safely casted to a u64
        if a > std::u32::MAX as i64 {
            panic!("cannot cast this to an unsaturated 28-bit limb safely")
        }
        Limb28::from(a as u32)
    }
}
// Trait implementations

/// From converts a u32 into a Limb28
/// There is no guarantee that the value is less than 2^28
impl From<u32> for Limb28 {
    fn from(a: u32) -> Limb28 {
        Limb28(a)
    }
}
/// Multiplying two limbs together produces a u64
/// So that we can manipulate the carry bit
impl Mul<Limb28> for Limb28 {
    type Output = u64;
    fn mul(self, rhs: Limb28) -> Self::Output {
        (self.0 as u64) * (rhs.0 as u64)
    }
}
/// Multiplying a u32 by a Limb is the same as
/// multiplying a Limb by a Limb
impl Mul<u32> for Limb28 {
    type Output = u64;
    fn mul(self, rhs: u32) -> Self::Output {
        self * Limb28::from(rhs)
    }
}

/// Adding two limbs produces another limb
/// XXX: In this case, we lose the carry bit
impl Add<Limb28> for Limb28 {
    type Output = Limb28;
    fn add(self, rhs: Limb28) -> Self::Output {
        Limb28::from(self.0 + rhs.0)
    }
}
/// Subtracting two limbs produces another limb]
/// XXX: In this case, we lose the carry bit
impl Sub<Limb28> for Limb28 {
    type Output = Limb28;
    fn sub(self, rhs: Limb28) -> Self::Output {
        Limb28::from(self.0 - rhs.0)
    }
}
/// Adding a u64 onto a limb produces a u64
impl Add<u64> for Limb28 {
    type Output = u64;
    fn add(self, rhs: u64) -> Self::Output {
        let a = self.0 as u64;
        a + rhs
    }
}
impl AddAssign<Limb28> for Limb28 {
    fn add_assign(&mut self, rhs: Limb28) {
        *self = Limb28(self.0 + rhs.0)
    }
}

impl Shr<usize> for Limb28 {
    type Output = Limb28;
    fn shr(self, shift_by: usize) -> Limb28 {
        let x = self.0 >> shift_by;
        Limb28::from(x)
    }
}
impl BitAnd<u32> for Limb28 {
    type Output = Limb28;
    fn bitand(self, operand: u32) -> Self::Output {
        let x = self.0 & operand;
        Limb28::from(x)
    }
}
impl BitXor<Limb28> for Limb28 {
    type Output = u32;
    fn bitxor(self, rhs: Limb28) -> Self::Output {
        self.0 ^ rhs.0
    }
}
