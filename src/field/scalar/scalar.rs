// This is the scalar field
// size = 4q = 2^446 - 0x8335dc163bb124b65129c96fde933d8d723a70aadc873d6d54a7bb0d
// We can therefore use 14 saturated 32-bit limbs
// in LE format
/// note: Montgomery reduction and barret are both used: montgomery when multiplying and barret when decoding
/// Probably easier to stick to montgomery I think
pub struct Scalar([u32; 14]);

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
}
