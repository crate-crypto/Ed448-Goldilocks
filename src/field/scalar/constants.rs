// Fow now, this will be specific to 32-bit saturated limbs

/// The number of limbs used to represent a scalar
pub(crate) const NUM_LIMBS: usize = 14;

/// The number of bits used to represent a single limb
/// This is also known as the number of bits in a word
pub(crate) const WORD_BITS: usize = 32;
