#[cfg(feature = "fiat_u64_backend")]
pub mod fiat_u64;

#[cfg(feature = "u32_backend")]
pub mod u32;

// XXX: Currently we only have one implementation for Scalar
mod scalar;
pub use crate::field::scalar::Scalar;

#[cfg(feature = "u32_backend")]
pub type FieldElement = crate::field::u32::FieldElement28;

#[cfg(feature = "fiat_u64_backend")]
pub type FieldElement = crate::field::fiat_u64::FieldElement56;

impl FieldElement {
    /// Inverts a field element
    pub(crate) fn invert(&self) -> FieldElement {
        let mut t1 = self.square();
        let (mut t2, _) = t1.inverse_square_root();
        t1 = t2.square();
        t2 = t1 * self;
        t2
    }
    /// Squares a field element  `n` times
    fn square_n(&self, mut n: u32) -> FieldElement {
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
    pub(crate) fn inverse_square_root(&self) -> (FieldElement, bool) {
        let (mut l0, mut l1, mut l2) = (
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
        );

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

        let is_residue = l0 == FieldElement::one();
        (l1, is_residue)
    }

    /// Computes the square root ratio of two elements
    pub(crate) fn sqrt_ratio(u: &FieldElement, v: &FieldElement) -> (FieldElement, bool) {
        let x = *u * v;
        let (inv_sqrt_x, is_res) = x.inverse_square_root();
        (inv_sqrt_x * u, is_res)
    }
}
