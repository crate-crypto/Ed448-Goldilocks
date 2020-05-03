use crate::constants::ONE_MINUS_D;
use crate::curve::edwards::extended::ExtendedPoint;
use crate::field::{FieldElement, Scalar};
use std::fmt;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

#[derive(Copy, Clone)]
pub struct MontgomeryPoint(pub [u8; 56]);

impl fmt::Debug for MontgomeryPoint {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        self.0[..].fmt(formatter)
    }
}

impl ConstantTimeEq for MontgomeryPoint {
    fn ct_eq(&self, other: &MontgomeryPoint) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl PartialEq for MontgomeryPoint {
    fn eq(&self, other: &MontgomeryPoint) -> bool {
        self.ct_eq(other).into()
    }
}
impl Eq for MontgomeryPoint {}

#[derive(Copy, Clone)]
pub struct ProjectiveMontgomeryPoint {
    U: FieldElement,
    W: FieldElement,
}

impl MontgomeryPoint {
    pub fn to_edwards(&self, sign: u8) -> Option<ExtendedPoint> {
        // We use the 4-isogeny to map to the Ed448.
        // This is different to Curve25519, where we use a birational map.
        todo!()
    }

    pub fn to_projective(&self) -> ProjectiveMontgomeryPoint {
        ProjectiveMontgomeryPoint {
            U: FieldElement::from_bytes(&self.0),
            W: FieldElement::one(),
        }
    }

    // Taken from Dalek
    pub fn mul(&self, scalar: &Scalar) -> MontgomeryPoint {
        // Algorithm 8 of Costello-Smith 2017
        let affine_u = FieldElement::from_bytes(&self.0);
        let mut x0 = ProjectiveMontgomeryPoint::identity();
        let mut x1 = ProjectiveMontgomeryPoint {
            U: affine_u,
            W: FieldElement::one(),
        };

        let bits = scalar.bits();
        let mut swap = 0;
        for s in (0..448).rev() {
            let bit = bits[s] as u8;
            let choice: u8 = (swap ^ bit) as u8;

            ProjectiveMontgomeryPoint::conditional_swap(&mut x0, &mut x1, Choice::from(choice));
            differential_add_and_double(&mut x0, &mut x1, &affine_u);

            swap = bit;
        }

        x0.to_affine()
    }
}

impl ConditionallySelectable for ProjectiveMontgomeryPoint {
    fn conditional_select(
        a: &ProjectiveMontgomeryPoint,
        b: &ProjectiveMontgomeryPoint,
        choice: Choice,
    ) -> ProjectiveMontgomeryPoint {
        ProjectiveMontgomeryPoint {
            U: FieldElement::conditional_select(&a.U, &b.U, choice),
            W: FieldElement::conditional_select(&a.W, &b.W, choice),
        }
    }
}

// Also taken from Dalek
fn differential_add_and_double(
    P: &mut ProjectiveMontgomeryPoint,
    Q: &mut ProjectiveMontgomeryPoint,
    affine_PmQ: &FieldElement,
) {
    let a24 = ONE_MINUS_D; //39082

    let t0 = P.U + P.W;
    let t1 = P.U - P.W;
    let t2 = Q.U + Q.W;
    let t3 = Q.U - Q.W;

    let t4 = t0.square(); // (U_P + W_P)^2 = U_P^2 + 2 U_P W_P + W_P^2
    let t5 = t1.square(); // (U_P - W_P)^2 = U_P^2 - 2 U_P W_P + W_P^2

    let t6 = t4 - t5; // 4 U_P W_P

    let t7 = t0 * t3; // (U_P + W_P) (U_Q - W_Q) = U_P U_Q + W_P U_Q - U_P W_Q - W_P W_Q
    let t8 = t1 * t2; // (U_P - W_P) (U_Q + W_Q) = U_P U_Q - W_P U_Q + U_P W_Q - W_P W_Q

    let t9 = t7 + t8; // 2 (U_P U_Q - W_P W_Q)
    let t10 = t7 - t8; // 2 (W_P U_Q - U_P W_Q)

    let t11 = t9.square(); // 4 (U_P U_Q - W_P W_Q)^2
    let t12 = t10.square(); // 4 (W_P U_Q - U_P W_Q)^2
    let t13 = a24 * t6; // (A + 2) U_P U_Q

    let t14 = t4 * t5; // ((U_P + W_P)(U_P - W_P))^2 = (U_P^2 - W_P^2)^2
    let t15 = t13 + t5; // (U_P - W_P)^2 + (A + 2) U_P W_P

    let t16 = t6 * t15; // 4 (U_P W_P) ((U_P - W_P)^2 + (A + 2) U_P W_P)
    let t17 = *affine_PmQ * t12; // U_D * 4 (W_P U_Q - U_P W_Q)^2
    let t18 = t11; // W_D * 4 (U_P U_Q - W_P W_Q)^2

    P.U = t14; // U_{P'} = (U_P + W_P)^2 (U_P - W_P)^2
    P.W = t16; // W_{P'} = (4 U_P W_P) ((U_P - W_P)^2 + ((A + 2)/4) 4 U_P W_P)
    Q.U = t18; // U_{Q'} = W_D * 4 (U_P U_Q - W_P W_Q)^2
    Q.W = t17; // W_{Q'} = U_D * 4 (W_P U_Q - U_P W_Q)^2
}

impl ProjectiveMontgomeryPoint {
    pub fn identity() -> ProjectiveMontgomeryPoint {
        ProjectiveMontgomeryPoint {
            U: FieldElement::one(),
            W: FieldElement::zero(),
        }
    }

    pub fn to_affine(&self) -> MontgomeryPoint {
        let x = self.U * self.W.invert();
        MontgomeryPoint(x.to_bytes())
    }
    pub fn double(&self) -> ProjectiveMontgomeryPoint {
        todo!()
    }
    pub fn add(&self, other: &ProjectiveMontgomeryPoint) -> ProjectiveMontgomeryPoint {
        todo!()
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use hex::decode as hex_decode;

    fn slice_to_fixed_array(b: &[u8]) -> [u8; 56] {
        let mut a: [u8; 56] = [0; 56];
        a.copy_from_slice(&b);
        a
    }

    fn hex_to_array(data: &str) -> [u8; 56] {
        let mut bytes = hex_decode(data).unwrap();
        slice_to_fixed_array(&bytes)
    }

    fn clamp(scalar: &mut [u8; 56]) {
        scalar[0] &= 252;
        scalar[55] |= 128;
    }

    // This will not stay here, it's only here so that we can be sure that the Ladder is being computed correctly
    #[test]
    fn test_rfc_vector() {
        // Load RFC basepoint
        let point_bytes = hex_to_array("0500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
        let point = MontgomeryPoint(point_bytes);
        // Clamp Scalar
        let mut scalar_bytes = hex_to_array("9a8f4925d1519f5775cf46b04b5800d4ee9ee8bae8bc5565d498c28dd9c9baf574a9419744897391006382a6f127ab1d9ac2d8c0a598726b");
        clamp(&mut scalar_bytes);
        // Load RFC scalar
        let scalar = Scalar::from_bytes(scalar_bytes);
        // Compute output
        let output = point.mul(&scalar);
        // Expected output
        let expected = hex_to_array("9b08f7cc31b7e3e67d22d5aea121074a273bd2b83de09c63faa73d2c22c5d9bbc836647241d953d40c5b12da88120d53177f80e532c41fa0");
        assert_eq!(&output.0[..], &expected[..]);
    }
}
