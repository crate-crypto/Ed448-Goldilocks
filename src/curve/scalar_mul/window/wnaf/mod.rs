use crate::curve::twedwards::extended::ExtendedPoint;
use crate::curve::twedwards::projective::ProjectiveNielsPoint;
use subtle::{ConditionallySelectable, ConstantTimeEq};

pub(crate) const WINDOW: usize = 5;
pub(crate) const WINDOW_MASK: usize = (1 << WINDOW) - 1;
pub(crate) const WINDOW_T_MASK: usize = WINDOW_MASK >> 1;
pub(crate) const TABLE_SIZE: usize = 16;

pub struct LookupTable([ProjectiveNielsPoint; TABLE_SIZE]);

/// Precomputes odd multiples of the point passed in
impl From<&ExtendedPoint> for LookupTable {
    fn from(point: &ExtendedPoint) -> LookupTable {
        let P = point.to_extensible();
        let P2 = P.double().to_projective_niels();

        let mut table = [P.to_projective_niels(); TABLE_SIZE];
        let mut p_original = point.to_extensible();

        for i in 1..TABLE_SIZE {
            p_original = p_original.add_projective_niels(&P2);
            table[i] = p_original.to_projective_niels();
        }

        LookupTable(table)
    }
}

impl LookupTable {
    /// Selects a projective niels point from a lookup table in constant time
    pub fn select(&self, index: u32) -> ProjectiveNielsPoint {
        let mut result = ProjectiveNielsPoint::identity();

        for i in 0..TABLE_SIZE {
            let swap = index.ct_eq(&(i as u32));
            result.conditional_assign(&self.0[i], swap);
        }
        result
    }
}
// XXX: Add back tests to ensure that select works correctly
