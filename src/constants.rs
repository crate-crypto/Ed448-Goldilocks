use crate::decaf::DecafPoint;

#[cfg(feature = "u32_backend")]
pub(crate) use crate::field::u32::constants::*;

pub(crate) const DECAF_BASEPOINT: DecafPoint = DecafPoint(TWISTED_EDWARDS_BASE_POINT);

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_base_point() {
        let goldilocks_base = GOLDILOCKS_BASE_POINT;
        // Use the isogeny to get the twisted edwards basepoint
        let twisted_edwards_base = goldilocks_base.to_twisted();
        assert!(twisted_edwards_base == TWISTED_EDWARDS_BASE_POINT)
    }
}
