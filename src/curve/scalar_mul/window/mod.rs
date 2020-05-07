pub mod wnaf;

// If they are the same return -1, if not return 0
fn select_mask(index: u32, current: u32) -> u32 {
    let equals = (index == current) as u32;
    0u32.wrapping_sub(equals)
}
