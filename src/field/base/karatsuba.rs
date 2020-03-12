use super::Fq;

// TODO: Fully comment code : (showing why it works-golden ratio equation)
// Code originally from Mike Hamburg.Referenced from otrv4/ed448
/// Karatsuba algorithm:
/// Given two field elements `a` and `b` in base `r` where `a` and `b` are length `n`
/// Divide each field element into halves `a_hi`, `a_lo`,`b_hi`,`b_lo` where `a_hi`denotes high bits of `a`
/// a = a_hi * r^n/2 + a_lo
/// b = b_hi * r^n/2 + b_lo
/// a*b = (a_hi * r^n/2 + a_lo)(b_hi * r^n/2 + b_lo) = (a_hi * b_hi * r^n) + (a_hi * b_lo + a_lo * b_hi) * r^n/2 + a_lo * b_lo
/// XXX: Write the full algorithm out for the case of golden-ratio primes

pub(crate) fn mul(a: &Fq, b: &Fq) -> Fq {
    let (mut accum0, mut accum1, mut accum2) = (0u64, 0u64, 0u64);
    const MASK: u32 = (1 << 28) - 1;

    // result
    let mut c = Fq::zero();

    let widemul = |a: u32, b: u32| -> u64 { (a as u64) * (b as u64) };

    let mut aa = [0u64; 8];
    let mut bb = [0u64; 8];
    for i in 0..8 {
        aa[i] = (a[i] as u64) + (a[i + 8] as u64);
        bb[i] = (b[i] as u64) + (b[i + 8] as u64);
    }

    for j in 0..8 {
        accum2 = 0;
        for i in 0..j + 1 {
            accum2 = accum2.wrapping_add(widemul(a[j - i], b[i]));
            accum1 = accum1.wrapping_add(aa[j - i] * bb[i]);
            accum0 = accum0.wrapping_add(widemul(a[8 + j - i], b[8 + i]));
        }
        accum1 = accum1.wrapping_sub(accum2);
        accum0 = accum0.wrapping_add(accum2);
        accum2 = 0;

        for i in j + 1..8 {
            accum0 = accum0.wrapping_sub(widemul(a[8 + j - i], b[i]));
            accum2 = accum2.wrapping_add(aa[8 + j - i] * bb[i]);
            accum1 = accum1.wrapping_add(widemul(a[16 + j - i], b[8 + i]));
        }
        accum1 = accum1.wrapping_add(accum2);
        accum0 = accum0.wrapping_add(accum2);

        c[j] = (accum0 as u32) & MASK;
        c[j + 8] = (accum1 as u32) & MASK;

        accum0 >>= 28;
        accum1 >>= 28;
    }

    accum0 += accum1;
    accum0 += c[8] as u64;
    accum1 += c[0] as u64;
    c[8] = (accum0 as u32) & MASK;
    c[0] = (accum1 as u32) & MASK;
    accum0 >>= 28;
    accum1 >>= 28;
    c[9] += accum0 as u32;
    c[1] += accum1 as u32;

    c
}

#[test]
fn test_mul() {
    let a = Fq::from(5u8);
    let b = Fq::from(210u32);
    let c = Fq::from(1050u32);

    let expected_c = mul(&a, &b);
    assert_eq!(expected_c, c);
}
#[test]
fn test_naive_square() {
    let a = Fq::from(5u8);
    let b = Fq::from(5u32);
    let c = Fq::from(25u32);

    let expected_c = mul(&a, &b);
    assert_eq!(expected_c, c);
}
