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

const MASK: u32 = (1 << 28) - 1;

pub(crate) fn mul(a: &Fq, b: &Fq) -> Fq {
    let (mut accum0, mut accum1, mut accum2) = (0u64, 0u64, 0u64);

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
            accum1 = accum1.wrapping_add(aa[j - i].wrapping_mul(bb[i]));
            accum0 = accum0.wrapping_add(widemul(a[8 + j - i], b[8 + i]));
        }
        accum1 = accum1.wrapping_sub(accum2);
        accum0 = accum0.wrapping_add(accum2);
        accum2 = 0;

        for i in j + 1..8 {
            accum0 = accum0.wrapping_sub(widemul(a[8 + j - i], b[i]));
            accum2 = accum2.wrapping_add(aa[8 + j - i].wrapping_mul(bb[i]));
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

pub(crate) fn square(a: &Fq) -> Fq {
    let (mut accum0, mut accum1, mut accum2) = (0u64, 0u64, 0u64);

    // result
    let mut c = Fq::zero();

    let widemul = |a: u32, b: u32| -> u64 { (a as u64) * (b as u64) };

    let mut aa = [0u64; 8];
    for i in 0..8 {
        aa[i] = (a[i] as u64) + (a[i + 8] as u64);
    }

    // start j = 0
    accum2 = 0;
    accum2 = accum2.wrapping_add(widemul(a[0], a[0]));
    accum1 = accum1.wrapping_add(aa[0].wrapping_mul(aa[0]));
    accum0 = accum0.wrapping_add(widemul(a[8], a[8]));
    accum1 = accum1.wrapping_sub(accum2);
    accum0 = accum0.wrapping_add(accum2);
    accum2 = 0;
    for i in 1..4 {
        accum2 = accum2.wrapping_add(aa[i].wrapping_mul(aa[8 - i]) << 1);
        accum1 = accum1.wrapping_add(widemul(a[8 + i], a[16 - i]) << 1);
        accum0 = accum0.wrapping_sub(widemul(a[i], a[8 - i]) << 1);
    }
    accum2 = accum2.wrapping_add(aa[4].wrapping_mul(aa[4]));
    accum1 = accum1
        .wrapping_add(widemul(a[12], a[12]))
        .wrapping_add(accum2);
    accum0 = accum0
        .wrapping_sub(widemul(a[4], a[4]))
        .wrapping_add(accum2);
    c[0] = (accum0 as u32) & MASK;
    c[8] = (accum1 as u32) & MASK;
    accum0 >>= 28;
    accum1 >>= 28;
    // end j = 0

    // start j = 1..7
    for j in 1..7 {
        // first round
        accum2 = 0;
        for i in 0..((j / 2) + 1) {
            if i != j - i {
                accum2 = accum2.wrapping_add(widemul(a[j - i], a[i]) << 1);
                accum1 = accum1.wrapping_add(aa[j - i].wrapping_mul(aa[i]) << 1);
                accum0 = accum0.wrapping_add(widemul(a[8 + (j - i)], a[8 + i]) << 1);
            } else {
                accum2 = accum2.wrapping_add(widemul(a[j / 2], a[j / 2]));
                accum1 = accum1.wrapping_add(aa[j / 2].wrapping_mul(aa[j / 2]));
                accum0 = accum0.wrapping_add(widemul(a[8 + (j / 2)], a[8 + (j / 2)]));
            }
        }
        accum1 = accum1.wrapping_sub(accum2);
        accum0 = accum0.wrapping_add(accum2);

        // second round
        accum2 = 0;
        for i in 1..(((8 - j) / 2) + 1) {
            if j + i != 8 - i {
                accum2 = accum2.wrapping_add(aa[8 - i].wrapping_mul(aa[j + i]) << 1);
                accum1 = accum1.wrapping_add(widemul(a[16 - i], a[8 + j + i]) << 1);
                accum0 = accum0.wrapping_sub(widemul(a[8 - i], a[j + i]) << 1);
            } else {
                accum2 = accum2.wrapping_add(aa[j + i].wrapping_mul(aa[j + i]));
                accum1 = accum1.wrapping_add(widemul(a[8 + j + i], a[8 + j + i]));
                accum0 = accum0.wrapping_sub(widemul(a[j + i], a[j + i]));
            }
        }
        accum1 = accum1.wrapping_add(accum2);
        accum0 = accum0.wrapping_add(accum2);

        c[j] = (accum0 as u32) & MASK;
        c[j + 8] = (accum1 as u32) & MASK;
        accum0 >>= 28;
        accum1 >>= 28;
    }
    // end j = 1..7

    // start j = 7
    accum2 = 0;
    for i in 0..4 {
        accum2 = accum2.wrapping_add(widemul(a[7 - i], a[i]) << 1);
        accum1 = accum1.wrapping_add(aa[7 - i].wrapping_mul(aa[i]) << 1);
        accum0 = accum0.wrapping_add(widemul(a[15 - i], a[8 + i]) << 1);
    }
    accum1 = accum1.wrapping_sub(accum2);
    accum0 = accum0.wrapping_add(accum2);
    c[7] = (accum0 as u32) & MASK;
    c[15] = (accum1 as u32) & MASK;
    accum0 >>= 28;
    accum1 >>= 28;
    // end j = 7

    // final
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
fn test_basic_mul() {
    let a = Fq::from(5u8);
    let b = Fq::from(210u32);
    let c = Fq::from(1050u32);

    let expected_c = mul(&a, &b);
    assert_eq!(expected_c, c);
}

#[test]
fn test_mul() {
    let x = Fq::from_bytes(&[
        0xf5, 0x81, 0x74, 0xd5, 0x7a, 0x33, 0x72, 0x36, 0x3c, 0x0d, 0x9f, 0xcf, 0xaa, 0x3d, 0xc1,
        0x8b, 0x1e, 0xff, 0x7e, 0x89, 0xbf, 0x76, 0x78, 0x63, 0x65, 0x80, 0xd1, 0x7d, 0xd8, 0x4a,
        0x87, 0x3b, 0x14, 0xb9, 0xc0, 0xe1, 0x68, 0x0b, 0xbd, 0xc8, 0x76, 0x47, 0xf3, 0xc3, 0x82,
        0x90, 0x2d, 0x2f, 0x58, 0xd2, 0x75, 0x4b, 0x39, 0xbc, 0xa8, 0x74,
    ]);

    let y = Fq::from_bytes(&[
        0x74, 0xa8, 0xbc, 0x39, 0x4b, 0x75, 0xd2, 0x58, 0x2f, 0x2d, 0x90, 0x82, 0xc3, 0xf3, 0x47,
        0x76, 0xc8, 0xbd, 0x0b, 0x68, 0xe1, 0xc0, 0xb9, 0x14, 0x3b, 0x87, 0x4a, 0xd8, 0x7d, 0xd1,
        0x80, 0x65, 0x63, 0x78, 0x76, 0xbf, 0x89, 0x7e, 0xff, 0x1e, 0x8b, 0xc1, 0x3d, 0xaa, 0xcf,
        0x9f, 0x0d, 0x3c, 0x36, 0x72, 0x33, 0x7a, 0xd5, 0x74, 0x81, 0xf5,
    ]);

    let expected = Fq::from_bytes(&[
        0x11, 0x95, 0x9c, 0x2e, 0x91, 0x78, 0x6f, 0xec, 0xff, 0x37, 0xe5, 0x8e, 0x2b, 0x50, 0x9e,
        0xf8, 0xfb, 0x41, 0x08, 0xc4, 0xa7, 0x02, 0x1c, 0xbf, 0x5a, 0x9f, 0x18, 0xa7, 0xec, 0x32,
        0x65, 0x7e, 0xed, 0xdc, 0x81, 0x81, 0x80, 0xa8, 0x4c, 0xdd, 0x95, 0x14, 0xe6, 0x67, 0x26,
        0xd3, 0xa1, 0x22, 0xdb, 0xb3, 0x9f, 0x17, 0x7a, 0x85, 0x16, 0x6c,
    ]);

    let result = mul(&x, &y);
    assert_eq!(result, expected)
}
#[test]
fn test_naive_square() {
    let a = Fq::from(5u8);
    let b = Fq::from(5u32);
    let c = Fq::from(25u32);

    let expected_c = mul(&a, &b);
    assert_eq!(expected_c, c);
}
#[test]
fn test_karatsuba_square() {
    let a = Fq::from(5u8);
    let expected_c = Fq::from(25u32);
    let c = square(&a);
    assert_eq!(expected_c, c);

    let d = Fq([
        268396374, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
        268435454, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455, 268435455,
    ]);
    let expected_d2 = mul(&d, &d);
    let d2 = square(&d);
    assert_eq!(expected_d2, d2);
}
