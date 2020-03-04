use super::{Fq, Limb28};

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
    const MASK: u64 = (1 << 28) - 1;

    let mut c = Fq::zero();

    let mut aa = [Limb28::zero(); 8];
    let mut bb = [Limb28::zero(); 8];
    for i in 0..8 {
        aa[i] = a[i] + a[i + 8];
        bb[i] = b[i] + b[i + 8];
    }

    let (mut z0, mut z1, mut z2) = (0u64, 0u64, 0u64);

    // j=0
    z2 = a[0] * b[0];
    z1 += aa[0] * bb[0];
    z1 -= z2;
    z0 += a[8] * b[8];
    z0 += z2;

    z2 = 0;
    z2 += aa[7] * bb[1];
    z2 += aa[6] * bb[2];
    z2 += aa[5] * bb[3];
    z2 += aa[4] * bb[4];
    z2 += aa[3] * bb[5];
    z2 += aa[2] * bb[6];
    z2 += aa[1] * bb[7];

    z1 += a[15] * b[9];
    z1 += a[14] * b[10];
    z1 += a[13] * b[11];
    z1 += a[12] * b[12];
    z1 += a[11] * b[13];
    z1 += a[10] * b[14];
    z1 += a[9] * b[15];
    z1 += z2;

    z0 -= a[7] * b[1];
    z0 -= a[6] * b[2];
    z0 -= a[5] * b[3];
    z0 -= a[4] * b[4];
    z0 -= a[3] * b[5];
    z0 -= a[2] * b[6];
    z0 -= a[1] * b[7];
    z0 += z2;

    c[0] = Limb28::from_checked_u64(z0 & MASK);
    c[8] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    //j=1
    z2 = 0;
    z2 += a[1] * b[0];
    z2 += a[0] * b[1];

    z1 += aa[1] * bb[0];
    z1 += aa[0] * bb[1];
    z1 -= z2;

    z0 += a[9] * b[8];
    z0 += a[8] * b[9];
    z0 += z2;

    z2 = 0;
    z2 += aa[7] * bb[2];
    z2 += aa[6] * bb[3];
    z2 += aa[5] * bb[4];
    z2 += aa[4] * bb[5];
    z2 += aa[3] * bb[6];
    z2 += aa[2] * bb[7];

    z1 += a[15] * b[10];
    z1 += a[14] * b[11];
    z1 += a[13] * b[12];
    z1 += a[12] * b[13];
    z1 += a[11] * b[14];
    z1 += a[10] * b[15];
    z1 += z2;

    z0 -= a[7] * b[2];
    z0 -= a[6] * b[3];
    z0 -= a[5] * b[4];
    z0 -= a[4] * b[5];
    z0 -= a[3] * b[6];
    z0 -= a[2] * b[7];
    z0 += z2;

    c[1] = Limb28::from_checked_u64(z0 & MASK);
    c[9] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    //j = 2
    z2 = 0;
    z2 += a[2] * b[0];
    z2 += a[1] * b[1];
    z2 += a[0] * b[2];

    z1 += aa[2] * bb[0];
    z1 += aa[1] * bb[1];
    z1 += aa[0] * bb[2];
    z1 -= z2;

    z0 += a[10] * b[8];
    z0 += a[9] * b[9];
    z0 += a[8] * b[10];
    z0 += z2;

    z2 = 0;
    z2 += aa[7] * bb[3];
    z2 += aa[6] * bb[4];
    z2 += aa[5] * bb[5];
    z2 += aa[4] * bb[6];
    z2 += aa[3] * bb[7];

    z1 += a[15] * b[11];
    z1 += a[14] * b[12];
    z1 += a[13] * b[13];
    z1 += a[12] * b[14];
    z1 += a[11] * b[15];
    z1 += z2;

    z0 -= a[7] * b[3];
    z0 -= a[6] * b[4];
    z0 -= a[5] * b[5];
    z0 -= a[4] * b[6];
    z0 -= a[3] * b[7];
    z0 += z2;

    c[2] = Limb28::from_checked_u64(z0 & MASK);
    c[10] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    //j = 3
    z2 = 0;
    z2 += a[3] * b[0];
    z2 += a[2] * b[1];
    z2 += a[1] * b[2];
    z2 += a[0] * b[3];

    z1 += aa[3] * bb[0];
    z1 += aa[2] * bb[1];
    z1 += aa[1] * bb[2];
    z1 += aa[0] * bb[3];
    z1 -= z2;

    z0 += a[11] * b[8];
    z0 += a[10] * b[9];
    z0 += a[9] * b[10];
    z0 += a[8] * b[11];
    z0 += z2;

    z2 = 0;
    z2 += aa[7] * bb[4];
    z2 += aa[6] * bb[5];
    z2 += aa[5] * bb[6];
    z2 += aa[4] * bb[7];

    z0 -= a[7] * b[4];
    z0 -= a[6] * b[5];
    z0 -= a[5] * b[6];
    z0 -= a[4] * b[7];
    z0 += z2;

    z1 += a[15] * b[12];
    z1 += a[14] * b[13];
    z1 += a[13] * b[14];
    z1 += a[12] * b[15];
    z1 += z2;

    c[3] = Limb28::from_checked_u64(z0 & MASK);
    c[11] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    //j = 4
    z2 = 0;
    z2 += a[4] * b[0];
    z2 += a[3] * b[1];
    z2 += a[2] * b[2];
    z2 += a[1] * b[3];
    z2 += a[0] * b[4];

    z1 += aa[4] * bb[0];
    z1 += aa[3] * bb[1];
    z1 += aa[2] * bb[2];
    z1 += aa[1] * bb[3];
    z1 += aa[0] * bb[4];
    z1 -= z2;

    z0 += a[12] * b[8];
    z0 += a[11] * b[9];
    z0 += a[10] * b[10];
    z0 += a[9] * b[11];
    z0 += a[8] * b[12];
    z0 += z2;

    z2 = 0;
    z2 += aa[7] * bb[5];
    z2 += aa[6] * bb[6];
    z2 += aa[5] * bb[7];

    z1 += a[15] * b[13];
    z1 += a[14] * b[14];
    z1 += a[13] * b[15];
    z1 += z2;

    z0 -= a[7] * b[5];
    z0 -= a[6] * b[6];
    z0 -= a[5] * b[7];
    z0 += z2;

    c[4] = Limb28::from_checked_u64(z0 & MASK);
    c[12] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    //j = 5
    z2 = 0;
    z2 += a[5] * b[0];
    z2 += a[4] * b[1];
    z2 += a[3] * b[2];
    z2 += a[2] * b[3];
    z2 += a[1] * b[4];
    z2 += a[0] * b[5];

    z1 += aa[5] * bb[0];
    z1 += aa[4] * bb[1];
    z1 += aa[3] * bb[2];
    z1 += aa[2] * bb[3];
    z1 += aa[1] * bb[4];
    z1 += aa[0] * bb[5];
    z1 -= z2;

    z0 += a[13] * b[8];
    z0 += a[12] * b[9];
    z0 += a[11] * b[10];
    z0 += a[10] * b[11];
    z0 += a[9] * b[12];
    z0 += a[8] * b[13];
    z0 += z2;

    z2 = 0;
    z2 += aa[7] * bb[6];
    z2 += aa[6] * bb[7];

    z1 += a[15] * b[14];
    z1 += a[14] * b[15];
    z1 += z2;

    z0 -= a[7] * b[6];
    z0 -= a[6] * b[7];
    z0 += z2;

    c[5] = Limb28::from_checked_u64(z0 & MASK);
    c[13] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    //j = 6
    z2 = 0;
    z2 += a[6] * b[0];
    z2 += a[5] * b[1];
    z2 += a[4] * b[2];
    z2 += a[3] * b[3];
    z2 += a[2] * b[4];
    z2 += a[1] * b[5];
    z2 += a[0] * b[6];

    z1 += aa[6] * bb[0];
    z1 += aa[5] * bb[1];
    z1 += aa[4] * bb[2];
    z1 += aa[3] * bb[3];
    z1 += aa[2] * bb[4];
    z1 += aa[1] * bb[5];
    z1 += aa[0] * bb[6];
    z1 -= z2;

    z0 += a[14] * b[8];
    z0 += a[13] * b[9];
    z0 += a[12] * b[10];
    z0 += a[11] * b[11];
    z0 += a[10] * b[12];
    z0 += a[9] * b[13];
    z0 += a[8] * b[14];
    z0 += z2;

    z2 = 0;
    z2 += aa[7] * bb[7];
    z1 += a[15] * b[15];
    z1 += z2;
    z0 -= a[7] * b[7];
    z0 += z2;

    c[6] = Limb28::from_checked_u64(z0 & MASK);
    c[14] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    //j = 7
    z2 = 0;
    z2 += a[7] * b[0];
    z2 += a[6] * b[1];
    z2 += a[5] * b[2];
    z2 += a[4] * b[3];
    z2 += a[3] * b[4];
    z2 += a[2] * b[5];
    z2 += a[1] * b[6];
    z2 += a[0] * b[7];

    z1 += aa[7] * bb[0];
    z1 += aa[6] * bb[1];
    z1 += aa[5] * bb[2];
    z1 += aa[4] * bb[3];
    z1 += aa[3] * bb[4];
    z1 += aa[2] * bb[5];
    z1 += aa[1] * bb[6];
    z1 += aa[0] * bb[7];
    z1 -= z2;

    z0 += a[15] * b[8];
    z0 += a[14] * b[9];
    z0 += a[13] * b[10];
    z0 += a[12] * b[11];
    z0 += a[11] * b[12];
    z0 += a[10] * b[13];
    z0 += a[9] * b[14];
    z0 += a[8] * b[15];
    z0 += z2;

    z2 = 0;
    z1 += z2;
    z0 += z2;

    c[7] = Limb28::from_checked_u64(z0 & MASK);
    c[15] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    // Final step

    z0 += z1;
    z0 += c[8].value() as u64;
    z1 += c[0].value() as u64;

    c[8] = Limb28::from_checked_u64(z0 & MASK);
    c[0] = Limb28::from_checked_u64(z1 & MASK);

    z0 >>= 28;
    z1 >>= 28;

    c[9] += Limb28::from_checked_u64(z0 & MASK);
    c[1] += Limb28::from_checked_u64(z1 & MASK);

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
