use crate::curve::edwards::projective::ProjectiveNielsPoint;
use crate::curve::edwards::ExtendedPoint;

pub(crate) const WINDOW: usize = 5;
pub(crate) const WINDOW_MASK: usize = (1 << WINDOW) - 1;
pub(crate) const WINDOW_T_MASK: usize = WINDOW_MASK >> 1;
pub(crate) const TABLE_SIZE: usize = 16;

pub fn prepare_fixed_window(point: &ExtendedPoint) -> [ProjectiveNielsPoint; TABLE_SIZE] {
    let P2 = point.double_internal(false).to_projective_niels();

    let mut table = [point.to_projective_niels(); TABLE_SIZE];

    let mut p_original = point.clone();
    for i in 1..TABLE_SIZE {
        p_original = p_original.add_projective_niels(&P2, false);
        table[i] = p_original.to_projective_niels();
    }
    table
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fq;
    #[test]
    fn test_prepare_fixed_window() {
        let X = Fq([
            0x0e0fbf9e, 0x0ba1bcd7, 0x01cc6d39, 0x053b56e8, 0x0635d142, 0x0383307a, 0x0f8a159b,
            0x097fd2cf, 0x0fa310f6, 0x05522bde, 0x0b981703, 0x0b095b1e, 0x042d4780, 0x05ae11df,
            0x0934fe80, 0x0dc6474d,
        ]);
        let Y = Fq([
            0x02c1149c, 0x0e72febf, 0x05259893, 0x0723e184, 0x0f7232ff, 0x019a5600, 0x05581d2c,
            0x07331444, 0x04e0124a, 0x09c3c5e5, 0x0945536e, 0x0b786a20, 0x0f75623f, 0x00ba30e8,
            0x0cc589a3, 0x04a2eea8,
        ]);

        let Z = Fq([
            0x02406c71, 0x0b2fdb67, 0x02591aa2, 0x085fc24e, 0x0dc50d09, 0x08692c5b, 0x0ba917d7,
            0x0aefea74, 0x037d0084, 0x04d5defa, 0x08bbe7ad, 0x050da977, 0x08adf827, 0x05425cdd,
            0x037d816d, 0x0d59cd0a,
        ]);

        let T = Fq([
            0x0baf8c30, 0x06686ad3, 0x0c149bac, 0x0f57f68d, 0x05cd321a, 0x02ff8d60, 0x09dcc4bd,
            0x0f731ec2, 0x0cd7ea75, 0x0be970e4, 0x043d30e0, 0x0dd64b9b, 0x04f78bf1, 0x0d1fde20,
            0x05c88e97, 0x026ce314,
        ]);

        let p = ExtendedPoint { X, Y, Z, T };
        let w = prepare_fixed_window(&p);
        assert_eq!(w.len(), TABLE_SIZE);

        let expected_window = [
            // 0
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x04b154fd, 0x02d141e7, 0x03592b5a, 0x01e88a9c, 0x093c61bd, 0x0e172586,
                    0x05ce0790, 0x0db34174, 0x053d0152, 0x04719a06, 0x0dad3c6b, 0x006f0f01,
                    0x0b481abf, 0x0b0c1f09, 0x03908b22, 0x06dca75b,
                ]),
                y_plus_x: Fq([
                    0x00d0d43b, 0x0a14bb97, 0x06f205cd, 0x0c5f386c, 0x05a80441, 0x051d867b,
                    0x04e232c7, 0x00b2e714, 0x04832342, 0x0f15f1c4, 0x04dd6a71, 0x0681c53f,
                    0x03a2a9c0, 0x066842c8, 0x05fa8823, 0x026935f6,
                ]),
                td: Fq([
                    0x0460a1f7, 0x0e76b0c4, 0x0bc48547, 0x0a643633, 0x0ff970aa, 0x0cb5cdc9,
                    0x0d2a0bc4, 0x06940a23, 0x0ad0577f, 0x07e65c18, 0x04b0332f, 0x059b353f,
                    0x010aebde, 0x09e69eb7, 0x084e54ff, 0x09ba3b12,
                ]),
                Z: Fq([
                    0x0480d8e3, 0x065fb6ce, 0x04b23545, 0x00bf849c, 0x0b8a1a13, 0x00d258b7,
                    0x07522faf, 0x05dfd4e9, 0x06fa010a, 0x09abbdf4, 0x0177cf5a, 0x0a1b52ef,
                    0x015bf04e, 0x0a84b9bb, 0x06fb02da, 0x0ab39a14,
                ]),
            },
            // 1
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x00a553ee, 0x00234003, 0x055c092c, 0x057bcbdd, 0x0e495a8c, 0x09c8b997,
                    0x0c143649, 0x0399ba66, 0x0df0f55f, 0x0c004b24, 0x0a7c0ab5, 0x0a95f91a,
                    0x005ee2fc, 0x011bb28d, 0x08ca5e5a, 0x0192776c,
                ]),
                y_plus_x: Fq([
                    0x0a2f247e, 0x02e7847f, 0x0a1e65d9, 0x0571128d, 0x0c710f20, 0x0d2d073c,
                    0x02f89f06, 0x01a67d7c, 0x0f1e3c52, 0x0e90f285, 0x0a97b2e7, 0x00f68012,
                    0x005f2af2, 0x0ddcb5d2, 0x0bd8a372, 0x0c15a881,
                ]),
                td: Fq([
                    0x0c7957ba, 0x02b7715c, 0x0951bac2, 0x08ceafd5, 0x0de55a10, 0x01edb51f,
                    0x019e372a, 0x0ef9da41, 0x00e4539d, 0x08f12616, 0x072899db, 0x04eabccc,
                    0x01c30552, 0x03c77348, 0x02e95364, 0x0d406ee6,
                ]),
                Z: Fq([
                    0x0f2830d8, 0x0056a35b, 0x08ce98ce, 0x0523ef98, 0x0922d7c0, 0x0dea9b3d,
                    0x04016a54, 0x0f7d6e3f, 0x0596053c, 0x0f26db88, 0x0e066eea, 0x0866a255,
                    0x0e8a62a7, 0x0dffb915, 0x084fa2b2, 0x04b35e84,
                ]),
            },
            // 2
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x0bdc0118, 0x02f3217f, 0x09ebc3b3, 0x07617f27, 0x0434daf9, 0x05c63130,
                    0x079ffbba, 0x0ed66a9f, 0x0ac6adc7, 0x0cea22df, 0x0ead4d20, 0x04d91dbe,
                    0x08aea04f, 0x0eff9bba, 0x08c6636d, 0x08fbf2fe,
                ]),
                y_plus_x: Fq([
                    0x0d6dae4e, 0x042ef7e0, 0x0e1974ba, 0x04c7f440, 0x0cc3f843, 0x00083ff4,
                    0x0490f26f, 0x0e95b6b5, 0x04b546d1, 0x059c8373, 0x03c0f841, 0x02ff119c,
                    0x087b2772, 0x0fc993a1, 0x02ddd9b6, 0x04a9c14d,
                ]),
                td: Fq([
                    0x0b9dbcd2, 0x0156cdf2, 0x02889448, 0x067306ba, 0x0bb5be76, 0x090eb070,
                    0x0a50ce79, 0x09b43e73, 0x080845ed, 0x09a2ecf7, 0x018e594a, 0x0b241a06,
                    0x051f02da, 0x031bbca7, 0x09b2848e, 0x0bbb2ec8,
                ]),
                Z: Fq([
                    0x08078601, 0x088d1850, 0x00e0769d, 0x08eb236c, 0x0130e14f, 0x09e4586f,
                    0x0911c3b0, 0x0ce7530f, 0x0c10cc39, 0x0da28bc9, 0x05845355, 0x0003ec5a,
                    0x07b7bf88, 0x0300bf50, 0x0653991f, 0x049db070,
                ]),
            },
            // 3
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x00496408, 0x0b651531, 0x067fd460, 0x085a3886, 0x04d21048, 0x0d621139,
                    0x071d8a8d, 0x03199572, 0x0e03e14c, 0x0d05c021, 0x0833d2b3, 0x0ef5a2c1,
                    0x0e1d24c6, 0x0027b3ae, 0x02ce6c8e, 0x08402aa3,
                ]),
                y_plus_x: Fq([
                    0x05ad9af1, 0x0a0c1d97, 0x0b2d81fb, 0x01ea8e2d, 0x0a22652b, 0x0efd1612,
                    0x082ba3ab, 0x0afad623, 0x0c37c3dc, 0x094a087f, 0x09442f18, 0x0813269f,
                    0x0a734c8d, 0x063a6cb0, 0x0010d46a, 0x06638790,
                ]),
                td: Fq([
                    0x0689ba7e, 0x081e4dc6, 0x0c1ef0d4, 0x068f65c6, 0x03cb2bd2, 0x09c5100c,
                    0x011feeda, 0x0f3d4306, 0x00af4e32, 0x0358eba7, 0x07994455, 0x03534ba3,
                    0x0ffa2ea6, 0x0621d409, 0x0b854ef5, 0x06997a18,
                ]),
                Z: Fq([
                    0x0a3dc259, 0x04b1a065, 0x03a2fa7c, 0x0099ce98, 0x0434976c, 0x04a0786d,
                    0x0e76a8b9, 0x0d291d23, 0x01dc0175, 0x05d1717b, 0x0c1d478a, 0x04ba7edf,
                    0x0b4b4324, 0x0c5f36a8, 0x0ed933aa, 0x0ac6d3f9,
                ]),
            },
            // 4
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x02df046b, 0x057e1a14, 0x0005dcd2, 0x035667ff, 0x0f83f892, 0x0d6b9dca,
                    0x0f8fb738, 0x025946da, 0x08263015, 0x01e7865e, 0x0c2cba18, 0x0e10546f,
                    0x0f091335, 0x0df8e689, 0x010a9906, 0x0830b70c,
                ]),
                y_plus_x: Fq([
                    0x0115bbe7, 0x0867aaae, 0x0ca910ef, 0x0c0ae406, 0x0a5a57ce, 0x03c2d826,
                    0x0457e63e, 0x0f34b795, 0x048c73f6, 0x0eb2313e, 0x0678f5b1, 0x05a6dbef,
                    0x0c78ba56, 0x094e3e5c, 0x0da22f13, 0x0678ecb3,
                ]),
                td: Fq([
                    0x0797d90b, 0x023eeee5, 0x06091221, 0x0ececc3a, 0x0ec1cca5, 0x041b1079,
                    0x0b4f4a87, 0x0a84bb5c, 0x0976bd53, 0x074c939e, 0x0b95416c, 0x05ea1867,
                    0x08ebff05, 0x0bc870f1, 0x0b309de3, 0x0010c34e,
                ]),
                Z: Fq([
                    0x0b38f0f1, 0x03bc6913, 0x061ba597, 0x0cf921ec, 0x0bd4a2a5, 0x02fac843,
                    0x083f5c15, 0x07d7877a, 0x09dfebf5, 0x015ee5a1, 0x09beba6e, 0x0a89a826,
                    0x08e66f3a, 0x063d031e, 0x09d450f7, 0x0a90a158,
                ]),
            },
            // 5
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x00235a46, 0x04d74adc, 0x0ddadd18, 0x061afbc5, 0x00b8d5e3, 0x04ad265e,
                    0x0b9914b0, 0x0e78fbb0, 0x09908a36, 0x0a2dde9e, 0x07ac5396, 0x0af6ff6a,
                    0x068ed7eb, 0x0dc11e8b, 0x02fec204, 0x0d5389db,
                ]),
                y_plus_x: Fq([
                    0x08ae1125, 0x0555a7b2, 0x0b489173, 0x07472588, 0x0a70d04c, 0x028ff109,
                    0x058ef8ab, 0x065f78fa, 0x0c88caf4, 0x026957ca, 0x02165c81, 0x05288d3d,
                    0x05ae2148, 0x0108fdde, 0x0d31cc57, 0x048b9c53,
                ]),
                td: Fq([
                    0x0be9e6a9, 0x0617ae56, 0x0759b26d, 0x0c3c93e8, 0x0779d05d, 0x0789edb3,
                    0x081cf9dd, 0x09aa41bb, 0x07e8f870, 0x054a61b5, 0x088f337d, 0x0bb738c2,
                    0x06bdc816, 0x0070143c, 0x0114a075, 0x08f49ba9,
                ]),
                Z: Fq([
                    0x048118f9, 0x0ec006fe, 0x05980ad8, 0x0e39107c, 0x0733036a, 0x0f5fc159,
                    0x07af1bad, 0x0e465b51, 0x0750076c, 0x0dd11964, 0x0b35aa73, 0x0a1fe8ae,
                    0x0294748c, 0x0f5150d7, 0x06ade020, 0x0659ffca,
                ]),
            },
            // 6
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x04e37267, 0x0df7e646, 0x002c1f13, 0x07145a47, 0x07f9da71, 0x0ef38604,
                    0x03b74cf6, 0x04f98272, 0x0b1e5de6, 0x0961e851, 0x066a9c11, 0x0637cd49,
                    0x0c40ca9e, 0x0ad4f30c, 0x0eddbe48, 0x0235fe80,
                ]),
                y_plus_x: Fq([
                    0x095ecd88, 0x007afcf3, 0x04645e96, 0x0929d281, 0x0c236979, 0x01487c7b,
                    0x01c93417, 0x0070157a, 0x000f5888, 0x0cd5c842, 0x09578536, 0x06bd9438,
                    0x0dde9e3f, 0x09d37fdc, 0x0c425dac, 0x081a7c61,
                ]),
                td: Fq([
                    0x0a01e5ed, 0x04ac3269, 0x0246ed1d, 0x01184f59, 0x086c9b60, 0x0b696eae,
                    0x0e544d50, 0x056e5ef4, 0x0e69db47, 0x0df87f7d, 0x08e76db7, 0x0f048dd2,
                    0x041b6fe1, 0x0d26dde7, 0x09177799, 0x0e5ebd77,
                ]),
                Z: Fq([
                    0x0f762e2c, 0x0e3d5341, 0x05ab9ead, 0x081acf6b, 0x068db9ea, 0x08cde1f2,
                    0x08cfc97f, 0x0df58221, 0x035f081e, 0x08b0a279, 0x0ee831ac, 0x0c40fd59,
                    0x07ae3f18, 0x086b3bed, 0x0ffcb0f3, 0x08492cfb,
                ]),
            },
            // 7
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x06baca70, 0x098c76c4, 0x06d43244, 0x063a0112, 0x0c712d42, 0x0c44c636,
                    0x0f7f2230, 0x040175a2, 0x05202d65, 0x08d900d7, 0x01da615d, 0x0894103b,
                    0x02fc86c6, 0x0b908128, 0x0ad5f030, 0x0a5d8188,
                ]),
                y_plus_x: Fq([
                    0x09c62e7a, 0x0194b55b, 0x00924732, 0x00baca1e, 0x0faab36a, 0x07815f7e,
                    0x0a75c250, 0x0e63bfa9, 0x03260e43, 0x0d0e684c, 0x06f30b44, 0x01e5fb64,
                    0x0de6297f, 0x03e7ae5a, 0x0fd9af92, 0x01d12c86,
                ]),
                td: Fq([
                    0x0076d1d3, 0x0e083ce0, 0x0140855d, 0x0efb8ca5, 0x0cff83d8, 0x0f899ccb,
                    0x0c461c57, 0x02bc586b, 0x0eb8ebf5, 0x03f8256a, 0x071b8d04, 0x0e48f5cb,
                    0x08a14607, 0x03882433, 0x0c2a92b6, 0x0562a977,
                ]),
                Z: Fq([
                    0x0028918f, 0x0f0c29c3, 0x0e1376f6, 0x09384b43, 0x0e8fc757, 0x00e5805c,
                    0x0455c0de, 0x055a162e, 0x0c4fd69a, 0x02330659, 0x061e4124, 0x094618ab,
                    0x040ef9be, 0x0a8f1974, 0x0b099f5c, 0x045c4715,
                ]),
            },
            // 8
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x05a5b500, 0x0985380b, 0x03849e60, 0x06582f04, 0x013402ba, 0x07c4a15b,
                    0x04afb6e5, 0x019465a6, 0x035d0446, 0x0df3f9e4, 0x050022f4, 0x08332f10,
                    0x082f2c46, 0x0f810973, 0x0e0216de, 0x0ebdd61f,
                ]),
                y_plus_x: Fq([
                    0x03aadbdc, 0x0f166b7e, 0x0a7619b8, 0x02528842, 0x008d0c13, 0x0f12020f,
                    0x06dd379b, 0x0736d3f4, 0x0ab1919d, 0x025ea6c2, 0x09bdc7e6, 0x0528572b,
                    0x0b5a9a26, 0x0f7e2e75, 0x0b326b12, 0x0ed1cd91,
                ]),
                td: Fq([
                    0x098e28fb, 0x03895751, 0x02d65b85, 0x077bc10a, 0x07b9a601, 0x066b9a42,
                    0x0dfabb87, 0x01bee158, 0x09b2d4e5, 0x06919086, 0x0169f374, 0x0942eb6b,
                    0x0bbe91e3, 0x047d5b20, 0x037523b6, 0x0473eb88,
                ]),
                Z: Fq([
                    0x077ca5f0, 0x0fdbec93, 0x0bf4673c, 0x0bcd691d, 0x0b4585cd, 0x0c082eb9,
                    0x03c5a97b, 0x058b3645, 0x00ef8e34, 0x0424ca79, 0x004655e4, 0x0d4a8fad,
                    0x08af0bb8, 0x000e75a7, 0x0b6bbaff, 0x0d40c5f3,
                ]),
            },
            // 9
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x06773538, 0x002e55cf, 0x09144496, 0x027e21a0, 0x0fa4fdbf, 0x01568a77,
                    0x087fe688, 0x05c6d9e6, 0x0508158c, 0x0b7433da, 0x080184db, 0x007477de,
                    0x00e3b950, 0x0f5b29b8, 0x0dafd838, 0x00336b81,
                ]),
                y_plus_x: Fq([
                    0x0faf83af, 0x08b8eac0, 0x045af9b9, 0x00fd0239, 0x025ab82f, 0x0a5a46c4,
                    0x005ab02b, 0x0548c499, 0x0af1167f, 0x0e3ff944, 0x0bbc8c50, 0x08df9171,
                    0x08421dcf, 0x05354fb7, 0x027fa656, 0x0d88cce5,
                ]),
                td: Fq([
                    0x06592c5f, 0x0758bc19, 0x0976a5c4, 0x080e3e50, 0x05637d60, 0x050266c4,
                    0x06cce0c2, 0x000f76ba, 0x0b39952b, 0x08118c1c, 0x0f903a05, 0x08dbaac8,
                    0x02bb60d0, 0x082a0d3a, 0x06b73845, 0x00f848b7,
                ]),
                Z: Fq([
                    0x0021161c, 0x0b667c0e, 0x08a9dc8c, 0x086fa0b9, 0x027b84c9, 0x0de272b8,
                    0x07dad62c, 0x0697d381, 0x03f9ac5a, 0x06b74d3b, 0x0416858d, 0x048dc187,
                    0x03cd18b2, 0x040eca1a, 0x0c49f066, 0x0ff05257,
                ]),
            },
            // 10
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x04992354, 0x0f624b93, 0x06c88150, 0x07f72b38, 0x02b76818, 0x01e09970,
                    0x07856dfb, 0x0b14c688, 0x02d8e364, 0x0d51048e, 0x0674d538, 0x0b822562,
                    0x0e4e17f9, 0x0b6bf11b, 0x0c10ee26, 0x070319cf,
                ]),
                y_plus_x: Fq([
                    0x036f3c2c, 0x054cefa5, 0x0cab8002, 0x0d1568b8, 0x0add5fe7, 0x0071a252,
                    0x04b3c490, 0x01e15013, 0x000e67bf, 0x01926c30, 0x0dc22cc5, 0x0d245c14,
                    0x03ed8abe, 0x0a7290c8, 0x058d2b11, 0x08aa1db2,
                ]),
                td: Fq([
                    0x0b43c9d5, 0x0d79d7fe, 0x040a1d18, 0x0a6369f1, 0x053e0d48, 0x09f1d213,
                    0x0bd51750, 0x0005ec64, 0x0297d428, 0x0acf4828, 0x090dd0e3, 0x0409de55,
                    0x0965d34c, 0x00b94f3a, 0x01dc5637, 0x085b74a4,
                ]),
                Z: Fq([
                    0x09ffe991, 0x04a4ae6a, 0x0851b2cb, 0x0dd241df, 0x0cb58ff7, 0x0e6f6cec,
                    0x00da8c9b, 0x09032938, 0x01a53998, 0x0d5fcd05, 0x066c843a, 0x0f153cee,
                    0x01ab32a1, 0x0abd6dc5, 0x06a28c1f, 0x07c2b400,
                ]),
            },
            // 11
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x054fad63, 0x0a6deab8, 0x03139dab, 0x08ec509b, 0x0f38948e, 0x02c721e1,
                    0x08dd641a, 0x0252cc09, 0x0c443f31, 0x0c8aaa66, 0x0c6458c1, 0x033b5b8b,
                    0x0089236e, 0x07b3592b, 0x09bf422c, 0x0261e002,
                ]),
                y_plus_x: Fq([
                    0x05817da6, 0x02cbe159, 0x0253d5b3, 0x004b79fb, 0x075cbe56, 0x059c447b,
                    0x00211c1b, 0x096e941d, 0x004f1774, 0x0c4b650a, 0x06cc4142, 0x018ce550,
                    0x01754ebb, 0x0afd200c, 0x0440c27e, 0x082d3511,
                ]),
                td: Fq([
                    0x08151172, 0x0bd650e1, 0x094aa489, 0x05afdd3a, 0x03e3c4a4, 0x0c54bade,
                    0x030e6eb2, 0x00ffc2c0, 0x0cd3eaa0, 0x0bab2965, 0x05e31ca3, 0x07dbd978,
                    0x01e50070, 0x01a31a70, 0x0e187e72, 0x04a6706d,
                ]),
                Z: Fq([
                    0x0a64587c, 0x09594b13, 0x07d06590, 0x0b62e1a6, 0x07a498ae, 0x040635bd,
                    0x04197574, 0x091c0d95, 0x09847940, 0x0be452a4, 0x0e2db46d, 0x077230bf,
                    0x0a8a3be9, 0x003e2ead, 0x090c8266, 0x01366b4e,
                ]),
            },
            // 12
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x0a3a6958, 0x0b7c6460, 0x0a26894a, 0x0a7f749b, 0x0d68a071, 0x08f5f7aa,
                    0x0129f200, 0x0be007b7, 0x0f8b0279, 0x0a1bc934, 0x0d0ece01, 0x00317338,
                    0x0f371340, 0x0f79d31d, 0x060ace29, 0x0ef4c1da,
                ]),
                y_plus_x: Fq([
                    0x0c1d0072, 0x0403e959, 0x0ec34cd9, 0x0c859e5e, 0x05fc7843, 0x0d850a1c,
                    0x0a50ccf9, 0x02dba347, 0x00d88f27, 0x08b61a79, 0x06a5f7e4, 0x0fe2f470,
                    0x0ecfb643, 0x05eae6fc, 0x0eeb0971, 0x0713d0d7,
                ]),
                td: Fq([
                    0x0db37a3f, 0x07c85d50, 0x04651c14, 0x084a7a43, 0x0c624b75, 0x09db7e57,
                    0x05b93aaf, 0x04aacb8d, 0x014b7e2c, 0x099cc5fe, 0x06c2374a, 0x04047476,
                    0x0c1b46d1, 0x073aed98, 0x0cb2a6e5, 0x045903bd,
                ]),
                Z: Fq([
                    0x0e1fcb85, 0x0ef06397, 0x01551f6a, 0x07ae3bcb, 0x08a026f5, 0x02398ab0,
                    0x0e39a35b, 0x05a9b180, 0x0df3cd2b, 0x0e7a9f02, 0x00d7f4b3, 0x070d723d,
                    0x0b84e97c, 0x0b1fa785, 0x0ba13ca1, 0x0a266cb0,
                ]),
            },
            // 13
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x0ec3db6b, 0x065f580b, 0x04fe9402, 0x07296ffe, 0x0a976860, 0x044b9376,
                    0x0060db4c, 0x08e2af69, 0x0a3c438a, 0x034a88bb, 0x06dd9552, 0x04c6d9ff,
                    0x0f880fef, 0x0c19a2ac, 0x003b6ea1, 0x0c5bc3a3,
                ]),
                y_plus_x: Fq([
                    0x026a4750, 0x0715decb, 0x0e648646, 0x0a122697, 0x0d0de5bf, 0x0b450ecb,
                    0x0eb9f9d6, 0x08d7a7cf, 0x04b5e294, 0x01a50226, 0x0eee5b89, 0x05b3390a,
                    0x0b64958e, 0x037eebed, 0x0041d7aa, 0x0b6ccce2,
                ]),
                td: Fq([
                    0x0824e359, 0x0a93715f, 0x043c899c, 0x0ac52d1d, 0x01c6c0d5, 0x0c28ad2d,
                    0x075514b9, 0x009a03c1, 0x00852ff9, 0x035a06ff, 0x02461b9c, 0x0451dfe8,
                    0x077838ab, 0x0d925573, 0x08016fb9, 0x0fe722ef,
                ]),
                Z: Fq([
                    0x0e50b079, 0x026c76da, 0x026b4fb0, 0x0661120a, 0x0d766626, 0x00f79967,
                    0x0fe99490, 0x0699c164, 0x085a39a8, 0x0ea7b9f1, 0x07103c79, 0x05046422,
                    0x06fe479f, 0x01827469, 0x044ee87b, 0x0cf58841,
                ]),
            },
            // 14
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x0ebfe2a0, 0x0eb5fbb2, 0x0dc6ec15, 0x0668e6e7, 0x0a2b46d6, 0x09716cff,
                    0x0b1f6161, 0x075f2e4b, 0x092bafec, 0x04a09fe7, 0x0d0693b5, 0x0181e9e6,
                    0x0ae2011c, 0x0339acaa, 0x05ab0851, 0x02e7f480,
                ]),
                y_plus_x: Fq([
                    0x00923994, 0x09e0df3c, 0x05d44abd, 0x0e120043, 0x0802cd52, 0x03d6dc68,
                    0x06848872, 0x0b6dc339, 0x0847737b, 0x0bc369ba, 0x02096051, 0x0bb48580,
                    0x07f1feb0, 0x01f3325e, 0x05347dc3, 0x0ac3e06b,
                ]),
                td: Fq([
                    0x0efcb0d5, 0x01516f87, 0x00c105be, 0x0c6c761b, 0x0b0e62e2, 0x0a170de3,
                    0x07e1d862, 0x0efee135, 0x0906ff9a, 0x097301d4, 0x0746a0ee, 0x0ebb2175,
                    0x0bfdd54e, 0x07828235, 0x0ed64735, 0x03c5ab55,
                ]),
                Z: Fq([
                    0x0cc630cf, 0x04be9e7e, 0x02806f73, 0x033000c5, 0x08901863, 0x06198942,
                    0x057694d2, 0x0c4ea1d8, 0x0c5c8ca4, 0x0ba40152, 0x0c694dd1, 0x085ae4e6,
                    0x05bc4f74, 0x0db70b9f, 0x09da2b77, 0x0a6633ef,
                ]),
            },
            // 15
            ProjectiveNielsPoint {
                y_minus_x: Fq([
                    0x077fc122, 0x08edd3d4, 0x0d81d8d0, 0x01ac9143, 0x08e4b336, 0x0f81211a,
                    0x0d1e5694, 0x008bd901, 0x04116e1f, 0x0fb16788, 0x092dd842, 0x0cbe52bb,
                    0x0ca568ee, 0x0f75e6bb, 0x03efc3bb, 0x00e50e2e,
                ]),
                y_plus_x: Fq([
                    0x05fbceed, 0x0eb215cb, 0x05eeff9f, 0x083a37a3, 0x0044f669, 0x0abb2e46,
                    0x0807c9b6, 0x0201d31e, 0x0ef4ff9c, 0x07157790, 0x087931f8, 0x0b23a524,
                    0x0b4c9a95, 0x012ae723, 0x0aa48f6e, 0x08f0775d,
                ]),
                td: Fq([
                    0x01160f9b, 0x08b66bfe, 0x0d9b8296, 0x0a7bbf2a, 0x0c1a10af, 0x0f3cc344,
                    0x02c0089c, 0x0f4d630f, 0x07901e60, 0x0c9ce1e7, 0x02f5380a, 0x0697ba06,
                    0x0e6065cf, 0x0e709b4e, 0x0de4fa1d, 0x043caab1,
                ]),
                Z: Fq([
                    0x012a04c3, 0x002a0678, 0x015d16e6, 0x04ec95aa, 0x0072c656, 0x01fee9b9,
                    0x0ea90f88, 0x04c44954, 0x0eada4ba, 0x0268f035, 0x048efc9e, 0x0905d82f,
                    0x0e8c31fc, 0x0f80f8f9, 0x0549be7b, 0x066b5184,
                ]),
            },
        ];

        for i in 0..TABLE_SIZE {
            let expected_i = expected_window[i];
            let w_i = w[i];

            assert!(expected_i.Z.equals(&w_i.Z));
            assert!(expected_i.y_plus_x.equals(&w_i.y_plus_x));
            assert!(expected_i.y_minus_x.equals(&w_i.y_minus_x));
            assert!(expected_i.td.equals(&w_i.td));
        }
    }
}
