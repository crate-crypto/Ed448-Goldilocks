use crate::curve::constants::EDWARDS_D;
use crate::field::Fq;

#[derive(Eq, Debug, PartialEq)]
pub struct Montgomery {
    z0: Fq,
    zd: Fq,
    za: Fq,
    xa: Fq,
    xd: Fq,
}

impl Montgomery {
    /// Montgomery Step
    pub fn step(&mut self) {
        let mut l0 = self.zd.add_no_reduce(&self.xd);
        let mut l1 = sub(&self.xd, &self.zd);
        self.zd = sub(&self.xa, &self.za);
        self.xd = l0 * self.zd;
        self.zd = self.za.add_no_reduce(&self.xa);
        self.za = l1 * self.zd;
        self.xa = self.za.add_no_reduce(&self.xd);
        self.zd = self.xa.square();
        self.xa = self.z0 * self.zd;
        self.zd = sub(&self.xd, &self.za);
        self.za = self.zd.square();
        self.xd = l0.square();
        l0 = l1.square();
        self.zd = Fq::vartime_mul_sdword_curve_constant(&self.xd, 1 - EDWARDS_D);
        l1 = sub(&self.xd, &l0);
        self.xd = l0 * self.zd;
        l0 = self.zd.sub_no_reduce(&l1);
        l0.bias(2);
        l0.weak_reduce();
        self.zd = l0 * l1;
    }
}

fn sub(x: &Fq, y: &Fq) -> Fq {
    let mut a = x.sub_no_reduce(y);
    a.bias(2);
    a.weak_reduce();
    a
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::num::ParseIntError;
    fn decode_hex(s: &str) -> Result<Vec<u8>, ParseIntError> {
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16))
            .collect()
    }

    fn slice_to_fixed_array(b: &[u8]) -> [u8; 56] {
        let mut a: [u8; 56] = [0; 56];
        a.copy_from_slice(&b);
        a
    }

    fn hex_to_fq(hex: &str) -> Fq {
        let mut bytes = decode_hex(hex).unwrap();
        bytes.reverse();
        Fq::from_bytes(&slice_to_fixed_array(&bytes))
    }

    #[test]
    fn test_step() {
        let z0 = hex_to_fq("e281b05e4051a52b331430897d9d950529a46637d3ca1f45e1d2dc4fbd164c956f25dd0cf30458b4129e900faa2ba9b8d305dc4ae1e1b343");
        let xd = hex_to_fq("dc7c2264cf2a3f6178ee7884793f2d0cfe98e602c32adfbec9a5fc225c904f5e1f45c614fc483aec252745e04a38f49e1a4cfc0e8bbf14c5");
        let zd = hex_to_fq("76ad01dbfd7dd72671ad1f827b762fe0c39c808084533b1e22ee18537b7e43c75b995f9e107ec055fbb3df4fb83ad78e69de76a188fb6db6");
        let xa = hex_to_fq("86032c9f990e2680726003f62a1ec5c01f18ad130ce0883b247d2ea9e8d591e6121e6007027d44d94d9659a05fb47e91c3c11b5552cb2185");
        let za = hex_to_fq("5c157ed45b60be2db18a494780b5b7ab79ae1afc3919c9b00c1879495ea079b73990eebf5f0def2897fe8ca78084c07ef89c5bfc336625fd");
        let mut montgomery = Montgomery { z0, xd, zd, xa, za };
        montgomery.step();
        let expected_z0 = hex_to_fq("e281b05e4051a52b331430897d9d950529a46637d3ca1f45e1d2dc4fbd164c956f25dd0cf30458b4129e900faa2ba9b8d305dc4ae1e1b343");
        let expected_xd = hex_to_fq("c88f896abf42ca2cbff1edf881d1246ee76abe7385932d7b54fb9d71307fdd8043d8a80c7d0363e7a45443d4e9a03bf3e0aab82fb4714c5f");
        let expected_zd = hex_to_fq("962fa8b019eeedd607eda6b44454e17b76b1536f6b336362257d72c3c1576339514f1f4d2d0ae7b0680469a432a2f54cb7f9dbc14473802d");
        let expected_xa = hex_to_fq("09e41fe2e74667a6676fb0492b496f7d69d45055601ec86839b95e9343407ed592ea357118e5568eea272e9349adf0efbe29307187cfff6e");
        let expected_za = hex_to_fq("b115a615745fc6f453a43d1466e12acd2215ac373cadcd633211235510c6a04c4f041006d07f543f2bd4b050ecdd472be4415ab7a3f79f95");

        let expected_montgomery = Montgomery {
            z0: expected_z0,
            xd: expected_xd,
            zd: expected_zd,
            xa: expected_xa,
            za: expected_za,
        };

        assert_eq!(montgomery, expected_montgomery);
    }
}
