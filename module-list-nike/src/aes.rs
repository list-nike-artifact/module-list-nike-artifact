//! AES-128 in counter mode, used as the XOF for all sampling.
//!
//! Follows \textsc{Swoosh}, which expands its public matrix with AES rather
//! than with a sponge; on AArch64 the round instructions make this several
//! times faster than SHAKE.  The NEON path keeps eight blocks in flight.  A
//! portable byte-oriented path exists so the crate builds and tests anywhere.

const SBOX: [u8; 256] = [
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
];

const RCON: [u8; 10] = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36];

/// The eleven AES-128 round keys.
#[derive(Clone)]
pub struct RoundKeys([[u8; 16]; 11]);

/// Expands a 128-bit key into the round keys.
pub fn key_expansion(key: &[u8; 16]) -> RoundKeys {
    let mut w = [[0u8; 16]; 11];
    w[0] = *key;
    for r in 1..11 {
        let prev = w[r - 1];
        let mut t = [prev[13], prev[14], prev[15], prev[12]]; // RotWord
        for b in t.iter_mut() {
            *b = SBOX[*b as usize]; // SubWord
        }
        t[0] ^= RCON[r - 1];
        for c in 0..4 {
            for b in 0..4 {
                let p = if c == 0 { t[b] } else { w[r][4 * (c - 1) + b] };
                w[r][4 * c + b] = prev[4 * c + b] ^ p;
            }
        }
    }
    RoundKeys(w)
}

fn xtime(x: u8) -> u8 {
    (x << 1) ^ (((x >> 7) & 1) * 0x1b)
}

/// Portable single-block encryption, also the oracle for the NEON path.
pub fn encrypt_block_ref(rk: &RoundKeys, block: &[u8; 16]) -> [u8; 16] {
    let mut s = *block;
    for i in 0..16 {
        s[i] ^= rk.0[0][i];
    }
    for r in 1..11 {
        for b in s.iter_mut() {
            *b = SBOX[*b as usize]; // SubBytes
        }
        let t = s;
        for c in 0..4 {
            for row in 0..4 {
                s[4 * c + row] = t[4 * ((c + row) % 4) + row]; // ShiftRows
            }
        }
        if r != 10 {
            let t = s;
            for c in 0..4 {
                let col = &t[4 * c..4 * c + 4];
                let x = col[0] ^ col[1] ^ col[2] ^ col[3];
                for row in 0..4 {
                    s[4 * c + row] = col[row] ^ x ^ xtime(col[row] ^ col[(row + 1) % 4]);
                }
            }
        }
        for i in 0..16 {
            s[i] ^= rk.0[r][i];
        }
    }
    s
}

/// AES-128-CTR keystream generator.
pub struct Ctr {
    rk: RoundKeys,
    #[cfg(target_arch = "aarch64")]
    rk_v: [core::arch::aarch64::uint8x16_t; 11],
    counter: u64,
    nonce: u64,
    buf: [u8; BLOCKS * 16],
    pos: usize,
}

/// Blocks produced per keystream refill.
const BLOCKS: usize = 8;

impl Ctr {
    /// A keystream seeded by `seed`, with `nonce` separating independent streams.
    pub fn new(seed: &[u8; 16], nonce: u64) -> Self {
        let rk = key_expansion(seed);
        #[cfg(target_arch = "aarch64")]
        let rk_v = unsafe {
            let mut v = [core::arch::aarch64::vdupq_n_u8(0); 11];
            for r in 0..11 {
                v[r] = core::arch::aarch64::vld1q_u8(rk.0[r].as_ptr());
            }
            v
        };
        Ctr {
            rk,
            #[cfg(target_arch = "aarch64")]
            rk_v,
            counter: 0,
            nonce,
            buf: [0u8; BLOCKS * 16],
            pos: BLOCKS * 16,
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "aes,neon")]
    unsafe fn refill_neon(&mut self) {
        use core::arch::aarch64::*;
        let mut st = [vdupq_n_u8(0); BLOCKS];
        for b in 0..BLOCKS {
            let mut blk = [0u8; 16];
            blk[0..8].copy_from_slice(&self.nonce.to_le_bytes());
            blk[8..16].copy_from_slice(&(self.counter + b as u64).to_le_bytes());
            st[b] = vld1q_u8(blk.as_ptr());
        }
        for r in 0..9 {
            for b in 0..BLOCKS {
                st[b] = vaesmcq_u8(vaeseq_u8(st[b], self.rk_v[r]));
            }
        }
        for b in 0..BLOCKS {
            st[b] = veorq_u8(vaeseq_u8(st[b], self.rk_v[9]), self.rk_v[10]);
            vst1q_u8(self.buf.as_mut_ptr().add(16 * b), st[b]);
        }
        self.counter += BLOCKS as u64;
        self.pos = 0;
    }

    fn refill_ref(&mut self) {
        for b in 0..BLOCKS {
            let mut blk = [0u8; 16];
            blk[0..8].copy_from_slice(&self.nonce.to_le_bytes());
            blk[8..16].copy_from_slice(&(self.counter + b as u64).to_le_bytes());
            self.buf[16 * b..16 * b + 16].copy_from_slice(&encrypt_block_ref(&self.rk, &blk));
        }
        self.counter += BLOCKS as u64;
        self.pos = 0;
    }

    #[inline(always)]
    fn refill(&mut self) {
        #[cfg(target_arch = "aarch64")]
        {
            if std::arch::is_aarch64_feature_detected!("aes") {
                unsafe { self.refill_neon() };
                return;
            }
        }
        self.refill_ref();
    }

    /// Writes the next `out.len()` keystream bytes.
    pub fn fill(&mut self, out: &mut [u8]) {
        let mut done = 0;
        while done < out.len() {
            if self.pos == self.buf.len() {
                self.refill();
            }
            let take = (out.len() - done).min(self.buf.len() - self.pos);
            out[done..done + take].copy_from_slice(&self.buf[self.pos..self.pos + take]);
            self.pos += take;
            done += take;
        }
    }

    /// The next keystream byte.
    #[inline(always)]
    pub fn byte(&mut self) -> u8 {
        if self.pos == self.buf.len() {
            self.refill();
        }
        let b = self.buf[self.pos];
        self.pos += 1;
        b
    }

    /// The next three keystream bytes as a little-endian integer.
    #[inline(always)]
    pub fn u24(&mut self) -> u32 {
        if self.pos + 3 > self.buf.len() {
            let mut t = [0u8; 3];
            self.fill(&mut t);
            return t[0] as u32 | (t[1] as u32) << 8 | (t[2] as u32) << 16;
        }
        let p = self.pos;
        self.pos += 3;
        self.buf[p] as u32 | (self.buf[p + 1] as u32) << 8 | (self.buf[p + 2] as u32) << 16
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// FIPS-197 Appendix C.1.
    #[test]
    fn fips197_vector() {
        let key: [u8; 16] = [
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
            0x0e, 0x0f,
        ];
        let pt: [u8; 16] = [
            0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd,
            0xee, 0xff,
        ];
        let ct: [u8; 16] = [
            0x69, 0xc4, 0xe0, 0xd8, 0x6a, 0x7b, 0x04, 0x30, 0xd8, 0xcd, 0xb7, 0x80, 0x70, 0xb4,
            0xc5, 0x5a,
        ];
        let rk = key_expansion(&key);
        assert_eq!(encrypt_block_ref(&rk, &pt), ct);
    }

    /// The NEON path must agree with the portable one byte for byte.
    #[test]
    fn neon_matches_reference() {
        let seed = [0x5au8; 16];
        let mut fast = Ctr::new(&seed, 0xdead_beef);
        let mut slow = Ctr::new(&seed, 0xdead_beef);
        let mut a = [0u8; 4096];
        let mut b = [0u8; 4096];
        fast.fill(&mut a);
        for _ in 0..(4096 / (BLOCKS * 16)) {
            slow.refill_ref();
            let off = slow.counter as usize * 16 - BLOCKS * 16;
            b[off..off + BLOCKS * 16].copy_from_slice(&slow.buf);
        }
        assert_eq!(a, b);
    }

    #[test]
    fn byte_and_fill_agree() {
        let seed = [7u8; 16];
        let mut x = Ctr::new(&seed, 1);
        let mut y = Ctr::new(&seed, 1);
        let mut buf = [0u8; 300];
        x.fill(&mut buf);
        for i in 0..300 {
            assert_eq!(buf[i], y.byte(), "byte {i}");
        }
    }
}
