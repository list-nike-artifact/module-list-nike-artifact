//! The module list-NIKE of `fig:mod-list-nike`.
//!
//! ```text
//!   Stp(1^lambda):  A <- R_q^{n x m}
//!   Gen(ID):        s <- R_q^n,  e <- CbdVec(eta),  u <- ([-beta,beta]^d)^m
//!                   sk = (s, u)
//!                   pk = (s^T A + e^T,  A u)
//!   Sdk(sk_1, pk_2): k_stat = s^T x,  k_comp = b^T u,  return ListRec(k_stat + k_comp)
//! ```
//!
//! # Where the work goes
//!
//! `Gen` is `2 n m = 1800` ring products: `n m` for `s^T A` and `n m` for `A u`.
//! `Sdk` is `n + m = 109` of them, and then `ListRec`.  Both raw-key halves are
//! accumulated in `i64` and reduced once per output coefficient, so a whole
//! inner product over `R_q` costs one reduction per coefficient rather than
//! `n` or `m` of them.
//!
//! # Which operand is expanded
//!
//! [`crate::poly::mac`] wants one plain and one expanded operand.  `A` is kept
//! expanded, which also gives its plain coefficients for free (the upper half of
//! the expanded array), so `s^T A` and `A u` both read the same table.  `u` is
//! kept expanded because it is the right operand of `b^T u`, and `A u`
//! multiplies it by `A`'s plain half.

use crate::aes::Ctr;
use crate::listrec::{self, Plan};
use crate::params::*;
use crate::poly::{self, Acc, Poly, PolyExp, PolyI32};
use crate::sample;

/// Domain separators for the keystreams.
const NONCE_A: u64 = 0x4100_0000_0000_0000;
const NONCE_S: u64 = 0x5300_0000_0000_0000;
const NONCE_E: u64 = 0x4500_0000_0000_0000;
const NONCE_U: u64 = 0x5500_0000_0000_0000;

/// The public parameter `A in R_q^{n x m}`, row-major, each entry expanded.
pub struct PublicParams {
    pub a: Vec<PolyExp>,
}

impl PublicParams {
    #[inline(always)]
    fn at(&self, i: usize, j: usize) -> &PolyExp {
        &self.a[i * M + j]
    }
}

/// `sk = (s, u)` with `s in R_q^n` and `u in ([-beta,beta]^d)^m`.
pub struct SecretKey {
    /// Plain coefficients of `s`, the left operand of `s^T A` and `s^T x`.
    pub s: Vec<PolyI32>,
    /// `u` in expanded form, the right operand of `A u` and `b^T u`.
    pub u: Vec<PolyExp>,
}

/// `pk = (s^T A + e^T, A u)`.
pub struct PublicKey {
    /// `b^T = s^T A + e^T in R_q^{1 x m}`, plain: it is the left operand of `b^T u`.
    pub l: Vec<PolyI32>,
    /// `x = A u in R_q^n`, expanded: it is the right operand of `s^T x`.
    pub r: Vec<PolyExp>,
}

/// `NIKE.Stp`: samples the public matrix `A`.
pub fn stp(seed: &[u8; 16]) -> PublicParams {
    let mut ctr = Ctr::new(seed, NONCE_A);
    let mut a = Vec::with_capacity(N * M);
    for _ in 0..N * M {
        a.push(sample::uniform_exp(&mut ctr));
    }
    PublicParams { a }
}

/// `NIKE.Gen`: one key pair.  `ID` plays no role in this construction, as the
/// paper notes, so the seed stands in for it.
pub fn gen(pp: &PublicParams, seed: &[u8; 16]) -> (SecretKey, PublicKey) {
    let mut cs = Ctr::new(seed, NONCE_S);
    let mut ce = Ctr::new(seed, NONCE_E);
    let mut cu = Ctr::new(seed, NONCE_U);

    let mut s: Vec<PolyI32> = Vec::with_capacity(N);
    for _ in 0..N {
        let p = sample::uniform(&mut cs);
        s.push(core::array::from_fn(|k| p[k] as i32));
    }
    let e: Vec<PolyI32> = (0..M).map(|_| sample::cbd(&mut ce)).collect();
    let u: Vec<PolyExp> = (0..M).map(|_| poly::expand(&sample::short(&mut cu))).collect();

    // b^T = s^T A + e^T:  m inner products of length n, each reduced once.
    let mut l: Vec<PolyI32> = Vec::with_capacity(M);
    for j in 0..M {
        let mut acc: Acc = poly::acc_zero();
        for i in 0..N {
            poly::mac(&mut acc, &s[i], pp.at(i, j));
        }
        poly::acc_add_centred(&mut acc, &e[j]);
        l.push(poly::to_i32(&poly::reduce_acc(&acc)));
    }

    // x = A u:  n inner products of length m, each reduced once.  Here the
    // right operand is short, so the accumulator stays far below its bound.
    let mut r: Vec<PolyExp> = Vec::with_capacity(N);
    for i in 0..N {
        let mut acc: Acc = poly::acc_zero();
        for j in 0..M {
            poly::mac(&mut acc, poly::plain(pp.at(i, j)), &u[j]);
        }
        r.push(poly::expand_u32(&poly::reduce_acc(&acc)));
    }

    (SecretKey { s, u }, PublicKey { l, r })
}

/// `k_stat = s^T x`: `n` full-width products, so at most `D n q^2 < 2^56.2`.
pub fn k_stat(sk: &SecretKey, pk: &PublicKey) -> Poly {
    let mut acc: Acc = poly::acc_zero();
    for i in 0..N {
        poly::mac(&mut acc, &sk.s[i], &pk.r[i]);
    }
    poly::reduce_acc(&acc)
}

/// `k_comp = b^T u`: `m` products against a short operand, so below `2^41`.
pub fn k_comp(sk: &SecretKey, pk: &PublicKey) -> Poly {
    let mut acc: Acc = poly::acc_zero();
    for j in 0..M {
        poly::mac(&mut acc, &pk.l[j], &sk.u[j]);
    }
    poly::reduce_acc(&acc)
}

/// The raw key `k_stat + k_comp in R_q`, before reconciliation.
pub fn raw_key(sk: &SecretKey, pk: &PublicKey) -> Poly {
    poly::add(&k_stat(sk, pk), &k_comp(sk, pk))
}

/// `NIKE.Sdk`, up to the enumeration: the raw key and its reconciliation plan.
/// The sender needs only [`Plan::first`]; the receiver enumerates.
pub fn sdk_plan(sk: &SecretKey, pk: &PublicKey) -> Plan {
    listrec::plan(&raw_key(sk, pk))
}

/// `NIKE.Sdk` returning only the first list entry, which is what a party that
/// sends a key-confirmation tag needs.
pub fn sdk_first(sk: &SecretKey, pk: &PublicKey) -> listrec::ShKey {
    sdk_plan(sk, pk).first()
}

/// Serialises a public key: `m + n` ring elements at `LOG_Q` bits per
/// coefficient.  The length is whatever that comes to; `bin/check` compares it
/// against the number in the table.
pub fn pk_bytes(pk: &PublicKey) -> Vec<u8> {
    let mut out = Vec::new();
    let mut bits = BitWriter::new(&mut out);
    for p in &pk.l {
        for &c in p.iter() {
            bits.put(c as u32, LOG_Q);
        }
    }
    for p in &pk.r {
        for &c in poly::plain(p).iter() {
            bits.put(c as u32, LOG_Q);
        }
    }
    bits.finish();
    out
}

/// Serialises a secret key: `s` at full width, `u` at `BETA_BITS` per
/// coefficient because its coefficients live in `[-beta, beta]`.
pub fn sk_bytes(sk: &SecretKey) -> Vec<u8> {
    let mut out = Vec::new();
    let mut bits = BitWriter::new(&mut out);
    for p in &sk.s {
        for &c in p.iter() {
            bits.put(c as u32, LOG_Q);
        }
    }
    for p in &sk.u {
        for &c in poly::plain(p).iter() {
            bits.put((c + BETA) as u32, BETA_BITS);
        }
    }
    bits.finish();
    out
}

struct BitWriter<'a> {
    out: &'a mut Vec<u8>,
    acc: u64,
    have: usize,
}

impl<'a> BitWriter<'a> {
    fn new(out: &'a mut Vec<u8>) -> Self {
        BitWriter { out, acc: 0, have: 0 }
    }
    fn put(&mut self, v: u32, bits: usize) {
        self.acc |= (v as u64) << self.have;
        self.have += bits;
        while self.have >= 8 {
            self.out.push(self.acc as u8);
            self.acc >>= 8;
            self.have -= 8;
        }
    }
    fn finish(&mut self) {
        if self.have > 0 {
            self.out.push(self.acc as u8);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair() -> (PublicParams, (SecretKey, PublicKey), (SecretKey, PublicKey)) {
        let pp = stp(&[1u8; 16]);
        let a = gen(&pp, &[2u8; 16]);
        let b = gen(&pp, &[3u8; 16]);
        (pp, a, b)
    }

    /// Serialise real keys and measure; the table's numbers are the expectation.
    #[test]
    fn serialised_sizes_match_the_table() {
        let (_, (ska, pka), _) = pair();
        assert_eq!(pk_bytes(&pka).len(), table::PK_BYTES);
        assert_eq!(sk_bytes(&ska).len(), table::SK_BYTES);
        let k = sdk_first(&ska, &pka);
        assert_eq!(k.len() * 8, table::SHK_BITS);
    }

    /// The exact identity `k_1 - k_2 = e_B^T u_A - e_A^T u_B` of `eq:mod-gap`
    /// says the two raw keys differ by a small amount; hypothesis `eq:mod-H`
    /// says that amount never exceeds `q / 2^(rho+1)`.
    #[test]
    fn raw_key_gap_is_inside_the_reconciliation_radius() {
        let pp = stp(&[9u8; 16]);
        let radius = Q / (1 << (RHO + 1));
        let bound = 2 * (M as u32) * (D as u32) * (ETA as u32) * (BETA as u32);
        assert!(bound <= radius, "eq:mod-H fails: {bound} > {radius}");
        for t in 0..8u8 {
            let (ska, pka) = gen(&pp, &[10 + t; 16]);
            let (skb, pkb) = gen(&pp, &[50 + t; 16]);
            let ka = raw_key(&ska, &pkb);
            let kb = raw_key(&skb, &pka);
            for i in 0..D {
                let g = crate::fq::to_centred(crate::fq::sub(ka[i], kb[i])).unsigned_abs();
                assert!(g <= bound, "gap {g} exceeds the worst case {bound}");
            }
        }
    }

    /// List correctness: `ListRec(k_1)[1]` is in `ListRec(k_2)`.
    #[test]
    fn list_correctness_holds_for_honest_pairs() {
        let pp = stp(&[4u8; 16]);
        for t in 0..24u8 {
            let (ska, pka) = gen(&pp, &[100 + t; 16]);
            let (skb, pkb) = gen(&pp, &[200 + t; 16]);
            let pa = sdk_plan(&ska, &pkb);
            let pb = sdk_plan(&skb, &pka);
            assert!(pb.contains(&pa.first()), "run {t}: A's key not in B's list");
            assert!(pa.contains(&pb.first()), "run {t}: B's key not in A's list");
        }
    }

    /// `s^T A u` computed the two ways must agree, which is the whole reason
    /// both halves are needed.
    #[test]
    fn the_two_halves_share_their_bilinear_part() {
        let pp = stp(&[6u8; 16]);
        let (ska, pka) = gen(&pp, &[7u8; 16]);
        let (skb, pkb) = gen(&pp, &[8u8; 16]);
        // s_A^T (A u_B) via pk_B.r
        let mut x = poly::acc_zero();
        for i in 0..N {
            poly::mac(&mut x, &ska.s[i], &pkb.r[i]);
        }
        // (s_A^T A) u_B via pk_A.l minus the error term
        let mut y = poly::acc_zero();
        for j in 0..M {
            poly::mac(&mut y, &pka.l[j], &skb.u[j]);
        }
        // The difference is e_A^T u_B, which is bounded by D * M * ETA * BETA.
        let xr = poly::reduce_acc(&x);
        let yr = poly::reduce_acc(&y);
        let bound = (M * D * ETA) as u32 * BETA as u32;
        for i in 0..D {
            let g = crate::fq::to_centred(crate::fq::sub(yr[i], xr[i])).unsigned_abs();
            assert!(g <= bound, "coefficient {i}: {g} > {bound}");
        }
    }
}
