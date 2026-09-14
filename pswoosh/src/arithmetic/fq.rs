//! Arithmetic in F_q for q = 2^214 - 255, ported to portable Rust.
//!
//! Upstream (`pswoosh/rust/ref0/src/arithmetic/fq.rs`) declares five functions
//! `extern` and links them against `fq.s`, 2718 lines of x86-64 assembly emitted
//! by the Jasmin compiler from libjbn.  This file replaces exactly those five,
//! `fp_add`, `fp_sub`, `fp_mul`, `fp_toM` and `fp_fromM`, with textbook 4-limb Montgomery arithmetic over `u64`/`u128`, so that the
//! scheme builds and runs natively on aarch64.  Everything else in this file
//! (`cmp`, `elem_frombytes`, `elem_tobytes`, the tests) is upstream's, unchanged.
//!
//! The Montgomery radix is R = 2^256, which is what the upstream assembly uses:
//! its own `.data` section stores R mod q and R^2 mod q, and both agree with the
//! constants below.  `ZETAS[0]` of `zetas.rs`, the Montgomery representation of
//! 1, is likewise R mod q.

pub const NLIMBS: usize = 4;
pub type Elem = [u64; NLIMBS];
const K: usize = 214; // bit size of q
pub const ELEM_BYTES: usize = K/8+1;
const RAD: usize = 64; //radix

/* Q = 2^214-255 */   /* 2^0              2^64               2^128              2^192   */
pub const Q: Elem =   [0xffffffffffffff01,0xffffffffffffffff,0xffffffffffffffff,0x3fffff];
/* HQ = Q/2 */
pub const HQ: Elem =  [0xffffffffffffff80,0xffffffffffffffff,0xffffffffffffffff,0x1fffff];
/* QQ = Q/4 */
pub const QQ: Elem =  [0xffffffffffffffc0,0xffffffffffffffff,0xffffffffffffffff,0xfffff];
/* TQQ = 3Q/4 */
pub const TQQ: Elem = [0xffffffffffffff40,0xffffffffffffffff,0xffffffffffffffff,0x2fffff];

/* -Q^{-1} mod 2^64, the Montgomery reduction constant */
const QINV: u64 = 0xfefefefefefefeff;
/* R^2 mod Q for R = 2^256 */
const R2: Elem = [0x0, 0x0000000fe0100000, 0x0, 0x0];

/* c = a + b - Q*carry, i.e. a + b reduced into [0,Q) */
pub fn add(c: &mut Elem, a: Elem, b: Elem) {
    let mut t: Elem = [0; NLIMBS];
    let mut carry = 0u64;

    for i in 0..NLIMBS {
        let s = (a[i] as u128) + (b[i] as u128) + (carry as u128);
        t[i] = s as u64;
        carry = (s >> RAD) as u64;
    }

    csub(c, t, carry);
}

/* c = a - b, i.e. a - b + Q if a < b */
pub fn sub(c: &mut Elem, a: Elem, b: Elem) {
    let mut borrow = 0u64;

    for i in 0..NLIMBS {
        let d = (a[i] as u128).wrapping_sub(b[i] as u128).wrapping_sub(borrow as u128);
        c[i] = d as u64;
        borrow = ((d >> RAD) as u64) & 1;
    }

    /* mask = all ones iff we borrowed out of the top limb */
    let mask = 0u64.wrapping_sub(borrow);
    let mut carry = 0u64;
    for i in 0..NLIMBS {
        let s = (c[i] as u128) + ((Q[i] & mask) as u128) + (carry as u128);
        c[i] = s as u64;
        carry = (s >> RAD) as u64;
    }
}

/* c = a*b*R^{-1} mod Q, Montgomery multiplication (CIOS).
 *
 * The accumulator never needs its top limb: entering each iteration t < 2Q, and
 * 2Q < 2^215, so the conditional subtraction at the end suffices to reduce. */
pub fn mul(c: &mut Elem, a: Elem, b: Elem) {
    #[cfg(all(target_arch = "aarch64", not(feature = "portable-fq")))]
    unsafe { mul_aarch64(c, a, b) }
    #[cfg(not(all(target_arch = "aarch64", not(feature = "portable-fq"))))]
    mul_portable(c, a, b)
}

/* c = a*b*R^{-1} mod Q, hand-written for aarch64.
 *
 * Two changes against `mul_portable`.  The carries run through the condition
 * flags as `adds`/`adcs`/`sbcs` chains, where the portable version has to
 * materialise each carry into a register, and the reduction uses the shape of
 * Q: since Q = 2^214 - 255, the product m*Q is (m << 214) - 255*m, a pair of
 * shifts and a small multiple, in place of four limb multiplications.  The two
 * agree on every input, which `test_mul_agrees_with_portable` checks.
 *
 * Bounds.  a, b < Q < 2^214, so the top limbs a[3], Q[3] are below 2^22 and the
 * high word of any product with them is below 2^22; hence t[4] stays small and
 * never carries out.  After the four rounds t < 2Q < 2^215, so one conditional
 * subtraction reduces, done here as t + 255 - 2^214 to avoid loading Q. */
#[cfg(all(target_arch = "aarch64", not(feature = "portable-fq")))]
unsafe fn mul_aarch64(c: &mut Elem, a: Elem, b: Elem) {
    /* a0-a3 hold a, t0-t4 the accumulator, l0-l3 and h0-h3 the low and the
     * high words of the four products of one round, bi is b[i] and m the
     * Montgomery multiplier; which machine register each becomes is left to
     * the compiler. */
    macro_rules! round { ($off:literal) => { concat!("
        ldr   {bi}, [{bp}, #", $off, "]
        mul   {l0}, {a0}, {bi}
        umulh {h0}, {a0}, {bi}
        mul   {l1}, {a1}, {bi}
        umulh {h1}, {a1}, {bi}
        mul   {l2}, {a2}, {bi}
        umulh {h2}, {a2}, {bi}
        mul   {l3}, {a3}, {bi}
        umulh {h3}, {a3}, {bi}
        adds  {l1}, {l1}, {h0}
        adcs  {l2}, {l2}, {h1}
        adcs  {l3}, {l3}, {h2}
        adc   {h3}, {h3}, xzr
        adds  {t0}, {t0}, {l0}
        adcs  {t1}, {t1}, {l1}
        adcs  {t2}, {t2}, {l2}
        adcs  {t3}, {t3}, {l3}
        adc   {t4}, {t4}, {h3}
        mul   {m},  {t0}, {qinv}
        lsl   {l0}, {m},  #22
        lsr   {l1}, {m},  #42
        adds  {t3}, {t3}, {l0}
        adc   {t4}, {t4}, {l1}
        lsl   {l2}, {m},  #8
        lsr   {l3}, {m},  #56
        subs  {l2}, {l2}, {m}
        sbc   {l3}, {l3}, xzr
        subs  {t0}, {t0}, {l2}
        sbcs  {t1}, {t1}, {l3}
        sbcs  {t2}, {t2}, xzr
        sbcs  {t3}, {t3}, xzr
        sbc   {t4}, {t4}, xzr
        mov   {t0}, {t1}
        mov   {t1}, {t2}
        mov   {t2}, {t3}
        mov   {t3}, {t4}
        mov   {t4}, xzr
    ") } }

    core::arch::asm!(
        "ldp  {a0}, {a1}, [{ap}]",
        "ldp  {a2}, {a3}, [{ap}, #16]",
        "mov  {t0}, xzr",
        "mov  {t1}, xzr",
        "mov  {t2}, xzr",
        "mov  {t3}, xzr",
        "mov  {t4}, xzr",
        round!("0"), round!("8"), round!("16"), round!("24"),
        /* t >= Q ?  form t + 255 - 2^214 and keep it if it did not borrow */
        "adds {l0}, {t0}, #255",
        "adcs {l1}, {t1}, xzr",
        "adcs {l2}, {t2}, xzr",
        "adcs {l3}, {t3}, xzr",
        "subs {l3}, {l3}, #0x400000",
        "csel {t0}, {l0}, {t0}, cs",
        "csel {t1}, {l1}, {t1}, cs",
        "csel {t2}, {l2}, {t2}, cs",
        "csel {t3}, {l3}, {t3}, cs",
        "stp  {t0}, {t1}, [{cp}]",
        "stp  {t2}, {t3}, [{cp}, #16]",
        ap = in(reg) a.as_ptr(),
        bp = in(reg) b.as_ptr(),
        cp = in(reg) c.as_mut_ptr(),
        qinv = in(reg) QINV,
        a0 = out(reg) _, a1 = out(reg) _, a2 = out(reg) _, a3 = out(reg) _,
        t0 = out(reg) _, t1 = out(reg) _, t2 = out(reg) _, t3 = out(reg) _,
        t4 = out(reg) _,
        l0 = out(reg) _, l1 = out(reg) _, l2 = out(reg) _, l3 = out(reg) _,
        h0 = out(reg) _, h1 = out(reg) _, h2 = out(reg) _, h3 = out(reg) _,
        bi = out(reg) _, m = out(reg) _,
        options(nostack),
    );
}

/* c = a*b*R^{-1} mod Q, Montgomery multiplication (CIOS), portable. */
pub fn mul_portable(c: &mut Elem, a: Elem, b: Elem) {
    let mut t: [u64; NLIMBS+2] = [0; NLIMBS+2];

    for i in 0..NLIMBS {
        /* t += a * b[i] */
        let mut carry: u128 = 0;
        for j in 0..NLIMBS {
            let s = (t[j] as u128) + (a[j] as u128)*(b[i] as u128) + carry;
            t[j] = s as u64;
            carry = s >> RAD;
        }
        let s = (t[NLIMBS] as u128) + carry;
        t[NLIMBS] = s as u64;
        t[NLIMBS+1] = (s >> RAD) as u64;

        /* t += Q * m, with m chosen so that the bottom limb of t becomes zero */
        let m = t[0].wrapping_mul(QINV);
        carry = 0;
        for j in 0..NLIMBS {
            let s = (t[j] as u128) + (Q[j] as u128)*(m as u128) + carry;
            t[j] = s as u64;
            carry = s >> RAD;
        }
        debug_assert_eq!(t[0], 0);
        let s = (t[NLIMBS] as u128) + carry;
        t[NLIMBS] = s as u64;
        t[NLIMBS+1] += (s >> RAD) as u64;

        /* t >>= 64 */
        for j in 0..=NLIMBS {
            t[j] = t[j+1];
        }
        t[NLIMBS+1] = 0;
    }

    csub(c, [t[0], t[1], t[2], t[3]], t[NLIMBS]);
}

pub fn toM(aM: &mut Elem, a: Elem) {
    mul(aM, a, R2);
}

pub fn fromM(a: &mut Elem, aM: Elem) {
    mul(a, aM, [1, 0, 0, 0]);
}

/* c = t - Q if t >= Q (as a 5-limb value t + high*2^256), else c = t */
fn csub(c: &mut Elem, t: Elem, high: u64) {
    let mut d: Elem = [0; NLIMBS];
    let mut borrow = 0u64;

    for i in 0..NLIMBS {
        let s = (t[i] as u128).wrapping_sub(Q[i] as u128).wrapping_sub(borrow as u128);
        d[i] = s as u64;
        borrow = ((s >> RAD) as u64) & 1;
    }

    /* keep the difference unless it went negative and nothing came in from above */
    let take = high | (1 - borrow);
    let mask = 0u64.wrapping_sub(take & 1);
    for i in 0..NLIMBS {
        c[i] = (d[i] & mask) | (t[i] & !mask);
    }
}

pub fn fp_init() -> Elem {
    [0; NLIMBS]
}
/* Description: Constant time comparison of two field elements
 *
 * Returns -1 if a < b; 0 if a = b; 1 if a > b
 */
pub fn cmp(a: Elem, b: Elem) -> u8 {
    let mut r: u8 = 0;
    let mut mask: u8 = 0xff;
    let mut s_ai: i64;
    let mut s_bi: i64;
    let mut gt: u8;
    let mut lt: u8;

    for i in (0..NLIMBS).rev() {
        s_ai = a[i] as i64;
        s_bi = b[i] as i64;

        lt = (((s_ai - s_bi) as u64) >> 63) as u8;
        gt = (((s_bi - s_ai) as u64) >> 63) as u8;

        //  high order limb comparisons take precedence
        lt &= mask;
        gt &= mask;

        mask ^= (lt | gt);
        r |= (lt << 7);
        r |= gt;
    }

    r
}

/*
 * Converts stream of bytes into value of type Elem
 */
pub fn elem_frombytes(ep: &[u8; ELEM_BYTES]) -> Elem {
    let mut e: Elem = fp_init();
    let mut t: [u8; 8] = [0; 8];

    for i in 0..NLIMBS-1 {
        e[i] = u64::from_le_bytes(ep[8*i..8*i+8].try_into().unwrap());
    }

    t[0..ELEM_BYTES-8*(NLIMBS-1)].copy_from_slice(&ep[8*(NLIMBS-1)..8*(NLIMBS-1)+3]);
    e[NLIMBS-1] = u64::from_le_bytes(t);

    e
}
/*
 * Converts field element into a byte buffer
 */
pub fn elem_tobytes(e: Elem) -> [u8; ELEM_BYTES] {
    let mut r: [u8; ELEM_BYTES] = [0; ELEM_BYTES];

    for i in 0..NLIMBS-1 {
        r[8*i..8*i+8].copy_from_slice(&e[i].to_le_bytes());
    }

    r[8*(NLIMBS-1)..8*(NLIMBS-1)+3].copy_from_slice(&e[NLIMBS-1].to_le_bytes()[0..3]); //remove trailing bytes

    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fp_add() {
        let a: Elem = QQ.clone();
        let b: Elem = QQ.clone();
        let rc: Elem = HQ.clone();
        let mut c: Elem = fp_init();

        add(&mut c, a, b); // q/4 + q/4

        assert_eq!(rc, c, "fp_add: Values don't match");
    }

    #[test]
    fn test_fp_sub() {
        let a: Elem = HQ.clone();
        let b: Elem = QQ.clone();
        let rc: Elem = QQ.clone();
        let mut c: Elem = fp_init();

        sub(&mut c, a, b); // q/2 - q/4

        assert_eq!(rc, c, "fp_sub: Values don't match");
    }

    /* `mul` is hand-written assembly on aarch64; check it against the portable
     * version on the boundary values and on inputs that span the range. */
    #[test]
    fn test_mul_agrees_with_portable() {
        let mut probes: Vec<Elem> = vec![
            [0, 0, 0, 0], [1, 0, 0, 0], QQ, HQ, TQQ, Q,
            [Q[0]-1, Q[1], Q[2], Q[3]],
        ];
        let mut x: u64 = 0x243f6a8885a308d3;
        for _ in 0..64 {
            let mut e: Elem = fp_init();
            for j in 0..NLIMBS {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                e[j] = x;
            }
            e[NLIMBS-1] &= 0x3fffff;
            probes.push(e);
        }

        for &a in probes.iter() {
            for &b in probes.iter() {
                let mut want: Elem = fp_init();
                let mut got: Elem = fp_init();
                mul_portable(&mut want, a, b);
                mul(&mut got, a, b);
                assert_eq!(want, got, "mul disagrees on {:x?} * {:x?}", a, b);
            }
        }
    }

    #[test]
    fn test_fp_mul() {
        let a: Elem = QQ.clone();
        let mut aM: Elem = fp_init();
        let b: Elem = [0x03, 0x0, 0x0, 0x0];
        let mut bM: Elem = fp_init();
        let rc: Elem = TQQ.clone();
        let mut c: Elem = fp_init();
        let mut cM: Elem = fp_init();

        toM(&mut aM, a);
        toM(&mut bM, b);
        mul(&mut cM, aM, bM); // q/4 * 3
        fromM(&mut c, cM);

        assert_eq!(rc, c, "fp_mul: Values don't match");
    }

    #[test]
    fn test_fp_bytes() {
        let a: Elem = QQ.clone();
        let mut r: Elem;
        let mut b: [u8; ELEM_BYTES] = [0; ELEM_BYTES];

        b = elem_tobytes(a);
        r = elem_frombytes(&b);
        assert_eq!(a, r, "fp_bytes: Values don't match");
    }

    #[test]
    fn test_cmp() {
        let a: Elem = HQ.clone();

        assert_eq!(cmp(a, HQ), 0x00);
        assert_eq!(cmp(a, QQ), 0x01);
        assert_eq!(cmp(a, TQQ), 0x80);
    }
}
