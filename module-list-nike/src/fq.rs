//! Arithmetic in `Z_q` for `q = 2^23 - 21`.
//!
//! The paper picks `q` as the largest prime below a power of two, so `q` is
//! pseudo-Mersenne and `2^23 = q + 21 == 21 (mod q)`.  A wide accumulator is
//! folded down by three rounds of shift, mask and multiply-add followed by one
//! conditional subtraction.  No Montgomery domain and no 128-bit multiply.
//!
//! Coefficients are held in `[0, q)` as `u32`.  Products of two of them are
//! below `2^46`, so `i64` accumulates more than `2^17` of them without
//! overflow; every inner product in the scheme stays well inside that.

use crate::params::{LOG_Q, Q};

/// `2^LOG_Q - q`, the fold-in constant.
const DELTA: i64 = (1i64 << LOG_Q) - Q as i64;

/// Low `LOG_Q` bits.
const MASK: i64 = (1i64 << LOG_Q) - 1;

/// A multiple of `q`, added to a signed accumulator to make it non-negative
/// without changing its class mod `q`.  It exceeds `2^57`, so it covers every
/// accumulator in this crate (the largest is below `2^56.2`), and stays below
/// `2^58`, so the sum still fits the input range of `reduce_u`.
pub const BIAS: i64 = (Q as i64) << 35;

const _: () = assert!(BIAS % Q as i64 == 0);
const _: () = assert!(BIAS > (1i64 << 57));
const _: () = assert!(BIAS < (1i64 << 58));

/// Reduces a non-negative accumulator below `2^59` to `[0, q)`.
#[inline(always)]
pub fn reduce_u(x: i64) -> u32 {
    debug_assert!(x >= 0 && x < (1i64 << 59));
    let x = (x >> LOG_Q) * DELTA + (x & MASK); //  < 2^41
    let x = (x >> LOG_Q) * DELTA + (x & MASK); //  < 2^24
    let x = (x >> LOG_Q) * DELTA + (x & MASK); //  < q + 42
    let x = x - Q as i64;
    (x + (((x >> 63) & 1) * Q as i64)) as u32
}

/// Reduces a signed accumulator with `|x| < 2^58` to `[0, q)`.
#[inline(always)]
pub fn reduce(x: i64) -> u32 {
    debug_assert!(x > -(1i64 << 58) && x < (1i64 << 58));
    reduce_u(x + BIAS)
}

/// `(a + b) mod q` for `a, b` in `[0, q)`.
#[inline(always)]
pub fn add(a: u32, b: u32) -> u32 {
    let t = a.wrapping_add(b).wrapping_sub(Q);
    t.wrapping_add((t >> 31).wrapping_neg() & Q)
}

/// `(a - b) mod q` for `a, b` in `[0, q)`.
#[inline(always)]
pub fn sub(a: u32, b: u32) -> u32 {
    let t = a.wrapping_sub(b);
    t.wrapping_add((t >> 31).wrapping_neg() & Q)
}

/// `(a * b) mod q` for `a, b` in `[0, q)`.
#[inline(always)]
pub fn mul(a: u32, b: u32) -> u32 {
    reduce_u(a as i64 * b as i64)
}

/// Lifts a centred representative in `(-q/2, q/2]` into `[0, q)`.
#[inline(always)]
pub fn from_centred(a: i32) -> u32 {
    (a + (((a >> 31) & 1) * Q as i32)) as u32
}

/// The centred representative of `a` in `(-q/2, q/2]`.
#[inline(always)]
pub fn to_centred(a: u32) -> i32 {
    let a = a as i32;
    a - (((Q as i32 / 2 - a) >> 31) & 1) * Q as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Textbook reduction, used as the oracle.
    fn slow(x: i64) -> u32 {
        x.rem_euclid(Q as i64) as u32
    }

    #[test]
    fn reduce_matches_rem_euclid() {
        let mut st = 0x243f_6a88_85a3_08d3u64;
        for _ in 0..200_000 {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let x = ((st >> 6) as i64) % (1i64 << 57);
            let x = if st & 1 == 0 { x } else { -x };
            assert_eq!(reduce(x), slow(x), "x = {x}");
        }
        for x in [0i64, 1, -1, Q as i64, -(Q as i64), (1i64 << 57) - 1, -((1i64 << 57) - 1)] {
            assert_eq!(reduce(x), slow(x), "x = {x}");
        }
    }

    #[test]
    fn ring_laws() {
        let mut st = 12345u64;
        let mut next = || {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((st >> 20) % Q as u64) as u32
        };
        for _ in 0..100_000 {
            let (a, b) = (next(), next());
            assert!(add(a, b) < Q && sub(a, b) < Q && mul(a, b) < Q);
            assert_eq!(add(a, b), ((a as u64 + b as u64) % Q as u64) as u32);
            assert_eq!(sub(a, b), ((a as u64 + Q as u64 - b as u64) % Q as u64) as u32);
            assert_eq!(mul(a, b), ((a as u64 * b as u64) % Q as u64) as u32);
            assert_eq!(from_centred(to_centred(a)), a);
        }
    }
}
