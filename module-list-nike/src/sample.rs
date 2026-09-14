//! The three distributions of the scheme.
//!
//! * uniform on `R_q`, for `A` and for `s`;
//! * `CbdVec(eta)` with `eta = 2`, for the MLWE error `e`;
//! * uniform on `[-beta, beta]^d` with `beta = 16`, for the short secret `u`.
//!
//! All three read from an AES-128-CTR keystream.  The uniform and the `beta`
//! sampler reject, so the number of keystream bytes consumed varies; the
//! rejection rate is `21/2^23` for the former and `31/64` for the latter.

use crate::aes::Ctr;
use crate::params::{BETA, D, ETA, Q};
use crate::poly::{Poly, PolyExp, PolyI32};

/// Uniform on `R_q`, coefficients in `[0, q)`.
pub fn uniform(ctr: &mut Ctr) -> Poly {
    let mut p = [0u32; D];
    let mut k = 0;
    while k < D {
        let v = ctr.u24() & ((1 << 23) - 1);
        if v < Q {
            p[k] = v;
            k += 1;
        }
    }
    p
}

/// Uniform on `R_q`, produced directly in negacyclic expanded form.
pub fn uniform_exp(ctr: &mut Ctr) -> PolyExp {
    let mut e = [0i32; 2 * D];
    let mut k = 0;
    while k < D {
        let v = ctr.u24() & ((1 << 23) - 1);
        if v < Q {
            e[k] = -(v as i32);
            e[D + k] = v as i32;
            k += 1;
        }
    }
    e
}

/// The centred binomial distribution of width `eta = 2`, coefficientwise.
///
/// Each coefficient consumes `2 eta = 4` keystream bits: `a` counts the low
/// `eta`, `b` the high `eta`, and the coefficient is `a - b` in `[-2, 2]`.
pub fn cbd(ctr: &mut Ctr) -> PolyI32 {
    const _: () = assert!(ETA == 2);
    let mut p = [0i32; D];
    for k in (0..D).step_by(2) {
        let byte = ctr.byte();
        for (h, slot) in [(0u8, k), (4u8, k + 1)] {
            let nib = (byte >> h) & 0xf;
            let a = (nib & 1) + ((nib >> 1) & 1);
            let b = ((nib >> 2) & 1) + ((nib >> 3) & 1);
            p[slot] = a as i32 - b as i32;
        }
    }
    p
}

/// Uniform on `[-beta, beta]^d`, the distribution of one component of `u`.
pub fn short(ctr: &mut Ctr) -> PolyI32 {
    let support = 2 * BETA + 1; // 33 values
    let mut p = [0i32; D];
    let mut k = 0;
    let mut bits: u32 = 0;
    let mut have = 0u32;
    while k < D {
        if have < 6 {
            bits |= (ctr.byte() as u32) << have;
            have += 8;
        }
        let v = (bits & 0x3f) as i32;
        bits >>= 6;
        have -= 6;
        if v < support {
            p[k] = v - BETA;
            k += 1;
        }
    }
    p
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ctr() -> Ctr {
        Ctr::new(&[0x11u8; 16], 0)
    }

    #[test]
    fn uniform_in_range() {
        let mut c = ctr();
        for _ in 0..200 {
            for v in uniform(&mut c) {
                assert!(v < Q);
            }
        }
    }

    #[test]
    fn uniform_exp_matches_uniform() {
        let mut a = ctr();
        let mut b = ctr();
        for _ in 0..50 {
            let p = uniform(&mut a);
            let e = uniform_exp(&mut b);
            assert_eq!(crate::poly::expand_u32(&p), e);
        }
    }

    #[test]
    fn cbd_range_and_mean() {
        let mut c = ctr();
        let mut hist = [0u64; 5];
        for _ in 0..2000 {
            for v in cbd(&mut c) {
                assert!((-2..=2).contains(&v));
                hist[(v + 2) as usize] += 1;
            }
        }
        // Binomial(4, 1/2) shifted: weights 1, 4, 6, 4, 1 out of 16.
        let n: u64 = hist.iter().sum();
        for (i, w) in [1.0, 4.0, 6.0, 4.0, 1.0].iter().enumerate() {
            let want = w / 16.0;
            let got = hist[i] as f64 / n as f64;
            assert!((got - want).abs() < 0.01, "bucket {i}: {got} vs {want}");
        }
    }

    #[test]
    fn short_range_and_uniformity() {
        let mut c = ctr();
        let mut hist = [0u64; 33];
        for _ in 0..2000 {
            for v in short(&mut c) {
                assert!((-BETA..=BETA).contains(&v));
                hist[(v + BETA) as usize] += 1;
            }
        }
        let n: u64 = hist.iter().sum();
        for (i, h) in hist.iter().enumerate() {
            let got = *h as f64 / n as f64;
            assert!((got - 1.0 / 33.0).abs() < 0.005, "bucket {i}: {got}");
        }
    }
}
