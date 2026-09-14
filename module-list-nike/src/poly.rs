//! Ring arithmetic in `R_q = Z_q[X]/(X^d + 1)` with `d = 128`.
//!
//! # Why there is no NTT
//!
//! The paper fixes `f = ord_{2d}(q) = d/2`, so `X^d + 1` splits into `g = 2`
//! irreducible factors of degree `64` and at most one NTT layer could ever
//! apply.  That layer would be `X^128 + 1 = (X^64 - r)(X^64 + r)` with
//! `r^2 = -1`, which needs `q == 1 (mod 4)`.  For the recommended modulus
//! `q = 2^23 - 21` we have `q == 3 (mod 4)`, so `-1` is a non-residue and the
//! layer does not exist: the two degree-64 factors are not of the form
//! `X^64 -+ r` and a CRT map onto them costs as much as the multiplication it
//! would replace.  Products are therefore schoolbook, which is cheap here
//! because `q < 2^23` lets an `i64` carry a whole inner product unreduced.
//!
//! # Representation
//!
//! A polynomial is `[u32; D]` in `[0, q)`.  For multiplication one operand is
//! held in the *negacyclic expanded* form [`PolyExp`],
//!
//! ```text
//!   exp[t] = -b[t]      for t in 0..D
//!   exp[D + t] = b[t]   for t in 0..D
//! ```
//!
//! so that the negacyclic product is a plain sliding dot product,
//! `c[k] = sum_i a[i] * exp[D + k - i]`, with the wrap-around sign already
//! folded into the data.  The upper half of `exp` is the ordinary coefficient
//! vector, so a `PolyExp` also serves wherever plain coefficients are wanted.

use crate::fq;
use crate::params::D;

/// A ring element with coefficients in `[0, q)`.
pub type Poly = [u32; D];

/// Plain signed coefficients, the shape consumed by [`mac`] as its left operand.
pub type PolyI32 = [i32; D];

/// Negacyclic expanded form of a ring element; see the module documentation.
pub type PolyExp = [i32; 2 * D];

/// An unreduced accumulator holding one ring element.
pub type Acc = [i64; D];

/// A fresh zero accumulator.
#[inline]
pub fn acc_zero() -> Acc {
    [0i64; D]
}

/// Builds the negacyclic expanded form of `b`, whose coefficients must be centred.
pub fn expand(b: &PolyI32) -> PolyExp {
    let mut e = [0i32; 2 * D];
    for t in 0..D {
        e[t] = -b[t];
        e[D + t] = b[t];
    }
    e
}

/// Builds the negacyclic expanded form of `b` given coefficients in `[0, q)`.
pub fn expand_u32(b: &Poly) -> PolyExp {
    let mut e = [0i32; 2 * D];
    for t in 0..D {
        let v = b[t] as i32;
        e[t] = -v;
        e[D + t] = v;
    }
    e
}

/// Reinterprets coefficients in `[0, q)` as the signed shape [`mac`] wants.
/// Both operands of the inner multiply must be *signed* 32-bit for the compiler
/// to contract it into a widening multiply-accumulate, so nothing in the hot
/// path is ever `u32`.
#[inline]
pub fn to_i32(a: &Poly) -> PolyI32 {
    core::array::from_fn(|k| a[k] as i32)
}

/// The plain coefficient vector inside an expanded polynomial.
#[inline(always)]
pub fn plain(e: &PolyExp) -> &[i32; D] {
    (&e[D..]).try_into().unwrap()
}

/// `acc += a * b` in `R_q`, unreduced.
///
/// `a` holds plain coefficients, `b` is in expanded form, so the product is the
/// sliding dot product `sum_i a[i] * b[D + k - i]`.  Eight rows of `a` are fused
/// into one pass over the accumulator: `acc` is 1 KiB and does not stay in
/// registers, and touching it once per eight rows instead of once per row is
/// worth about 1.5x.  There is no data-dependent branch, because `a` is the
/// secret `s` at one of the call sites and skipping zero coefficients would leak
/// it.
///
/// The caller owns the accumulator bound: one call adds at most
/// `D * max|a| * max|b|`, and [`reduce_acc`] needs the total below `2^58`.  The
/// two uses in this crate are `D q^2 < 2^53.1` repeated `n = 9` times, and
/// `D q beta < 2^34` repeated `m = 100` times.
#[inline]
pub fn mac(acc: &mut Acc, a: &[i32; D], b: &PolyExp) {
    /// Rows of `a` fused into one pass over `acc`.
    const B: usize = 8;
    const _: () = assert!(D % B == 0);

    for i in (0..D).step_by(B) {
        // Row r is a[i+r] * X^(i+r); in expanded form its window of `b` starts
        // at D - i - r.
        let coef: [i64; B] = core::array::from_fn(|r| a[i + r] as i64);
        let src: [&[i32]; B] = core::array::from_fn(|r| &b[D - i - r..2 * D - i - r]);
        for k in 0..D {
            let mut t = acc[k];
            for r in 0..B {
                t += coef[r] * src[r][k] as i64;
            }
            acc[k] = t;
        }
    }
}

/// Reduces an accumulator to a ring element with coefficients in `[0, q)`.
#[inline]
pub fn reduce_acc(acc: &Acc) -> Poly {
    let mut r = [0u32; D];
    for k in 0..D {
        r[k] = fq::reduce(acc[k]);
    }
    r
}

/// Adds a centred polynomial into an accumulator.
#[inline]
pub fn acc_add_centred(acc: &mut Acc, a: &PolyI32) {
    for k in 0..D {
        acc[k] += a[k] as i64;
    }
}

/// Coefficientwise sum in `R_q`.
pub fn add(a: &Poly, b: &Poly) -> Poly {
    let mut c = [0u32; D];
    for k in 0..D {
        c[k] = fq::add(a[k], b[k]);
    }
    c
}

/// Schoolbook negacyclic product, reduced.  Reference path for the tests.
pub fn mul_ref(a: &Poly, b: &Poly) -> Poly {
    let mut acc = acc_zero();
    for i in 0..D {
        for j in 0..D {
            let p = a[i] as i64 * b[j] as i64;
            if i + j < D {
                acc[i + j] += p;
            } else {
                acc[i + j - D] -= p;
            }
        }
    }
    // The double loop can reach D * q^2 < 2^53.1, inside the reduce bound.
    reduce_acc(&acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Q;

    fn rnd_poly(st: &mut u64) -> Poly {
        let mut p = [0u32; D];
        for k in 0..D {
            *st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            p[k] = ((*st >> 24) % Q as u64) as u32;
        }
        p
    }

    #[test]
    fn mac_matches_schoolbook() {
        let mut st = 0xdeadbeefu64;
        for _ in 0..40 {
            let a = rnd_poly(&mut st);
            let b = rnd_poly(&mut st);
            let want = mul_ref(&a, &b);
            let mut acc = acc_zero();
            mac(&mut acc, &to_i32(&a), &expand_u32(&b));
            assert_eq!(reduce_acc(&acc), want);
        }
    }

    #[test]
    fn mac_is_commutative() {
        let mut st = 7u64;
        for _ in 0..40 {
            let a = rnd_poly(&mut st);
            let b = rnd_poly(&mut st);
            let mut x = acc_zero();
            let mut y = acc_zero();
            mac(&mut x, &to_i32(&a), &expand_u32(&b));
            mac(&mut y, &to_i32(&b), &expand_u32(&a));
            assert_eq!(reduce_acc(&x), reduce_acc(&y));
        }
    }

    /// `X^d = -1`, the defining relation of the ring.
    #[test]
    fn negacyclic_relation() {
        let mut a = [0u32; D];
        a[1] = 1; // X
        let mut b = [0u32; D];
        b[D - 1] = 1; // X^(d-1)
        let c = mul_ref(&a, &b); // X^d = -1
        let mut want = [0u32; D];
        want[0] = Q - 1;
        assert_eq!(c, want);
    }

    #[test]
    fn nine_term_inner_product_stays_in_range() {
        let mut st = 99u64;
        let mut acc = acc_zero();
        let mut want = [0u32; D];
        for _ in 0..crate::params::N {
            let a = rnd_poly(&mut st);
            let b = rnd_poly(&mut st);
            mac(&mut acc, &to_i32(&a), &expand_u32(&b));
            want = add(&want, &mul_ref(&a, &b));
        }
        assert_eq!(reduce_acc(&acc), want);
    }
}
