//! The recommended parameter set of the module list-NIKE (`tab:mod-params`, shaded row).
//!
//! ```text
//!   d = 128   rho = 2   c = 128   beta = 16   n = 9   m = 100
//!   q = 2^23 - 21       |pk| = 40112 B        |sk| = 12912 B
//!   MLWE bit security 133,  correctness 2^-103 at ell = 2^24
//! ```

/// Degree of the ring `R_q = Z_q[X]/(X^d + 1)`.
pub const D: usize = 128;

/// Modulus, the largest prime below `2^23` with multiplicative order `d/2` modulo `2d`.
pub const Q: u32 = (1 << 23) - 21;

/// Bit length of a `Z_q` coefficient, so `ceil(log2 q)`.
pub const LOG_Q: usize = 23;

/// Module rank: `s <- R_q^n`, the height of `A`.
pub const N: usize = 9;

/// MLWE samples: `e in R_q^m`, the width of `A`.
pub const M: usize = 100;

/// Binomial width of the error distribution `Cbd(eta)`.
pub const ETA: usize = 2;

/// Secret norm: `u <- ([-beta, beta]^d)^m`.
pub const BETA: i32 = 16;

/// Bits extracted per coefficient of the raw key.
pub const RHO: usize = 2;

/// Coefficients of the raw key that are reconciled, of the `d` available.
pub const C: usize = 128;

/// `log2` of the list length. The list has `ELL = 2^LOG_ELL` entries.
pub const LOG_ELL: usize = 24;

/// List length enumerated by `ListRec`.
pub const ELL: u64 = 1 << LOG_ELL;

/// How many hedged coefficients the receiver deviates on before it gives up on
/// the light entries and walks the whole list.  Only the coefficients whose raw
/// values straddle a rounding boundary have to deviate, and over 1000 honest
/// pairs that was at most two of them, so 3 covers every case we have observed
/// at a cost of `1 + L + L(L-1)/2 + L(L-1)(L-2)/6 = 2325` entries.
pub const WEIGHT: usize = 3;

/// Bits per coefficient of `u` when a secret key is serialised, enough for the
/// `2 beta + 1 = 33` values it takes.
pub const BETA_BITS: usize = 6;
const _: () = assert!((1 << BETA_BITS) >= 2 * BETA + 1);

/// Shared key length in bytes: `rho` bits out of each of the `c` coefficients.
pub const SHK_BYTES: usize = RHO * C / 8;

/// `ListRec` hedges `log ell` of the `c` reconciled coefficients.
const _: () = assert!(LOG_ELL <= C);

/// The sizes `tab:mod-params` prints for this row.
///
/// Nothing in the implementation is built from these.  They are the paper's
/// claim; `bin/check` serialises a real key pair and a real shared key and
/// compares the lengths it measures against them.
pub mod table {
    pub const PK_BYTES: usize = 40112;
    pub const SK_BYTES: usize = 12912;
    pub const SHK_BITS: usize = 256;
    pub const MLWE_BITS: u32 = 133;
    pub const CORRECTNESS_LOG2: i32 = -103;
}
