//! Module list-NIKE over `R_q = Z_q[X]/(X^d+1)`, at the recommended parameter
//! set of `tab:mod-params`:
//!
//! ```text
//!   d = 128   rho = 2   c = 128   beta = 16   n = 9   m = 100
//!   q = 2^23 - 21       |pk| = 40112 B      |sk| = 12912 B
//!   MLWE hardness 133 bits          correctness 2^-103
//! ```
//!
//! The scheme is [`nike`]; [`listrec`] is the list reconciliation of
//! `fig:adv-list-rec`; [`fq`] and [`poly`] are the arithmetic; [`aes`] is the
//! expansion function and [`sample`] the three distributions.

pub mod aes;
pub mod cycles;
pub mod fq;
pub mod listrec;
pub mod nike;
pub mod params;
pub mod poly;
pub mod sample;
