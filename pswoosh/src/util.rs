//! PORT: upstream's `util.rs` read the cycle counter with `core::arch::x86_64::_rdtsc`.
//! Apple silicon has no unprivileged cycle counter, so timing lives in
//! `crate::cycles` instead and this file keeps only the statistics helpers.

pub const NRUNS: usize = 1000;

pub fn median(t: &mut Vec<f64>) -> f64 {
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if t.len() % 2 == 1 {
        t[t.len()/2]
    } else {
        (t[t.len()/2 - 1] + t[t.len()/2]) / 2.0
    }
}
