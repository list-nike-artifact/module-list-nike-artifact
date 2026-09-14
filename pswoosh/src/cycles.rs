//! Cycle counting on a machine without `rdtsc`.
//!
//! Apple silicon exposes no unprivileged cycle counter: `PMCR0`/`PMC0` need a
//! kernel entry point, and `CNTVCT_EL0` is a fixed 24 MHz timebase, not a cycle
//! count.  So the numbers here are wall-clock nanoseconds multiplied by a
//! frequency that the program measures for itself.
//!
//! [`calibrate_ghz`] times a chain of dependent `add` instructions.  Every one
//! of them depends on the previous result, so the core retires exactly one per
//! cycle no matter how wide it is; the elapsed time therefore divides out to the
//! clock period directly.  Before measuring anything, [`pin_to_p_core`] raises
//! the thread's QoS class so that the scheduler keeps it on a performance core
//! at its top frequency.

use std::time::Instant;

/// `qos_class_t::QOS_CLASS_USER_INTERACTIVE`, from `<sys/qos.h>`.
const QOS_CLASS_USER_INTERACTIVE: u32 = 0x21;

extern "C" {
    fn pthread_set_qos_class_self_np(qos_class: u32, relative_priority: i32) -> i32;
}

/// Asks the scheduler for a performance core at full frequency.
pub fn pin_to_p_core() {
    // Safe: the call only sets a scheduling attribute of the calling thread.
    unsafe {
        pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
    }
}

/// Eight dependent increments; each one needs the previous result, so the core
/// retires exactly one per cycle however wide it is.
macro_rules! add8 {
    () => {
        concat!(
            "add {x}, {x}, #1\n", "add {x}, {x}, #1\n", "add {x}, {x}, #1\n",
            "add {x}, {x}, #1\n", "add {x}, {x}, #1\n", "add {x}, {x}, #1\n",
            "add {x}, {x}, #1\n", "add {x}, {x}, #1\n",
        )
    };
}

/// Defines `fn $name(iters) -> nanoseconds`, a loop whose body is `$body`.
macro_rules! chain_fn {
    ($name:ident, $body:expr) => {
        #[cfg(target_arch = "aarch64")]
        #[inline(never)]
        fn $name(iters: u64) -> f64 {
            let mut x: u64 = 0;
            let t = Instant::now();
            // Safe: no memory is touched, and both operands are scratch.
            unsafe {
                core::arch::asm!(
                    "2:",
                    $body,
                    "subs {i}, {i}, #1",
                    "b.ne 2b",
                    x = inout(reg) x,
                    i = inout(reg) iters => _,
                    options(nostack),
                );
            }
            let ns = t.elapsed().as_nanos() as f64;
            std::hint::black_box(x);
            ns
        }
    };
}

chain_fn!(chain8, add8!());
chain_fn!(chain24, concat!(add8!(), add8!(), add8!()));

/// Portable stand-in: `black_box` keeps the chain dependent and unoptimised.
#[cfg(not(target_arch = "aarch64"))]
fn chain_generic(iters: u64, adds: usize) -> f64 {
    let mut x: u64 = 0;
    let t = Instant::now();
    for _ in 0..iters {
        for _ in 0..adds {
            x = std::hint::black_box(x + 1);
        }
    }
    let ns = t.elapsed().as_nanos() as f64;
    std::hint::black_box(x);
    ns
}
#[cfg(not(target_arch = "aarch64"))]
fn chain8(iters: u64) -> f64 {
    chain_generic(iters, 8)
}
#[cfg(not(target_arch = "aarch64"))]
fn chain24(iters: u64) -> f64 {
    chain_generic(iters, 24)
}

/// The core frequency in GHz, measured rather than assumed.
///
/// Two loops are timed: one with 8 increments in the body, one with 24.
/// Everything else about them is identical, and the loop counter and branch sit
/// off the dependency chain, so their cost cancels in the difference.  What is
/// left is 16 increments, hence 16 cycles, per iteration.  The fastest of
/// several attempts is kept, since interference can only slow a run down.
pub fn calibrate_ghz() -> f64 {
    pin_to_p_core();
    const ITERS: u64 = 1 << 21;
    const EXTRA_ADDS: f64 = 16.0;
    // The clock ramps up over the first few milliseconds of load, so spin once
    // before believing anything.
    for _ in 0..4 {
        std::hint::black_box(chain24(ITERS));
    }
    let mut best = 0.0f64;
    for _ in 0..24 {
        let d = chain24(ITERS) - chain8(ITERS);
        if d > 0.0 {
            let ghz = ITERS as f64 * EXTRA_ADDS / d;
            if ghz > best {
                best = ghz;
            }
        }
    }
    best
}

/// One measured operation: sample count, median and mean, in cycles.
pub struct Stat {
    pub label: String,
    pub runs: usize,
    pub median: f64,
    pub mean: f64,
}

impl Stat {
    /// Median cycles divided by a count, for "per ring product" figures.
    pub fn per(&self, n: usize) -> f64 {
        self.median / n as f64
    }

    /// Microseconds, given the clock the cycles were computed with.
    pub fn micros(&self, ghz: f64) -> f64 {
        self.median / ghz / 1000.0
    }
}

/// Times `f` `runs` times and converts to cycles at `ghz`.
///
/// The closure's result is passed through `black_box` so that nothing it
/// computes can be optimised away, and the whole body runs once before the
/// timed loop so that the first call's page faults and branch mispredictions do
/// not land in the sample.
pub fn measure<T, F: FnMut() -> T>(label: &str, ghz: f64, runs: usize, mut f: F) -> Stat {
    pin_to_p_core();
    std::hint::black_box(f());
    let mut samples = Vec::with_capacity(runs);
    for _ in 0..runs {
        let t = Instant::now();
        let out = f();
        let ns = t.elapsed().as_nanos() as f64;
        std::hint::black_box(out);
        samples.push(ns * ghz);
    }
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if runs % 2 == 1 {
        samples[runs / 2]
    } else {
        0.5 * (samples[runs / 2 - 1] + samples[runs / 2])
    };
    let mean = samples.iter().sum::<f64>() / runs as f64;
    Stat { label: label.to_string(), runs, median, mean }
}
