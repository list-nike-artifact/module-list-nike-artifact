//! Passive-SWOOSH key generation and shared-key derivation, timed the same way
//! as `mlnike`: one pinned performance core, median of 1000 calls, wall clock
//! converted at a fixed clock rate.

use pswoosh::*;
use pswoosh::cycles::{calibrate_ghz, measure, pin_to_p_core};

const GHZ: f64 = 4.05;
const RUNS: usize = 1000;

fn main() {
    pswoosh::with_stack(run);
}

fn run() {
    pin_to_p_core();
    /* Let the clock ramp before anything is timed. */
    let warm = std::time::Instant::now();
    let mut acc: u64 = 0;
    while warm.elapsed().as_millis() < 300 {
        acc = acc.wrapping_mul(6364136223846793005).wrapping_add(1);
    }
    std::hint::black_box(acc);

    println!("calibrated clock: {:.2} GHz (using {:.2})", calibrate_ghz(), GHZ);
    println!("runs per figure:  {}", RUNS);

    /* A is a system parameter, generated once and outside every timed region. */
    let a: Matrix = genmatrix(&[0x42; SYMBYTES], true);

    let (sk, pk) = pswoosh_keygen(&a, true);

    let stats = [
        measure("Keygen", GHZ, RUNS, || pswoosh_keygen(&a, true)),
        measure("SdK", GHZ, RUNS, || pswoosh_skey_deriv(&pk, &pk, &sk, true)),
        measure("of which expand_seed_aes (x2 per Keygen)", GHZ, RUNS, || {
            let mut buf = [0u8; NOISE_BYTES];
            expand_seed_aes(&[0x5a; SYMBYTES], 0, &mut buf);
            buf[0]
        }),
    ];

    println!();
    for s in &stats {
        println!("{:<42} {:>14.0} cycles  {:>10.1} us", s.label, s.median, s.micros(GHZ));
    }
    println!();
    println!("|pk| {} bytes, |sk| {} bytes, |shk| {} bytes",
             PUBLICKEY_BYTES, SECRETKEY_BYTES, SYMBYTES);
}
